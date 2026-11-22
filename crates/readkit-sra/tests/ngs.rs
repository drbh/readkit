//! The NCBI VDB backend on real runs (feature `ngs`).
//!
//! Runs live in `<workspace>/.data/sra/<run>/<run>` or `READKIT_SRA_DIR`.
//!
//! - Yeast (single-end, ~15–25 MB each): run by default when the files
//!   exist; otherwise each test prints why it skipped.
//! - Human (paired, 713 MB): only with `READKIT_SRA_HUMAN=1` (in the runs
//!   directory) or `READKIT_SRA_HUMAN=<path>`.
//! - Exact equality with `fasterq-dump --split-spot --stdout` (bases,
//!   qualities and names), streamed and compared in memory.

#![cfg(feature = "ngs")]

use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::Instant;

use readkit_batch::Pool;
use readkit_io::BatchSource;
use readkit_sra::{
    DecodeOptions, ErrorKind, Names, PairedDecoder, QualityKind, SegmentSelection, SegmentsDecoder,
    SingleDecoder, run_info,
};

/// Where local runs live: `READKIT_SRA_DIR`, else `<workspace>/.data/sra`.
fn sra_dir() -> PathBuf {
    std::env::var_os("READKIT_SRA_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../.data/sra"))
}

/// (run, spots) for the yeast runs.
const YEAST: [(&str, u64); 4] = [
    ("SRR12999707", 89_234),
    ("SRR12999708", 109_901),
    ("SRR12999711", 134_282),
    ("SRR12999712", 138_347),
];

fn yeast(run: &str) -> Option<PathBuf> {
    let p = sra_dir().join(run).join(run);
    if p.exists() {
        Some(p)
    } else {
        eprintln!("skipping: {} not found", p.display());
        None
    }
}

fn human() -> Option<PathBuf> {
    let Ok(v) = std::env::var("READKIT_SRA_HUMAN") else {
        eprintln!("skipping: set READKIT_SRA_HUMAN=1 (or a path) for the paired human run");
        return None;
    };
    let p = if v == "1" {
        sra_dir().join("SRR13638690/SRR13638690")
    } else {
        PathBuf::from(v)
    };
    if p.exists() {
        Some(p)
    } else {
        eprintln!("skipping: {} not found", p.display());
        None
    }
}

const FASTERQ_DUMP: &str = "fasterq-dump";

/// Whether `program` is on `PATH`.
fn on_path(program: &str) -> bool {
    std::env::var_os("PATH")
        .is_some_and(|path| std::env::split_paths(&path).any(|dir| dir.join(program).is_file()))
}

/// One FASTQ record: (header without `@`, seq, qual).
type Record = (Vec<u8>, Vec<u8>, Vec<u8>);

/// Test-only FASTQ reader for fasterq-dump's 4-line records.
struct TinyFastq<R> {
    r: R,
    line: Vec<u8>,
}

impl<R: BufRead> TinyFastq<R> {
    fn line(&mut self) -> Option<Vec<u8>> {
        self.line.clear();
        if self.r.read_until(b'\n', &mut self.line).unwrap() == 0 {
            return None;
        }
        assert_eq!(self.line.pop(), Some(b'\n'), "unterminated line");
        Some(self.line.clone())
    }

    fn next(&mut self) -> Option<Record> {
        let head = self.line()?;
        let seq = self.line().expect("truncated record");
        let plus = self.line().expect("truncated record");
        let qual = self.line().expect("truncated record");
        assert_eq!(head[0], b'@');
        assert_eq!(plus[0], b'+');
        Some((head[1..].to_vec(), seq, qual))
    }
}

/// Stream `fasterq-dump <mode> --stdout` of `path` and compare it record by
/// record with `decoder` (see [`compare_with`]). Returns (records, seconds).
fn compare(path: &Path, mode: &[&str], decoder: SegmentsDecoder) -> (u64, f64) {
    static SEQ: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let seq = SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let tmp = std::env::temp_dir().join(format!("readkit-sra-fqd-{}-{seq}", std::process::id()));
    std::fs::create_dir_all(&tmp).unwrap();
    let _cleanup = Cleanup(tmp.clone());
    let mut cmd = Command::new(FASTERQ_DUMP);
    cmd.args(mode)
        .args(["--stdout", "-e", "4", "-t"])
        .arg(&tmp)
        .arg(path)
        .current_dir(&tmp);
    compare_with(cmd, decoder)
}

struct Cleanup(PathBuf);

impl Drop for Cleanup {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// Run `cmd`, read its FASTQ stdout and compare it record by record with
/// `decoder`: names, bases and qualities must be identical. Empty segments
/// are skipped, because the NCBI tools do not print them.
fn compare_with(mut cmd: Command, mut decoder: SegmentsDecoder) -> (u64, f64) {
    let t = Instant::now();
    let mut child = cmd
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .expect("starting the reference dumper");
    let mut fq = TinyFastq {
        r: BufReader::with_capacity(1 << 20, child.stdout.take().unwrap()),
        line: Vec::new(),
    };
    let mut n = 0u64;
    while let Some(batch) = decoder.next_batch().unwrap() {
        for spot in batch.spots() {
            for read in spot.iter().filter(|r| !r.is_empty()) {
                let (name, seq, qual) = fq
                    .next()
                    .unwrap_or_else(|| panic!("reference ended early after {n} records"));
                if name != read.name().unwrap() || seq != read.seq() || qual != read.qual() {
                    panic!(
                        "record {n} differs:\n reference {}\n readkit   {}",
                        String::from_utf8_lossy(&name),
                        String::from_utf8_lossy(read.name().unwrap())
                    );
                }
                n += 1;
            }
        }
    }
    assert!(fq.next().is_none(), "reference has records beyond {n}");
    assert!(child.wait().unwrap().success());
    (n, t.elapsed().as_secs_f64())
}

fn have_fasterq_dump() -> bool {
    let ok = on_path(FASTERQ_DUMP);
    if !ok {
        eprintln!("skipping: {FASTERQ_DUMP} not on PATH");
    }
    ok
}

#[test]
fn yeast_spot_counts_and_quality() {
    for (run, spots) in YEAST {
        let Some(p) = yeast(run) else { continue };
        let info = run_info(&p).unwrap();
        assert_eq!(info.spots, spots, "{run}");
        assert_eq!(info.first_spot, 1);
        assert_eq!(info.quality, QualityKind::Original, "{run}");
        let mut d = SingleDecoder::open(&p, DecodeOptions::default()).unwrap();
        assert_eq!(d.backend_name(), "ncbi-vdb");
        assert_eq!(d.info().spots, spots);
        let mut records = 0;
        while let Some(b) = d.next_batch().unwrap() {
            records += b.len() as u64;
        }
        assert_eq!(records, spots, "{run}");
        assert_eq!(d.stats().spots, spots);
    }
}

#[test]
fn yeast_equals_fasterq_dump_split_spot() {
    if !have_fasterq_dump() {
        return;
    }
    for (run, spots) in YEAST {
        let Some(p) = yeast(run) else { continue };
        let opts = DecodeOptions {
            names: Names::FasterqSplitSpot,
            ..Default::default()
        };
        let t = Instant::now();
        let d = SegmentsDecoder::open(&p, opts).unwrap();
        let (n, secs) = compare(&p, &["--split-spot"], d);
        assert_eq!(n, spots, "{run}: one biological read per spot");
        eprintln!(
            "{run}: {n} records identical to fasterq-dump ({:.2} s for both)",
            secs.max(t.elapsed().as_secs_f64())
        );
    }
}

#[test]
fn yeast_default_names_equal_fasterq_dump_default_defline() {
    if !have_fasterq_dump() {
        return;
    }
    let Some(p) = yeast("SRR12999707") else {
        return;
    };
    let opts = DecodeOptions {
        names: Names::Fasterq,
        ..Default::default()
    };
    // fasterq-dump refuses --split-files with --stdout, so ask --split-spot
    // for its default (split-3 / split-files) defline explicitly.
    let mode = ["--split-spot", "--seq-defline", "@$ac.$si $sn length=$rl"];
    let (n, _) = compare(&p, &mode, SegmentsDecoder::open(&p, opts).unwrap());
    assert_eq!(n, 89_234);
}

#[test]
fn yeast_read_structure_and_technical_reads() {
    let Some(p) = yeast("SRR12999707") else {
        return;
    };
    // The run stores two reads per spot: a biological one and an empty
    // technical one.
    let all = SegmentsDecoder::open(
        &p,
        DecodeOptions {
            segments: SegmentSelection::All,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(all.segments(), 2);
    assert!(all.read_types()[0].is_biological());
    assert!(all.read_types()[1].is_technical());
    let mut all = all;
    let b = all.next_batch().unwrap().unwrap();
    assert!(b.segment(1).seqs().data().is_empty());
    assert_eq!(b.segment(0).seqs().get(0).len(), 301);

    let e = PairedDecoder::open(&p, DecodeOptions::default()).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Layout);
}

#[test]
fn yeast_threads_match_sequential() {
    let Some(p) = yeast("SRR12999708") else {
        return;
    };
    let digest = |threads| {
        let mut d = SingleDecoder::open(
            &p,
            DecodeOptions {
                threads,
                batch_spots: 5_000,
                ..Default::default()
            },
        )
        .unwrap();
        let mut h = std::collections::hash_map::DefaultHasher::new();
        let mut n = 0;
        while let Some(b) = d.next_batch().unwrap() {
            use std::hash::Hash;
            let c = b.column();
            c.seqs().offsets().hash(&mut h);
            c.seqs().data().hash(&mut h);
            c.quals().data().hash(&mut h);
            n += b.len();
        }
        (n, std::hash::Hasher::finish(&h))
    };
    assert_eq!(digest(1), digest(6));
}

#[test]
fn yeast_pool_reuse_across_batches() {
    let Some(p) = yeast("SRR12999711") else {
        return;
    };
    let pool = Pool::new();
    let mut d = SingleDecoder::open(
        &p,
        DecodeOptions {
            threads: 1,
            batch_spots: 4_096,
            ..Default::default()
        },
    )
    .unwrap()
    .with_pool(&pool);
    drop(d.next_batch().unwrap().unwrap());
    let fresh = pool.stats().fresh;
    let mut batches = 1;
    while let Some(b) = d.next_batch().unwrap() {
        batches += 1;
        drop(b);
    }
    assert!(batches > 30);
    assert_eq!(
        pool.stats().fresh,
        fresh,
        "no new buffers after the first batch"
    );
    assert!(pool.stats().reused as usize >= (batches - 1) * 3);
}

#[test]
fn bad_inputs_are_typed_errors() {
    let e = SingleDecoder::open("/nonexistent/SRR1", DecodeOptions::default()).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::NotFound);
    // Not an SRA file at all.
    let e = SingleDecoder::open(
        concat!(env!("CARGO_MANIFEST_DIR"), "/Cargo.toml"),
        DecodeOptions::default(),
    )
    .unwrap_err();
    assert!(e.code().is_some(), "{e}");
    eprintln!("non-SRA file: {e}");
}

const FASTQ_DUMP: &str = "fastq-dump";

/// The paired human run: spot count, a full decode (throughput), and exact
/// equality with `fastq-dump --split-spot` on spot windows at the start,
/// middle and end. `fastq-dump` streams a spot range directly; a full
/// `fasterq-dump --stdout` of this run needs several GB of temporary files,
/// more than this machine has free.
#[test]
fn human_paired_equals_ncbi_on_windows() {
    let Some(p) = human() else { return };
    let info = run_info(&p).unwrap();
    assert_eq!(info.spots, 4_451_399);
    assert_eq!(info.quality, QualityKind::Original);

    for threads in [1, 4, 8] {
        let t = Instant::now();
        let mut d = PairedDecoder::open(
            &p,
            DecodeOptions {
                threads,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(d.read_types().iter().all(|t| t.is_biological()));
        let mut bases = 0u64;
        while let Some(b) = d.next_batch().unwrap() {
            bases += (b.r1().seqs().data().len() + b.r2().seqs().data().len()) as u64;
        }
        let secs = t.elapsed().as_secs_f64();
        assert_eq!(d.stats().spots, info.spots);
        eprintln!(
            "decode, {threads} threads: {} spots, {bases} bases in {secs:.2} s ({:.0} Mbases/s)",
            info.spots,
            bases as f64 / secs / 1e6
        );
    }

    if !on_path(FASTQ_DUMP) {
        eprintln!("skipping comparison: {FASTQ_DUMP} not on PATH");
        return;
    }
    for (first, count) in [(1u64, 100_000u64), (2_200_000, 50_000), (4_401_400, 50_000)] {
        let last = first + count - 1;
        let mut cmd = Command::new(FASTQ_DUMP);
        cmd.args(["--split-spot", "--stdout", "-N"])
            .arg(first.to_string())
            .arg("-X")
            .arg(last.to_string())
            .arg(&p);
        let opts = DecodeOptions {
            names: Names::Fasterq,
            spot_range: Some(first..last + 1),
            ..Default::default()
        };
        let (n, secs) = compare_with(cmd, SegmentsDecoder::open(&p, opts).unwrap());
        eprintln!("spots {first}..={last}: {n} records identical to fastq-dump ({secs:.1} s)");
        assert!(n > count, "{n} records for {count} paired spots");
    }
}
