//! The `fasterq-dump` backend.
//!
//! - Offline (default run, Unix): fake `fasterq-dump` / `vdb-dump` shell
//!   scripts written to a temp dir print fixture FASTQ or fail, covering
//!   spot grouping, the command line, empty spots, spot ranges, error
//!   mapping and child cleanup.
//! - Real tool: yeast runs (when present) against the VDB backend (feature
//!   `ngs`), and early-stop cleanup. Skipped with a message when
//!   `fasterq-dump` or the runs are missing.

#![cfg(unix)]

use std::os::unix::fs::PermissionsExt as _;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{Duration, Instant};

use readkit_io::BatchSource;
use readkit_sra::{
    DecodeOptions, ErrorKind, FasterqDumpBackend, LitePolicy, Names, OpenRequest, PairedDecoder,
    SingleDecoder, SpotReader as _,
};

/// A scratch dir removed on drop.
struct Dir(PathBuf);

impl Dir {
    fn new(tag: &str) -> Self {
        static SEQ: AtomicUsize = AtomicUsize::new(0);
        let n = SEQ.fetch_add(1, Ordering::Relaxed);
        let p =
            std::env::temp_dir().join(format!("readkit-sra-fq-{tag}-{}-{n}", std::process::id()));
        std::fs::create_dir_all(&p).unwrap();
        Self(p)
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// Fake SRA Toolkit: `bin/fasterq-dump` and `bin/vdb-dump` scripts, a temp
/// root, and a dummy run file.
struct Fake {
    dir: Dir,
}

impl Fake {
    /// `vdb` is the body of the vdb-dump script, `fqd` that of the
    /// fasterq-dump script. `$D` is the fake's dir in both.
    fn new(vdb: &str, fqd: &str) -> Self {
        let dir = Dir::new("fake");
        let d = &dir.0;
        std::fs::create_dir_all(d.join("bin")).unwrap();
        std::fs::create_dir_all(d.join("tmp")).unwrap();
        std::fs::write(d.join("SRR1"), b"not really a run").unwrap();
        let write = |name: &str, body: &str| {
            let p = d.join("bin").join(name);
            let text = format!(
                "#!/bin/sh\nD='{}'\nfor a in \"$@\"; do printf '%s\\n' \"$a\"; done > \"$D/{name}.args\"\n{body}\n",
                d.display()
            );
            std::fs::write(&p, text).unwrap();
            std::fs::set_permissions(&p, std::fs::Permissions::from_mode(0o755)).unwrap();
        };
        write("vdb-dump", vdb);
        write("fasterq-dump", fqd);
        Self { dir }
    }

    /// A fake whose vdb-dump reports `spots` spots and whose fasterq-dump
    /// prints `fastq` and exits 0.
    fn printing(spots: u64, fastq: &str) -> Self {
        let f = Self::new(&id_range(spots), "cat \"$D/out.fq\"");
        std::fs::write(f.dir.0.join("out.fq"), fastq).unwrap();
        f
    }

    fn backend(&self) -> FasterqDumpBackend {
        FasterqDumpBackend::new()
            .program(self.dir.0.join("bin/fasterq-dump"))
            .temp_dir(self.dir.0.join("tmp"))
    }

    fn run(&self) -> PathBuf {
        self.dir.0.join("SRR1")
    }

    fn args(&self, tool: &str) -> Vec<String> {
        std::fs::read_to_string(self.dir.0.join(format!("{tool}.args")))
            .unwrap_or_default()
            .lines()
            .map(str::to_owned)
            .collect()
    }

    /// Entries left in the temp root.
    fn leftovers(&self) -> Vec<String> {
        std::fs::read_dir(self.dir.0.join("tmp"))
            .unwrap()
            .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
            .collect()
    }
}

fn id_range(spots: u64) -> String {
    format!("echo 'id-range: first-row = 1, row-count = {spots}'")
}

/// FASTQ text; qualities vary so the Lite probe sees original qualities.
fn fq(records: &[(&str, &str)]) -> String {
    let quals = b"I5?+#AE";
    let mut s = String::new();
    for (head, seq) in records {
        let q: String = (0..seq.len())
            .map(|i| quals[i % quals.len()] as char)
            .collect();
        s += &format!("@{head}\n{seq}\n+\n{q}\n");
    }
    s
}

fn qual(len: usize) -> Vec<u8> {
    let quals = b"I5?+#AE";
    (0..len).map(|i| quals[i % quals.len()]).collect()
}

fn opts(names: Names) -> DecodeOptions {
    DecodeOptions {
        names,
        accession: Some("ACC".into()),
        ..Default::default()
    }
}

#[test]
fn paired_records_are_grouped_into_spots() {
    let fake = Fake::printing(
        3,
        &fq(&[
            ("1 n1", "ACGT"),
            ("1 n1", "TT"),
            ("2 n2", "G"),
            ("2 n2", "CCCAA"),
            ("3 n3", "A"),
            ("3 n3", "C"),
        ]),
    );
    let backend = fake
        .backend()
        .threads(3)
        .include_technical(true)
        .args(["-m", "50MB"]);
    let mut d =
        PairedDecoder::open_with(Arc::new(backend), fake.run(), opts(Names::FasterqSplitSpot))
            .unwrap();
    assert_eq!(d.backend_name(), "fasterq-dump");
    assert_eq!(d.info().spots, 3);
    assert!(d.read_types().iter().all(|t| t.is_biological()));
    let b = d.next_batch().unwrap().unwrap();
    assert!(d.next_batch().unwrap().is_none());
    let got: Vec<_> = b
        .pairs()
        .map(|(r1, r2)| {
            (
                String::from_utf8(r1.name().unwrap().to_vec()).unwrap(),
                r1.seq().to_vec(),
                r2.seq().to_vec(),
                r2.qual().to_vec(),
            )
        })
        .collect();
    assert_eq!(
        got,
        vec![
            (
                "ACC.1/1 n1 length=4".into(),
                b"ACGT".to_vec(),
                b"TT".to_vec(),
                qual(2)
            ),
            (
                "ACC.2/1 n2 length=1".into(),
                b"G".to_vec(),
                b"CCCAA".to_vec(),
                qual(5)
            ),
            (
                "ACC.3/1 n3 length=1".into(),
                b"A".to_vec(),
                b"C".to_vec(),
                qual(1)
            ),
        ]
    );
    assert_eq!(b.r2().names().unwrap().get(1), b"ACC.2/2 n2 length=5");

    let tmp = fake.dir.0.join("tmp");
    let args = fake.args("fasterq-dump");
    let t = args.iter().position(|a| a == "-t").unwrap();
    assert!(
        Path::new(&args[t + 1]).starts_with(&tmp) && args[t + 1].contains("readkit-fasterq-"),
        "{args:?}"
    );
    let mut want: Vec<String> = [
        "--split-spot",
        "--stdout",
        "--seq-defline",
        "@$si $sn",
        "--qual-defline",
        "+",
        "-t",
    ]
    .map(String::from)
    .to_vec();
    want.push(args[t + 1].clone());
    want.extend(["-e", "3", "--include-technical", "-m", "50MB"].map(String::from));
    want.push(fake.run().display().to_string());
    assert_eq!(args, want);
    assert_eq!(
        fake.args("vdb-dump"),
        vec!["--id_range".to_string(), fake.run().display().to_string()]
    );
    drop(d);
    assert!(fake.leftovers().is_empty(), "{:?}", fake.leftovers());
}

#[test]
fn without_names_the_defline_is_the_spot_id() {
    let fake = Fake::printing(2, &fq(&[("1", "ACGT"), ("2", "GTCA")]));
    let mut d =
        SingleDecoder::open_with(Arc::new(fake.backend()), fake.run(), opts(Names::None)).unwrap();
    let b = d.next_batch().unwrap().unwrap();
    assert!(!b.has_names());
    assert_eq!(b.column().seqs().data(), b"ACGTGTCA");
    let args = fake.args("fasterq-dump");
    assert!(args.windows(2).any(|w| w == ["--seq-defline", "@$si"]));
    assert!(!args.iter().any(|a| a == "-e" || a == "--include-technical"));
}

/// The tool runs in its temp dir, where a relative run path does not exist
/// and a bare one is taken for an accession to look up.
#[test]
fn a_relative_path_reaches_the_tool_absolute() {
    let fake = Fake::printing(1, &fq(&[("1", "ACGT")]));
    // Tests run in the package directory; the fake tools never read the run.
    let mut d = SingleDecoder::open_with(Arc::new(fake.backend()), "Cargo.toml", opts(Names::None))
        .unwrap();
    assert_eq!(d.next_batch().unwrap().unwrap().len(), 1);
    let expected = std::env::current_dir().unwrap().join("Cargo.toml");
    for tool in ["vdb-dump", "fasterq-dump"] {
        let args = fake.args(tool);
        assert_eq!(
            args.last().map(Path::new),
            Some(expected.as_path()),
            "{tool}"
        );
    }
}

#[test]
fn spots_without_records_are_empty_reads() {
    // Spot 2 and the last spot 4 print nothing (empty reads).
    let fake = Fake::printing(4, &fq(&[("1 a", "ACGT"), ("3 c", "GG")]));
    let mut d =
        SingleDecoder::open_with(Arc::new(fake.backend()), fake.run(), opts(Names::Fasterq))
            .unwrap();
    let b = d.next_batch().unwrap().unwrap();
    let seqs: Vec<&[u8]> = b.column().seqs().iter().collect();
    assert_eq!(seqs, vec![&b"ACGT"[..], b"", b"GG", b""]);
    let names: Vec<&[u8]> = b.column().names().unwrap().iter().collect();
    assert_eq!(names[1], b"ACC.2 2 length=0");
    assert_eq!(d.stats().spots, 4);
}

#[test]
fn spot_ranges_skip_printed_spots() {
    let fake = Fake::printing(
        4,
        &fq(&[("1", "A"), ("2", "CC"), ("3", "GGG"), ("4", "TTTT")]),
    );
    let o = DecodeOptions {
        spot_range: Some(2..4),
        ..opts(Names::None)
    };
    let mut d = SingleDecoder::open_with(Arc::new(fake.backend()), fake.run(), o).unwrap();
    let b = d.next_batch().unwrap().unwrap();
    assert_eq!(b.column().seqs().data(), b"CCGGG");
    assert!(d.next_batch().unwrap().is_none());
}

#[test]
fn quality_is_unknown_so_lite_is_sampled() {
    let records: Vec<(String, String)> = (1..=3)
        .map(|i| (i.to_string(), "ACGT".to_string()))
        .collect();
    let text: String = records
        .iter()
        .map(|(h, s)| format!("@{h}\n{s}\n+\n????\n"))
        .collect();
    let fake = Fake::printing(3, &text);
    let mut d =
        SingleDecoder::open_with(Arc::new(fake.backend()), fake.run(), opts(Names::None)).unwrap();
    assert_eq!(d.quality(), readkit_sra::QualityKind::Unknown);
    let e = d.next_batch().unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Lite, "{e}");

    let o = DecodeOptions {
        lite: LitePolicy::Allow,
        ..opts(Names::None)
    };
    let mut d = SingleDecoder::open_with(Arc::new(fake.backend()), fake.run(), o).unwrap();
    assert_eq!(d.next_batch().unwrap().unwrap().len(), 3);
}

/// Open a single-end decoder on `fake` and decode to the end.
fn decode_all(fake: &Fake) -> Result<u64, readkit_sra::Error> {
    let mut d = SingleDecoder::open_with(Arc::new(fake.backend()), fake.run(), opts(Names::None))?;
    while d.next_batch()?.is_some() {}
    Ok(d.stats().spots)
}

#[test]
fn tool_failures_are_typed() {
    let cases = [
        (
            "echo \"$(date) fasterq-dump err: cannot find 'x' -> RC(rcVFS,rcQuery,rcResolving,rcName,rcNotFound)\" >&2; exit 3",
            ErrorKind::NotFound,
            "rcNotFound",
        ),
        (
            "echo 'err: file invalid' >&2; echo 'err: VCursorOpen -> RC(rcFS,rcDirectory,rcOpening,rcFile,rcInvalid)' >&2; exit 3",
            ErrorKind::Corrupt,
            "file invalid",
        ),
        (
            "echo 'err: the input data is missing the QUALITY-column' >&2; exit 3",
            ErrorKind::Unsupported,
            "QUALITY-column",
        ),
        ("echo boom >&2; exit 2", ErrorKind::Backend, "boom"),
        ("kill -9 $$", ErrorKind::Backend, "signal"),
    ];
    for (body, kind, needle) in cases {
        let fake = Fake::new(&id_range(2), body);
        let e = decode_all(&fake).unwrap_err();
        assert_eq!(e.kind(), kind, "{body}: {e}");
        assert!(e.to_string().contains(needle), "{e}");
        assert!(e.to_string().contains("fasterq-dump"), "{e}");
        assert!(fake.leftovers().is_empty());
    }
}

#[test]
fn a_failure_after_all_records_is_still_an_error() {
    let fake = Fake::new(
        &id_range(2),
        "cat \"$D/out.fq\"; echo 'late failure' >&2; exit 3",
    );
    std::fs::write(fake.dir.0.join("out.fq"), fq(&[("1", "A"), ("2", "C")])).unwrap();
    let e = decode_all(&fake).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Backend);
    assert!(e.message().contains("late failure"), "{e}");
}

#[test]
fn spot_counting_failures_are_typed() {
    let cases = [
        (
            "echo \"vdb-dump int: item not found while constructing - the path 'x' cannot be opened as vdb-database or vdb-table\" >&2; exit 3",
            ErrorKind::Unsupported,
        ),
        ("echo 'something else'", ErrorKind::Backend),
    ];
    for (body, kind) in cases {
        let fake = Fake::new(body, "exit 0");
        let e = decode_all(&fake).unwrap_err();
        assert_eq!(e.kind(), kind, "{e}");
        assert!(e.to_string().contains("vdb-dump"), "{e}");
        // fasterq-dump never ran.
        assert!(fake.args("fasterq-dump").is_empty());
    }
}

#[test]
fn bad_output_is_corrupt() {
    let cases = [
        // Sequence and quality lengths differ.
        ("@1\nACGT\n+\nII\n".to_string(), 1, "length"),
        // Spots out of order.
        (fq(&[("2", "A"), ("1", "C")]), 2, "after spot 2"),
        // Reads of one spot not consecutive.
        (fq(&[("1", "A"), ("2", "C"), ("1", "G")]), 2, "after spot 2"),
        // More spots than vdb-dump counted.
        (fq(&[("1", "A"), ("2", "C")]), 1, "outside the run"),
        // No spot id in the header.
        (fq(&[("x1", "A")]), 1, "spot id"),
        // Truncated record from a tool that exited 0.
        ("@1\nACGT\n+\n".to_string(), 1, "truncated"),
    ];
    for (text, spots, needle) in cases {
        let fake = Fake::printing(spots, &text);
        let e = decode_all(&fake).unwrap_err();
        assert_eq!(e.kind(), ErrorKind::Corrupt, "{text:?}: {e}");
        assert!(e.to_string().contains(needle), "{needle}: {e}");
    }
}

#[test]
fn a_missing_program_is_unsupported() {
    let fake = Fake::printing(1, &fq(&[("1", "A")]));
    let backend = fake
        .backend()
        .program(fake.dir.0.join("bin/no-such-dumper"));
    // vdb-dump is not next to the missing program: name it explicitly.
    let backend = backend.vdb_dump(fake.dir.0.join("bin/vdb-dump"));
    let e = SingleDecoder::open_with(Arc::new(backend), fake.run(), opts(Names::None)).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Unsupported, "{e}");
    assert!(
        !FasterqDumpBackend::new()
            .program(fake.dir.0.join("bin/no-such-dumper"))
            .is_available()
    );
    assert!(fake.backend().is_available());

    let e = SingleDecoder::open_with(
        Arc::new(fake.backend()),
        fake.dir.0.join("no-such-run"),
        opts(Names::None),
    )
    .unwrap_err();
    assert_eq!(e.kind(), ErrorKind::NotFound);
}

/// Whether process `pid` exists (zombies included).
fn alive(pid: u32) -> bool {
    let out = Command::new("ps")
        .args(["-o", "stat=", "-p", &pid.to_string()])
        .output()
        .unwrap();
    !String::from_utf8_lossy(&out.stdout).trim().is_empty()
}

/// Wait up to 5 s for `pid` to disappear (an orphaned grandchild is reaped
/// by init, not by us).
fn gone(pid: u32) -> bool {
    let t = Instant::now();
    while alive(pid) {
        if t.elapsed() > Duration::from_secs(5) {
            return false;
        }
        std::thread::sleep(Duration::from_millis(20));
    }
    true
}

/// A fasterq-dump that forks a long-lived grandchild (like the SRA
/// Toolkit's driver) and prints spots forever.
fn endless() -> Fake {
    Fake::new(
        &id_range(1_000_000_000),
        "sleep 300 & echo $! > \"$D/grandchild.pid\"\n\
         i=1; while :; do printf '@%d\\nACGT\\n+\\nI5?+\\n' $i; i=$((i+1)); done",
    )
}

#[test]
fn dropping_a_reader_kills_the_process_group() {
    let fake = endless();
    let backend = fake.backend();
    let mut r = backend
        .open_reader(&fake.run(), &OpenRequest::new(false))
        .unwrap();
    assert!(r.child_pid().is_none(), "started lazily");
    assert!(!r.random_access());
    for spot in 1..=5000 {
        assert_eq!(r.read_spot(spot).unwrap().bases, b"ACGT");
    }
    // Asking again for the current spot is allowed; going back is not.
    assert_eq!(r.read_spot(5000).unwrap().lens, &[4]);
    assert_eq!(r.read_spot(10).unwrap_err().kind(), ErrorKind::Unsupported);
    let pid = r.child_pid().unwrap();
    let grandchild: u32 = std::fs::read_to_string(fake.dir.0.join("grandchild.pid"))
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    assert!(alive(pid) && alive(grandchild));
    let t = Instant::now();
    drop(r);
    assert!(t.elapsed() < Duration::from_secs(5));
    assert!(!alive(pid), "child {pid} not reaped");
    assert!(gone(grandchild), "grandchild {grandchild} survived");
    assert!(fake.leftovers().is_empty(), "{:?}", fake.leftovers());
}

#[test]
fn dropping_a_decoder_early_stops_the_child() {
    let fake = endless();
    let o = DecodeOptions {
        batch_spots: 1000,
        ..opts(Names::Fasterq)
    };
    let mut d = SingleDecoder::open_with(Arc::new(fake.backend()), fake.run(), o).unwrap();
    assert_eq!(d.next_batch().unwrap().unwrap().len(), 1000);
    let grandchild: u32 = std::fs::read_to_string(fake.dir.0.join("grandchild.pid"))
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    drop(d);
    assert!(gone(grandchild));
    assert!(fake.leftovers().is_empty());
    assert!(!ps_mentions(&fake.dir.0));
}

/// Whether any process's command line mentions `needle`.
fn ps_mentions(needle: &Path) -> bool {
    let out = Command::new("ps")
        .args(["-axo", "command"])
        .output()
        .unwrap();
    String::from_utf8_lossy(&out.stdout)
        .lines()
        .any(|l| l.contains(&*needle.to_string_lossy()))
}

// ---------------------------------------------------------------------------
// The real tool on the yeast runs.

/// Where local runs live: `READKIT_SRA_DIR`, else `<workspace>/.data/sra`.
fn sra_dir() -> PathBuf {
    std::env::var_os("READKIT_SRA_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../.data/sra"))
}

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

/// The real backend with its temp files under a fresh dir: inside
/// `READKIT_SRA_TMP`, else the system temp dir.
fn real() -> Option<(FasterqDumpBackend, Dir)> {
    let backend = FasterqDumpBackend::new();
    if !backend.is_available() {
        eprintln!("skipping: fasterq-dump / vdb-dump not on PATH");
        return None;
    }
    let dir = match std::env::var_os("READKIT_SRA_TMP") {
        Some(root) => {
            static SEQ: AtomicUsize = AtomicUsize::new(0);
            let n = SEQ.fetch_add(1, Ordering::Relaxed);
            let d =
                Path::new(&root).join(format!("readkit-sra-fq-real-{}-{n}", std::process::id()));
            std::fs::create_dir_all(&d).unwrap();
            Dir(d)
        }
        None => Dir::new("real"),
    };
    Some((backend.threads(4).temp_dir(&dir.0), dir))
}

#[test]
fn real_yeast_spot_counts() {
    let Some((backend, _dir)) = real() else {
        return;
    };
    for (run, spots) in YEAST {
        let Some(p) = yeast(run) else { continue };
        let r = backend.open_reader(&p, &OpenRequest::new(false)).unwrap();
        assert_eq!((r.info().first_spot, r.info().spots), (1, spots), "{run}");
        assert!(r.child_pid().is_none());
    }
}

#[test]
fn real_early_stop_leaves_no_child() {
    let Some((backend, dir)) = real() else {
        return;
    };
    let Some(p) = yeast("SRR12999712") else {
        return;
    };
    let mut r = backend.open_reader(&p, &OpenRequest::new(true)).unwrap();
    assert_eq!(r.read_spot(1).unwrap().lens.len(), 1);
    let pid = r.child_pid().unwrap();
    drop(r);
    assert!(!alive(pid), "fasterq-dump {pid} left behind");

    let o = DecodeOptions {
        batch_spots: 1000,
        ..Default::default()
    };
    let mut d = SingleDecoder::open_with(Arc::new(backend), &p, o).unwrap();
    assert_eq!(d.next_batch().unwrap().unwrap().len(), 1000);
    drop(d);
    let t = Instant::now();
    while ps_mentions(&dir.0) {
        assert!(
            t.elapsed() < Duration::from_secs(5),
            "a fasterq-dump using {} is still running",
            dir.0.display()
        );
        std::thread::sleep(Duration::from_millis(20));
    }
    let left: Vec<_> = std::fs::read_dir(&dir.0).unwrap().collect();
    assert!(left.is_empty(), "{left:?}");
}

/// Exact equality with the VDB backend: batch boundaries, bases,
/// qualities, names in both fasterq-dump styles, and spot counts.
#[cfg(feature = "ngs")]
#[test]
fn real_yeast_equals_ngs() {
    use readkit_sra::{NgsBackend, SegmentsDecoder};
    let Some((backend, dir)) = real() else {
        return;
    };
    let backend = Arc::new(backend);
    for (run, spots) in YEAST {
        let Some(p) = yeast(run) else { continue };
        for names in [Names::None, Names::FasterqSplitSpot, Names::Fasterq] {
            let o = DecodeOptions {
                names,
                ..Default::default()
            };
            let mut a = SegmentsDecoder::open_with(backend.clone(), &p, o.clone()).unwrap();
            let mut b = SegmentsDecoder::open_with(Arc::new(NgsBackend::new()), &p, o).unwrap();
            assert_eq!(a.segments(), b.segments());
            let mut batches = 0;
            loop {
                let (x, y) = (a.next_batch().unwrap(), b.next_batch().unwrap());
                let (x, y) = match (x, y) {
                    (Some(x), Some(y)) => (x, y),
                    (None, None) => break,
                    _ => panic!("{run}: batch counts differ"),
                };
                assert_eq!(x.len(), y.len());
                for (cx, cy) in x.columns().zip(y.columns()) {
                    assert_eq!(cx.seqs().offsets(), cy.seqs().offsets(), "{run}");
                    assert!(cx.seqs().data() == cy.seqs().data(), "{run}: bases differ");
                    assert!(
                        cx.quals().data() == cy.quals().data(),
                        "{run}: quals differ"
                    );
                    match (cx.names(), cy.names()) {
                        (None, None) => assert_eq!(names, Names::None),
                        (Some(nx), Some(ny)) => {
                            assert_eq!(nx.offsets(), ny.offsets());
                            assert!(nx.data() == ny.data(), "{run}: names differ");
                        }
                        _ => panic!("{run}: names present on one side only"),
                    }
                }
                batches += 1;
            }
            assert_eq!(a.stats(), b.stats());
            assert_eq!(a.stats().spots, spots);
            assert_eq!(a.quality(), readkit_sra::QualityKind::Original);
            eprintln!("{run} {names:?}: {batches} batches, {spots} spots identical");
        }
    }
    let left: Vec<_> = std::fs::read_dir(&dir.0).unwrap().collect();
    assert!(left.is_empty(), "{left:?}");
}
