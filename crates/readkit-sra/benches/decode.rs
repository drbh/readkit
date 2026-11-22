//! `.sra` decoding benchmarks on pinned real runs (feature `ngs`).
//!
//! ```text
//! cargo bench -p readkit-sra --features ngs --bench decode -- --fetch     # download missing runs once
//! cargo bench -p readkit-sra --features ngs --bench decode                # every run and mode
//! cargo bench -p readkit-sra --features ngs --bench decode -- --run SRR10971381 --mode t1
//! cargo bench -p readkit-sra --features ngs --bench decode -- --save main
//! cargo bench -p readkit-sra --features ngs --bench decode -- --compare main
//! cargo bench -p readkit-sra --features ngs --bench decode -- --reference # also fastq-dump
//! ```
//!
//! Runs live in `<workspace>/.data/sra/<run>/<run>` (or `READKIT_SRA_DIR`).
//! They are pinned by size and S3 ETag; `--fetch` downloads missing ones from
//! NCBI's bucket and checks the pin. A run with more than `--spots` spots
//! (default 1 000 000) is decoded over its first `--spots` only, so runs of
//! any size cost about the same.
//!
//! Modes:
//! - `read`: the backend alone, one reader serving every spot
//!   (`SpotReader::read_spot`), nothing copied: the floor for one thread.
//! - `t1`, `t4`, `tN` (N = all cores): the full decoder into batches.
//! - `t4-names`: four threads with fasterq-dump names.
//! - `vdb-…`: the same mode with `NgsBackend::vdb_only()` (no fast path).
//! - `open`, `vdb-open`: opening a reader and reading its first spot, 20
//!   times (a cost every decoder thread pays once).
//! - `fastq-dump` (with `--reference`): `fastq-dump --split-spot
//!   --stdout -N -X` over the same spots, to /dev/null. It is the toolkit's
//!   single-threaded dumper; `fasterq-dump` cannot select a spot range.
//!
//! Results: a table (seconds, Mspots/s, Mbases/s) and JSON at
//! `target/sra-bench/latest.json`, also `bench/baselines/sra-NAME.json`
//! with `--save NAME`.

use readkit_io::BatchSource;
use readkit_sra::{DecodeOptions, Names, NgsBackend, OpenRequest, SegmentsDecoder, SraBackend};
use serde_json::{Value, json};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::Arc;
use std::time::Instant;

struct Run {
    acc: &'static str,
    about: &'static str,
    size: u64,
    etag: &'static str,
}

const RUNS: &[Run] = &[
    Run {
        acc: "SRR12999707",
        about: "yeast, 1×301 + empty technical read, 16 MB",
        size: 16_143_261,
        etag: "7279b48a19275f5644569ab0073dc0b1-2",
    },
    Run {
        acc: "SRR13638690",
        about: "human, 2×301, full qualities, 713 MB",
        size: 713_058_679,
        etag: "de419e86e94fc1f72632919e46b2cf20-86",
    },
    Run {
        acc: "SRR10971381",
        about: "2×151, binned qualities (NovaSeq-style), 2.7 GB",
        size: 2_743_429_212,
        etag: "bc4c1ce8d7e0d929c6e3ef4db0d48a28-328",
    },
];

const S3: &str = "https://sra-pub-run-odp.s3.amazonaws.com/sra";

struct Args {
    runs: Vec<String>,
    modes: Vec<String>,
    spots: u64,
    repeat: usize,
    fetch: bool,
    reference: bool,
    save: Option<String>,
    compare: Option<String>,
}

fn args() -> Args {
    let mut a = Args {
        runs: Vec::new(),
        modes: Vec::new(),
        spots: 1_000_000,
        repeat: 3,
        fetch: false,
        reference: false,
        save: None,
        compare: None,
    };
    let mut it = std::env::args().skip(1);
    while let Some(arg) = it.next() {
        let mut val = || it.next().unwrap_or_else(|| panic!("{arg} needs a value"));
        match arg.as_str() {
            "--run" => a.runs.push(val()),
            "--mode" => a.modes.push(val()),
            "--spots" => a.spots = val().parse().expect("--spots N"),
            "--repeat" => a.repeat = val().parse().expect("--repeat N"),
            "--fetch" => a.fetch = true,
            "--reference" => a.reference = true,
            "--save" => a.save = Some(val()),
            "--compare" => a.compare = Some(val()),
            "--bench" => {}
            _ => panic!("unknown argument {arg} (see the file's docs)"),
        }
    }
    a
}

fn workspace() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn data_dir() -> PathBuf {
    std::env::var_os("READKIT_SRA_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| workspace().join(".data/sra"))
}

fn run_path(r: &Run) -> PathBuf {
    data_dir().join(r.acc).join(r.acc)
}

fn fetch(r: &Run) {
    use readkit_fetch::{FetchOptions, Fetcher, FileSink};
    let fetcher = Fetcher::new(FetchOptions::default()).expect("fetcher");
    let object = fetcher
        .probe_s3(&format!("{S3}/{}/{}", r.acc, r.acc))
        .expect("probe");
    let etag = match &object.checksum {
        readkit_fetch::Checksum::S3MultipartEtag { etag, .. } => etag.clone(),
        c => format!("{c:?}"),
    };
    assert!(
        object.size == r.size && etag == r.etag,
        "{}: the bucket now has {} bytes, ETag {etag}; the bench pins {} bytes, {}",
        r.acc,
        object.size,
        r.size,
        r.etag
    );
    println!("fetching {} ({} MB)", r.acc, r.size / 1_000_000);
    let done = fetcher
        .fetch(object, FileSink::new(run_path(r)))
        .join()
        .expect("download");
    println!("  {:.1}s", done.report.elapsed.as_secs_f64());
}

struct Measure {
    secs: Vec<f64>,
    spots: u64,
    bases: u64,
}

fn main() {
    let a = args();
    let cores = std::thread::available_parallelism().map_or(4, |n| n.get());
    let mut modes = vec![
        "read".to_string(),
        "vdb-read".into(),
        "t1".into(),
        "vdb-t1".into(),
        "t4".into(),
        "vdb-t4".into(),
        format!("t{cores}"),
        format!("vdb-t{cores}"),
        "t4-names".into(),
        "open".into(),
        "vdb-open".into(),
    ];
    if a.reference {
        modes.push("fastq-dump".into());
    }
    let mut results = Vec::new();
    for r in RUNS {
        if !a.runs.is_empty() && !a.runs.iter().any(|x| x == r.acc) {
            continue;
        }
        let path = run_path(r);
        if !path.exists() {
            if a.fetch {
                fetch(r);
            } else {
                println!(
                    "skipping {}: {} not found (pass --fetch)",
                    r.acc,
                    path.display()
                );
                continue;
            }
        }
        let info = readkit_sra::run_info(&path).expect("run info");
        let range = info.first_spot..info.end_spot().min(info.first_spot + a.spots);
        println!(
            "== {}: {} — spots {}..{} of {}",
            r.acc, r.about, range.start, range.end, info.spots
        );
        for mode in &modes {
            if !a.modes.is_empty() && !a.modes.iter().any(|m| m == mode) {
                continue;
            }
            let m = Measure::of(a.repeat, || measure(mode, &path, range.clone()));
            let med = median(&m.secs);
            println!(
                "  {mode:10} {med:7.3}s {:7.2} Mspots/s {:7.0} Mbases/s | spread {:4.1}%",
                m.spots as f64 / med / 1e6,
                m.bases as f64 / med / 1e6,
                (max(&m.secs) - min(&m.secs)) / med * 100.0
            );
            results.push(json!({
                "run": r.acc,
                "mode": mode,
                "spots": m.spots,
                "bases": m.bases,
                "runs_s": m.secs,
                "median_s": med,
                "mspots_s": m.spots as f64 / med / 1e6,
                "mbases_s": m.bases as f64 / med / 1e6,
            }));
        }
    }
    report(&a, results);
}

impl Measure {
    fn of(repeat: usize, mut f: impl FnMut() -> (f64, u64, u64)) -> Measure {
        let mut m = Measure {
            secs: Vec::new(),
            spots: 0,
            bases: 0,
        };
        for _ in 0..repeat.max(1) {
            let (s, spots, bases) = f();
            m.secs.push(s);
            (m.spots, m.bases) = (spots, bases);
        }
        m
    }
}

/// One timed pass: (seconds, spots, bases).
fn measure(mode: &str, path: &Path, range: std::ops::Range<u64>) -> (f64, u64, u64) {
    let (backend, mode) = match mode.strip_prefix("vdb-") {
        Some(m) => (NgsBackend::vdb_only(), m),
        None => (NgsBackend::new(), mode),
    };
    let t = Instant::now();
    match mode {
        "open" => {
            for _ in 0..20 {
                let mut r = backend.open(path, &OpenRequest::new(false)).expect("open");
                std::hint::black_box(r.read_spot(range.start).expect("spot").bases.len());
            }
            (t.elapsed().as_secs_f64(), 20, 0)
        }
        "read" => {
            let mut r = backend.open(path, &OpenRequest::new(false)).expect("open");
            let mut bases = 0u64;
            let mut touch = 0u8;
            for spot in range.clone() {
                let s = r.read_spot(spot).expect("spot");
                bases += s.bases.len() as u64;
                touch ^= s.quals.first().copied().unwrap_or(0);
            }
            std::hint::black_box(touch);
            (t.elapsed().as_secs_f64(), range.end - range.start, bases)
        }
        "fastq-dump" => {
            let ok = Command::new("fastq-dump")
                .args(["--split-spot", "--stdout", "-N"])
                .arg(range.start.to_string())
                .arg("-X")
                .arg((range.end - 1).to_string())
                .arg(path)
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .status()
                .expect("fastq-dump on PATH")
                .success();
            assert!(ok, "fastq-dump failed");
            (t.elapsed().as_secs_f64(), range.end - range.start, 0)
        }
        _ => {
            let (threads, names) = match mode.strip_suffix("-names") {
                Some(m) => (m, Names::Fasterq),
                None => (mode, Names::None),
            };
            let threads = threads[1..].parse().expect("mode tN");
            let mut d = SegmentsDecoder::open_with(
                Arc::new(backend),
                path,
                DecodeOptions {
                    threads,
                    names,
                    spot_range: Some(range),
                    ..Default::default()
                },
            )
            .expect("open");
            while let Some(b) = d.next_batch().expect("batch") {
                std::hint::black_box(&b);
            }
            let s = d.stats();
            (t.elapsed().as_secs_f64(), s.spots, s.bases)
        }
    }
}

fn median(v: &[f64]) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(f64::total_cmp);
    s[s.len() / 2]
}

fn max(v: &[f64]) -> f64 {
    v.iter().cloned().fold(f64::MIN, f64::max)
}

fn min(v: &[f64]) -> f64 {
    v.iter().cloned().fold(f64::MAX, f64::min)
}

fn git(args: &[&str]) -> String {
    Command::new("git")
        .args(args)
        .current_dir(workspace())
        .output()
        .ok()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_default()
}

fn report(a: &Args, results: Vec<Value>) {
    let dirty = !git(&["status", "--porcelain", "--untracked-files=no"]).is_empty();
    let doc = json!({
        "schema": 1,
        "git": format!("{}{}", git(&["rev-parse", "--short", "HEAD"]), if dirty { "-dirty" } else { "" }),
        "unix_time": std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| d.as_secs()),
        "machine": {
            "os": std::env::consts::OS,
            "arch": std::env::consts::ARCH,
            "cpus": std::thread::available_parallelism().map_or(0, |n| n.get()),
        },
        "spots_limit": a.spots,
        "repeat": a.repeat,
        "results": results,
    });
    let text = serde_json::to_string_pretty(&doc).unwrap();
    let out = workspace().join("target/sra-bench");
    std::fs::create_dir_all(&out).unwrap();
    std::fs::write(out.join("latest.json"), &text).unwrap();
    if let Some(name) = &a.save {
        let dir = workspace().join("bench/baselines");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(format!("sra-{name}.json"));
        std::fs::write(&path, &text).unwrap();
        println!("saved {}", path.display());
    }
    if let Some(name) = &a.compare {
        let path = workspace().join(format!("bench/baselines/sra-{name}.json"));
        let base: Value =
            serde_json::from_slice(&std::fs::read(&path).expect("baseline")).expect("JSON");
        println!(
            "== against {name} ({})",
            base["git"].as_str().unwrap_or("?")
        );
        for r in doc["results"].as_array().unwrap() {
            let old = base["results"].as_array().into_iter().flatten().find(|b| {
                b["run"] == r["run"] && b["mode"] == r["mode"] && b["spots"] == r["spots"]
            });
            let new = r["median_s"].as_f64().unwrap();
            match old.and_then(|o| o["median_s"].as_f64()) {
                Some(old) => println!(
                    "  {:12} {:10} {:7.3}s -> {:7.3}s  {:+6.1}%",
                    r["run"].as_str().unwrap(),
                    r["mode"].as_str().unwrap(),
                    old,
                    new,
                    (new / old - 1.0) * 100.0
                ),
                None => println!(
                    "  {:12} {:10} new",
                    r["run"].as_str().unwrap(),
                    r["mode"].as_str().unwrap()
                ),
            }
        }
    }
}
