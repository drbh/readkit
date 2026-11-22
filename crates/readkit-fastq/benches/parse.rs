//! FASTQ parsing benchmarks: the parser alone and gzip → parser, on a
//! deterministic synthetic paired dataset (or your own files).
//!
//! ```text
//! cargo bench -p readkit-fastq --bench parse                    # every mode
//! cargo bench -p readkit-fastq --bench parse -- --mode text-t4 --mode gz-t4
//! cargo bench -p readkit-fastq --bench parse -- --save main
//! cargo bench -p readkit-fastq --bench parse -- --compare main
//! cargo bench -p readkit-fastq --bench parse -- --input r1.fq.gz r2.fq.gz --label srr
//! ```
//!
//! The synthetic data is `--pairs` (default 2 000 000) NovaSeq-like pairs:
//! Illumina names with a comment, 2×150 bases (0.1% N), four-level binned
//! qualities, from a fixed seed. It is written once to
//! `target/fastq-bench/data/` as text and as single-member gzip (libdeflate
//! level 6, like most `.fastq.gz` files: no recorded member sizes, so
//! readkit-gz decodes it region-parallel).
//!
//! Modes (`tN`: `FastqOptions::threads(N)`; `t0` parses on the caller's
//! thread, the default):
//! - `text-tN`: the parser alone over both files held in memory.
//! - `gz-only`: readkit-gz decoding both files at once (default threads),
//!   nothing parsed: the floor the next mode can reach.
//! - `gz-tN`: decoding and parsing, as a tool reads `.fastq.gz` input.
//!
//! N runs over 0, 1, 2, 4 and all cores. On the synthetic data decoding is
//! the slower stage (random bases), so `gz-tN` mostly measures readkit-gz;
//! real files of short reads (`--input`) show the parser's share there too.
//! Results: a table and JSON at
//! `target/fastq-bench/latest.json`, also `bench/baselines/fastq-NAME.json`
//! with `--save NAME`; `--compare NAME` prints the change per mode.

use readkit_batch::Paired;
use readkit_fastq::{FastqInput, FastqOptions, FastqSource};
use readkit_io::BatchSource;
use serde_json::{Value, json};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

/// Bump when the generator changes, so cached data is regenerated.
const DATA_VERSION: u32 = 1;
const BUFFER: usize = 1 << 20;

struct Args {
    modes: Vec<String>,
    pairs: u64,
    repeat: usize,
    input: Option<(PathBuf, PathBuf)>,
    label: Option<String>,
    save: Option<String>,
    compare: Option<String>,
}

fn args() -> Args {
    let mut a = Args {
        modes: Vec::new(),
        pairs: 2_000_000,
        repeat: 5,
        input: None,
        label: None,
        save: None,
        compare: None,
    };
    let mut it = std::env::args().skip(1).filter(|x| x != "--bench");
    while let Some(x) = it.next() {
        let mut val = || it.next().unwrap_or_else(|| panic!("{x} needs a value"));
        match x.as_str() {
            "--mode" => a.modes.push(val()),
            "--pairs" => a.pairs = val().parse().expect("--pairs N"),
            "--repeat" => a.repeat = val().parse().expect("--repeat N"),
            "--input" => {
                let r1 = val();
                let r2 = val();
                a.input = Some((r1.into(), r2.into()));
            }
            "--label" => a.label = Some(val()),
            "--save" => a.save = Some(val()),
            "--compare" => a.compare = Some(val()),
            other => panic!("unknown argument {other}"),
        }
    }
    a
}

fn workspace() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// splitmix64: small, fast, the same everywhere.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
}

/// Writes `pairs` synthetic pairs to `r1` and `r2`.
fn generate(pairs: u64, r1: &Path, r2: &Path) {
    let mut rng = Rng(0x5eed_f00d);
    let mut out = [
        std::io::BufWriter::with_capacity(BUFFER, File::create(r1).unwrap()),
        std::io::BufWriter::with_capacity(BUFFER, File::create(r2).unwrap()),
    ];
    let mut seq = Vec::with_capacity(150);
    let mut qual = Vec::with_capacity(150);
    for i in 0..pairs {
        let (tile, x, y) = (
            1101 + (i / 400_000) % 80,
            rng.next() % 32_000,
            rng.next() % 40_000,
        );
        for (mate, o) in out.iter_mut().enumerate() {
            seq.clear();
            qual.clear();
            for _ in 0..150 {
                let r = rng.next();
                if r.is_multiple_of(1000) {
                    seq.push(b'N');
                    qual.push(b'#');
                    continue;
                }
                seq.push(b"ACGT"[(r >> 10) as usize & 3]);
                qual.push(match (r >> 20) % 100 {
                    0..=84 => b'F',
                    85..=94 => b':',
                    _ => b',',
                });
            }
            writeln!(
                o,
                "@A00123:8:H7KJ2DSXX:1:{tile}:{x}:{y} {}:N:0:ACGTACGT+TGCATGCA",
                mate + 1
            )
            .unwrap();
            o.write_all(&seq).unwrap();
            o.write_all(b"\n+\n").unwrap();
            o.write_all(&qual).unwrap();
            o.write_all(b"\n").unwrap();
        }
    }
    for o in &mut out {
        o.flush().unwrap();
    }
}

fn gzip(src: &Path, dst: &Path) {
    let text = std::fs::read(src).unwrap();
    let mut c = libdeflater::Compressor::new(libdeflater::CompressionLvl::new(6).unwrap());
    let mut out = vec![0; c.gzip_compress_bound(text.len())];
    let n = c.gzip_compress(&text, &mut out).unwrap();
    std::fs::write(dst, &out[..n]).unwrap();
}

/// The synthetic text and gzip files, generated once per size.
fn dataset(pairs: u64) -> [PathBuf; 4] {
    let dir = workspace().join(format!("target/fastq-bench/data/v{DATA_VERSION}-{pairs}"));
    let files = ["r1.fq", "r2.fq", "r1.fq.gz", "r2.fq.gz"].map(|f| dir.join(f));
    if !files.iter().all(|f| f.exists()) {
        std::fs::create_dir_all(&dir).unwrap();
        println!("generating {pairs} pairs in {}", dir.display());
        generate(pairs, &files[0], &files[1]);
        gzip(&files[0], &files[2]);
        gzip(&files[1], &files[3]);
    }
    files
}

fn is_gz(p: &Path) -> bool {
    p.extension().is_some_and(|e| e == "gz")
}

fn main() {
    let a = args();
    let cores = std::thread::available_parallelism().map_or(4, |n| n.get());
    let (label, text, gz) = match &a.input {
        Some((r1, r2)) => {
            let label = a.label.clone().unwrap_or_else(|| "input".into());
            if is_gz(r1) {
                (label, None, Some([r1.clone(), r2.clone()]))
            } else {
                (label, Some([r1.clone(), r2.clone()]), None)
            }
        }
        None => {
            let [t1, t2, g1, g2] = dataset(a.pairs);
            (
                format!("synthetic-{}", a.pairs),
                Some([t1, t2]),
                Some([g1, g2]),
            )
        }
    };
    let mut threads = vec![0, 1, 2, 4, cores];
    threads.dedup();
    let mut modes = Vec::new();
    if text.is_some() {
        modes.extend(threads.iter().map(|t| format!("text-t{t}")));
    }
    if gz.is_some() {
        modes.push("gz-only".to_string());
        modes.extend(threads.iter().map(|t| format!("gz-t{t}")));
    }
    // Text held in memory for `text-…` modes (read once).
    let held: Option<[&'static [u8]; 2]> = text
        .as_ref()
        .filter(|_| {
            modes
                .iter()
                .any(|m| m.starts_with("text-") && wanted(&a, m))
        })
        .map(|[r1, r2]| [r1, r2].map(|p| &*Vec::leak(std::fs::read(p).unwrap())));
    println!("== {label}");
    let mut results = Vec::new();
    for mode in &modes {
        if !wanted(&a, mode) {
            continue;
        }
        let mut secs = Vec::new();
        let (mut bytes, mut pairs) = (0, 0);
        for _ in 0..a.repeat.max(1) {
            let t = Instant::now();
            (bytes, pairs) = measure(mode, held, gz.as_ref());
            secs.push(t.elapsed().as_secs_f64());
        }
        let med = median(&secs);
        let rate = match pairs {
            0 => "     -".to_string(),
            p => format!("{:6.2}", p as f64 / med / 1e6),
        };
        println!(
            "  {mode:10} {med:7.3}s {:6.0} MB/s {rate} Mpairs/s | spread {:4.1}%",
            bytes as f64 / med / 1e6,
            (max(&secs) - min(&secs)) / med * 100.0
        );
        results.push(json!({
            "data": label,
            "mode": mode,
            "pairs": pairs,
            "bytes": bytes,
            "runs_s": secs,
            "median_s": med,
            "mb_s": bytes as f64 / med / 1e6,
            "mpairs_s": pairs as f64 / med / 1e6,
        }));
    }
    report(&a, results);
}

fn wanted(a: &Args, mode: &str) -> bool {
    a.modes.is_empty() || a.modes.iter().any(|m| m == mode)
}

/// One pass: (FASTQ text bytes, pairs).
fn measure(mode: &str, held: Option<[&'static [u8]; 2]>, gz: Option<&[PathBuf; 2]>) -> (u64, u64) {
    if mode == "gz-only" {
        let gz = gz.unwrap();
        return std::thread::scope(|s| {
            let h: Vec<_> = gz
                .iter()
                .map(|p| {
                    s.spawn(move || {
                        let mut r = readkit_gz::open(p, &Default::default()).unwrap();
                        let mut buf = vec![0; BUFFER];
                        let mut n = 0u64;
                        loop {
                            match r.read(&mut buf).unwrap() {
                                0 => return n,
                                k => n += k as u64,
                            }
                        }
                    })
                })
                .collect();
            (h.into_iter().map(|h| h.join().unwrap()).sum(), 0)
        });
    }
    let (kind, t) = mode.split_once("-t").expect("mode KIND-tN");
    let opts = FastqOptions::default().threads(t.parse().expect("mode KIND-tN"));
    let inputs: [Box<dyn BufRead + Send>; 2] = match kind {
        "text" => held
            .unwrap()
            .map(|b| Box::new(b) as Box<dyn BufRead + Send>),
        "gz" => gz.unwrap().clone().map(|p| {
            let r = readkit_gz::open(&p, &Default::default()).unwrap();
            Box::new(BufReader::with_capacity(BUFFER, r)) as Box<dyn BufRead + Send>
        }),
        _ => panic!("unknown mode {mode}"),
    };
    let [r1, r2] = inputs;
    let mut src = FastqSource::<Paired>::paired(FastqInput::new(r1), FastqInput::new(r2), opts);
    let (mut bytes, mut pairs) = (0u64, 0u64);
    while let Some(b) = src.next_batch().unwrap() {
        pairs += b.len() as u64;
        for s in 0..2 {
            let c = b.segment(s);
            // Text bytes: name, sequence, quality, `@`, `+` and four `\n`s.
            bytes += (c.names().map_or(0, |n| n.data().len()) + 2 * c.seqs().data().len()) as u64
                + 6 * b.len() as u64;
        }
        std::hint::black_box(&b);
    }
    (bytes, pairs)
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
        "repeat": a.repeat,
        "results": results,
    });
    let text = serde_json::to_string_pretty(&doc).unwrap();
    let out = workspace().join("target/fastq-bench");
    std::fs::create_dir_all(&out).unwrap();
    std::fs::write(out.join("latest.json"), &text).unwrap();
    if let Some(name) = &a.save {
        let dir = workspace().join("bench/baselines");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(format!("fastq-{name}.json"));
        std::fs::write(&path, &text).unwrap();
        println!("saved {}", path.display());
    }
    if let Some(name) = &a.compare {
        let path = workspace().join(format!("bench/baselines/fastq-{name}.json"));
        let base: Value =
            serde_json::from_slice(&std::fs::read(&path).expect("baseline")).expect("JSON");
        println!(
            "== against {name} ({})",
            base["git"].as_str().unwrap_or("?")
        );
        for r in doc["results"].as_array().unwrap() {
            let old = base["results"].as_array().into_iter().flatten().find(|b| {
                b["data"] == r["data"] && b["mode"] == r["mode"] && b["pairs"] == r["pairs"]
            });
            let new = r["median_s"].as_f64().unwrap();
            match old.and_then(|o| o["median_s"].as_f64()) {
                Some(old) => println!(
                    "  {:10} {:7.3}s -> {:7.3}s  {:+6.1}%",
                    r["mode"].as_str().unwrap(),
                    old,
                    new,
                    (new / old - 1.0) * 100.0
                ),
                None => println!("  {:10} new", r["mode"].as_str().unwrap()),
            }
        }
    }
}
