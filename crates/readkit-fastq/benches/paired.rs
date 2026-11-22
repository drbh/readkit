//! Paired-read parsing throughput: `FastqSource<Paired>` against fasterp's
//! batch reader (vendored in `tests/support/fasterp.rs`) on the same files.
//!
//! ```text
//! cargo bench -p readkit-fastq --bench paired -- r1.fq r2.fq [runs]
//! ```
//!
//! Without files it reads the `parse` benchmark's synthetic pairs from
//! `target/fastq-bench/data/`, and does nothing if that has not run yet.
//!
//! Each variant reads both files from the start through 1 MiB `BufReader`s
//! with an 8 MiB batch budget, on one thread, and the best of `runs`
//! (default 5) is reported. fasterp validates records lazily on its worker
//! threads, so it is measured twice: splitting only, and splitting plus
//! validating every record (the work `FastqSource` does before it returns a
//! batch). `FastqSource` also copies names, sequences and qualities into
//! separate columns, which fasterp's raw blocks do not.

#[path = "../tests/support/fasterp.rs"]
mod fasterp;

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::time::{Duration, Instant};

use readkit_batch::Paired;
use readkit_fastq::{FastqInput, FastqOptions, FastqSource};
use readkit_io::BatchSource;

const BUFFER: usize = 1 << 20;
const BUDGET: usize = 8 << 20;

type Variant<'a> = (&'static str, Box<dyn FnMut() -> u64 + 'a>);

fn open(path: &str) -> BufReader<File> {
    BufReader::with_capacity(BUFFER, File::open(path).unwrap())
}

fn readkit(r1: &str, r2: &str, keep_names: bool) -> u64 {
    let mut source = FastqSource::<Paired>::paired(
        FastqInput::named(open(r1), r1),
        FastqInput::named(open(r2), r2),
        FastqOptions::default()
            .batch_bytes(BUDGET)
            .keep_names(keep_names),
    );
    let mut pairs = 0;
    while let Some(b) = source.next_batch().unwrap() {
        pairs += b.len() as u64;
    }
    pairs
}

fn fasterp_blocks(r1: &str, r2: &str, validate: bool) -> u64 {
    let mut first = fasterp::Reader::new(Box::new(open(r1)), false);
    let mut second = fasterp::Reader::new(Box::new(open(r2)), false);
    let mut pairs = 0;
    while let Some(blocks) =
        fasterp::read_batch(&mut first, Some(&mut second), BUDGET, usize::MAX).unwrap()
    {
        if validate {
            for mate in &blocks {
                for i in 0..mate.len() {
                    std::hint::black_box(mate.read(i, pairs as usize + i, false).unwrap());
                }
            }
        }
        pairs += blocks[0].len() as u64;
    }
    pairs
}

fn fasterp_owned(r1: &str, r2: &str) -> u64 {
    let mut first = fasterp::Reader::new(Box::new(open(r1)), false);
    let mut second = fasterp::Reader::new(Box::new(open(r2)), false);
    let mut pairs = 0;
    while fasterp::read_pair(&mut first, Some(&mut second), false)
        .unwrap()
        .is_some()
    {
        pairs += 1;
    }
    pairs
}

fn best(runs: usize, mut f: impl FnMut() -> u64) -> (Duration, u64) {
    let mut best = Duration::MAX;
    let mut pairs = 0;
    for _ in 0..runs {
        let t = Instant::now();
        pairs = f();
        best = best.min(t.elapsed());
    }
    (best, pairs)
}

/// The text files the `parse` benchmark generated, if any.
fn synthetic() -> Option<[String; 2]> {
    let data = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../target/fastq-bench/data");
    let dir = std::fs::read_dir(data)
        .ok()?
        .flatten()
        .map(|e| e.path())
        .max()?;
    let file = |name: &str| {
        let path = dir.join(name);
        path.exists().then(|| path.to_string_lossy().into_owned())
    };
    Some([file("r1.fq")?, file("r2.fq")?])
}

fn main() {
    let mut args: Vec<String> = std::env::args()
        .skip(1)
        .filter(|a| a != "--bench")
        .collect();
    if args.len() < 2 {
        let Some(found) = synthetic() else {
            eprintln!("no input: pass r1.fq r2.fq, or run `cargo bench --bench parse` first");
            return;
        };
        args.splice(0..0, found);
    }
    let (r1, r2) = (&args[0], &args[1]);
    let runs = args.get(2).map_or(5, |r| r.parse().unwrap());
    let bytes = std::fs::metadata(r1).unwrap().len() + std::fs::metadata(r2).unwrap().len();
    // Warm the page cache.
    fasterp_blocks(r1, r2, false);
    let variants: [Variant<'_>; 5] = [
        (
            "readkit FastqSource (names)",
            Box::new(|| readkit(r1, r2, true)),
        ),
        (
            "readkit FastqSource (no names)",
            Box::new(|| readkit(r1, r2, false)),
        ),
        (
            "fasterp blocks, split only",
            Box::new(|| fasterp_blocks(r1, r2, false)),
        ),
        (
            "fasterp blocks + validation",
            Box::new(|| fasterp_blocks(r1, r2, true)),
        ),
        ("fasterp owned Reader", Box::new(|| fasterp_owned(r1, r2))),
    ];
    println!("{bytes} bytes, best of {runs}");
    for (name, mut f) in variants {
        let (t, pairs) = best(runs, &mut f);
        println!(
            "{name:32} {pairs:>9} pairs {:>8.1} ms {:>8.0} MB/s",
            t.as_secs_f64() * 1e3,
            bytes as f64 / t.as_secs_f64() / 1e6
        );
    }
}
