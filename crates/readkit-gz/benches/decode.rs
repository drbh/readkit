//! Decode throughput, in MB/s of decompressed output.
//!
//! ```text
//! cargo bench -p readkit-gz --bench decode                      # every method
//! cargo bench -p readkit-gz --bench decode -- --method regions --threads 8
//! cargo bench -p readkit-gz --bench decode -- --input reads.fastq.gz
//! ```
//!
//! Without `--input`, `--reads` (default 1 000 000) reads of synthetic FASTQ
//! are compressed once to `target/gz-bench/` as a single libdeflate level 6
//! member, like most `.fastq.gz` files.
//!
//! Methods (`--method`, repeatable; all by default):
//! - `flate2`: flate2's `MultiGzDecoder`, one thread;
//! - `sequential`: [`readkit_gz::GzReader`] on the calling thread;
//! - `auto`: [`readkit_gz::open`] with default options and `--threads`;
//! - `regions`: region decoding forced ([`readkit_gz::Split::Always`]);
//! - `never`: [`readkit_gz::Split::Never`] (sequential on a background
//!   thread unless members record their sizes).
//!
//! Prints the best of `--repeat` (default 3) and the CRC32 of the output, so
//! methods can be checked against each other. `--threads` defaults to 4, and
//! `--in-flight` sets [`readkit_gz::DecodeOptions::in_flight`].

#[path = "../tests/common/mod.rs"]
mod common;

use std::io::{BufReader, Read};
use std::path::PathBuf;
use std::time::Instant;

const METHODS: [&str; 5] = ["flate2", "sequential", "auto", "regions", "never"];

fn drain(mut reader: impl Read) -> (u64, u32) {
    let mut buf = vec![0u8; 1 << 20];
    let mut crc = crc32fast::Hasher::new();
    let mut total = 0u64;
    loop {
        let n = reader.read(&mut buf).expect("decoding failed");
        if n == 0 {
            return (total, crc.finalize());
        }
        crc.update(&buf[..n]);
        total += n as u64;
    }
}

/// The synthetic input, written on first use.
fn synthetic(reads: usize) -> PathBuf {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../target/gz-bench");
    let path = dir.join(format!("synthetic-{reads}.fastq.gz"));
    if !path.exists() {
        eprintln!("writing {}", path.display());
        std::fs::create_dir_all(&dir).unwrap();
        let gz = common::libdeflate_gz(&common::fastq(reads, 0x5eed_f00d), 6);
        std::fs::write(&path, gz).unwrap();
    }
    path
}

fn main() {
    let mut input = None;
    let mut methods = Vec::new();
    let mut reads = 1_000_000usize;
    let mut threads = 4usize;
    let mut repeat = 3usize;
    let mut in_flight = None;
    let mut it = std::env::args().skip(1).filter(|x| x != "--bench");
    while let Some(x) = it.next() {
        let mut val = || it.next().unwrap_or_else(|| panic!("{x} needs a value"));
        match x.as_str() {
            "--input" => input = Some(PathBuf::from(val())),
            "--method" => methods.push(val()),
            "--reads" => reads = val().parse().expect("--reads N"),
            "--threads" => threads = val().parse().expect("--threads N"),
            "--repeat" => repeat = val().parse().expect("--repeat N"),
            "--in-flight" => in_flight = Some(val().parse().expect("--in-flight N")),
            other => panic!("unknown argument {other}"),
        }
    }
    if methods.is_empty() {
        methods = METHODS.map(String::from).to_vec();
    }
    let path = input.unwrap_or_else(|| synthetic(reads));
    let options = |split| readkit_gz::DecodeOptions {
        threads: Some(threads),
        split,
        in_flight,
        ..Default::default()
    };
    println!("{}, best of {repeat}", path.display());
    for method in &methods {
        let mut best = f64::INFINITY;
        let mut result = (0, 0);
        for _ in 0..repeat {
            let file = std::fs::File::open(&path).unwrap();
            let t = Instant::now();
            result = match method.as_str() {
                "flate2" => drain(flate2::read::MultiGzDecoder::new(BufReader::with_capacity(
                    1 << 20,
                    file,
                ))),
                "sequential" => drain(readkit_gz::GzReader::from_reader(file)),
                "auto" | "regions" | "never" => {
                    let split = match method.as_str() {
                        "auto" => readkit_gz::Split::Auto,
                        "regions" => readkit_gz::Split::Always,
                        _ => readkit_gz::Split::Never,
                    };
                    drain(readkit_gz::open(&path, &options(split)).unwrap())
                }
                other => panic!("unknown method {other}"),
            };
            best = best.min(t.elapsed().as_secs_f64());
        }
        let (bytes, crc) = result;
        println!(
            "{method:<10} threads {threads:>2}: {:>7.0} MB/s  ({} MB in {best:.3} s, crc {crc:08x})",
            bytes as f64 / 1e6 / best,
            bytes / 1_000_000,
        );
    }
}
