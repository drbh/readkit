//! End-to-end throughput through the facade: gzip → FASTQ parser → batches,
//! and on to the FASTQ writer and gzip encoder, in MB/s of FASTQ text.
//!
//! ```text
//! cargo bench -p readkit --bench pipeline
//! cargo bench -p readkit --bench pipeline -- --reads 4000000 --repeat 5
//! ```
//!
//! The input is `--reads` (default 1 000 000) synthetic 150-base reads with
//! binned qualities, compressed in memory by `GzWriter`, so it decodes
//! member-parallel.
//!
//! Cases:
//! - `gz > fastq`: decoding and parsing, counting bases.
//! - `gz > fastq > fastq`: also formatting FASTQ text again.
//! - `gz > fastq > fastq > gz`: also compressing it, as a tool that rewrites
//!   its input does.
//!
//! The best of `--repeat` (default 3) runs is reported.

use std::io::Write;
use std::sync::Arc;
use std::time::{Duration, Instant};

use readkit::fastq::{FastqInput, FastqOptions, FastqSink, FastqSource};
use readkit::gz::{DecodeOptions, EncodeOptions, GzWriter, ParGzReader};
use readkit::{BatchSource, ReadAhead, Single, pump};

const LEN: usize = 150;

/// `reads` records of FASTQ text from a fixed seed.
fn fastq(reads: usize) -> Vec<u8> {
    let mut state = 0x5eed_f00d_u64;
    let mut next = |m: u64| {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) % m
    };
    let mut out = Vec::with_capacity(reads * (2 * LEN + 48));
    for i in 0..reads {
        writeln!(out, "@A00123:45:HXXXXXXXX:1:1101:{i}:2000 1:N:0:ACGT").unwrap();
        out.extend((0..LEN).map(|_| b"ACGT"[next(4) as usize]));
        out.extend_from_slice(b"\n+\n");
        out.extend((0..LEN).map(|_| b"#,:F"[next(4) as usize]));
        out.push(b'\n');
    }
    out
}

fn source(gz: &Arc<Vec<u8>>) -> ReadAhead<Single, FastqSource<Single>> {
    let reader = ParGzReader::new(gz.clone(), gz.len() as u64, &DecodeOptions::default())
        .expect("valid gzip");
    ReadAhead::new(FastqSource::single(
        FastqInput::new(reader),
        FastqOptions::default(),
    ))
}

/// Counts the bytes written to it.
struct Count(u64);

impl Write for Count {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0 += buf.len() as u64;
        Ok(buf.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

fn main() {
    let mut reads = 1_000_000usize;
    let mut repeat = 3usize;
    let mut it = std::env::args().skip(1).filter(|x| x != "--bench");
    while let Some(x) = it.next() {
        let mut val = || it.next().unwrap_or_else(|| panic!("{x} needs a value"));
        match x.as_str() {
            "--reads" => reads = val().parse().expect("--reads N"),
            "--repeat" => repeat = val().parse().expect("--repeat N"),
            other => panic!("unknown argument {other}"),
        }
    }
    let text = fastq(reads);
    let mut writer = GzWriter::new(Vec::new(), &EncodeOptions::default()).unwrap();
    writer.write_all(&text).unwrap();
    let gz = Arc::new(writer.finish().unwrap());
    let bytes = text.len();
    drop(text);

    type Case<'a> = (&'static str, Box<dyn Fn() -> u64 + 'a>);
    let cases: [Case<'_>; 3] = [
        (
            "gz > fastq",
            Box::new(|| {
                let mut source = source(&gz);
                let mut bases = 0;
                while let Some(batch) = source.next_batch().unwrap() {
                    bases += batch.column().seqs().data().len() as u64;
                }
                assert_eq!(bases, (reads * LEN) as u64);
                bytes as u64
            }),
        ),
        (
            "gz > fastq > fastq",
            Box::new(|| {
                let mut sink = FastqSink::new(Count(0));
                pump(source(&gz), &mut sink).unwrap();
                sink.into_inner().unwrap()[0].0
            }),
        ),
        (
            "gz > fastq > fastq > gz",
            Box::new(|| {
                let writer = GzWriter::new(Count(0), &EncodeOptions::default()).unwrap();
                let mut sink = FastqSink::new(writer);
                pump(source(&gz), &mut sink).unwrap();
                let writer = sink.into_inner().unwrap().pop().unwrap();
                writer.finish().unwrap();
                bytes as u64
            }),
        ),
    ];
    println!(
        "{reads} reads, {} MB of FASTQ in {} MB of gzip, best of {repeat}",
        bytes / 1_000_000,
        gz.len() / 1_000_000
    );
    for (name, run) in cases {
        let mut best = Duration::MAX;
        for _ in 0..repeat {
            let t = Instant::now();
            assert_eq!(run(), bytes as u64, "{name}");
            best = best.min(t.elapsed());
        }
        let s = best.as_secs_f64();
        println!(
            "{name:24} {:>8.1} ms {:>7.0} MB/s",
            s * 1e3,
            bytes as f64 / s / 1e6
        );
    }
}
