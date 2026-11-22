//! Batch building and reading throughput, on synthetic 2×150 pairs.
//!
//! ```text
//! cargo bench -p readkit-batch --bench build
//! cargo bench -p readkit-batch --bench build -- --records 8000000 --repeat 10
//! ```
//!
//! Cases:
//! - `push_pair`: copying named reads into batches from one recycling pool.
//! - `push_pair, no reuse`: the same with a new pool per batch, so every
//!   buffer is allocated.
//! - `push_with`: filling sequence and quality slices in place, without names.
//! - `pairs()`: counting G and C through per-record views.
//! - `columns`: the same count over each segment's sequence buffer.
//!
//! The best of `--repeat` (default 5) runs over `--records` (default
//! 2 000 000) pairs is reported.

use std::hint::black_box;
use std::time::{Duration, Instant};

use readkit_batch::{Batch, BatchBuilder, Paired, Pool, Read};

const LEN: usize = 150;
const BATCH: usize = 16_384;
const NAME: &[u8] = b"A00123:45:HXXXXXXXX:1:1101:1000:2000";

fn best(repeat: usize, mut f: impl FnMut() -> u64) -> (Duration, u64) {
    let mut best = Duration::MAX;
    let mut bytes = 0;
    for _ in 0..repeat {
        let t = Instant::now();
        bytes = black_box(f());
        best = best.min(t.elapsed());
    }
    (best, bytes)
}

fn push_pair(pool: &Pool, seq: &[u8], qual: &[u8]) -> Batch<Paired> {
    let mut b = BatchBuilder::<Paired>::with_pool(pool, BATCH);
    for _ in 0..BATCH {
        b.push_pair(Read::named(NAME, seq, qual), Read::named(NAME, seq, qual))
            .expect("valid pair");
    }
    b.finish()
}

fn push_with(pool: &Pool, seq: &[u8], qual: &[u8]) -> Batch<Paired> {
    let mut b = BatchBuilder::<Paired>::with_pool(pool, BATCH);
    for _ in 0..BATCH {
        b.push_with(&[LEN, LEN], None, |_, s, q| {
            s.copy_from_slice(seq);
            q.copy_from_slice(qual);
        })
        .expect("valid pair");
    }
    b.finish()
}

fn is_gc(base: &u8) -> bool {
    matches!(base, b'G' | b'C')
}

fn main() {
    let mut records = 2_000_000usize;
    let mut repeat = 5usize;
    let mut it = std::env::args().skip(1).filter(|x| x != "--bench");
    while let Some(x) = it.next() {
        let mut val = || it.next().unwrap_or_else(|| panic!("{x} needs a value"));
        match x.as_str() {
            "--records" => records = val().parse().expect("--records N"),
            "--repeat" => repeat = val().parse().expect("--repeat N"),
            other => panic!("unknown argument {other}"),
        }
    }
    let batches = records.div_ceil(BATCH);
    let seq: Vec<u8> = (0..LEN).map(|i| b"ACGT"[i * 7 % 4]).collect();
    let qual = vec![b'F'; LEN];
    let pool = Pool::new();
    let built = push_pair(&pool, &seq, &qual);

    type Case<'a> = (&'static str, Box<dyn FnMut() -> u64 + 'a>);
    let cases: [Case<'_>; 5] = [
        (
            "push_pair",
            Box::new(|| {
                (0..batches)
                    .map(|_| push_pair(&pool, &seq, &qual).byte_len() as u64)
                    .sum()
            }),
        ),
        (
            "push_pair, no reuse",
            Box::new(|| {
                (0..batches)
                    .map(|_| push_pair(&Pool::new(), &seq, &qual).byte_len() as u64)
                    .sum()
            }),
        ),
        (
            "push_with",
            Box::new(|| {
                (0..batches)
                    .map(|_| push_with(&pool, &seq, &qual).byte_len() as u64)
                    .sum()
            }),
        ),
        (
            "pairs()",
            Box::new(|| {
                for _ in 0..batches {
                    let gc: usize = built
                        .pairs()
                        .map(|(a, b)| a.seq().iter().chain(b.seq()).filter(|x| is_gc(x)).count())
                        .sum();
                    black_box(gc);
                }
                (batches * built.byte_len()) as u64
            }),
        ),
        (
            "columns",
            Box::new(|| {
                for _ in 0..batches {
                    let gc: usize = built
                        .columns()
                        .map(|c| c.seqs().data().iter().filter(|x| is_gc(x)).count())
                        .sum();
                    black_box(gc);
                }
                (batches * built.byte_len()) as u64
            }),
        ),
    ];
    println!(
        "{} pairs in batches of {BATCH}, best of {repeat}",
        batches * BATCH
    );
    for (name, mut f) in cases {
        let (t, bytes) = best(repeat, &mut f);
        let s = t.as_secs_f64();
        println!(
            "{name:22} {:>8.1} ms {:>7.1} M pairs/s {:>7.0} MB/s",
            s * 1e3,
            (batches * BATCH) as f64 / s / 1e6,
            bytes as f64 / s / 1e6
        );
    }
}
