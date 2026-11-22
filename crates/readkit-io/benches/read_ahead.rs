//! What `ReadAhead` buys a two-stage chain: a producer that fills batches
//! and a consumer that scans them, each doing real work.
//!
//! ```text
//! cargo bench -p readkit-io --bench read_ahead
//! cargo bench -p readkit-io --bench read_ahead -- --batches 400 --repeat 10
//! ```
//!
//! Cases:
//! - `direct`: the consumer pulls from the producer on one thread.
//! - `sequential`: through `ReadAhead::sequential`, also one thread; the
//!   wrapper's own cost.
//! - `read-ahead`: through `ReadAhead::new`, the producer one batch ahead on
//!   its own thread.
//!
//! The best of `--repeat` (default 5) runs over `--batches` (default 200)
//! batches of 4096 pairs of 2×100 bases is reported.

use std::convert::Infallible;
use std::time::{Duration, Instant};

use readkit_batch::{Batch, BatchBuilder, Paired, Pool};
use readkit_io::{BatchSink, BatchSource, ReadAhead, pump};

const RECORDS: usize = 4096;
const LEN: usize = 100;

/// Fills pairs with pseudo-random bases, as a decoder would.
struct Producer {
    pool: Pool,
    batches_left: usize,
    rng: u64,
}

impl BatchSource<Paired> for Producer {
    type Error = Infallible;

    fn next_batch(&mut self) -> Result<Option<Batch<Paired>>, Infallible> {
        if self.batches_left == 0 {
            return Ok(None);
        }
        self.batches_left -= 1;
        let mut b = BatchBuilder::<Paired>::with_pool(&self.pool, RECORDS);
        for _ in 0..RECORDS {
            let rng = &mut self.rng;
            b.push_with(&[LEN, LEN], None, |_, seq, qual| {
                for base in seq.iter_mut() {
                    *rng ^= *rng << 13;
                    *rng ^= *rng >> 7;
                    *rng ^= *rng << 17;
                    *base = b"ACGT"[(*rng & 3) as usize];
                }
                qual.fill(b'F');
            })
            .expect("valid pair");
        }
        Ok(Some(b.finish()))
    }
}

/// Counts G and C over the sequence columns.
#[derive(Default)]
struct Consumer {
    gc: u64,
}

impl BatchSink<Paired> for &mut Consumer {
    type Error = Infallible;

    fn accept(&mut self, batch: Batch<Paired>) -> Result<(), Infallible> {
        for col in batch.columns() {
            self.gc += col
                .seqs()
                .data()
                .iter()
                .filter(|x| matches!(x, b'G' | b'C'))
                .count() as u64;
        }
        Ok(())
    }

    fn finish(self) -> Result<(), Infallible> {
        Ok(())
    }
}

fn main() {
    let mut batches = 200usize;
    let mut repeat = 5usize;
    let mut it = std::env::args().skip(1).filter(|x| x != "--bench");
    while let Some(x) = it.next() {
        let mut val = || it.next().unwrap_or_else(|| panic!("{x} needs a value"));
        match x.as_str() {
            "--batches" => batches = val().parse().expect("--batches N"),
            "--repeat" => repeat = val().parse().expect("--repeat N"),
            other => panic!("unknown argument {other}"),
        }
    }
    let pool = Pool::new();
    let producer = || Producer {
        pool: pool.clone(),
        batches_left: batches,
        rng: 0x9E37_79B9_7F4A_7C15,
    };
    type Case<'a> = (&'static str, Box<dyn Fn(&mut Consumer) -> u64 + 'a>);
    let cases: [Case<'_>; 3] = [
        ("direct", Box::new(|c| pump(producer(), c).unwrap())),
        (
            "sequential",
            Box::new(|c| pump(ReadAhead::sequential(producer()), c).unwrap()),
        ),
        (
            "read-ahead",
            Box::new(|c| pump(ReadAhead::new(producer()), c).unwrap()),
        ),
    ];
    println!("{batches} batches of {RECORDS} pairs, best of {repeat}");
    let mut expected = None;
    for (name, run) in cases {
        let mut best = Duration::MAX;
        let mut records = 0;
        for _ in 0..repeat {
            let mut consumer = Consumer::default();
            let t = Instant::now();
            records = run(&mut consumer);
            best = best.min(t.elapsed());
            assert_eq!(*expected.get_or_insert(consumer.gc), consumer.gc, "{name}");
        }
        let s = best.as_secs_f64();
        println!(
            "{name:12} {:>8.1} ms {:>7.1} M pairs/s {:>7.0} MB/s",
            s * 1e3,
            records as f64 / s / 1e6,
            (records as usize * 4 * LEN) as f64 / s / 1e6
        );
    }
}
