//! A minimal decoder → tool → consumer chain.
//!
//! - The decoder fills paired batches straight from base and quality arrays
//!   (as an SRA decoder would), from its own pool.
//! - The tool tags each pair with its GC percentage bucket, reading ahead of
//!   the decoder on its own thread.
//! - The consumer counts k-mers per tag from the sequence columns, without
//!   per-record views.
//!
//! Run with `cargo run -p readkit-io --example chain`.

use std::convert::Infallible;

use readkit_batch::{Batch, BatchBuilder, Paired, Pool};
use readkit_io::{BatchSink, BatchSource, Error, ReadAhead, pump};

/// Decoder: deterministic pseudo-random spots, 2 x 100 bp.
struct Decoder {
    pool: Pool,
    batches_left: usize,
    rng: u64,
}

impl BatchSource<Paired> for Decoder {
    type Error = Infallible;

    fn next_batch(&mut self) -> Result<Option<Batch<Paired>>, Infallible> {
        if self.batches_left == 0 {
            return Ok(None);
        }
        self.batches_left -= 1;
        let mut b = BatchBuilder::<Paired>::with_pool(&self.pool, 4096);
        for _ in 0..4096 {
            // bases as 2-bit codes and phred values, as VDB stores them
            let mut codes = [0u8; 200];
            for c in &mut codes {
                self.rng ^= self.rng << 13;
                self.rng ^= self.rng >> 7;
                self.rng ^= self.rng << 17;
                *c = (self.rng & 3) as u8;
            }
            b.push_with(&[100, 100], None, |s, seq, qual| {
                for (i, (base, q)) in seq.iter_mut().zip(qual.iter_mut()).enumerate() {
                    *base = b"ACGT"[codes[s * 100 + i] as usize];
                    *q = 30 + 33;
                }
            })
            .expect("valid spot");
        }
        Ok(Some(b.finish()))
    }
}

/// Tool: tags each pair with GC% / 10 over both mates.
struct GcTagger<S: BatchSource<Paired>> {
    input: ReadAhead<Paired, S>,
    pool: Pool,
}

impl<S: BatchSource<Paired> + Send + 'static> GcTagger<S> {
    fn new(input: S) -> Self {
        Self {
            input: ReadAhead::new(input),
            pool: Pool::new(),
        }
    }
}

impl<S: BatchSource<Paired>> BatchSource<Paired> for GcTagger<S> {
    type Error = Error;

    fn next_batch(&mut self) -> Result<Option<Batch<Paired>>, Error> {
        let Some(mut batch) = self.input.next_batch().map_err(Error::upstream)? else {
            return Ok(None);
        };
        let gc: Vec<u16> = batch
            .pairs()
            .map(|(a, b)| {
                let n = a.len() + b.len();
                let gc = a
                    .seq()
                    .iter()
                    .chain(b.seq())
                    .filter(|&&x| x == b'G' || x == b'C')
                    .count();
                (gc * 10 / n.max(1)) as u16
            })
            .collect();
        batch.tags_or_insert(&self.pool, 0).copy_from_slice(&gc);
        Ok(Some(batch))
    }
}

/// Consumer: k-mers per GC bucket, read straight from the columns.
#[derive(Default)]
struct KmerCounter {
    per_bucket: [u64; 11],
    records: u64,
}

impl BatchSink<Paired> for KmerCounter {
    type Error = Infallible;

    fn accept(&mut self, batch: Batch<Paired>) -> Result<(), Infallible> {
        const K: usize = 31;
        let tags = batch.tags().unwrap_or(&[]);
        for col in batch.columns() {
            for (w, &tag) in col.seqs().offsets().windows(2).zip(tags) {
                let len = (w[1] - w[0]) as usize;
                self.per_bucket[tag as usize] += len.saturating_sub(K - 1) as u64;
            }
        }
        self.records += batch.len() as u64;
        Ok(())
    }

    fn finish(self) -> Result<(), Infallible> {
        println!("records: {}", self.records);
        for (bucket, n) in self.per_bucket.iter().enumerate().filter(|(_, n)| **n > 0) {
            println!("GC {:>3}%: {n} 31-mers", bucket * 10);
        }
        Ok(())
    }
}

fn main() -> Result<(), Error> {
    let pool = Pool::new();
    let decoder = Decoder {
        pool: pool.clone(),
        batches_left: 50,
        rng: 0x9E37_79B9_7F4A_7C15,
    };
    let tool = GcTagger::new(decoder);
    let records = pump(tool, KmerCounter::default())?;
    let s = pool.stats();
    println!(
        "{records} records; decoder pool: {} buffers allocated, {} reused",
        s.fresh, s.reused
    );
    Ok(())
}
