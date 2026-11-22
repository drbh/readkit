//! Builds a paired batch and reads it back by record and by column.
//!
//! Run with `cargo run -p readkit-batch --example columns`.

use readkit_batch::{BatchBuilder, Error, Paired, Pool, Read};

fn main() -> Result<(), Error> {
    let pool = Pool::new();
    for round in 0..2 {
        let mut b = BatchBuilder::<Paired>::with_pool(&pool, 2);
        b.push_pair(
            Read::named(b"a/1", b"ACGT", b"IIII"),
            Read::named(b"a/2", b"TTGGC", b"IIIII"),
        )?;
        b.set_tag(7);
        b.push_pair(
            Read::named(b"b/1", b"GG", b"##"),
            Read::named(b"b/2", b"", b""),
        )?;
        let batch = b.finish();

        if round == 0 {
            // By record.
            for (r1, r2) in batch.pairs() {
                println!(
                    "{} {} | {} {}",
                    text(r1.name().unwrap_or_default()),
                    text(r1.seq()),
                    text(r2.name().unwrap_or_default()),
                    text(r2.seq()),
                );
            }
            // By column: one byte buffer and n + 1 offsets per segment.
            let seqs = batch.r1().seqs();
            println!(
                "r1 offsets {:?} over {:?}",
                seqs.offsets(),
                text(seqs.data())
            );
            println!("tags {:?}", batch.tags());
        }
        // Dropping the batch returns its buffers to the pool, so the second
        // round allocates nothing.
    }
    let s = pool.stats();
    println!("pool: {} buffers allocated, {} reused", s.fresh, s.reused);
    Ok(())
}

fn text(bytes: &[u8]) -> &str {
    std::str::from_utf8(bytes).unwrap_or("<binary>")
}
