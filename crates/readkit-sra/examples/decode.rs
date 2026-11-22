//! Decodes a local `.sra` file into batches and prints what it held.
//!
//! ```text
//! cargo run --release -p readkit-sra --example decode -- SRR000001.sra
//! cargo run --release -p readkit-sra --features ngs --example decode -- SRR000001.sra
//! ```
//!
//! Without `ngs` the decoder runs `fasterq-dump`, which must be on `PATH`.
//! SRA Lite runs are rejected by default (`LitePolicy`).

use readkit_io::BatchSource;
use readkit_sra::{DecodeOptions, SegmentsDecoder};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args().nth(1).ok_or("usage: decode FILE.sra")?;
    let mut decoder = SegmentsDecoder::open(&path, DecodeOptions::default())?;
    eprintln!(
        "{} segments per spot, {:?} qualities, backend {}",
        decoder.segments(),
        decoder.quality(),
        decoder.backend_name(),
    );

    let mut gc = 0u64;
    while let Some(batch) = decoder.next_batch()? {
        // Bases arrive as one buffer per segment; no FASTQ text is formatted.
        for col in batch.columns() {
            gc += col
                .seqs()
                .data()
                .iter()
                .filter(|x| matches!(x, b'G' | b'C'))
                .count() as u64;
        }
    }
    let s = decoder.stats();
    println!(
        "{} spots, {} reads, {} bases ({:.1}% GC), longest read {}",
        s.spots,
        s.reads,
        s.bases,
        gc as f64 * 100.0 / s.bases.max(1) as f64,
        s.max_read_len,
    );
    Ok(())
}
