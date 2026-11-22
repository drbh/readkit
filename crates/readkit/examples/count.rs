//! Counts the reads and bases of a FASTQ file, gzip or plain.
//!
//! ```text
//! cargo run --release -p readkit --example count -- reads.fastq.gz
//! ```
//!
//! gzip input is decoded in parallel where the file allows it, and the
//! parser reads ahead of the counting loop on its own thread.

use std::fs::File;
use std::io::BufReader;

use readkit::fastq::{FastqInput, FastqOptions, FastqSource};
use readkit::gz::DecodeOptions;
use readkit::{BatchSource, ReadAhead};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args()
        .nth(1)
        .ok_or("usage: count FILE.fastq[.gz]")?;
    let input = if path.ends_with(".gz") {
        FastqInput::named(
            readkit::gz::open(&path, &DecodeOptions::default())?,
            path.as_str(),
        )
    } else {
        FastqInput::named(
            BufReader::with_capacity(1 << 20, File::open(&path)?),
            path.as_str(),
        )
    };
    // Names are not needed to count, so their columns are not filled.
    let options = FastqOptions::default().keep_names(false);
    let mut source = ReadAhead::new(FastqSource::single(input, options));

    let (mut reads, mut bases) = (0u64, 0u64);
    while let Some(batch) = source.next_batch()? {
        reads += batch.len() as u64;
        bases += batch.column().seqs().data().len() as u64;
    }
    println!("{reads} reads, {bases} bases");
    Ok(())
}
