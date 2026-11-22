//! Recompresses a gzip file: the parallel decoder feeding the parallel
//! encoder.
//!
//! ```text
//! cargo run --release -p readkit-gz --example recompress -- IN.gz OUT.gz [LEVEL]
//! ```
//!
//! The output is one gzip member per 256 KiB of input, each recording its
//! compressed size, so readkit-gz decodes it member-parallel. LEVEL defaults
//! to 6.

use std::fs::File;
use std::io::BufWriter;

use readkit_gz::{DecodeOptions, EncodeOptions, GzWriter};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [input, output, rest @ ..] = args.as_slice() else {
        return Err("usage: recompress IN.gz OUT.gz [LEVEL]".into());
    };
    let level = match rest.first() {
        Some(level) => level.parse()?,
        None => 6,
    };

    let mut reader = readkit_gz::open(input, &DecodeOptions::default())?;
    let strategy = reader.strategy();
    let options = EncodeOptions {
        level,
        ..Default::default()
    };
    let mut writer = GzWriter::new(BufWriter::new(File::create(output)?), &options)?;
    let bytes = std::io::copy(&mut reader, &mut writer)?;
    writer.finish()?;
    eprintln!("{bytes} bytes, decoded with {strategy:?}, written at level {level}");
    Ok(())
}
