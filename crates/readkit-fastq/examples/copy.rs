//! Copies FASTQ through batches: `k` inputs in lockstep to `k` outputs.
//!
//! ```text
//! cargo run --release -p readkit-fastq --example copy -- [--plus-name] [--phred64] in1 [in2 ...] -- out1 [out2 ...]
//! ```
//!
//! With all filtering disabled, `fasterp -A -G -Q -L` writes the same bytes
//! for inputs without empty reads (`--plus-name` when the input's `+` lines
//! repeat the name).

use std::fs::File;
use std::io::{BufReader, BufWriter};

use readkit_batch::Segments;
use readkit_fastq::{FastqInput, FastqOptions, FastqSink, FastqSource, PlusLine};
use readkit_io::{ReadAhead, pump};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let mut flag = |name: &str| {
        let found = args.iter().any(|a| a == name);
        args.retain(|a| a != name);
        found
    };
    let plus = if flag("--plus-name") {
        PlusLine::Name
    } else {
        PlusLine::Bare
    };
    let phred64 = flag("--phred64");
    let split = args.iter().position(|a| a == "--").ok_or("missing --")?;
    let (ins, outs) = (&args[..split], &args[split + 1..]);
    let inputs = ins
        .iter()
        .map(|p| {
            Ok(FastqInput::named(
                BufReader::with_capacity(1 << 20, File::open(p)?),
                p.as_str(),
            ))
        })
        .collect::<std::io::Result<Vec<_>>>()?;
    let writers = outs
        .iter()
        .map(|p| Ok(BufWriter::with_capacity(1 << 20, File::create(p)?)))
        .collect::<std::io::Result<Vec<_>>>()?;
    let source =
        FastqSource::<Segments>::lockstep(inputs, FastqOptions::default().phred64(phred64))?;
    let spots = pump(
        ReadAhead::new(source),
        FastqSink::split(writers).with_plus_line(plus),
    )?;
    eprintln!("{spots} spots");
    Ok(())
}
