//! Shared test helpers: inputs over byte slices, collection of batches, and
//! fasterp's reference reader as the oracle.
#![allow(dead_code)]

pub mod fasterp;

use std::io::{BufReader, Cursor};

use readkit_batch::{Batch, Segments};
use readkit_fastq::{
    Error, ErrorKind, FastqInput, FastqLayout, FastqOptions, FastqSink, FastqSource, PlusLine,
};
use readkit_io::{BatchSink, BatchSource};

/// An input over a copy of `data`, buffered with `capacity` bytes.
pub fn input(data: &[u8], capacity: usize) -> FastqInput {
    FastqInput::new(BufReader::with_capacity(
        capacity,
        Cursor::new(data.to_vec()),
    ))
}

/// A named input over a copy of `data`.
pub fn named(data: &[u8], capacity: usize, name: &str) -> FastqInput {
    FastqInput::named(
        BufReader::with_capacity(capacity, Cursor::new(data.to_vec())),
        name,
    )
}

/// A fasterp reference reader over a copy of `data`.
pub fn reference(data: &[u8], capacity: usize, phred64: bool) -> fasterp::Reader {
    fasterp::Reader::new(
        Box::new(BufReader::with_capacity(
            capacity,
            Cursor::new(data.to_vec()),
        )),
        phred64,
    )
}

/// All batches of `source` as runtime-`k` batches, or the first error.
pub fn drain<L: FastqLayout>(mut source: FastqSource<L>) -> Result<Vec<Batch<Segments>>, Error> {
    let mut out = Vec::new();
    while let Some(b) = source.next_batch()? {
        assert!(!b.is_empty(), "empty batch");
        out.push(b.into_segments());
    }
    // After the end, still the end.
    assert!(source.next_batch()?.is_none());
    Ok(out)
}

/// Like `drain`, but returns the batches read before an error too.
pub fn drain_partial<L: FastqLayout>(
    mut source: FastqSource<L>,
) -> (Vec<Batch<Segments>>, Option<Error>) {
    let mut out = Vec::new();
    loop {
        match source.next_batch() {
            Ok(Some(b)) => out.push(b.into_segments()),
            Ok(None) => return (out, None),
            Err(e) => {
                assert!(
                    source.next_batch().unwrap().is_none(),
                    "source not done after error"
                );
                return (out, Some(e));
            }
        }
    }
}

/// Names of segment `s` over all batches.
pub fn names(batches: &[Batch<Segments>], s: usize) -> Vec<Vec<u8>> {
    batches
        .iter()
        .flat_map(|b| {
            b.segment(s)
                .names()
                .unwrap()
                .iter()
                .map(<[u8]>::to_vec)
                .collect::<Vec<_>>()
        })
        .collect()
}

/// FASTQ text of `batches`, one output per segment (`split`) or one
/// interleaved output.
pub fn write(batches: Vec<Batch<Segments>>, split: bool, plus: PlusLine) -> Vec<Vec<u8>> {
    let k = batches.first().map_or(1, |b| b.segments());
    let writers = if split {
        vec![Vec::new(); k]
    } else {
        vec![Vec::new()]
    };
    let mut sink = FastqSink::split(writers).with_plus_line(plus);
    for b in batches {
        sink.accept(b).unwrap();
    }
    sink.into_inner().unwrap()
}

/// fasterp's message for an error, for comparison with the reference.
pub fn fasterp_message(e: &Error) -> String {
    let r = e.record().unwrap_or(0);
    match e.kind() {
        ErrorKind::Header => format!("Invalid FASTQ header at record {r}"),
        ErrorKind::Separator => format!("Invalid FASTQ separator at record {r}"),
        ErrorKind::LengthMismatch => format!("Sequence/quality length mismatch at record {r}"),
        ErrorKind::Quality => format!("Invalid quality encoding at record {r}"),
        ErrorKind::Truncated => format!("Truncated FASTQ record {r}"),
        ErrorKind::MateCount => "Paired input files have different read counts".into(),
        ErrorKind::Interleaved => "Odd read count in interleaved FASTQ".into(),
        other => format!("{other:?}: {e}"),
    }
}

/// What fasterp's owned reader does with these inputs: the FASTQ text it
/// writes per input (interleaved input: one text), or its error.
pub struct Expected {
    pub out: Vec<Vec<u8>>,
    pub error: Option<String>,
}

/// Runs fasterp's reference reader. `inputs.len()` is 1 (single or
/// interleaved) or 2 (paired).
pub fn fasterp_expected(inputs: &[&[u8]], interleaved: bool, phred64: bool) -> Expected {
    let mut first = reference(inputs[0], 7, phred64);
    let mut second = inputs.get(1).map(|d| reference(d, 13, phred64));
    let mut out = vec![Vec::new(); inputs.len()];
    loop {
        match fasterp::read_pair(&mut first, second.as_mut(), interleaved) {
            Ok(Some((r1, r2))) => {
                r1.write(&mut out[0]).unwrap();
                if let Some(r2) = r2 {
                    let o = if interleaved { 0 } else { 1 };
                    r2.write(&mut out[o]).unwrap();
                }
            }
            Ok(None) => return Expected { out, error: None },
            Err(e) => {
                return Expected {
                    out,
                    error: Some(e),
                };
            }
        }
    }
}

/// Default source options with `phred64`.
pub fn options(phred64: bool) -> FastqOptions {
    FastqOptions::default().phred64(phred64)
}

/// One record as (name without `@`, sequence, quality).
pub type Triple = (Vec<u8>, Vec<u8>, Vec<u8>);

/// Every record fasterp's owned reader yields from `data`, which must parse.
pub fn reference_records(data: &[u8], phred64: bool) -> Vec<Triple> {
    let mut r = reference(data, 17, phred64);
    let mut out = Vec::new();
    while let Some(read) = r.next().unwrap() {
        out.push((
            read.name[1..].to_vec(),
            read.seq.to_vec(),
            read.qual.to_vec(),
        ));
    }
    out
}

/// Records of segment `s` over all batches, as triples.
pub fn records(batches: &[Batch<Segments>], s: usize) -> Vec<Triple> {
    batches
        .iter()
        .flat_map(|b| {
            let c = b.segment(s);
            (0..b.len())
                .map(|i| {
                    let r = c.read(i);
                    (
                        r.name().unwrap().to_vec(),
                        r.seq().to_vec(),
                        r.qual().to_vec(),
                    )
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

/// Records of all segments in spot order (interleaved order).
pub fn spot_records(batches: &[Batch<Segments>]) -> Vec<Triple> {
    let mut out = Vec::new();
    for b in batches {
        for spot in b.spots() {
            for r in spot.iter() {
                out.push((
                    r.name().unwrap().to_vec(),
                    r.seq().to_vec(),
                    r.qual().to_vec(),
                ));
            }
        }
    }
    out
}
