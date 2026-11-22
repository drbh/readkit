//! readkit: zero-copy batches of sequencing reads, and the stages that produce
//! and consume them.
//!
//! This crate re-exports the readkit crates behind features. The core is
//! always available:
//!
//! - [`batch`] (`readkit-batch`): the columnar, Arrow-compatible [`Batch`] of
//!   reads with one or more segments per record, built from a recycling [`Pool`].
//! - [`io`] (`readkit-io`): [`BatchSource`] and [`BatchSink`], and the
//!   [`ReadAhead`] helper tools use to read ahead of their input.
//!
//! Formats and transport are features:
//!
//! | Feature | Module | Crate |
//! | --- | --- | --- |
//! | `fastq` (default) | `fastq` | strict FASTQ parsing and writing |
//! | `gz` (default) | `gz` | gzip decoding and encoding |
//! | `parallel` (default) | — | parallel gzip (threads, libdeflate) |
//! | `fetch` | `fetch` | ranged, verified HTTP downloads |
//! | `sra` | `sra` | `.sra` decoding into batches (no HTTP stack) |
//! | `resolve` | `sra::Resolver` | SRA accession resolution via ENA and NCBI S3 (implies `sra`) |
//! | `ngs` | — | in-process `.sra` decoding via the NCBI VDB library (implies `sra`) |
//! | `arrow` | `batch::arrow` | zero-copy export to arrow-rs |
//!
//! Stages compose by pulling batches from each other. Copying paired FASTQ
//! through a batch source and sink:
//!
//! ```
//! # #[cfg(feature = "fastq")]
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! use readkit::fastq::{FastqInput, FastqOptions, FastqSink, FastqSource};
//!
//! let r1 = b"@a/1\nACGT\n+\nIIII\n".as_slice();
//! let r2 = b"@a/2\nTTGG\n+\nIIII\n".as_slice();
//! let source = FastqSource::paired(FastqInput::new(r1), FastqInput::new(r2), FastqOptions::default());
//! let mut sink = FastqSink::split(vec![Vec::new(), Vec::new()]);
//! let records = readkit::pump(source, &mut sink)?;
//!
//! assert_eq!(records, 1);
//! let out = sink.into_inner()?;
//! assert_eq!(out[0], b"@a/1\nACGT\n+\nIIII\n");
//! assert_eq!(out[1], b"@a/2\nTTGG\n+\nIIII\n");
//! # Ok(())
//! # }
//! # #[cfg(not(feature = "fastq"))]
//! # fn main() {}
//! ```
#![deny(missing_docs)]

pub use readkit_batch as batch;
pub use readkit_io as io;

#[cfg(feature = "fastq")]
pub use readkit_fastq as fastq;

#[cfg(feature = "gz")]
pub use readkit_gz as gz;

#[cfg(feature = "fetch")]
pub use readkit_fetch as fetch;

#[cfg(feature = "sra")]
pub use readkit_sra as sra;

pub use readkit_batch::{Batch, BatchBuilder, Paired, Pool, Read, Segments, Single};
pub use readkit_io::{BatchSink, BatchSource, ReadAhead, pump};
