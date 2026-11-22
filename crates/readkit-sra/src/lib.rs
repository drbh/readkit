//! SRA for readkit: `.sra` decoding into [`Batch`](readkit_batch::Batch)es,
//! typed SRA Lite detection, and (feature `resolve`) accession resolution.
//!
//! | Feature | Adds | Dependencies |
//! | --- | --- | --- |
//! | (none) | decoding through the `fasterq-dump` backend | no HTTP stack |
//! | `ngs` | the in-process NCBI VDB backend (`NgsBackend`) | `libncbi-ngs` at build time |
//! | `resolve` | `Resolver`, `Transport`, `HttpTransport` | readkit-fetch (reqwest, rustls, tokio) |
//!
//! `ngs` and `resolve` are independent; neither is on by default.
//!
//! # Decoding
//!
//! A [`Decoder<L>`] opens a local `.sra` file and is a
//! [`BatchSource<L>`](readkit_io::BatchSource). Bases and qualities go from
//! the backend straight into the batch's sequence and quality columns; no
//! FASTQ text is formatted. Names are optional ([`Names`]).
//!
//! - Layouts: [`SingleDecoder`], [`PairedDecoder`], or [`SegmentsDecoder`]
//!   for any number of segments (single-cell).
//! - Which reads of a spot become segments: [`SegmentSelection`] (biological
//!   only by default; all reads, technical included; or chosen indices).
//!   Read types keep their orientation bits ([`ReadType`]).
//! - SRA Lite: [`QualityKind`] in [`RunInfo`], and [`LitePolicy`]
//!   (rejecting with [`ErrorKind::Lite`] by default).
//! - Backends sit behind [`SraBackend`]; this crate ships the NCBI VDB
//!   backend (feature `ngs`, `NgsBackend`; it decodes bases and qualities of
//!   common unaligned runs from the stored blobs itself, about 2× faster
//!   than through VDB) and a `fasterq-dump` subprocess
//!   backend ([`FasterqDumpBackend`], always built). [`default_backend`]
//!   picks VDB when built with `ngs`, otherwise `fasterq-dump` if it is
//!   installed.
//!
//! ```no_run
//! use readkit_io::BatchSource;
//! use readkit_sra::{DecodeOptions, ErrorKind, SingleDecoder};
//!
//! match SingleDecoder::open("SRR12999707", DecodeOptions::default()) {
//!     Ok(mut d) => {
//!         while let Some(batch) = d.next_batch()? {
//!             let seqs = batch.column().seqs(); // one buffer + offsets
//!             let _ = (seqs.offsets(), seqs.data());
//!         }
//!     }
//!     Err(e) if e.kind() == ErrorKind::Lite => { /* use ENA FASTQ instead */ }
//!     Err(e) => return Err(e),
//! }
//! # Ok::<(), readkit_sra::Error>(())
//! ```
//!
//! # Resolution (feature `resolve`)
//!
//! `Resolver` turns run, experiment, sample, study and project accessions
//! into `RunRecord`s: the `.sra` object on NCBI S3 (size and ETag, as a
//! `readkit_fetch::Object`) and ENA's FASTQ files (MD5s, mates classified
//! by [`classify_mate`]). All HTTP goes through readkit-fetch, behind the
//! `Transport` trait so that tests and mirrors can stay offline.

#![deny(missing_docs)]

mod backend;
mod decode;
mod error;
mod fasterq;
mod lite;
mod mates;
#[cfg(feature = "ngs")]
mod ngs;
mod rc;
mod read_type;
#[cfg(feature = "resolve")]
mod resolve;
#[cfg(feature = "ngs")]
mod vdbfmt;

pub use backend::{OpenRequest, RunInfo, Spot, SpotReader, SraBackend};
pub use decode::{
    DecodeLayout, DecodeOptions, DecodeStats, Decoder, Names, PairedDecoder, SegmentSelection,
    SegmentsDecoder, SingleDecoder, default_backend, run_info,
};
pub use error::{Error, ErrorKind};
pub use fasterq::{FasterqDumpBackend, FasterqDumpReader};
pub use lite::{LITE_SAMPLE, LitePolicy, QualityKind, QualityProbe};
pub use mates::{Mate, classify_mate};
#[cfg(feature = "ngs")]
pub use ngs::{NgsBackend, NgsReader};
pub use read_type::{Orientation, ReadType};
#[cfg(feature = "resolve")]
pub use resolve::{
    AccessionKind, ENA_FIELDS, ENA_PORTAL, FastqFile, FastqGroup, HttpTransport, LibraryLayout,
    NCBI_S3, Resolver, RunRecord, SraStatus, Transport, parse_ena_rows, run_from_row,
};
