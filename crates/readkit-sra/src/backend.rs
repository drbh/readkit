//! The backend boundary: what a decoder needs from an SRA reader.
//!
//! A backend opens a local `.sra` run and serves spots one at a time as
//! borrowed slices: the spot's bases and qualities (all reads
//! concatenated) plus each read's start, length and [`ReadType`]. The
//! [`Decoder`](crate::Decoder) does everything else (validation, segment
//! selection, names, batching, threads, Lite policy), so every backend gets
//! the same checks and the same output.
//!
//! - Random-access backends (VDB) report `random_access() == true`; the
//!   decoder then opens one reader per worker thread and decodes spot
//!   ranges in parallel.
//! - Streaming backends (a `fasterq-dump` subprocess, a future sracha-rs
//!   iterator) return `false` and are read sequentially, spots in order.

use std::path::Path;

use crate::error::Error;
use crate::lite::QualityKind;
use crate::read_type::ReadType;

/// An SRA decoding backend. Shared by all worker threads of a decoder.
pub trait SraBackend: Send + Sync {
    /// A short name for diagnostics (`"vdb"`, `"fasterq-dump"`, ...).
    fn name(&self) -> &'static str;

    /// Open the run at `path` for reading.
    fn open(&self, path: &Path, request: &OpenRequest) -> Result<Box<dyn SpotReader>, Error>;
}

/// What the decoder asks of a reader when opening it.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct OpenRequest {
    /// Whether spot names (VDB `NAME`) are wanted. Backends may skip
    /// reading them otherwise.
    pub names: bool,
}

impl OpenRequest {
    /// A request; `names` says whether spot names are wanted.
    pub fn new(names: bool) -> Self {
        Self { names }
    }
}

/// Facts about an opened run, reported by its backend.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub struct RunInfo {
    /// First spot id (VDB row id; 1 for SRA runs).
    pub first_spot: u64,
    /// Number of spots.
    pub spots: u64,
    /// Whether qualities are original or synthetic (SRA Lite).
    pub quality: QualityKind,
    /// Whether the backend can serve spot names (`NAME`).
    pub has_names: bool,
    /// Whether the run is a VDB database (cSRA and most newer runs) rather
    /// than a bare table. Informational.
    pub is_database: bool,
}

impl RunInfo {
    /// Info for `spots` spots starting at `first_spot`; other fields default
    /// to unknown quality, no names, not a database.
    pub fn new(first_spot: u64, spots: u64) -> Self {
        Self {
            first_spot,
            spots,
            quality: QualityKind::Unknown,
            has_names: false,
            is_database: false,
        }
    }

    /// One past the last spot id.
    pub fn end_spot(&self) -> u64 {
        self.first_spot + self.spots
    }
}

/// One spot as served by a backend, borrowed until the next read.
///
/// `starts`, `lens` and `types` have one entry per read in the spot.
/// Nothing is trusted: the decoder checks every length and range.
#[derive(Clone, Copy, Debug)]
pub struct Spot<'a> {
    /// All reads' bases, concatenated (`INSDC:dna:text`, bytes as stored).
    pub bases: &'a [u8],
    /// All reads' qualities, phred+33 text, same layout as `bases`.
    pub quals: &'a [u8],
    /// Each read's start in `bases`.
    pub starts: &'a [i32],
    /// Each read's length.
    pub lens: &'a [u32],
    /// Each read's type.
    pub types: &'a [ReadType],
    /// The spot's name (`NAME`), if requested and stored.
    pub name: Option<&'a [u8]>,
}

/// Reads spots of one opened run. One per thread.
pub trait SpotReader: Send {
    /// Facts about the run.
    fn info(&self) -> &RunInfo;

    /// Whether [`read_spot`](Self::read_spot) accepts any spot id in any
    /// order. If `false`, the decoder uses a single reader and asks for
    /// spots in increasing order, except that the most recent spot may be
    /// asked for again (the decoder reads the first spot once at open to
    /// learn the read structure, then again when decoding), so a streaming
    /// reader must keep its last spot.
    fn random_access(&self) -> bool {
        true
    }

    /// The spot with id `spot` (from `info().first_spot` to
    /// `info().end_spot()`, exclusive).
    fn read_spot(&mut self, spot: u64) -> Result<Spot<'_>, Error>;
}
