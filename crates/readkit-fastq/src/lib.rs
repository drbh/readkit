//! Strict FASTQ parsing into [`Batch`]es and FASTQ writing out of them: the
//! text edge of readkit.
//!
//! [`FastqSource`] reads one or more [`BufRead`](std::io::BufRead) inputs:
//! `k` files in lockstep (single-end, R1/R2, single-cell R1/R2/I1/I2, ...)
//! or one interleaved file of `k`-record spots. It implements
//! [`BatchSource`](readkit_io::BatchSource) and fills batches from a
//! [`Pool`] under one byte budget across all inputs.
//! [`FastqSink`] writes batches back as FASTQ and implements
//! [`BatchSink`](readkit_io::BatchSink).
//!
//! The parser is fasterp's strict reader, generalized from two segments to
//! `k`:
//!
//! - every record has an `@` header and a `+` separator; the `+` line's
//!   content is not compared to the name, and is not stored (the sink writes
//!   a bare `+` or repeats the name, see [`PlusLine`]);
//! - sequence and quality lengths are equal, quality bytes lie in
//!   `offset..=126` (offset 33, or 64 with [`FastqOptions::phred64`], which
//!   converts qualities to Phred+33);
//! - one trailing `\r` per line is stripped; multi-line records are rejected;
//!   empty sequences are allowed;
//! - blank and whitespace-only lines where a header is expected are skipped;
//!   a blank line inside a record is an error;
//! - truncated records, record-count mismatches between lockstep inputs and
//!   interleaved counts that are not a multiple of `k` are errors.
//!
//! Errors name the record number, the segment and the input (see [`Error`]).
//! No decompression happens here: wrap a decoder (readkit-gz) in a
//! `BufReader` and pass that.
//!
//! ```
//! use readkit_batch::Paired;
//! use readkit_fastq::{FastqInput, FastqOptions, FastqSink, FastqSource};
//! use readkit_io::pump;
//!
//! let r1 = b"@a/1\nACGT\n+\nIIII\n@b/1\nGG\n+\nII\n";
//! let r2 = b"@a/2\nTT\n+\nII\n@b/2\nCCA\n+\nIII\n";
//! let source = FastqSource::<Paired>::paired(
//!     FastqInput::named(&r1[..], "r1.fq"),
//!     FastqInput::named(&r2[..], "r2.fq"),
//!     FastqOptions::default(),
//! );
//! let mut out = [Vec::new(), Vec::new()];
//! let [o1, o2] = &mut out;
//! let records = pump(source, FastqSink::split(vec![o1, o2]))?;
//! assert_eq!(records, 2);
//! assert_eq!(out[0], r1);
//! assert_eq!(out[1], r2);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

#![deny(missing_docs)]

mod error;
#[cfg(not(target_arch = "wasm32"))]
mod par;
mod sink;
mod source;

pub use error::{Error, ErrorKind};
pub use sink::{FastqSink, PlusLine};
pub use source::{FastqInput, FastqSource};

use readkit_batch::{Batch, Layout, Paired, Pool, Segments, Single};

/// Default byte budget of one batch: 8 MiB of FASTQ text across all inputs
/// (fasterp's effective batch size).
pub const DEFAULT_BATCH_BYTES: usize = 8 << 20;

/// Options for [`FastqSource`].
///
/// ```
/// use readkit_fastq::FastqOptions;
/// let opts = FastqOptions::default().phred64(true).keep_names(false);
/// assert!(opts.phred64 && !opts.keep_names);
/// ```
#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct FastqOptions {
    /// Input qualities are Phred+64: accept bytes in `64..=126` and convert
    /// them to Phred+33 in the batch (fasterp's `--phred64`). Default false.
    pub phred64: bool,
    /// Store read names (the header line without `@`). Default true.
    pub keep_names: bool,
    /// At most this many bytes of FASTQ text per batch, summed over all
    /// inputs: a batch stops before the spot that would take it over, but
    /// always takes at least one. Batches therefore depend only on the text
    /// and these options, not on input buffering or
    /// [`threads`](Self::threads). Default [`DEFAULT_BATCH_BYTES`].
    pub batch_bytes: usize,
    /// At most this many spots per batch. Default `usize::MAX`.
    pub batch_records: usize,
    /// The pool batches draw their buffers from. Default: a new pool per
    /// source.
    pub pool: Option<Pool>,
    /// Worker threads that validate records and copy them into batches.
    /// Default 0: everything on the caller's thread. With `n > 0`, a thread
    /// per input reads it and finds record ends, one more cuts batches, and
    /// `n` workers fill them: on short reads in memory, 2.5× the
    /// single-thread rate at 4 workers and 3× at 14. Batches, records and
    /// errors are the same for every value. Ignored on wasm32.
    ///
    /// The threads take plain four-line records only. From a blank line
    /// between records, and for the last batch of the input, the source
    /// reads on the caller's thread to the end.
    pub threads: usize,
}

impl Default for FastqOptions {
    fn default() -> Self {
        Self {
            phred64: false,
            keep_names: true,
            batch_bytes: DEFAULT_BATCH_BYTES,
            batch_records: usize::MAX,
            pool: None,
            threads: 0,
        }
    }
}

impl FastqOptions {
    /// Sets [`phred64`](Self::phred64).
    pub fn phred64(mut self, on: bool) -> Self {
        self.phred64 = on;
        self
    }

    /// Sets [`keep_names`](Self::keep_names).
    pub fn keep_names(mut self, on: bool) -> Self {
        self.keep_names = on;
        self
    }

    /// Sets [`batch_bytes`](Self::batch_bytes).
    pub fn batch_bytes(mut self, bytes: usize) -> Self {
        self.batch_bytes = bytes;
        self
    }

    /// Sets [`batch_records`](Self::batch_records).
    pub fn batch_records(mut self, records: usize) -> Self {
        self.batch_records = records;
        self
    }

    /// Sets [`threads`](Self::threads).
    pub fn threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    /// Sets [`pool`](Self::pool).
    pub fn pool(mut self, pool: &Pool) -> Self {
        self.pool = Some(pool.clone());
        self
    }
}

mod sealed {
    pub trait Sealed {}
    impl Sealed for readkit_batch::Single {}
    impl Sealed for readkit_batch::Paired {}
    impl Sealed for readkit_batch::Segments {}
}

/// The batch layouts a [`FastqSource`] produces: [`Single`], [`Paired`] and
/// [`Segments`]. Sealed.
pub trait FastqLayout: Layout + sealed::Sealed {
    #[doc(hidden)]
    fn from_segments(batch: Batch<Segments>) -> Batch<Self>
    where
        Self: Sized;
}

impl FastqLayout for Single {
    fn from_segments(batch: Batch<Segments>) -> Batch<Self> {
        batch
            .try_into_layout()
            .unwrap_or_else(|_| unreachable!("segment count checked on construction"))
    }
}

impl FastqLayout for Paired {
    fn from_segments(batch: Batch<Segments>) -> Batch<Self> {
        batch
            .try_into_layout()
            .unwrap_or_else(|_| unreachable!("segment count checked on construction"))
    }
}

impl FastqLayout for Segments {
    fn from_segments(batch: Batch<Segments>) -> Batch<Self> {
        batch
    }
}
