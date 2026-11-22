//! The error type.

use std::fmt;
use std::sync::Arc;

type BoxError = Box<dyn std::error::Error + Send + Sync + 'static>;

/// What went wrong. Match on this rather than on messages.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum ErrorKind {
    /// A record does not start with `@`, or a non-blank line stands where a
    /// header is expected.
    Header,
    /// The third line of a record does not start with `+` (this is also how
    /// multi-line records and blank lines inside a record surface).
    Separator,
    /// Sequence and quality lengths differ.
    LengthMismatch,
    /// A quality byte lies outside `offset..=126` (offset 33, or 64 with
    /// `phred64`).
    Quality,
    /// The input ends inside a record.
    Truncated,
    /// Inputs read in lockstep have different record counts.
    MateCount,
    /// An interleaved input's record count is not a multiple of the segment
    /// count (an odd count for paired input).
    Interleaved,
    /// Reading or writing failed; the cause is `source()`.
    Io,
    /// Building a batch failed (a column over the `i32` offset limit); the
    /// cause is a `readkit_batch::Error`.
    Batch,
    /// Invalid options or arguments (for example a segment count that does
    /// not match the layout, or the wrong number of output writers).
    Options,
}

impl ErrorKind {
    fn as_str(self) -> &'static str {
        match self {
            ErrorKind::Header => "invalid FASTQ header",
            ErrorKind::Separator => "invalid FASTQ separator",
            ErrorKind::LengthMismatch => "sequence/quality length mismatch",
            ErrorKind::Quality => "invalid quality encoding",
            ErrorKind::Truncated => "truncated FASTQ record",
            ErrorKind::MateCount => "inputs read in lockstep have different record counts",
            ErrorKind::Interleaved => {
                "interleaved record count is not a multiple of the segment count"
            }
            ErrorKind::Io => "I/O error",
            ErrorKind::Batch => "building a batch failed",
            ErrorKind::Options => "invalid options",
        }
    }
}

impl fmt::Display for ErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// An error from [`FastqSource`](crate::FastqSource) or
/// [`FastqSink`](crate::FastqSink).
///
/// Input errors carry their position: the 1-based record number within the
/// input (blank lines skipped between records are not counted; in an
/// interleaved input every FASTQ record counts), the 0-based segment slot
/// when records have more than one segment, and the input's name when one
/// was given.
pub struct Error {
    kind: ErrorKind,
    record: Option<u64>,
    segment: Option<usize>,
    input: Option<Arc<str>>,
    detail: String,
    cause: Option<BoxError>,
}

impl Error {
    pub(crate) fn new(kind: ErrorKind, detail: impl Into<String>) -> Self {
        Self {
            kind,
            record: None,
            segment: None,
            input: None,
            detail: detail.into(),
            cause: None,
        }
    }

    pub(crate) fn with_cause(kind: ErrorKind, cause: impl Into<BoxError>) -> Self {
        Self {
            cause: Some(cause.into()),
            ..Self::new(kind, "")
        }
    }

    pub(crate) fn at(
        mut self,
        record: Option<u64>,
        segment: Option<usize>,
        input: Option<&Arc<str>>,
    ) -> Self {
        self.record = record;
        self.segment = segment;
        self.input = input.cloned();
        self
    }

    /// The error kind.
    pub fn kind(&self) -> ErrorKind {
        self.kind
    }

    /// The 1-based record number within the input, if the error concerns
    /// one record.
    pub fn record(&self) -> Option<u64> {
        self.record
    }

    /// The 0-based segment slot, for records with more than one segment.
    pub fn segment(&self) -> Option<usize> {
        self.segment
    }

    /// The input's name, if one was given.
    pub fn input(&self) -> Option<&str> {
        self.input.as_deref()
    }
}

impl fmt::Debug for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Error")
            .field("kind", &self.kind)
            .field("record", &self.record)
            .field("segment", &self.segment)
            .field("input", &self.input)
            .field("detail", &self.detail)
            .field("cause", &self.cause)
            .finish()
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.kind.as_str())?;
        if let Some(r) = self.record {
            write!(f, " at record {r}")?;
        }
        if let Some(s) = self.segment {
            write!(f, ", segment {s}")?;
        }
        if let Some(i) = &self.input {
            write!(f, " in {i}")?;
        }
        if !self.detail.is_empty() {
            write!(f, ": {}", self.detail)?;
        }
        Ok(())
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.cause.as_deref().map(|e| e as _)
    }
}
