//! Errors raised while building or converting batches.

use std::fmt;

/// What went wrong. Match on this rather than on messages.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum ErrorKind {
    /// A segment's sequence and quality lengths differ.
    LengthMismatch,
    /// A record has a different number of segments than the batch.
    SegmentCount,
    /// Some segments or records have names and others do not.
    Names,
    /// A column would exceed the `i32` offset limit (2 GiB).
    Overflow,
    /// A batch's segment count does not match the requested layout.
    Layout,
    /// Offsets given to [`Batch::from_parts`](crate::Batch::from_parts) are
    /// the wrong count, negative, decreasing or past the end of their bytes.
    Offsets,
}

impl ErrorKind {
    fn as_str(self) -> &'static str {
        match self {
            ErrorKind::LengthMismatch => "sequence and quality lengths differ",
            ErrorKind::SegmentCount => "wrong number of segments",
            ErrorKind::Names => "names must be present for all records or none",
            ErrorKind::Overflow => "column exceeds the i32 offset limit",
            ErrorKind::Layout => "segment count does not match layout",
            ErrorKind::Offsets => "invalid offsets",
        }
    }
}

impl fmt::Display for ErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Error from [`BatchBuilder`](crate::BatchBuilder) or
/// [`Batch::from_parts`](crate::Batch::from_parts).
///
/// Carries the record number within the batch (0-based) and the segment
/// slot when they apply.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Error {
    kind: ErrorKind,
    record: Option<usize>,
    segment: Option<usize>,
    detail: String,
}

impl Error {
    pub(crate) fn new(kind: ErrorKind, detail: impl Into<String>) -> Self {
        Self {
            kind,
            record: None,
            segment: None,
            detail: detail.into(),
        }
    }

    pub(crate) fn at(mut self, record: usize, segment: Option<usize>) -> Self {
        self.record = Some(record);
        self.segment = segment;
        self
    }

    pub(crate) fn at_segment(mut self, segment: usize) -> Self {
        self.segment = Some(segment);
        self
    }

    /// The error kind.
    pub fn kind(&self) -> ErrorKind {
        self.kind
    }

    /// 0-based record number within the batch being built, if known.
    pub fn record(&self) -> Option<usize> {
        self.record
    }

    /// 0-based segment slot, if the error concerns one segment.
    pub fn segment(&self) -> Option<usize> {
        self.segment
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.kind.as_str())?;
        match (self.record, self.segment) {
            (Some(r), Some(s)) => write!(f, " at record {r}, segment {s}")?,
            (Some(r), None) => write!(f, " at record {r}")?,
            (None, Some(s)) => write!(f, " in segment {s}")?,
            (None, None) => {}
        }
        if !self.detail.is_empty() {
            write!(f, ": {}", self.detail)?;
        }
        Ok(())
    }
}

impl std::error::Error for Error {}
