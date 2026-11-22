//! Typed errors.

use std::fmt;
use std::sync::Arc;

type Cause = Arc<dyn std::error::Error + Send + Sync + 'static>;

/// What went wrong. Match on this rather than on messages.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum ErrorKind {
    /// The run is SRA Lite: it stores synthetic qualities only, and the
    /// [`LitePolicy`](crate::LitePolicy) rejects it.
    Lite,
    /// The file, accession or table does not exist.
    NotFound,
    /// The data is inconsistent: sequence and quality lengths differ, a read
    /// lies outside its spot, the spot count is wrong, ...
    Corrupt,
    /// The run or request uses something this crate or backend cannot do
    /// (an accession type, a VDB schema, a feature not compiled in).
    Unsupported,
    /// The run's read structure does not fit the requested layout or
    /// segment selection, or changes between spots.
    Layout,
    /// Repository metadata (ENA, S3) is malformed or incomplete.
    Metadata,
    /// A network request failed. `source()` is the transport's error, for
    /// example a `readkit_fetch::Error`.
    Network,
    /// Local I/O failed.
    Io,
    /// The backend failed for another reason (an unclassified VDB error).
    Backend,
    /// A bug: a worker thread panicked.
    Internal,
}

impl ErrorKind {
    fn as_str(self) -> &'static str {
        match self {
            ErrorKind::Lite => "SRA Lite run (synthetic qualities)",
            ErrorKind::NotFound => "not found",
            ErrorKind::Corrupt => "corrupt data",
            ErrorKind::Unsupported => "unsupported",
            ErrorKind::Layout => "read structure does not fit the layout",
            ErrorKind::Metadata => "bad repository metadata",
            ErrorKind::Network => "network request failed",
            ErrorKind::Io => "I/O error",
            ErrorKind::Backend => "backend error",
            ErrorKind::Internal => "internal error",
        }
    }
}

impl fmt::Display for ErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Error returned by every fallible readkit-sra operation.
///
/// Carries the input's name (path or accession), the spot (VDB row id) and
/// segment when they apply, the backend's own status code when there is one
/// (a VDB `rc_t`), and the underlying cause. Cheap to clone.
#[derive(Clone)]
pub struct Error {
    kind: ErrorKind,
    message: String,
    input: Option<String>,
    spot: Option<u64>,
    segment: Option<usize>,
    code: Option<u32>,
    cause: Option<Cause>,
}

impl Error {
    /// An error of `kind` with a message.
    pub fn new(kind: ErrorKind, message: impl Into<String>) -> Self {
        Self {
            kind,
            message: message.into(),
            input: None,
            spot: None,
            segment: None,
            code: None,
            cause: None,
        }
    }

    /// Attach the underlying cause, returned by `source()`.
    pub fn with_source(mut self, cause: impl std::error::Error + Send + Sync + 'static) -> Self {
        self.cause = Some(Arc::new(cause));
        self
    }

    /// Attach the input's name (a path or accession).
    pub fn with_input(mut self, input: impl Into<String>) -> Self {
        self.input = Some(input.into());
        self
    }

    /// Attach the spot (row id) the error concerns.
    pub fn at_spot(mut self, spot: u64) -> Self {
        self.spot = Some(spot);
        self
    }

    /// Attach the segment slot the error concerns.
    pub fn at_segment(mut self, segment: usize) -> Self {
        self.segment = Some(segment);
        self
    }

    /// Attach the backend's status code (for VDB, an `rc_t`; for HTTP, the
    /// response status).
    pub fn with_code(mut self, code: u32) -> Self {
        self.code = Some(code);
        self
    }

    /// The kind of failure.
    pub fn kind(&self) -> ErrorKind {
        self.kind
    }

    /// The message, without kind, position or cause.
    pub fn message(&self) -> &str {
        &self.message
    }

    /// The input (path or accession), if known.
    pub fn input(&self) -> Option<&str> {
        self.input.as_deref()
    }

    /// The spot (VDB row id, 1-based for SRA runs), if the error concerns one.
    pub fn spot(&self) -> Option<u64> {
        self.spot
    }

    /// The segment slot, if the error concerns one.
    pub fn segment(&self) -> Option<usize> {
        self.segment
    }

    /// The backend's status code, if any.
    pub fn code(&self) -> Option<u32> {
        self.code
    }

    /// Whether this is [`ErrorKind::Lite`].
    pub fn is_lite(&self) -> bool {
        self.kind == ErrorKind::Lite
    }
}

impl fmt::Debug for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut d = f.debug_struct("Error");
        d.field("kind", &self.kind).field("message", &self.message);
        if let Some(i) = &self.input {
            d.field("input", i);
        }
        if let Some(s) = self.spot {
            d.field("spot", &s);
        }
        if let Some(s) = self.segment {
            d.field("segment", &s);
        }
        if let Some(c) = self.code {
            d.field("code", &c);
        }
        if let Some(c) = &self.cause {
            d.field("cause", &c.to_string());
        }
        d.finish()
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.kind.as_str())?;
        if let Some(i) = &self.input {
            write!(f, " in {i}")?;
        }
        if let Some(s) = self.spot {
            write!(f, " at spot {s}")?;
        }
        if let Some(s) = self.segment {
            write!(f, ", segment {s}")?;
        }
        if !self.message.is_empty() {
            write!(f, ": {}", self.message)?;
        }
        Ok(())
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.cause.as_deref().map(|c| c as _)
    }
}

impl From<readkit_batch::Error> for Error {
    fn from(e: readkit_batch::Error) -> Self {
        let kind = match e.kind() {
            readkit_batch::ErrorKind::LengthMismatch => ErrorKind::Corrupt,
            readkit_batch::ErrorKind::SegmentCount | readkit_batch::ErrorKind::Layout => {
                ErrorKind::Layout
            }
            _ => ErrorKind::Internal,
        };
        Error::new(kind, "building a batch").with_source(e)
    }
}

#[cfg(feature = "resolve")]
impl From<readkit_fetch::Error> for Error {
    fn from(e: readkit_fetch::Error) -> Self {
        let kind = match e.kind() {
            readkit_fetch::ErrorKind::NotFound => ErrorKind::NotFound,
            readkit_fetch::ErrorKind::Io => ErrorKind::Io,
            _ => ErrorKind::Network,
        };
        let status = e.status();
        let out = Error::new(kind, e.message().to_string()).with_source(e);
        match status {
            Some(s) => out.with_code(u32::from(s)),
            None => out,
        }
    }
}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        let kind = match e.kind() {
            std::io::ErrorKind::NotFound => ErrorKind::NotFound,
            _ => ErrorKind::Io,
        };
        Error::new(kind, e.to_string()).with_source(e)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_and_accessors() {
        let e = Error::new(ErrorKind::Corrupt, "seq 3 bytes, qual 2 bytes")
            .with_input("SRR1")
            .at_spot(7)
            .at_segment(1);
        assert_eq!(
            e.to_string(),
            "corrupt data in SRR1 at spot 7, segment 1: seq 3 bytes, qual 2 bytes"
        );
        assert_eq!(
            (e.spot(), e.segment(), e.input()),
            (Some(7), Some(1), Some("SRR1"))
        );
        assert!(!e.is_lite());
        assert!(Error::new(ErrorKind::Lite, "").is_lite());
    }

    #[test]
    fn lite_survives_a_readkit_io_chain() {
        let lite = Error::new(ErrorKind::Lite, "synthetic qualities");
        let wrapped = readkit_io::Error::upstream(lite);
        let found = wrapped.find_cause::<Error>().expect("cause kept");
        assert_eq!(found.kind(), ErrorKind::Lite);
    }
}
