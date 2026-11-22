//! Typed errors.

use std::fmt;
use std::sync::Arc;

/// What went wrong, for callers to match on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum ErrorKind {
    /// The server says the object does not exist (HTTP 404 or 410).
    NotFound,
    /// A client error that retrying cannot fix (4xx other than 404, 408, 410
    /// and 429), for example 403.
    Permanent,
    /// The server asked us to slow down (HTTP 429 or 503).
    Throttled,
    /// A failure that may succeed on retry: a 5xx, a dropped connection, a
    /// short body, or an unexpected status for a ranged request.
    Transient,
    /// A request timed out, or the whole transfer made no forward progress for
    /// [`FetchOptions::no_progress_timeout`](crate::FetchOptions::no_progress_timeout).
    Timeout,
    /// The downloaded bytes do not match the expected checksum, even after
    /// repair.
    Verification,
    /// Local I/O failed (disk, sink).
    Io,
    /// The download was cancelled, or a stream reader was dropped.
    Cancelled,
    /// The request or options are inconsistent (bad URL, impossible part
    /// layout, malformed ETag).
    InvalidInput,
    /// A response is larger than the caller's limit
    /// ([`GetOptions::max_bytes`](crate::GetOptions::max_bytes)).
    TooLarge,
    /// A bug: a worker panicked.
    Internal,
}

impl ErrorKind {
    /// Whether a chunk failing with this kind is retried.
    pub fn is_retryable(self) -> bool {
        matches!(self, Self::Throttled | Self::Transient | Self::Timeout)
    }
}

/// Error returned by every fallible readkit-fetch operation.
///
/// Cheap to clone; the source, if any, is shared.
#[derive(Clone)]
pub struct Error {
    kind: ErrorKind,
    message: String,
    status: Option<u16>,
    source: Option<Arc<dyn std::error::Error + Send + Sync + 'static>>,
}

impl Error {
    /// A new error of `kind` with a human-readable `message`.
    pub fn new(kind: ErrorKind, message: impl Into<String>) -> Self {
        Self {
            kind,
            message: message.into(),
            status: None,
            source: None,
        }
    }

    /// Attach the underlying cause, returned by [`std::error::Error::source`].
    pub fn with_source(mut self, source: impl std::error::Error + Send + Sync + 'static) -> Self {
        self.source = Some(Arc::new(source));
        self
    }

    pub(crate) fn with_status(mut self, status: u16) -> Self {
        self.status = Some(status);
        self
    }

    /// The kind of failure.
    pub fn kind(&self) -> ErrorKind {
        self.kind
    }

    /// The HTTP status that caused this error, if any.
    pub fn status(&self) -> Option<u16> {
        self.status
    }

    /// The message, without the source chain.
    pub fn message(&self) -> &str {
        &self.message
    }

    pub(crate) fn cancelled() -> Self {
        Self::new(ErrorKind::Cancelled, "download cancelled")
    }
}

impl fmt::Debug for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut d = f.debug_struct("Error");
        d.field("kind", &self.kind).field("message", &self.message);
        if let Some(s) = self.status {
            d.field("status", &s);
        }
        if let Some(s) = &self.source {
            d.field("source", &s.to_string());
        }
        d.finish()
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.source
            .as_deref()
            .map(|s| s as &(dyn std::error::Error + 'static))
    }
}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Self::new(ErrorKind::Io, e.to_string()).with_source(e)
    }
}

impl From<Error> for std::io::Error {
    fn from(e: Error) -> Self {
        use std::io::ErrorKind as Io;
        let kind = match e.kind {
            ErrorKind::NotFound => Io::NotFound,
            ErrorKind::Timeout => Io::TimedOut,
            ErrorKind::Verification => Io::InvalidData,
            ErrorKind::InvalidInput => Io::InvalidInput,
            _ => Io::Other,
        };
        std::io::Error::new(kind, e)
    }
}
