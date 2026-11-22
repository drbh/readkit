//! Typed errors.

use std::fmt;
use std::io;
use std::sync::Arc;

/// What went wrong, for callers to match on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum ErrorKind {
    /// The input is not valid gzip: a bad member header, invalid deflate
    /// data, or bytes after a member that do not start another member.
    Corrupt,
    /// A member decoded, but its CRC32 or length (ISIZE) does not match its
    /// trailer, or its optional header CRC16 does not match the header.
    Checksum,
    /// The input ended inside a member (including an empty input).
    Truncated,
    /// Reading the compressed input failed. The I/O error is the
    /// [`source`](std::error::Error::source).
    Io,
    /// The options are unusable (for example a zero buffer size).
    InvalidInput,
    /// A bug: a worker thread panicked or stopped without a result.
    Internal,
}

/// Error returned by every fallible readkit-gz operation.
///
/// Cheap to clone; the source, if any, is shared. Readers keep returning the
/// same error after the first failure.
///
/// At the [`std::io::Read`] and [`std::io::Write`] boundary it is converted
/// into an [`io::Error`] that carries it as its inner error, so callers can
/// recover it:
///
/// ```
/// # fn f(e: std::io::Error) {
/// if let Some(gz) = readkit_gz::Error::from_io(&e) {
///     println!("{:?} in the member at byte {:?}", gz.kind(), gz.offset());
/// }
/// # }
/// ```
#[derive(Clone)]
pub struct Error {
    kind: ErrorKind,
    message: String,
    offset: Option<u64>,
    source: Option<Arc<dyn std::error::Error + Send + Sync + 'static>>,
}

impl Error {
    /// A new error of `kind` with a human-readable `message`.
    pub fn new(kind: ErrorKind, message: impl Into<String>) -> Self {
        Self {
            kind,
            message: message.into(),
            offset: None,
            source: None,
        }
    }

    /// Attach the underlying cause, returned by [`std::error::Error::source`].
    pub fn with_source(mut self, source: impl std::error::Error + Send + Sync + 'static) -> Self {
        self.source = Some(Arc::new(source));
        self
    }

    /// Set the compressed byte offset of the gzip member the error is in,
    /// unless one is already set.
    pub(crate) fn at(mut self, offset: u64) -> Self {
        self.offset.get_or_insert(offset);
        self
    }

    /// The kind of failure.
    pub fn kind(&self) -> ErrorKind {
        self.kind
    }

    /// The message, without the source chain.
    pub fn message(&self) -> &str {
        &self.message
    }

    /// Compressed byte offset at which the gzip member that failed starts,
    /// when known.
    pub fn offset(&self) -> Option<u64> {
        self.offset
    }

    /// The readkit-gz error inside an [`io::Error`] returned by one of this
    /// crate's readers or writers, if it carries one.
    pub fn from_io(error: &io::Error) -> Option<&Error> {
        error.get_ref()?.downcast_ref::<Error>()
    }

    pub(crate) fn corrupt(message: impl Into<String>) -> Self {
        Self::new(ErrorKind::Corrupt, message)
    }

    pub(crate) fn checksum(message: impl Into<String>) -> Self {
        Self::new(ErrorKind::Checksum, message)
    }

    pub(crate) fn truncated() -> Self {
        Self::new(ErrorKind::Truncated, "gzip input ends inside a member")
    }

    #[cfg(feature = "parallel")]
    pub(crate) fn internal(message: impl Into<String>) -> Self {
        Self::new(ErrorKind::Internal, message)
    }

    /// An error from reading compressed input: end of file means truncated
    /// input; a readkit-gz error passed through an `io::Error` is unwrapped.
    pub(crate) fn input(e: io::Error) -> Self {
        if e.kind() == io::ErrorKind::UnexpectedEof {
            return Self::truncated().with_source(e);
        }
        match e.get_ref().and_then(|i| i.downcast_ref::<Error>()) {
            Some(inner) => inner.clone(),
            None => Self::new(ErrorKind::Io, format!("reading gzip input: {e}")).with_source(e),
        }
    }
}

impl fmt::Debug for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut d = f.debug_struct("Error");
        d.field("kind", &self.kind).field("message", &self.message);
        if let Some(o) = self.offset {
            d.field("offset", &o);
        }
        if let Some(s) = &self.source {
            d.field("source", &s.to_string());
        }
        d.finish()
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)?;
        if let Some(o) = self.offset {
            write!(f, " (gzip member at byte {o})")?;
        }
        Ok(())
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.source
            .as_deref()
            .map(|s| s as &(dyn std::error::Error + 'static))
    }
}

impl From<Error> for io::Error {
    fn from(e: Error) -> Self {
        let kind = match e.kind {
            ErrorKind::Corrupt | ErrorKind::Checksum => io::ErrorKind::InvalidData,
            ErrorKind::Truncated => io::ErrorKind::UnexpectedEof,
            ErrorKind::InvalidInput => io::ErrorKind::InvalidInput,
            ErrorKind::Io => e
                .source
                .as_deref()
                .and_then(|s| s.downcast_ref::<io::Error>())
                .map(io::Error::kind)
                // A retryable kind would make callers retry a reader that
                // keeps failing.
                .filter(|k| !matches!(k, io::ErrorKind::Interrupted | io::ErrorKind::WouldBlock))
                .unwrap_or(io::ErrorKind::Other),
            ErrorKind::Internal => io::ErrorKind::Other,
        };
        io::Error::new(kind, e)
    }
}
