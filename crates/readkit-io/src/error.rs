//! The error type for stages built on readkit-io.

use std::borrow::Cow;
use std::fmt;

type BoxError = Box<dyn std::error::Error + Send + Sync + 'static>;

/// Where in a chain an error arose. Match on this rather than on messages.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum ErrorKind {
    /// The upstream source failed; the cause is [`source()`](std::error::Error::source).
    Upstream,
    /// The downstream sink failed; the cause is `source()`.
    Sink,
    /// The stage itself failed (bad input, invalid options, ...).
    Stage,
    /// Building a batch failed; `source()` is a `readkit_batch::Error`.
    Batch,
}

impl ErrorKind {
    fn as_str(self) -> &'static str {
        match self {
            ErrorKind::Upstream => "upstream source failed",
            ErrorKind::Sink => "sink failed",
            ErrorKind::Stage => "stage failed",
            ErrorKind::Batch => "building a batch failed",
        }
    }
}

impl fmt::Display for ErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// An error with a kind, an optional message and an optional cause.
///
/// Stages that wrap a source keep its error as `source()`, so a pipeline can
/// find the root cause through several stages with [`find_cause`].
pub struct Error {
    kind: ErrorKind,
    message: Cow<'static, str>,
    cause: Option<BoxError>,
}

impl Error {
    /// An error of `kind` with a message and no cause.
    pub fn new(kind: ErrorKind, message: impl Into<Cow<'static, str>>) -> Self {
        Self {
            kind,
            message: message.into(),
            cause: None,
        }
    }

    /// An error of `kind` caused by `cause`.
    pub fn with_cause(
        kind: ErrorKind,
        message: impl Into<Cow<'static, str>>,
        cause: impl Into<BoxError>,
    ) -> Self {
        Self {
            kind,
            message: message.into(),
            cause: Some(cause.into()),
        }
    }

    /// The upstream source failed with `cause`.
    pub fn upstream(cause: impl Into<BoxError>) -> Self {
        Self::with_cause(ErrorKind::Upstream, "", cause)
    }

    /// The sink failed with `cause`.
    pub fn sink(cause: impl Into<BoxError>) -> Self {
        Self::with_cause(ErrorKind::Sink, "", cause)
    }

    /// The error kind.
    pub fn kind(&self) -> ErrorKind {
        self.kind
    }

    /// The message, possibly empty.
    pub fn message(&self) -> &str {
        &self.message
    }

    /// The first error of type `E` in this error's `source()` chain,
    /// including itself.
    pub fn find_cause<E: std::error::Error + 'static>(&self) -> Option<&E> {
        find_cause(self)
    }
}

impl From<readkit_batch::Error> for Error {
    fn from(e: readkit_batch::Error) -> Self {
        Self::with_cause(ErrorKind::Batch, "", e)
    }
}

impl fmt::Debug for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Error")
            .field("kind", &self.kind)
            .field("message", &self.message)
            .field("cause", &self.cause)
            .finish()
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.kind.as_str())?;
        if !self.message.is_empty() {
            write!(f, ": {}", self.message)?;
        }
        Ok(())
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.cause.as_deref().map(|e| e as _)
    }
}

/// The first error of type `E` in `err`'s `source()` chain, including `err`.
pub fn find_cause<'a, E: std::error::Error + 'static>(
    err: &'a (dyn std::error::Error + 'static),
) -> Option<&'a E> {
    let mut cur = Some(err);
    while let Some(e) = cur {
        if let Some(hit) = e.downcast_ref::<E>() {
            return Some(hit);
        }
        cur = e.source();
    }
    None
}
