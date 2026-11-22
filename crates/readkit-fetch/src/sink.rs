//! Where downloaded bytes go.

use crate::error::{Error, ErrorKind};
use crate::object::Object;

mod file;
mod memory;
mod stream;

pub use file::FileSink;
pub use memory::MemorySink;
pub use stream::{StreamReader, StreamSink};

/// What a sink already holds when a download starts.
#[derive(Debug)]
#[non_exhaustive]
pub enum Opened {
    /// Nothing: fetch every chunk.
    Empty,
    /// A checkpoint previously passed to [`Sink::checkpoint`], whose bytes the
    /// sink still holds. The fetcher validates it and fetches only what is
    /// missing; an invalid checkpoint is ignored. Only honoured for
    /// [`Access::Random`] sinks.
    Resume(Vec<u8>),
    /// The object is already complete and verified; nothing is fetched and
    /// [`Sink::commit`] is called directly.
    Complete,
}

/// How the fetcher may write to a sink.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Access {
    /// Chunks are written in any order, concurrently, and may be read back
    /// (for verification and repair) and rewritten (repair).
    Random,
    /// Chunks are written once, in file order. The fetcher fetches at most
    /// `window` bytes ahead of the last byte written, and verifies by hashing
    /// in order; a mismatch cannot be repaired and fails the download.
    Sequential {
        /// Read-ahead in bytes beyond the write position.
        window: u64,
    },
}

/// A destination for one download.
///
/// The fetcher calls [`open`](Sink::open) once, then
/// [`write_at`](Sink::write_at) (concurrently, from several threads) and,
/// for resumable sinks, [`checkpoint`](Sink::checkpoint), and finally exactly
/// one of [`commit`](Sink::commit) (verified) or [`abort`](Sink::abort).
pub trait Sink: Send + Sync + 'static {
    /// What a successful download yields.
    type Output: Send + 'static;

    /// Prepare for `object` and report what is already present.
    fn open(&mut self, object: &Object) -> Result<Opened, Error>;

    /// Write access pattern; [`Access::Random`] by default.
    fn access(&self) -> Access {
        Access::Random
    }

    /// Store `data` at `offset`. For [`Access::Sequential`] sinks calls arrive
    /// in order and may block (backpressure) until [`interrupt`](Sink::interrupt).
    fn write_at(&self, offset: u64, data: Vec<u8>) -> Result<(), Error>;

    /// Fill `buf` from `offset`. Required for [`Access::Random`] sinks.
    fn read_at(&self, offset: u64, buf: &mut [u8]) -> Result<(), Error> {
        let _ = (offset, buf);
        Err(Error::new(ErrorKind::Io, "sink does not support reads"))
    }

    /// Persist resume state (opaque bytes, at most about once per second).
    /// Sinks that cannot resume ignore it.
    fn checkpoint(&self, state: &[u8]) -> Result<(), Error> {
        let _ = state;
        Ok(())
    }

    /// Called from another thread when the download fails or is cancelled:
    /// a blocked [`write_at`](Sink::write_at) should return promptly.
    fn interrupt(&self) {}

    /// Every byte is written and verified: finish and produce the output.
    fn commit(self) -> Result<Self::Output, Error>;

    /// The download failed with `error`.
    fn abort(self, error: &Error);
}
