//! Small sources and sinks, and a driver.

use std::convert::Infallible;
use std::fmt;

use readkit_batch::{Batch, Layout};

use crate::{BatchSink, BatchSource, Error};

/// A source over an iterator of `Result<Batch<L>, E>`. See [`iter_source`].
pub struct IterSource<I>(I);

/// Turns anything iterable over `Result<Batch<L>, E>` into a source.
pub fn iter_source<I: IntoIterator>(iter: I) -> IterSource<I::IntoIter> {
    IterSource(iter.into_iter())
}

impl<L, I, E> BatchSource<L> for IterSource<I>
where
    I: Iterator<Item = Result<Batch<L>, E>>,
    E: std::error::Error + Send + Sync + 'static,
{
    type Error = E;
    fn next_batch(&mut self) -> Result<Option<Batch<L>>, E> {
        self.0.next().transpose()
    }
}

impl<I> fmt::Debug for IterSource<I> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("IterSource")
    }
}

/// A source calling a closure for each batch. See [`source_fn`].
pub struct FnSource<F>(F);

/// A source whose `next_batch` is `f` (for decoders written as a loop).
pub fn source_fn<F>(f: F) -> FnSource<F> {
    FnSource(f)
}

impl<L, F, E> BatchSource<L> for FnSource<F>
where
    F: FnMut() -> Result<Option<Batch<L>>, E>,
    E: std::error::Error + Send + Sync + 'static,
{
    type Error = E;
    fn next_batch(&mut self) -> Result<Option<Batch<L>>, E> {
        (self.0)()
    }
}

impl<F> fmt::Debug for FnSource<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("FnSource")
    }
}

/// A sink calling a closure for each batch. See [`sink_fn`].
pub struct FnSink<F>(F);

/// A sink whose `accept` is `f`; `finish` does nothing.
pub fn sink_fn<F>(f: F) -> FnSink<F> {
    FnSink(f)
}

impl<L, F, E> BatchSink<L> for FnSink<F>
where
    F: FnMut(Batch<L>) -> Result<(), E>,
    E: std::error::Error + Send + Sync + 'static,
{
    type Error = E;
    fn accept(&mut self, batch: Batch<L>) -> Result<(), E> {
        (self.0)(batch)
    }
    fn finish(self) -> Result<(), E> {
        Ok(())
    }
}

impl<F> fmt::Debug for FnSink<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("FnSink")
    }
}

/// A sink that keeps every batch (tests, small inputs). Holding batches
/// keeps their buffers out of the pool.
pub struct VecSink<L> {
    batches: Vec<Batch<L>>,
}

impl<L> VecSink<L> {
    /// An empty sink.
    pub fn new() -> Self {
        Self {
            batches: Vec::new(),
        }
    }

    /// The batches accepted so far.
    pub fn batches(&self) -> &[Batch<L>] {
        &self.batches
    }

    /// Takes the batches.
    pub fn into_batches(self) -> Vec<Batch<L>> {
        self.batches
    }
}

impl<L> Default for VecSink<L> {
    fn default() -> Self {
        Self::new()
    }
}

impl<L> fmt::Debug for VecSink<L> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("VecSink")
            .field("batches", &self.batches.len())
            .finish()
    }
}

impl<L> BatchSink<L> for VecSink<L> {
    type Error = Infallible;
    fn accept(&mut self, batch: Batch<L>) -> Result<(), Infallible> {
        self.batches.push(batch);
        Ok(())
    }
    fn finish(self) -> Result<(), Infallible> {
        Ok(())
    }
}

/// Drives a chain: pulls every batch from `source` into `sink`, then
/// finishes the sink. Returns the number of records moved.
///
/// Runs on the caller's thread; the stages inside `source` do their own
/// read-ahead. Source errors come back as [`ErrorKind::Upstream`](crate::ErrorKind::Upstream),
/// sink errors as [`ErrorKind::Sink`](crate::ErrorKind::Sink), each with the
/// original error as `source()`.
pub fn pump<L, S, K>(mut source: S, mut sink: K) -> Result<u64, Error>
where
    L: Layout,
    S: BatchSource<L>,
    K: BatchSink<L>,
{
    let mut records = 0u64;
    while let Some(batch) = source.next_batch().map_err(Error::upstream)? {
        records += batch.len() as u64;
        sink.accept(batch).map_err(Error::sink)?;
    }
    sink.finish().map_err(Error::sink)?;
    Ok(records)
}
