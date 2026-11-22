//! Stages that pass [`Batch`]es: the [`BatchSource`] and [`BatchSink`]
//! traits, the [`ReadAhead`] helper, and small adapters.
//!
//! Chains are pulled from the end: the last consumer drives, and each stage
//! pulls from the one before it. A stage that consumes a source wraps it in
//! [`ReadAhead`], which reads at most one batch ahead on its own thread
//! (sequentially on wasm32), so stages overlap without pipeline authors
//! placing thread boundaries. Cancellation is drop. Stages keep upstream
//! errors as their `source()`.
//!
//! ```
//! use readkit_batch::{Batch, BatchBuilder, Pool, Read, Single};
//! use readkit_io::{BatchSource, ReadAhead, VecSink, iter_source, pump};
//!
//! let pool = Pool::new();
//! let batches = (0..3).map(|i| {
//!     let mut b = BatchBuilder::<Single>::new(&pool);
//!     b.push_read(Read::new(b"ACGT", b"IIII"))?;
//!     b.set_tag(i);
//!     Ok::<_, readkit_batch::Error>(b.finish())
//! });
//! let source = ReadAhead::new(iter_source(batches.collect::<Vec<_>>()));
//! let mut sink = VecSink::new();
//! let records = pump(source, &mut sink)?;
//! assert_eq!(records, 3);
//! assert_eq!(sink.batches()[2].tags(), Some(&[2][..]));
//! # Ok::<(), readkit_io::Error>(())
//! ```

#![deny(missing_docs)]

mod adapters;
mod error;
mod read_ahead;

pub use adapters::{FnSink, FnSource, IterSource, VecSink, iter_source, pump, sink_fn, source_fn};
pub use error::{Error, ErrorKind, find_cause};
pub use read_ahead::ReadAhead;

use readkit_batch::Batch;

/// Produces batches. Pulled by the next stage.
pub trait BatchSource<L> {
    /// The source's error.
    type Error: std::error::Error + Send + Sync + 'static;

    /// The next batch, or `None` at the end of input.
    fn next_batch(&mut self) -> Result<Option<Batch<L>>, Self::Error>;
}

/// Consumes batches.
pub trait BatchSink<L> {
    /// The sink's error.
    type Error: std::error::Error + Send + Sync + 'static;

    /// Takes one batch.
    fn accept(&mut self, batch: Batch<L>) -> Result<(), Self::Error>;

    /// Flushes and closes the sink.
    fn finish(self) -> Result<(), Self::Error>
    where
        Self: Sized;
}

impl<L, S: BatchSource<L> + ?Sized> BatchSource<L> for &mut S {
    type Error = S::Error;
    fn next_batch(&mut self) -> Result<Option<Batch<L>>, S::Error> {
        (**self).next_batch()
    }
}

impl<L, S: BatchSource<L> + ?Sized> BatchSource<L> for Box<S> {
    type Error = S::Error;
    fn next_batch(&mut self) -> Result<Option<Batch<L>>, S::Error> {
        (**self).next_batch()
    }
}

/// `finish` on a borrowed sink does nothing; the owner finishes it.
impl<L, K: BatchSink<L> + ?Sized> BatchSink<L> for &mut K {
    type Error = K::Error;
    fn accept(&mut self, batch: Batch<L>) -> Result<(), K::Error> {
        (**self).accept(batch)
    }
    fn finish(self) -> Result<(), K::Error> {
        Ok(())
    }
}
