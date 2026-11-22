//! Columnar, pooled batches of sequencing reads.
//!
//! A [`Batch<L>`] holds `n` records of `k` segments each. Every segment slot
//! stores names (optional for the whole batch), sequences and qualities in
//! Arrow `Binary` layout: `n + 1` `i32` offsets plus a byte buffer. Sequence
//! and quality share one offsets array. Buffers come from a [`Pool`] and go
//! back to it when the batch drops. An optional `u16` tag per record lets
//! tools carry facts across a boundary; readkit gives tags no meaning.
//!
//! The layout `L` is [`Single`] (k = 1), [`Paired`] (k = 2) or [`Segments`]
//! (k chosen at runtime). Every batch is stored as runtime-`k` segments; the
//! typed layouts only add methods that make sense for them, so
//! `pairs()` on a single-end batch does not compile.
//!
//! ```
//! use readkit_batch::{BatchBuilder, Pool, Read, Single};
//! let pool = Pool::new();
//! let mut b = BatchBuilder::<Single>::with_pool(&pool, 2);
//! b.push_read(Read::named(b"a", b"ACGT", b"IIII"))?;
//! b.push_read(Read::named(b"b", b"", b""))?;
//! let batch = b.finish();
//! let seqs = batch.column().seqs();
//! assert_eq!(seqs.offsets(), &[0, 4, 4]);
//! assert_eq!(seqs.data(), b"ACGT");
//! # Ok::<(), readkit_batch::Error>(())
//! ```
//!
//! Feature `arrow` adds zero-copy export to arrow-rs (see the `arrow`
//! module).

#![deny(missing_docs)]

mod batch;
mod builder;
mod error;
mod layout;
mod pool;

#[cfg(feature = "arrow")]
pub mod arrow;

pub use batch::{Batch, BatchParts, Binary, Column, LayoutError, Read, SegmentParts, Spot};
pub use builder::BatchBuilder;
pub use error::{Error, ErrorKind};
pub use layout::{FixedLayout, Layout, Paired, Segments, Single};
pub use pool::{DEFAULT_MAX_IDLE, Element, Pool, PoolStats, PooledVec};
