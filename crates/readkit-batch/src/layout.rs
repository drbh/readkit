//! Layout markers: how many segments each record has.

mod sealed {
    pub trait Sealed {}
    impl Sealed for super::Single {}
    impl Sealed for super::Paired {}
    impl Sealed for super::Segments {}
}

/// Single-end reads: one segment per record (`k = 1`).
///
/// Methods that need another segment count do not exist on single-end
/// batches and builders:
///
/// ```compile_fail
/// use readkit_batch::{BatchBuilder, Pool, Single};
/// let batch = BatchBuilder::<Single>::new(&Pool::new()).finish();
/// let _ = batch.pairs(); // no `pairs` on Batch<Single>
/// ```
///
/// ```compile_fail
/// use readkit_batch::{BatchBuilder, Pool, Read, Single};
/// let mut b = BatchBuilder::<Single>::new(&Pool::new());
/// let r = Read::new(b"A", b"I");
/// b.push_pair(r, r); // no `push_pair` on a single-end builder
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Single;

/// Paired reads: two segments per record (`k = 2`).
///
/// ```
/// use readkit_batch::{BatchBuilder, Paired, Pool};
/// let batch = BatchBuilder::<Paired>::new(&Pool::new()).finish();
/// assert_eq!(batch.pairs().count(), 0);
/// ```
///
/// ```compile_fail
/// use readkit_batch::{BatchBuilder, Paired, Pool};
/// let batch = BatchBuilder::<Paired>::new(&Pool::new()).finish();
/// let _ = batch.reads(); // no `reads` on Batch<Paired>
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Paired;

/// Any number of segments chosen at runtime (for example single-cell
/// barcode + cDNA + index reads). This is the internal representation of
/// every batch; [`Single`] and [`Paired`] are typed views over it.
///
/// A runtime-`k` builder needs `k`; it has no `with_pool`:
///
/// ```compile_fail
/// use readkit_batch::{BatchBuilder, Pool, Segments};
/// let b = BatchBuilder::<Segments>::with_pool(&Pool::new(), 8);
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Segments;

/// A batch layout marker. Sealed: implemented by [`Single`], [`Paired`] and
/// [`Segments`] only.
pub trait Layout: sealed::Sealed + Send + Sync + 'static {
    /// The segment count fixed by the type, or `None` for [`Segments`].
    const SEGMENTS: Option<usize>;
}

/// A layout whose segment count is known at compile time.
pub trait FixedLayout: Layout {
    /// The segment count.
    const K: usize;
}

impl Layout for Single {
    const SEGMENTS: Option<usize> = Some(1);
}
impl Layout for Paired {
    const SEGMENTS: Option<usize> = Some(2);
}
impl Layout for Segments {
    const SEGMENTS: Option<usize> = None;
}
impl FixedLayout for Single {
    const K: usize = 1;
}
impl FixedLayout for Paired {
    const K: usize = 2;
}
