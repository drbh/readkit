//! The batch type and its borrowed views.

use std::fmt;
use std::marker::PhantomData;

use crate::error::{Error, ErrorKind};
use crate::layout::{Layout, Paired, Segments, Single};
use crate::pool::{Pool, PooledVec};

/// The pooled columns of one segment slot, for building a batch from
/// columns filled elsewhere ([`Batch::from_parts`]) or taking one apart
/// ([`Batch::into_parts`]).
///
/// All three columns are in Arrow `Binary` layout: `n + 1` offsets into a
/// byte buffer. Sequence and quality share `offsets`, so record `i` has
/// `seqs[offsets[i]..offsets[i + 1]]` and `quals` over the same range, and
/// their lengths are equal by construction. Offsets need not start at 0,
/// and the byte buffers may extend past the last offset.
pub struct SegmentParts {
    /// Name offsets and bytes, or `None` when the batch has no names.
    pub names: Option<(PooledVec<i32>, PooledVec<u8>)>,
    /// Offsets shared by `seqs` and `quals`.
    pub offsets: PooledVec<i32>,
    /// Sequence bytes.
    pub seqs: PooledVec<u8>,
    /// Quality bytes.
    pub quals: PooledVec<u8>,
}

impl fmt::Debug for SegmentParts {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SegmentParts")
            .field("records", &self.offsets.len().saturating_sub(1))
            .field("names", &self.names.is_some())
            .field("seqs", &self.seqs.len())
            .field("quals", &self.quals.len())
            .finish()
    }
}

/// A batch as its pooled buffers: one [`SegmentParts`] per segment slot
/// plus the optional tags column (one `u16` per record).
///
/// See [`Batch::from_parts`] and [`Batch::into_parts`].
#[derive(Debug)]
pub struct BatchParts {
    /// The segment slots' columns, in segment order.
    pub segments: Vec<SegmentParts>,
    /// Per-record tags, if any.
    pub tags: Option<PooledVec<u16>>,
}

/// `n` records of `k` segments each, stored column-wise in pooled buffers.
///
/// `L` is a [`Layout`] marker. Dropping the batch returns its buffers to the
/// pools they came from. A batch is `Send` and is moved between stages,
/// never shared mutably.
pub struct Batch<L> {
    pub(crate) len: usize,
    pub(crate) segs: Vec<SegmentParts>,
    pub(crate) tags: Option<PooledVec<u16>>,
    pub(crate) _layout: PhantomData<fn() -> L>,
}

impl<L: Layout> Batch<L> {
    /// A batch from already validated columns.
    pub(crate) fn assemble(
        len: usize,
        segs: Vec<SegmentParts>,
        tags: Option<PooledVec<u16>>,
    ) -> Self {
        debug_assert!(L::SEGMENTS.is_none_or(|k| k == segs.len()));
        Self {
            len,
            segs,
            tags,
            _layout: PhantomData,
        }
    }

    /// A batch from columns filled elsewhere, for example on worker
    /// threads, without copying them.
    ///
    /// The record count `n` is `offsets.len() - 1` of the first segment.
    /// Checks, each failing with the [`ErrorKind`] shown:
    ///
    /// - at least one segment, and as many as `L` fixes
    ///   ([`SegmentCount`](ErrorKind::SegmentCount));
    /// - every offsets array (sequence and name) has `n + 1` entries, is
    ///   non-negative and non-decreasing, and its last offset lies within
    ///   the bytes it indexes: `seqs` and `quals` for the shared offsets
    ///   ([`Offsets`](ErrorKind::Offsets));
    /// - names in every segment or in none ([`Names`](ErrorKind::Names));
    /// - tags, if present, one per record
    ///   ([`LengthMismatch`](ErrorKind::LengthMismatch)).
    ///
    /// Sequence and quality lengths are equal per record because they share
    /// one offsets array, and `i32` offsets that pass these checks keep every
    /// column within the 2 GiB limit. Every buffer is a [`PooledVec`], so it
    /// came from a pool and returns there when the batch drops; buffers may
    /// come from different pools. On error the parts are dropped, returning
    /// their buffers.
    ///
    /// ```
    /// use readkit_batch::{Batch, BatchParts, Pool, SegmentParts, Single};
    /// let pool = Pool::new();
    /// let mut offsets = pool.take::<i32>(3);
    /// offsets.vec_mut().extend([0, 4, 6]);
    /// let mut seqs = pool.take::<u8>(6);
    /// seqs.vec_mut().extend_from_slice(b"ACGTNN");
    /// let mut quals = pool.take::<u8>(6);
    /// quals.vec_mut().extend_from_slice(b"IIII##");
    /// let segment = SegmentParts { names: None, offsets, seqs, quals };
    /// let batch = Batch::<Single>::from_parts(BatchParts { segments: vec![segment], tags: None })?;
    /// assert_eq!(batch.reads().nth(1).unwrap().seq(), b"NN");
    /// # Ok::<(), readkit_batch::Error>(())
    /// ```
    pub fn from_parts(parts: BatchParts) -> Result<Self, Error> {
        let BatchParts { segments, tags } = parts;
        let k = segments.len();
        if k == 0 || L::SEGMENTS.is_some_and(|want| want != k) {
            let want = L::SEGMENTS.map_or("at least 1".to_string(), |w| w.to_string());
            return Err(Error::new(
                ErrorKind::SegmentCount,
                format!("got {k} segments, layout needs {want}"),
            ));
        }
        let Some(len) = segments[0].offsets.len().checked_sub(1) else {
            return Err(Error::new(ErrorKind::Offsets, "no offsets (need n + 1)").at_segment(0));
        };
        let named = segments[0].names.is_some();
        for (s, seg) in segments.iter().enumerate() {
            let bytes = seg.seqs.len().min(seg.quals.len());
            check_offsets(&seg.offsets, len, bytes, "sequence").map_err(|e| e.at_segment(s))?;
            match &seg.names {
                Some((o, v)) => {
                    check_offsets(o, len, v.len(), "name").map_err(|e| e.at_segment(s))?
                }
                None if named => {
                    return Err(Error::new(ErrorKind::Names, "segment has no names").at_segment(s));
                }
                None => {}
            }
            if seg.names.is_some() && !named {
                return Err(Error::new(ErrorKind::Names, "segment has names").at_segment(s));
            }
        }
        if let Some(t) = &tags
            && t.len() != len
        {
            return Err(Error::new(
                ErrorKind::LengthMismatch,
                format!("{} tags for {len} records", t.len()),
            ));
        }
        Ok(Self::assemble(len, segments, tags))
    }

    /// Takes the batch apart into its pooled buffers, for example to refill
    /// them and build the next batch with [`from_parts`](Self::from_parts).
    /// Buffers keep their contents; clear them with
    /// [`PooledVec::vec_mut`] before refilling.
    pub fn into_parts(self) -> BatchParts {
        BatchParts {
            segments: self.segs,
            tags: self.tags,
        }
    }

    /// Number of records (spots).
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether the batch has no records.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Number of segments per record (`k`).
    pub fn segments(&self) -> usize {
        self.segs.len()
    }

    /// Whether the batch carries read names.
    pub fn has_names(&self) -> bool {
        self.segs.first().is_some_and(|s| s.names.is_some())
    }

    /// The columns of segment slot `i`.
    ///
    /// # Panics
    /// If `i >= self.segments()`.
    pub fn segment(&self, i: usize) -> Column<'_> {
        Column {
            seg: &self.segs[i],
            len: self.len,
        }
    }

    /// Iterates over all segment slots' columns.
    pub fn columns(&self) -> impl ExactSizeIterator<Item = Column<'_>> + '_ {
        self.segs
            .iter()
            .map(move |seg| Column { seg, len: self.len })
    }

    /// The per-record tags, if present. readkit assigns them no meaning.
    pub fn tags(&self) -> Option<&[u16]> {
        self.tags.as_deref()
    }

    /// Mutable per-record tags, if present.
    pub fn tags_mut(&mut self) -> Option<&mut [u16]> {
        self.tags.as_mut().map(|t| t.vec_mut().as_mut_slice())
    }

    /// The per-record tags, adding a column filled with `fill` (taken from
    /// `pool`) if the batch has none.
    pub fn tags_or_insert(&mut self, pool: &Pool, fill: u16) -> &mut [u16] {
        let len = self.len;
        let tags = self.tags.get_or_insert_with(|| {
            let mut t = pool.take::<u16>(len);
            t.vec_mut().resize(len, fill);
            t
        });
        tags.vec_mut().as_mut_slice()
    }

    /// Removes the tags column, returning its buffer to its pool.
    pub fn clear_tags(&mut self) {
        self.tags = None;
    }

    /// Total bytes of names, sequences and qualities.
    pub fn byte_len(&self) -> usize {
        self.columns()
            .map(|c| c.seqs().data().len() * 2 + c.names().map_or(0, |n| n.data().len()))
            .sum()
    }

    /// Erases the layout type, keeping the data.
    pub fn into_segments(self) -> Batch<Segments> {
        self.relabel()
    }

    /// Converts to layout `M` when the segment count matches: any batch
    /// converts to [`Segments`], and to [`Single`] or [`Paired`] when it has
    /// one or two segments. Code generic over `M: Layout` can name its
    /// target layout without matching on it.
    ///
    /// ```
    /// use readkit_batch::{Batch, BatchBuilder, Layout, Paired, Pool, Segments};
    /// fn typed<M: Layout>(batch: Batch<Segments>) -> Option<Batch<M>> {
    ///     batch.try_into_layout::<M>().ok()
    /// }
    /// let batch = BatchBuilder::with_segments(&Pool::new(), 2, 0).finish();
    /// let batch: Batch<Paired> = typed(batch).unwrap();
    /// assert!(typed::<readkit_batch::Single>(batch.into_segments()).is_none());
    /// ```
    pub fn try_into_layout<M: Layout>(self) -> Result<Batch<M>, LayoutError> {
        match M::SEGMENTS {
            Some(k) if k != self.segs.len() => Err(LayoutError {
                batch: self.relabel(),
                expected: k,
            }),
            _ => Ok(self.relabel()),
        }
    }

    fn relabel<M: Layout>(self) -> Batch<M> {
        Batch {
            len: self.len,
            segs: self.segs,
            tags: self.tags,
            _layout: PhantomData,
        }
    }
}

impl Batch<Segments> {
    /// Iterates over records as [`Spot`]s.
    pub fn spots(&self) -> impl ExactSizeIterator<Item = Spot<'_>> + '_ {
        (0..self.len).map(move |index| Spot {
            segs: &self.segs,
            index,
        })
    }

    /// The record at `index` as a [`Spot`].
    ///
    /// # Panics
    /// If `index >= self.len()`.
    pub fn spot(&self, index: usize) -> Spot<'_> {
        assert!(index < self.len, "spot index out of range");
        Spot {
            segs: &self.segs,
            index,
        }
    }
}

impl Batch<Single> {
    /// The only segment slot's columns.
    pub fn column(&self) -> Column<'_> {
        self.segment(0)
    }

    /// Iterates over the reads.
    pub fn reads(&self) -> impl ExactSizeIterator<Item = Read<'_>> + '_ {
        let c = self.column();
        (0..self.len).map(move |i| c.read(i))
    }
}

impl Batch<Paired> {
    /// Columns of the first mate.
    pub fn r1(&self) -> Column<'_> {
        self.segment(0)
    }

    /// Columns of the second mate.
    pub fn r2(&self) -> Column<'_> {
        self.segment(1)
    }

    /// Iterates over mate pairs.
    pub fn pairs(&self) -> impl ExactSizeIterator<Item = (Read<'_>, Read<'_>)> + '_ {
        let (a, b) = (self.r1(), self.r2());
        (0..self.len).map(move |i| (a.read(i), b.read(i)))
    }
}

impl From<Batch<Single>> for Batch<Segments> {
    fn from(b: Batch<Single>) -> Self {
        b.into_segments()
    }
}

impl From<Batch<Paired>> for Batch<Segments> {
    fn from(b: Batch<Paired>) -> Self {
        b.into_segments()
    }
}

impl<L: Layout> fmt::Debug for Batch<L> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Batch")
            .field("len", &self.len)
            .field("segments", &self.segs.len())
            .field("names", &self.has_names())
            .field("tags", &self.tags.is_some())
            .finish()
    }
}

/// A failed layout conversion; gives the batch back.
pub struct LayoutError {
    batch: Batch<Segments>,
    expected: usize,
}

impl LayoutError {
    /// Always [`ErrorKind::Layout`].
    pub fn kind(&self) -> ErrorKind {
        ErrorKind::Layout
    }

    /// The segment count the requested layout needs.
    pub fn expected(&self) -> usize {
        self.expected
    }

    /// The batch that could not be converted.
    pub fn into_batch(self) -> Batch<Segments> {
        self.batch
    }
}

impl fmt::Debug for LayoutError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("LayoutError")
            .field("segments", &self.batch.segments())
            .field("expected", &self.expected)
            .finish()
    }
}

impl fmt::Display for LayoutError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}: batch has {} segments, layout needs {}",
            ErrorKind::Layout,
            self.batch.segments(),
            self.expected
        )
    }
}

impl std::error::Error for LayoutError {}

/// One variable-length binary column in Arrow `Binary` layout: `len + 1`
/// offsets into a byte buffer. Offsets need not start at 0.
#[derive(Clone, Copy)]
pub struct Binary<'a> {
    offsets: &'a [i32],
    values: &'a [u8],
}

impl<'a> Binary<'a> {
    /// Number of values.
    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    /// Whether there are no values.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The `len + 1` offsets.
    pub fn offsets(&self) -> &'a [i32] {
        self.offsets
    }

    /// The whole value buffer the offsets index into.
    pub fn values(&self) -> &'a [u8] {
        self.values
    }

    /// The bytes covered by this column (`offsets[0]..offsets[len]`): all
    /// values back to back.
    pub fn data(&self) -> &'a [u8] {
        let (a, b) = (self.offsets[0], self.offsets[self.len()]);
        &self.values[a as usize..b as usize]
    }

    /// Value `i`.
    ///
    /// # Panics
    /// If `i >= self.len()`.
    pub fn get(&self, i: usize) -> &'a [u8] {
        let (a, b) = (self.offsets[i], self.offsets[i + 1]);
        &self.values[a as usize..b as usize]
    }

    /// Iterates over the values.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = &'a [u8]> + 'a {
        let values = self.values;
        self.offsets
            .windows(2)
            .map(move |w| &values[w[0] as usize..w[1] as usize])
    }
}

impl fmt::Debug for Binary<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Binary")
            .field("len", &self.len())
            .field("bytes", &self.data().len())
            .finish()
    }
}

/// The columns of one segment slot, for bulk consumers.
#[derive(Clone, Copy)]
pub struct Column<'a> {
    seg: &'a SegmentParts,
    len: usize,
}

impl<'a> Column<'a> {
    /// Number of records.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether there are no records.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Read names, if the batch has them.
    pub fn names(&self) -> Option<Binary<'a>> {
        self.seg.names.as_ref().map(|(o, v)| Binary {
            offsets: o.as_slice(),
            values: v.as_slice(),
        })
    }

    /// Sequences.
    pub fn seqs(&self) -> Binary<'a> {
        Binary {
            offsets: self.seg.offsets.as_slice(),
            values: self.seg.seqs.as_slice(),
        }
    }

    /// Qualities. Shares its offsets with [`seqs`](Self::seqs).
    pub fn quals(&self) -> Binary<'a> {
        Binary {
            offsets: self.seg.offsets.as_slice(),
            values: self.seg.quals.as_slice(),
        }
    }

    /// Record `i` of this segment slot.
    ///
    /// # Panics
    /// If `i >= self.len()`.
    pub fn read(&self, i: usize) -> Read<'a> {
        assert!(i < self.len, "record index out of range");
        read_at(self.seg, i)
    }
}

impl fmt::Debug for Column<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Column")
            .field("len", &self.len)
            .field("names", &self.seg.names.is_some())
            .finish()
    }
}

fn read_at(seg: &SegmentParts, i: usize) -> Read<'_> {
    let span = |o: &[i32]| o[i] as usize..o[i + 1] as usize;
    let r = span(&seg.offsets);
    Read {
        name: seg.names.as_ref().map(|(o, v)| &v[span(o)]),
        seq: &seg.seqs[r.clone()],
        qual: &seg.quals[r],
    }
}

/// One segment of one record: optional name, sequence and quality bytes.
///
/// Also the input type of [`BatchBuilder`](crate::BatchBuilder) pushes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Read<'a> {
    name: Option<&'a [u8]>,
    seq: &'a [u8],
    qual: &'a [u8],
}

impl<'a> Read<'a> {
    /// A read without a name.
    pub fn new(seq: &'a [u8], qual: &'a [u8]) -> Self {
        Self {
            name: None,
            seq,
            qual,
        }
    }

    /// A read with a name.
    pub fn named(name: &'a [u8], seq: &'a [u8], qual: &'a [u8]) -> Self {
        Self {
            name: Some(name),
            seq,
            qual,
        }
    }

    /// The name, without FASTQ's `@`.
    pub fn name(&self) -> Option<&'a [u8]> {
        self.name
    }

    /// The bases, exactly as given.
    pub fn seq(&self) -> &'a [u8] {
        self.seq
    }

    /// The quality bytes, exactly as given.
    pub fn qual(&self) -> &'a [u8] {
        self.qual
    }

    /// Sequence length.
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    /// Whether the sequence is empty.
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
}

/// One record of a [`Batch<Segments>`]: `k` segments.
#[derive(Clone, Copy)]
pub struct Spot<'a> {
    segs: &'a [SegmentParts],
    index: usize,
}

impl<'a> Spot<'a> {
    /// The record number within the batch.
    pub fn index(&self) -> usize {
        self.index
    }

    /// Number of segments.
    pub fn segments(&self) -> usize {
        self.segs.len()
    }

    /// Segment `i`.
    ///
    /// # Panics
    /// If `i >= self.segments()`.
    pub fn get(&self, i: usize) -> Read<'a> {
        read_at(&self.segs[i], self.index)
    }

    /// Iterates over the segments.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = Read<'a>> + 'a {
        let (segs, index) = (self.segs, self.index);
        segs.iter().map(move |seg| read_at(seg, index))
    }
}

impl fmt::Debug for Spot<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

/// Checks one offsets array for `len` records over `bytes` bytes.
fn check_offsets(offsets: &[i32], len: usize, bytes: usize, what: &str) -> Result<(), Error> {
    if offsets.len() != len + 1 {
        return Err(Error::new(
            ErrorKind::Offsets,
            format!(
                "{what} offsets have {} entries for {len} records",
                offsets.len()
            ),
        ));
    }
    if offsets[0] < 0 {
        return Err(Error::new(
            ErrorKind::Offsets,
            format!("{what} offsets start at {}", offsets[0]),
        ));
    }
    if let Some(i) = offsets.windows(2).position(|w| w[1] < w[0]) {
        return Err(Error::new(
            ErrorKind::Offsets,
            format!(
                "{what} offsets decrease ({} to {})",
                offsets[i],
                offsets[i + 1]
            ),
        )
        .at(i, None));
    }
    let end = offsets[len] as usize;
    if end > bytes {
        return Err(Error::new(
            ErrorKind::Offsets,
            format!("{what} offsets end at {end}, past {bytes} bytes"),
        ));
    }
    Ok(())
}
