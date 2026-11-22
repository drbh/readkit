//! Building batches.

use std::fmt;
use std::marker::PhantomData;

use crate::batch::{Batch, Read, SegmentParts};
use crate::error::{Error, ErrorKind};
use crate::layout::{FixedLayout, Layout, Paired, Segments, Single};
use crate::pool::{Pool, PooledVec};

/// Builds a [`Batch`] record by record, drawing buffers from a [`Pool`].
///
/// Every push is validated before anything is written, so a failed push
/// leaves the builder unchanged. Checks: segment count, sequence and quality
/// lengths equal, names present for every segment of every record or for
/// none (decided by the first record), and the `i32` offset limit. Bytes
/// are stored exactly as given.
///
/// ```
/// use readkit_batch::{BatchBuilder, Paired, Pool, Read};
/// let pool = Pool::new();
/// let mut b = BatchBuilder::<Paired>::with_pool(&pool, 1024);
/// b.push_pair(Read::named(b"r1", b"ACGT", b"IIII"), Read::named(b"r1", b"GG", b"II"))?;
/// let batch = b.finish();
/// assert_eq!(batch.len(), 1);
/// assert_eq!(batch.r2().seqs().get(0), b"GG");
/// # Ok::<(), readkit_batch::Error>(())
/// ```
pub struct BatchBuilder<L> {
    pool: Pool,
    hint: usize,
    len: usize,
    names: Option<bool>,
    segs: Vec<SegmentParts>,
    tags: Option<PooledVec<u16>>,
    _layout: PhantomData<fn() -> L>,
}

impl<L: FixedLayout> BatchBuilder<L> {
    /// A builder for about `records` records.
    pub fn with_pool(pool: &Pool, records: usize) -> Self {
        Self::make(pool, L::K, records)
    }

    /// A builder with no capacity hint.
    pub fn new(pool: &Pool) -> Self {
        Self::with_pool(pool, 0)
    }
}

impl BatchBuilder<Segments> {
    /// A builder for records of `k` segments, about `records` of them.
    ///
    /// # Panics
    /// If `k == 0`.
    pub fn with_segments(pool: &Pool, k: usize, records: usize) -> Self {
        assert!(k > 0, "a batch needs at least one segment");
        Self::make(pool, k, records)
    }
}

impl BatchBuilder<Single> {
    /// Appends a single-end read.
    pub fn push_read(&mut self, read: Read<'_>) -> Result<(), Error> {
        self.push(&[read])
    }
}

impl BatchBuilder<Paired> {
    /// Appends a mate pair.
    pub fn push_pair(&mut self, r1: Read<'_>, r2: Read<'_>) -> Result<(), Error> {
        self.push(&[r1, r2])
    }
}

impl<L: Layout> BatchBuilder<L> {
    fn make(pool: &Pool, k: usize, records: usize) -> Self {
        let segs = (0..k)
            .map(|_| {
                let mut offsets = pool.take::<i32>(records + 1);
                offsets.vec_mut().push(0);
                SegmentParts {
                    names: None,
                    offsets,
                    seqs: pool.take(0),
                    quals: pool.take(0),
                }
            })
            .collect();
        Self {
            pool: pool.clone(),
            hint: records,
            len: 0,
            names: None,
            segs,
            tags: None,
            _layout: PhantomData,
        }
    }

    /// Records pushed so far.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether no record has been pushed.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Segments per record.
    pub fn segments(&self) -> usize {
        self.segs.len()
    }

    /// Bytes of names, sequences and qualities so far (for batch sizing).
    pub fn byte_len(&self) -> usize {
        self.segs
            .iter()
            .map(|s| s.seqs.len() + s.quals.len() + s.names.as_ref().map_or(0, |(_, v)| v.len()))
            .sum()
    }

    /// Appends one record given as one [`Read`] per segment.
    pub fn push(&mut self, record: &[Read<'_>]) -> Result<(), Error> {
        let rec = self.len;
        if record.len() != self.segs.len() {
            return Err(self.count_error(record.len()));
        }
        let named = record[0].name().is_some();
        for (s, r) in record.iter().enumerate() {
            if r.seq().len() != r.qual().len() {
                return Err(Error::new(
                    ErrorKind::LengthMismatch,
                    format!("seq {} bytes, qual {} bytes", r.seq().len(), r.qual().len()),
                )
                .at(rec, Some(s)));
            }
            if r.name().is_some() != named {
                return Err(Error::new(ErrorKind::Names, "").at(rec, Some(s)));
            }
        }
        self.prepare(
            named,
            |s| record[s].seq().len(),
            |s| record[s].name().map_or(0, <[u8]>::len),
        )?;
        for (seg, r) in self.segs.iter_mut().zip(record) {
            push_name(seg, r.name());
            seg.seqs.vec_mut().extend_from_slice(r.seq());
            seg.quals.vec_mut().extend_from_slice(r.qual());
            seg.offsets.vec_mut().push(seg.seqs.len() as i32);
        }
        self.end_record();
        Ok(())
    }

    /// Appends one record by writing directly into the columns, for
    /// producers that already hold bases and qualities in separate arrays
    /// (a VDB decoder, say) and should not format text.
    ///
    /// `lens[s]` is segment `s`'s length; `names`, if given, has one name per
    /// segment. `fill(s, seq, qual)` is called once per segment with
    /// zero-initialized slices of exactly `lens[s]` bytes each, so sequence
    /// and quality lengths are equal by construction.
    pub fn push_with<F>(
        &mut self,
        lens: &[usize],
        names: Option<&[&[u8]]>,
        mut fill: F,
    ) -> Result<(), Error>
    where
        F: FnMut(usize, &mut [u8], &mut [u8]),
    {
        if lens.len() != self.segs.len() {
            return Err(self.count_error(lens.len()));
        }
        if let Some(n) = names
            && n.len() != self.segs.len()
        {
            return Err(Error::new(
                ErrorKind::Names,
                format!("{} names for {} segments", n.len(), self.segs.len()),
            )
            .at(self.len, None));
        }
        self.prepare(
            names.is_some(),
            |s| lens[s],
            |s| names.map_or(0, |n| n[s].len()),
        )?;
        for (s, seg) in self.segs.iter_mut().enumerate() {
            append(seg, names.map(|n| n[s]), lens[s], |sq, ql| fill(s, sq, ql));
        }
        self.end_record();
        Ok(())
    }

    /// Sets the tag of the most recently pushed record. The first call adds
    /// a tags column, with tag 0 for earlier records.
    ///
    /// # Panics
    /// If no record has been pushed.
    pub fn set_tag(&mut self, tag: u16) {
        assert!(self.len > 0, "set_tag before any record");
        let (len, hint) = (self.len, self.hint);
        let pool = &self.pool;
        let tags = self.tags.get_or_insert_with(|| {
            let mut t = pool.take::<u16>(hint.max(len));
            t.vec_mut().resize(len, 0);
            t
        });
        tags.vec_mut()[len - 1] = tag;
    }

    /// Finishes the batch. A batch with no records has no names.
    pub fn finish(self) -> Batch<L> {
        Batch::assemble(self.len, self.segs, self.tags)
    }

    fn count_error(&self, got: usize) -> Error {
        Error::new(
            ErrorKind::SegmentCount,
            format!("got {got}, batch has {}", self.segs.len()),
        )
        .at(self.len, None)
    }

    /// Checks names consistency and the offset limit for a record whose
    /// segment `s` has `seq_len(s)` bases and a `name_len(s)`-byte name.
    fn prepare(
        &mut self,
        named: bool,
        seq_len: impl Fn(usize) -> usize,
        name_len: impl Fn(usize) -> usize,
    ) -> Result<(), Error> {
        let rec = self.len;
        if self.names.is_some_and(|n| n != named) {
            return Err(Error::new(
                ErrorKind::Names,
                if named {
                    "named record in a batch without names"
                } else {
                    "unnamed record in a batch with names"
                },
            )
            .at(rec, None));
        }
        for (s, seg) in self.segs.iter().enumerate() {
            let fits = |end: i32, add: usize| {
                (end as usize)
                    .checked_add(add)
                    .is_some_and(|t| t <= i32::MAX as usize)
            };
            if !fits(*seg.offsets.last().unwrap_or(&0), seq_len(s)) {
                return Err(Error::new(ErrorKind::Overflow, "sequence column").at(rec, Some(s)));
            }
            if let Some((o, _)) = &seg.names
                && !fits(*o.last().unwrap_or(&0), name_len(s))
            {
                return Err(Error::new(ErrorKind::Overflow, "name column").at(rec, Some(s)));
            }
        }
        if self.names.is_none() {
            self.names = Some(named);
            if named {
                for seg in &mut self.segs {
                    let mut o = self.pool.take::<i32>(self.hint + 1);
                    o.vec_mut().push(0);
                    seg.names = Some((o, self.pool.take(0)));
                }
            }
        }
        Ok(())
    }

    fn end_record(&mut self) {
        self.len += 1;
        if let Some(t) = &mut self.tags {
            t.vec_mut().push(0);
        }
    }
}

fn push_name(seg: &mut SegmentParts, name: Option<&[u8]>) {
    if let (Some((o, v)), Some(name)) = (&mut seg.names, name) {
        v.vec_mut().extend_from_slice(name);
        o.vec_mut().push(v.len() as i32);
    }
}

/// Appends one segment through `fill`; limits already checked.
fn append(
    seg: &mut SegmentParts,
    name: Option<&[u8]>,
    len: usize,
    fill: impl FnOnce(&mut [u8], &mut [u8]),
) {
    push_name(seg, name);
    let start = seg.seqs.len();
    seg.seqs.vec_mut().resize(start + len, 0);
    seg.quals.vec_mut().resize(start + len, 0);
    fill(
        &mut seg.seqs.vec_mut()[start..],
        &mut seg.quals.vec_mut()[start..],
    );
    seg.offsets.vec_mut().push((start + len) as i32);
}

impl<L: Layout> fmt::Debug for BatchBuilder<L> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("BatchBuilder")
            .field("len", &self.len)
            .field("segments", &self.segs.len())
            .field("names", &self.names)
            .finish()
    }
}
