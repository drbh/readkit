//! Building batches from pooled columns filled elsewhere, taking them apart
//! again, and converting between layouts generically.

use readkit_batch::{
    Batch, BatchBuilder, BatchParts, ErrorKind, Layout, Paired, Pool, PooledVec, Read,
    SegmentParts, Segments, Single,
};

fn filled<T: readkit_batch::Element>(pool: &Pool, values: &[T]) -> PooledVec<T> {
    let mut v = pool.take::<T>(values.len());
    v.vec_mut().extend_from_slice(values);
    v
}

/// One segment slot: `seqs` back to back, quality `I` for every base, and
/// names when given.
fn segment(pool: &Pool, seqs: &[&[u8]], names: Option<&[&[u8]]>) -> SegmentParts {
    let binary = |values: &[&[u8]]| {
        let mut offsets = vec![0i32];
        let mut bytes = Vec::new();
        for v in values {
            bytes.extend_from_slice(v);
            offsets.push(bytes.len() as i32);
        }
        (filled(pool, &offsets), filled(pool, &bytes))
    };
    let (offsets, seqs) = binary(seqs);
    let quals = filled(pool, &vec![b'I'; seqs.len()]);
    SegmentParts {
        names: names.map(binary),
        offsets,
        seqs,
        quals,
    }
}

fn parts(segments: Vec<SegmentParts>) -> BatchParts {
    BatchParts {
        segments,
        tags: None,
    }
}

#[test]
fn paired_batch_from_columns_filled_on_worker_threads() {
    let pool = Pool::new();
    // Two "workers" each fill one mate's columns on their own thread.
    let [r1, r2] = std::thread::scope(|s| {
        let a = s.spawn(|| segment(&pool, &[b"ACGT", b"", b"GG"], Some(&[b"a", b"b", b"c"])));
        let b = s.spawn(|| segment(&pool, &[b"T", b"CC", b"AAA"], Some(&[b"a", b"b", b"c"])));
        [a.join().unwrap(), b.join().unwrap()]
    });
    let fresh = pool.stats().fresh;
    let mut tags = pool.take::<u16>(3);
    tags.vec_mut().extend([7, 8, 9]);
    let batch = Batch::<Paired>::from_parts(BatchParts {
        segments: vec![r1, r2],
        tags: Some(tags),
    })
    .unwrap();
    assert_eq!(pool.stats().fresh, fresh + 1, "only the tags were taken");
    assert_eq!((batch.len(), batch.segments()), (3, 2));
    assert!(batch.has_names());
    assert_eq!(batch.tags(), Some(&[7, 8, 9][..]));
    let pairs: Vec<_> = batch.pairs().collect();
    assert_eq!(pairs[0].0, Read::named(b"a", b"ACGT", b"IIII"));
    assert_eq!(pairs[1].0.seq(), b"");
    assert_eq!(pairs[2].1, Read::named(b"c", b"AAA", b"III"));
}

#[test]
fn offsets_need_not_start_at_zero_and_bytes_may_run_past_them() {
    let pool = Pool::new();
    let seg = SegmentParts {
        names: None,
        offsets: filled(&pool, &[2, 4, 7]),
        seqs: filled(&pool, b"xxACGTTxxxx"),
        quals: filled(&pool, b"xxIIIIIxx"),
    };
    let batch = Batch::<Single>::from_parts(parts(vec![seg])).unwrap();
    let seqs = batch.column().seqs();
    assert_eq!(seqs.data(), b"ACGTT");
    assert_eq!(batch.reads().nth(1).unwrap(), Read::new(b"GTT", b"III"));
}

#[test]
fn into_parts_and_back_reuses_the_buffers() {
    let pool = Pool::new();
    let mut b = BatchBuilder::<Single>::with_pool(&pool, 2);
    b.push_read(Read::named(b"r", b"ACGT", b"IIII")).unwrap();
    b.set_tag(3);
    let batch = b.finish();
    let before = pool.stats();
    let mut parts = batch.into_parts();
    assert_eq!(pool.stats(), before, "taking apart returns nothing");
    assert_eq!(parts.segments[0].seqs.as_slice(), b"ACGT");

    // Refill the same buffers with a different record.
    let seg = &mut parts.segments[0];
    let (no, nv) = seg.names.as_mut().unwrap();
    no.vec_mut().truncate(1);
    nv.vec_mut().clear();
    nv.vec_mut().extend_from_slice(b"s");
    no.vec_mut().push(1);
    for (v, text) in [(&mut seg.seqs, b"GG"), (&mut seg.quals, b"#I")] {
        v.vec_mut().clear();
        v.vec_mut().extend_from_slice(text);
    }
    seg.offsets.vec_mut().truncate(1);
    seg.offsets.vec_mut().push(2);
    let batch = Batch::<Single>::from_parts(parts).unwrap();
    assert_eq!(
        batch.reads().next().unwrap(),
        Read::named(b"s", b"GG", b"#I")
    );
    assert_eq!(batch.tags(), Some(&[3][..]));
    assert_eq!(pool.stats().fresh, before.fresh);
}

#[test]
fn empty_batch_from_parts() {
    let pool = Pool::new();
    let batch = Batch::<Segments>::from_parts(parts(vec![
        segment(&pool, &[], None),
        segment(&pool, &[], None),
        segment(&pool, &[], None),
    ]))
    .unwrap();
    assert_eq!((batch.len(), batch.segments()), (0, 3));
}

#[test]
fn from_parts_rejects_invalid_columns() {
    let pool = Pool::new();
    let err = |p: BatchParts| Batch::<Segments>::from_parts(p).unwrap_err();
    let seg = |seqs: &[&[u8]]| segment(&pool, seqs, None);

    // no segments, or not the layout's count
    assert_eq!(err(parts(vec![])).kind(), ErrorKind::SegmentCount);
    let e = Batch::<Paired>::from_parts(parts(vec![seg(&[b"A"])])).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::SegmentCount);

    // no offsets at all
    let mut p = parts(vec![seg(&[b"A"])]);
    p.segments[0].offsets.vec_mut().clear();
    assert_eq!(err(p).kind(), ErrorKind::Offsets);

    // segments with different record counts
    let e = err(parts(vec![seg(&[b"A", b"C"]), seg(&[b"A"])]));
    assert_eq!((e.kind(), e.segment()), (ErrorKind::Offsets, Some(1)));

    // decreasing offsets
    let mut p = parts(vec![seg(&[b"AC", b"G", b"T"])]);
    p.segments[0].offsets.vec_mut()[2] = 1;
    let e = err(p);
    assert_eq!((e.kind(), e.record()), (ErrorKind::Offsets, Some(1)));

    // negative first offset
    let mut p = parts(vec![seg(&[b"AC"])]);
    p.segments[0].offsets.vec_mut()[0] = -1;
    assert_eq!(err(p).kind(), ErrorKind::Offsets);

    // offsets past the sequence bytes, or past the quality bytes: a
    // record's sequence and quality would differ in length
    let mut p = parts(vec![seg(&[b"ACGT"])]);
    p.segments[0].seqs.vec_mut().pop();
    assert_eq!(err(p).kind(), ErrorKind::Offsets);
    let mut p = parts(vec![seg(&[b"ACGT"])]);
    p.segments[0].quals.vec_mut().pop();
    let e = err(p);
    assert_eq!((e.kind(), e.segment()), (ErrorKind::Offsets, Some(0)));

    // names in some segments only, either way round
    let named = || segment(&pool, &[b"A"], Some(&[b"n"]));
    let e = err(parts(vec![named(), seg(&[b"A"])]));
    assert_eq!((e.kind(), e.segment()), (ErrorKind::Names, Some(1)));
    let e = err(parts(vec![seg(&[b"A"]), named()]));
    assert_eq!((e.kind(), e.segment()), (ErrorKind::Names, Some(1)));

    // name offsets checked like sequence offsets
    let mut p = parts(vec![named()]);
    p.segments[0].names.as_mut().unwrap().1.vec_mut().clear();
    assert_eq!(err(p).kind(), ErrorKind::Offsets);
    let mut p = parts(vec![named()]);
    p.segments[0].names.as_mut().unwrap().0.vec_mut().push(1);
    assert_eq!(err(p).kind(), ErrorKind::Offsets);

    // one tag per record
    let mut p = parts(vec![seg(&[b"A", b"C"])]);
    p.tags = Some(filled(&pool, &[1u16]));
    assert_eq!(err(p).kind(), ErrorKind::LengthMismatch);

    // the message names what is wrong and where
    let mut p = parts(vec![seg(&[b"A"]), seg(&[b"ACGT"])]);
    p.segments[1].quals.vec_mut().truncate(2);
    assert_eq!(
        err(p).to_string(),
        "invalid offsets in segment 1: sequence offsets end at 4, past 2 bytes"
    );
}

#[test]
fn rejected_parts_return_their_buffers() {
    let pool = Pool::new();
    let mut p = parts(vec![segment(&pool, &[b"AC"], None)]);
    p.segments[0].quals.vec_mut().clear();
    let taken = pool.stats().fresh;
    assert!(Batch::<Single>::from_parts(p).is_err());
    assert_eq!(pool.stats().returned, taken);
}

/// Generic over the target layout, with no match on it.
fn typed<M: Layout>(batch: Batch<Segments>) -> Result<Batch<M>, usize> {
    batch.try_into_layout::<M>().map_err(|e| e.expected())
}

#[test]
fn generic_layout_conversion() {
    let pool = Pool::new();
    let two = || BatchBuilder::with_segments(&pool, 2, 0).finish();
    let paired: Batch<Paired> = typed(two()).unwrap();
    assert_eq!(paired.pairs().count(), 0);
    let segments: Batch<Segments> = typed(two()).unwrap();
    assert_eq!(segments.segments(), 2);
    assert_eq!(typed::<Single>(two()).unwrap_err(), 1);

    // typed layouts convert too: to Segments always, between fixed
    // layouts never (the counts differ)
    let mut b = BatchBuilder::<Single>::new(&pool);
    b.push_read(Read::new(b"A", b"I")).unwrap();
    let single = b.finish();
    let err = single.try_into_layout::<Paired>().unwrap_err();
    assert_eq!(err.expected(), 2);
    let single: Batch<Single> = err.into_batch().try_into_layout().unwrap();
    let back: Batch<Segments> = single.try_into_layout().unwrap();
    assert_eq!(back.spot(0).get(0).seq(), b"A");
}
