//! Producer sketches: how a FASTQ parser, an SRA/VDB decoder and a
//! single-cell decoder fill batches.

mod common;

use readkit_batch::{BatchBuilder, ErrorKind, Paired, Pool, Read, Segments, Single};

/// (a) A FASTQ parser pushes mate pairs read in lockstep from two files.
#[test]
fn paired_from_fastq_records() {
    let pool = Pool::new();
    let (r1, r2) = (common::records(common::R1), common::records(common::R2));
    let mut b = BatchBuilder::<Paired>::with_pool(&pool, r1.len());
    for (a, m) in r1.iter().zip(&r2) {
        b.push_pair(*a, *m).unwrap();
    }
    let batch = b.finish();
    assert_eq!((batch.len(), batch.segments()), (3, 2));
    assert!(batch.has_names());
    let pairs: Vec<_> = batch.pairs().collect();
    assert_eq!(pairs[0].0.seq(), b"ACGTNacgtRYKM"); // case and IUPAC kept
    assert_eq!(pairs[0].1.name(), Some(&b"SRR1.1 1/2"[..]));
    assert!(pairs[1].0.is_empty()); // zero-length read
    assert_eq!(pairs[1].1.len(), 12);
    assert_eq!(pairs[2].0.qual(), b"!\"#$%&'");
    assert_eq!(batch.r1().seqs().offsets(), &[0, 13, 13, 20]);
    assert_eq!(batch.r2().seqs().offsets(), &[0, 4, 16, 17]);
}

/// A VDB-like spot: one `READ` array of all bases, numeric phred qualities,
/// and per-segment starts and lengths.
struct VdbSpot {
    read: &'static [u8],
    phred: &'static [u8],
    read_start: [usize; 2],
    read_len: [usize; 2],
}

/// (b) An SRA decoder writes bases and qualities straight into the columns,
/// converting phred to ASCII in place, with no FASTQ text in between.
#[test]
fn columns_filled_from_base_and_quality_arrays() {
    let spots = [
        VdbSpot {
            read: b"ACGTTGCA",
            phred: &[30, 30, 2, 40, 10, 20, 30, 41],
            read_start: [0, 4],
            read_len: [4, 4],
        },
        VdbSpot {
            read: b"NNACG",
            phred: &[0, 0, 35, 35, 35],
            read_start: [0, 2],
            read_len: [2, 3],
        },
    ];
    let pool = Pool::new();
    let mut b = BatchBuilder::<Paired>::with_pool(&pool, spots.len());
    for spot in &spots {
        b.push_with(&spot.read_len, None, |s, seq, qual| {
            let r = spot.read_start[s]..spot.read_start[s] + spot.read_len[s];
            seq.copy_from_slice(&spot.read[r.clone()]);
            for (q, p) in qual.iter_mut().zip(&spot.phred[r]) {
                *q = p + 33;
            }
        })
        .unwrap();
    }
    let batch = b.finish();
    assert!(!batch.has_names());
    assert_eq!(batch.r1().seqs().data(), b"ACGTNN");
    assert_eq!(batch.r2().seqs().data(), b"TGCAACG");
    assert_eq!(batch.r1().quals().get(0), &[63, 63, 35, 73]);
    assert_eq!(batch.r2().quals().get(1), b"DDD");
}

/// (c) A 10x-style spot: barcode+UMI 28 bp, cDNA 90 bp, I1 10 bp, I2 10 bp.
#[test]
fn four_segment_single_cell_spot() {
    let pool = Pool::new();
    let lens = [28usize, 90, 10, 10];
    let mut b = BatchBuilder::<Segments>::with_segments(&pool, 4, 8);
    for i in 0..8u8 {
        // the cDNA read is trimmed on some spots
        let lens = [lens[0], lens[1] - usize::from(i) * 5, lens[2], lens[3]];
        let name = format!("spot{i}");
        let names = [name.as_bytes(); 4];
        b.push_with(&lens, Some(&names), |s, seq, qual| {
            seq.fill(b"ACGT"[s]);
            qual.fill(b'F');
        })
        .unwrap();
    }
    let batch = b.finish();
    assert_eq!(batch.segments(), 4);
    for spot in batch.spots() {
        let l: Vec<usize> = spot.iter().map(|r| r.len()).collect();
        assert_eq!(l, [28, 90 - spot.index() * 5, 10, 10]);
        assert_eq!(
            spot.get(3).name(),
            Some(format!("spot{}", spot.index()).as_bytes())
        );
    }
    assert_eq!(
        batch.segment(1).seqs().data().len(),
        (0..8).map(|i| 90 - i * 5).sum()
    );
    assert!(batch.segment(2).seqs().data().iter().all(|&b| b == b'G'));

    // runtime k converts to a typed layout only when it matches
    let err = batch.try_into_layout::<Paired>().unwrap_err();
    assert_eq!((err.kind(), err.expected()), (ErrorKind::Layout, 2));
    let batch = err.into_batch();
    assert_eq!(batch.len(), 8);
}

#[test]
fn typed_layouts_convert_through_segments() {
    let pool = Pool::new();
    let mut b = BatchBuilder::<Single>::new(&pool);
    b.push_read(Read::new(b"AC", b"II")).unwrap();
    let seg: readkit_batch::Batch<Segments> = b.finish().into();
    assert_eq!(seg.spot(0).get(0).seq(), b"AC");
    let single = seg.try_into_layout::<Single>().unwrap();
    assert_eq!(single.reads().next().unwrap().qual(), b"II");
}

#[test]
fn builder_rejects_invalid_records_without_changing_state() {
    let pool = Pool::new();
    let mut b = BatchBuilder::<Paired>::new(&pool);
    b.push_pair(Read::named(b"a", b"A", b"I"), Read::named(b"a", b"C", b"I"))
        .unwrap();

    let e = b
        .push_pair(
            Read::named(b"b", b"AC", b"I"),
            Read::named(b"b", b"C", b"I"),
        )
        .unwrap_err();
    assert_eq!(e.kind(), ErrorKind::LengthMismatch);
    assert_eq!((e.record(), e.segment()), (Some(1), Some(0)));

    let e = b
        .push_pair(Read::named(b"b", b"A", b"I"), Read::new(b"C", b"I"))
        .unwrap_err();
    assert_eq!((e.kind(), e.segment()), (ErrorKind::Names, Some(1)));

    let e = b
        .push_pair(Read::new(b"A", b"I"), Read::new(b"C", b"I"))
        .unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Names);

    let e = b.push(&[Read::named(b"x", b"A", b"I")]).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::SegmentCount);
    let e = b.push_with(&[1, 1, 1], None, |_, _, _| {}).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::SegmentCount);
    assert!(e.to_string().contains("record 1"), "{e}");

    let batch = b.finish();
    assert_eq!(batch.len(), 1);
    assert_eq!(batch.r1().seqs().offsets(), &[0, 1]);
    assert_eq!(batch.r1().names().unwrap().offsets(), &[0, 1]);
}

#[test]
fn builder_tags() {
    let pool = Pool::new();
    let mut b = BatchBuilder::<Single>::new(&pool);
    for i in 0..4u16 {
        b.push_read(Read::new(b"A", b"I")).unwrap();
        if i % 2 == 1 {
            b.set_tag(i);
        }
    }
    b.push_read(Read::new(b"A", b"I")).unwrap();
    let mut batch = b.finish();
    assert_eq!(batch.tags(), Some(&[0, 1, 0, 3, 0][..]));
    batch.tags_mut().unwrap()[0] = 9;
    assert_eq!(batch.tags().unwrap()[0], 9);
    batch.clear_tags();
    assert_eq!(batch.tags_or_insert(&pool, 7), &[7; 5]);
}

#[test]
fn empty_batch() {
    let pool = Pool::new();
    let batch = BatchBuilder::<Paired>::new(&pool).finish();
    assert!(batch.is_empty());
    assert!(!batch.has_names());
    assert_eq!(batch.r1().seqs().offsets(), &[0]);
    assert_eq!(batch.pairs().count(), 0);
    assert_eq!(batch.byte_len(), 0);
}

#[test]
fn batches_are_send_and_sync() {
    fn check<T: Send + Sync>() {}
    check::<readkit_batch::Batch<Single>>();
    check::<readkit_batch::Batch<Segments>>();
    check::<BatchBuilder<Paired>>();
    check::<Pool>();
}
