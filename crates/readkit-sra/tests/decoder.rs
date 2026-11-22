//! Decoder logic against an in-memory backend: segment selection, read
//! types, names, Lite policy, validation, threads and pools. Offline.

use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use readkit_batch::{Batch, Paired, Pool, Segments, Single};
use readkit_io::BatchSource;
use readkit_sra::{
    DecodeOptions, Decoder, Error, ErrorKind, LitePolicy, Names, OpenRequest, Orientation,
    QualityKind, ReadType, RunInfo, SegmentSelection, Spot, SpotReader, SraBackend,
};

#[derive(Clone)]
struct MockSpot {
    bases: Vec<u8>,
    quals: Vec<u8>,
    starts: Vec<i32>,
    lens: Vec<u32>,
    types: Vec<ReadType>,
    name: Option<Vec<u8>>,
}

/// A spot from `(seq, type bits)` reads; qualities derived from the spot id.
fn spot(id: u64, reads: &[(&str, u8)]) -> MockSpot {
    let mut s = MockSpot {
        bases: Vec::new(),
        quals: Vec::new(),
        starts: Vec::new(),
        lens: Vec::new(),
        types: Vec::new(),
        name: None,
    };
    for (seq, t) in reads {
        s.starts.push(s.bases.len() as i32);
        s.lens.push(seq.len() as u32);
        s.types.push(ReadType::from_bits(*t));
        s.bases.extend_from_slice(seq.as_bytes());
        s.quals
            .extend((0..seq.len()).map(|i| b'#' + ((id as usize + i) % 40) as u8));
    }
    s
}

struct Mock {
    spots: Arc<Vec<MockSpot>>,
    quality: QualityKind,
    random: bool,
    opens: Arc<AtomicUsize>,
}

impl Mock {
    fn new(spots: Vec<MockSpot>) -> Self {
        Self {
            spots: Arc::new(spots),
            quality: QualityKind::Original,
            random: true,
            opens: Arc::new(AtomicUsize::new(0)),
        }
    }
}

struct MockReader {
    spots: Arc<Vec<MockSpot>>,
    info: RunInfo,
    random: bool,
    names: bool,
    last: Option<u64>,
}

impl SraBackend for Mock {
    fn name(&self) -> &'static str {
        "mock"
    }
    fn open(&self, _: &Path, request: &OpenRequest) -> Result<Box<dyn SpotReader>, Error> {
        self.opens.fetch_add(1, Ordering::SeqCst);
        let mut info = RunInfo::new(1, self.spots.len() as u64);
        info.quality = self.quality;
        info.has_names = true;
        Ok(Box::new(MockReader {
            spots: self.spots.clone(),
            info,
            random: self.random,
            names: request.names,
            last: None,
        }))
    }
}

impl SpotReader for MockReader {
    fn info(&self) -> &RunInfo {
        &self.info
    }
    fn random_access(&self) -> bool {
        self.random
    }
    fn read_spot(&mut self, id: u64) -> Result<Spot<'_>, Error> {
        if !self.random {
            // A streaming backend: the first spot again (the decoder peeks
            // it at open), then strictly increasing ids.
            assert!(self.last.is_none_or(|l| id == l + 1 || id == 1 && l == 1));
        }
        self.last = Some(id);
        let s = self
            .spots
            .get((id - 1) as usize)
            .ok_or_else(|| Error::new(ErrorKind::NotFound, "no such spot").at_spot(id))?;
        Ok(Spot {
            bases: &s.bases,
            quals: &s.quals,
            starts: &s.starts,
            lens: &s.lens,
            types: &s.types,
            name: if self.names { s.name.as_deref() } else { None },
        })
    }
}

fn path() -> std::path::PathBuf {
    // Any existing path: the decoder checks the file exists.
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml")
}

fn open<L: readkit_sra::DecodeLayout>(
    mock: Mock,
    opts: DecodeOptions,
) -> Result<Decoder<L>, Error> {
    Decoder::<L>::open_with(Arc::new(mock), path(), opts)
}

fn drain<L: readkit_sra::DecodeLayout>(d: &mut Decoder<L>) -> Result<Vec<Batch<L>>, Error> {
    let mut out = Vec::new();
    while let Some(b) = d.next_batch()? {
        out.push(b);
    }
    Ok(out)
}

/// Single-cell-like spots: barcode (technical, forward), cDNA (biological),
/// index (technical, reverse); lengths vary with the spot id.
fn sc_spots(n: u64) -> Vec<MockSpot> {
    (1..=n)
        .map(|i| {
            let cdna = "ACGTN".repeat(1 + (i % 7) as usize);
            spot(
                i,
                &[("TTTTGGGGCCCC", 0b010), (&cdna, 0b001), ("AAC", 0b100)],
            )
        })
        .collect()
}

#[test]
fn biological_selection_skips_technical_reads() {
    let mut d = open::<Single>(Mock::new(sc_spots(50)), DecodeOptions::default()).unwrap();
    assert_eq!(d.segments(), 1);
    assert_eq!(d.read_types(), &[ReadType::from_bits(1)]);
    let batches = drain(&mut d).unwrap();
    let reads: Vec<Vec<u8>> = batches
        .iter()
        .flat_map(|b| b.reads().map(|r| r.seq().to_vec()).collect::<Vec<_>>())
        .collect();
    assert_eq!(reads.len(), 50);
    for (i, r) in reads.iter().enumerate() {
        assert_eq!(r, &"ACGTN".repeat(1 + ((i + 1) % 7)).into_bytes());
    }
    assert_eq!(d.stats().spots, 50);
}

#[test]
fn all_segments_keep_technical_reads_and_orientation() {
    let mut d = open::<Segments>(
        Mock::new(sc_spots(20)),
        DecodeOptions {
            segments: SegmentSelection::All,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(d.segments(), 3);
    let t = d.read_types();
    assert!(t[0].is_technical() && t[0].orientation() == Orientation::Forward);
    assert!(t[1].is_biological());
    assert!(t[2].is_technical() && t[2].orientation() == Orientation::Reverse);
    let batches = drain(&mut d).unwrap();
    let spot = batches[0].spot(3);
    assert_eq!(spot.get(0).seq(), b"TTTTGGGGCCCC");
    assert_eq!(spot.get(2).seq(), b"AAC");
    let expected = spot_ref(4);
    assert_eq!(
        spot.get(1).qual(),
        &expected.quals[12..expected.quals.len() - 3]
    );
}

fn spot_ref(id: u64) -> MockSpot {
    sc_spots(id).pop().unwrap()
}

#[test]
fn selected_reads_in_given_order() {
    let mut d = open::<Paired>(
        Mock::new(sc_spots(10)),
        DecodeOptions {
            segments: SegmentSelection::Reads(vec![1, 0]),
            ..Default::default()
        },
    )
    .unwrap();
    let b = drain(&mut d).unwrap();
    let (r1, r2) = b[0].pairs().next().unwrap();
    assert_eq!(r1.seq(), "ACGTN".repeat(2).as_bytes());
    assert_eq!(r2.seq(), b"TTTTGGGGCCCC");

    let e = open::<Single>(
        Mock::new(sc_spots(3)),
        DecodeOptions {
            segments: SegmentSelection::Reads(vec![5]),
            ..Default::default()
        },
    )
    .unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Layout);
}

#[test]
fn layout_must_fit() {
    let e = open::<Paired>(Mock::new(sc_spots(3)), DecodeOptions::default()).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Layout);
    let only_technical = vec![spot(1, &[("AC", 0b010)])];
    let e = open::<Segments>(Mock::new(only_technical), DecodeOptions::default()).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Layout);
}

#[test]
fn changing_read_structure_is_a_layout_error_at_its_spot() {
    let mut spots = sc_spots(30);
    spots[19] = spot(20, &[("AC", 1), ("GT", 1)]);
    let mut d = open::<Single>(
        Mock::new(spots),
        DecodeOptions {
            batch_spots: 8,
            threads: 1,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(d.next_batch().unwrap().unwrap().len(), 8);
    assert_eq!(d.next_batch().unwrap().unwrap().len(), 8);
    let e = d.next_batch().unwrap_err();
    assert_eq!((e.kind(), e.spot()), (ErrorKind::Layout, Some(20)));
    assert!(d.next_batch().unwrap().is_none(), "nothing after an error");
}

#[test]
fn corrupt_spots_are_rejected_not_misread() {
    let mut short_qual = spot(2, &[("ACGT", 1)]);
    short_qual.quals.pop();
    let mut out_of_range = spot(2, &[("ACGT", 1)]);
    out_of_range.lens[0] = 5;
    let mut negative = spot(2, &[("ACGT", 1)]);
    negative.starts[0] = -1;
    let mut counts = spot(2, &[("ACGT", 1)]);
    counts.starts.push(0);
    for bad in [short_qual, out_of_range, negative, counts] {
        let spots = vec![spot(1, &[("ACGT", 1)]), bad, spot(3, &[("ACGT", 1)])];
        let mut d = open::<Single>(
            Mock::new(spots),
            DecodeOptions {
                threads: 1,
                ..Default::default()
            },
        )
        .unwrap();
        let e = d.next_batch().unwrap_err();
        assert_eq!((e.kind(), e.spot()), (ErrorKind::Corrupt, Some(2)), "{e}");
    }
}

#[test]
fn names_follow_fasterq_dump() {
    let mut spots = sc_spots(3);
    spots[1].name = Some(b"HWI-1:2:3".to_vec());
    let named = |names| {
        let mut d = open::<Segments>(
            Mock::new(spots.clone()),
            DecodeOptions {
                names,
                segments: SegmentSelection::All,
                accession: Some("SRR9".into()),
                ..Default::default()
            },
        )
        .unwrap();
        let b = drain(&mut d).unwrap().remove(0);
        assert!(b.has_names());
        (0..3)
            .map(|i| {
                let s = b.spot(i);
                (0..3)
                    .map(|j| String::from_utf8(s.get(j).name().unwrap().to_vec()).unwrap())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    };
    let n = named(Names::Fasterq);
    assert_eq!(
        n[0],
        [
            "SRR9.1 1 length=12",
            "SRR9.1 1 length=10",
            "SRR9.1 1 length=3"
        ]
    );
    assert_eq!(n[1][1], "SRR9.2 HWI-1:2:3 length=15");
    let n = named(Names::FasterqSplitSpot);
    assert_eq!(
        n[2],
        [
            "SRR9.3/1 3 length=12",
            "SRR9.3/2 3 length=20",
            "SRR9.3/3 3 length=3"
        ]
    );

    let mut d = open::<Single>(Mock::new(sc_spots(3)), DecodeOptions::default()).unwrap();
    assert!(!d.next_batch().unwrap().unwrap().has_names());
}

#[test]
fn declared_lite_is_a_typed_error_or_a_flag() {
    let mut mock = Mock::new(sc_spots(5));
    mock.quality = QualityKind::Synthetic;
    let e = open::<Single>(mock, DecodeOptions::default()).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Lite);
    assert!(e.is_lite());

    let mut mock = Mock::new(sc_spots(5));
    mock.quality = QualityKind::Synthetic;
    let mut d = open::<Single>(
        mock,
        DecodeOptions {
            lite: LitePolicy::Allow,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(d.quality(), QualityKind::Synthetic);
    assert_eq!(drain(&mut d).unwrap().len(), 1);
}

fn lite_spots(n: u64) -> Vec<MockSpot> {
    (1..=n)
        .map(|i| {
            let mut s = spot(i, &[("ACGTACGT", 1)]);
            let q = if i % 10 == 0 { b'$' } else { b'?' };
            s.quals.iter_mut().for_each(|b| *b = q);
            s
        })
        .collect()
}

#[test]
fn undeclared_lite_is_detected_by_sampling() {
    let mut mock = Mock::new(lite_spots(1000));
    mock.quality = QualityKind::Unknown;
    let mut d = open::<Single>(
        mock,
        DecodeOptions {
            batch_spots: 100,
            ..Default::default()
        },
    )
    .unwrap();
    let e = drain(&mut d).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Lite);
    assert_eq!(d.quality(), QualityKind::Synthetic);

    let mut mock = Mock::new(sc_spots(1000));
    mock.quality = QualityKind::Unknown;
    let mut d = open::<Single>(mock, DecodeOptions::default()).unwrap();
    drain(&mut d).unwrap();
    assert_eq!(d.quality(), QualityKind::Original);

    // Lite errors survive wrapping by a downstream stage.
    let mut mock = Mock::new(lite_spots(10));
    mock.quality = QualityKind::Unknown;
    let mut d = open::<Single>(mock, DecodeOptions::default()).unwrap();
    let wrapped = readkit_io::Error::upstream(drain(&mut d).unwrap_err());
    assert!(wrapped.find_cause::<Error>().unwrap().is_lite());
}

fn flatten(batches: &[Batch<Single>]) -> Vec<(Vec<u8>, Vec<u8>)> {
    batches
        .iter()
        .flat_map(|b| {
            b.reads()
                .map(|r| (r.seq().to_vec(), r.qual().to_vec()))
                .collect::<Vec<_>>()
        })
        .collect()
}

#[test]
fn threads_preserve_spot_order() {
    let mut seq = open::<Single>(
        Mock::new(sc_spots(1000)),
        DecodeOptions {
            threads: 1,
            batch_spots: 7,
            ..Default::default()
        },
    )
    .unwrap();
    let a = drain(&mut seq).unwrap();
    let mock = Mock::new(sc_spots(1000));
    let opens = mock.opens.clone();
    let mut par = open::<Single>(
        mock,
        DecodeOptions {
            threads: 4,
            batch_spots: 7,
            ..Default::default()
        },
    )
    .unwrap();
    let b = drain(&mut par).unwrap();
    assert_eq!(a.len(), 143);
    assert_eq!(b.len(), 143);
    assert_eq!(flatten(&a), flatten(&b));
    assert_eq!(opens.load(Ordering::SeqCst), 4, "one reader per worker");
    assert_eq!(par.stats().spots, 1000);
    assert_eq!(par.stats().batches, 143);
}

#[test]
fn streaming_backends_are_read_sequentially() {
    let mut mock = Mock::new(sc_spots(100));
    mock.random = false;
    let opens = mock.opens.clone();
    let mut d = open::<Single>(
        mock,
        DecodeOptions {
            threads: 8,
            batch_spots: 9,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(drain(&mut d).unwrap().len(), 12);
    assert_eq!(opens.load(Ordering::SeqCst), 1);
}

#[test]
fn worker_errors_arrive_in_order() {
    let mut spots = sc_spots(500);
    spots[300].quals.pop();
    let mut d = open::<Single>(
        Mock::new(spots),
        DecodeOptions {
            threads: 4,
            batch_spots: 50,
            ..Default::default()
        },
    )
    .unwrap();
    for _ in 0..6 {
        assert_eq!(d.next_batch().unwrap().unwrap().len(), 50);
    }
    let e = d.next_batch().unwrap_err();
    assert_eq!((e.kind(), e.spot()), (ErrorKind::Corrupt, Some(301)));
    assert!(d.next_batch().unwrap().is_none());
}

#[test]
fn dropping_mid_stream_stops_workers() {
    let mut d = open::<Single>(
        Mock::new(sc_spots(10_000)),
        DecodeOptions {
            threads: 4,
            batch_spots: 10,
            ..Default::default()
        },
    )
    .unwrap();
    d.next_batch().unwrap().unwrap();
    drop(d); // joins; must not hang
}

#[test]
fn pool_buffers_are_reused() {
    let pool = Pool::new();
    let mut d = open::<Single>(
        Mock::new(sc_spots(2000)),
        DecodeOptions {
            threads: 1,
            batch_spots: 100,
            ..Default::default()
        },
    )
    .unwrap()
    .with_pool(&pool);
    d.next_batch().unwrap().unwrap();
    let after_first = pool.stats().fresh;
    while let Some(b) = d.next_batch().unwrap() {
        drop(b);
    }
    assert_eq!(
        pool.stats().fresh,
        after_first,
        "steady state allocates nothing"
    );
    assert!(pool.stats().reused > 0);
}

#[test]
fn missing_files_are_not_found() {
    let e = Decoder::<Single>::open_with(
        Arc::new(Mock::new(Vec::new())),
        "/nonexistent/SRR1",
        DecodeOptions::default(),
    )
    .unwrap_err();
    assert_eq!(e.kind(), ErrorKind::NotFound);
    assert_eq!(e.input(), Some("/nonexistent/SRR1"));
}

#[test]
fn spot_ranges_are_clamped_and_named_by_spot_id() {
    for threads in [1, 3] {
        let mut d = open::<Single>(
            Mock::new(sc_spots(100)),
            DecodeOptions {
                spot_range: Some(41..1000),
                names: Names::Fasterq,
                accession: Some("R".into()),
                batch_spots: 16,
                threads,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(d.spot_range(), 41..101);
        let batches = drain(&mut d).unwrap();
        let first = batches[0].reads().next().unwrap();
        assert_eq!(first.name().unwrap(), b"R.41 41 length=35");
        assert_eq!(batches.iter().map(Batch::len).sum::<usize>(), 60);
    }
    let mut d = open::<Single>(
        Mock::new(sc_spots(10)),
        DecodeOptions {
            spot_range: Some(50..60),
            ..Default::default()
        },
    )
    .unwrap();
    assert!(d.next_batch().unwrap().is_none());
}

#[test]
fn empty_runs_end_at_once() {
    let mut d = open::<Single>(Mock::new(Vec::new()), DecodeOptions::default()).unwrap();
    assert!(d.next_batch().unwrap().is_none());
}

#[cfg(not(feature = "ngs"))]
#[test]
fn without_a_backend_open_is_unsupported() {
    let e = Decoder::<Single>::open(path(), DecodeOptions::default()).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Unsupported);
}
