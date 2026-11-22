//! The fast path against VDB, spot by spot (feature `ngs`).
//!
//! Runs live in `<workspace>/.data/sra/<run>/<run>` or `READKIT_SRA_DIR`
//! (`cargo bench -p readkit-sra --features ngs --bench decode -- --fetch`
//! downloads them); missing runs are skipped with a message.
//! `READKIT_SRA_FAST_SPOTS=N` limits each run to its first N spots (default:
//! every spot).
#![cfg(feature = "ngs")]

use std::path::{Path, PathBuf};
use std::sync::Arc;

use readkit_io::BatchSource;
use readkit_sra::{
    DecodeOptions, Names, NgsBackend, NgsReader, SegmentSelection, SegmentsDecoder, SpotReader,
};

const RUNS: [&str; 3] = ["SRR12999707", "SRR13638690", "SRR10971381"];

fn run(acc: &str) -> Option<PathBuf> {
    let dir = std::env::var_os("READKIT_SRA_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../.data/sra"));
    let p = dir.join(acc).join(acc);
    if p.exists() {
        Some(p)
    } else {
        eprintln!("skipping {acc}: {} not found", p.display());
        None
    }
}

fn spot_limit() -> u64 {
    std::env::var("READKIT_SRA_FAST_SPOTS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(u64::MAX)
}

/// Compare spots `[start, end)` read both ways; returns (spots, fallbacks).
fn compare(path: &Path, start: u64, end: u64) -> (u64, u64) {
    let mut fast = NgsReader::open(path, true).unwrap();
    assert!(fast.fast_path(), "{}: fast path not used", path.display());
    let mut vdb = NgsReader::open_vdb(path, true).unwrap();
    for spot in start..end {
        let a = fast.read_spot(spot).unwrap();
        let b = vdb.read_spot(spot).unwrap();
        let same = a.bases == b.bases
            && a.quals == b.quals
            && a.starts == b.starts
            && a.lens == b.lens
            && a.types == b.types
            && a.name == b.name;
        if !same {
            panic!(
                "{} spot {spot} differs:\n fast bases {:?}\n  vdb bases {:?}\n fast quals {:?}\n  vdb quals {:?}\n lens {:?} / {:?}",
                path.display(),
                String::from_utf8_lossy(a.bases),
                String::from_utf8_lossy(b.bases),
                String::from_utf8_lossy(a.quals),
                String::from_utf8_lossy(b.quals),
                a.lens,
                b.lens,
            );
        }
    }
    if let Some(r) = fast.fallback_reason() {
        eprintln!(
            "{} [{start}, {end}): through VDB because {r}",
            path.display()
        );
    }
    (end - start, fast.vdb_fallback_spots())
}

/// Every spot of every run, across all cores.
#[test]
fn fast_path_matches_vdb_on_every_spot() {
    for acc in RUNS {
        let Some(path) = run(acc) else { continue };
        let info = NgsReader::open_vdb(&path, false).unwrap().info().clone();
        let end = info
            .end_spot()
            .min(info.first_spot.saturating_add(spot_limit()));
        let threads = std::thread::available_parallelism().map_or(4, |n| n.get()) as u64;
        let per = (end - info.first_spot).div_ceil(threads).max(1);
        let t = std::time::Instant::now();
        let (spots, fallbacks) = std::thread::scope(|s| {
            let handles: Vec<_> = (0..threads)
                .map(|i| {
                    let a = info.first_spot + i * per;
                    let b = (a + per).min(end);
                    let path = &path;
                    s.spawn(move || if a < b { compare(path, a, b) } else { (0, 0) })
                })
                .collect();
            handles
                .into_iter()
                .map(|h| h.join().unwrap())
                .fold((0, 0), |x, y| (x.0 + y.0, x.1 + y.1))
        });
        eprintln!(
            "{acc}: {spots} spots identical ({fallbacks} through VDB) in {:.1}s",
            t.elapsed().as_secs_f64()
        );
        assert_eq!(spots, end - info.first_spot);
        assert_eq!(fallbacks, 0, "{acc}: blobs the fast path did not decode");
    }
}

/// Through the decoder: every read (technical too), names, threads.
#[test]
fn decoder_batches_match_vdb() {
    for acc in RUNS {
        let Some(path) = run(acc) else { continue };
        let opts = DecodeOptions {
            segments: SegmentSelection::All,
            names: Names::FasterqSplitSpot,
            threads: 4,
            batch_spots: 5_000,
            spot_range: Some(1..200_001),
            ..Default::default()
        };
        let mut fast =
            SegmentsDecoder::open_with(Arc::new(NgsBackend::new()), &path, opts.clone()).unwrap();
        let mut vdb =
            SegmentsDecoder::open_with(Arc::new(NgsBackend::vdb_only()), &path, opts).unwrap();
        let mut batches = 0;
        loop {
            match (fast.next_batch().unwrap(), vdb.next_batch().unwrap()) {
                (Some(a), Some(b)) => {
                    assert_eq!(a.len(), b.len(), "{acc} batch {batches}");
                    assert_eq!(a.segments(), b.segments(), "{acc}");
                    for s in 0..a.segments() {
                        let (x, y) = (a.segment(s), b.segment(s));
                        assert_eq!(x.seqs().offsets(), y.seqs().offsets(), "{acc}");
                        assert_eq!(x.seqs().data(), y.seqs().data(), "{acc} bases");
                        assert_eq!(x.quals().data(), y.quals().data(), "{acc} qualities");
                        assert_eq!(
                            x.names().map(|n| n.data()),
                            y.names().map(|n| n.data()),
                            "{acc} names"
                        );
                    }
                    batches += 1;
                }
                (None, None) => break,
                (a, b) => panic!(
                    "{acc}: one decoder ended first ({} {})",
                    a.is_some(),
                    b.is_some()
                ),
            }
        }
        assert_eq!(fast.stats(), vdb.stats(), "{acc}");
        eprintln!("{acc}: {batches} batches identical");
    }
}

/// A run the fast path does not take is decoded through VDB unchanged.
#[test]
fn vdb_only_backend_never_uses_the_fast_path() {
    let Some(path) = RUNS.iter().find_map(|a| run(a)) else {
        return;
    };
    let r = NgsReader::open_vdb(&path, false).unwrap();
    assert!(!r.fast_path());
    assert_eq!(r.vdb_fallback_spots(), 0);
}

/// Readers opened and closed on many threads at once. ncbi-vdb frees a
/// process-global string in every manager teardown without a lock; with a
/// manager per reader this double-freed (shim keeps one shared manager).
#[test]
fn readers_close_concurrently() {
    let Some(path) = RUNS.iter().find_map(|a| run(a)) else {
        return;
    };
    std::thread::scope(|s| {
        for _ in 0..16 {
            s.spawn(|| {
                for i in 0..40 {
                    let mut r = if i % 2 == 0 {
                        NgsReader::open_vdb(&path, false).unwrap()
                    } else {
                        NgsReader::open(&path, false).unwrap()
                    };
                    r.read_spot(1).unwrap();
                }
            });
        }
    });
}
