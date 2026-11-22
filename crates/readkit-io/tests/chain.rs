//! A decoder → stage → sink chain using read-ahead: order, buffer reuse,
//! cancellation by drop, and error chaining.

use std::fmt;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::mpsc;
use std::time::{Duration, Instant};

use readkit_batch::{Batch, BatchBuilder, Pool, Read, Single};
use readkit_io::{BatchSink, BatchSource, Error, ErrorKind, ReadAhead, VecSink, pump, sink_fn};

const PER_BATCH: usize = 16;

/// A decoder-specific error with its own kind, as readkit-sra would have.
#[derive(Debug)]
struct DecodeError {
    kind: DecodeErrorKind,
}

#[derive(Debug, PartialEq, Eq)]
enum DecodeErrorKind {
    Lite,
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "decode failed: {:?}", self.kind)
    }
}
impl std::error::Error for DecodeError {}

/// A synthetic decoder with its own pool.
struct Decoder {
    pool: Pool,
    next: usize,
    batches: Option<usize>, // None = endless
    fail_at: Option<usize>,
    panic_at: Option<usize>,
    delay: Duration,
    reads: Arc<AtomicUsize>,
    dropped: Arc<AtomicBool>,
}

impl Decoder {
    fn new(pool: &Pool, batches: Option<usize>) -> Self {
        Self {
            pool: pool.clone(),
            next: 0,
            batches,
            fail_at: None,
            panic_at: None,
            delay: Duration::ZERO,
            reads: Arc::default(),
            dropped: Arc::default(),
        }
    }
}

impl Drop for Decoder {
    fn drop(&mut self) {
        self.dropped.store(true, Ordering::SeqCst);
    }
}

impl BatchSource<Single> for Decoder {
    type Error = DecodeError;
    fn next_batch(&mut self) -> Result<Option<Batch<Single>>, DecodeError> {
        let i = self.next;
        if self.batches.is_some_and(|n| i >= n) {
            return Ok(None);
        }
        if self.fail_at == Some(i) {
            return Err(DecodeError {
                kind: DecodeErrorKind::Lite,
            });
        }
        assert_ne!(self.panic_at, Some(i), "decoder exploded");
        std::thread::sleep(self.delay);
        self.reads.fetch_add(1, Ordering::SeqCst);
        self.next += 1;
        let mut b = BatchBuilder::<Single>::with_pool(&self.pool, PER_BATCH);
        for r in 0..PER_BATCH {
            let name = (i * PER_BATCH + r).to_string();
            b.push_read(Read::named(name.as_bytes(), b"ACGTNacgt", b"IIIIIIIII"))
                .unwrap();
        }
        Ok(Some(b.finish()))
    }
}

/// A tool stage: passes reads through unchanged and re-tags them, reading
/// ahead of its input. Its tag buffers come from its own pool.
struct Retag<S: BatchSource<Single>> {
    input: ReadAhead<Single, S>,
    pool: Pool,
}

impl<S: BatchSource<Single> + Send + 'static> Retag<S> {
    fn new(input: S, read_ahead: bool) -> Self {
        Self {
            input: ReadAhead::with_read_ahead(input, read_ahead),
            pool: Pool::new(),
        }
    }
}

impl<S: BatchSource<Single>> BatchSource<Single> for Retag<S> {
    type Error = Error;
    fn next_batch(&mut self) -> Result<Option<Batch<Single>>, Error> {
        let Some(mut batch) = self.input.next_batch().map_err(Error::upstream)? else {
            return Ok(None);
        };
        let seqs: Vec<u16> = batch
            .column()
            .seqs()
            .iter()
            .map(|s| s.len() as u16)
            .collect();
        let tags = batch.tags_or_insert(&self.pool, 0);
        for (t, len) in tags.iter_mut().zip(seqs) {
            *t = len;
        }
        Ok(Some(batch))
    }
}

fn names(batch: &Batch<Single>) -> Vec<usize> {
    let names = batch.column().names().unwrap();
    names
        .iter()
        .map(|n| std::str::from_utf8(n).unwrap().parse().unwrap())
        .collect()
}

fn run_ordered(read_ahead: bool) {
    let pool = Pool::new();
    let decoder = Decoder::new(&pool, Some(200));
    let stage = Retag::new(decoder, read_ahead);
    assert_eq!(stage.input.is_threaded(), read_ahead);
    let stage_pool = stage.pool.clone();
    let mut expected = 0usize;
    let sink = sink_fn(|batch: Batch<Single>| {
        for n in names(&batch) {
            assert_eq!(n, expected);
            expected += 1;
        }
        assert_eq!(batch.tags().unwrap(), &[9; PER_BATCH]);
        assert_eq!(batch.column().seqs().get(0), b"ACGTNacgt");
        Ok::<_, std::convert::Infallible>(())
    });
    assert_eq!(pump(stage, sink).unwrap(), 200 * PER_BATCH as u64);

    // Each decoder batch takes 5 buffers. In flight at once: one being
    // built, one waiting in the read-ahead hand-off, one in the stage/sink.
    let s = pool.stats();
    assert!(
        s.fresh <= 5 * 4,
        "decoder pool allocated {} buffers",
        s.fresh
    );
    assert_eq!(s.fresh + s.reused, 200 * 5);
    assert_eq!(s.idle as u64, s.fresh); // all back after the chain ends
    let t = stage_pool.stats();
    assert!(t.fresh <= 2, "stage pool allocated {} buffers", t.fresh);
    assert_eq!(t.fresh + t.reused, 200);
}

#[test]
fn chain_preserves_order_and_recycles_buffers() {
    run_ordered(true);
}

#[test]
fn chain_without_read_ahead() {
    run_ordered(false);
}

#[test]
fn reads_at_most_one_batch_ahead() {
    let pool = Pool::new();
    let decoder = Decoder::new(&pool, None);
    let reads = Arc::clone(&decoder.reads);
    let mut ahead = ReadAhead::new(decoder);
    let first = ahead.next_batch().unwrap().unwrap();
    std::thread::sleep(Duration::from_millis(100));
    // one delivered, one finished and waiting; no third read
    assert_eq!(reads.load(Ordering::SeqCst), 2);
    drop(first);
    let _second = ahead.next_batch().unwrap().unwrap();
    std::thread::sleep(Duration::from_millis(100));
    assert_eq!(reads.load(Ordering::SeqCst), 3);
}

fn within<T: Send + 'static>(limit: Duration, f: impl FnOnce() -> T + Send + 'static) -> T {
    let (tx, rx) = mpsc::channel();
    std::thread::spawn(move || {
        let _ = tx.send(f());
    });
    rx.recv_timeout(limit).expect("chain hung")
}

#[test]
fn dropping_mid_chain_does_not_hang() {
    let pool = Pool::new();
    let decoder = Decoder::new(&pool, None); // endless
    let dropped = Arc::clone(&decoder.dropped);
    within(Duration::from_secs(10), move || {
        let mut stage = Retag::new(decoder, true);
        for _ in 0..3 {
            stage.next_batch().unwrap().unwrap();
        }
        drop(stage);
    });
    // the read-ahead thread was joined and dropped the decoder
    assert!(dropped.load(Ordering::SeqCst));
    // everything that was in flight went back to the pool
    let s = pool.stats();
    assert_eq!(s.idle as u64, s.fresh);
}

#[test]
fn drop_waits_only_for_the_read_in_progress() {
    let pool = Pool::new();
    let mut decoder = Decoder::new(&pool, None);
    decoder.delay = Duration::from_millis(200);
    let dropped = Arc::clone(&decoder.dropped);
    let elapsed = within(Duration::from_secs(10), move || {
        let mut ahead = ReadAhead::new(decoder);
        ahead.next_batch().unwrap().unwrap();
        let t = Instant::now();
        drop(ahead); // the thread is 0..200 ms into the next read
        t.elapsed()
    });
    assert!(dropped.load(Ordering::SeqCst));
    assert!(
        elapsed < Duration::from_millis(1000),
        "drop took {elapsed:?}"
    );
}

#[test]
fn source_error_surfaces_through_the_stage() {
    let pool = Pool::new();
    let mut decoder = Decoder::new(&pool, None);
    decoder.fail_at = Some(3);
    let stage = Retag::new(decoder, true);
    let mut sink = VecSink::new();
    let err = pump(stage, &mut sink).unwrap_err();
    assert_eq!(sink.batches().len(), 3);

    // pump → Retag's error → DecodeError
    assert_eq!(err.kind(), ErrorKind::Upstream);
    let stage_err = std::error::Error::source(&err)
        .unwrap()
        .downcast_ref::<Error>()
        .unwrap();
    assert_eq!(stage_err.kind(), ErrorKind::Upstream);
    let root = std::error::Error::source(stage_err)
        .unwrap()
        .downcast_ref::<DecodeError>()
        .unwrap();
    assert_eq!(root.kind, DecodeErrorKind::Lite);
    // or without knowing the depth
    assert_eq!(
        err.find_cause::<DecodeError>().unwrap().kind,
        DecodeErrorKind::Lite
    );
}

#[test]
fn read_ahead_is_fused_after_an_error() {
    let pool = Pool::new();
    for threaded in [true, false] {
        let mut decoder = Decoder::new(&pool, None);
        decoder.fail_at = Some(1);
        let mut ahead = ReadAhead::with_read_ahead(decoder, threaded);
        assert!(ahead.next_batch().unwrap().is_some());
        assert!(ahead.next_batch().is_err());
        assert!(ahead.next_batch().unwrap().is_none());
    }
}

#[test]
#[should_panic(expected = "decoder exploded")]
fn source_panic_resumes_on_the_consumer() {
    let pool = Pool::new();
    let mut decoder = Decoder::new(&pool, None);
    decoder.panic_at = Some(1);
    let mut ahead = ReadAhead::new(decoder);
    while ahead.next_batch().unwrap().is_some() {}
}

#[test]
fn sink_error_is_reported_as_sink() {
    #[derive(Debug)]
    struct Full;
    impl fmt::Display for Full {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            f.write_str("disk full")
        }
    }
    impl std::error::Error for Full {}

    struct Failing;
    impl BatchSink<Single> for Failing {
        type Error = Full;
        fn accept(&mut self, _: Batch<Single>) -> Result<(), Full> {
            Err(Full)
        }
        fn finish(self) -> Result<(), Full> {
            Ok(())
        }
    }
    let pool = Pool::new();
    let err = pump(Retag::new(Decoder::new(&pool, Some(5)), true), Failing).unwrap_err();
    assert_eq!(err.kind(), ErrorKind::Sink);
    assert!(err.find_cause::<Full>().is_some());
}
