//! The transfer engine.
//!
//! An object is split into transfer *chunks* (at most `chunk_size`) grouped
//! into verification *units*: one per S3 multipart part, or the whole object
//! for MD5. Workers pull chunks in order, share a per-host connection budget
//! that adapts (AIMD), retry failing chunks with a per-chunk delay at the back
//! of the queue, and hedge the oldest slow chunk once nothing else is queued.
//! A unit that fails verification is repaired by re-fetching its chunks,
//! comparing per-chunk BLAKE3 digests, and rewriting only the chunks that
//! differ.

use crate::error::{Error, ErrorKind};
use crate::http;
use crate::ledger;
use crate::object::{self, Checksum, Layout, Object};
use crate::options::{FetchOptions, GetOptions};
use crate::sink::{Access, Opened, Sink, StreamReader, StreamSink};
use md5::{Digest as _, Md5};
use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};
use std::io::Read;
use std::sync::atomic::{AtomicBool, AtomicU8, AtomicU32, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex, MutexGuard};
use std::thread::JoinHandle;
use std::time::{Duration, Instant};
use tokio::sync::watch;

const PENDING: u8 = 0;
const IN_FLIGHT: u8 = 1;
const DONE: u8 = 2;
/// Fetched by the request that won the chunk, not yet on the sink.
const CLAIMED: u8 = 3;

/// User-visible events are batched to at least this many received bytes.
const PROGRESS_STEP: u64 = 64 * 1024;
const AIMD_INTERVAL: Duration = Duration::from_secs(5);
const LEDGER_INTERVAL: Duration = Duration::from_secs(1);
const TICK: Duration = Duration::from_millis(100);
const PROBE_TIMEOUT: Duration = Duration::from_secs(120);
const REPAIR_TIMEOUT: Duration = Duration::from_secs(600);

fn lock<T>(m: &Mutex<T>) -> MutexGuard<'_, T> {
    m.lock().unwrap_or_else(|e| e.into_inner())
}

fn try_lock<T>(m: &Mutex<T>) -> Option<MutexGuard<'_, T>> {
    match m.try_lock() {
        Ok(g) => Some(g),
        Err(std::sync::TryLockError::Poisoned(e)) => Some(e.into_inner()),
        Err(std::sync::TryLockError::WouldBlock) => None,
    }
}

/// Progress of a running download, delivered to the observer passed to
/// [`Fetcher::fetch_with`].
///
/// Observers run on worker threads and must be quick.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum Event {
    /// Body bytes received, in batches of at least 64 KiB (the last batch of a
    /// request may be smaller). Hedged duplicates count too.
    Received {
        /// Bytes since the previous `Received` event of the same request.
        bytes: u64,
    },
    /// A chunk was written to the sink.
    ChunkDone {
        /// Chunk start.
        offset: u64,
        /// Chunk length.
        len: u64,
    },
    /// A chunk request failed and will be retried.
    Retrying {
        /// Chunk start.
        offset: u64,
        /// Why it failed.
        error: Error,
    },
    /// A slow chunk was re-issued on another connection.
    Hedged {
        /// Chunk start.
        offset: u64,
    },
    /// Verification failed; the object is being re-fetched to find and
    /// rewrite corrupt chunks.
    Repairing,
}

type Observer = Arc<dyn Fn(&Event) + Send + Sync>;

/// Statistics of a finished download.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct Report {
    /// Object size in bytes.
    pub size: u64,
    /// The sink already held the verified object; nothing was fetched.
    pub reused: bool,
    /// Chunks taken from a resume checkpoint instead of being fetched.
    pub resumed_chunks: usize,
    /// Wall time.
    pub elapsed: Duration,
    /// Failed chunk requests that were retried.
    pub chunk_retries: u64,
    /// Chunks re-issued because they were slow.
    pub hedged: u64,
    /// Chunks rewritten by repair after a verification failure.
    pub repaired_chunks: u64,
}

/// A successful download: the sink's output and statistics.
#[derive(Debug)]
#[non_exhaustive]
pub struct Download<T> {
    /// What the sink produced (a path, bytes, `()`).
    pub output: T,
    /// Statistics.
    pub report: Report,
}

/// Downloads objects. Cheap to clone; clones share connections and per-host
/// connection budgets.
#[derive(Clone)]
pub struct Fetcher {
    shared: Arc<Shared>,
}

impl std::fmt::Debug for Fetcher {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Fetcher")
            .field("options", &self.shared.options)
            .finish_non_exhaustive()
    }
}

struct Shared {
    options: FetchOptions,
    client: reqwest::Client,
    runtime: Option<tokio::runtime::Runtime>,
    hosts: Mutex<HashMap<String, Arc<Host>>>,
    seq: AtomicU64,
}

impl Drop for Shared {
    fn drop(&mut self) {
        if let Some(rt) = self.runtime.take() {
            rt.shutdown_background();
        }
    }
}

impl Shared {
    fn rt(&self) -> &tokio::runtime::Runtime {
        self.runtime.as_ref().expect("runtime present until drop")
    }

    /// Run `fut` on the I/O runtime and wait for it; safe to call from any
    /// thread, including inside another async runtime.
    fn block<T: Send + 'static>(
        &self,
        fut: impl std::future::Future<Output = T> + Send + 'static,
    ) -> Result<T, Error> {
        let (tx, rx) = std::sync::mpsc::channel();
        self.rt().spawn(async move {
            let _ = tx.send(fut.await);
        });
        rx.recv()
            .map_err(|_| Error::new(ErrorKind::Internal, "I/O runtime stopped"))
    }

    fn host(&self, key: &str, name: &str) -> Arc<Host> {
        lock(&self.hosts)
            .entry(key.to_string())
            .or_insert_with(|| {
                let max = self.options.host_max_connections(name);
                Arc::new(Host::new(
                    max,
                    self.options.initial_connections.clamp(1, max),
                ))
            })
            .clone()
    }
}

impl Fetcher {
    /// A fetcher with `options`. Starts a small I/O runtime (two threads).
    pub fn new(options: FetchOptions) -> Result<Self, Error> {
        let client = reqwest::Client::builder()
            .user_agent(options.user_agent.clone())
            .connect_timeout(options.connect_timeout)
            .build()
            .map_err(|e| http::request_error(e, "building HTTP client"))?;
        let runtime = tokio::runtime::Builder::new_multi_thread()
            .worker_threads(2)
            .thread_name("readkit-fetch-io")
            .enable_all()
            .build()?;
        Ok(Self {
            shared: Arc::new(Shared {
                options,
                client,
                runtime: Some(runtime),
                hosts: Mutex::new(HashMap::new()),
                seq: AtomicU64::new(0),
            }),
        })
    }

    /// The options this fetcher was built with.
    pub fn options(&self) -> &FetchOptions {
        &self.shared.options
    }

    /// Start downloading `object` into `sink`. Returns at once; the download
    /// runs on its own threads until the handle is joined, cancelled or
    /// dropped. Downloads started earlier get connections first.
    pub fn fetch<S: Sink>(&self, object: Object, sink: S) -> FetchHandle<S::Output> {
        self.start(object, sink, None)
    }

    /// Like [`fetch`](Self::fetch), calling `observer` with progress events.
    pub fn fetch_with<S: Sink>(
        &self,
        object: Object,
        sink: S,
        observer: impl Fn(&Event) + Send + Sync + 'static,
    ) -> FetchHandle<S::Output> {
        self.start(object, sink, Some(Arc::new(observer)))
    }

    /// Download `object` as an in-order byte stream with about `buffer` bytes
    /// of read-ahead (see [`StreamSink`]).
    ///
    /// Only chunks within `buffer` of the reader are fetched, so `buffer`
    /// bounds the parallelism: about `buffer / chunk_size` connections. To use
    /// every connection to S3 (32 by default, 16 MiB chunks) that is 512 MiB.
    pub fn stream(&self, object: Object, buffer: usize) -> FetchReader {
        let (sink, reader) = StreamSink::new(buffer);
        FetchReader {
            reader,
            handle: Some(self.fetch(object, sink)),
        }
    }

    /// GET `url` and return its body, with [`GetOptions::default`]: at most
    /// 64 MiB, 120 s per attempt, 3 retries. Blocks the calling thread (any
    /// thread, including one inside an async runtime).
    ///
    /// For small bodies of unknown size, such as API responses; use
    /// [`fetch`](Self::fetch) for objects of known size and checksum.
    /// Errors: [`ErrorKind::NotFound`] (404, 410), [`ErrorKind::Permanent`]
    /// (other 4xx), [`ErrorKind::TooLarge`], and, once retries are spent,
    /// [`ErrorKind::Throttled`] (429, 503), [`ErrorKind::Transient`] (5xx,
    /// network) or [`ErrorKind::Timeout`]. [`Error::status`] holds the HTTP
    /// status when there was one.
    pub fn get(&self, url: &str) -> Result<Vec<u8>, Error> {
        self.get_with(url, &GetOptions::default())
    }

    /// Like [`get`](Self::get), with a size limit, deadline and retry count
    /// from `options`.
    ///
    /// ```no_run
    /// use readkit_fetch::{FetchOptions, Fetcher, GetOptions};
    ///
    /// let fetcher = Fetcher::new(FetchOptions::default())?;
    /// let limit = GetOptions { max_bytes: 1 << 20, ..Default::default() };
    /// let body = fetcher.get_with("https://example.org/manifest.json", &limit)?;
    /// # Ok::<(), readkit_fetch::Error>(())
    /// ```
    pub fn get_with(&self, url: &str, options: &GetOptions) -> Result<Vec<u8>, Error> {
        let mut attempt = 0u32;
        loop {
            let client = self.shared.client.clone();
            let owned = url.to_string();
            let (max, timeout) = (options.max_bytes, options.timeout);
            let result = self
                .shared
                .block(async move { http::get(&client, &owned, max, timeout).await })?;
            match result {
                Ok(body) => return Ok(body),
                Err((e, retry_after)) if e.kind().is_retryable() && attempt < options.retries => {
                    let o = &self.shared.options;
                    let backoff = o
                        .retry_delay
                        .saturating_mul(1u32 << attempt.min(16))
                        .min(o.max_retry_delay);
                    let wait = retry_after
                        .map_or(backoff, |r| r.max(backoff))
                        .min(o.max_retry_delay);
                    attempt += 1;
                    std::thread::sleep(wait);
                }
                Err((e, _)) => return Err(e),
            }
        }
    }

    /// Describe an S3 object from a `HEAD`: its size and ETag. For a
    /// multipart ETag (`-N` suffix, the part count) the part size comes from a
    /// second `HEAD` with `?partNumber=1`.
    pub fn probe_s3(&self, url: &str) -> Result<Object, Error> {
        let client = self.shared.client.clone();
        let owned = url.to_string();
        let (size, etag) = self
            .shared
            .block(async move { http::head(&client, &owned, PROBE_TIMEOUT).await })??;
        let etag = etag
            .ok_or_else(|| Error::new(ErrorKind::Transient, format!("HEAD {url} without ETag")))?;
        let checksum = match Checksum::from_s3_etag(&etag, size)? {
            Checksum::S3MultipartEtag { etag, parts, .. } if parts > 1 => {
                let sep = if url.contains('?') { '&' } else { '?' };
                let part_url = format!("{url}{sep}partNumber=1");
                let client = self.shared.client.clone();
                let (part_size, _) = self
                    .shared
                    .block(async move { http::head(&client, &part_url, PROBE_TIMEOUT).await })??;
                Checksum::S3MultipartEtag {
                    etag,
                    parts,
                    part_size,
                }
            }
            c => c,
        };
        Ok(Object::new(url, size, checksum))
    }

    fn start<S: Sink>(
        &self,
        object: Object,
        sink: S,
        observer: Option<Observer>,
    ) -> FetchHandle<S::Output> {
        let (cancel_tx, _) = watch::channel(false);
        let control = Arc::new(Control {
            cancelled: AtomicBool::new(false),
            cancel_tx,
            wake: Mutex::new(()),
            cv: Condvar::new(),
        });
        let seq = self.shared.seq.fetch_add(1, Ordering::SeqCst);
        let shared = self.shared.clone();
        let ctl = control.clone();
        let spawned = std::thread::Builder::new()
            .name("readkit-fetch".into())
            .spawn(move || drive(&shared, &ctl, seq, object, sink, observer));
        FetchHandle {
            control,
            inner: Some(match spawned {
                Ok(t) => Inner::Running(t),
                Err(e) => Inner::Failed(e.into()),
            }),
        }
    }
}

/// Re-verify a local file against `object`'s checksum.
pub fn verify_file(object: &Object, path: &std::path::Path) -> Result<(), Error> {
    let mut file = std::fs::File::open(path).map_err(|e| {
        Error::new(ErrorKind::Io, format!("opening {}: {e}", path.display())).with_source(e)
    })?;
    let len = file.metadata()?.len();
    if len != object.size {
        return Err(Error::new(
            ErrorKind::Verification,
            format!("{}: size {len} != expected {}", path.display(), object.size),
        ));
    }
    let layout = object::plan(object, 8 << 20, 1)?;
    let mut buf = vec![0u8; 4 << 20];
    let mut digests = Vec::with_capacity(layout.units.len());
    // Units are contiguous and in order, so one sequential pass suffices.
    for unit in &layout.units {
        let mut h = Md5::new();
        let mut left = unit.len;
        while left > 0 {
            let n = left.min(buf.len() as u64) as usize;
            file.read_exact(&mut buf[..n])?;
            h.update(&buf[..n]);
            left -= n as u64;
        }
        digests.push(h.finalize().into());
    }
    object::verify(object, &digests)
}

/// A running download. Join it for the result; cancel or drop it to stop the
/// download promptly (in-flight requests are aborted and every worker thread
/// is joined).
pub struct FetchHandle<T> {
    control: Arc<Control>,
    inner: Option<Inner<T>>,
}

enum Inner<T> {
    Running(JoinHandle<Result<Download<T>, Error>>),
    Failed(Error),
}

impl<T> std::fmt::Debug for FetchHandle<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FetchHandle")
            .field("finished", &self.is_finished())
            .finish()
    }
}

impl<T> FetchHandle<T> {
    /// Ask the download to stop; [`join`](Self::join) then returns
    /// [`ErrorKind::Cancelled`] (unless it had already finished). File
    /// downloads keep their resume state.
    pub fn cancel(&self) {
        self.control.cancel();
    }

    /// Whether the download has finished (successfully or not).
    pub fn is_finished(&self) -> bool {
        match &self.inner {
            Some(Inner::Running(t)) => t.is_finished(),
            _ => true,
        }
    }

    /// Wait for the download to finish.
    pub fn join(mut self) -> Result<Download<T>, Error> {
        match self.inner.take() {
            Some(Inner::Running(t)) => t
                .join()
                .unwrap_or_else(|_| Err(Error::new(ErrorKind::Internal, "fetch thread panicked"))),
            Some(Inner::Failed(e)) => Err(e),
            None => Err(Error::new(ErrorKind::Internal, "handle already joined")),
        }
    }
}

impl<T> Drop for FetchHandle<T> {
    fn drop(&mut self) {
        if let Some(Inner::Running(t)) = self.inner.take() {
            self.control.cancel();
            let _ = t.join();
        }
    }
}

/// A download as an in-order `impl Read` (see [`Fetcher::stream`]).
///
/// Reads yield the object's bytes in order, then `Ok(0)` once it is complete
/// and verified, or an error if it failed (including verification). Dropping
/// it cancels the download.
#[derive(Debug)]
pub struct FetchReader {
    reader: StreamReader,
    handle: Option<FetchHandle<()>>,
}

impl FetchReader {
    /// Cancel the download; pending and later reads fail.
    pub fn cancel(&self) {
        if let Some(h) = &self.handle {
            h.cancel();
        }
    }

    /// Wait for the download and return its statistics. Cancels it first if
    /// the stream was not read to the end.
    pub fn finish(mut self) -> Result<Report, Error> {
        let handle = self.handle.take().expect("handle present until finish");
        if !self.reader.ended() {
            handle.cancel();
        }
        handle.join().map(|d| d.report)
    }
}

impl Read for FetchReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.reader.read(buf)
    }
}

struct Control {
    cancelled: AtomicBool,
    cancel_tx: watch::Sender<bool>,
    wake: Mutex<()>,
    cv: Condvar,
}

impl Control {
    fn cancel(&self) {
        self.cancelled.store(true, Ordering::SeqCst);
        self.cancel_tx.send_replace(true);
        self.notify();
    }

    fn notify(&self) {
        let _g = lock(&self.wake);
        self.cv.notify_all();
    }
}

/// Per-host connection budget shared by every download to that host.
struct Host {
    max: usize,
    slots: Mutex<Slots>,
    cv: Condvar,
    bytes: AtomicU64,
    errors: AtomicU64,
    last_progress: Mutex<Instant>,
    /// EWMA bytes/s per connection, used for chunk deadlines.
    per_conn_rate: Mutex<f64>,
    aimd: Mutex<Aimd>,
}

struct Slots {
    limit: usize,
    active: usize,
    jobs: usize,
    /// Waiting workers per download sequence number; the oldest download goes
    /// first.
    waiting: BTreeMap<u64, usize>,
}

struct Aimd {
    at: Instant,
    bytes: u64,
    rate: f64,
    errors: u64,
}

impl Host {
    fn new(max: usize, initial: usize) -> Self {
        Self {
            max,
            slots: Mutex::new(Slots {
                limit: initial,
                active: 0,
                jobs: 0,
                waiting: BTreeMap::new(),
            }),
            cv: Condvar::new(),
            bytes: AtomicU64::new(0),
            errors: AtomicU64::new(0),
            last_progress: Mutex::new(Instant::now()),
            per_conn_rate: Mutex::new(256.0 * 1024.0),
            aimd: Mutex::new(Aimd {
                at: Instant::now(),
                bytes: 0,
                rate: 0.0,
                errors: 0,
            }),
        }
    }

    fn begin_job(&self) {
        let first = {
            let mut s = lock(&self.slots);
            s.jobs += 1;
            s.jobs == 1
        };
        if first {
            // An idle host's clocks are stale.
            self.touch();
            let mut a = lock(&self.aimd);
            a.at = Instant::now();
            a.bytes = self.bytes.load(Ordering::SeqCst);
            a.errors = self.errors.load(Ordering::SeqCst);
        }
    }

    fn end_job(&self) {
        lock(&self.slots).jobs -= 1;
    }

    fn touch(&self) {
        *lock(&self.last_progress) = Instant::now();
    }

    /// Take a connection slot, waiting for one; `false` if `halt` is set first.
    fn acquire(&self, seq: u64, halt: &AtomicBool) -> bool {
        let mut s = lock(&self.slots);
        *s.waiting.entry(seq).or_insert(0) += 1;
        let ok = loop {
            if halt.load(Ordering::SeqCst) {
                break false;
            }
            if s.active < s.limit && s.waiting.keys().next() == Some(&seq) {
                s.active += 1;
                break true;
            }
            s = self
                .cv
                .wait_timeout(s, TICK)
                .unwrap_or_else(|e| e.into_inner())
                .0;
        };
        let w = s.waiting.get_mut(&seq).expect("registered");
        *w -= 1;
        if *w == 0 {
            s.waiting.remove(&seq);
        }
        self.cv.notify_all();
        ok
    }

    fn release(&self) {
        lock(&self.slots).active -= 1;
        self.cv.notify_all();
    }

    fn wake_all(&self) {
        let _s = lock(&self.slots);
        self.cv.notify_all();
    }

    /// AIMD, at most every 5 s: halve the limit after throttling, add two
    /// connections while throughput still grows by 5 %.
    fn maybe_adjust(&self) {
        let mut a = lock(&self.aimd);
        let elapsed = a.at.elapsed();
        if elapsed < AIMD_INTERVAL {
            return;
        }
        let bytes = self.bytes.load(Ordering::SeqCst);
        let rate = bytes.saturating_sub(a.bytes) as f64 / elapsed.as_secs_f64();
        let errors = self.errors.load(Ordering::SeqCst);
        let active = {
            let mut s = lock(&self.slots);
            if errors > a.errors {
                s.limit = (s.limit / 2).max(1);
            } else if rate >= a.rate * 1.05 && s.limit < self.max {
                s.limit = (s.limit + 2).min(self.max);
            }
            s.active.max(1) as f64
        };
        a.at = Instant::now();
        a.bytes = bytes;
        a.rate = rate;
        a.errors = errors;
        if rate > 0.0 {
            let mut r = lock(&self.per_conn_rate);
            *r = 0.5 * *r + 0.5 * (rate / active);
        }
        self.cv.notify_all();
    }
}

fn host_key(url: &str) -> Result<(String, String), Error> {
    let parsed = reqwest::Url::parse(url).map_err(|e| {
        Error::new(ErrorKind::InvalidInput, format!("invalid URL {url:?}: {e}")).with_source(e)
    })?;
    let host = parsed
        .host_str()
        .ok_or_else(|| Error::new(ErrorKind::InvalidInput, format!("URL {url:?} has no host")))?
        .to_ascii_lowercase();
    let port = parsed.port_or_known_default().unwrap_or(0);
    Ok((format!("{}://{host}:{port}", parsed.scheme()), host))
}

/// State of an in-order (sequential) sink.
struct Seq {
    window: u64,
    written: AtomicU64,
    reorder: Mutex<BTreeMap<usize, Vec<u8>>>,
    drain: Mutex<Drain>,
    consumer_blocked: AtomicBool,
}

struct Drain {
    next: usize,
    hasher: Option<Md5>,
}

/// A unit's MD5 over its chunks `first_chunk..next`.
struct UnitHash {
    next: usize,
    hasher: Md5,
}

struct Job<'a, S: Sink> {
    shared: &'a Shared,
    control: &'a Control,
    host: &'a Host,
    object: &'a Object,
    observer: Option<Observer>,
    seq_no: u64,
    layout: Layout,
    sink: S,
    seq: Option<Seq>,
    verifying: bool,
    state: Vec<AtomicU8>,
    started: Vec<Mutex<Option<Instant>>>,
    hedges: Vec<AtomicU8>,
    /// Seconds per byte of chunks fetched without a hedge, and their median.
    chunk_spb: Mutex<(Vec<f64>, Option<f64>)>,
    attempts: Vec<AtomicU32>,
    retry_at: Vec<Mutex<Option<Instant>>>,
    blake3: Vec<Mutex<Option<[u8; 32]>>>,
    unit_md5: Vec<Mutex<Option<[u8; 16]>>>,
    /// In-order MD5 of each unit, advanced while the download runs
    /// ([`Access::Random`] sinks).
    unit_hash: Vec<Mutex<UnitHash>>,
    /// Set by a worker that found a unit's hasher busy: the holder looks again.
    unit_kick: Vec<AtomicBool>,
    queue: Mutex<VecDeque<usize>>,
    remaining: AtomicUsize,
    resumed: usize,
    dirty: AtomicBool,
    retries: AtomicU64,
    hedged: AtomicU64,
    failed: Mutex<Option<Error>>,
    last_error: Mutex<Option<Error>>,
    halt: AtomicBool,
    halt_tx: watch::Sender<bool>,
    last_progress: Mutex<Instant>,
}

fn drive<S: Sink>(
    shared: &Shared,
    control: &Control,
    seq_no: u64,
    object: Object,
    mut sink: S,
    observer: Option<Observer>,
) -> Result<Download<S::Output>, Error> {
    let t0 = Instant::now();
    // At least one chunk per connection the host allows (and 16), so small
    // objects still spread over every connection.
    let prepared = host_key(&object.url)
        .and_then(|(key, name)| {
            let min_chunks = shared.options.host_max_connections(&name).max(16) as u64;
            let layout = object::plan(&object, shared.options.chunk_size, min_chunks)?;
            Ok((layout, (key, name)))
        })
        .and_then(|prepared| Ok((prepared, sink.open(&object)?)));
    let ((layout, (key, name)), opened) = match prepared {
        Ok(p) => p,
        Err(e) => {
            sink.abort(&e);
            return Err(e);
        }
    };
    if let Opened::Complete = opened {
        let output = sink.commit()?;
        return Ok(Download {
            output,
            report: Report {
                size: object.size,
                reused: true,
                resumed_chunks: 0,
                elapsed: t0.elapsed(),
                chunk_retries: 0,
                hedged: 0,
                repaired_chunks: 0,
            },
        });
    }
    let access = sink.access();
    let resume = match (opened, access) {
        (Opened::Resume(bytes), Access::Random) => ledger::parse(&bytes, &object, &layout),
        _ => None,
    };
    let host = shared.host(&key, &name);
    host.begin_job();
    struct EndJob<'h>(&'h Host);
    impl Drop for EndJob<'_> {
        fn drop(&mut self) {
            self.0.end_job();
        }
    }
    let _end = EndJob(&host);
    let job = Job::new(
        shared, control, &host, &object, observer, seq_no, layout, sink, access, resume,
    );
    if job.remaining.load(Ordering::SeqCst) > 0 && !control.cancelled.load(Ordering::SeqCst) {
        // One worker beyond the connection cap: a worker that is hashing has
        // released its slot, and the spare keeps that slot busy.
        let workers = (host.max + 1)
            .min(job.remaining.load(Ordering::SeqCst))
            .max(1);
        std::thread::scope(|scope| {
            for _ in 0..workers {
                let spawned = std::thread::Builder::new()
                    .name("readkit-fetch-worker".into())
                    .spawn_scoped(scope, || job.worker());
                if let Err(e) = spawned {
                    job.fail(e.into());
                    break;
                }
            }
            job.supervise();
        });
    }
    if control.cancelled.load(Ordering::SeqCst) {
        lock(&job.failed).get_or_insert_with(Error::cancelled);
    }
    job.finish(t0)
}

impl<'a, S: Sink> Job<'a, S> {
    #[allow(clippy::too_many_arguments)]
    fn new(
        shared: &'a Shared,
        control: &'a Control,
        host: &'a Host,
        object: &'a Object,
        observer: Option<Observer>,
        seq_no: u64,
        layout: Layout,
        sink: S,
        access: Access,
        resume: Option<ledger::Resume>,
    ) -> Self {
        let n = layout.chunks.len();
        let done = |i: usize| resume.as_ref().is_some_and(|r| r.done[i]);
        let state: Vec<AtomicU8> = (0..n)
            .map(|i| AtomicU8::new(if done(i) { DONE } else { PENDING }))
            .collect();
        let resumed = (0..n).filter(|i| done(*i)).count();
        let queue = (0..n).filter(|i| !done(*i)).collect::<VecDeque<_>>();
        let blake3 = (0..n)
            .map(|i| Mutex::new(resume.as_ref().and_then(|r| r.blake3[i])))
            .collect();
        let seq = match access {
            Access::Sequential { window } => Some(Seq {
                window: window.max(1),
                written: AtomicU64::new(0),
                reorder: Mutex::new(BTreeMap::new()),
                drain: Mutex::new(Drain {
                    next: 0,
                    hasher: None,
                }),
                consumer_blocked: AtomicBool::new(false),
            }),
            _ => None,
        };
        let units = layout.units.len();
        Self {
            shared,
            control,
            host,
            object,
            observer,
            seq_no,
            verifying: object.checksum.verifies(),
            state,
            started: (0..n).map(|_| Mutex::new(None)).collect(),
            hedges: (0..n).map(|_| AtomicU8::new(0)).collect(),
            chunk_spb: Mutex::new((Vec::new(), None)),
            attempts: (0..n).map(|_| AtomicU32::new(0)).collect(),
            retry_at: (0..n).map(|_| Mutex::new(None)).collect(),
            blake3,
            unit_md5: (0..units).map(|_| Mutex::new(None)).collect(),
            unit_hash: layout
                .units
                .iter()
                .map(|u| {
                    Mutex::new(UnitHash {
                        next: u.first_chunk,
                        hasher: Md5::new(),
                    })
                })
                .collect(),
            unit_kick: (0..units).map(|_| AtomicBool::new(false)).collect(),
            remaining: AtomicUsize::new(n - resumed),
            resumed,
            queue: Mutex::new(queue),
            dirty: AtomicBool::new(true),
            retries: AtomicU64::new(0),
            hedged: AtomicU64::new(0),
            failed: Mutex::new(None),
            last_error: Mutex::new(None),
            halt: AtomicBool::new(false),
            halt_tx: watch::channel(false).0,
            last_progress: Mutex::new(Instant::now()),
            layout,
            sink,
            seq,
        }
    }

    fn emit(&self, event: impl FnOnce() -> Event) {
        if let Some(o) = &self.observer {
            o(&event());
        }
    }

    fn opts(&self) -> &FetchOptions {
        &self.shared.options
    }

    fn halted(&self) -> bool {
        self.halt.load(Ordering::SeqCst)
    }

    fn fail(&self, error: Error) {
        lock(&self.failed).get_or_insert(error);
        self.control.notify();
    }

    fn touch(&self) {
        *lock(&self.last_progress) = Instant::now();
    }

    /// Wait up to `d` for a notification (or halt).
    fn idle(&self, d: Duration) {
        let g = lock(&self.control.wake);
        if !self.halted() {
            let _ = self.control.cv.wait_timeout(g, d);
        }
    }

    /// Sleep `d` unless halted.
    fn pause(&self, d: Duration) {
        let end = Instant::now() + d;
        let mut g = lock(&self.control.wake);
        while !self.halted() {
            let now = Instant::now();
            if now >= end {
                break;
            }
            g = self
                .control
                .cv
                .wait_timeout(g, end - now)
                .unwrap_or_else(|e| e.into_inner())
                .0;
        }
    }

    fn supervise(&self) {
        let mut last_flush = Instant::now();
        loop {
            {
                let g = lock(&self.control.wake);
                if !self.should_stop() {
                    let _ = self.control.cv.wait_timeout(g, TICK);
                }
            }
            if self.control.cancelled.load(Ordering::SeqCst) {
                self.fail(Error::cancelled());
            }
            if self.should_stop() {
                break;
            }
            if last_flush.elapsed() >= LEDGER_INTERVAL {
                self.flush_ledger(false);
                last_flush = Instant::now();
            }
            self.host.maybe_adjust();
            if let Some(seq) = &self.seq {
                // Waiting on a slow consumer is not a stall.
                if seq.consumer_blocked.load(Ordering::SeqCst) {
                    self.touch();
                }
            }
            let t = self.opts().no_progress_timeout;
            if lock(&self.last_progress).elapsed() > t
                && lock(&self.host.last_progress).elapsed() > t
            {
                let mut e = Error::new(
                    ErrorKind::Timeout,
                    format!(
                        "{}: no forward progress for {:.1}s",
                        self.object.url,
                        t.as_secs_f64()
                    ),
                );
                if let Some(last) = lock(&self.last_error).clone() {
                    e = e.with_source(last);
                }
                self.fail(e);
                break;
            }
        }
        // Stop workers: abort in-flight requests and wake every waiter.
        self.halt.store(true, Ordering::SeqCst);
        self.halt_tx.send_replace(true);
        self.control.notify();
        self.host.wake_all();
        if lock(&self.failed).is_some() {
            self.sink.interrupt();
        }
    }

    fn should_stop(&self) -> bool {
        self.control.cancelled.load(Ordering::SeqCst)
            || lock(&self.failed).is_some()
            || self.remaining.load(Ordering::SeqCst) == 0
    }

    fn worker(&self) {
        let run = std::panic::AssertUnwindSafe(|| self.worker_loop());
        if std::panic::catch_unwind(run).is_err() {
            // Never hang the caller on a worker bug.
            self.fail(Error::new(ErrorKind::Internal, "fetch worker panicked"));
        }
    }

    fn worker_loop(&self) {
        let mut halt_rx = self.halt_tx.subscribe();
        loop {
            if self.halted() || lock(&self.failed).is_some() {
                return;
            }
            let Some((ci, hedge)) = self.next_task() else {
                self.idle(TICK);
                continue;
            };
            if !self.host.acquire(self.seq_no, &self.halt) {
                if !hedge {
                    self.unclaim(ci);
                }
                return;
            }
            if hedge {
                self.emit(|| Event::Hedged {
                    offset: self.layout.chunks[ci].offset,
                });
            } else {
                *lock(&self.started[ci]) = Some(Instant::now());
            }
            let result = self.fetch_chunk(ci, &mut halt_rx);
            self.host.release();
            match result {
                Ok(Some(buffer)) => {
                    if !hedge {
                        self.record_chunk_time(ci);
                    }
                    if let Err(e) = self.complete_chunk(ci, buffer) {
                        self.fail(e);
                    }
                }
                Ok(None) => {} // lost a hedge race
                Err(e) if e.kind() == ErrorKind::Cancelled => {
                    if !hedge {
                        self.unclaim(ci);
                    }
                    return;
                }
                Err(e) if !e.kind().is_retryable() => self.fail(e),
                Err(e) => self.retry(ci, hedge, e),
            }
        }
    }

    fn unclaim(&self, ci: usize) {
        if self.state[ci]
            .compare_exchange(IN_FLIGHT, PENDING, Ordering::SeqCst, Ordering::SeqCst)
            .is_ok()
        {
            lock(&self.queue).push_front(ci);
        }
    }

    fn retry(&self, ci: usize, hedge: bool, error: Error) {
        self.retries.fetch_add(1, Ordering::SeqCst);
        let throttled = error.kind() == ErrorKind::Throttled;
        if throttled {
            self.host.errors.fetch_add(1, Ordering::SeqCst);
        }
        self.emit(|| Event::Retrying {
            offset: self.layout.chunks[ci].offset,
            error: error.clone(),
        });
        *lock(&self.last_error) = Some(error);
        if !hedge {
            // Per-chunk exponential delay; the chunk goes to the back of the
            // queue so other chunks keep flowing and no worker sleeps on a bad
            // chunk.
            let attempts = self.attempts[ci].fetch_add(1, Ordering::SeqCst).min(6);
            let delay =
                (self.opts().retry_delay * (1u32 << attempts)).min(self.opts().max_retry_delay);
            *lock(&self.retry_at[ci]) = Some(Instant::now() + delay);
            if self.state[ci]
                .compare_exchange(IN_FLIGHT, PENDING, Ordering::SeqCst, Ordering::SeqCst)
                .is_ok()
            {
                lock(&self.queue).push_back(ci);
            }
        }
        if throttled {
            // Host-level throttling: back this connection off briefly.
            self.pause(self.opts().retry_delay);
        }
    }

    fn admissible(&self, ci: usize) -> bool {
        match &self.seq {
            Some(seq) => {
                self.layout.chunks[ci].offset < seq.written.load(Ordering::SeqCst) + seq.window
            }
            None => true,
        }
    }

    fn record_chunk_time(&self, ci: usize) {
        let len = self.layout.chunks[ci].len;
        let Some(started) = *lock(&self.started[ci]) else {
            return;
        };
        if len == 0 {
            return;
        }
        let mut g = lock(&self.chunk_spb);
        let (samples, median) = &mut *g;
        samples.push(started.elapsed().as_secs_f64() / len as f64);
        let mut sorted = samples.clone();
        sorted.sort_by(f64::total_cmp);
        *median = Some(sorted[sorted.len() / 2]);
    }

    /// Age past which in-flight chunk `ci` is hedged: `hedge_after`, or three
    /// times as long as this download's chunks usually take, if longer. On a
    /// saturated link every chunk is slow; duplicating them only takes
    /// bandwidth from the originals.
    fn hedge_threshold(&self, ci: usize) -> Duration {
        let typical = lock(&self.chunk_spb)
            .1
            .map_or(0.0, |spb| 3.0 * spb * self.layout.chunks[ci].len as f64);
        self.opts()
            .hedge_after
            .max(Duration::from_secs_f64(typical))
    }

    /// Pop the next pending chunk, or hedge the oldest in-flight chunk once
    /// nothing is ready.
    fn next_task(&self) -> Option<(usize, bool)> {
        {
            let mut queue = lock(&self.queue);
            for _ in 0..queue.len() {
                let ci = queue.pop_front()?;
                if self.state[ci].load(Ordering::SeqCst) != PENDING {
                    continue; // stale entry
                }
                let waiting = lock(&self.retry_at[ci]).is_some_and(|t| t > Instant::now());
                if waiting || !self.admissible(ci) {
                    queue.push_back(ci);
                    continue;
                }
                if self.state[ci]
                    .compare_exchange(PENDING, IN_FLIGHT, Ordering::SeqCst, Ordering::SeqCst)
                    .is_ok()
                {
                    return Some((ci, false));
                }
            }
        }
        let mut best: Option<(usize, Duration)> = None;
        for ci in 0..self.layout.chunks.len() {
            if self.state[ci].load(Ordering::SeqCst) != IN_FLIGHT
                || self.hedges[ci].load(Ordering::SeqCst) > 0
            {
                continue;
            }
            if let Some(started) = *lock(&self.started[ci]) {
                let age = started.elapsed();
                if age > self.hedge_threshold(ci) && best.is_none_or(|b| age > b.1) {
                    best = Some((ci, age));
                }
            }
        }
        let (ci, _) = best?;
        if self.hedges[ci].swap(1, Ordering::SeqCst) > 0 {
            return None;
        }
        self.hedged.fetch_add(1, Ordering::SeqCst);
        Some((ci, true))
    }

    fn fetch_chunk(
        &self,
        ci: usize,
        halt_rx: &mut watch::Receiver<bool>,
    ) -> Result<Option<Vec<u8>>, Error> {
        let spec = self.layout.chunks[ci];
        let rate = *lock(&self.host.per_conn_rate);
        let deadline = Duration::from_secs_f64(spec.len as f64 / (0.25 * rate.max(1.0)))
            .max(self.opts().min_chunk_timeout);
        let mut unreported = 0u64;
        let result = {
            let mut progress = |n: u64| {
                self.host.touch();
                self.touch();
                unreported += n;
                if unreported >= PROGRESS_STEP {
                    let bytes = std::mem::take(&mut unreported);
                    self.emit(|| Event::Received { bytes });
                }
            };
            let fut = http::fetch_range(
                &self.shared.client,
                &self.object.url,
                spec.offset,
                spec.len,
                deadline,
                &mut progress,
            );
            self.shared
                .rt()
                .handle()
                .block_on(http::or_cancel(fut, halt_rx))
        };
        if unreported > 0 {
            self.emit(|| Event::Received { bytes: unreported });
        }
        let buffer = result.ok_or_else(Error::cancelled)??;
        self.host.bytes.fetch_add(spec.len, Ordering::SeqCst);
        self.host.touch();
        self.touch();
        // Winner takes the chunk; a losing duplicate discards its copy.
        let won = [IN_FLIGHT, PENDING].into_iter().any(|from| {
            self.state[ci]
                .compare_exchange(from, CLAIMED, Ordering::SeqCst, Ordering::SeqCst)
                .is_ok()
        });
        Ok(won.then_some(buffer))
    }

    fn complete_chunk(&self, ci: usize, buffer: Vec<u8>) -> Result<(), Error> {
        let spec = self.layout.chunks[ci];
        let unit = self.layout.units[spec.unit];
        match &self.seq {
            None => {
                if self.verifying {
                    *lock(&self.blake3[ci]) = Some(*blake3::hash(&buffer).as_bytes());
                    // The unit's next chunk is hashed from memory, before the
                    // write; later chunks are read back as the cursor reaches them.
                    if let Some(mut h) = try_lock(&self.unit_hash[spec.unit])
                        && h.next == ci
                    {
                        h.hasher.update(&buffer);
                        h.next += 1;
                        self.finish_unit_hash(spec.unit, &mut h);
                    }
                }
                self.sink.write_at(spec.offset, buffer)?;
                self.state[ci].store(DONE, Ordering::SeqCst);
                self.dirty.store(true, Ordering::SeqCst);
                if self.verifying {
                    self.advance_unit_hash(spec.unit)?;
                }
            }
            Some(seq) => {
                if self.verifying && unit.chunks == 1 {
                    // Whole-unit chunks hash here, in parallel, not in drain.
                    *lock(&self.unit_md5[spec.unit]) = Some(Md5::digest(&buffer).into());
                }
                lock(&seq.reorder).insert(ci, buffer);
                self.drain(seq)?;
                self.state[ci].store(DONE, Ordering::SeqCst);
            }
        }
        self.emit(|| Event::ChunkDone {
            offset: spec.offset,
            len: spec.len,
        });
        self.remaining.fetch_sub(1, Ordering::SeqCst);
        self.control.notify();
        Ok(())
    }

    /// Write every chunk that is next in file order, hashing units as they
    /// pass. Only one thread drains at a time; a thread that finds nothing to
    /// do knows its chunk was written by the previous drainer.
    fn drain(&self, seq: &Seq) -> Result<(), Error> {
        let mut d = lock(&seq.drain);
        loop {
            if self.halted() {
                return Ok(());
            }
            let ci = d.next;
            let Some(data) = lock(&seq.reorder).remove(&ci) else {
                return Ok(());
            };
            d.next += 1;
            let spec = self.layout.chunks[ci];
            if self.verifying && self.layout.units[spec.unit].chunks > 1 {
                let unit = self.layout.units[spec.unit];
                if ci == unit.first_chunk {
                    d.hasher = Some(Md5::new());
                }
                if let Some(h) = d.hasher.as_mut() {
                    h.update(&data);
                }
                if ci + 1 == unit.first_chunk + unit.chunks {
                    let digest = d.hasher.take().map(|h| h.finalize().into());
                    *lock(&self.unit_md5[spec.unit]) = digest;
                }
            }
            seq.consumer_blocked.store(true, Ordering::SeqCst);
            let written = self.sink.write_at(spec.offset, data);
            seq.consumer_blocked.store(false, Ordering::SeqCst);
            written?;
            seq.written.store(spec.offset + spec.len, Ordering::SeqCst);
            self.control.notify();
        }
    }

    /// Advance unit `u`'s MD5 over chunks already on the sink, unless another
    /// worker is doing so (it is told to look again).
    fn advance_unit_hash(&self, u: usize) -> Result<(), Error> {
        self.unit_kick[u].store(true, Ordering::SeqCst);
        loop {
            let Some(mut h) = try_lock(&self.unit_hash[u]) else {
                return Ok(());
            };
            self.unit_kick[u].store(false, Ordering::SeqCst);
            self.hash_done_chunks(u, &mut h)?;
            drop(h);
            if !self.unit_kick[u].load(Ordering::SeqCst) {
                return Ok(());
            }
        }
    }

    /// Hash chunks of unit `u` from the cursor while they are on the sink;
    /// record the digest when the unit is complete.
    fn hash_done_chunks(&self, u: usize, h: &mut UnitHash) -> Result<(), Error> {
        let unit = self.layout.units[u];
        let end = unit.first_chunk + unit.chunks;
        let mut buf = Vec::new();
        while h.next < end && self.state[h.next].load(Ordering::SeqCst) == DONE {
            let spec = self.layout.chunks[h.next];
            buf.resize(spec.len as usize, 0);
            self.sink.read_at(spec.offset, &mut buf)?;
            h.hasher.update(&buf);
            h.next += 1;
        }
        self.finish_unit_hash(u, h);
        Ok(())
    }

    fn finish_unit_hash(&self, u: usize, h: &mut UnitHash) {
        let unit = self.layout.units[u];
        if h.next == unit.first_chunk + unit.chunks {
            let mut slot = lock(&self.unit_md5[u]);
            if slot.is_none() {
                *slot = Some(std::mem::take(&mut h.hasher).finalize().into());
            }
        }
    }

    /// MD5 of a range read back from the sink with positioned reads (handles
    /// that share a cursor would interleave under concurrent hashing).
    fn md5_of_range(&self, offset: u64, len: u64) -> Result<[u8; 16], Error> {
        let mut h = Md5::new();
        let mut buf = vec![0u8; (len.min(4 << 20)) as usize];
        let mut at = offset;
        while at < offset + len {
            let n = ((offset + len - at) as usize).min(buf.len());
            self.sink.read_at(at, &mut buf[..n])?;
            h.update(&buf[..n]);
            at += n as u64;
        }
        Ok(h.finalize().into())
    }

    fn flush_ledger(&self, force: bool) {
        if self.seq.is_some() || !(self.dirty.swap(false, Ordering::SeqCst) || force) {
            return;
        }
        let state = ledger::encode(
            self.object,
            &self.layout,
            self.state.iter().map(|s| s.load(Ordering::SeqCst) == DONE),
            self.blake3.iter().map(|b| *lock(b)),
        );
        let _ = self.sink.checkpoint(&state);
    }

    fn digests(&self) -> Result<Vec<[u8; 16]>, Error> {
        let mut out = Vec::with_capacity(self.layout.units.len());
        for (u, unit) in self.layout.units.iter().enumerate() {
            let have = *lock(&self.unit_md5[u]);
            let d = match (have, &self.seq) {
                (Some(d), _) => d,
                (None, None) => {
                    // Finish the cursor; every chunk is on the sink by now.
                    self.hash_done_chunks(u, &mut lock(&self.unit_hash[u]))?;
                    match *lock(&self.unit_md5[u]) {
                        Some(d) => d,
                        None => {
                            let d = self.md5_of_range(unit.offset, unit.len)?;
                            *lock(&self.unit_md5[u]) = Some(d);
                            d
                        }
                    }
                }
                (None, Some(_)) => {
                    return Err(Error::new(
                        ErrorKind::Internal,
                        "stream unit was not hashed",
                    ));
                }
            };
            out.push(d);
        }
        Ok(out)
    }

    fn report(&self, t0: Instant, repaired: u64) -> Report {
        Report {
            size: self.object.size,
            reused: false,
            resumed_chunks: self.resumed,
            elapsed: t0.elapsed(),
            chunk_retries: self.retries.load(Ordering::SeqCst),
            hedged: self.hedged.load(Ordering::SeqCst),
            repaired_chunks: repaired,
        }
    }

    fn finish(self, t0: Instant) -> Result<Download<S::Output>, Error> {
        let failed = lock(&self.failed).clone();
        if let Some(e) = failed {
            return Err(self.abort(e));
        }
        let mut repaired = 0;
        if self.verifying {
            let first = match self.digests() {
                Ok(d) => object::verify(self.object, &d),
                Err(e) => return Err(self.abort(e)),
            };
            if let Err(first) = first {
                if self.seq.is_some() {
                    return Err(self.abort(first));
                }
                // Localize corruption: re-fetch every chunk, rewrite those
                // whose BLAKE3 differs, re-hash affected units.
                self.emit(|| Event::Repairing);
                repaired = match self.repair() {
                    Ok(n) => n,
                    Err(e) => {
                        let e = Error::new(e.kind(), format!("{first}; repair failed: {e}"))
                            .with_source(e);
                        return Err(self.abort(e));
                    }
                };
                let second = self.digests().and_then(|d| object::verify(self.object, &d));
                if let Err(second) = second {
                    let e = Error::new(
                        second.kind(),
                        format!("{second} after repairing {repaired} chunk(s)"),
                    );
                    return Err(self.abort(e));
                }
            }
        }
        let report = self.report(t0, repaired);
        let output = self.sink.commit()?;
        Ok(Download { output, report })
    }

    fn abort(self, error: Error) -> Error {
        self.flush_ledger(true);
        self.sink.abort(&error);
        error
    }

    /// Re-fetch every chunk in parallel, compare BLAKE3 with what was stored
    /// and with what the sink holds, and rewrite only differing chunks; then
    /// re-hash the affected units. Returns the number of rewritten chunks.
    fn repair(&self) -> Result<u64, Error> {
        let n = self.layout.chunks.len();
        let next = AtomicUsize::new(0);
        let repaired = AtomicU64::new(0);
        let touched = Mutex::new(BTreeSet::new());
        let failure = Mutex::new(None::<Error>);
        let workers = self.opts().repair_connections.clamp(1, n.max(1));
        std::thread::scope(|scope| {
            for _ in 0..workers {
                scope.spawn(|| {
                    let mut cancel = self.control.cancel_tx.subscribe();
                    loop {
                        let ci = next.fetch_add(1, Ordering::SeqCst);
                        if ci >= n || lock(&failure).is_some() {
                            return;
                        }
                        if let Err(e) = self.repair_chunk(ci, &mut cancel, &repaired, &touched) {
                            lock(&failure).get_or_insert(e);
                        }
                    }
                });
            }
        });
        if let Some(e) = failure.into_inner().unwrap_or_else(|e| e.into_inner()) {
            return Err(e);
        }
        for u in touched.into_inner().unwrap_or_else(|e| e.into_inner()) {
            let unit = self.layout.units[u];
            let d = self.md5_of_range(unit.offset, unit.len)?;
            *lock(&self.unit_md5[u]) = Some(d);
        }
        self.dirty.store(true, Ordering::SeqCst);
        Ok(repaired.into_inner())
    }

    fn repair_chunk(
        &self,
        ci: usize,
        cancel: &mut watch::Receiver<bool>,
        repaired: &AtomicU64,
        touched: &Mutex<BTreeSet<usize>>,
    ) -> Result<(), Error> {
        let spec = self.layout.chunks[ci];
        let mut fetch = || {
            let mut nop = |_: u64| {};
            let fut = http::fetch_range(
                &self.shared.client,
                &self.object.url,
                spec.offset,
                spec.len,
                REPAIR_TIMEOUT,
                &mut nop,
            );
            self.shared
                .rt()
                .handle()
                .block_on(http::or_cancel(fut, cancel))
                .unwrap_or_else(|| Err(Error::cancelled()))
        };
        let fresh = match fetch() {
            Err(e) if e.kind() != ErrorKind::Cancelled => fetch(),
            r => r,
        }?;
        let digest = *blake3::hash(&fresh).as_bytes();
        let stored = *lock(&self.blake3[ci]);
        let mut on_disk = vec![0u8; spec.len as usize];
        self.sink.read_at(spec.offset, &mut on_disk)?;
        if stored != Some(digest) || *blake3::hash(&on_disk).as_bytes() != digest {
            self.sink.write_at(spec.offset, fresh)?;
            *lock(&self.blake3[ci]) = Some(digest);
            repaired.fetch_add(1, Ordering::SeqCst);
            lock(touched).insert(spec.unit);
        }
        Ok(())
    }
}
