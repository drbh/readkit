//! `.sra` → [`Batch`]es: options, segment selection, names, threads.

use std::collections::BTreeMap;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Condvar, Mutex, MutexGuard};
use std::thread::JoinHandle;

use readkit_batch::{Batch, BatchBuilder, Layout, Paired, Pool, Segments, Single};
use readkit_io::BatchSource;

use crate::backend::{OpenRequest, RunInfo, Spot, SpotReader, SraBackend};
use crate::error::{Error, ErrorKind};
use crate::lite::{LitePolicy, QualityKind, QualityProbe};
use crate::read_type::ReadType;

/// Which reads of each spot become batch segments.
///
/// SRA spots hold every read of a sequencing cluster: biological reads
/// (cDNA, genomic mates) and technical ones (cell barcodes, UMIs, sample
/// indexes, linkers). Each selected read becomes one segment slot, in spot
/// order. The selection must give the same number of segments for every
/// spot; otherwise decoding fails with [`ErrorKind::Layout`].
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum SegmentSelection {
    /// Biological reads only, as `fasterq-dump` does by default. The
    /// default: single-end and paired RNA-seq.
    #[default]
    Biological,
    /// Every read, technical ones included (`fasterq-dump
    /// --include-technical`). Single-cell needs this: the barcode and UMI
    /// reads are technical.
    All,
    /// These reads by 0-based index within the spot, in this order (for
    /// example `[1, 0]` for cDNA first, barcode second). Duplicates are
    /// allowed.
    Reads(Vec<usize>),
}

/// Whether and how to produce read names.
///
/// Names cost a formatted string per segment; quantifiers do not need them.
/// With names on, segment `j` of spot `s` is named like `fasterq-dump`
/// names it (without the leading `@`), where `<acc>` is
/// [`DecodeOptions::accession`] or the file name, `<name>` is the spot's
/// stored `NAME` or, when the run stores none, the spot id, and `<len>` is
/// the segment's length.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Names {
    /// No names (the batch has no name columns).
    #[default]
    None,
    /// `<acc>.<spot> <name> length=<len>`: `fasterq-dump`'s default defline
    /// (split-3 and `--split-files` modes).
    Fasterq,
    /// `<acc>.<spot>/<j+1> <name> length=<len>`: `fasterq-dump
    /// --split-spot`'s defline, `<j+1>` being the 1-based position among
    /// the selected segments.
    FasterqSplitSpot,
}

/// Options for [`Decoder`].
///
/// Construct with `DecodeOptions { threads: 2, ..Default::default() }`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DecodeOptions {
    /// Which reads become segments (default [`SegmentSelection::Biological`]).
    pub segments: SegmentSelection,
    /// Read names (default [`Names::None`]).
    pub names: Names,
    /// What to do with SRA Lite runs (default [`LitePolicy::Reject`]).
    pub lite: LitePolicy,
    /// Spots per batch (default 16 384; at least 1).
    pub batch_spots: u64,
    /// Decoding threads (default 4). With 1, or a backend without random
    /// access, decoding runs on the caller's thread inside
    /// [`next_batch`](BatchSource::next_batch). Otherwise worker threads
    /// decode spot ranges ahead, at most `2 × threads` batches beyond the
    /// one the consumer is waiting for, and batches come out in spot order.
    pub threads: usize,
    /// Accession used in names (default: the file name).
    pub accession: Option<String>,
    /// Decode only these spot ids (default: the whole run). Clamped to the
    /// run; useful for sharding a run across processes.
    pub spot_range: Option<std::ops::Range<u64>>,
}

impl Default for DecodeOptions {
    fn default() -> Self {
        Self {
            segments: SegmentSelection::Biological,
            names: Names::None,
            lite: LitePolicy::Reject,
            batch_spots: 16_384,
            threads: 4,
            accession: None,
            spot_range: None,
        }
    }
}

/// Counts of what a decoder has emitted so far.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct DecodeStats {
    /// Spots (records) emitted.
    pub spots: u64,
    /// Segments emitted (`spots × segments`, empty ones included).
    pub reads: u64,
    /// Bases emitted.
    pub bases: u64,
    /// Longest segment emitted.
    pub max_read_len: u64,
    /// Batches emitted.
    pub batches: u64,
}

mod sealed {
    pub trait Sealed {}
    impl Sealed for readkit_batch::Single {}
    impl Sealed for readkit_batch::Paired {}
    impl Sealed for readkit_batch::Segments {}
}

/// Batch layouts a [`Decoder`] can produce: [`Single`], [`Paired`] and
/// [`Segments`]. Sealed.
pub trait DecodeLayout: Layout + sealed::Sealed {
    #[doc(hidden)]
    fn builder(pool: &Pool, k: usize, records: usize) -> BatchBuilder<Self>
    where
        Self: Sized;
}

impl DecodeLayout for Single {
    fn builder(pool: &Pool, _k: usize, records: usize) -> BatchBuilder<Self> {
        BatchBuilder::with_pool(pool, records)
    }
}
impl DecodeLayout for Paired {
    fn builder(pool: &Pool, _k: usize, records: usize) -> BatchBuilder<Self> {
        BatchBuilder::with_pool(pool, records)
    }
}
impl DecodeLayout for Segments {
    fn builder(pool: &Pool, k: usize, records: usize) -> BatchBuilder<Self> {
        BatchBuilder::with_segments(pool, k, records)
    }
}

/// The default backend: VDB when built with feature `ngs`, otherwise
/// [`FasterqDumpBackend`](crate::FasterqDumpBackend) when `fasterq-dump`
/// and `vdb-dump` are on `PATH`.
///
/// Fails with [`ErrorKind::Unsupported`] when neither is available.
pub fn default_backend() -> Result<Arc<dyn SraBackend>, Error> {
    #[cfg(feature = "ngs")]
    {
        Ok(Arc::new(crate::ngs::NgsBackend::new()))
    }
    #[cfg(not(feature = "ngs"))]
    {
        let fasterq = crate::fasterq::FasterqDumpBackend::new();
        if fasterq.is_available() {
            return Ok(Arc::new(fasterq));
        }
        Err(Error::new(
            ErrorKind::Unsupported,
            "readkit-sra was built without the `ngs` feature and `fasterq-dump` \
             (with `vdb-dump`) is not on PATH",
        ))
    }
}

/// Open `path` with the default backend and return its [`RunInfo`] (spot
/// count, quality kind) without decoding.
pub fn run_info(path: impl AsRef<Path>) -> Result<RunInfo, Error> {
    let path = path.as_ref();
    check_path(path)?;
    let reader = default_backend()?.open(path, &OpenRequest::new(false))?;
    Ok(reader.info().clone())
}

fn check_path(path: &Path) -> Result<(), Error> {
    match std::fs::metadata(path) {
        Ok(_) => Ok(()),
        Err(e) => Err(Error::from(e).with_input(path.display().to_string())),
    }
}

/// What every chunk decode needs; shared by workers.
struct Plan {
    path: PathBuf,
    accession: String,
    names: Names,
    selection: SegmentSelection,
    k: usize,
    first: u64,
    end: u64,
    batch_spots: u64,
}

impl Plan {
    fn chunks(&self) -> u64 {
        (self.end - self.first).div_ceil(self.batch_spots)
    }

    fn range(&self, chunk: u64) -> (u64, u64) {
        let start = self.first + chunk * self.batch_spots;
        (start, (start + self.batch_spots).min(self.end))
    }

    fn error(&self, kind: ErrorKind, msg: impl Into<String>, spot: u64) -> Error {
        Error::new(kind, msg)
            .with_input(self.path.display().to_string())
            .at_spot(spot)
    }

    /// Validate `s` and write the selected read indices into `sel`.
    fn select(&self, s: &Spot<'_>, spot: u64, sel: &mut Vec<usize>) -> Result<(), Error> {
        let n = s.types.len();
        if s.starts.len() != n || s.lens.len() != n {
            return Err(self.error(
                ErrorKind::Corrupt,
                format!(
                    "{} read types, {} starts, {} lengths",
                    n,
                    s.starts.len(),
                    s.lens.len()
                ),
                spot,
            ));
        }
        if s.bases.len() != s.quals.len() {
            return Err(self.error(
                ErrorKind::Corrupt,
                format!(
                    "sequence {} bytes, quality {} bytes",
                    s.bases.len(),
                    s.quals.len()
                ),
                spot,
            ));
        }
        for r in 0..n {
            let ok = usize::try_from(s.starts[r])
                .ok()
                .and_then(|st| st.checked_add(s.lens[r] as usize))
                .is_some_and(|end| end <= s.bases.len());
            if !ok {
                return Err(self.error(
                    ErrorKind::Corrupt,
                    format!(
                        "read {r} at {}+{} outside a {}-base spot",
                        s.starts[r],
                        s.lens[r],
                        s.bases.len()
                    ),
                    spot,
                ));
            }
        }
        sel.clear();
        match &self.selection {
            SegmentSelection::Biological => {
                sel.extend((0..n).filter(|&r| s.types[r].is_biological()));
            }
            SegmentSelection::All => sel.extend(0..n),
            SegmentSelection::Reads(idx) => {
                if let Some(bad) = idx.iter().find(|&&r| r >= n) {
                    return Err(self.error(
                        ErrorKind::Layout,
                        format!("selected read {bad} but the spot has {n} reads"),
                        spot,
                    ));
                }
                sel.extend_from_slice(idx);
            }
        }
        if sel.len() != self.k {
            return Err(self.error(
                ErrorKind::Layout,
                format!(
                    "spot gives {} segments, earlier spots gave {}",
                    sel.len(),
                    self.k
                ),
                spot,
            ));
        }
        Ok(())
    }

    /// Decode spots `[start, end)` into one batch.
    fn decode<L: DecodeLayout>(
        &self,
        reader: &mut dyn SpotReader,
        pool: &Pool,
        (start, end): (u64, u64),
    ) -> Result<Batch<L>, Error> {
        let mut builder = L::builder(pool, self.k, (end - start) as usize);
        let mut sel = Vec::with_capacity(self.k);
        let mut lens = vec![0usize; self.k];
        let mut names: Vec<Vec<u8>> = vec![Vec::new(); self.k];
        for spot in start..end {
            let s = reader.read_spot(spot)?;
            self.select(&s, spot, &mut sel)?;
            for (j, &r) in sel.iter().enumerate() {
                lens[j] = s.lens[r] as usize;
            }
            let named = self.names != Names::None;
            if named {
                for (j, buf) in names.iter_mut().enumerate() {
                    buf.clear();
                    let _ = write!(buf, "{}.{spot}", self.accession);
                    if self.names == Names::FasterqSplitSpot {
                        let _ = write!(buf, "/{}", j + 1);
                    }
                    buf.push(b' ');
                    match s.name {
                        Some(n) if !n.is_empty() => buf.extend_from_slice(n),
                        _ => {
                            let _ = write!(buf, "{spot}");
                        }
                    }
                    let _ = write!(buf, " length={}", lens[j]);
                }
            }
            let refs: Vec<&[u8]>;
            let name_arg = if named {
                refs = names.iter().map(Vec::as_slice).collect();
                Some(refs.as_slice())
            } else {
                None
            };
            let (bases, quals) = (s.bases, s.quals);
            builder
                .push_with(&lens, name_arg, |j, seq, qual| {
                    let r = sel[j];
                    let st = s.starts[r] as usize;
                    seq.copy_from_slice(&bases[st..st + seq.len()]);
                    qual.copy_from_slice(&quals[st..st + qual.len()]);
                })
                .map_err(|e| {
                    Error::from(e)
                        .with_input(self.path.display().to_string())
                        .at_spot(spot)
                })?;
        }
        Ok(builder.finish())
    }
}

struct State<L> {
    ready: BTreeMap<u64, Result<Batch<L>, Error>>,
    next_claim: u64,
    next_emit: u64,
    stop: bool,
    fatal: Option<Error>,
}

struct Shared<L> {
    state: Mutex<State<L>>,
    cv: Condvar,
}

fn lock<T>(m: &Mutex<T>) -> MutexGuard<'_, T> {
    m.lock().unwrap_or_else(|e| e.into_inner())
}

/// Reports a panicking worker so the consumer does not wait forever.
struct PanicGuard<L>(Arc<Shared<L>>);

impl<L> Drop for PanicGuard<L> {
    fn drop(&mut self) {
        if std::thread::panicking() {
            let mut st = lock(&self.0.state);
            st.fatal
                .get_or_insert_with(|| Error::new(ErrorKind::Internal, "decode worker panicked"));
            drop(st);
            self.0.cv.notify_all();
        }
    }
}

enum Engine<L> {
    /// Not started yet; holds the reader opened by `open`.
    Idle(Box<dyn SpotReader>),
    /// Decoding on the caller's thread.
    Inline {
        reader: Box<dyn SpotReader>,
        next: u64,
    },
    /// Worker threads.
    Threads {
        shared: Arc<Shared<L>>,
        workers: Vec<JoinHandle<()>>,
    },
    Done,
}

/// Decodes a local `.sra` run into [`Batch<L>`]s; a
/// [`BatchSource`](readkit_io::BatchSource).
///
/// Bases and qualities are copied from the backend straight into the
/// batch's sequence and quality columns ([`BatchBuilder::push_with`]): no
/// FASTQ text is produced. `L` is [`Single`] (one segment per spot after
/// [`SegmentSelection`]), [`Paired`] (two) or [`Segments`] (any number,
/// for single-cell). Opening checks the first spot's read structure against
/// `L`, and every later spot must have the same number of selected reads.
///
/// ```no_run
/// use readkit_io::BatchSource;
/// use readkit_sra::{DecodeOptions, PairedDecoder};
///
/// let mut decoder = PairedDecoder::open("SRR13638690/SRR13638690", DecodeOptions::default())?;
/// println!("{} spots", decoder.info().spots);
/// while let Some(batch) = decoder.next_batch()? {
///     for (r1, r2) in batch.pairs() {
///         let _ = (r1.seq(), r2.qual());
///     }
/// }
/// # Ok::<(), readkit_sra::Error>(())
/// ```
pub struct Decoder<L: DecodeLayout> {
    backend: Arc<dyn SraBackend>,
    plan: Arc<Plan>,
    info: RunInfo,
    read_types: Vec<ReadType>,
    pool: Pool,
    threads: usize,
    lite: LitePolicy,
    probe: Option<QualityProbe>,
    stats: DecodeStats,
    engine: Engine<L>,
}

/// A [`Decoder`] for single-end runs.
pub type SingleDecoder = Decoder<Single>;
/// A [`Decoder`] for paired runs.
pub type PairedDecoder = Decoder<Paired>;
/// A [`Decoder`] with a runtime segment count (single-cell).
pub type SegmentsDecoder = Decoder<Segments>;

impl<L: DecodeLayout> Decoder<L> {
    /// Open the run at `path` with the [`default_backend`].
    pub fn open(path: impl AsRef<Path>, options: DecodeOptions) -> Result<Self, Error> {
        Self::open_with(default_backend()?, path, options)
    }

    /// Open the run at `path` with `backend`.
    ///
    /// Fails with [`ErrorKind::NotFound`] if the file does not exist,
    /// [`ErrorKind::Lite`] if the run declares synthetic qualities and
    /// `options.lite` is [`LitePolicy::Reject`], and [`ErrorKind::Layout`]
    /// if the first spot's selected reads do not fit `L`.
    pub fn open_with(
        backend: Arc<dyn SraBackend>,
        path: impl AsRef<Path>,
        options: DecodeOptions,
    ) -> Result<Self, Error> {
        let path = path.as_ref().to_path_buf();
        let input = path.display().to_string();
        check_path(&path)?;
        let request = OpenRequest::new(options.names != Names::None);
        let mut reader = backend.open(&path, &request)?;
        let info = reader.info().clone();
        if info.quality == QualityKind::Synthetic && options.lite == LitePolicy::Reject {
            return Err(Error::new(
                ErrorKind::Lite,
                "the run declares synthetic qualities only (SRA Lite)",
            )
            .with_input(input));
        }
        let accession = options.accession.clone().unwrap_or_else(|| {
            path.file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_default()
        });
        let (first, end) = match &options.spot_range {
            Some(r) => {
                let first = r.start.clamp(info.first_spot, info.end_spot());
                (first, r.end.clamp(first, info.end_spot()))
            }
            None => (info.first_spot, info.end_spot()),
        };
        let mut plan = Plan {
            path,
            accession,
            names: options.names,
            selection: options.segments.clone(),
            k: 0,
            first,
            end,
            batch_spots: options.batch_spots.max(1),
        };
        // The first spot decides the segment count.
        let mut read_types = Vec::new();
        if end > first {
            let s = reader.read_spot(first)?;
            plan.k = match &plan.selection {
                SegmentSelection::Biological => {
                    s.types.iter().filter(|t| t.is_biological()).count()
                }
                SegmentSelection::All => s.types.len(),
                SegmentSelection::Reads(r) => r.len(),
            };
            let mut sel = Vec::new();
            plan.select(&s, first, &mut sel)?;
            read_types = sel.iter().map(|&r| s.types[r]).collect();
        } else {
            plan.k = L::SEGMENTS.unwrap_or(1);
        }
        if plan.k == 0 {
            return Err(Error::new(
                ErrorKind::Layout,
                format!(
                    "the selection gives no segments ({} reads in the first spot)",
                    read_types.len()
                ),
            )
            .with_input(input));
        }
        if let Some(want) = L::SEGMENTS
            && want != plan.k
        {
            return Err(Error::new(
                ErrorKind::Layout,
                format!(
                    "the run has {} selected segments per spot, the layout needs {want}",
                    plan.k
                ),
            )
            .with_input(input));
        }
        let probe = (info.quality == QualityKind::Unknown && options.lite == LitePolicy::Reject)
            .then(QualityProbe::new);
        let threads = if reader.random_access() {
            options.threads.max(1)
        } else {
            1
        };
        Ok(Self {
            backend,
            plan: Arc::new(plan),
            info,
            read_types,
            pool: Pool::new(),
            threads,
            lite: options.lite,
            probe,
            stats: DecodeStats::default(),
            engine: Engine::Idle(reader),
        })
    }

    /// Use `pool` for batch buffers instead of the decoder's own. Takes
    /// effect only before the first batch.
    pub fn with_pool(mut self, pool: &Pool) -> Self {
        if matches!(self.engine, Engine::Idle(_)) {
            self.pool = pool.clone();
        }
        self
    }

    /// The pool batches are drawn from.
    pub fn pool(&self) -> &Pool {
        &self.pool
    }

    /// Facts about the run (spot count, quality kind).
    pub fn info(&self) -> &RunInfo {
        &self.info
    }

    /// Segments per record.
    pub fn segments(&self) -> usize {
        self.plan.k
    }

    /// The spot ids this decoder emits.
    pub fn spot_range(&self) -> std::ops::Range<u64> {
        self.plan.first..self.plan.end
    }

    /// The read types of the selected segments, from the first spot.
    pub fn read_types(&self) -> &[ReadType] {
        &self.read_types
    }

    /// The run's quality kind. [`QualityKind::Unknown`] runs become
    /// `Original` or `Synthetic` once sampling has decided.
    pub fn quality(&self) -> QualityKind {
        self.info.quality
    }

    /// The Lite policy in force.
    pub fn lite_policy(&self) -> LitePolicy {
        self.lite
    }

    /// The backend's name.
    pub fn backend_name(&self) -> &'static str {
        self.backend.name()
    }

    /// What has been emitted so far.
    pub fn stats(&self) -> DecodeStats {
        self.stats
    }

    fn start(&mut self) {
        let Engine::Idle(reader) = std::mem::replace(&mut self.engine, Engine::Done) else {
            return;
        };
        if self.threads <= 1 || self.plan.chunks() <= 1 {
            self.engine = Engine::Inline { reader, next: 0 };
            return;
        }
        let shared = Arc::new(Shared {
            state: Mutex::new(State {
                ready: BTreeMap::new(),
                next_claim: 0,
                next_emit: 0,
                stop: false,
                fatal: None,
            }),
            cv: Condvar::new(),
        });
        let workers_n = self.threads.min(self.plan.chunks() as usize);
        let window = 2 * workers_n as u64;
        let mut first_reader = Some(reader);
        let mut workers = Vec::with_capacity(workers_n);
        for i in 0..workers_n {
            let reader = first_reader.take();
            let (shared2, plan, pool, backend) = (
                shared.clone(),
                self.plan.clone(),
                self.pool.clone(),
                self.backend.clone(),
            );
            let names = plan.names != Names::None;
            let spawned = std::thread::Builder::new()
                .name(format!("readkit-sra-{i}"))
                .spawn(move || {
                    let _guard = PanicGuard(shared2.clone());
                    let mut reader = match reader {
                        Some(r) => r,
                        None => match backend.open(&plan.path, &OpenRequest::new(names)) {
                            Ok(r) => r,
                            Err(e) => {
                                lock(&shared2.state).fatal.get_or_insert(e);
                                shared2.cv.notify_all();
                                return;
                            }
                        },
                    };
                    worker(&shared2, &plan, &pool, reader.as_mut(), window);
                });
            match spawned {
                Ok(h) => workers.push(h),
                Err(e) => {
                    lock(&shared.state).fatal.get_or_insert(
                        Error::from(e).with_input(self.plan.path.display().to_string()),
                    );
                    break;
                }
            }
        }
        self.engine = Engine::Threads { shared, workers };
    }

    fn next_raw(&mut self) -> Result<Option<Batch<L>>, Error> {
        if matches!(self.engine, Engine::Idle(_)) {
            self.start();
        }
        let chunks = self.plan.chunks();
        match &mut self.engine {
            Engine::Idle(_) | Engine::Done => Ok(None),
            Engine::Inline { reader, next } => {
                if *next >= chunks {
                    return Ok(None);
                }
                let range = self.plan.range(*next);
                *next += 1;
                self.plan
                    .decode::<L>(reader.as_mut(), &self.pool, range)
                    .map(Some)
            }
            Engine::Threads { shared, .. } => {
                let mut st = lock(&shared.state);
                loop {
                    let k = st.next_emit;
                    if let Some(r) = st.ready.remove(&k) {
                        st.next_emit += 1;
                        drop(st);
                        shared.cv.notify_all();
                        return r.map(Some);
                    }
                    if let Some(e) = st.fatal.take() {
                        return Err(e);
                    }
                    if k >= chunks {
                        return Ok(None);
                    }
                    st = shared.cv.wait(st).unwrap_or_else(|e| e.into_inner());
                }
            }
        }
    }

    fn shutdown(&mut self) {
        if let Engine::Threads { shared, workers } =
            std::mem::replace(&mut self.engine, Engine::Done)
        {
            lock(&shared.state).stop = true;
            shared.cv.notify_all();
            for w in workers {
                let _ = w.join();
            }
        } else {
            self.engine = Engine::Done;
        }
    }

    fn account(&mut self, batch: &Batch<L>) -> Result<(), Error> {
        self.stats.batches += 1;
        self.stats.spots += batch.len() as u64;
        for col in batch.columns() {
            let seqs = col.seqs();
            self.stats.reads += seqs.len() as u64;
            self.stats.bases += seqs.data().len() as u64;
            for s in seqs.iter() {
                self.stats.max_read_len = self.stats.max_read_len.max(s.len() as u64);
            }
        }
        if let Some(probe) = &mut self.probe {
            'outer: for i in 0..batch.len() {
                for col in batch.columns() {
                    if probe.is_full() {
                        break 'outer;
                    }
                    probe.observe(col.quals().get(i));
                }
            }
            let finished = self.stats.spots >= self.plan.end - self.plan.first;
            if probe.is_full() || (finished && probe.sampled() > 0) {
                let synthetic = probe.looks_synthetic();
                let distinct = probe.distinct_values();
                self.probe = None;
                if synthetic {
                    self.info.quality = QualityKind::Synthetic;
                    return Err(Error::new(
                        ErrorKind::Lite,
                        format!(
                            "sampled qualities look synthetic ({distinct} distinct values; SRA Lite)"
                        ),
                    )
                    .with_input(self.plan.path.display().to_string()));
                }
                self.info.quality = QualityKind::Original;
            }
        }
        Ok(())
    }
}

fn worker<L: DecodeLayout>(
    shared: &Shared<L>,
    plan: &Plan,
    pool: &Pool,
    reader: &mut dyn SpotReader,
    window: u64,
) {
    let chunks = plan.chunks();
    loop {
        let k = {
            let mut st = lock(&shared.state);
            loop {
                if st.stop || st.next_claim >= chunks {
                    return;
                }
                if st.next_claim < st.next_emit + window {
                    let k = st.next_claim;
                    st.next_claim += 1;
                    break k;
                }
                st = shared.cv.wait(st).unwrap_or_else(|e| e.into_inner());
            }
        };
        let result = plan.decode::<L>(reader, pool, plan.range(k));
        let failed = result.is_err();
        let mut st = lock(&shared.state);
        if failed {
            // Nothing after a failed chunk will be emitted.
            st.next_claim = chunks;
        }
        st.ready.insert(k, result);
        drop(st);
        shared.cv.notify_all();
        if failed {
            return;
        }
    }
}

impl<L: DecodeLayout> BatchSource<L> for Decoder<L> {
    type Error = Error;

    /// The next batch in spot order. After an error or the end, returns
    /// `Ok(None)`.
    fn next_batch(&mut self) -> Result<Option<Batch<L>>, Error> {
        let result = self.next_raw().and_then(|b| match b {
            Some(batch) => self.account(&batch).map(|()| Some(batch)),
            None => Ok(None),
        });
        match &result {
            Ok(Some(_)) => {}
            Ok(None) => {
                let expected = self.plan.end - self.plan.first;
                if !matches!(self.engine, Engine::Done) && self.stats.spots != expected {
                    self.shutdown();
                    return Err(Error::new(
                        ErrorKind::Corrupt,
                        format!("decoded {} of {expected} spots", self.stats.spots),
                    )
                    .with_input(self.plan.path.display().to_string()));
                }
                self.shutdown();
            }
            Err(_) => self.shutdown(),
        }
        result
    }
}

impl<L: DecodeLayout> Drop for Decoder<L> {
    fn drop(&mut self) {
        self.shutdown();
    }
}

impl<L: DecodeLayout> std::fmt::Debug for Decoder<L> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Decoder")
            .field("backend", &self.backend.name())
            .field("path", &self.plan.path)
            .field("info", &self.info)
            .field("segments", &self.plan.k)
            .field("threads", &self.threads)
            .field("stats", &self.stats)
            .finish()
    }
}
