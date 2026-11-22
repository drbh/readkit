//! Parallel parsing ([`FastqOptions::threads`](crate::FastqOptions::threads)).
//!
//! Three kinds of thread:
//!
//! - a scanner per input reads it into large buffers and finds every
//!   record's end (every fourth `\n`, the record starting with `@`), then
//!   ships each buffer cut after its last whole record ("pieces");
//! - one splitter takes record ends from the pieces spot by spot and cuts
//!   batches under the same byte budget and record limit as the sequential
//!   path; a batch's text is a list of ranges in pieces, so nothing is copied;
//! - workers validate the records and copy them into columns, exactly as the
//!   sequential fast path does.
//!
//! Batches and errors are delivered in input order.
//!
//! Only plain text is cut this way: records that start with `@` and whose
//! four lines end in `\n`. At anything else (a blank line, a record that does
//! not start with `@`) and at the end of any input, the splitter stops and
//! hands the inputs back with the text not yet in a batch (from the start of
//! the batch it was cutting) in front of them, and the sequential reader
//! takes over from that point. So blank lines, the last batch, truncation and
//! record-count errors all go through the code that already handles them.

use std::collections::{BTreeMap, VecDeque};
use std::io::{BufRead, Cursor, Read};
use std::ops::Range;
use std::panic::{AssertUnwindSafe, catch_unwind, resume_unwind};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{Receiver, Sender, SyncSender, channel, sync_channel};
use std::sync::{Arc, Mutex};
use std::thread::JoinHandle;

use readkit_batch::{Batch, BatchBuilder, Pool, Segments};

use crate::error::{Error, ErrorKind};
use crate::source::{Lines, batch_error, check, push_spot, split_record};

type Panic = Box<dyn std::any::Any + Send + 'static>;
type Reader = Box<dyn BufRead + Send>;

/// Bytes a scanner reads into one piece (more for a record that does not
/// fit).
const PIECE: usize = if cfg!(test) { 64 } else { 4 << 20 };
/// Pieces a scanner may have shipped that the splitter has not taken.
const PIECES_AHEAD: usize = 4;

/// What the parse needs besides the text.
#[derive(Clone)]
pub(crate) struct Config {
    pub k: usize,
    pub interleaved: bool,
    pub phred64: bool,
    pub keep_names: bool,
    pub batch_bytes: usize,
    pub batch_records: usize,
    pub pool: Pool,
    pub names: Arc<[Option<Arc<str>>]>,
}

impl Config {
    fn reader_of(&self, s: usize) -> usize {
        if self.interleaved { 0 } else { s }
    }

    fn error(&self, kind: ErrorKind, record: u64, s: usize) -> Error {
        Error::new(kind, "").at(
            Some(record),
            (self.k > 1).then_some(s),
            self.names[self.reader_of(s)].as_ref(),
        )
    }
}

/// The inputs handed back to the sequential reader.
pub(crate) struct Tail {
    /// Each input with the text read ahead of it chained in front.
    pub readers: Vec<Reader>,
    /// FASTQ records consumed per input before that text.
    pub records: Vec<u64>,
}

/// A read error a scanner met, replayed where it happened: after the text
/// read before it, so the sequential reader reports it as its own.
struct Failed(Option<std::io::Error>);

impl Read for Failed {
    fn read(&mut self, _: &mut [u8]) -> std::io::Result<usize> {
        self.fill_buf().map(|_| 0)
    }
}

impl BufRead for Failed {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        match self.0.take() {
            Some(e) => Err(e),
            None => Ok(&[]),
        }
    }

    fn consume(&mut self, _: usize) {}
}

/// A buffer of whole records: `text[..ends.last()]`, records ending at
/// `ends`.
struct Piece {
    text: Vec<u8>,
    ends: Vec<usize>,
}

/// What a scanner sends: pieces, then once its input stops being plain
/// text (or it is asked to stop), the rest.
enum Feed {
    Piece(Piece),
    End {
        /// Text read after the last shipped record.
        rest: Vec<u8>,
        reader: Reader,
        error: Option<std::io::Error>,
    },
}

/// One batch: per input, the ranges of pieces holding its records.
struct Job {
    seq: u64,
    text: Vec<Vec<(Arc<Piece>, Range<usize>)>>,
    spots: usize,
    /// Records of each input before this batch.
    records: Vec<u64>,
}

enum Msg {
    Batch(u64, Result<Batch<Segments>, Error>),
    Tail(u64, Tail),
    Panic(Panic),
}

/// The running threads, seen from the consumer.
pub(crate) struct Par {
    results: Option<Receiver<Msg>>,
    credit: Option<Receiver<()>>,
    pending: BTreeMap<u64, Result<Batch<Segments>, Error>>,
    tail: Option<(u64, Tail)>,
    next: u64,
    stop: Arc<AtomicBool>,
    threads: Vec<JoinHandle<()>>,
}

pub(crate) enum Next {
    Batch(Result<Batch<Segments>, Error>),
    /// Every batch cut in parallel has been delivered: continue sequentially
    /// from here.
    Sequential(Tail),
}

fn spawn(name: String, f: impl FnOnce() + Send + 'static) -> Option<JoinHandle<()>> {
    std::thread::Builder::new().name(name).spawn(f).ok()
}

impl Par {
    /// Starts the scanners, the splitter and `workers` workers over
    /// `readers`. Gives the readers back if the OS refuses a thread.
    pub fn start(
        readers: Vec<Reader>,
        records: Vec<u64>,
        cfg: Config,
        workers: usize,
    ) -> Result<Par, Vec<Reader>> {
        let workers = workers.max(1);
        let n = readers.len();
        let stop = Arc::new(AtomicBool::new(false));
        let (results_tx, results) = channel::<Msg>();
        let (jobs_tx, jobs_rx) = channel::<Job>();
        // At most this many batches cut and not yet delivered.
        let (credit_tx, credit) = sync_channel::<()>(workers + 2);
        let jobs_rx = Arc::new(Mutex::new(jobs_rx));
        let mut par = Par {
            results: Some(results),
            credit: Some(credit),
            pending: BTreeMap::new(),
            tail: None,
            next: 0,
            stop: stop.clone(),
            threads: Vec::new(),
        };
        // Empty buffers go back to the scanner of their input.
        let (spare_tx, spare_rx): (Vec<_>, Vec<_>) = (0..n).map(|_| channel::<Vec<u8>>()).unzip();
        for i in 0..workers {
            let (rx, tx, spare) = (jobs_rx.clone(), results_tx.clone(), spare_tx.clone());
            let cfg = cfg.clone();
            match spawn(format!("readkit-fastq-parse-{i}"), move || {
                worker(rx, tx, spare, cfg)
            }) {
                Some(h) => par.threads.push(h),
                None if i > 0 => break,
                None => return Err(readers),
            }
        }
        drop(jobs_rx);
        // Readers travel over channels so a refused spawn does not lose them.
        let mut feeds = Vec::with_capacity(n);
        let mut handoff = Vec::with_capacity(n);
        for (r, spare) in spare_rx.into_iter().enumerate() {
            let (feed_tx, feed_rx) = sync_channel::<Feed>(PIECES_AHEAD);
            let (reader_tx, reader_rx) = channel::<Reader>();
            let stop = stop.clone();
            let results = results_tx.clone();
            let h = spawn(format!("readkit-fastq-scan-{r}"), move || {
                let Ok(reader) = reader_rx.recv() else { return };
                if let Err(p) = catch_unwind(AssertUnwindSafe(|| {
                    scan_input(reader, &feed_tx, &spare, &stop)
                })) {
                    let _ = results.send(Msg::Panic(p));
                }
            });
            match h {
                Some(h) => par.threads.push(h),
                None => {
                    // Workers wait on the job channel: close it before the
                    // drop of `par` joins them.
                    drop((jobs_tx, handoff, par));
                    return Err(readers);
                }
            }
            feeds.push(feed_rx);
            handoff.push(reader_tx);
        }
        let splitter = Splitter {
            lanes: feeds.into_iter().map(Lane::new).collect(),
            records,
            cfg,
            jobs: jobs_tx,
            results: results_tx.clone(),
            credit: credit_tx,
            spare: spare_tx,
            stop,
            seq: 0,
        };
        let h = spawn("readkit-fastq-split".into(), move || {
            let mut s = splitter;
            if let Err(p) = catch_unwind(AssertUnwindSafe(|| s.run())) {
                let _ = results_tx.send(Msg::Panic(p));
            }
        });
        match h {
            Some(h) => par.threads.push(h),
            None => return Err(readers),
        }
        for (tx, reader) in handoff.into_iter().zip(readers) {
            tx.send(reader).expect("scanners wait for their inputs");
        }
        Ok(par)
    }

    /// The next batch in input order, or the hand-back to the sequential
    /// reader.
    pub fn next(&mut self) -> Next {
        loop {
            if let Some(r) = self.pending.remove(&self.next) {
                self.next += 1;
                if let Some(c) = &self.credit {
                    let _ = c.recv();
                }
                return Next::Batch(r);
            }
            if self.tail.as_ref().is_some_and(|(n, _)| *n == self.next) {
                let (_, t) = self.tail.take().unwrap();
                return Next::Sequential(t);
            }
            match self.results.as_ref().map(|r| r.recv()) {
                Some(Ok(Msg::Batch(seq, r))) => {
                    self.pending.insert(seq, r);
                }
                Some(Ok(Msg::Tail(n, t))) => self.tail = Some((n, t)),
                Some(Ok(Msg::Panic(p))) => resume_unwind(p),
                _ => {
                    return Next::Batch(Err(Error::new(
                        ErrorKind::Io,
                        "a parser thread ended without its batches",
                    )));
                }
            }
        }
    }
}

impl Drop for Par {
    fn drop(&mut self) {
        // Closing the result and credit channels stops the splitter at its
        // next send; scanners stop at their next read or send, workers when
        // the job channel closes. The join waits at most for a read in
        // progress.
        self.stop.store(true, Ordering::Relaxed);
        self.results.take();
        self.credit.take();
        self.pending.clear();
        self.tail.take();
        for h in self.threads.drain(..) {
            let _ = h.join();
        }
    }
}

/// The high bit of byte `i` set where byte `i` (little-endian) of `word` is
/// `\n`, and no other bit. Eight bytes per step without `unsafe`: about 2.5
/// times `memchr` per line on short-read FASTQ, where lines are ~40 bytes.
#[inline]
fn newlines(word: u64) -> u64 {
    const LOW7: u64 = 0x7f7f_7f7f_7f7f_7f7f;
    let x = word ^ 0x0a0a_0a0a_0a0a_0a0a;
    !(((x & LOW7) + LOW7) | x | LOW7)
}

/// Where a scanner is inside its buffer.
#[derive(Default)]
struct Scan {
    /// Bytes searched for line ends.
    scanned: usize,
    /// Start of the record being scanned and its line ends so far.
    start: usize,
    lines: u8,
    /// A record start that is not `@`: no record ends past it.
    odd: bool,
}

impl Scan {
    /// Finds record ends in `text[self.scanned..]`: every fourth line end,
    /// where the record's first byte is `@` (what `split_record` accepts).
    fn run(&mut self, text: &[u8], ends: &mut Vec<usize>) {
        let len = text.len();
        let mut base = self.scanned;
        self.scanned = len;
        while base < len {
            let mut mask = match text.get(base..base + 8) {
                Some(w) => newlines(u64::from_le_bytes(w.try_into().unwrap())),
                None => {
                    // The last few bytes; zero padding never matches `\n`.
                    let mut w = [0u8; 8];
                    w[..len - base].copy_from_slice(&text[base..]);
                    newlines(u64::from_le_bytes(w))
                }
            };
            while mask != 0 {
                let p = base + (mask.trailing_zeros() / 8) as usize;
                mask &= mask - 1;
                if self.lines == 0 && text[self.start] != b'@' {
                    self.odd = true;
                    return;
                }
                self.lines += 1;
                if self.lines == 4 {
                    self.lines = 0;
                    self.start = p + 1;
                    ends.push(self.start);
                }
            }
            base += 8;
        }
    }
}

/// A scanner: reads `reader` into pieces of whole records until its text
/// stops being plain, it ends, a read fails or `stop` is set, then sends the
/// rest.
fn scan_input(
    mut reader: Reader,
    feed: &SyncSender<Feed>,
    spare: &Receiver<Vec<u8>>,
    stop: &AtomicBool,
) {
    let fresh = |spare: &Receiver<Vec<u8>>| match spare.try_recv() {
        Ok(b) => b,
        Err(_) => vec![0; PIECE],
    };
    let mut buf = fresh(spare);
    let mut filled = 0;
    let mut ends = Vec::new();
    let mut scan = Scan::default();
    // Ships the whole records in `buf`; the partial one after them moves to
    // a new buffer.
    let ship =
        |buf: &mut Vec<u8>, filled: &mut usize, ends: &mut Vec<usize>, scan: &mut Scan| -> bool {
            let Some(&cut) = ends.last() else { return true };
            let mut next = fresh(spare);
            let rest = *filled - cut;
            if next.len() < rest.max(PIECE) {
                next.resize(rest.max(PIECE), 0);
            }
            next[..rest].copy_from_slice(&buf[cut..*filled]);
            let text = std::mem::replace(buf, next);
            *filled = rest;
            scan.scanned -= cut;
            scan.start -= cut;
            feed.send(Feed::Piece(Piece {
                text,
                ends: std::mem::take(ends),
            }))
            .is_ok()
        };
    let error = loop {
        if stop.load(Ordering::Relaxed) || scan.odd {
            break None;
        }
        if filled == buf.len() {
            if ends.is_empty() {
                // One record larger than the buffer.
                buf.resize(buf.len() * 2, 0);
            } else if !ship(&mut buf, &mut filled, &mut ends, &mut scan) {
                return;
            }
        }
        // `read` on a `BufReader` with nothing buffered reads straight into
        // `buf` when the request is at least its capacity: no copy.
        match reader.read(&mut buf[filled..]) {
            Ok(0) => break None,
            Ok(got) => {
                filled += got;
                scan.run(&buf[..filled], &mut ends);
            }
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
            Err(e) => break Some(e),
        }
    };
    if !ship(&mut buf, &mut filled, &mut ends, &mut scan) {
        return;
    }
    buf.truncate(filled);
    let _ = feed.send(Feed::End {
        rest: buf,
        reader,
        error,
    });
}

/// The splitter's view of one input: pieces from the one holding the start
/// of the batch being cut.
struct Lane {
    feed: Receiver<Feed>,
    pieces: VecDeque<Arc<Piece>>,
    /// Offset of the batch's start in `pieces[0]`.
    from: usize,
    /// The next record: `pieces[cur].ends[rec]`.
    cur: usize,
    rec: usize,
    end: Option<(Vec<u8>, Reader, Option<std::io::Error>)>,
}

impl Lane {
    fn new(feed: Receiver<Feed>) -> Self {
        Lane {
            feed,
            pieces: VecDeque::new(),
            from: 0,
            cur: 0,
            rec: 0,
            end: None,
        }
    }

    /// The next record's length; `None` where the parallel cut stops.
    fn next_record(&mut self) -> Result<Option<usize>, Stop> {
        loop {
            if let Some(p) = self.pieces.get(self.cur) {
                if let Some(&end) = p.ends.get(self.rec) {
                    let start = match self.rec {
                        0 => 0,
                        r => p.ends[r - 1],
                    };
                    self.rec += 1;
                    return Ok(Some(end - start));
                }
                self.cur += 1;
                self.rec = 0;
                continue;
            }
            if self.end.is_some() {
                return Ok(None);
            }
            match self.feed.recv() {
                Ok(Feed::Piece(p)) => self.pieces.push_back(Arc::new(p)),
                Ok(Feed::End {
                    rest,
                    reader,
                    error,
                }) => self.end = Some((rest, reader, error)),
                // The scanner panicked; it reported that.
                Err(_) => return Err(Stop::Closed),
            }
        }
    }

    /// The end of the last record taken: (piece, offset).
    fn taken(&self) -> (usize, usize) {
        match self.rec {
            0 => (self.cur, if self.cur == 0 { self.from } else { 0 }),
            r => (self.cur, self.pieces[self.cur].ends[r - 1]),
        }
    }

    /// Ranges from the batch start to the last record taken; that becomes
    /// the new batch start.
    fn cut(&mut self, spare: &Sender<Vec<u8>>) -> Vec<(Arc<Piece>, Range<usize>)> {
        let (last, end) = self.taken();
        let mut out = Vec::with_capacity(last + 1);
        for (i, p) in self.pieces.iter().enumerate().take(last + 1) {
            let a = if i == 0 { self.from } else { 0 };
            let b = if i == last {
                end
            } else {
                p.ends.last().copied().unwrap_or(0)
            };
            if a < b {
                out.push((p.clone(), a..b));
            }
        }
        for _ in 0..last {
            recycle(self.pieces.pop_front().unwrap(), spare);
        }
        self.cur -= last;
        self.from = end;
        out
    }

    /// Everything not yet in a batch, then the reader (or the read error the
    /// scanner met): drains the scanner, which `stop` (or its own end) ends.
    fn hand_back(&mut self) -> Reader {
        while self.end.is_none() {
            match self.feed.recv() {
                Ok(Feed::Piece(p)) => self.pieces.push_back(Arc::new(p)),
                Ok(Feed::End {
                    rest,
                    reader,
                    error,
                }) => self.end = Some((rest, reader, error)),
                Err(_) => unreachable!("a scanner always ends with End unless it panicked"),
            }
        }
        let (rest, reader, error) = self.end.take().unwrap();
        let mut text = Vec::new();
        for (i, p) in self.pieces.iter().enumerate() {
            let a = if i == 0 { self.from } else { 0 };
            let b = p.ends.last().copied().unwrap_or(0);
            text.extend_from_slice(&p.text[a.min(b)..b]);
        }
        text.extend_from_slice(&rest);
        self.pieces.clear();
        let reader: Reader = match error {
            Some(e) => Box::new(Failed(Some(e))),
            None => reader,
        };
        if text.is_empty() {
            reader
        } else {
            Box::new(Cursor::new(text).chain(reader))
        }
    }
}

/// Gives a piece's buffer back to its scanner once nothing else holds it.
fn recycle(p: Arc<Piece>, spare: &Sender<Vec<u8>>) {
    if let Ok(p) = Arc::try_unwrap(p) {
        let _ = spare.send(p.text);
    }
}

struct Splitter {
    lanes: Vec<Lane>,
    records: Vec<u64>,
    cfg: Config,
    jobs: Sender<Job>,
    results: Sender<Msg>,
    credit: SyncSender<()>,
    spare: Vec<Sender<Vec<u8>>>,
    stop: Arc<AtomicBool>,
    seq: u64,
}

/// Why the splitter stopped cutting.
enum Stop {
    /// Hand the inputs back.
    Sequential,
    /// The consumer, or a scanner, is gone.
    Closed,
}

impl Splitter {
    fn run(&mut self) {
        loop {
            match self.cut() {
                Ok(()) => {}
                Err(Stop::Closed) => return,
                Err(Stop::Sequential) => break,
            }
        }
        self.stop.store(true, Ordering::Relaxed);
        let readers = self.lanes.iter_mut().map(Lane::hand_back).collect();
        let tail = Tail {
            readers,
            records: self.records.clone(),
        };
        let _ = self.results.send(Msg::Tail(self.seq, tail));
    }

    /// Cuts one batch and sends it to the workers.
    fn cut(&mut self) -> Result<(), Stop> {
        let k = self.cfg.k;
        let mut used = 0usize;
        let mut spots = 0usize;
        let mut saved = vec![(0, 0); self.lanes.len()];
        while spots < self.cfg.batch_records {
            for (s, l) in saved.iter_mut().zip(&self.lanes) {
                *s = (l.cur, l.rec);
            }
            let mut bytes = 0;
            for s in 0..k {
                let r = self.cfg.reader_of(s);
                match self.lanes[r].next_record()? {
                    Some(len) => bytes += len,
                    None => return Err(Stop::Sequential),
                }
            }
            if used + bytes > self.cfg.batch_bytes && spots > 0 {
                // The spot goes into the next batch.
                for (l, &(cur, rec)) in self.lanes.iter_mut().zip(&saved) {
                    (l.cur, l.rec) = (cur, rec);
                }
                break;
            }
            used += bytes;
            spots += 1;
        }
        let text = self
            .lanes
            .iter_mut()
            .zip(&self.spare)
            .map(|(l, spare)| l.cut(spare))
            .collect();
        let records = self.records.clone();
        let stride = if self.cfg.interleaved { k as u64 } else { 1 };
        for r in &mut self.records {
            *r += spots as u64 * stride;
        }
        let job = Job {
            seq: self.seq,
            text,
            spots,
            records,
        };
        self.credit.send(()).map_err(|_| Stop::Closed)?;
        self.jobs.send(job).map_err(|_| Stop::Closed)?;
        self.seq += 1;
        Ok(())
    }
}

fn worker(
    jobs: Arc<Mutex<Receiver<Job>>>,
    results: Sender<Msg>,
    spare: Vec<Sender<Vec<u8>>>,
    cfg: Config,
) {
    let mut qual33 = vec![Vec::new(); cfg.k];
    loop {
        let job = match jobs.lock() {
            Ok(rx) => match rx.recv() {
                Ok(job) => job,
                Err(_) => return,
            },
            Err(_) => return,
        };
        let msg = match catch_unwind(AssertUnwindSafe(|| parse(&job, &cfg, &mut qual33))) {
            Ok(r) => Msg::Batch(job.seq, r),
            Err(p) => Msg::Panic(p),
        };
        for (ranges, spare) in job.text.into_iter().zip(&spare) {
            for (p, _) in ranges {
                recycle(p, spare);
            }
        }
        if results.send(msg).is_err() {
            return;
        }
    }
}

/// Validates and copies one batch, as the sequential fast path does.
fn parse(job: &Job, cfg: &Config, qual33: &mut [Vec<u8>]) -> Result<Batch<Segments>, Error> {
    let k = cfg.k;
    let offset = if cfg.phred64 { 64 } else { 33 };
    let mut b = BatchBuilder::with_segments(&cfg.pool, k, job.spots);
    // Per input: the range being read and the position in it.
    let mut at: Vec<(usize, usize)> = job
        .text
        .iter()
        .map(|r| (0, r.first().map_or(0, |(_, x)| x.start)))
        .collect();
    let mut count = vec![0u64; job.text.len()];
    let mut fields: Vec<Lines<'_>> = vec![[&[][..]; 4]; k];
    for _ in 0..job.spots {
        for (s, f) in fields.iter_mut().enumerate() {
            let r = cfg.reader_of(s);
            let ranges = &job.text[r];
            let (mut i, mut pos) = at[r];
            if pos == ranges[i].1.end {
                i += 1;
                pos = ranges[i].1.start;
            }
            let (piece, range) = &ranges[i];
            let (lines, end) = split_record(&piece.text[..range.end], pos)
                .expect("the splitter cut whole records");
            *f = lines;
            at[r] = (i, end);
        }
        for (s, f) in fields.iter().enumerate() {
            let r = cfg.reader_of(s);
            count[r] += 1;
            if let Err(kind) = check(f, offset) {
                return Err(cfg.error(kind, job.records[r] + count[r], s));
            }
        }
        push_spot(
            &mut b,
            &fields,
            cfg.keep_names,
            cfg.phred64.then_some(&mut qual33[..]),
        )
        .map_err(batch_error)?;
    }
    Ok(b.finish())
}

/// The parallel path against the sequential one on generated inputs: every
/// layout, option, reader shape and kind of damage. Pieces are 64 bytes in
/// tests, so records span pieces and buffers grow.
#[cfg(test)]
mod tests {
    use std::io::{BufReader, Cursor};

    use readkit_batch::Segments;
    use readkit_io::BatchSource;

    use crate::{ErrorKind, FastqInput, FastqOptions, FastqSource};

    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            z ^ (z >> 31)
        }

        fn below(&mut self, n: u64) -> u64 {
            self.next() % n.max(1)
        }
    }

    /// Reads at most a few bytes at a time, and fails after `fail_at` bytes.
    struct Chunky {
        data: Cursor<Vec<u8>>,
        rng: Rng,
        max: usize,
        fail_at: Option<u64>,
    }

    impl std::io::Read for Chunky {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            let pos = self.data.position();
            if self.fail_at.is_some_and(|f| pos >= f) {
                self.fail_at = None;
                return Err(std::io::Error::other("injected failure"));
            }
            let mut n = 1 + self.rng.below(self.max as u64) as usize;
            if let Some(f) = self.fail_at {
                n = n.min((f - pos) as usize).max(1);
            }
            let n = n.min(buf.len());
            std::io::Read::read(&mut self.data, &mut buf[..n])
        }
    }

    #[derive(Clone, Copy, Debug)]
    enum Shape {
        /// A slice: `fill_buf` returns all of it.
        Slice,
        /// A `BufReader` of this capacity over a slice.
        Buffered(usize),
        /// A `BufReader` over short random reads, failing after some bytes.
        Chunky(usize, Option<u64>),
    }

    fn reader(data: &[u8], shape: Shape, seed: u64) -> FastqInput {
        let data = data.to_vec();
        match shape {
            Shape::Slice => FastqInput::named(Cursor::new(data), "in"),
            Shape::Buffered(c) => {
                FastqInput::named(BufReader::with_capacity(c, Cursor::new(data)), "in")
            }
            Shape::Chunky(max, fail_at) => FastqInput::named(
                BufReader::with_capacity(
                    16,
                    Chunky {
                        data: Cursor::new(data),
                        rng: Rng(seed),
                        max,
                        fail_at,
                    },
                ),
                "in",
            ),
        }
    }

    /// FASTQ records with the given damage.
    fn text(rng: &mut Rng, records: usize, phred64: bool, damage: Damage) -> Vec<u8> {
        let mut out = Vec::new();
        let base = if phred64 { 64 } else { 33 };
        let bad_at = rng.below(records as u64) as usize;
        let eol: &[u8] = if matches!(damage, Damage::Crlf) {
            b"\r\n"
        } else {
            b"\n"
        };
        for i in 0..records {
            let len = match rng.below(20) {
                0 => 0,
                1 => 100 + rng.below(400) as usize,
                _ => 1 + rng.below(30) as usize,
            };
            let seq: Vec<u8> = (0..len).map(|_| b"ACGTN"[rng.below(5) as usize]).collect();
            let mut qual: Vec<u8> = (0..len)
                .map(|_| base + rng.below(127 - base as u64) as u8)
                .collect();
            let mut header = format!("@r{i} c{}", rng.below(1000)).into_bytes();
            let mut plus = b"+".to_vec();
            if i == bad_at {
                match damage {
                    Damage::Quality if len > 0 => qual[0] = base - 1,
                    Damage::Length => qual.push(base),
                    Damage::Separator => plus = b"-".to_vec(),
                    Damage::Header => header[0] = b'X',
                    Damage::Blank => out.extend_from_slice(b"\n \t\n"),
                    _ => {}
                }
            }
            if rng.below(4) == 0 {
                plus.extend_from_slice(&header[1..]);
            }
            for line in [&header[..], &seq, &plus, &qual] {
                out.extend_from_slice(line);
                out.extend_from_slice(eol);
            }
        }
        match damage {
            Damage::Trailing => out.extend_from_slice(b"\n\n  \n"),
            Damage::NoFinalNewline => {
                out.pop();
            }
            Damage::Truncated if !out.is_empty() => {
                let cut = rng.below(out.len() as u64) as usize;
                out.truncate(cut);
            }
            _ => {}
        }
        out
    }

    #[derive(Clone, Copy, Debug, PartialEq)]
    enum Damage {
        None,
        Crlf,
        Quality,
        Length,
        Separator,
        Header,
        Blank,
        Trailing,
        NoFinalNewline,
        Truncated,
        /// One input of a lockstep set has one record fewer.
        MateCount,
    }

    const DAMAGE: [Damage; 11] = [
        Damage::None,
        Damage::Crlf,
        Damage::Quality,
        Damage::Length,
        Damage::Separator,
        Damage::Header,
        Damage::Blank,
        Damage::Trailing,
        Damage::NoFinalNewline,
        Damage::Truncated,
        Damage::MateCount,
    ];

    type Record = (Option<Vec<u8>>, Vec<u8>, Vec<u8>);
    /// Kind, record, segment, input and message.
    type Failure = (
        ErrorKind,
        Option<u64>,
        Option<usize>,
        Option<String>,
        String,
    );

    #[derive(Debug, PartialEq)]
    struct Outcome {
        /// Per spot, per segment.
        spots: Vec<Vec<Record>>,
        /// Spots per batch.
        batches: Vec<usize>,
        error: Option<Failure>,
    }

    fn collect(mut src: FastqSource<Segments>, limit: usize) -> Outcome {
        let mut o = Outcome {
            spots: Vec::new(),
            batches: Vec::new(),
            error: None,
        };
        loop {
            match src.next_batch() {
                Ok(Some(b)) => {
                    assert!(!b.is_empty() && b.len() <= limit);
                    o.batches.push(b.len());
                    for i in 0..b.len() {
                        o.spots.push(
                            (0..b.segments())
                                .map(|s| {
                                    let c = b.segment(s);
                                    (
                                        c.names().map(|n| n.get(i).to_vec()),
                                        c.seqs().get(i).to_vec(),
                                        c.quals().get(i).to_vec(),
                                    )
                                })
                                .collect(),
                        );
                    }
                }
                Ok(None) => break,
                Err(e) => {
                    o.error = Some((
                        e.kind(),
                        e.record(),
                        e.segment(),
                        e.input().map(str::to_owned),
                        e.to_string(),
                    ));
                    assert!(src.next_batch().unwrap().is_none());
                    break;
                }
            }
        }
        assert!(src.next_batch().unwrap().is_none());
        o
    }

    struct Case {
        inputs: Vec<Vec<u8>>,
        k: usize,
        interleaved: bool,
        opts: FastqOptions,
        shape: Shape,
    }

    fn source(c: &Case, threads: usize) -> FastqSource<Segments> {
        let opts = c.opts.clone().threads(threads);
        let inputs: Vec<_> = c
            .inputs
            .iter()
            .enumerate()
            .map(|(i, d)| reader(d, c.shape, i as u64 + 7))
            .collect();
        if c.interleaved {
            FastqSource::interleaved(inputs.into_iter().next().unwrap(), c.k, opts).unwrap()
        } else {
            FastqSource::lockstep(inputs, opts).unwrap()
        }
    }

    fn case(rng: &mut Rng, damage: Damage) -> Case {
        let interleaved = rng.below(4) == 0;
        let k = if interleaved {
            2
        } else {
            1 + rng.below(3) as usize
        };
        let phred64 = rng.below(5) == 0;
        let records = rng.below(60) as usize;
        let n = if interleaved { 1 } else { k };
        let spot_records = if interleaved { records * k } else { records };
        let mut inputs: Vec<Vec<u8>> = (0..n)
            .map(|i| {
                // Damage in one input only (MateCount: the last one).
                let d = match damage {
                    Damage::MateCount => Damage::None,
                    d if i == n - 1 => d,
                    Damage::Crlf => Damage::Crlf,
                    _ => Damage::None,
                };
                let mut r = spot_records;
                if damage == Damage::MateCount && i == n - 1 && r > 0 {
                    r -= 1;
                }
                text(rng, r, phred64, d)
            })
            .collect();
        if damage == Damage::MateCount && interleaved {
            // An odd record count.
            let extra = text(rng, 1, phred64, Damage::None);
            inputs[0].extend_from_slice(&extra);
        }
        let budget = [1, 50, 300, 2000, 1 << 20][rng.below(5) as usize];
        let limit = [1, 3, usize::MAX][rng.below(3) as usize];
        let shape = match rng.below(4) {
            0 => Shape::Slice,
            1 => Shape::Buffered(1 + rng.below(200) as usize),
            2 => Shape::Chunky(1 + rng.below(100) as usize, None),
            _ => {
                let len = inputs[0].len() as u64;
                Shape::Chunky(1 + rng.below(100) as usize, Some(rng.below(len + 1)))
            }
        };
        let opts = FastqOptions::default()
            .phred64(phred64)
            .keep_names(rng.below(3) != 0)
            .batch_bytes(budget)
            .batch_records(limit);
        Case {
            inputs,
            k,
            interleaved,
            opts,
            shape,
        }
    }

    fn env(name: &str, default: u64) -> u64 {
        std::env::var(name)
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(default)
    }

    /// `READKIT_FASTQ_PAR_ROUNDS` and `READKIT_FASTQ_PAR_SEED` run more or
    /// other cases.
    #[test]
    fn parallel_matches_sequential_on_generated_inputs() {
        let rounds = env("READKIT_FASTQ_PAR_ROUNDS", 400) as usize;
        let mut rng = Rng(env("READKIT_FASTQ_PAR_SEED", 42));
        let mut checked = 0;
        for round in 0..rounds {
            for damage in DAMAGE {
                let c = case(&mut rng, damage);
                let limit = c.opts.batch_records;
                let budget = c.opts.batch_bytes;
                let want = collect(source(&c, 0), limit);
                for threads in [1, 3] {
                    let got = collect(source(&c, threads), limit);
                    let ctx = format!(
                        "round {round} {damage:?} k={} interleaved={} shape={:?} budget={budget} limit={limit} threads={threads}",
                        c.k, c.interleaved, c.shape
                    );
                    assert_eq!(got.error, want.error, "{ctx}");
                    assert!(
                        got.spots == want.spots,
                        "{ctx}: records differ ({} spots vs {})",
                        got.spots.len(),
                        want.spots.len()
                    );
                    assert_eq!(got.batches, want.batches, "{ctx}");
                    checked += 1;
                }
            }
        }
        assert_eq!(checked, rounds * DAMAGE.len() * 2);
    }

    /// Every batch the splitter cuts stays within the byte budget unless it
    /// holds one spot.
    #[test]
    fn parallel_batches_respect_the_budget() {
        let mut rng = Rng(7);
        let data = text(&mut rng, 3000, false, Damage::None);
        for budget in [1, 100, 1000, 10_000] {
            let opts = FastqOptions::default().batch_bytes(budget).threads(2);
            let mut src = FastqSource::<Segments>::lockstep(
                vec![FastqInput::new(Cursor::new(data.clone()))],
                opts,
            )
            .unwrap();
            let mut spots = 0;
            while let Some(b) = src.next_batch().unwrap() {
                let c = b.segment(0);
                let text: usize = (0..b.len())
                    .map(|i| {
                        // `@name\nseq\n+…\nqual\n`, the `+` line as written.
                        c.names().unwrap().get(i).len() + 2 * c.seqs().get(i).len() + 6
                    })
                    .sum();
                assert!(
                    b.len() == 1 || text <= budget + b.len() * 1100,
                    "budget {budget}"
                );
                spots += b.len();
            }
            assert_eq!(spots, 3000);
        }
    }

    /// Dropping a source mid-stream stops its threads (the test would hang
    /// otherwise), on every reader shape.
    #[test]
    fn dropping_mid_stream_stops_the_threads() {
        let mut rng = Rng(9);
        let data = text(&mut rng, 2000, false, Damage::None);
        for shape in [Shape::Slice, Shape::Buffered(100), Shape::Chunky(50, None)] {
            for take in [0, 1, 5] {
                let opts = FastqOptions::default().batch_bytes(500).threads(3);
                let mut src = FastqSource::<Segments>::lockstep(
                    vec![reader(&data, shape, 1), reader(&data, shape, 2)],
                    opts,
                )
                .unwrap();
                for _ in 0..take {
                    src.next_batch().unwrap().unwrap();
                }
                drop(src);
            }
        }
    }
}
