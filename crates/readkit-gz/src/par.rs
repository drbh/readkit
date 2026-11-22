//! Parallel decoding of positioned input, behind `Read` / `BufRead`.

use crate::error::{Error, ErrorKind};
use crate::header;
use crate::pgzip::{self, Output, Plan};
use crate::read_at::{Cursor, ReadAt};
use crate::sequential::GzReader;
use crossbeam_channel::{Receiver, Sender};
use std::io::{self, BufRead, Read};
use std::path::Path;
use std::sync::Arc;

/// When to cut one large gzip member into regions that are decoded in
/// parallel (the pugz / rapidgzip method).
///
/// Region decoding costs about 2.5 times zlib's CPU per byte, so it pays off
/// only when decoding limits the whole run. fasterp, for example, splits
/// only for a single input whose output is not gzip or whose qualities are
/// binned ([`binned_qualities`](crate::binned_qualities)), and uses
/// [`Split::Never`] otherwise.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Split {
    /// Split when there are at least two threads and the member spans at
    /// least three regions. Members that record their sizes (BGZF, fasterp
    /// and [`GzWriter`](crate::GzWriter) output) are decoded member-parallel
    /// instead.
    #[default]
    Auto,
    /// Never split a member. Sized members are still decoded in parallel;
    /// anything else is decoded sequentially on one background thread.
    Never,
    /// Always decode by regions, whatever the size, even for sized members.
    /// For tests and fuzzing, with a small
    /// [`region_size`](DecodeOptions::region_size).
    Always,
}

/// Options for [`ParGzReader`].
///
/// Construct with `DecodeOptions { threads: Some(4), ..Default::default() }`.
///
/// Memory: region decoding has at most `in_flight` regions between being
/// queued and being delivered. Each needs up to three times its decoded
/// size, `region_size` times the compression ratio (about 9 MiB for FASTQ
/// at the default 2 MiB): two bytes per decoded byte until references to
/// the unknown window before it are resolved, which for FASTQ is usually
/// the whole region, then the bytes themselves. Three more decoded regions
/// are queued to or held by the reader, and each thread has a compressed
/// buffer of `region_size` plus 1 MiB. Region buffers are reserved for 9/8
/// of the largest region seen, and for five times `region_size` before the
/// first; whatever a region does not use is address space, not resident.
/// Buffers are reused rather than freed, so resident memory stays near the
/// live heap: on macOS a 1.6 GB FASTQ decode with 2 MiB regions peaked at
/// 75 MB (1 thread), 210 MB (4) and 300 MB (8). Sized members hold at most
/// `2 × in_flight` batches of `buffer_size` decoded bytes; sequential
/// decoding holds a few buffers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DecodeOptions {
    /// Decoding threads. `None` uses fasterp's defaults: half the available
    /// cores clamped to 2..=8 for region decoding, a quarter clamped to
    /// 1..=4 for sized members (libdeflate is fast). One more thread always
    /// coordinates and delivers output.
    pub threads: Option<usize>,
    /// When to split a member into regions.
    pub split: Split,
    /// Compressed bytes per region (default 2 MiB).
    pub region_size: usize,
    /// Regions, or batches of sized members, decoded ahead of the reader.
    /// `None` means `threads + 2` regions or `2 × threads` batches.
    pub in_flight: Option<usize>,
    /// Decoded bytes per output buffer, and per batch of sized members
    /// (default 1 MiB).
    pub buffer_size: usize,
}

impl Default for DecodeOptions {
    fn default() -> Self {
        Self {
            threads: None,
            split: Split::Auto,
            region_size: pgzip::DEFAULT_REGION,
            in_flight: None,
            buffer_size: 1 << 20,
        }
    }
}

impl DecodeOptions {
    fn check(&self) -> Result<(), Error> {
        let bad = |what: &str| Err(Error::new(ErrorKind::InvalidInput, what.to_string()));
        if self.threads == Some(0) {
            return bad("threads must be at least 1");
        }
        if self.in_flight == Some(0) {
            return bad("in_flight must be at least 1");
        }
        if self.region_size == 0 {
            return bad("region_size must be at least 1");
        }
        if self.buffer_size == 0 {
            return bad("buffer_size must be at least 1");
        }
        Ok(())
    }

    fn cores() -> usize {
        std::thread::available_parallelism().map_or(1, |n| n.get())
    }

    fn region_threads(&self) -> usize {
        self.threads
            .unwrap_or_else(|| Self::cores().div_ceil(2).clamp(2, 8))
    }

    fn member_threads(&self) -> usize {
        self.threads
            .unwrap_or_else(|| Self::cores().div_ceil(4).clamp(1, 4))
    }
}

/// How a [`ParGzReader`] decodes its input, chosen when it is opened.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Strategy {
    /// Sequentially, on one background thread.
    Sequential,
    /// The first member by regions in parallel; later members too while they
    /// stay large, the rest sequentially.
    Regions,
    /// Members that record their compressed size (BGZF `BC`, fasterp `FP`)
    /// in parallel with libdeflate; from the first member without a size,
    /// sequentially.
    Members,
}

/// gzip decoder that decodes on background threads and is read like a file.
///
/// - [`ParGzReader::new`] and [`open`](crate::open) take positioned input
///   ([`ReadAt`] plus its length) and decode in parallel where the input
///   allows ([`Strategy`]).
/// - [`ParGzReader::from_reader`] takes any `Read`, decoded sequentially
///   on one background thread so decoding overlaps with the consumer.
///
/// Output is identical to [`GzReader`]'s for every strategy and thread
/// count, with the same checks: every member's header, CRC32 and ISIZE.
/// Reaching end of file means everything was verified; errors arrive in
/// order, after all the output that precedes them. It implements
/// [`BufRead`], so wrapping it in a `BufReader` is unnecessary.
///
/// Dropping the reader stops its threads and, for positioned input, waits
/// for them (at most the reads and region decodes in progress).
pub struct ParGzReader {
    receiver: Option<Receiver<Result<Vec<u8>, Error>>>,
    recycle: Sender<Vec<u8>>,
    worker: Option<std::thread::JoinHandle<()>>,
    /// Join the worker on drop; not done for streams, whose read may block.
    join: bool,
    buffer: Vec<u8>,
    position: usize,
    finished: bool,
    error: Option<Error>,
    strategy: Strategy,
}

impl std::fmt::Debug for ParGzReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ParGzReader")
            .field("strategy", &self.strategy)
            .field("finished", &self.finished)
            .finish_non_exhaustive()
    }
}

/// Open the gzip file at `path` for parallel decoding.
///
/// A regular file is decoded with [`ParGzReader::new`]. Anything else that
/// opens (a FIFO, `/dev/stdin`) is decoded sequentially, as
/// [`ParGzReader::from_reader`] does.
pub fn open(path: impl AsRef<Path>, options: &DecodeOptions) -> Result<ParGzReader, Error> {
    let path = path.as_ref();
    let io = |e: io::Error| {
        Error::new(ErrorKind::Io, format!("opening {}: {e}", path.display())).with_source(e)
    };
    let file = std::fs::File::open(path).map_err(io)?;
    let meta = file.metadata().map_err(io)?;
    if meta.is_file() {
        ParGzReader::new(file, meta.len(), options)
    } else {
        ParGzReader::from_reader(file, options)
    }
}

impl ParGzReader {
    /// Decode the `len` bytes of `src` in parallel where possible.
    ///
    /// Reads the start of the input to choose a [`Strategy`]; fails only if
    /// that read fails or the options are invalid.
    pub fn new<S: ReadAt + 'static>(
        src: S,
        len: u64,
        options: &DecodeOptions,
    ) -> Result<Self, Error> {
        options.check()?;
        let len = usize::try_from(len)
            .map_err(|_| Error::new(ErrorKind::InvalidInput, "input too large"))?;
        let mut head = vec![0; len.min((64 << 10) + 12)];
        src.read_exact_at(&mut head, 0).map_err(Error::input)?;
        let sized = header::member_size(&head).is_some();
        let strategy = match options.split {
            Split::Always => Strategy::Regions,
            _ if sized => Strategy::Members,
            Split::Auto
                if options.region_threads() >= 2
                    && pgzip::worth_it(len as u64, options.region_size) =>
            {
                Strategy::Regions
            }
            _ => Strategy::Sequential,
        };
        let src: Arc<dyn ReadAt> = Arc::new(src);
        let buffer = options.buffer_size;
        Ok(match strategy {
            Strategy::Regions => {
                let threads = options.region_threads();
                let plan = Plan {
                    region: options.region_size,
                    threads,
                    ahead: options.in_flight.unwrap_or(threads + 2),
                    buffer,
                };
                // Finished buffers go back to the decoder for reuse; the
                // channel holds every buffer that can exist, so none is freed.
                let (recycle, returned) =
                    crossbeam_channel::bounded(pgzip::output_buffers(plan.ahead));
                let pool = (recycle.clone(), returned);
                Self::spawn(
                    strategy,
                    pgzip::OUTPUT_DEPTH,
                    recycle,
                    true,
                    move |sender| pgzip::decompress(&*src, len, plan, &sender, pool),
                )?
            }
            Strategy::Members => {
                let threads = options.member_threads();
                let ahead = options.in_flight.unwrap_or(2 * threads);
                let (recycle, _) = crossbeam_channel::bounded(1);
                Self::spawn(strategy, 2, recycle, true, move |sender| {
                    members(&*src, len, threads, ahead, buffer, &sender)
                })?
            }
            Strategy::Sequential => {
                let (recycle, reused) = crossbeam_channel::bounded(2);
                Self::spawn(strategy, 2, recycle, true, move |sender| {
                    let cursor = Cursor {
                        src: &*src,
                        pos: 0,
                        end: len as u64,
                    };
                    stream(cursor, buffer, &sender, &reused)
                })?
            }
        })
    }

    /// Decode any reader sequentially on one background thread, buffering
    /// up to two output buffers ahead of the consumer.
    ///
    /// Dropping the reader does not wait for that thread, which may be
    /// blocked reading `reader`; it exits at its next read.
    pub fn from_reader<R: Read + Send + 'static>(
        reader: R,
        options: &DecodeOptions,
    ) -> Result<Self, Error> {
        options.check()?;
        let buffer = options.buffer_size;
        let (recycle, reused) = crossbeam_channel::bounded(2);
        Self::spawn(Strategy::Sequential, 2, recycle, false, move |sender| {
            stream(reader, buffer, &sender, &reused)
        })
    }

    fn spawn(
        strategy: Strategy,
        depth: usize,
        recycle: Sender<Vec<u8>>,
        join: bool,
        work: impl FnOnce(Output) + Send + 'static,
    ) -> Result<Self, Error> {
        let (sender, receiver) = crossbeam_channel::bounded(depth);
        let worker = std::thread::Builder::new()
            .name("readkit-gz".into())
            .spawn(move || work(sender))
            .map_err(|e| {
                Error::new(ErrorKind::Io, "cannot start a decoding thread").with_source(e)
            })?;
        Ok(Self {
            receiver: Some(receiver),
            recycle,
            worker: Some(worker),
            join,
            buffer: Vec::new(),
            position: 0,
            finished: false,
            error: None,
            strategy,
        })
    }

    /// The strategy chosen for this input.
    pub fn strategy(&self) -> Strategy {
        self.strategy
    }

    fn next_buffer(&mut self) -> Result<(), Error> {
        let next = self
            .receiver
            .as_ref()
            .expect("receiver present until drop")
            .recv();
        let next = match next {
            Ok(Ok(next)) => next,
            Ok(Err(e)) => return Err(e),
            Err(_) => {
                return Err(Error::internal(
                    "decoding thread stopped before the end of the stream",
                ));
            }
        };
        let old = std::mem::replace(&mut self.buffer, next);
        let _ = self.recycle.try_send(old);
        self.position = 0;
        self.finished = self.buffer.is_empty();
        Ok(())
    }
}

impl BufRead for ParGzReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if let Some(e) = &self.error {
            return Err(e.clone().into());
        }
        while self.position == self.buffer.len() && !self.finished {
            if let Err(e) = self.next_buffer() {
                self.error = Some(e.clone());
                return Err(e.into());
            }
        }
        Ok(&self.buffer[self.position..])
    }

    fn consume(&mut self, amount: usize) {
        self.position = (self.position + amount).min(self.buffer.len());
    }
}

impl Read for ParGzReader {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        if output.is_empty() {
            return Ok(0);
        }
        let available = self.fill_buf()?;
        let count = output.len().min(available.len());
        output[..count].copy_from_slice(&available[..count]);
        self.consume(count);
        Ok(count)
    }
}

impl Drop for ParGzReader {
    fn drop(&mut self) {
        // Disconnect first, so a decoder blocked sending output stops.
        drop(self.receiver.take());
        if let Some(worker) = self.worker.take()
            && self.join
        {
            let _ = worker.join();
        }
    }
}

/// Decode `input` sequentially into buffers of `buffer` bytes.
fn stream<R: Read>(input: R, buffer: usize, sender: &Output, reused: &Receiver<Vec<u8>>) {
    let mut decoder = GzReader::from_reader(input);
    loop {
        let mut chunk = reused.try_recv().unwrap_or_default();
        chunk.resize(buffer, 0);
        let mut filled = 0;
        let mut failed = None;
        while filled < chunk.len() {
            match decoder.read_typed(&mut chunk[filled..]) {
                Ok(0) => break,
                Ok(n) => filled += n,
                Err(e) => {
                    failed = Some(e);
                    break;
                }
            }
        }
        chunk.truncate(filled);
        // Output before an error is delivered first.
        if filled > 0 && sender.send(Ok(chunk)).is_err() {
            return;
        }
        if let Some(e) = failed {
            let _ = sender.send(Err(e));
            return;
        }
        if filled == 0 {
            let _ = sender.send(Ok(Vec::new()));
            return;
        }
    }
}

/// libdeflate's decompressor, for whole gzip members.
struct Libdeflate(std::ptr::NonNull<libdeflate_sys::libdeflate_decompressor>);

// SAFETY: a libdeflate decompressor holds no thread-local state; it is used
// by one thread at a time through `&mut self`.
unsafe impl Send for Libdeflate {}

impl Libdeflate {
    fn new() -> Self {
        // SAFETY: plain allocation, checked for null below.
        let p = unsafe { libdeflate_sys::libdeflate_alloc_decompressor() };
        Self(std::ptr::NonNull::new(p).expect("libdeflate could not allocate a decompressor"))
    }

    /// Decode the gzip member that is exactly `member` into exactly `out`.
    /// libdeflate checks the trailer's CRC32 and ISIZE; this also requires
    /// that the member ends where its size field says and fills `out`.
    fn member(&mut self, member: &[u8], out: &mut [u8]) -> bool {
        let (mut used, mut made) = (0usize, 0usize);
        // SAFETY: the pointers and lengths come from live slices, and the
        // decompressor is exclusively borrowed.
        let result = unsafe {
            libdeflate_sys::libdeflate_gzip_decompress_ex(
                self.0.as_ptr(),
                member.as_ptr().cast(),
                member.len(),
                out.as_mut_ptr().cast(),
                out.len(),
                &mut used,
                &mut made,
            )
        };
        result == libdeflate_sys::libdeflate_result_LIBDEFLATE_SUCCESS
            && used == member.len()
            && made == out.len()
    }
}

impl Drop for Libdeflate {
    fn drop(&mut self) {
        // SAFETY: allocated in `new`, freed once.
        unsafe { libdeflate_sys::libdeflate_free_decompressor(self.0.as_ptr()) }
    }
}

/// A batch of whole sized members, back to back.
struct Batch {
    /// Offset of the first member in the input.
    start: usize,
    bytes: Vec<u8>,
    /// Total decompressed size, from the members' trailers.
    size: usize,
}

enum Job {
    Members(Batch),
    /// Stop here and decode the rest sequentially from this offset.
    Rest(usize),
}

/// The members of `src` from `at` on: each one's total size, while the
/// members record it and lie inside the input.
fn next_batch(
    src: &dyn ReadAt,
    len: usize,
    at: &mut usize,
    window: &mut Vec<u8>,
    batch_bytes: usize,
) -> Result<Option<Batch>, Error> {
    let start = *at;
    let mut bytes = Vec::new();
    let mut size = 0usize;
    while size < batch_bytes && *at < len {
        // The fixed header, then the extra field holding the size.
        if len - *at < 12 {
            break;
        }
        window.resize(12, 0);
        src.read_exact_at(window, *at as u64)
            .map_err(Error::input)?;
        let xlen = usize::from(u16::from_le_bytes([window[10], window[11]]));
        if window[3] & 4 == 0 || len - *at < 12 + xlen {
            break;
        }
        window.resize(12 + xlen, 0);
        src.read_exact_at(&mut window[12..], *at as u64 + 12)
            .map_err(Error::input)?;
        let Some(total) = header::member_size(window).filter(|&t| t >= 18 + 10 && t <= len - *at)
        else {
            break;
        };
        let from = bytes.len();
        bytes.resize(from + total, 0);
        src.read_exact_at(&mut bytes[from..], *at as u64)
            .map_err(Error::input)?;
        let isize = u32::from_le_bytes(bytes[from + total - 4..].try_into().unwrap());
        size += isize as usize;
        *at += total;
    }
    Ok((!bytes.is_empty()).then_some(Batch { start, bytes, size }))
}

/// Decode a batch; `None` if any member fails, which the sequential decoder
/// then explains.
fn decode_batch(inflater: &mut Libdeflate, batch: &Batch) -> Option<Vec<u8>> {
    let mut out = vec![0; batch.size];
    let (mut at, mut written) = (0, 0);
    while at < batch.bytes.len() {
        let rest = &batch.bytes[at..];
        let total = header::member_size(rest)?;
        let member = &rest[..total];
        // libdeflate skips the optional header CRC; check the header here.
        header::parse(member).ok()??;
        let isize = u32::from_le_bytes(member[total - 4..].try_into().unwrap()) as usize;
        let slot = out.get_mut(written..written + isize)?;
        if !inflater.member(member, slot) {
            return None;
        }
        written += isize;
        at += total;
    }
    Some(out)
}

/// Decompress members that record their sizes on `threads` threads with
/// libdeflate, sending output in order. From the first member without a
/// size, or the first batch that fails, the rest is decoded sequentially,
/// which also produces the exact error.
fn members(
    src: &dyn ReadAt,
    len: usize,
    threads: usize,
    ahead: usize,
    batch_bytes: usize,
    sender: &Output,
) {
    use std::collections::BTreeMap;
    let rest = std::thread::scope(|scope| -> Option<usize> {
        let (jobs, queue) = crossbeam_channel::bounded::<(usize, Job)>(ahead);
        let (done, results) = crossbeam_channel::bounded::<(usize, Result<Vec<u8>, usize>)>(ahead);
        for _ in 0..threads {
            let (queue, done) = (queue.clone(), done.clone());
            scope.spawn(move || {
                let mut inflater = Libdeflate::new();
                for (seq, job) in queue {
                    let result = match job {
                        Job::Rest(at) => Err(at),
                        Job::Members(batch) => {
                            decode_batch(&mut inflater, &batch).ok_or(batch.start)
                        }
                    };
                    if done.send((seq, result)).is_err() {
                        break;
                    }
                }
            });
        }
        drop((queue, done));
        // Cut the input into batches of whole members.
        scope.spawn(move || {
            let mut window = Vec::new();
            let mut at = 0;
            let mut seq = 0;
            loop {
                let job = match next_batch(src, len, &mut at, &mut window, batch_bytes) {
                    Ok(Some(batch)) => Job::Members(batch),
                    Ok(None) | Err(_) => Job::Rest(at),
                };
                let last = matches!(job, Job::Rest(_));
                if jobs.send((seq, job)).is_err() || last {
                    return;
                }
                seq += 1;
            }
        });
        // Deliver in order.
        let mut pending = BTreeMap::new();
        let mut next = 0;
        for (seq, result) in results {
            pending.insert(seq, result);
            while let Some(result) = pending.remove(&next) {
                match result {
                    Ok(out) => {
                        if !out.is_empty() && sender.send(Ok(out)).is_err() {
                            return None;
                        }
                    }
                    // Dropping `results` on return stops the workers and
                    // then the splitter.
                    Err(at) => return Some(at),
                }
                next += 1;
            }
        }
        None
    });
    let Some(at) = rest else {
        return;
    };
    let ok = if at == len {
        Ok(true)
    } else {
        pgzip::sequential_rest(
            src,
            len,
            at,
            batch_bytes,
            sender,
            &crossbeam_channel::never(),
        )
    };
    match ok {
        Ok(true) => {
            let _ = sender.send(Ok(Vec::new()));
        }
        Ok(false) => {}
        Err(e) => {
            let _ = sender.send(Err(e));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{EncodeOptions, GzWriter};
    use std::io::Write;

    fn read_all(gz: &[u8], options: &DecodeOptions) -> Result<Vec<u8>, Error> {
        let mut out = Vec::new();
        let mut reader = ParGzReader::new(gz.to_vec(), gz.len() as u64, options)?;
        reader
            .read_to_end(&mut out)
            .map_err(|e| Error::from_io(&e).cloned().unwrap())?;
        Ok(out)
    }

    /// fasterp's own gzip output, written with the parallel compressor.
    fn fasterp_gzip(data: &[u8]) -> Vec<u8> {
        let options = EncodeOptions {
            threads: 3,
            ..Default::default()
        };
        let mut writer = GzWriter::new(Vec::new(), &options).unwrap();
        writer.write_all(data).unwrap();
        writer.finish().unwrap()
    }

    fn plain_gzip(data: &[u8]) -> Vec<u8> {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(6));
        encoder.write_all(data).unwrap();
        encoder.finish().unwrap()
    }

    /// A BGZF member: `BC` subfield holding the member size minus one.
    fn bgzf_member(data: &[u8]) -> Vec<u8> {
        let plain = plain_gzip(data);
        let total = plain.len() + 8;
        let mut out = plain[..3].to_vec();
        out.push(plain[3] | 4);
        out.extend_from_slice(&plain[4..10]);
        out.extend_from_slice(&6u16.to_le_bytes());
        out.extend_from_slice(b"BC");
        out.extend_from_slice(&2u16.to_le_bytes());
        out.extend_from_slice(&((total - 1) as u16).to_le_bytes());
        out.extend_from_slice(&plain[10..]);
        out
    }

    #[test]
    fn sized_members_decompress_in_parallel_and_in_order() {
        let data: Vec<u8> = (0..5_000_000u32)
            .map(|i| b"ACGTN\n@+"[(i.wrapping_mul(2654435761) >> 29) as usize])
            .collect();
        let ours = fasterp_gzip(&data);
        assert!(
            header::member_size(&ours).is_some(),
            "output records member sizes"
        );
        let options = DecodeOptions::default();
        for (name, bytes, expected) in [
            ("fasterp", ours.clone(), data.clone()),
            // Sized members, then a member without a size: the rest is read
            // sequentially.
            (
                "mixed",
                [ours.clone(), plain_gzip(b"tail\n"), fasterp_gzip(b"more\n")].concat(),
                [&data[..], b"tail\n", b"more\n"].concat(),
            ),
            (
                "bgzf",
                data.chunks(60_000).flat_map(bgzf_member).collect(),
                data.clone(),
            ),
            ("empty", fasterp_gzip(b""), Vec::new()),
        ] {
            let reader = ParGzReader::new(bytes.clone(), bytes.len() as u64, &options).unwrap();
            assert_eq!(reader.strategy(), Strategy::Members, "{name}");
            assert!(read_all(&bytes, &options).unwrap() == expected, "{name}");
        }
        // A corrupt member and a truncated file are errors, not short output.
        let mut corrupt = ours.clone();
        let middle = corrupt.len() / 2;
        corrupt[middle] ^= 0xff;
        assert!(read_all(&corrupt, &options).is_err(), "corrupt member");
        let e = read_all(&ours[..ours.len() - 100], &options).unwrap_err();
        assert_eq!(e.kind(), ErrorKind::Truncated, "truncated file");
    }

    #[test]
    fn sequential_members_check_every_trailer() {
        // Member ends at, just before and just after the decoder's 1 MiB
        // buffer size, plus an empty member, all checked against their CRC.
        const BUFFER: usize = 1 << 20;
        let sizes = [BUFFER, 2 * BUFFER, BUFFER - 1, BUFFER + 1, 0, 1, 3 * BUFFER];
        let parts: Vec<Vec<u8>> = sizes
            .iter()
            .enumerate()
            .map(|(m, &n)| (0..n).map(|i| b"ACGT\n"[(i * 7 + m) % 5]).collect())
            .collect();
        let members: Vec<Vec<u8>> = parts.iter().map(|p| plain_gzip(p)).collect();
        let options = DecodeOptions {
            split: Split::Never,
            ..Default::default()
        };
        assert!(read_all(&members.concat(), &options).unwrap() == parts.concat());
        for bad in 0..members.len() {
            let mut corrupt = members.clone();
            let crc = corrupt[bad].len() - 8;
            corrupt[bad][crc] ^= 1;
            let e = read_all(&corrupt.concat(), &options).unwrap_err();
            assert_eq!(e.kind(), ErrorKind::Checksum, "bad CRC in member {bad}");
            let start: usize = members[..bad].iter().map(Vec::len).sum();
            assert_eq!(e.offset(), Some(start as u64), "offset of member {bad}");
        }
    }

    #[test]
    fn gzip_headers_are_checked_like_zlib() {
        // Found by fuzzing: reserved flag bits and the optional header CRC
        // were ignored, so damaged headers were accepted.
        let data = b"@r\nACGT\n+\nFFFF\n".repeat(1000);
        let plain = plain_gzip(&data);
        let with_flags = |flags: u8, crc: Option<u16>| {
            let mut out = plain[..10].to_vec();
            out[3] = flags;
            let crc = crc.unwrap_or(crc32fast::hash(&out) as u16);
            out.extend_from_slice(&crc.to_le_bytes());
            out.extend_from_slice(&plain[10..]);
            out
        };
        let good = with_flags(2, None);
        let bad_crc = with_flags(2, Some(crc32fast::hash(&good[..10]) as u16 ^ 1));
        let mut reserved = plain.clone();
        reserved[3] = 0x20;
        let sequential = DecodeOptions {
            split: Split::Never,
            ..Default::default()
        };
        let parallel = DecodeOptions {
            split: Split::Always,
            region_size: 1024,
            threads: Some(2),
            ..Default::default()
        };
        for (name, bytes, ok) in [
            ("header crc", &good, true),
            ("bad header crc", &bad_crc, false),
            ("reserved flag", &reserved, false),
        ] {
            for (how, options) in [("sequential", &sequential), ("parallel", &parallel)] {
                let result = read_all(bytes, options);
                assert_eq!(result.is_ok(), ok, "{name}, {how}");
                if ok {
                    assert!(result.unwrap() == data);
                }
            }
        }
    }
}
