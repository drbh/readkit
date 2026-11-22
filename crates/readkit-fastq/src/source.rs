//! The strict FASTQ reader (fasterp's `compat::Reader` and `compat/block.rs`,
//! generalized to `k` segments).

use std::fmt;
use std::io::BufRead;
use std::marker::PhantomData;
use std::sync::Arc;

use readkit_batch::{Batch, BatchBuilder, Paired, Pool, Read, Segments, Single};
use readkit_io::BatchSource;

use crate::error::{Error, ErrorKind};
use crate::{FastqLayout, FastqOptions};

/// One FASTQ text input: a buffered reader and an optional name (a path or
/// URL) used in error messages.
///
/// The reader's buffer size matters for speed: records are parsed in place
/// from the buffer, and only records crossing a buffer boundary are copied
/// line by line. 1 MiB (`BufReader::with_capacity(1 << 20, ..)`, as fasterp
/// uses) is a good choice for files.
pub struct FastqInput {
    reader: Box<dyn BufRead + Send>,
    name: Option<Arc<str>>,
}

impl FastqInput {
    /// An unnamed input.
    pub fn new(reader: impl BufRead + Send + 'static) -> Self {
        Self {
            reader: Box::new(reader),
            name: None,
        }
    }

    /// An input named `name` in error messages.
    pub fn named(reader: impl BufRead + Send + 'static, name: impl Into<Arc<str>>) -> Self {
        Self {
            reader: Box::new(reader),
            name: Some(name.into()),
        }
    }

    /// The input's name, if any.
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }
}

impl fmt::Debug for FastqInput {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FastqInput")
            .field("name", &self.name)
            .finish()
    }
}

/// Reads FASTQ text into [`Batch`]es: `k` inputs in lockstep, or one
/// interleaved input of `k`-record spots.
///
/// Each call to [`next_batch`](BatchSource::next_batch) fills one batch
/// until it holds [`FastqOptions::batch_bytes`] of text or
/// [`FastqOptions::batch_records`] spots. Records are validated as they are
/// read (see the crate docs for the rules), so every batch it returns is
/// valid. After an error, or at the end of input, it returns `Ok(None)`.
///
/// By default the source reads on the caller's thread, and a consuming stage
/// wraps it in [`ReadAhead`](readkit_io::ReadAhead) to overlap parsing with
/// its own work. With [`FastqOptions::threads`] set, it parses on threads of
/// its own; batches, records and errors are the same either way.
pub struct FastqSource<L> {
    readers: Vec<Box<dyn BufRead + Send>>,
    names: Vec<Option<Arc<str>>>,
    /// FASTQ records consumed per reader.
    records: Vec<u64>,
    k: usize,
    interleaved: bool,
    phred64: bool,
    keep_names: bool,
    batch_bytes: usize,
    batch_records: usize,
    pool: Pool,
    /// Per segment: slow-path record text and its four line ends.
    slow: Vec<(Vec<u8>, [usize; 4])>,
    /// A slow-path spot read while filling a batch it did not fit: it starts
    /// the next one.
    held: Option<Held>,
    /// Per segment: Phred+33 conversion of the current quality line.
    qual33: Vec<Vec<u8>>,
    spots: u64,
    hint: usize,
    done: bool,
    /// Worker threads for [`FastqOptions::threads`]; 0 once started or off.
    #[cfg_attr(target_arch = "wasm32", allow(dead_code))]
    threads: usize,
    #[cfg(not(target_arch = "wasm32"))]
    par: Option<crate::par::Par>,
    _layout: PhantomData<fn() -> L>,
}

impl FastqSource<Single> {
    /// Single-end reads from one input.
    pub fn single(input: FastqInput, opts: FastqOptions) -> Self {
        Self::make(vec![input], 1, false, opts)
    }
}

impl FastqSource<Paired> {
    /// Mate pairs from two inputs read in lockstep.
    pub fn paired(r1: FastqInput, r2: FastqInput, opts: FastqOptions) -> Self {
        Self::make(vec![r1, r2], 2, false, opts)
    }
}

impl<L: FastqLayout> FastqSource<L> {
    /// Spots of `inputs.len()` segments, one input per segment, read in
    /// lockstep (for example R1, R2, I1, I2).
    ///
    /// Fails with [`ErrorKind::Options`] if there are no inputs or the
    /// layout fixes another segment count.
    pub fn lockstep(inputs: Vec<FastqInput>, opts: FastqOptions) -> Result<Self, Error> {
        let k = inputs.len();
        check_k::<L>(k)?;
        Ok(Self::make(inputs, k, false, opts))
    }

    /// Spots of `k` consecutive records each from one interleaved input
    /// (`k = 2` for ordinary interleaved paired-end FASTQ).
    ///
    /// Fails with [`ErrorKind::Options`] if `k` is 0 or the layout fixes
    /// another segment count.
    pub fn interleaved(input: FastqInput, k: usize, opts: FastqOptions) -> Result<Self, Error> {
        check_k::<L>(k)?;
        Ok(Self::make(vec![input], k, true, opts))
    }
}

fn check_k<L: FastqLayout>(k: usize) -> Result<(), Error> {
    if k == 0 {
        return Err(Error::new(ErrorKind::Options, "no segments"));
    }
    if let Some(want) = L::SEGMENTS
        && want != k
    {
        return Err(Error::new(
            ErrorKind::Options,
            format!("{k} segments for a layout of {want}"),
        ));
    }
    Ok(())
}

impl<L> FastqSource<L> {
    fn make(inputs: Vec<FastqInput>, k: usize, interleaved: bool, opts: FastqOptions) -> Self {
        let (readers, names): (Vec<_>, Vec<_>) =
            inputs.into_iter().map(|i| (i.reader, i.name)).unzip();
        let n = readers.len();
        Self {
            readers,
            names,
            records: vec![0; n],
            k,
            interleaved,
            phred64: opts.phred64,
            keep_names: opts.keep_names,
            batch_bytes: opts.batch_bytes,
            batch_records: opts.batch_records.max(1),
            pool: opts.pool.unwrap_or_default(),
            slow: vec![(Vec::new(), [0; 4]); k],
            held: None,
            qual33: vec![Vec::new(); k],
            spots: 0,
            hint: 0,
            done: false,
            threads: opts.threads,
            #[cfg(not(target_arch = "wasm32"))]
            par: None,
            _layout: PhantomData,
        }
    }

    /// Segments per spot (`k`).
    pub fn segments(&self) -> usize {
        self.k
    }

    /// Whether the source reads one interleaved input.
    pub fn is_interleaved(&self) -> bool {
        self.interleaved
    }

    /// Spots returned so far.
    pub fn spots(&self) -> u64 {
        self.spots
    }

    /// The pool batches draw from.
    pub fn pool(&self) -> &Pool {
        &self.pool
    }
}

impl<L: FastqLayout> BatchSource<L> for FastqSource<L> {
    type Error = Error;

    fn next_batch(&mut self) -> Result<Option<Batch<L>>, Error> {
        if self.done {
            return Ok(None);
        }
        match self.next_segments() {
            Ok(Some(b)) => Ok(Some(L::from_segments(b))),
            Ok(None) => {
                self.done = true;
                Ok(None)
            }
            Err(e) => {
                self.done = true;
                Err(e)
            }
        }
    }
}

impl<L> fmt::Debug for FastqSource<L> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FastqSource")
            .field("segments", &self.k)
            .field("interleaved", &self.interleaved)
            .field("inputs", &self.names)
            .field("spots", &self.spots)
            .finish()
    }
}

/// The four lines of one record, without line endings.
pub(crate) type Lines<'a> = [&'a [u8]; 4];

/// Spots with at most this many segments are pushed without allocating.
const INLINE_SEGMENTS: usize = 8;

impl<L> FastqSource<L> {
    fn next_segments(&mut self) -> Result<Option<Batch<Segments>>, Error> {
        #[cfg(not(target_arch = "wasm32"))]
        {
            use crate::par::{Config, Next, Par};
            if self.threads > 0 {
                let workers = std::mem::take(&mut self.threads);
                let cfg = Config {
                    k: self.k,
                    interleaved: self.interleaved,
                    phred64: self.phred64,
                    keep_names: self.keep_names,
                    batch_bytes: self.batch_bytes,
                    batch_records: self.batch_records,
                    pool: self.pool.clone(),
                    names: self.names.clone().into(),
                };
                let readers = std::mem::take(&mut self.readers);
                match Par::start(readers, self.records.clone(), cfg, workers) {
                    Ok(p) => self.par = Some(p),
                    Err(readers) => self.readers = readers,
                }
            }
            if let Some(par) = &mut self.par {
                match par.next() {
                    Next::Batch(r) => {
                        let b = r?;
                        self.spots += b.len() as u64;
                        return Ok(Some(b));
                    }
                    Next::Sequential(tail) => {
                        self.par = None;
                        self.readers = tail.readers;
                        self.records = tail.records;
                    }
                }
            }
        }
        self.fill()
    }

    fn fill(&mut self) -> Result<Option<Batch<Segments>>, Error> {
        let mut b = BatchBuilder::with_segments(&self.pool, self.k, self.hint);
        let mut used = 0usize;
        while b.len() < self.batch_records && (b.is_empty() || used < self.batch_bytes) {
            if self.held.is_none() {
                let (spots, full) = self.fast_round(&mut b, &mut used)?;
                if full {
                    break;
                }
                if spots > 0 {
                    continue;
                }
            }
            match self.slow_spot(&mut b, &mut used)? {
                Slow::Pushed => {}
                Slow::Held | Slow::End => break,
            }
        }
        if b.is_empty() {
            return Ok(None);
        }
        self.spots += b.len() as u64;
        self.hint = self.hint.max(b.len());
        Ok(Some(b.finish()))
    }

    /// Reader index of segment `s`.
    fn reader_of(&self, s: usize) -> usize {
        if self.interleaved { 0 } else { s }
    }

    fn segment_label(&self, s: usize) -> Option<usize> {
        (self.k > 1).then_some(s)
    }

    fn input_error(&self, kind: ErrorKind, record: u64, s: usize) -> Error {
        Error::new(kind, "").at(
            Some(record),
            self.segment_label(s),
            self.names[self.reader_of(s)].as_ref(),
        )
    }

    fn io_error(&self, e: std::io::Error, s: usize) -> Error {
        let r = self.reader_of(s);
        Error::with_cause(ErrorKind::Io, e).at(
            Some(self.records[r] + 1),
            self.segment_label(s),
            self.names[r].as_ref(),
        )
    }

    /// Parses every spot that lies complete in the readers' current buffers,
    /// then consumes them. Returns the spot count and whether the byte budget
    /// stopped it.
    fn fast_round(
        &mut self,
        b: &mut BatchBuilder<Segments>,
        used: &mut usize,
    ) -> Result<(usize, bool), Error> {
        let n = self.readers.len();
        let mut bufs: Vec<&[u8]> = Vec::with_capacity(n);
        let mut io_failure = None;
        for (r, reader) in self.readers.iter_mut().enumerate() {
            match reader.fill_buf() {
                Ok(buf) => bufs.push(buf),
                Err(e) => {
                    io_failure = Some((e, r));
                    break;
                }
            }
        }
        if let Some((e, r)) = io_failure {
            drop(bufs);
            return Err(self.io_error(e, r));
        }
        // Positions after the last accepted spot, and the next spot's.
        let mut pos = vec![0usize; n];
        let mut next = vec![0usize; n];
        let mut count = vec![0u64; n];
        let mut fields: Vec<Lines<'_>> = vec![[&[][..]; 4]; self.k];
        let (k, offset) = (self.k, if self.phred64 { 64 } else { 33 });
        let mut spots = 0;
        let mut full = false;
        let mut failure = None;
        'spots: while b.len() < self.batch_records {
            next.copy_from_slice(&pos);
            for (s, f) in fields.iter_mut().enumerate() {
                let r = if self.interleaved { 0 } else { s };
                // Leave anything that is not a header (a blank line, a
                // record crossing the buffer end) to the slow path.
                match split_record(bufs[r], next[r]) {
                    Some((lines, end)) => {
                        *f = lines;
                        next[r] = end;
                    }
                    None => break 'spots,
                }
            }
            let bytes: usize = next.iter().zip(&pos).map(|(a, b)| a - b).sum();
            if *used + bytes > self.batch_bytes && !b.is_empty() {
                full = true;
                break;
            }
            for (s, f) in fields.iter().enumerate() {
                let r = if self.interleaved { 0 } else { s };
                count[r] += 1;
                if let Err(kind) = check(f, offset) {
                    failure = Some((kind, self.records[r] + count[r], s));
                    break 'spots;
                }
            }
            if let Err(e) = push_spot(
                b,
                &fields,
                self.keep_names,
                self.phred64.then_some(&mut self.qual33[..k]),
            ) {
                return Err(batch_error(e));
            }
            pos.copy_from_slice(&next);
            *used += bytes;
            spots += 1;
        }
        drop(bufs);
        if let Some((kind, record, s)) = failure {
            return Err(self.input_error(kind, record, s));
        }
        for (r, reader) in self.readers.iter_mut().enumerate() {
            reader.consume(pos[r]);
        }
        let stride = if self.interleaved { k as u64 } else { 1 };
        for rec in &mut self.records {
            *rec += spots as u64 * stride;
        }
        Ok((spots, full))
    }

    /// Reads one spot line by line (records crossing a buffer boundary,
    /// blank lines, end of input), or pushes the held one.
    ///
    /// The byte budget applies to the whole spot before its checks, as on
    /// the fast path: a spot that would take a non-empty batch over it is
    /// held for the next batch along with the outcome of its checks, so
    /// batches do not depend on where input buffers end. The first failure
    /// in segment order wins: a failed check of a segment comes before a
    /// read error (truncation, I/O) in a later one, and both before a
    /// record-count error.
    fn slow_spot(
        &mut self,
        b: &mut BatchBuilder<Segments>,
        used: &mut usize,
    ) -> Result<Slow, Error> {
        let bytes = match self.held.take() {
            Some(Held { error: Some(e), .. }) => return Err(e),
            Some(Held { bytes, .. }) => bytes,
            None => {
                let (bytes, error) = match self.read_spot()? {
                    Some(r) => r,
                    None => return Ok(Slow::End),
                };
                if !b.is_empty() && *used + bytes > self.batch_bytes {
                    self.held = Some(Held { bytes, error });
                    return Ok(Slow::Held);
                }
                if let Some(e) = error {
                    return Err(e);
                }
                bytes
            }
        };
        let k = self.k;
        let fields: Vec<Lines<'_>> = self.slow.iter().map(|(t, e)| lines_of(t, e)).collect();
        *used += bytes;
        push_spot(
            b,
            &fields,
            self.keep_names,
            self.phred64.then_some(&mut self.qual33[..k]),
        )
        .map_err(batch_error)?;
        Ok(Slow::Pushed)
    }

    /// Reads one spot's records into `slow`: `None` at the end of input,
    /// else its text bytes and the first failed check or record-count
    /// error. A read error is returned as `Err`, unless an earlier segment
    /// failed its check, which then is.
    fn read_spot(&mut self) -> Result<Option<(usize, Option<Error>)>, Error> {
        let offset = if self.phred64 { 64 } else { 33 };
        let (mut present, mut missing, mut bytes) = (0, None, 0);
        let mut failed: Option<Error> = None;
        for s in 0..self.k {
            let r = self.reader_of(s);
            let record = self.records[r] + 1;
            let (text, ends) = &mut self.slow[s];
            let got = read_record(&mut *self.readers[r], text, ends)
                .map_err(|e| match e {
                    SlowError::Io(e) => Error::with_cause(ErrorKind::Io, e),
                    SlowError::Kind(kind) => Error::new(kind, ""),
                })
                .map_err(|e| {
                    e.at(
                        Some(record),
                        (self.k > 1).then_some(s),
                        self.names[r].as_ref(),
                    )
                });
            let got = match (got, failed.take()) {
                (Ok(got), f) => {
                    failed = f;
                    got
                }
                (Err(_), Some(f)) => return Err(f),
                (Err(e), None) => return Err(e),
            };
            if !got {
                if self.interleaved {
                    if s == 0 {
                        return Ok(None);
                    }
                    let e = self.input_error(ErrorKind::Interleaved, record, s);
                    return Err(failed.unwrap_or(e));
                }
                missing.get_or_insert(s);
                continue;
            }
            self.records[r] = record;
            present += 1;
            let (text, ends) = &self.slow[s];
            bytes += text.len();
            if failed.is_none()
                && let Err(kind) = check(&lines_of(text, ends), offset)
            {
                failed = Some(self.input_error(kind, record, s));
            }
        }
        if present == 0 {
            return Ok(None);
        }
        if let (None, Some(s)) = (&failed, missing) {
            let r = self.reader_of(s);
            failed = Some(self.input_error(ErrorKind::MateCount, self.records[r] + 1, s));
        }
        Ok(Some((bytes, failed)))
    }
}

/// A slow-path spot read while filling a batch it did not fit; its records
/// are in `slow`.
struct Held {
    bytes: usize,
    /// Its first failure, reported when the next batch is filled.
    error: Option<Error>,
}

enum Slow {
    Pushed,
    Held,
    End,
}

pub(crate) fn batch_error(e: readkit_batch::Error) -> Error {
    Error::with_cause(ErrorKind::Batch, e)
}

/// Splits the record starting at `pos` if it begins with `@` and its four
/// lines all end in `buf`. Returns the lines and the position after it.
#[inline]
pub(crate) fn split_record(buf: &[u8], pos: usize) -> Option<(Lines<'_>, usize)> {
    if buf.get(pos) != Some(&b'@') {
        return None;
    }
    let mut start = pos;
    let mut lines: Lines<'_> = [&[]; 4];
    for line in &mut lines {
        let end = start + memchr::memchr(b'\n', &buf[start..])?;
        *line = strip_cr(&buf[start..end]);
        start = end + 1;
    }
    Some((lines, start))
}

#[inline]
fn strip_cr(line: &[u8]) -> &[u8] {
    line.strip_suffix(b"\r").unwrap_or(line)
}

/// The four lines of a slow-path record: `text` ends at `ends[i]` for line
/// `i`, each line with its `\n` (except possibly the last) and `\r`.
fn lines_of<'a>(text: &'a [u8], ends: &[usize; 4]) -> Lines<'a> {
    let mut start = 0;
    ends.map(|end| {
        let line = &text[start..end];
        start = end;
        strip_cr(line.strip_suffix(b"\n").unwrap_or(line))
    })
}

/// The record checks of fasterp's `Read::validate`, in its order.
#[inline]
pub(crate) fn check(lines: &Lines<'_>, offset: u8) -> Result<(), ErrorKind> {
    let [name, seq, plus, qual] = lines;
    if !name.starts_with(b"@") {
        return Err(ErrorKind::Header);
    }
    if !plus.starts_with(b"+") {
        return Err(ErrorKind::Separator);
    }
    if seq.len() != qual.len() {
        return Err(ErrorKind::LengthMismatch);
    }
    if !valid_quality(qual, offset) {
        return Err(ErrorKind::Quality);
    }
    Ok(())
}

/// Whether every byte lies in `offset..=126`. A branch-free reduction, so
/// the compiler vectorizes it; wrapping subtraction rejects bytes below
/// `offset`.
#[inline]
fn valid_quality(qual: &[u8], offset: u8) -> bool {
    qual.iter().fold(0u8, |bad, &q| {
        bad | u8::from(q.wrapping_sub(offset) > 126 - offset)
    }) == 0
}

/// Pushes one validated spot. `qual33`, given for Phred+64 input, holds one
/// scratch buffer per segment for the converted qualities.
pub(crate) fn push_spot(
    b: &mut BatchBuilder<Segments>,
    fields: &[Lines<'_>],
    keep_names: bool,
    qual33: Option<&mut [Vec<u8>]>,
) -> Result<(), readkit_batch::Error> {
    let quals: Option<&[Vec<u8>]> = match qual33 {
        Some(bufs) => {
            for (buf, f) in bufs.iter_mut().zip(fields) {
                buf.clear();
                buf.extend(f[3].iter().map(|q| q - 31));
            }
            Some(bufs)
        }
        None => None,
    };
    let read = |s: usize| {
        let [name, seq, _, qual] = fields[s];
        let qual = quals.map_or(qual, |q| &q[s][..]);
        if keep_names {
            Read::named(&name[1..], seq, qual)
        } else {
            Read::new(seq, qual)
        }
    };
    let k = fields.len();
    if k <= INLINE_SEGMENTS {
        let mut reads = [Read::new(b"", b""); INLINE_SEGMENTS];
        for (s, slot) in reads[..k].iter_mut().enumerate() {
            *slot = read(s);
        }
        b.push(&reads[..k])
    } else {
        let reads: Vec<Read<'_>> = (0..k).map(read).collect();
        b.push(&reads)
    }
}

enum SlowError {
    Io(std::io::Error),
    Kind(ErrorKind),
}

impl From<std::io::Error> for SlowError {
    fn from(e: std::io::Error) -> Self {
        SlowError::Io(e)
    }
}

/// Whitespace other than LF, as C `isspace`.
fn is_space(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\r' | b'\x0b' | b'\x0c')
}

/// fastp skips lines that cannot start a record. Skip only blank and
/// whitespace-only ones, and only where a header is expected (fasterp's
/// `Reader::skip_blank_lines`).
fn skip_blank_lines(reader: &mut dyn BufRead) -> Result<(), SlowError> {
    loop {
        let buffer = reader.fill_buf()?;
        if buffer.first().is_none_or(|&b| b == b'@') {
            return Ok(());
        }
        match buffer.iter().position(|&b| !is_space(b)) {
            Some(i) if buffer[i] == b'\n' => reader.consume(i + 1),
            Some(0) => return Ok(()),
            Some(_) => return Err(SlowError::Kind(ErrorKind::Header)),
            None => {
                // A whitespace run reaching the buffer end: read the rest of
                // the line to see whether it is blank.
                let mut line = Vec::new();
                reader.read_until(b'\n', &mut line)?;
                if !line.iter().all(|&b| is_space(b) || b == b'\n') {
                    return Err(SlowError::Kind(ErrorKind::Header));
                }
            }
        }
    }
}

/// Reads one record's four lines into `text` (fasterp's `Records::read_one`).
/// Returns false at the end of input; a final line without `\n` is kept.
fn read_record(
    reader: &mut dyn BufRead,
    text: &mut Vec<u8>,
    ends: &mut [usize; 4],
) -> Result<bool, SlowError> {
    skip_blank_lines(reader)?;
    text.clear();
    for (line, end) in ends.iter_mut().enumerate() {
        let n = reader.read_until(b'\n', text)?;
        if n == 0 {
            if line == 0 {
                return Ok(false);
            }
            return Err(SlowError::Kind(ErrorKind::Truncated));
        }
        *end = text.len();
    }
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quality_validation_matches_scalar_at_every_byte_and_vector_boundary() {
        for offset in [33, 64] {
            for len in [0, 1, 15, 16, 17, 31, 32, 33, 51, 64, 65, 148, 513] {
                let mut quality = vec![offset; len];
                assert!(valid_quality(&quality, offset));
                for i in 0..len {
                    for byte in 0..=255 {
                        quality[i] = byte;
                        assert_eq!(
                            valid_quality(&quality, offset),
                            byte >= offset && byte <= 126,
                            "offset={offset}, len={len}, index={i}, byte={byte}"
                        );
                    }
                    quality[i] = offset;
                }
            }
        }
    }
}
