//! Parallel gzip encoding into independent members.

use crate::error::{Error, ErrorKind};
use crate::fqdeflate;
use crate::header::MEMBER_SIZE_ID;
use crossbeam_channel::{Receiver, Sender};
use libdeflater::{CompressionLvl, Compressor};
use std::collections::BTreeMap;
use std::io::{self, Write};

/// When [`GzWriter`] uses the FASTQ-aware encoder for a block instead of
/// libdeflate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Fastq {
    /// fasterp's rule: at levels 2 to 4, for blocks whose first records have
    /// binned qualities and are not sorted by position
    /// ([`binned_qualities`](crate::binned_qualities)). There the FASTQ
    /// encoder is as small as libdeflate and about 2.4 times faster.
    #[default]
    Auto,
    /// Every block, whatever the level and contents. Output is valid gzip
    /// for any input, but larger than libdeflate's on data that is not
    /// binned FASTQ.
    Always,
    /// Never: libdeflate only.
    Never,
}

/// Options for [`GzWriter`].
///
/// Construct with `EncodeOptions { level: 4, ..Default::default() }`.
///
/// Memory: at most `in_flight` blocks of input and their compressed members
/// are held at once, plus the block being filled.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EncodeOptions {
    /// libdeflate level, 0 (stored) to 12 (default 6).
    pub level: u32,
    /// Compression threads (default: half the available cores, clamped to
    /// 1..=8, as fasterp uses for one output; divide it among outputs that
    /// are written at once). 0 compresses on the calling thread. Output does
    /// not depend on this.
    pub threads: usize,
    /// Input bytes per gzip member (default 256 KiB, at most 1 GiB).
    pub block_size: usize,
    /// When to use the FASTQ-aware encoder.
    pub fastq: Fastq,
    /// Blocks compressed or waiting to be written at once. `None` means
    /// `2 × threads`.
    pub in_flight: Option<usize>,
}

impl Default for EncodeOptions {
    fn default() -> Self {
        Self {
            level: 6,
            threads: std::thread::available_parallelism()
                .map_or(1, |n| n.get())
                .div_ceil(2)
                .clamp(1, 8),
            block_size: 256 << 10,
            fastq: Fastq::Auto,
            in_flight: None,
        }
    }
}

/// Add the member-size subfield to a gzip member with no extra field.
fn with_member_size(member: Vec<u8>) -> Vec<u8> {
    const FEXTRA: u8 = 4;
    if member.len() < 18 || member[3] & !0x08 != 0 {
        return member;
    }
    let total = (member.len() + 10) as u32;
    let mut out = Vec::with_capacity(member.len() + 10);
    out.extend_from_slice(&member[..3]);
    out.push(member[3] | FEXTRA);
    out.extend_from_slice(&member[4..10]);
    out.extend_from_slice(&8u16.to_le_bytes());
    out.extend_from_slice(&MEMBER_SIZE_ID);
    out.extend_from_slice(&4u16.to_le_bytes());
    out.extend_from_slice(&total.to_le_bytes());
    out.extend_from_slice(&member[10..]);
    out
}

/// Compresses blocks into members; one per thread.
struct Block {
    level: u32,
    fastq: Fastq,
    compressor: Compressor,
}

impl Block {
    fn new(level: u32, fastq: Fastq) -> Self {
        let lvl = CompressionLvl::new(level as i32).expect("level checked when the writer opened");
        Self {
            level,
            fastq,
            compressor: Compressor::new(lvl),
        }
    }

    fn compress(&mut self, input: &[u8]) -> Result<Vec<u8>, Error> {
        let fastq = match self.fastq {
            Fastq::Always => true,
            Fastq::Never => false,
            Fastq::Auto => (2..=4).contains(&self.level) && fqdeflate::suits(input),
        };
        if fastq {
            let mut output = Vec::new();
            fqdeflate::gzip_member(input, &mut output);
            return Ok(output);
        }
        let mut output = vec![0; self.compressor.gzip_compress_bound(input.len())];
        let size = self
            .compressor
            .gzip_compress(input, &mut output)
            .map_err(|e| Error::internal(format!("libdeflate compression failed: {e:?}")))?;
        output.truncate(size);
        Ok(with_member_size(output))
    }
}

struct Workers {
    jobs: Option<Sender<(u64, Vec<u8>)>>,
    results: Receiver<(u64, Result<Vec<u8>, Error>)>,
    handles: Vec<std::thread::JoinHandle<()>>,
}

impl Drop for Workers {
    fn drop(&mut self) {
        drop(self.jobs.take());
        for h in self.handles.drain(..) {
            let _ = h.join();
        }
    }
}

/// gzip encoder writing independent members of
/// [`block_size`](EncodeOptions::block_size) input bytes each, compressed
/// on worker threads and written in order to the inner writer.
///
/// Each member carries fasterp's `FP` extra subfield with its compressed
/// size, so [`ParGzReader`](crate::ParGzReader) decodes the output member
/// by member in parallel; every other gzip reader skips the field. The
/// output depends only on the input, the options other than `threads`, and
/// where [`flush`](Write::flush) was called (a flush ends the current member).
///
/// The inner writer is written on the calling thread, inside `write`,
/// `flush` and [`finish`](Self::finish), so it need not be `Send`. Empty
/// input still produces one (empty) member.
///
/// Call [`finish`](Self::finish) to see errors; dropping the writer finishes
/// it and ignores them.
///
/// ```
/// use std::io::{Read, Write};
/// let mut writer = readkit_gz::GzWriter::new(Vec::new(), &Default::default())?;
/// writer.write_all(b"@r1\nACGT\n+\nFFFF\n")?;
/// let gz = writer.finish()?;
/// let mut text = String::new();
/// readkit_gz::GzReader::from_reader(&gz[..]).read_to_string(&mut text)?;
/// assert_eq!(text, "@r1\nACGT\n+\nFFFF\n");
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub struct GzWriter<W: Write> {
    inner: Option<W>,
    block: Vec<u8>,
    block_size: usize,
    /// Compression on the calling thread (`threads == 0`).
    local: Option<Block>,
    workers: Option<Workers>,
    in_flight: usize,
    next_submit: u64,
    next_write: u64,
    done: BTreeMap<u64, Vec<u8>>,
    members: u64,
    error: Option<Error>,
}

impl<W: Write> std::fmt::Debug for GzWriter<W> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GzWriter")
            .field("members", &self.members)
            .finish_non_exhaustive()
    }
}

impl<W: Write> GzWriter<W> {
    /// A writer compressing into `inner`.
    pub fn new(inner: W, options: &EncodeOptions) -> Result<Self, Error> {
        let bad = |what: &str| Err(Error::new(ErrorKind::InvalidInput, what.to_string()));
        if options.level > 12 {
            return bad("level must be 0 to 12");
        }
        if options.block_size == 0 || options.block_size > 1 << 30 {
            return bad("block_size must be 1 byte to 1 GiB");
        }
        if options.in_flight == Some(0) {
            return bad("in_flight must be at least 1");
        }
        let (level, fastq) = (options.level, options.fastq);
        let (local, workers) = if options.threads == 0 {
            (Some(Block::new(level, fastq)), None)
        } else {
            let (jobs, queue) = crossbeam_channel::unbounded::<(u64, Vec<u8>)>();
            let (done, results) = crossbeam_channel::unbounded();
            let mut handles = Vec::new();
            for _ in 0..options.threads {
                let (queue, done) = (queue.clone(), done.clone());
                let handle = std::thread::Builder::new()
                    .name("readkit-gz-encode".into())
                    .spawn(move || {
                        let mut block = Block::new(level, fastq);
                        for (seq, input) in queue {
                            if done.send((seq, block.compress(&input))).is_err() {
                                break;
                            }
                        }
                    })
                    .map_err(|e| {
                        Error::new(ErrorKind::Io, "cannot start a compression thread")
                            .with_source(e)
                    })?;
                handles.push(handle);
            }
            let workers = Workers {
                jobs: Some(jobs),
                results,
                handles,
            };
            (None, Some(workers))
        };
        Ok(Self {
            inner: Some(inner),
            block: Vec::with_capacity(options.block_size),
            block_size: options.block_size,
            local,
            workers,
            in_flight: options.in_flight.unwrap_or(2 * options.threads).max(1),
            next_submit: 0,
            next_write: 0,
            done: BTreeMap::new(),
            members: 0,
            error: None,
        })
    }

    /// The inner writer.
    pub fn get_ref(&self) -> &W {
        self.inner.as_ref().expect("present until finish")
    }

    /// Members written so far.
    pub fn members(&self) -> u64 {
        self.members
    }

    fn put(&mut self, member: &[u8]) -> Result<(), Error> {
        let inner = self.inner.as_mut().expect("present until finish");
        inner.write_all(member).map_err(|e| {
            Error::new(ErrorKind::Io, format!("writing gzip output: {e}")).with_source(e)
        })?;
        self.members += 1;
        Ok(())
    }

    /// Write finished members in order; with `wait`, block until the oldest
    /// outstanding one is written.
    fn drain(&mut self, wait: bool) -> Result<(), Error> {
        let Some(workers) = &self.workers else {
            return Ok(());
        };
        let results = workers.results.clone();
        let mut need = wait && self.next_write < self.next_submit;
        loop {
            let got = if need {
                Some(
                    results
                        .recv()
                        .map_err(|_| Error::internal("compression workers stopped"))?,
                )
            } else {
                results.try_recv().ok()
            };
            let Some((seq, result)) = got else {
                return Ok(());
            };
            self.done.insert(seq, result?);
            while let Some(member) = self.done.remove(&self.next_write) {
                self.put(&member)?;
                self.next_write += 1;
                need = false;
            }
        }
    }

    fn submit(&mut self, input: Vec<u8>) -> Result<(), Error> {
        if let Some(block) = &mut self.local {
            let member = block.compress(&input)?;
            self.next_submit += 1;
            self.next_write += 1;
            return self.put(&member);
        }
        while self.next_submit - self.next_write >= self.in_flight as u64 {
            self.drain(true)?;
        }
        let workers = self.workers.as_ref().expect("workers when not local");
        workers
            .jobs
            .as_ref()
            .expect("present until drop")
            .send((self.next_submit, input))
            .map_err(|_| Error::internal("compression workers stopped"))?;
        self.next_submit += 1;
        self.drain(false)
    }

    /// Compress the partial block and write everything submitted.
    fn complete(&mut self, at_end: bool) -> Result<(), Error> {
        if !self.block.is_empty() || (at_end && self.next_submit == 0) {
            let block = std::mem::replace(&mut self.block, Vec::with_capacity(self.block_size));
            self.submit(block)?;
        }
        while self.next_write < self.next_submit {
            self.drain(true)?;
        }
        Ok(())
    }

    fn guard<T>(&mut self, f: impl FnOnce(&mut Self) -> Result<T, Error>) -> Result<T, Error> {
        if let Some(e) = &self.error {
            return Err(e.clone());
        }
        f(self).inspect_err(|e| self.error = Some(e.clone()))
    }

    /// Write the last member, flush, stop the workers and return the inner
    /// writer.
    pub fn finish(mut self) -> Result<W, Error> {
        self.guard(|w| w.complete(true))?;
        let mut inner = self.inner.take().expect("present until finish");
        inner.flush().map_err(|e| {
            Error::new(ErrorKind::Io, format!("flushing gzip output: {e}")).with_source(e)
        })?;
        Ok(inner)
    }
}

impl<W: Write> Write for GzWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.guard(|w| {
            let mut rest = buf;
            while !rest.is_empty() {
                let take = rest.len().min(w.block_size - w.block.len());
                w.block.extend_from_slice(&rest[..take]);
                rest = &rest[take..];
                if w.block.len() == w.block_size {
                    let block = std::mem::replace(&mut w.block, Vec::with_capacity(w.block_size));
                    w.submit(block)?;
                }
            }
            Ok(buf.len())
        })
        .map_err(io::Error::from)
    }

    /// Ends the current member, writes every member and flushes the inner
    /// writer.
    fn flush(&mut self) -> io::Result<()> {
        self.guard(|w| w.complete(false))?;
        self.inner.as_mut().expect("present until finish").flush()
    }
}

impl<W: Write> Drop for GzWriter<W> {
    fn drop(&mut self) {
        if self.inner.is_some() && self.error.is_none() {
            let _ = self.complete(true);
            if let Some(inner) = &mut self.inner {
                let _ = inner.flush();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    #[test]
    fn parallel_members_roundtrip_levels_boundaries_and_empty_output() {
        const BLOCK: usize = 256 << 10;
        let data: Vec<_> = (0..BLOCK * 3 + 31)
            .map(|i| ((i * 71 + i / 239) % 256) as u8)
            .collect();
        for level in [1, 4, 9] {
            for len in [0, 1, BLOCK - 1, BLOCK, data.len()] {
                let options = EncodeOptions {
                    level,
                    threads: 2,
                    ..Default::default()
                };
                let mut writer = GzWriter::new(Vec::new(), &options).unwrap();
                for block in data[..len].chunks(27013) {
                    writer.write_all(block).unwrap();
                }
                writer.flush().unwrap();
                let bytes = writer.finish().unwrap();
                assert!(
                    bytes.starts_with(&[31, 139]),
                    "empty output must still be gzip"
                );
                let mut actual = Vec::new();
                flate2::read::MultiGzDecoder::new(&bytes[..])
                    .read_to_end(&mut actual)
                    .unwrap();
                assert_eq!(actual, data[..len]);
            }
        }
    }

    #[test]
    fn compressed_output_errors_are_returned_by_finish() {
        #[derive(Debug)]
        struct Broken;
        impl Write for Broken {
            fn write(&mut self, _: &[u8]) -> io::Result<usize> {
                Err(io::Error::new(io::ErrorKind::BrokenPipe, "closed output"))
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        for threads in [0, 2] {
            let options = EncodeOptions {
                threads,
                ..Default::default()
            };
            let mut writer = GzWriter::new(Broken, &options).unwrap();
            writer.write_all(b"some input").unwrap();
            let e = writer.finish().unwrap_err();
            assert_eq!(e.kind(), ErrorKind::Io);
        }
    }
}
