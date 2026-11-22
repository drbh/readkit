//! The FASTQ writer (fasterp's `Read::write`).

use std::fmt;
use std::io::Write;
use std::marker::PhantomData;

use readkit_batch::{Batch, Layout};
use readkit_io::BatchSink;

use crate::error::{Error, ErrorKind};

/// What a written record's third line holds.
///
/// Batches do not store the input's `+` line (see the crate docs), so the
/// sink chooses. fasterp copies the input's `+` line; for the two
/// conventions found in practice, bare `+` (current Illumina output) and
/// `+` repeating the name (older Illumina and SRA dumps), one of these
/// styles reproduces fasterp's output byte for byte.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum PlusLine {
    /// `+` alone.
    #[default]
    Bare,
    /// `+` followed by the record's name.
    Name,
}

/// Writes batches as FASTQ text.
///
/// [`new`](Self::new) writes every segment of a spot, in order, to one
/// writer (interleaved output for `k > 1`); [`split`](Self::split) writes
/// segment `s` to writer `s` (R1/R2 files, ...). Records come out as
///
/// ```text
/// @<name>
/// <sequence>
/// +[<name>]
/// <quality>
/// ```
///
/// with `\n` line endings. Batches without names get the 1-based spot
/// number, counted across all batches this sink has written, as name.
/// Qualities are written as stored (Phred+33 when read by
/// [`FastqSource`](crate::FastqSource), which converts Phred+64 input).
///
/// Each batch is formatted into a buffer per writer and written with one
/// `write_all` per writer, so `W` need not be buffered.
pub struct FastqSink<L, W: Write> {
    writers: Vec<W>,
    bufs: Vec<Vec<u8>>,
    plus: PlusLine,
    spots: u64,
    _layout: PhantomData<fn(L)>,
}

impl<L, W: Write> FastqSink<L, W> {
    /// Writes all segments to `writer`.
    pub fn new(writer: W) -> Self {
        Self::make(vec![writer])
    }

    /// Writes segment `s` of every spot to `writers[s]`. Batches must have
    /// `writers.len()` segments ([`ErrorKind::Options`] otherwise). With one
    /// writer this is [`new`](Self::new).
    pub fn split(writers: Vec<W>) -> Self {
        Self::make(writers)
    }

    fn make(writers: Vec<W>) -> Self {
        let bufs = writers.iter().map(|_| Vec::new()).collect();
        Self {
            writers,
            bufs,
            plus: PlusLine::Bare,
            spots: 0,
            _layout: PhantomData,
        }
    }

    /// Sets the `+` line style (default [`PlusLine::Bare`]).
    pub fn with_plus_line(mut self, plus: PlusLine) -> Self {
        self.plus = plus;
        self
    }

    /// Spots written so far.
    pub fn spots(&self) -> u64 {
        self.spots
    }

    /// Flushes and returns the writers.
    pub fn into_inner(mut self) -> Result<Vec<W>, Error> {
        self.flush()?;
        Ok(self.writers)
    }

    fn flush(&mut self) -> Result<(), Error> {
        for w in &mut self.writers {
            w.flush().map_err(|e| Error::with_cause(ErrorKind::Io, e))?;
        }
        Ok(())
    }

    fn format<M: Layout>(&mut self, batch: &Batch<M>) -> Result<(), Error> {
        let k = batch.segments();
        let split = self.writers.len() > 1;
        if split && self.writers.len() != k {
            return Err(Error::new(
                ErrorKind::Options,
                format!("{k} segments for {} writers", self.writers.len()),
            ));
        }
        if self.writers.is_empty() {
            return Err(Error::new(ErrorKind::Options, "no writers"));
        }
        let columns: Vec<_> = batch.columns().collect();
        for (s, c) in columns.iter().enumerate() {
            let buf = &mut self.bufs[if split { s } else { 0 }];
            let seqs = c.seqs().data().len();
            let names = c.names().map_or(batch.len() * 12, |n| n.data().len());
            let plus = if self.plus == PlusLine::Name {
                names
            } else {
                0
            };
            buf.reserve(2 * seqs + names + plus + 5 * batch.len());
        }
        let mut number = itoa::Buffer::default();
        for i in 0..batch.len() {
            let spot = self.spots + i as u64 + 1;
            for (s, c) in columns.iter().enumerate() {
                let buf = &mut self.bufs[if split { s } else { 0 }];
                let name: &[u8] = match c.names() {
                    Some(n) => n.get(i),
                    None => number.format(spot),
                };
                let r = c.read(i);
                buf.push(b'@');
                buf.extend_from_slice(name);
                buf.push(b'\n');
                buf.extend_from_slice(r.seq());
                buf.extend_from_slice(b"\n+");
                if self.plus == PlusLine::Name {
                    buf.extend_from_slice(name);
                }
                buf.push(b'\n');
                buf.extend_from_slice(r.qual());
                buf.push(b'\n');
            }
        }
        self.spots += batch.len() as u64;
        Ok(())
    }
}

impl<L: Layout, W: Write> BatchSink<L> for FastqSink<L, W> {
    type Error = Error;

    fn accept(&mut self, batch: Batch<L>) -> Result<(), Error> {
        self.format(&batch)?;
        drop(batch);
        for (w, buf) in self.writers.iter_mut().zip(&mut self.bufs) {
            let r = w.write_all(buf);
            buf.clear();
            r.map_err(|e| Error::with_cause(ErrorKind::Io, e))?;
        }
        Ok(())
    }

    fn finish(mut self) -> Result<(), Error> {
        self.flush()
    }
}

impl<L, W: Write> fmt::Debug for FastqSink<L, W> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FastqSink")
            .field("writers", &self.writers.len())
            .field("plus", &self.plus)
            .field("spots", &self.spots)
            .finish()
    }
}

/// Decimal formatting of spot numbers without allocating.
mod itoa {
    #[derive(Default)]
    pub(super) struct Buffer([u8; 20]);

    impl Buffer {
        pub(super) fn format(&mut self, mut n: u64) -> &[u8] {
            let mut i = self.0.len();
            loop {
                i -= 1;
                self.0[i] = b'0' + (n % 10) as u8;
                n /= 10;
                if n == 0 {
                    break;
                }
            }
            &self.0[i..]
        }
    }
}
