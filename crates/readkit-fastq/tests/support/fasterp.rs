//! fasterp's strict FASTQ reader and writer, vendored for tests and
//! benchmarks as the reference behaviour (fasterp commit 84ba478, "fix: skip
//! blank lines between fastq records"; `src/compat.rs` `Read`, `Reader`,
//! `is_space`, `valid_quality` and `src/compat/block.rs`; MIT).
//!
//! Changes: `anyhow` replaced by `String` errors with the same messages;
//! the reader is generic-free (`Box<dyn BufRead + Send>`, as in fasterp);
//! trimming and statistics fields dropped; the sample batch variant dropped.
#![allow(dead_code)]

use std::borrow::Cow;
use std::io::{BufRead, Write};

pub type Result<T> = std::result::Result<T, String>;

macro_rules! ensure {
    ($cond:expr, $($arg:tt)*) => {
        if !$cond {
            return Err(format!($($arg)*));
        }
    };
}

fn io(e: std::io::Error) -> String {
    e.to_string()
}

#[derive(Clone, Debug, Default)]
pub struct Read<'a> {
    pub name: Cow<'a, [u8]>,
    pub seq: Cow<'a, [u8]>,
    pub plus: Cow<'a, [u8]>,
    pub qual: Cow<'a, [u8]>,
}

impl Read<'_> {
    pub fn validate(&mut self, record: usize, phred64: bool) -> Result<()> {
        let [name, seq, plus, qual] = [&self.name, &self.seq, &self.plus, &self.qual];
        ensure!(
            name.starts_with(b"@"),
            "Invalid FASTQ header at record {}",
            record
        );
        ensure!(
            plus.starts_with(b"+"),
            "Invalid FASTQ separator at record {}",
            record
        );
        ensure!(
            seq.len() == qual.len(),
            "Sequence/quality length mismatch at record {}",
            record
        );
        let offset = if phred64 { 64 } else { 33 };
        ensure!(
            valid_quality(qual, offset),
            "Invalid quality encoding at record {}",
            record
        );
        if phred64 {
            self.qual.to_mut().iter_mut().for_each(|q| *q -= 31);
        }
        Ok(())
    }

    pub fn write(&self, writer: &mut impl Write) -> Result<()> {
        writer.write_all(&self.name).map_err(io)?;
        for line in [&b""[..], &self.seq, &self.plus, &self.qual] {
            writer.write_all(line).map_err(io)?;
            writer.write_all(b"\n").map_err(io)?;
        }
        Ok(())
    }
}

pub struct Reader {
    pub inner: Box<dyn BufRead + Send>,
    pub record: usize,
    pub phred64: bool,
}

impl Reader {
    pub fn new(inner: Box<dyn BufRead + Send>, phred64: bool) -> Self {
        Self {
            inner,
            record: 0,
            phred64,
        }
    }

    pub fn next(&mut self) -> Result<Option<Read<'static>>> {
        let mut read = Read::default();
        Ok(self.read_into(&mut read)?.then_some(read))
    }

    fn skip_blank_lines(&mut self) -> Result<()> {
        loop {
            let buffer = self.inner.fill_buf().map_err(io)?;
            if buffer.first().is_none_or(|&b| b == b'@') {
                return Ok(());
            }
            match buffer.iter().position(|&b| !is_space(b)) {
                Some(i) if buffer[i] == b'\n' => self.inner.consume(i + 1),
                Some(0) => return Ok(()),
                Some(_) => {
                    return Err(format!(
                        "Invalid FASTQ header at record {}",
                        self.record + 1
                    ));
                }
                None => {
                    let mut line = Vec::new();
                    self.inner.read_until(b'\n', &mut line).map_err(io)?;
                    ensure!(
                        line.iter().all(|&b| is_space(b) || b == b'\n'),
                        "Invalid FASTQ header at record {}",
                        self.record + 1
                    );
                }
            }
        }
    }

    pub fn read_into(&mut self, read: &mut Read) -> Result<bool> {
        let mut lines = [
            read.name.to_mut(),
            read.seq.to_mut(),
            read.plus.to_mut(),
            read.qual.to_mut(),
        ];
        let mut buffer = self.inner.fill_buf().map_err(io)?;
        if buffer.first() != Some(&b'@') {
            self.skip_blank_lines()?;
            buffer = self.inner.fill_buf().map_err(io)?;
        }
        let mut newlines = memchr::memchr_iter(b'\n', buffer);
        let ends = std::array::from_fn::<_, 4, _>(|_| newlines.next());
        if let [Some(a), Some(b), Some(c), Some(d)] = ends {
            let mut start = 0;
            for (line, end) in lines.iter_mut().zip([a, b, c, d]) {
                let content = &buffer[start..end];
                line.clear();
                line.extend_from_slice(content.strip_suffix(b"\r").unwrap_or(content));
                start = end + 1;
            }
            self.inner.consume(d + 1);
        } else {
            for (i, line) in lines.iter_mut().enumerate() {
                line.clear();
                let n = self.inner.read_until(b'\n', line).map_err(io)?;
                if n == 0 && i == 0 {
                    return Ok(false);
                }
                ensure!(n > 0, "Truncated FASTQ record {}", self.record + 1);
                if line.last() == Some(&b'\n') {
                    line.pop();
                }
                if line.last() == Some(&b'\r') {
                    line.pop();
                }
            }
        }
        self.record += 1;
        read.validate(self.record, self.phred64)?;
        Ok(true)
    }
}

fn is_space(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\r' | b'\x0b' | b'\x0c')
}

pub fn valid_quality(qual: &[u8], offset: u8) -> bool {
    qual.iter().fold(0u8, |bad, &q| {
        bad | u8::from(q.wrapping_sub(offset) > 126 - offset)
    }) == 0
}

/// fasterp's owned `BatchReader` loop (one record or pair at a time).
pub fn read_pair(
    first: &mut Reader,
    second: Option<&mut Reader>,
    interleaved: bool,
) -> Result<Option<(Read<'static>, Option<Read<'static>>)>> {
    let Some(r1) = first.next()? else {
        if let Some(r2) = second {
            ensure!(
                r2.next()?.is_none(),
                "Paired input files have different read counts"
            );
        }
        return Ok(None);
    };
    let r2 = if interleaved {
        let Some(r) = first.next()? else {
            return Err("Odd read count in interleaved FASTQ".into());
        };
        Some(r)
    } else if let Some(reader) = second {
        let Some(r) = reader.next()? else {
            return Err("Paired input files have different read counts".into());
        };
        Some(r)
    } else {
        None
    };
    Ok(Some((r1, r2)))
}

// ---- compat/block.rs ----

#[derive(Default)]
pub struct Records {
    data: Vec<u8>,
    ends: Vec<[usize; 4]>,
}

impl Records {
    pub fn len(&self) -> usize {
        self.ends.len()
    }

    pub fn bytes(&self) -> usize {
        self.data.len()
    }

    pub fn read(&self, index: usize, ordinal: usize, phred64: bool) -> Result<Read<'_>> {
        let mut start = if index == 0 {
            0
        } else {
            self.ends[index - 1][3]
        };
        let fields = self.ends[index].map(|end| {
            let line = &self.data[start..end];
            start = end;
            let line = line.strip_suffix(b"\n").unwrap_or(line);
            Cow::Borrowed(line.strip_suffix(b"\r").unwrap_or(line))
        });
        let [name, seq, plus, qual] = fields;
        let mut read = Read {
            name,
            seq,
            plus,
            qual,
        };
        read.validate(ordinal + 1, phred64)?;
        Ok(read)
    }

    fn read_one(&mut self, reader: &mut Reader) -> Result<bool> {
        reader.skip_blank_lines()?;
        let mut ends = [0; 4];
        for (line, end) in ends.iter_mut().enumerate() {
            let n = reader.inner.read_until(b'\n', &mut self.data).map_err(io)?;
            if n == 0 && line == 0 {
                return Ok(false);
            }
            ensure!(n > 0, "Truncated FASTQ record {}", reader.record + 1);
            *end = self.data.len();
        }
        self.ends.push(ends);
        reader.record += 1;
        Ok(true)
    }
}

#[inline]
fn next_ends(newlines: &mut impl Iterator<Item = usize>) -> Option<[usize; 4]> {
    Some([
        newlines.next()? + 1,
        newlines.next()? + 1,
        newlines.next()? + 1,
        newlines.next()? + 1,
    ])
}

pub fn read_batch(
    first: &mut Reader,
    second: Option<&mut Reader>,
    budget: usize,
    limit: usize,
) -> Result<Option<[Records; 2]>> {
    read_batch_inner::<false>(first, second, budget, limit)
}

pub fn read_interleaved_batch(
    first: &mut Reader,
    budget: usize,
    pairs: usize,
) -> Result<Option<[Records; 2]>> {
    read_batch_inner::<true>(first, None, budget, pairs)
}

fn read_batch_inner<const INTERLEAVED: bool>(
    first: &mut Reader,
    mut second: Option<&mut Reader>,
    budget: usize,
    limit: usize,
) -> Result<Option<[Records; 2]>> {
    let mut blocks: [Records; 2] = Default::default();
    let stride = if INTERLEAVED { 2 } else { 1 };
    while blocks[0].len() / stride < limit
        && blocks.iter().map(|b| b.data.len()).sum::<usize>() < budget
    {
        let remaining = limit - blocks[0].len() / stride;
        let b1 = first.inner.fill_buf().map_err(io)?;
        let b2 = match second.as_mut() {
            Some(r) => Some(r.inner.fill_buf().map_err(io)?),
            None => None,
        };
        let mut newlines1 = memchr::memchr_iter(b'\n', b1);
        let mut newlines2 = b2.map(|b| memchr::memchr_iter(b'\n', b));
        let base = [blocks[0].data.len(), blocks[1].data.len()];
        let used = base[0] + base[1];
        let mut count = 0;
        let mut consumed = [0; 2];
        let mut full = false;
        while count < remaining {
            if b1.get(consumed[0]) != Some(&b'@')
                || b2.is_some_and(|b| b.get(consumed[1]) != Some(&b'@'))
            {
                break;
            }
            let Some(ends1) = next_ends(&mut newlines1) else {
                break;
            };
            let ends2 = if INTERLEAVED {
                if b1.get(ends1[3]) != Some(&b'@') {
                    break;
                }
                let Some(ends) = next_ends(&mut newlines1) else {
                    break;
                };
                Some(ends)
            } else if let Some(lines) = &mut newlines2 {
                let Some(ends) = next_ends(lines) else {
                    break;
                };
                Some(ends)
            } else {
                None
            };
            let next_bytes = if INTERLEAVED {
                ends2.unwrap()[3]
            } else {
                ends1[3] + ends2.map_or(0, |e| e[3])
            };
            if used + next_bytes > budget && blocks[0].len() > 0 {
                full = true;
                break;
            }
            consumed[0] = ends1[3];
            blocks[0].ends.push(ends1.map(|n| n + base[0]));
            if let Some(ends) = ends2 {
                let mate = usize::from(!INTERLEAVED);
                consumed[mate] = ends[3];
                blocks[mate].ends.push(ends.map(|n| n + base[mate]));
            }
            count += 1;
        }
        if count > 0 {
            blocks[0].data.extend_from_slice(&b1[..consumed[0]]);
            if let Some(b) = b2 {
                blocks[1].data.extend_from_slice(&b[..consumed[1]]);
            }
            first.inner.consume(consumed[0]);
            first.record += count * stride;
            if let Some(r) = second.as_mut() {
                r.inner.consume(consumed[1]);
                r.record += count;
            }
            if full {
                break;
            }
        } else if full {
            break;
        } else {
            let has_first = read_slow::<INTERLEAVED>(&mut blocks, first, second.as_deref_mut())
                .map_err(|error| {
                    let records = [first.record, second.as_ref().map_or(0, |r| r.record)];
                    earliest_error(&blocks, records, first.phred64).unwrap_or(error)
                })?;
            if !has_first {
                break;
            }
        }
    }
    Ok((blocks[0].len() > 0).then_some(blocks))
}

fn read_slow<const INTERLEAVED: bool>(
    blocks: &mut [Records; 2],
    first: &mut Reader,
    second: Option<&mut Reader>,
) -> Result<bool> {
    let has_first = blocks[0].read_one(first)?;
    if INTERLEAVED && has_first {
        ensure!(
            blocks[0].read_one(first)?,
            "Odd read count in interleaved FASTQ"
        );
    }
    if let Some(r) = second {
        let has_second = blocks[1].read_one(r)?;
        ensure!(
            has_first == has_second,
            "Paired input files have different read counts"
        );
    }
    Ok(has_first)
}

fn earliest_error(blocks: &[Records; 2], records: [usize; 2], phred64: bool) -> Option<String> {
    let first = records[0] - blocks[0].len();
    let second = records[1] - blocks[1].len();
    (0..blocks[0].len().max(blocks[1].len())).find_map(|i| {
        let check = |mate: usize, start: usize| {
            (i < blocks[mate].len())
                .then(|| blocks[mate].read(i, start + i, phred64).err())
                .flatten()
        };
        check(0, first).or_else(|| check(1, second))
    })
}
