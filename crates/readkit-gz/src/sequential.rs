//! Sequential multi-member gzip decoding of any `BufRead`.

use crate::error::Error;
use crate::header;
use flate2::{Decompress, FlushDecompress, Status};
use std::io::{self, BufRead, BufReader, Read};

/// Bytes of compressed input buffered by [`GzReader::from_reader`].
pub(crate) const INPUT_BUFFER: usize = 1 << 20;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum State {
    /// At a member boundary: a header follows, or the end of input.
    Header,
    /// Inside a member's deflate stream.
    Body,
    Done,
}

/// Sequential gzip decoder for any input, including streams that cannot
/// seek (pipes, sockets, a download in progress).
///
/// Decodes every member of a multi-member stream (concatenated `.gz` files,
/// BGZF) and checks each member's header, CRC32 and ISIZE. It accepts exactly
/// what flate2's `MultiGzDecoder` accepts: an empty input is
/// [`Truncated`](crate::ErrorKind::Truncated), and bytes after a member must
/// start another member ([`Corrupt`](crate::ErrorKind::Corrupt) otherwise).
///
/// Output is produced on the calling thread, so dropping the reader stops
/// decoding at once. Reaching end of file (`Ok(0)`) means every member was
/// verified. After an error every later read returns the same error.
///
/// The inner reader must block: `WouldBlock` is treated as a failure.
///
/// ```
/// use std::io::{BufReader, Read};
/// # let gz: Vec<u8> = {
/// #     use std::io::Write;
/// #     let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
/// #     e.write_all(b"@r\nACGT\n+\nFFFF\n").unwrap();
/// #     e.finish().unwrap()
/// # };
/// let mut text = String::new();
/// readkit_gz::GzReader::from_reader(&gz[..]).read_to_string(&mut text)?;
/// assert_eq!(text, "@r\nACGT\n+\nFFFF\n");
/// # Ok::<(), std::io::Error>(())
/// ```
pub struct GzReader<R> {
    inner: R,
    state: State,
    inflate: Decompress,
    crc: crc32fast::Hasher,
    /// Decoded bytes of the current member.
    produced: u64,
    /// Compressed bytes consumed so far.
    consumed: u64,
    /// Where the current member starts.
    member_start: u64,
    /// Offset of the input's first byte in a larger file, for errors.
    base: u64,
    members: u64,
    error: Option<Error>,
}

impl<R: Read> GzReader<BufReader<R>> {
    /// Decode `reader`, buffering 1 MiB of compressed input.
    pub fn from_reader(reader: R) -> Self {
        GzReader::new(BufReader::with_capacity(INPUT_BUFFER, reader))
    }
}

impl<R: BufRead> GzReader<R> {
    /// Decode `inner`, reading through its own buffer. The decoder consumes
    /// exactly the gzip bytes it uses, so `inner` can be recovered with
    /// [`into_inner`](Self::into_inner) after end of file.
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            state: State::Header,
            inflate: Decompress::new(false),
            crc: crc32fast::Hasher::new(),
            produced: 0,
            consumed: 0,
            member_start: 0,
            base: 0,
            members: 0,
            error: None,
        }
    }

    /// Report error offsets relative to `base`, the input's position in a
    /// larger file.
    #[cfg_attr(not(feature = "parallel"), allow(dead_code))]
    pub(crate) fn with_base(mut self, base: u64) -> Self {
        self.base = base;
        self
    }

    /// Members decoded and verified so far.
    pub fn members(&self) -> u64 {
        self.members
    }

    /// Compressed bytes consumed so far.
    pub fn compressed_bytes(&self) -> u64 {
        self.consumed
    }

    /// The inner reader.
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// The inner reader, positioned after the last byte consumed.
    pub fn into_inner(self) -> R {
        self.inner
    }

    /// Like `read`, but with the typed error.
    pub(crate) fn read_typed(&mut self, out: &mut [u8]) -> Result<usize, Error> {
        if let Some(e) = &self.error {
            return Err(e.clone());
        }
        if out.is_empty() {
            return Ok(0);
        }
        match self.step(out) {
            Ok(n) => Ok(n),
            Err(e) => {
                let e = e.at(self.base + self.member_start);
                self.error = Some(e.clone());
                Err(e)
            }
        }
    }

    fn fill(inner: &mut R) -> Result<&[u8], Error> {
        let available = loop {
            match inner.fill_buf() {
                Ok(buf) => break buf.len(),
                Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                Err(e) => return Err(Error::input(e)),
            }
        };
        if available == 0 {
            return Ok(&[]);
        }
        // Buffered bytes are returned again without another read; the loop
        // above cannot return the borrow itself.
        inner.fill_buf().map_err(Error::input)
    }

    fn step(&mut self, out: &mut [u8]) -> Result<usize, Error> {
        let mut stalls = 0;
        loop {
            match self.state {
                State::Done => return Ok(0),
                State::Header => {
                    if Self::fill(&mut self.inner)?.is_empty() {
                        if self.members == 0 {
                            return Err(Error::truncated());
                        }
                        self.state = State::Done;
                        return Ok(0);
                    }
                    self.member_start = self.consumed;
                    self.consumed += header::read(&mut self.inner)?;
                    self.inflate.reset(false);
                    self.crc.reset();
                    self.produced = 0;
                    self.state = State::Body;
                }
                State::Body => {
                    let input = Self::fill(&mut self.inner)?;
                    let eof = input.is_empty();
                    let (in0, out0) = (self.inflate.total_in(), self.inflate.total_out());
                    let status = self
                        .inflate
                        .decompress(input, out, FlushDecompress::None)
                        .map_err(|e| Error::corrupt("invalid deflate data").with_source(e))?;
                    let used = (self.inflate.total_in() - in0) as usize;
                    let made = (self.inflate.total_out() - out0) as usize;
                    self.inner.consume(used);
                    self.consumed += used as u64;
                    self.crc.update(&out[..made]);
                    self.produced += made as u64;
                    if status == Status::StreamEnd {
                        self.trailer()?;
                        if made > 0 {
                            return Ok(made);
                        }
                        continue;
                    }
                    if made > 0 {
                        return Ok(made);
                    }
                    if eof {
                        return Err(Error::truncated());
                    }
                    if used == 0 {
                        stalls += 1;
                        if stalls > 2 {
                            return Err(Error::corrupt("deflate decoder made no progress"));
                        }
                    }
                }
            }
        }
    }

    fn trailer(&mut self) -> Result<(), Error> {
        let mut trailer = [0u8; 8];
        self.inner.read_exact(&mut trailer).map_err(Error::input)?;
        self.consumed += 8;
        let crc = u32::from_le_bytes(trailer[..4].try_into().unwrap());
        let isize = u32::from_le_bytes(trailer[4..].try_into().unwrap());
        let actual = std::mem::take(&mut self.crc).finalize();
        if crc != actual || isize != self.produced as u32 {
            return Err(Error::checksum(
                "gzip member CRC32 or length does not match its trailer",
            ));
        }
        self.members += 1;
        self.state = State::Header;
        Ok(())
    }
}

impl<R: BufRead> Read for GzReader<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        self.read_typed(out).map_err(io::Error::from)
    }
}

impl<R> std::fmt::Debug for GzReader<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GzReader")
            .field("state", &self.state)
            .field("members", &self.members)
            .field("compressed_bytes", &self.consumed)
            .finish_non_exhaustive()
    }
}
