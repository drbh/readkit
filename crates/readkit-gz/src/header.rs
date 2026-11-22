//! gzip member headers (RFC 1952), checked as zlib and flate2 check them:
//! magic, method 8, reserved flag bits clear, and the optional header CRC16.

use crate::error::Error;
use std::io::{self, BufRead};

const FHCRC: u8 = 2;
const FEXTRA: u8 = 4;
const FNAME: u8 = 8;
const FCOMMENT: u8 = 16;
const RESERVED: u8 = 0xe0;

/// Extra-field subfield that records a member's total compressed size, like
/// BGZF's `BC` but 32 bits so members can exceed 64 KiB. fasterp writes it
/// (and so does [`GzWriter`](crate::GzWriter)); other gzip readers skip it.
pub(crate) const MEMBER_SIZE_ID: [u8; 2] = *b"FP";

fn bad_header() -> Error {
    Error::corrupt("invalid gzip header")
}

fn check_fixed(fixed: &[u8]) -> Result<(), Error> {
    if fixed[..3] != [0x1f, 0x8b, 8] || fixed[3] & RESERVED != 0 {
        return Err(bad_header());
    }
    Ok(())
}

fn check_crc16(stored: [u8; 2], crc: u32) -> Result<(), Error> {
    if u16::from_le_bytes(stored) != crc as u16 {
        return Err(Error::checksum("gzip header checksum mismatch"));
    }
    Ok(())
}

/// Length in bytes of the member header at the start of `data`, or `None`
/// if `data` ends before the header does.
#[cfg_attr(not(feature = "parallel"), allow(dead_code))]
pub(crate) fn parse(data: &[u8]) -> Result<Option<usize>, Error> {
    if data.len() < 10 {
        return Ok(None);
    }
    check_fixed(data)?;
    let flags = data[3];
    let mut at = 10;
    if flags & FEXTRA != 0 {
        let Some(xlen) = data.get(at..at + 2) else {
            return Ok(None);
        };
        at += 2 + usize::from(u16::from_le_bytes([xlen[0], xlen[1]]));
    }
    for flag in [FNAME, FCOMMENT] {
        if flags & flag != 0 {
            let Some(rest) = data.get(at..) else {
                return Ok(None);
            };
            match memchr::memchr(0, rest) {
                Some(i) => at += i + 1,
                None => return Ok(None),
            }
        }
    }
    if flags & FHCRC != 0 {
        let Some(crc) = data.get(at..at + 2) else {
            return Ok(None);
        };
        check_crc16([crc[0], crc[1]], crc32fast::hash(&data[..at]))?;
        at += 2;
    }
    if at > data.len() {
        return Ok(None);
    }
    Ok(Some(at))
}

/// Read and check a member header from `input`; returns its length in bytes.
/// End of input inside the header is an `UnexpectedEof` I/O error.
pub(crate) fn read<R: BufRead>(input: &mut R) -> Result<u64, Error> {
    let mut fixed = [0u8; 10];
    input.read_exact(&mut fixed).map_err(Error::input)?;
    check_fixed(&fixed)?;
    let mut crc = crc32fast::Hasher::new();
    crc.update(&fixed);
    let flags = fixed[3];
    let mut len = 10u64;
    if flags & FEXTRA != 0 {
        let mut xlen = [0u8; 2];
        input.read_exact(&mut xlen).map_err(Error::input)?;
        crc.update(&xlen);
        let xlen = u64::from(u16::from_le_bytes(xlen));
        skip(input, &mut crc, Some(xlen)).map_err(Error::input)?;
        len += 2 + xlen;
    }
    for flag in [FNAME, FCOMMENT] {
        if flags & flag != 0 {
            len += skip(input, &mut crc, None).map_err(Error::input)?;
        }
    }
    if flags & FHCRC != 0 {
        let mut stored = [0u8; 2];
        input.read_exact(&mut stored).map_err(Error::input)?;
        check_crc16(stored, crc.finalize())?;
        len += 2;
    }
    Ok(len)
}

/// Consume `count` bytes, or through the next zero byte when `None`, adding
/// them to `crc`. Returns how many bytes were consumed.
fn skip<R: BufRead>(
    input: &mut R,
    crc: &mut crc32fast::Hasher,
    mut count: Option<u64>,
) -> io::Result<u64> {
    let mut total = 0u64;
    loop {
        if count == Some(0) {
            return Ok(total);
        }
        let buf = match input.fill_buf() {
            Ok(b) => b,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e),
        };
        if buf.is_empty() {
            return Err(io::ErrorKind::UnexpectedEof.into());
        }
        let (take, done) = match &mut count {
            Some(left) => {
                let take = buf.len().min(usize::try_from(*left).unwrap_or(usize::MAX));
                *left -= take as u64;
                (take, false)
            }
            None => match memchr::memchr(0, buf) {
                Some(i) => (i + 1, true),
                None => (buf.len(), false),
            },
        };
        crc.update(&buf[..take]);
        input.consume(take);
        total += take as u64;
        if done {
            return Ok(total);
        }
    }
}

/// Total compressed size of the member starting at `data`, from a size
/// subfield in its extra field (fasterp's `FP`, or BGZF's `BC`), or `None`
/// when it has none. Needs the member's first 12 bytes plus its extra field.
#[cfg_attr(not(feature = "parallel"), allow(dead_code))]
pub(crate) fn member_size(data: &[u8]) -> Option<usize> {
    if data.len() < 12 || data[..3] != [0x1f, 0x8b, 8] || data[3] & FEXTRA == 0 {
        return None;
    }
    let xlen = usize::from(u16::from_le_bytes([data[10], data[11]]));
    let mut extra = data.get(12..12 + xlen)?;
    while extra.len() >= 4 {
        let len = usize::from(u16::from_le_bytes([extra[2], extra[3]]));
        let field = extra.get(4..4 + len)?;
        match (&extra[..2], len) {
            (id, 4) if id == MEMBER_SIZE_ID => {
                return Some(u32::from_le_bytes(field.try_into().unwrap()) as usize);
            }
            (b"BC", 2) => return Some(usize::from(u16::from_le_bytes([field[0], field[1]])) + 1),
            _ => extra = &extra[4 + len..],
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slice_and_stream_parsers_agree() {
        let mut h = vec![
            0x1f,
            0x8b,
            8,
            FEXTRA | FNAME | FCOMMENT | FHCRC,
            0,
            0,
            0,
            0,
            0,
            3,
        ];
        h.extend_from_slice(&3u16.to_le_bytes());
        h.extend_from_slice(b"xyz");
        h.extend_from_slice(b"name\0comment\0");
        let crc = crc32fast::hash(&h) as u16;
        h.extend_from_slice(&crc.to_le_bytes());
        assert_eq!(parse(&h).unwrap(), Some(h.len()));
        for cut in 0..h.len() {
            assert_eq!(parse(&h[..cut]).unwrap(), None, "cut at {cut}");
        }
        let mut input = &h[..];
        assert_eq!(read(&mut input).unwrap(), h.len() as u64);
        let mut bad = h.clone();
        let n = bad.len();
        bad[n - 1] ^= 1;
        assert_eq!(parse(&bad).unwrap_err().kind(), crate::ErrorKind::Checksum);
        assert_eq!(
            read(&mut &bad[..]).unwrap_err().kind(),
            crate::ErrorKind::Checksum
        );
        bad = h.clone();
        bad[3] |= 0x20;
        assert_eq!(parse(&bad).unwrap_err().kind(), crate::ErrorKind::Corrupt);
        assert_eq!(
            read(&mut &h[..5]).unwrap_err().kind(),
            crate::ErrorKind::Truncated
        );
    }
}
