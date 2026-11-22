//! Decoding of VDB physical blobs for the fast path (feature `ngs`).
//!
//! A KDB blob of a version 2 VDB column is an envelope: one header byte
//! (version 2, a size-field variant, byte order, trailing pad bits), the
//! sizes of the encoding header and of the page map, then the encoding
//! header (`BlobHeaders`), the page map (row lengths, and which rows repeat
//! the previous row's data) and the payload. Formats follow ncbi-vdb:
//! `libs/vdb/blob.c` (envelope), `libs/vdb/blob-headers.c` (encoding
//! header), `libs/kdb/page-map.c` (page map), `libs/vxf/unzip.c` (zip),
//! `libs/klib/vlen-encode.c` (integers).
//!
//! Only what the fast path has been checked against is accepted. Anything
//! else is a [`FormatError`], and the caller decodes that blob through VDB.

use std::fmt;

/// A blob this module does not decode. Not a corruption verdict: VDB may
/// still decode it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct FormatError(pub &'static str);

impl fmt::Display for FormatError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.0)
    }
}

type Result<T> = std::result::Result<T, FormatError>;

fn err<T>(what: &'static str) -> Result<T> {
    Err(FormatError(what))
}

/// Big-endian base-128 unsigned integer (`vlen_decodeU1`, and the page
/// map's `deserialize_lengths`): bit 7 marks a continuation byte.
pub(crate) fn varint_u(src: &[u8], at: &mut usize) -> Result<u64> {
    let mut v: u64 = 0;
    for _ in 0..10 {
        let Some(&c) = src.get(*at) else {
            return err("truncated integer");
        };
        *at += 1;
        v = (v << 7) | u64::from(c & 0x7f);
        if c & 0x80 == 0 {
            return Ok(v);
        }
    }
    err("integer longer than 10 bytes")
}

/// Signed integer (`vlen_decode1`): the first byte holds a continuation
/// bit, a sign bit (0x40) and 6 value bits; later bytes 7 value bits each.
pub(crate) fn varint_i(src: &[u8], at: &mut usize) -> Result<i64> {
    let Some(&c) = src.get(*at) else {
        return err("truncated integer");
    };
    *at += 1;
    let negative = c & 0x40 != 0;
    let mut v = u64::from(c & 0x3f);
    let mut more = c & 0x80 != 0;
    let mut n = 1;
    while more {
        let Some(&c) = src.get(*at) else {
            return err("truncated integer");
        };
        *at += 1;
        n += 1;
        if n > 10 {
            return err("integer longer than 10 bytes");
        }
        v = (v << 7) | u64::from(c & 0x7f);
        more = c & 0x80 != 0;
    }
    let v = i64::try_from(v).map_err(|_| FormatError("integer out of range"))?;
    Ok(if negative { -v } else { v })
}

/// The parts of a version 2 blob.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Envelope<'a> {
    /// Pad bits at the end of the payload's last byte.
    pub adjust: u8,
    /// Serialized `BlobHeaders`; empty for a column without an encoding.
    pub header: &'a [u8],
    /// Serialized page map.
    pub page_map: &'a [u8],
    pub payload: &'a [u8],
}

/// Split a blob into its envelope (`decode_header_v2`).
pub(crate) fn envelope(blob: &[u8]) -> Result<Envelope<'_>> {
    let Some(&h) = blob.first() else {
        return err("empty blob");
    };
    if h >> 6 != 2 {
        return err("blob is not version 2");
    }
    let adjust = (8 - (h & 7)) & 7;
    let le = |b: &[u8]| {
        b.iter()
            .rev()
            .fold(0usize, |v, &x| (v << 8) | usize::from(x))
    };
    let field = |a: usize, b: usize| {
        blob.get(a..b)
            .map(le)
            .ok_or(FormatError("truncated envelope"))
    };
    let (hdr, map, at): (usize, usize, usize) = match (h >> 4) & 3 {
        0 => (field(1, 2)?, field(2, 3)?, 3),
        1 => (field(1, 2)?, field(2, 4)?, 4),
        2 => (field(1, 2)?, field(2, 6)?, 6),
        _ => (field(1, 5)?, field(5, 9)?, 9),
    };
    let end = at
        .checked_add(hdr)
        .and_then(|x| x.checked_add(map))
        .filter(|&e| e <= blob.len())
        .ok_or(FormatError("envelope sizes exceed the blob"))?;
    Ok(Envelope {
        adjust,
        header: &blob[at..at + hdr],
        page_map: &blob[at + hdr..end],
        payload: &blob[end..],
    })
}

/// A zip-encoded payload's header: zip format version and the size of the
/// data before compression, in bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ZipHeader {
    pub osize: u64,
}

/// Parse the encoding header of a `zip_encoding` column of bytes
/// (`deserialize_v0`; `zip_func` writes version 1 for byte-aligned input,
/// version 2 with one "trailing bits" argument otherwise). A blob with no
/// data at all has version 0 and size 0: zip never ran on it.
pub(crate) fn zip_header(h: &[u8]) -> Result<ZipHeader> {
    if h.first() != Some(&0) {
        return err("encoding header is not format 0");
    }
    let mut at = 1;
    let _flags = h.get(at).ok_or(FormatError("truncated encoding header"))?;
    let version = *h
        .get(at + 1)
        .ok_or(FormatError("truncated encoding header"))?;
    at += 2;
    let fmt = varint_i(h, &mut at)?;
    let osize = varint_i(h, &mut at)?;
    let ops = varint_i(h, &mut at)?;
    let args = varint_i(h, &mut at)?;
    if fmt != 0 || ops != 0 || osize < 0 {
        return err("unexpected encoding header fields");
    }
    match (version, args) {
        (0, 0) if osize == 0 => {}
        (1, 0) => {}
        (2, 1) => {
            // Bits beyond the last whole byte; bytes have none.
            if varint_i(h, &mut at)? != 0 {
                return err("zip of data that is not byte-aligned");
            }
        }
        _ => return err("unexpected zip version"),
    }
    if at != h.len() {
        // A further header means another encoding below zip.
        return err("chained encoding headers");
    }
    Ok(ZipHeader {
        osize: osize as u64,
    })
}

/// Where each row of a blob finds its data: element offset and length.
pub(crate) type Rows = Vec<(u64, u32)>;

/// Raw-deflate inflation of a page map body, given an upper bound on its size.
pub(crate) type Inflate<'a> = dyn FnMut(&[u8], usize) -> Result<Vec<u8>> + 'a;

/// Decode a page map for `row_count` rows (`PageMapDeserialize`). Random
/// access maps (format version 2) are not accepted.
pub(crate) fn page_map(src: &[u8], row_count: u32, inflate: &mut Inflate<'_>) -> Result<Rows> {
    let Some(&b0) = src.first() else {
        return err("empty page map");
    };
    match b0 >> 2 {
        0 => page_map_v0(src, row_count),
        1 => {
            // Header integers stay as they are; the rest is deflated.
            let variant = b0 & 3;
            let mut at = 1;
            let mut bound = 0usize;
            match variant {
                0 => return page_map_v0(src, row_count),
                1 => {
                    varint_u(src, &mut at)?;
                    bound = 5 * varint_u(src, &mut at)? as usize;
                }
                2 => bound = 10 * varint_u(src, &mut at)? as usize,
                _ => {
                    bound += 10 * varint_u(src, &mut at)? as usize;
                    bound += 5 * varint_u(src, &mut at)? as usize;
                }
            }
            let body = inflate(&src[at..], bound)?;
            let mut whole = Vec::with_capacity(at + body.len());
            whole.extend_from_slice(&src[..at]);
            whole.extend_from_slice(&body);
            page_map_v0(&whole, row_count)
        }
        2 => err("random-access page map"),
        _ => err("unknown page map version"),
    }
}

fn lengths(src: &[u8], at: &mut usize, n: u64) -> Result<Vec<u64>> {
    if n > src.len() as u64 {
        return err("page map lists more values than it has bytes");
    }
    (0..n).map(|_| varint_u(src, at)).collect()
}

fn page_map_v0(src: &[u8], row_count: u32) -> Result<Rows> {
    let variant = src[0] & 3;
    let rows = u64::from(row_count);
    let mut at = 1;
    // (length, rows with that length) runs, and data records: how many
    // consecutive rows share one copy of the data.
    let (length_runs, data_runs): (Vec<(u64, u64)>, Option<Vec<u64>>) = match variant {
        0 => (vec![(varint_u(src, &mut at)?, rows)], None),
        1 => {
            let len = varint_u(src, &mut at)?;
            let recs = varint_u(src, &mut at)?;
            (vec![(len, rows)], Some(lengths(src, &mut at, recs)?))
        }
        2 => {
            let recs = varint_u(src, &mut at)?;
            let v = lengths(src, &mut at, 2 * recs)?;
            let (len, run) = v.split_at(recs as usize);
            (len.iter().copied().zip(run.iter().copied()).collect(), None)
        }
        _ => {
            let lrecs = varint_u(src, &mut at)?;
            let drecs = varint_u(src, &mut at)?;
            let v = lengths(src, &mut at, 2 * lrecs + drecs)?;
            let (len, rest) = v.split_at(lrecs as usize);
            let (run, data) = rest.split_at(lrecs as usize);
            (
                len.iter().copied().zip(run.iter().copied()).collect(),
                Some(data.to_vec()),
            )
        }
    };
    if length_runs.iter().map(|r| r.1).sum::<u64>() != rows {
        return err("page map length runs do not cover the blob's rows");
    }
    if let Some(d) = &data_runs
        && (d.iter().sum::<u64>() != rows || d.contains(&0))
    {
        return err("page map data runs do not cover the blob's rows");
    }
    let mut out = Vec::with_capacity(row_count as usize);
    for &(len, run) in &length_runs {
        let len = u32::try_from(len).map_err(|_| FormatError("row longer than 4 G elements"))?;
        out.extend(std::iter::repeat_n((0u64, len), run as usize));
    }
    // Offsets: each data record starts where the previous record's data
    // ended; rows of one record share its offset (and must share its length).
    let mut offset = 0u64;
    match data_runs {
        None => {
            for row in &mut out {
                row.0 = offset;
                offset += u64::from(row.1);
            }
        }
        Some(runs) => {
            let mut r = 0;
            for run in runs {
                let len = out[r].1;
                for row in &mut out[r..r + run as usize] {
                    if row.1 != len {
                        return err("rows sharing data differ in length");
                    }
                    row.0 = offset;
                }
                offset += u64::from(len);
                r += run as usize;
            }
        }
    }
    Ok(out)
}

/// 4na codes as text (`out_dna_text`).
pub(crate) const TEXT_4NA: &[u8; 16] = b".ACMGRSVTWYHKDBN";
/// 2na code to 4na (`out_2na_4na_bin`).
pub(crate) const FOUR_OF_TWO: [u8; 4] = [1, 2, 4, 8];

/// Write `len` bases starting at base `first` of a 2na packed buffer
/// (four bases per byte, first base in the high bits) as text.
pub(crate) fn unpack_2na(packed: &[u8], first: u64, len: usize, out: &mut Vec<u8>) {
    static LUT: [[u8; 4]; 256] = {
        let mut t = [[0u8; 4]; 256];
        let mut i = 0;
        while i < 256 {
            let mut j = 0;
            while j < 4 {
                t[i][j] = b"ACGT"[(i >> (6 - 2 * j)) & 3];
                j += 1;
            }
            i += 1;
        }
        t
    };
    out.clear();
    out.reserve(len);
    let mut b = first as usize;
    let end = b + len;
    while b < end && !b.is_multiple_of(4) {
        out.push(b"ACGT"[(packed[b / 4] >> (6 - 2 * (b % 4))) as usize & 3]);
        b += 1;
    }
    while b + 4 <= end {
        out.extend_from_slice(&LUT[packed[b / 4] as usize]);
        b += 4;
    }
    while b < end {
        out.push(b"ACGT"[(packed[b / 4] >> (6 - 2 * (b % 4))) as usize & 3]);
        b += 1;
    }
}

/// Apply a row of `ALTREAD` to text bases (`bit_or #1 <1>` of the 2na-as-4na
/// bases with `.ALTREAD`, whose rows are trimmed of leading zeros, so they
/// align with the end of the read), then text. `bases` holds the read's
/// 2na bases as `ACGT` text.
pub(crate) fn apply_alt(bases: &mut [u8], alt: &[u8]) -> Result<()> {
    if alt.len() > bases.len() {
        return err("ALTREAD row longer than its read");
    }
    let start = bases.len() - alt.len();
    for (b, &a) in bases[start..].iter_mut().zip(alt) {
        if a == 0 {
            continue;
        }
        if a > 15 {
            return err("ALTREAD value above 15");
        }
        let two = match *b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            _ => 3,
        };
        *b = TEXT_4NA[usize::from(FOUR_OF_TWO[two] | a)];
    }
    Ok(())
}

/// Element count of a 2na packed payload: 2 bits per base, `adjust` pad
/// bits in the last byte.
pub(crate) fn bases_2na(payload: &[u8], adjust: u8) -> Result<u64> {
    let bits = (payload.len() as u64 * 8)
        .checked_sub(u64::from(adjust))
        .ok_or(FormatError("pad bits exceed the payload"))?;
    if bits % 2 != 0 {
        return err("2na payload with an odd bit count");
    }
    Ok(bits / 2)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn enc_u(v: u64, out: &mut Vec<u8>) {
        let mut groups = vec![(v & 0x7f) as u8];
        let mut v = v >> 7;
        while v > 0 {
            groups.push((v & 0x7f) as u8 | 0x80);
            v >>= 7;
        }
        out.extend(groups.iter().rev());
    }

    fn no_inflate(_: &[u8], _: usize) -> Result<Vec<u8>> {
        panic!("no inflate expected")
    }

    #[test]
    fn unsigned_integers() {
        for v in [0u64, 1, 127, 128, 131_156, u32::MAX as u64, 1 << 62] {
            let mut b = Vec::new();
            enc_u(v, &mut b);
            let mut at = 0;
            assert_eq!(varint_u(&b, &mut at).unwrap(), v);
            assert_eq!(at, b.len());
        }
        assert!(varint_u(&[0x80], &mut 0).is_err());
        assert!(varint_u(&[0xff; 11], &mut 0).is_err());
    }

    #[test]
    fn signed_integers() {
        // From a real QUALITY header: osize 131 156 is 88 80 54.
        let mut at = 0;
        assert_eq!(varint_i(&[0x88, 0x80, 0x54], &mut at).unwrap(), 131_156);
        assert_eq!(at, 3);
        assert_eq!(varint_i(&[0x00], &mut 0).unwrap(), 0);
        assert_eq!(varint_i(&[0x41], &mut 0).unwrap(), -1);
        assert_eq!(varint_i(&[0x3f], &mut 0).unwrap(), 63);
        assert_eq!(varint_i(&[0x80, 0x40], &mut 0).unwrap(), 64);
        assert!(varint_i(&[0x80], &mut 0).is_err());
    }

    #[test]
    fn envelopes() {
        // Variant 1 (as in the NovaSeq QUALITY blobs): 9-byte header, map of 2.
        let mut b = vec![0x90, 9, 2, 0];
        b.extend([0xaa; 9]);
        b.extend([0xbb; 2]);
        b.extend([1, 2, 3]);
        let e = envelope(&b).unwrap();
        assert_eq!(
            (e.adjust, e.header.len(), e.page_map, e.payload),
            (0, 9, &[0xbb, 0xbb][..], &[1, 2, 3][..])
        );
        // Variant 0 with 2 pad bits (low 3 bits = 6).
        let e = envelope(&[0x86, 0, 1, 0x00, 7]).unwrap();
        assert_eq!((e.adjust, e.page_map, e.payload), (2, &[0][..], &[7][..]));
        // Variant 3: 32-bit sizes.
        let mut b = vec![0xb0, 1, 0, 0, 0, 0, 0, 0, 0, 0x55, 9];
        let e = envelope(&b).unwrap();
        assert_eq!((e.header, e.payload), (&[0x55][..], &[9][..]));
        b[1] = 9;
        assert!(envelope(&b).is_err(), "sizes beyond the blob");
        assert!(envelope(&[0x50, 0, 0]).is_err(), "version 1");
        assert!(envelope(&[]).is_err());
    }

    #[test]
    fn zip_headers() {
        // Real NovaSeq QUALITY blob header.
        let h = [0x00, 0x00, 0x01, 0x00, 0x88, 0x80, 0x54, 0x00, 0x00];
        assert_eq!(zip_header(&h).unwrap(), ZipHeader { osize: 131_156 });
        // An empty ALTREAD blob, zipped (NovaSeq run) and not (human run:
        // flags 1, version 0).
        assert_eq!(zip_header(&[0, 0, 1, 0, 0, 0, 0]).unwrap().osize, 0);
        assert_eq!(zip_header(&[0, 1, 0, 0, 0, 0, 0]).unwrap().osize, 0);
        assert!(
            zip_header(&[0, 1, 0, 0, 0x10, 0, 0]).is_err(),
            "version 0 with data"
        );
        // Version 2 with 0 trailing bits is fine; with 4 it is not bytes.
        assert!(zip_header(&[0, 0, 2, 0, 0x10, 0, 1, 0]).is_ok());
        assert!(zip_header(&[0, 0, 2, 0, 0x10, 0, 1, 4]).is_err());
        assert!(
            zip_header(&[0, 0, 3, 0, 0x10, 0, 0]).is_err(),
            "unknown version"
        );
        assert!(zip_header(&[0, 0, 1, 0, 0x10, 1, 0, 7]).is_err(), "an op");
        assert!(
            zip_header(&[0, 0, 1, 0, 0x10, 0, 0, 0]).is_err(),
            "a chained header"
        );
        assert!(
            zip_header(&[1, 0, 1, 0, 0x10, 0, 0]).is_err(),
            "header format 1"
        );
    }

    #[test]
    fn page_map_fixed_length() {
        // v0 variant 0: every row 151 long, each with its own data.
        let mut m = vec![0x00];
        enc_u(151, &mut m);
        let rows = page_map(&m, 3, &mut no_inflate).unwrap();
        assert_eq!(rows, vec![(0, 151), (151, 151), (302, 151)]);
    }

    #[test]
    fn page_map_fixed_length_shared_data() {
        // v0 variant 1: length 10, data runs 2,1 over 3 rows.
        let mut m = vec![0x01];
        enc_u(10, &mut m);
        enc_u(2, &mut m);
        enc_u(2, &mut m);
        enc_u(1, &mut m);
        let rows = page_map(&m, 3, &mut no_inflate).unwrap();
        assert_eq!(rows, vec![(0, 10), (0, 10), (10, 10)]);
    }

    #[test]
    fn page_map_variable_lengths() {
        // v0 variant 2: lengths 5 (2 rows), 300 (1 row), 0 (1 row).
        let mut m = vec![0x02];
        enc_u(3, &mut m);
        for v in [5, 300, 0, 2, 1, 1] {
            enc_u(v, &mut m);
        }
        let rows = page_map(&m, 4, &mut no_inflate).unwrap();
        assert_eq!(rows, vec![(0, 5), (5, 5), (10, 300), (310, 0)]);
    }

    #[test]
    fn page_map_variable_lengths_shared_data() {
        // v0 variant 3: lengths 4 (3 rows), 7 (1 row); data runs 1, 2, 1:
        // rows 1 and 2 share data.
        let mut m = vec![0x03];
        enc_u(2, &mut m);
        enc_u(3, &mut m);
        for v in [4, 7, 3, 1, 1, 2, 1] {
            enc_u(v, &mut m);
        }
        let rows = page_map(&m, 4, &mut no_inflate).unwrap();
        assert_eq!(rows, vec![(0, 4), (4, 4), (4, 4), (8, 7)]);
    }

    #[test]
    fn page_map_deflated_body() {
        // v1 variant 2: the integers after the header count are deflated.
        let mut body = Vec::new();
        for v in [8, 2] {
            enc_u(v, &mut body);
        }
        let m = [0x06, 1, 0xde, 0xad];
        let mut inflate = |src: &[u8], bound: usize| {
            assert_eq!((src, bound), (&[0xde, 0xad][..], 10));
            Ok(body.clone())
        };
        assert_eq!(page_map(&m, 2, &mut inflate).unwrap(), vec![(0, 8), (8, 8)]);
    }

    #[test]
    fn page_map_rejects_what_it_cannot_check() {
        // Runs that do not cover the rows.
        let mut m = vec![0x02];
        enc_u(1, &mut m);
        enc_u(5, &mut m);
        enc_u(2, &mut m);
        assert!(page_map(&m, 3, &mut no_inflate).is_err());
        // Shared data with differing lengths.
        let mut m = vec![0x03];
        for v in [2, 1, 4, 7, 1, 1, 2] {
            enc_u(v, &mut m);
        }
        assert!(page_map(&m, 2, &mut no_inflate).is_err());
        // Random access (version 2).
        assert!(page_map(&[0x08, 1], 1, &mut no_inflate).is_err());
        // More values than bytes.
        assert!(page_map(&[0x02, 0x7f], 1, &mut no_inflate).is_err());
        assert!(page_map(&[], 1, &mut no_inflate).is_err());
    }

    #[test]
    fn unpacking_2na() {
        // 0b00_01_10_11 = ACGT, 0b11_10_01_00 = TGCA
        let packed = [0x1b, 0xe4, 0x1b];
        let mut out = Vec::new();
        for first in 0..12u64 {
            for len in 0..(12 - first as usize) {
                unpack_2na(&packed, first, len, &mut out);
                let all = b"ACGTTGCAACGT";
                assert_eq!(
                    out,
                    &all[first as usize..first as usize + len],
                    "{first}+{len}"
                );
            }
        }
        assert_eq!(bases_2na(&packed, 0).unwrap(), 12);
        assert_eq!(bases_2na(&packed, 2).unwrap(), 11);
        assert!(bases_2na(&packed, 1).is_err());
        assert!(bases_2na(&[], 2).is_err());
    }

    #[test]
    fn ambiguity_codes() {
        // N over A (1 | 15 = 15), M over C (2 | 3 = 3), aligned to the end.
        let mut b = b"ACGTAC".to_vec();
        apply_alt(&mut b, &[0, 0, 15, 3]).unwrap();
        assert_eq!(b, b"ACGTNM");
        let mut b = b"AC".to_vec();
        apply_alt(&mut b, &[]).unwrap();
        assert_eq!(b, b"AC");
        assert!(apply_alt(&mut b"A".to_vec(), &[1, 2]).is_err());
        assert!(apply_alt(&mut b"A".to_vec(), &[16]).is_err());
    }
}
