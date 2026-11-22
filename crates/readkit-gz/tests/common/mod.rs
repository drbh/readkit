//! Inputs and reference decoders shared by the integration tests.
#![allow(dead_code)]

use readkit_gz::{DecodeOptions, Error, GzReader, ParGzReader, Split, Strategy};
use std::io::{self, Read, Write};

/// FASTQ text with variable read lengths, occasional N and binned qualities.
pub fn fastq(n: usize, seed: u64) -> Vec<u8> {
    let mut state = seed;
    let mut next = |m: u64| {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) % m
    };
    let mut out = Vec::new();
    for i in 0..n {
        let len = 50 + next(100) as usize;
        writeln!(out, "@SRR5808766.{i} {i} length={len}").unwrap();
        out.extend((0..len).map(|_| b"ACGTN"[next(if i % 50 == 0 { 5 } else { 4 }) as usize]));
        out.extend_from_slice(b"\n+\n");
        out.extend((0..len).map(|_| b"#,:FFFF"[next(7) as usize]));
        out.push(b'\n');
    }
    out
}

pub fn flate2_gz(data: &[u8], level: u32) -> Vec<u8> {
    let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(level));
    e.write_all(data).unwrap();
    e.finish().unwrap()
}

/// A flate2 member whose header carries a file name, as `gzip file` writes.
pub fn named_gz(data: &[u8], name: &str) -> Vec<u8> {
    let mut e = flate2::GzBuilder::new()
        .filename(name)
        .write(Vec::new(), flate2::Compression::new(6));
    e.write_all(data).unwrap();
    e.finish().unwrap()
}

pub fn libdeflate_gz(data: &[u8], level: i32) -> Vec<u8> {
    let mut c = libdeflater::Compressor::new(libdeflater::CompressionLvl::new(level).unwrap());
    let mut out = vec![0; c.gzip_compress_bound(data.len())];
    let n = c.gzip_compress(data, &mut out).unwrap();
    out.truncate(n);
    out
}

pub fn fq_gz(data: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    readkit_gz::fastq_gzip_member(data, &mut out);
    out
}

/// One member from each kind of compressor, so decoders meet stored, fixed
/// and dynamic blocks.
pub fn every_compressor(data: &[u8]) -> Vec<(&'static str, Vec<u8>)> {
    vec![
        ("flate2-0", flate2_gz(data, 0)),
        ("flate2-1", flate2_gz(data, 1)),
        ("flate2-6", flate2_gz(data, 6)),
        ("flate2-9", flate2_gz(data, 9)),
        ("libdeflate-1", libdeflate_gz(data, 1)),
        ("libdeflate-12", libdeflate_gz(data, 12)),
        ("fqdeflate", fq_gz(data)),
    ]
}

/// The BGZF end-of-file marker block.
pub const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0, 0x1b, 0, 3, 0, 0, 0, 0, 0, 0, 0,
    0, 0,
];

/// A BGZF member: `BC` subfield holding the member size minus one.
pub fn bgzf_member(data: &[u8]) -> Vec<u8> {
    let plain = libdeflate_gz(data, 6);
    let total = plain.len() + 8;
    assert!(total <= 1 << 16);
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

/// `data` as BGZF: 64 KiB blocks and the EOF marker.
pub fn bgzf(data: &[u8]) -> Vec<u8> {
    let mut out: Vec<u8> = data.chunks(0xff00).flat_map(bgzf_member).collect();
    out.extend_from_slice(&BGZF_EOF);
    out
}

pub fn flate2_decode(gz: &[u8]) -> io::Result<Vec<u8>> {
    let mut out = Vec::new();
    flate2::read::MultiGzDecoder::new(gz).read_to_end(&mut out)?;
    Ok(out)
}

fn typed(e: io::Error) -> Error {
    Error::from_io(&e)
        .cloned()
        .unwrap_or_else(|| panic!("untyped error: {e}"))
}

pub fn sequential(gz: &[u8]) -> Result<Vec<u8>, Error> {
    let mut out = Vec::new();
    GzReader::from_reader(gz)
        .read_to_end(&mut out)
        .map_err(typed)?;
    Ok(out)
}

pub fn parallel_with(gz: &[u8], options: &DecodeOptions) -> Result<(Vec<u8>, Strategy), Error> {
    let mut reader = ParGzReader::new(gz.to_vec(), gz.len() as u64, options)?;
    let strategy = reader.strategy();
    let mut out = Vec::new();
    reader.read_to_end(&mut out).map_err(typed)?;
    Ok((out, strategy))
}

/// Region decoding of every member, in regions of `region` bytes.
pub fn regions(gz: &[u8], region: usize, threads: usize) -> Result<Vec<u8>, Error> {
    let options = DecodeOptions {
        threads: Some(threads),
        split: Split::Always,
        region_size: region,
        ..Default::default()
    };
    parallel_with(gz, &options).map(|(out, _)| out)
}

/// Every decoding path of this crate: the sequential reader, the
/// background sequential reader, sized-member or region decoding as chosen,
/// and forced region decoding at several region sizes and thread counts.
/// All must agree; returns their common result.
pub fn all_paths(gz: &[u8]) -> Result<Vec<u8>, Error> {
    let expected = sequential(gz);
    let mut results = Vec::new();
    let never = DecodeOptions {
        split: Split::Never,
        threads: Some(2),
        ..Default::default()
    };
    results.push(("never", parallel_with(gz, &never).map(|r| r.0)));
    let auto = DecodeOptions {
        threads: Some(3),
        region_size: 64 << 10,
        ..Default::default()
    };
    results.push(("auto", parallel_with(gz, &auto).map(|r| r.0)));
    // Region sizes relative to the input, so boundaries land in different
    // places without decoding thousands of regions in a debug build.
    let n = gz.len();
    for (region, threads) in [(64 + n / 24, 1), (100 + n / 7, 2), (n / 2 + 1, 4)] {
        results.push(("regions", regions(gz, region, threads)));
    }
    let mut piped = Vec::new();
    let r = ParGzReader::from_reader(Trickle::new(gz.to_vec(), 7), &DecodeOptions::default())
        .unwrap()
        .read_to_end(&mut piped)
        .map(|_| piped)
        .map_err(typed);
    results.push(("from_reader", r));
    for (name, result) in results {
        match (&expected, &result) {
            (Ok(a), Ok(b)) => assert!(a == b, "{name}: different bytes"),
            (Err(a), Err(b)) => assert_eq!(a.kind(), b.kind(), "{name}: {a} vs {b}"),
            (a, b) => panic!("{name}: sequential {a:?} but {b:?}"),
        }
    }
    expected
}

/// A pipe-like reader: no seeking, and reads of varying small sizes.
pub struct Trickle {
    data: io::Cursor<Vec<u8>>,
    max: usize,
    turn: usize,
}

impl Trickle {
    pub fn new(data: Vec<u8>, max: usize) -> Self {
        Self {
            data: io::Cursor::new(data),
            max,
            turn: 0,
        }
    }
}

impl Read for Trickle {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.turn += 1;
        if self.turn.is_multiple_of(5) {
            return Err(io::ErrorKind::Interrupted.into());
        }
        let n = buf.len().min(1 + self.turn * 31 % self.max);
        self.data.read(&mut buf[..n])
    }
}
