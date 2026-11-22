//! gzip for sequencing reads: decoding behind [`std::io::Read`], in parallel
//! where the input allows it, and a FASTQ-aware parallel encoder behind
//! [`std::io::Write`].
//!
//! Extracted from fasterp's `pgzip`, `io` and `fqdeflate` modules.
//!
//! # Decoding
//!
//! - [`GzReader`] decodes any reader sequentially: pipes, sockets, a
//!   download in progress. Every member of a multi-member stream is decoded,
//!   with header, CRC32 and ISIZE checks. No thread is started. This is the
//!   only decoder on wasm32.
//! - [`ParGzReader`] (feature `parallel`) decodes positioned input (a local
//!   file, bytes in memory, anything implementing [`ReadAt`]) on background
//!   threads, choosing a [`Strategy`] from the input:
//!   - members that record their size (BGZF, fasterp and [`GzWriter`]
//!     output) are decoded member-parallel with libdeflate;
//!   - a large single member is cut into regions decoded in parallel
//!     (pugz / rapidgzip style; see [`Split`]), and so are later members
//!     while they stay large (concatenated files);
//!   - anything else is decoded sequentially on one background thread.
//!
//!   Output and errors are identical to [`GzReader`]'s for every strategy
//!   and thread count. [`open`] opens a path.
//!
//! ```no_run
//! use std::io::BufRead;
//! let options = readkit_gz::DecodeOptions { threads: Some(4), ..Default::default() };
//! let reader = readkit_gz::open("reads.fastq.gz", &options)?;
//! let records = reader.lines().count() / 4;
//! # Ok::<(), readkit_gz::Error>(())
//! ```
//!
//! # Encoding
//!
//! [`GzWriter`] (feature `parallel`) compresses blocks into independent
//! gzip members on worker threads, with libdeflate or, for binned FASTQ
//! qualities at levels 2 to 4, the FASTQ-aware encoder ([`Fastq`]).
//! [`fastq_gzip_member`] is that encoder for one block, on every target.
//!
//! # Errors
//!
//! Fallible functions return [`Error`]; its [`ErrorKind`] separates corrupt
//! data, checksum mismatches, truncated input and I/O failures. Readers and
//! writers return it inside an [`std::io::Error`]; [`Error::from_io`]
//! recovers it.
//!
//! # Features
//!
//! - `parallel` (default): threads and libdeflate (C). Without it the crate
//!   is pure Rust (flate2 with the zlib-rs backend) and builds for
//!   `wasm32-unknown-unknown`.

#![deny(missing_docs)]

mod error;
mod fqdeflate;
mod header;
mod sequential;

#[cfg(feature = "parallel")]
mod encode;
#[cfg(feature = "parallel")]
mod par;
#[cfg(feature = "parallel")]
mod pgzip;
#[cfg(feature = "parallel")]
mod read_at;

pub use error::{Error, ErrorKind};
pub use sequential::GzReader;

#[cfg(feature = "parallel")]
pub use encode::{EncodeOptions, Fastq, GzWriter};
#[cfg(feature = "parallel")]
pub use par::{DecodeOptions, ParGzReader, Split, Strategy, open};
#[cfg(feature = "parallel")]
pub use read_at::ReadAt;

use std::io::Read;

/// Whether a block of FASTQ text suits the FASTQ-aware encoder: its first
/// records have binned qualities (at most six distinct values, as NovaSeq
/// and NextSeq 2000 write them) and do not overlap like reads sorted by
/// position.
///
/// Decoding a large member in parallel costs extra CPU; fasterp uses this
/// on the first 256 KiB of input to decide whether gzip output will be cheap
/// enough to leave that CPU to decoding ([`binned_gzip`]).
pub fn binned_qualities(fastq: &[u8]) -> bool {
    fqdeflate::suits(fastq)
}

/// [`binned_qualities`] of the first 256 KiB of the gzip stream `gz`
/// (fasterp's `binned_gzip_input`).
pub fn binned_gzip<R: Read>(gz: R) -> Result<bool, Error> {
    let mut head = Vec::with_capacity(256 << 10);
    let mut reader = GzReader::new(std::io::BufReader::new(gz));
    let mut buf = vec![0; 64 << 10];
    while head.len() < 256 << 10 {
        match reader.read_typed(&mut buf)? {
            0 => break,
            n => head.extend_from_slice(&buf[..n]),
        }
    }
    head.truncate(256 << 10);
    Ok(binned_qualities(&head))
}

/// Append one complete gzip member for `input` to `out`, written by the
/// FASTQ-aware deflate encoder.
///
/// The encoder exploits FASTQ structure instead of a general-purpose match
/// finder: sequence lines become Huffman literals (plus one probe per read
/// for duplicates), headers are matched against the previous header, and
/// `+` lines and qualities go through a small hash matcher. Line
/// classification only affects size: any input gives valid gzip that every
/// inflater decodes back to `input`. The member carries an `FP` extra
/// subfield with its compressed size, which [`ParGzReader`] uses to decode
/// members in parallel.
///
/// One dynamic Huffman block per member, so blocks of 64 KiB to 1 MiB suit
/// it best ([`GzWriter`] uses 256 KiB).
///
/// # Panics
///
/// If `input` is larger than 1 GiB.
pub fn fastq_gzip_member(input: &[u8], out: &mut Vec<u8>) {
    assert!(
        input.len() <= 1 << 30,
        "fastq_gzip_member input above 1 GiB"
    );
    fqdeflate::gzip_member(input, out);
}
