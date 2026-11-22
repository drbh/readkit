//! Ranged, verified, resumable HTTP downloads.
//!
//! A [`Fetcher`] downloads an [`Object`] (URL, size, [`Checksum`]) into a
//! [`Sink`] over parallel range requests:
//!
//! - per-host connection budgets that adapt (AIMD) and are shared by every
//!   download to that host ([`FetchOptions::host_limits`]);
//! - per-chunk retries with exponential delay at the back of the queue;
//!   permanent 4xx errors fail at once;
//! - hedging: the oldest slow chunk is re-issued once nothing else is queued;
//! - verification by MD5 or S3 multipart ETag, and BLAKE3-based repair that
//!   rewrites only corrupt chunks;
//! - resumable file downloads (a ledger next to the `.partial` file).
//!
//! Sinks: [`FileSink`] (resumable, renamed into place when verified),
//! [`MemorySink`] (`Vec<u8>`) and [`StreamSink`] / [`Fetcher::stream`] (an
//! in-order `impl Read` with bounded read-ahead, usable while the download is
//! still running).
//!
//! ```no_run
//! use readkit_fetch::{Checksum, FetchOptions, Fetcher, FileSink, Object};
//!
//! let fetcher = Fetcher::new(FetchOptions::default())?;
//! let object = Object::new(
//!     "https://example.org/reads.fastq.gz",
//!     1_234_567,
//!     Checksum::Md5 { hex: "0123456789abcdef0123456789abcdef".into() },
//! );
//! let done = fetcher.fetch(object, FileSink::new("reads.fastq.gz")).join()?;
//! println!("{} in {:?}", done.output.display(), done.report.elapsed);
//! # Ok::<(), readkit_fetch::Error>(())
//! ```
//!
//! [`Fetcher::get`] / [`Fetcher::get_with`] read a small body of unknown size
//! (an API response, a manifest) into memory, with a size limit
//! ([`GetOptions`]) and the same retry rules, client and user agent.
//!
//! Every operation returns [`Error`], whose [`ErrorKind`] tells not-found,
//! permanent, throttled, verification and other failures apart.
//!
//! Platforms: positioned file I/O uses `pread`/`pwrite` on unix; elsewhere
//! [`FileSink`] serializes file access behind a lock (correct, slower).

#![deny(missing_docs)]

mod engine;
mod error;
mod http;
mod ledger;
mod object;
mod options;
mod sink;

pub use engine::{Download, Event, FetchHandle, FetchReader, Fetcher, Report, verify_file};
pub use error::{Error, ErrorKind};
pub use object::{Checksum, Object};
pub use options::{FetchOptions, GetOptions, HostLimit};
pub use sink::{Access, FileSink, MemorySink, Opened, Sink, StreamReader, StreamSink};
