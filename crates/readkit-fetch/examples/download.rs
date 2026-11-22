//! Downloads an S3 object to a file over parallel ranges, verified against
//! its ETag.
//!
//! ```text
//! cargo run --release -p readkit-fetch --example download -- URL [OUT]
//! cargo run --release -p readkit-fetch --example download -- \
//!     https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR000001/SRR000001
//! ```
//!
//! OUT defaults to the last segment of the URL. An interrupted download
//! resumes from the `.partial` file next to OUT.

use readkit_fetch::{FetchOptions, Fetcher, FileSink};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let url = args.next().ok_or("usage: download URL [OUT]")?;
    let out = match args.next() {
        Some(out) => out,
        None => url
            .rsplit('/')
            .next()
            .filter(|name| !name.is_empty())
            .ok_or("no file name in the URL; pass OUT")?
            .to_string(),
    };

    let fetcher = Fetcher::new(FetchOptions::default())?;
    // Size and ETag from a HEAD request.
    let object = fetcher.probe_s3(&url)?;
    eprintln!("{} bytes, {:?}", object.size, object.checksum);

    let done = fetcher.fetch(object, FileSink::new(&out)).join()?;
    let r = &done.report;
    let secs = r.elapsed.as_secs_f64();
    eprintln!(
        "{}: {} bytes in {secs:.1} s ({:.0} MB/s), {} chunks resumed, {} retried",
        done.output.display(),
        r.size,
        r.size as f64 / secs.max(1e-9) / 1e6,
        r.resumed_chunks,
        r.chunk_retries,
    );
    Ok(())
}
