//! Fault-injection tests against a local HTTP range server (offline).

mod server;

use md5::{Digest as _, Md5};
use readkit_fetch::{
    Checksum, ErrorKind, Event, FetchOptions, Fetcher, FileSink, GetOptions, HostLimit, MemorySink,
    Object, StreamSink, verify_file,
};
use server::{data, serve};
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

const MIB: usize = 1024 * 1024;

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

fn options() -> FetchOptions {
    FetchOptions {
        chunk_size: MIB as u64,
        initial_connections: 4,
        no_progress_timeout: Duration::from_secs(6),
        ..Default::default()
    }
}

fn fetcher(options: FetchOptions) -> Fetcher {
    Fetcher::new(options).unwrap()
}

fn md5_object(url: &str, bytes: &[u8]) -> Object {
    Object::new(
        url,
        bytes.len() as u64,
        Checksum::Md5 {
            hex: hex(&Md5::digest(bytes)),
        },
    )
}

fn s3_etag(bytes: &[u8], part: usize) -> (String, usize) {
    let mut h = Md5::new();
    let parts = bytes.chunks(part).collect::<Vec<_>>();
    for p in &parts {
        h.update(Md5::digest(p));
    }
    (
        format!("{}-{}", hex(&h.finalize()), parts.len()),
        parts.len(),
    )
}

fn s3_object(url: &str, bytes: &[u8], part: usize) -> Object {
    let (etag, parts) = s3_etag(bytes, part);
    Object::new(
        url,
        bytes.len() as u64,
        Checksum::S3MultipartEtag {
            etag,
            parts: parts as u32,
            part_size: part as u64,
        },
    )
}

/// A per-test scratch directory, removed when the test ends (even on panic).
struct Scratch(PathBuf);

impl std::ops::Deref for Scratch {
    type Target = Path;
    fn deref(&self) -> &Path {
        &self.0
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn tmp(name: &str) -> Scratch {
    let dir = std::env::temp_dir().join(format!("readkit-fetch-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    Scratch(dir)
}

fn to_file(
    object: Object,
    dir: &Path,
    options: FetchOptions,
) -> Result<readkit_fetch::Download<PathBuf>, readkit_fetch::Error> {
    fetcher(options)
        .fetch(object, FileSink::new(dir.join("obj.bin")))
        .join()
}

fn served(server: &server::Server) -> u64 {
    server.served.load(Ordering::SeqCst)
}

// ---- downloads ----

#[test]
fn md5_object_parallel_ranges() {
    let bytes = data(10 * MIB + 123);
    let server = serve(bytes.clone());
    let dir = tmp("md5");
    let done = to_file(md5_object(&server.url, &bytes), &dir, options()).unwrap();
    assert_eq!(done.output, dir.join("obj.bin"));
    assert_eq!(fs::read(dir.join("obj.bin")).unwrap(), *bytes);
    assert_eq!(served(&server), bytes.len() as u64);
    assert!(FileSink::marker_path(&dir.join("obj.bin")).exists());
    assert!(!dir.join("obj.bin.partial").exists());
    assert!(!dir.join("obj.bin.ledger.json").exists());
    // Second call reuses the verified file without any transfer.
    let again = to_file(md5_object(&server.url, &bytes), &dir, options()).unwrap();
    assert!(again.report.reused);
    assert_eq!(served(&server), bytes.len() as u64);
}

#[test]
fn s3_multipart_etag_and_corrupt_part_repair() {
    let part = 2 * MIB;
    let bytes = data(3 * part + 777);
    let server = serve(bytes.clone());
    server.faults.corrupt_once.lock().unwrap().push(part as u64);
    let dir = tmp("s3");
    let opts = FetchOptions {
        chunk_size: 2 * MIB as u64,
        ..options()
    };
    let done = to_file(s3_object(&server.url, &bytes, part), &dir, opts).unwrap();
    // A hedge may deliver the clean copy first; otherwise exactly one part is
    // rewritten. Never more.
    assert!(done.report.repaired_chunks <= 1);
    assert_eq!(fs::read(dir.join("obj.bin")).unwrap(), *bytes);
}

#[test]
fn md5_corrupt_chunk_rewrites_one_chunk() {
    let bytes = data(6 * MIB);
    let server = serve(bytes.clone());
    server
        .faults
        .corrupt_once
        .lock()
        .unwrap()
        .push(3 * MIB as u64);
    let dir = tmp("md5-corrupt");
    let events = Arc::new(Mutex::new(Vec::new()));
    let seen = events.clone();
    let done = fetcher(options())
        .fetch_with(
            md5_object(&server.url, &bytes),
            FileSink::new(dir.join("obj.bin")),
            move |e| {
                if matches!(e, Event::Repairing) {
                    seen.lock().unwrap().push("repairing");
                }
            },
        )
        .join()
        .unwrap();
    assert_eq!(done.report.repaired_chunks, 1);
    assert_eq!(*events.lock().unwrap(), vec!["repairing"]);
    assert_eq!(fs::read(dir.join("obj.bin")).unwrap(), *bytes);
}

#[test]
fn throttle_and_stall_are_retried_or_hedged() {
    let bytes = data(4 * MIB);
    let server = serve(bytes.clone());
    server.faults.throttle_once.lock().unwrap().push(0);
    server
        .faults
        .stall_once
        .lock()
        .unwrap()
        .push((2 * MIB as u64, Duration::from_secs(30)));
    let dir = tmp("stall");
    let t = Instant::now();
    let done = to_file(md5_object(&server.url, &bytes), &dir, options()).unwrap();
    assert!(done.report.chunk_retries >= 1);
    assert!(done.report.hedged >= 1, "stalled chunk should be hedged");
    // The stalled original request is aborted, not waited for.
    assert!(
        t.elapsed() < Duration::from_secs(15),
        "hedge should beat the stall"
    );
    assert_eq!(fs::read(dir.join("obj.bin")).unwrap(), *bytes);
}

#[test]
fn not_found_fails_fast_without_retrying() {
    let bytes = data(MIB);
    let server = serve(bytes.clone());
    let dir = tmp("404");
    let mut object = md5_object(&server.url, &bytes);
    object.url = object.url.replace("/obj", "/missing");
    let t = Instant::now();
    let err = to_file(object, &dir, options()).unwrap_err();
    assert_eq!(err.kind(), ErrorKind::NotFound, "{err}");
    assert_eq!(err.status(), Some(404));
    assert!(t.elapsed() < Duration::from_secs(3), "no retry loop on 404");
}

#[test]
fn resume_refetches_only_missing_chunks() {
    let bytes = data(8 * MIB);
    let server = serve(bytes.clone());
    let bad = 5 * MIB as u64;
    server.faults.always_fail.lock().unwrap().push(bad);
    let dir = tmp("resume");
    let err = to_file(md5_object(&server.url, &bytes), &dir, options()).unwrap_err();
    assert_eq!(err.kind(), ErrorKind::Timeout, "{err}");
    let cause = std::error::Error::source(&err)
        .and_then(|s| s.downcast_ref::<readkit_fetch::Error>())
        .expect("last chunk error as source");
    assert_eq!(cause.status(), Some(500));
    let before = served(&server);
    assert!(before >= 7 * MIB as u64, "all good chunks were fetched");
    assert!(dir.join("obj.bin.ledger.json").exists());
    server.faults.always_fail.lock().unwrap().clear();
    let done = to_file(md5_object(&server.url, &bytes), &dir, options()).unwrap();
    assert_eq!(done.report.resumed_chunks, 7);
    // Only the missing chunk (plus at most one hedged duplicate) is re-fetched.
    let again = served(&server) - before;
    assert!(
        (MIB as u64..=2 * MIB as u64).contains(&again),
        "re-fetched {again} bytes"
    );
    assert_eq!(fs::read(dir.join("obj.bin")).unwrap(), *bytes);
}

// ---- sinks ----

#[test]
fn memory_sink_returns_verified_bytes() {
    let bytes = data(3 * MIB + 5);
    let server = serve(bytes.clone());
    let done = fetcher(options())
        .fetch(md5_object(&server.url, &bytes), MemorySink::new())
        .join()
        .unwrap();
    assert_eq!(done.output, *bytes);
    assert!(!done.report.reused);
}

#[test]
fn memory_sink_without_checksum_and_empty_object() {
    let bytes = data(2 * MIB);
    let server = serve(bytes.clone());
    let f = fetcher(options());
    let object = Object::new(&server.url, bytes.len() as u64, Checksum::None);
    assert_eq!(
        f.fetch(object, MemorySink::new()).join().unwrap().output,
        *bytes
    );
    let empty = Object::new(
        &server.url,
        0,
        Checksum::Md5 {
            hex: hex(&Md5::digest(b"")),
        },
    );
    assert!(
        f.fetch(empty, MemorySink::new())
            .join()
            .unwrap()
            .output
            .is_empty()
    );
}

#[test]
fn memory_sink_repairs_corruption() {
    let bytes = data(4 * MIB);
    let server = serve(bytes.clone());
    server.faults.corrupt_once.lock().unwrap().push(MIB as u64);
    let done = fetcher(options())
        .fetch(md5_object(&server.url, &bytes), MemorySink::new())
        .join()
        .unwrap();
    assert_eq!(done.report.repaired_chunks, 1);
    assert_eq!(done.output, *bytes);
}

#[test]
fn stream_delivers_in_order_while_slow_chunk_outstanding() {
    let bytes = data(6 * MIB);
    let server = serve(bytes.clone());
    // Chunk 2 is slow; hedging is off so it stays outstanding.
    server
        .faults
        .stall_once
        .lock()
        .unwrap()
        .push((2 * MIB as u64, Duration::from_secs(2)));
    let opts = FetchOptions {
        hedge_after: Duration::from_secs(60),
        ..options()
    };
    let f = fetcher(opts);
    let (sink, mut reader) = StreamSink::new(4 * MIB);
    let t = Instant::now();
    let handle = f.fetch(md5_object(&server.url, &bytes), sink);
    let mut head = vec![0u8; 2 * MIB];
    reader.read_exact(&mut head).unwrap();
    assert!(
        t.elapsed() < Duration::from_millis(1500),
        "first chunks arrive before the slow one: {:?}",
        t.elapsed()
    );
    assert!(!handle.is_finished(), "slow chunk still outstanding");
    assert_eq!(head, bytes[..2 * MIB]);
    let mut rest = Vec::new();
    reader.read_to_end(&mut rest).unwrap();
    assert!(t.elapsed() >= Duration::from_millis(1900));
    assert_eq!(rest, bytes[2 * MIB..]);
    assert!(reader.ended());
    handle.join().unwrap();
}

#[test]
fn stream_verification_failure_is_a_read_error_at_the_end() {
    let bytes = data(3 * MIB);
    let server = serve(bytes.clone());
    let object = Object::new(
        &server.url,
        bytes.len() as u64,
        Checksum::Md5 {
            hex: "00000000000000000000000000000000".into(),
        },
    );
    let mut reader = fetcher(options()).stream(object, 2 * MIB);
    let mut got = Vec::new();
    let err = reader.read_to_end(&mut got).unwrap_err();
    // Every byte was delivered before the failure.
    assert_eq!(got, *bytes);
    assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    let inner = err
        .get_ref()
        .and_then(|e| e.downcast_ref::<readkit_fetch::Error>())
        .unwrap();
    assert_eq!(inner.kind(), ErrorKind::Verification);
    assert_eq!(reader.finish().unwrap_err().kind(), ErrorKind::Verification);
}

#[test]
fn stream_s3_multipart_verifies_in_order() {
    let part = 2 * MIB;
    let bytes = data(2 * part + 999);
    let server = serve(bytes.clone());
    let mut reader = fetcher(options()).stream(s3_object(&server.url, &bytes, part), MIB);
    let mut got = Vec::new();
    reader.read_to_end(&mut got).unwrap();
    assert_eq!(got, *bytes);
    let report = reader.finish().unwrap();
    assert_eq!(report.size, bytes.len() as u64);
}

#[test]
fn dropping_the_stream_reader_cancels_the_download() {
    let bytes = data(8 * MIB);
    let server = serve(bytes.clone());
    let f = fetcher(options());
    let (sink, mut reader) = StreamSink::new(MIB);
    let handle = f.fetch(md5_object(&server.url, &bytes), sink);
    let mut first = vec![0u8; 1000];
    reader.read_exact(&mut first).unwrap();
    drop(reader);
    let t = Instant::now();
    let err = handle.join().unwrap_err();
    assert_eq!(err.kind(), ErrorKind::Cancelled);
    assert!(t.elapsed() < Duration::from_secs(2));
}

// ---- cancellation ----

#[test]
fn cancel_stops_promptly() {
    let bytes = data(4 * MIB);
    let server = serve(bytes.clone());
    *server.faults.stall_all.lock().unwrap() = Some(Duration::from_secs(30));
    let dir = tmp("cancel");
    let f = fetcher(options());
    let handle = f.fetch(
        md5_object(&server.url, &bytes),
        FileSink::new(dir.join("obj.bin")),
    );
    std::thread::sleep(Duration::from_millis(300));
    let t = Instant::now();
    handle.cancel();
    let err = handle.join().unwrap_err();
    assert_eq!(err.kind(), ErrorKind::Cancelled);
    assert!(t.elapsed() < Duration::from_secs(1), "{:?}", t.elapsed());

    // Dropping the handle stops it too, without hanging.
    let handle = f.fetch(md5_object(&server.url, &bytes), MemorySink::new());
    std::thread::sleep(Duration::from_millis(300));
    let t = Instant::now();
    drop(handle);
    assert!(t.elapsed() < Duration::from_secs(1), "{:?}", t.elapsed());
}

#[test]
fn cancelled_file_download_resumes() {
    let bytes = data(8 * MIB);
    let server = serve(bytes.clone());
    let slow = 7 * MIB as u64;
    server
        .faults
        .stall_always
        .lock()
        .unwrap()
        .push((slow, Duration::from_secs(30)));
    let dir = tmp("cancel-resume");
    let opts = FetchOptions {
        hedge_after: Duration::from_secs(60),
        ..options()
    };
    let done = Arc::new(AtomicUsize::new(0));
    let counter = done.clone();
    let f = fetcher(opts.clone());
    let handle = f.fetch_with(
        md5_object(&server.url, &bytes),
        FileSink::new(dir.join("obj.bin")),
        move |e| {
            if matches!(e, Event::ChunkDone { .. }) {
                counter.fetch_add(1, Ordering::SeqCst);
            }
        },
    );
    let t = Instant::now();
    while done.load(Ordering::SeqCst) < 7 {
        assert!(t.elapsed() < Duration::from_secs(10));
        std::thread::sleep(Duration::from_millis(20));
    }
    handle.cancel();
    assert_eq!(handle.join().unwrap_err().kind(), ErrorKind::Cancelled);
    assert!(!dir.join("obj.bin").exists());
    server.faults.stall_always.lock().unwrap().clear();
    let before = served(&server);
    let done = fetcher(opts)
        .fetch(
            md5_object(&server.url, &bytes),
            FileSink::new(dir.join("obj.bin")),
        )
        .join()
        .unwrap();
    assert_eq!(done.report.resumed_chunks, 7);
    assert_eq!(served(&server) - before, MIB as u64);
    assert_eq!(fs::read(dir.join("obj.bin")).unwrap(), *bytes);
}

// ---- typed errors ----

#[test]
fn forbidden_is_permanent_and_fails_fast() {
    let bytes = data(MIB);
    let server = serve(bytes.clone());
    let object = Object::new(
        format!("{}/forbidden", server.base),
        bytes.len() as u64,
        Checksum::None,
    );
    let t = Instant::now();
    let err = fetcher(options())
        .fetch(object, MemorySink::new())
        .join()
        .unwrap_err();
    assert_eq!(err.kind(), ErrorKind::Permanent, "{err}");
    assert_eq!(err.status(), Some(403));
    assert!(t.elapsed() < Duration::from_secs(3));
}

#[test]
fn throttling_is_typed_and_retried() {
    let bytes = data(2 * MIB);
    let server = serve(bytes.clone());
    server.faults.throttle_once.lock().unwrap().push(0);
    server.faults.too_many_once.lock().unwrap().push(MIB as u64);
    let seen = Arc::new(Mutex::new(Vec::new()));
    let log = seen.clone();
    let done = fetcher(options())
        .fetch_with(
            md5_object(&server.url, &bytes),
            MemorySink::new(),
            move |e| {
                if let Event::Retrying { error, .. } = e {
                    log.lock().unwrap().push((error.kind(), error.status()));
                }
            },
        )
        .join()
        .unwrap();
    assert_eq!(done.output, *bytes);
    let mut seen = seen.lock().unwrap().clone();
    seen.sort_by_key(|s| s.1);
    assert_eq!(
        seen,
        vec![
            (ErrorKind::Throttled, Some(429)),
            (ErrorKind::Throttled, Some(503))
        ]
    );
}

#[test]
fn persistent_throttling_times_out_with_throttled_cause() {
    let bytes = data(MIB);
    let server = serve(bytes.clone());
    let opts = FetchOptions {
        no_progress_timeout: Duration::from_secs(1),
        ..options()
    };
    for (path, status) in [("throttle503", 503), ("throttle429", 429)] {
        let object = Object::new(
            format!("{}/{path}", server.base),
            bytes.len() as u64,
            Checksum::None,
        );
        let err = fetcher(opts.clone())
            .fetch(object, MemorySink::new())
            .join()
            .unwrap_err();
        assert_eq!(err.kind(), ErrorKind::Timeout, "{err}");
        let cause = std::error::Error::source(&err)
            .and_then(|s| s.downcast_ref::<readkit_fetch::Error>())
            .unwrap();
        assert_eq!(cause.kind(), ErrorKind::Throttled);
        assert_eq!(cause.status(), Some(status));
    }
}

// ---- S3 ----

#[test]
fn s3_wrong_etag_fails_verification_after_repair() {
    let part = MIB;
    let bytes = data(3 * part);
    let server = serve(bytes.clone());
    let object = Object::new(
        &server.url,
        bytes.len() as u64,
        Checksum::S3MultipartEtag {
            etag: "00000000000000000000000000000000-3".into(),
            parts: 3,
            part_size: part as u64,
        },
    );
    let err = fetcher(options())
        .fetch(object, MemorySink::new())
        .join()
        .unwrap_err();
    assert_eq!(err.kind(), ErrorKind::Verification, "{err}");
    assert!(
        err.to_string().contains("after repairing 0 chunk(s)"),
        "{err}"
    );
}

#[test]
fn probe_s3_discovers_parts_and_part_size() {
    let part = MIB;
    let bytes = data(2 * part + 10);
    let server = serve(bytes.clone());
    let (etag, _) = s3_etag(&bytes, part);
    *server.s3.lock().unwrap() = Some((etag.clone(), part as u64));
    let f = fetcher(options());
    let object = f.probe_s3(&server.url).unwrap();
    assert_eq!(object.size, bytes.len() as u64);
    assert_eq!(
        object.checksum,
        Checksum::S3MultipartEtag {
            etag,
            parts: 3,
            part_size: part as u64
        }
    );
    let done = f.fetch(object, MemorySink::new()).join().unwrap();
    assert_eq!(done.output, *bytes);

    // Single-part ETag: the plain MD5.
    let md5 = hex(&Md5::digest(&*bytes));
    *server.s3.lock().unwrap() = Some((md5.clone(), 0));
    let object = f.probe_s3(&server.url).unwrap();
    assert!(
        matches!(&object.checksum, Checksum::S3MultipartEtag { parts: 1, etag, .. } if *etag == md5)
    );
    assert_eq!(
        f.fetch(object, MemorySink::new()).join().unwrap().output,
        *bytes
    );

    let missing = f.probe_s3(&format!("{}/missing", server.base)).unwrap_err();
    assert_eq!(missing.kind(), ErrorKind::NotFound);
}

// ---- small-body GET ----

fn quick_retries() -> Fetcher {
    fetcher(FetchOptions {
        retry_delay: Duration::from_millis(10),
        ..options()
    })
}

fn get_opts(max_bytes: u64, retries: u32) -> GetOptions {
    GetOptions {
        max_bytes,
        retries,
        ..Default::default()
    }
}

#[test]
fn get_returns_the_whole_body_with_the_user_agent() {
    let bytes = data(300 * 1024);
    let server = serve(bytes.clone());
    let f = fetcher(FetchOptions {
        user_agent: "readkit-get-test/1".into(),
        ..options()
    });
    assert_eq!(f.get(&server.url).unwrap(), *bytes);
    assert_eq!(
        server.user_agent.lock().unwrap().as_deref(),
        Some("readkit-get-test/1")
    );
    // An empty body is fine.
    let empty = serve(Arc::new(Vec::new()));
    assert!(f.get(&empty.url).unwrap().is_empty());
}

#[test]
fn get_enforces_the_size_limit() {
    let bytes = data(100 * 1024);
    let server = serve(bytes.clone());
    let f = quick_retries();
    let n = bytes.len() as u64;
    // Declared length over the limit: refused before reading, not retried.
    let e = f.get_with(&server.url, &get_opts(n - 1, 3)).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::TooLarge);
    assert!(!e.kind().is_retryable());
    assert_eq!(server.faults.requests.load(Ordering::SeqCst), 1);
    assert_eq!(f.get_with(&server.url, &get_opts(n, 0)).unwrap(), *bytes);
    // No Content-Length: the limit applies while the body arrives.
    let no_length = format!("{}/no-length", server.base);
    assert_eq!(
        f.get_with(&no_length, &get_opts(n - 1, 0))
            .unwrap_err()
            .kind(),
        ErrorKind::TooLarge
    );
    assert_eq!(f.get_with(&no_length, &get_opts(n, 0)).unwrap(), *bytes);
}

#[test]
fn get_error_kinds_and_statuses() {
    let server = serve(data(1024));
    let f = quick_retries();
    let at = |path: &str| format!("{}/{path}", server.base);
    for (path, kind, status) in [
        ("missing", ErrorKind::NotFound, 404),
        ("forbidden", ErrorKind::Permanent, 403),
        ("throttle503", ErrorKind::Throttled, 503),
        ("throttle429", ErrorKind::Throttled, 429),
    ] {
        let before = server.faults.requests.load(Ordering::SeqCst);
        let e = f.get_with(&at(path), &get_opts(1 << 20, 2)).unwrap_err();
        assert_eq!((e.kind(), e.status()), (kind, Some(status)), "{path}");
        // Throttling is retried; not-found and permanent errors are not.
        let tries = if kind == ErrorKind::Throttled { 3 } else { 1 };
        let after = server.faults.requests.load(Ordering::SeqCst);
        assert_eq!(after - before, tries, "{path}");
    }
    let e = f.get("not a url").unwrap_err();
    assert_eq!(e.kind(), ErrorKind::InvalidInput);
    // Nothing listens on port 1.
    let e = f
        .get_with("http://127.0.0.1:1/x", &get_opts(1 << 20, 0))
        .unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Transient);
}

#[test]
fn get_retries_throttling_and_server_errors() {
    let bytes = data(4096);
    let server = serve(bytes.clone());
    *server.faults.flaky.lock().unwrap() =
        vec!["503 Slow Down", "500 Internal", "429 Too Many Requests"];
    let f = quick_retries();
    let url = format!("{}/flaky", server.base);
    assert_eq!(f.get_with(&url, &get_opts(1 << 20, 3)).unwrap(), *bytes);
    assert_eq!(server.faults.requests.load(Ordering::SeqCst), 4);

    // One retry too few: the last failure is returned.
    *server.faults.flaky.lock().unwrap() = vec!["500 Internal", "502 Bad Gateway"];
    let e = f.get_with(&url, &get_opts(1 << 20, 1)).unwrap_err();
    assert_eq!((e.kind(), e.status()), (ErrorKind::Transient, Some(502)));
}

#[test]
fn get_honours_retry_after() {
    let server = serve(data(16));
    let f = quick_retries();
    let started = Instant::now();
    let e = f
        .get_with(
            &format!("{}/retry-after", server.base),
            &get_opts(1 << 20, 1),
        )
        .unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Throttled);
    // `Retry-After: 1` beats the 10 ms backoff.
    assert!(started.elapsed() >= Duration::from_secs(1));
    assert_eq!(server.faults.requests.load(Ordering::SeqCst), 2);
}

#[test]
fn get_times_out() {
    let server = serve(data(16));
    let f = quick_retries();
    let opts = GetOptions {
        timeout: Duration::from_millis(200),
        retries: 0,
        ..Default::default()
    };
    let started = Instant::now();
    let e = f
        .get_with(&format!("{}/stall", server.base), &opts)
        .unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Timeout);
    assert!(started.elapsed() < Duration::from_millis(1500));
}

#[test]
fn get_works_inside_an_async_runtime() {
    let bytes = data(2048);
    let server = serve(bytes.clone());
    let f = quick_retries();
    let rt = tokio::runtime::Builder::new_current_thread()
        .build()
        .unwrap();
    let got = rt.block_on(async { f.get(&server.url) }).unwrap();
    assert_eq!(got, *bytes);
}

// ---- options ----

#[test]
fn host_limits_cap_connections() {
    let bytes = data(6 * MIB);
    let server = serve(bytes.clone());
    let opts = FetchOptions {
        host_limits: vec![HostLimit::new("127.0.0.1", 1)],
        ..options()
    };
    let done = fetcher(opts)
        .fetch(md5_object(&server.url, &bytes), MemorySink::new())
        .join()
        .unwrap();
    assert_eq!(done.output, *bytes);
    assert_eq!(server.max_active.load(Ordering::SeqCst), 1);
}

#[test]
fn user_agent_is_configurable() {
    let bytes = data(MIB);
    let server = serve(bytes.clone());
    let object = Object::new(&server.url, bytes.len() as u64, Checksum::None);
    fetcher(options())
        .fetch(object.clone(), MemorySink::new())
        .join()
        .unwrap();
    let ua = server.user_agent.lock().unwrap().clone().unwrap();
    assert!(ua.starts_with("readkit-fetch/"), "{ua}");
    let opts = FetchOptions {
        user_agent: "custom/1".into(),
        ..options()
    };
    fetcher(opts)
        .fetch(object, MemorySink::new())
        .join()
        .unwrap();
    assert_eq!(
        server.user_agent.lock().unwrap().as_deref(),
        Some("custom/1")
    );
}

#[test]
fn verify_file_checks_size_and_checksum() {
    let bytes = data(3 * MIB);
    let dir = tmp("verify");
    let path = dir.join("f");
    fs::write(&path, &*bytes).unwrap();
    verify_file(&md5_object("http://x/o", &bytes), &path).unwrap();
    verify_file(&s3_object("http://x/o", &bytes, MIB), &path).unwrap();
    let mut wrong = md5_object("http://x/o", &bytes);
    wrong.checksum = Checksum::Md5 {
        hex: "00".repeat(16),
    };
    assert_eq!(
        verify_file(&wrong, &path).unwrap_err().kind(),
        ErrorKind::Verification
    );
}

#[test]
fn invalid_url_is_invalid_input_and_reaches_the_stream() {
    let object = Object::new("not a url", 10, Checksum::None);
    let mut reader = fetcher(options()).stream(object, MIB);
    let err = reader.read_to_end(&mut Vec::new()).unwrap_err();
    assert_eq!(err.kind(), std::io::ErrorKind::InvalidInput);
    assert_eq!(reader.finish().unwrap_err().kind(), ErrorKind::InvalidInput);
}
