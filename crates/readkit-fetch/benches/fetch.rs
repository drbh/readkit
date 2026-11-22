//! Download benchmarks against an in-process, shaped S3 stand-in, or live URLs.
//!
//! ```text
//! cargo bench -p readkit-fetch --bench fetch                         # every profile
//! cargo bench -p readkit-fetch --bench fetch -- --profile s3-wan     # one profile
//! cargo bench -p readkit-fetch --bench fetch -- --save main          # bench/baselines/main.json
//! cargo bench -p readkit-fetch --bench fetch -- --compare main       # medians against it
//! cargo bench -p readkit-fetch --bench fetch -- --list               # profiles and cases
//! cargo bench -p readkit-fetch --bench fetch -- \
//!     --url https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR000001/SRR000001 --label aws-batch
//! ```
//!
//! No network: the server runs in this process on 127.0.0.1 and paces its
//! responses. Every request waits one round trip before its first byte, a new
//! connection waits two more (TCP and TLS setup), each connection is capped
//! at `conn` bytes/s and all of them share a `link` of bytes/s. Objects are
//! deterministic bytes served under S3 ETags: a whole-object MD5, and
//! multipart with 16 MiB parts and with 8 MiB parts (as NCBI's SRA bucket
//! has). The fetcher uses its defaults, with the bench host given the S3
//! connection cap.
//!
//! Results: a table, and JSON at `target/fetch-bench/latest.json` (with
//! `--save NAME` also at `bench/baselines/NAME.json`). `ceiling` is the
//! fastest the profile allows, min(link, S3 cap × conn), and `efficiency`
//! the median throughput over it. Live runs (`--url`) write the same records
//! with `profile` set to `--label`, so in-region numbers compare directly.
//!
//! Other options: `--runs N` (default 3), `--case SUBSTR`, `--chunk MiB`,
//! `--conns N` (start and cap), `--hedge-ms N`.

use md5::{Digest, Md5};
use readkit_fetch::{Event, FetchOptions, Fetcher, FileSink, HostLimit, MemorySink, Object};
use serde_json::{Value, json};
use std::io::{BufRead, BufReader, Read, Write};
use std::net::{TcpListener, TcpStream};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

const MIB: u64 = 1 << 20;
/// Read-ahead of stream cases.
const STREAM_BUFFER: usize = 512 << 20;

struct Profile {
    name: &'static str,
    about: &'static str,
    rtt: Duration,
    /// Bytes/s per connection.
    conn: Option<f64>,
    /// Bytes/s shared by all connections.
    link: Option<f64>,
    size: u64,
}

const fn mb(n: f64) -> Option<f64> {
    Some(n * 1e6)
}

const PROFILES: &[Profile] = &[
    Profile {
        name: "loopback",
        about: "no shaping: the client's CPU ceiling",
        rtt: Duration::ZERO,
        conn: None,
        link: None,
        size: 512 * MIB,
    },
    Profile {
        name: "s3-region",
        about: "S3 from EC2 in its region: 2 ms, 90 MB/s per connection, 10 Gbit/s",
        rtt: Duration::from_millis(2),
        conn: mb(90.0),
        link: mb(1250.0),
        size: 512 * MIB,
    },
    Profile {
        name: "conn-capped",
        about: "per-connection limit only: 40 ms, 40 MB/s per connection",
        rtt: Duration::from_millis(40),
        conn: mb(40.0),
        link: None,
        size: 512 * MIB,
    },
    Profile {
        name: "s3-wan",
        about: "S3 over a congested route (measured): 30 ms, 20 MB/s per connection, 47 MB/s",
        rtt: Duration::from_millis(30),
        conn: mb(20.0),
        link: mb(47.0),
        // Large enough for 8 MiB chunks that take longer than `hedge_after`
        // on a saturated link (a fixed hedge threshold then duplicates them).
        size: 256 * MIB,
    },
];

#[derive(Clone, Copy, PartialEq)]
enum SinkKind {
    Memory,
    File,
    Stream,
}

impl SinkKind {
    fn name(self) -> &'static str {
        match self {
            SinkKind::Memory => "memory",
            SinkKind::File => "file",
            SinkKind::Stream => "stream",
        }
    }
}

/// (case name, object path, sink)
const CASES: &[(&str, &str, SinkKind)] = &[
    ("mpu16/memory", "mpu16", SinkKind::Memory),
    ("md5/memory", "md5", SinkKind::Memory),
    ("mpu8/file", "mpu8", SinkKind::File),
    ("mpu16/stream", "mpu16", SinkKind::Stream),
];

struct Args {
    profiles: Vec<String>,
    cases: Option<String>,
    urls: Vec<String>,
    label: String,
    runs: usize,
    save: Option<String>,
    compare: Option<String>,
    chunk: Option<u64>,
    conns: Option<usize>,
    hedge_ms: Option<u64>,
}

fn args() -> Args {
    let mut a = Args {
        profiles: Vec::new(),
        cases: None,
        urls: Vec::new(),
        label: "live".into(),
        runs: 3,
        save: None,
        compare: None,
        chunk: None,
        conns: None,
        hedge_ms: None,
    };
    let mut it = std::env::args().skip(1);
    while let Some(arg) = it.next() {
        let mut val = || it.next().unwrap_or_else(|| panic!("{arg} needs a value"));
        match arg.as_str() {
            "--profile" => a.profiles.push(val()),
            "--case" => a.cases = Some(val()),
            "--url" => a.urls.push(val()),
            "--label" => a.label = val(),
            "--runs" => a.runs = val().parse().expect("--runs N"),
            "--save" => a.save = Some(val()),
            "--compare" => a.compare = Some(val()),
            "--chunk" => a.chunk = Some(val().parse::<u64>().expect("--chunk MiB") * MIB),
            "--conns" => a.conns = Some(val().parse().expect("--conns N")),
            "--hedge-ms" => a.hedge_ms = Some(val().parse().expect("--hedge-ms N")),
            "--list" => {
                for p in PROFILES {
                    println!("{:12} {} ({} MiB)", p.name, p.about, p.size / MIB);
                }
                for (c, _, _) in CASES {
                    println!("case {c}");
                }
                std::process::exit(0);
            }
            // `cargo bench` passes `--bench`.
            "--bench" => {}
            _ => panic!("unknown argument {arg} (see --list and the file's docs)"),
        }
    }
    a
}

fn options(a: &Args, host: &str) -> FetchOptions {
    let mut o = FetchOptions::default();
    // The bench host gets what S3 gets.
    let cap = a
        .conns
        .unwrap_or_else(|| o.host_max_connections("bucket.s3.amazonaws.com"));
    o.host_limits.push(HostLimit::new(host, cap));
    if let Some(c) = a.conns {
        o.initial_connections = c;
        o.max_connections = o.max_connections.max(c);
    }
    if let Some(c) = a.chunk {
        o.chunk_size = c;
    }
    if let Some(h) = a.hedge_ms {
        o.hedge_after = Duration::from_millis(h);
    }
    o
}

fn main() {
    let a = args();
    let mut results = Vec::new();
    if a.urls.is_empty() {
        let selected: Vec<&Profile> = PROFILES
            .iter()
            .filter(|p| a.profiles.is_empty() || a.profiles.iter().any(|n| n == p.name))
            .collect();
        assert!(!selected.is_empty(), "no such profile (see --list)");
        for p in selected {
            println!("== {}: {}", p.name, p.about);
            let server = Server::start(p);
            let opts = options(&a, "127.0.0.1");
            let cap = opts.host_max_connections("127.0.0.1") as f64;
            let ceiling = match (p.conn, p.link) {
                (Some(c), Some(l)) => Some(l.min(c * cap)),
                (Some(c), None) => Some(c * cap),
                (None, l) => l,
            };
            for &(case, path, sink) in CASES {
                if a.cases.as_ref().is_some_and(|f| !case.contains(f.as_str())) {
                    continue;
                }
                let url = format!("{}/{path}", server.base);
                results.push(bench(&a, p.name, case, &url, sink, &opts, ceiling));
            }
        }
    } else {
        println!("== {} (live)", a.label);
        for url in &a.urls {
            let host = reqwest_host(url);
            let opts = options(&a, &host);
            let name = url.rsplit('/').next().unwrap_or(url);
            for sink in [SinkKind::Memory, SinkKind::File, SinkKind::Stream] {
                let case = format!("{name}/{}", sink.name());
                if a.cases.as_ref().is_some_and(|f| !case.contains(f.as_str())) {
                    continue;
                }
                results.push(bench(&a, &a.label, &case, url, sink, &opts, None));
            }
        }
    }
    report(&a, results);
}

fn reqwest_host(url: &str) -> String {
    let rest = url.split("://").nth(1).unwrap_or(url);
    let authority = rest.split('/').next().unwrap_or(rest);
    authority.split(':').next().unwrap_or(authority).to_string()
}

#[derive(Default)]
struct Trace {
    first_byte: Option<Instant>,
    last_chunk: Option<Instant>,
}

fn bench(
    a: &Args,
    profile: &str,
    case: &str,
    url: &str,
    sink: SinkKind,
    opts: &FetchOptions,
    ceiling: Option<f64>,
) -> Value {
    let mut runs = Vec::new();
    let (mut ttfb, mut finish, mut retries, mut hedged) = (Vec::new(), Vec::new(), 0, 0);
    let mut size = 0;
    for _ in 0..a.runs {
        // A fresh fetcher per run: every run pays connection setup, as the
        // first download of a job does.
        let fetcher = Fetcher::new(opts.clone()).expect("fetcher");
        let object = fetcher.probe_s3(url).expect("probe");
        size = object.size;
        let (secs, fb, fin, report) = run(&fetcher, object, sink);
        runs.push(secs);
        ttfb.push(fb);
        finish.push(fin);
        retries += report.chunk_retries;
        hedged += report.hedged;
    }
    let med = median(&runs);
    let mib_s = size as f64 / MIB as f64 / med;
    let ceiling_mib_s = ceiling.map(|c| c / MIB as f64);
    let efficiency = ceiling_mib_s.map(|c| mib_s / c);
    println!(
        "  {case:18} {:7.3}s {:8.1} MiB/s{} | ttfb {:5.0}ms finish {:5.3}s | spread {:4.1}% | retries {retries} hedged {hedged}",
        med,
        mib_s,
        efficiency.map_or(String::new(), |e| format!(
            " ({:3.0}% of ceiling)",
            e * 100.0
        )),
        median(&ttfb) * 1e3,
        median(&finish),
        (runs.iter().cloned().fold(0.0, f64::max) - runs.iter().cloned().fold(f64::MAX, f64::min))
            / med
            * 100.0,
    );
    json!({
        "profile": profile,
        "case": case,
        "sink": sink.name(),
        "size": size,
        "runs_s": runs,
        "median_s": med,
        "mib_s": mib_s,
        "ceiling_mib_s": ceiling_mib_s,
        "efficiency": efficiency,
        "ttfb_ms": median(&ttfb) * 1e3,
        "finish_s": median(&finish),
        "retries": retries,
        "hedged": hedged,
    })
}

/// One download: (total s, first byte s, finish s, report).
fn run(
    fetcher: &Fetcher,
    object: Object,
    sink: SinkKind,
) -> (f64, f64, f64, readkit_fetch::Report) {
    let trace = Arc::new(Mutex::new(Trace::default()));
    let tr = trace.clone();
    let observer = move |e: &Event| {
        let now = Instant::now();
        let mut t = tr.lock().unwrap();
        match e {
            Event::Received { .. } => {
                t.first_byte.get_or_insert(now);
            }
            Event::ChunkDone { .. } => t.last_chunk = Some(now),
            _ => {}
        }
    };
    let t0 = Instant::now();
    let report = match sink {
        SinkKind::Memory => {
            let d = fetcher
                .fetch_with(object, MemorySink::new(), observer)
                .join()
                .expect("download");
            d.report
        }
        SinkKind::File => {
            let dest =
                std::env::temp_dir().join(format!("readkit-fetch-bench-{}", std::process::id()));
            let s = FileSink::new(&dest).verified_marker(false);
            let d = fetcher.fetch_with(object, s, observer).join();
            let _ = std::fs::remove_file(&dest);
            d.expect("download").report
        }
        SinkKind::Stream => {
            let mut reader = fetcher.stream(object, STREAM_BUFFER);
            let mut buf = vec![0u8; 1 << 20];
            loop {
                let n = reader.read(&mut buf).expect("read");
                let now = Instant::now();
                let mut t = trace.lock().unwrap();
                t.first_byte.get_or_insert(now);
                if n == 0 {
                    break;
                }
                t.last_chunk = Some(now);
            }
            reader.finish().expect("download")
        }
    };
    let total = t0.elapsed();
    let t = trace.lock().unwrap();
    let fb = t.first_byte.map_or(total, |f| f - t0);
    let lc = t.last_chunk.map_or(total, |l| l - t0);
    (
        total.as_secs_f64(),
        fb.as_secs_f64(),
        (total - lc).as_secs_f64(),
        report,
    )
}

fn median(v: &[f64]) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(f64::total_cmp);
    s[s.len() / 2]
}

fn workspace() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn git(args: &[&str]) -> String {
    std::process::Command::new("git")
        .args(args)
        .current_dir(workspace())
        .output()
        .ok()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_default()
}

fn report(a: &Args, results: Vec<Value>) {
    let dirty = !git(&["status", "--porcelain", "--untracked-files=no"]).is_empty();
    let doc = json!({
        "schema": 1,
        "git": format!("{}{}", git(&["rev-parse", "--short", "HEAD"]), if dirty { "-dirty" } else { "" }),
        "unix_time": std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| d.as_secs()),
        "machine": {
            "os": std::env::consts::OS,
            "arch": std::env::consts::ARCH,
            "cpus": std::thread::available_parallelism().map_or(0, |n| n.get()),
        },
        "runs": a.runs,
        "results": results,
    });
    let text = serde_json::to_string_pretty(&doc).unwrap();
    let out = workspace().join("target/fetch-bench");
    std::fs::create_dir_all(&out).unwrap();
    std::fs::write(out.join("latest.json"), &text).unwrap();
    if let Some(name) = &a.save {
        let dir = workspace().join("bench/baselines");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(format!("{name}.json"));
        std::fs::write(&path, &text).unwrap();
        println!("saved {}", path.display());
    }
    if let Some(name) = &a.compare {
        let path = workspace().join(format!("bench/baselines/{name}.json"));
        let base: Value = serde_json::from_slice(&std::fs::read(&path).expect("baseline"))
            .expect("baseline JSON");
        println!(
            "== against {name} ({})",
            base["git"].as_str().unwrap_or("?")
        );
        for r in doc["results"].as_array().unwrap() {
            let old = base["results"].as_array().into_iter().flatten().find(|b| {
                b["profile"] == r["profile"] && b["case"] == r["case"] && b["size"] == r["size"]
            });
            let new = r["median_s"].as_f64().unwrap();
            // Wasted requests (hedges, retries) change before wall time does:
            // on the stand-in, duplicates share the link without TCP contention.
            let waste =
                |v: &Value| v["hedged"].as_u64().unwrap_or(0) + v["retries"].as_u64().unwrap_or(0);
            match old.and_then(|o| Some((o["median_s"].as_f64()?, waste(o)))) {
                Some((old, old_waste)) => println!(
                    "  {:12} {:18} {:7.3}s -> {:7.3}s  {:+6.1}%{}",
                    r["profile"].as_str().unwrap(),
                    r["case"].as_str().unwrap(),
                    old,
                    new,
                    (new / old - 1.0) * 100.0,
                    if waste(r) != old_waste {
                        format!("  hedges+retries {old_waste} -> {}", waste(r))
                    } else {
                        String::new()
                    }
                ),
                None => println!(
                    "  {:12} {:18} new",
                    r["profile"].as_str().unwrap(),
                    r["case"].as_str().unwrap()
                ),
            }
        }
    }
}

// ---------------------------------------------------------------------------
// The shaped server.

struct Served {
    data: Vec<u8>,
    /// path -> (ETag, part size; 0 for a single part)
    etags: Vec<(&'static str, String, u64)>,
}

struct Shaping {
    rtt: Duration,
    conn: Option<f64>,
    link: Option<f64>,
    /// When the link is free next.
    link_clock: Mutex<Instant>,
}

struct Server {
    base: String,
}

impl Server {
    fn start(p: &Profile) -> Server {
        let data = bytes(p.size as usize);
        let hex = |d: &[u8]| d.iter().map(|b| format!("{b:02x}")).collect::<String>();
        let multipart = |part: u64| {
            let mut all = Md5::new();
            let parts = data.chunks(part as usize);
            let n = parts.len();
            for c in parts {
                all.update(Md5::digest(c));
            }
            format!("{}-{n}", hex(&all.finalize()))
        };
        let served = Arc::new(Served {
            etags: vec![
                ("md5", hex(&Md5::digest(&data)), 0),
                ("mpu16", multipart(16 * MIB), 16 * MIB),
                ("mpu8", multipart(8 * MIB), 8 * MIB),
            ],
            data,
        });
        let shaping = Arc::new(Shaping {
            rtt: p.rtt,
            conn: p.conn,
            link: p.link,
            link_clock: Mutex::new(Instant::now()),
        });
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let base = format!("http://{}", listener.local_addr().unwrap());
        std::thread::spawn(move || {
            for stream in listener.incoming().flatten() {
                let (served, shaping) = (served.clone(), shaping.clone());
                std::thread::spawn(move || connection(stream, &served, &shaping));
            }
        });
        Server { base }
    }
}

/// Deterministic, incompressible bytes.
fn bytes(len: usize) -> Vec<u8> {
    let mut v = Vec::with_capacity(len + 8);
    let mut x: u64 = 0x9e3779b97f4a7c15;
    while v.len() < len {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        v.extend_from_slice(&x.to_le_bytes());
    }
    v.truncate(len);
    v
}

fn connection(stream: TcpStream, served: &Served, shaping: &Shaping) {
    let _ = stream.set_nodelay(true);
    let mut reader = BufReader::new(stream.try_clone().unwrap());
    let mut out = stream;
    // TCP and TLS setup: two round trips before the first request lands.
    std::thread::sleep(shaping.rtt * 2);
    let mut conn_clock = Instant::now();
    loop {
        let mut request = String::new();
        if reader.read_line(&mut request).unwrap_or(0) == 0 {
            return;
        }
        let mut range = None;
        loop {
            let mut line = String::new();
            if reader.read_line(&mut line).unwrap_or(0) == 0 {
                return;
            }
            if line == "\r\n" {
                break;
            }
            if let Some(r) = line.to_ascii_lowercase().strip_prefix("range: bytes=") {
                let (a, b) = r.trim().split_once('-').unwrap();
                range = Some((a.parse::<u64>().unwrap(), b.parse::<u64>().unwrap()));
            }
        }
        // The request travels half a round trip, the response the other half.
        std::thread::sleep(shaping.rtt);
        let path = request.split_whitespace().nth(1).unwrap_or("/");
        let (name, query) = path
            .trim_start_matches('/')
            .split_once('?')
            .unwrap_or((path.trim_start_matches('/'), ""));
        let Some((_, etag, part)) = served.etags.iter().find(|(n, _, _)| *n == name) else {
            let _ = out.write_all(b"HTTP/1.1 404 Not Found\r\nContent-Length: 0\r\n\r\n");
            continue;
        };
        let size = served.data.len() as u64;
        if request.starts_with("HEAD ") {
            let len = if query.contains("partNumber=1") && *part > 0 {
                *part
            } else {
                size
            };
            let _ = write!(
                out,
                "HTTP/1.1 200 OK\r\nContent-Length: {len}\r\nETag: \"{etag}\"\r\n\r\n"
            );
            continue;
        }
        let (a, b) = range.unwrap_or((0, size - 1));
        let head = format!(
            "HTTP/1.1 206 Partial Content\r\nContent-Length: {}\r\nContent-Range: bytes {a}-{b}/{size}\r\n\r\n",
            b - a + 1
        );
        if out.write_all(head.as_bytes()).is_err()
            || send(
                &mut out,
                &served.data[a as usize..=b as usize],
                shaping,
                &mut conn_clock,
            )
            .is_err()
        {
            return;
        }
    }
}

/// Write `body` paced by the connection's and the link's rates.
fn send(
    out: &mut TcpStream,
    body: &[u8],
    s: &Shaping,
    conn_clock: &mut Instant,
) -> std::io::Result<()> {
    if s.conn.is_none() && s.link.is_none() {
        return out.write_all(body);
    }
    // Slices of about 1 ms at the connection's rate, at least 16 KiB.
    let slice = s
        .conn
        .map_or(256 << 10, |r| ((r / 1000.0) as usize).max(16 << 10));
    for piece in body.chunks(slice) {
        let n = piece.len() as f64;
        if let Some(rate) = s.conn {
            let now = Instant::now();
            *conn_clock = (*conn_clock).max(now) + Duration::from_secs_f64(n / rate);
            sleep_until(*conn_clock - Duration::from_secs_f64(n / rate));
        }
        if let Some(rate) = s.link {
            let start = {
                let mut c = s.link_clock.lock().unwrap();
                let start = (*c).max(Instant::now());
                *c = start + Duration::from_secs_f64(n / rate);
                start
            };
            sleep_until(start);
        }
        out.write_all(piece)?;
    }
    Ok(())
}

fn sleep_until(t: Instant) {
    let now = Instant::now();
    if t > now {
        std::thread::sleep(t - now);
    }
}
