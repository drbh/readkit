//! A local HTTP/1.1 range server with fault injection. Offline only.

use std::io::{BufRead, BufReader, Write};
use std::net::{TcpListener, TcpStream};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Duration;

#[derive(Default)]
pub struct Faults {
    /// Offsets whose first response is corrupted.
    pub corrupt_once: Mutex<Vec<u64>>,
    /// Offsets whose first response is a 503.
    pub throttle_once: Mutex<Vec<u64>>,
    /// Offsets whose first response is a 429.
    pub too_many_once: Mutex<Vec<u64>>,
    /// Offsets whose first response stalls for this long before sending.
    pub stall_once: Mutex<Vec<(u64, Duration)>>,
    /// Offsets that stall for this long on every request.
    pub stall_always: Mutex<Vec<(u64, Duration)>>,
    /// Every ranged response stalls this long.
    pub stall_all: Mutex<Option<Duration>>,
    /// Offsets that always fail with 500.
    pub always_fail: Mutex<Vec<u64>>,
    /// Status lines served, in order, to `GET /flaky` before it succeeds.
    pub flaky: Mutex<Vec<&'static str>>,
    /// Requests received, of any kind.
    pub requests: AtomicUsize,
}

pub struct Server {
    pub base: String,
    pub url: String,
    pub served: Arc<AtomicU64>,
    pub faults: Arc<Faults>,
    pub max_active: Arc<AtomicUsize>,
    pub user_agent: Arc<Mutex<Option<String>>>,
    /// `(etag, part_size)` for HEAD responses.
    pub s3: Arc<Mutex<Option<(String, u64)>>>,
}

struct Ctx {
    data: Arc<Vec<u8>>,
    served: Arc<AtomicU64>,
    faults: Arc<Faults>,
    active: Arc<AtomicUsize>,
    max_active: Arc<AtomicUsize>,
    user_agent: Arc<Mutex<Option<String>>>,
    s3: Arc<Mutex<Option<(String, u64)>>>,
}

pub fn serve(data: Arc<Vec<u8>>) -> Server {
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let base = format!("http://{}", listener.local_addr().unwrap());
    let ctx = Arc::new(Ctx {
        data,
        served: Default::default(),
        faults: Default::default(),
        active: Default::default(),
        max_active: Default::default(),
        user_agent: Default::default(),
        s3: Default::default(),
    });
    let server = Server {
        url: format!("{base}/obj"),
        base,
        served: ctx.served.clone(),
        faults: ctx.faults.clone(),
        max_active: ctx.max_active.clone(),
        user_agent: ctx.user_agent.clone(),
        s3: ctx.s3.clone(),
    };
    std::thread::spawn(move || {
        for stream in listener.incoming().flatten() {
            let ctx = ctx.clone();
            std::thread::spawn(move || handle(stream, &ctx));
        }
    });
    server
}

fn take(list: &Mutex<Vec<u64>>, offset: u64) -> bool {
    let mut l = list.lock().unwrap();
    if let Some(i) = l.iter().position(|o| *o == offset) {
        l.remove(i);
        true
    } else {
        false
    }
}

fn status(out: &mut TcpStream, line: &str) {
    let _ = write!(out, "HTTP/1.1 {line}\r\nContent-Length: 0\r\n\r\n");
}

fn handle(stream: TcpStream, ctx: &Ctx) {
    let mut reader = BufReader::new(stream.try_clone().unwrap());
    let mut out = stream;
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
            let lower = line.to_ascii_lowercase();
            if let Some(r) = lower.strip_prefix("range: bytes=") {
                let (a, b) = r.trim().split_once('-').unwrap();
                range = Some((a.parse::<u64>().unwrap(), b.parse::<u64>().unwrap()));
            }
            if lower.starts_with("user-agent:") {
                *ctx.user_agent.lock().unwrap() = Some(line[11..].trim().to_string());
            }
        }
        let n = ctx.active.fetch_add(1, Ordering::SeqCst) + 1;
        ctx.max_active.fetch_max(n, Ordering::SeqCst);
        let keep = respond(&mut out, &request, range, ctx);
        ctx.active.fetch_sub(1, Ordering::SeqCst);
        if !keep {
            return;
        }
    }
}

/// Returns whether the connection is still usable.
fn respond(out: &mut TcpStream, request: &str, range: Option<(u64, u64)>, ctx: &Ctx) -> bool {
    let data = &ctx.data;
    let faults = &ctx.faults;
    faults.requests.fetch_add(1, Ordering::SeqCst);
    if request.contains("/flaky") {
        let next = {
            let mut l = faults.flaky.lock().unwrap();
            (!l.is_empty()).then(|| l.remove(0))
        };
        if let Some(line) = next {
            status(out, line);
            return true;
        }
    }
    if request.contains("/retry-after") {
        let _ = write!(
            out,
            "HTTP/1.1 503 Slow Down\r\nRetry-After: 1\r\nContent-Length: 0\r\n\r\n"
        );
        return true;
    }
    if request.contains("/stall") {
        std::thread::sleep(Duration::from_secs(2));
    }
    if request.contains("/no-length") {
        // No Content-Length: the body ends when the connection closes.
        let _ = write!(out, "HTTP/1.1 200 OK\r\nConnection: close\r\n\r\n");
        let _ = out.write_all(data);
        return false;
    }
    if request.contains("/missing") {
        status(out, "404 Not Found");
        return true;
    }
    if request.contains("/forbidden") {
        status(out, "403 Forbidden");
        return true;
    }
    if request.contains("/throttle503") {
        status(out, "503 Slow Down");
        return true;
    }
    if request.contains("/throttle429") {
        status(out, "429 Too Many Requests");
        return true;
    }
    if request.starts_with("HEAD ") {
        let s3 = ctx.s3.lock().unwrap().clone();
        let (etag, part) = s3.unwrap_or_default();
        let len = if request.contains("partNumber=1") {
            part
        } else {
            data.len() as u64
        };
        let _ = write!(
            out,
            "HTTP/1.1 200 OK\r\nContent-Length: {len}\r\nETag: \"{etag}\"\r\n\r\n"
        );
        return true;
    }
    let Some((a, b)) = range else {
        let _ = write!(
            out,
            "HTTP/1.1 200 OK\r\nContent-Length: {}\r\n\r\n",
            data.len()
        );
        return out.write_all(data).is_ok();
    };
    if faults.always_fail.lock().unwrap().contains(&a) {
        status(out, "500 Internal");
        return true;
    }
    if take(&faults.throttle_once, a) {
        status(out, "503 Slow Down");
        return true;
    }
    if take(&faults.too_many_once, a) {
        status(out, "429 Too Many Requests");
        return true;
    }
    let stall = {
        let mut l = faults.stall_once.lock().unwrap();
        l.iter().position(|(o, _)| *o == a).map(|i| l.remove(i).1)
    };
    let stall = stall
        .or_else(|| {
            let l = faults.stall_always.lock().unwrap();
            l.iter().find(|(o, _)| *o == a).map(|(_, d)| *d)
        })
        .or(*faults.stall_all.lock().unwrap());
    if let Some(d) = stall {
        std::thread::sleep(d);
    }
    let mut body = data[a as usize..=b as usize].to_vec();
    if take(&faults.corrupt_once, a) {
        body[0] ^= 0xff;
    }
    let head = format!(
        "HTTP/1.1 206 Partial Content\r\nContent-Length: {}\r\nContent-Range: bytes {a}-{b}/{}\r\n\r\n",
        body.len(),
        data.len()
    );
    if out.write_all(head.as_bytes()).is_err() || out.write_all(&body).is_err() {
        return false;
    }
    ctx.served.fetch_add(body.len() as u64, Ordering::SeqCst);
    true
}

pub fn data(len: usize) -> Arc<Vec<u8>> {
    let mut v = Vec::with_capacity(len);
    let mut x: u64 = 0x9e3779b97f4a7c15;
    while v.len() < len {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        v.extend_from_slice(&x.to_le_bytes());
    }
    v.truncate(len);
    Arc::new(v)
}
