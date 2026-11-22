//! Ranged GETs, small whole-body GETs and HEAD probes on reqwest's async client, driven from
//! blocking threads so every request can be cancelled at once.

use crate::error::{Error, ErrorKind};
use std::future::Future;
use std::pin::pin;
use std::task::Poll;
use std::time::Duration;
use tokio::sync::watch;

/// Map an HTTP status to an error kind, or `None` if it is a success.
pub(crate) fn classify(status: u16) -> Option<ErrorKind> {
    match status {
        200..=299 => None,
        404 | 410 => Some(ErrorKind::NotFound),
        429 | 503 => Some(ErrorKind::Throttled),
        408 => Some(ErrorKind::Transient),
        400..=499 => Some(ErrorKind::Permanent),
        _ => Some(ErrorKind::Transient),
    }
}

fn status_error(status: reqwest::StatusCode, url: &str) -> Error {
    let kind = classify(status.as_u16()).unwrap_or(ErrorKind::Transient);
    let what = match kind {
        ErrorKind::NotFound => "not found",
        ErrorKind::Permanent => "permanent",
        ErrorKind::Throttled => "throttled",
        _ => "transient",
    };
    Error::new(kind, format!("{what}: HTTP {status} for {url}")).with_status(status.as_u16())
}

pub(crate) fn request_error(e: reqwest::Error, what: &str) -> Error {
    let kind = if e.is_timeout() {
        ErrorKind::Timeout
    } else if e.is_builder() {
        ErrorKind::InvalidInput
    } else {
        ErrorKind::Transient
    };
    Error::new(kind, format!("{what}: {e}")).with_source(e)
}

/// Run `fut` unless `cancel` becomes true first.
pub(crate) async fn or_cancel<F: Future>(
    fut: F,
    cancel: &mut watch::Receiver<bool>,
) -> Option<F::Output> {
    let mut fut = pin!(fut);
    let mut stop = pin!(cancel.wait_for(|v| *v));
    std::future::poll_fn(|cx| {
        if let Poll::Ready(out) = fut.as_mut().poll(cx) {
            return Poll::Ready(Some(out));
        }
        match stop.as_mut().poll(cx) {
            // Sender dropped counts as cancelled too.
            Poll::Ready(_) => Poll::Ready(None),
            Poll::Pending => Poll::Pending,
        }
    })
    .await
}

/// Ranged GET of `len` bytes at `offset`. `progress` sees every piece of the
/// body as it arrives.
pub(crate) async fn fetch_range(
    client: &reqwest::Client,
    url: &str,
    offset: u64,
    len: u64,
    timeout: Duration,
    progress: &mut (dyn FnMut(u64) + Send),
) -> Result<Vec<u8>, Error> {
    if len == 0 {
        return Ok(Vec::new());
    }
    let mut response = client
        .get(url)
        .header("Range", format!("bytes={}-{}", offset, offset + len - 1))
        .timeout(timeout)
        .send()
        .await
        .map_err(|e| request_error(e, "request failed"))?;
    let status = response.status();
    if status.as_u16() != 206 {
        if classify(status.as_u16()).is_some() {
            return Err(status_error(status, url));
        }
        return Err(Error::new(
            ErrorKind::Transient,
            format!("unexpected HTTP {status} for ranged GET of {url}"),
        )
        .with_status(status.as_u16()));
    }
    let mut buffer = Vec::with_capacity(len as usize);
    while let Some(piece) = response
        .chunk()
        .await
        .map_err(|e| request_error(e, "body interrupted"))?
    {
        if buffer.len() as u64 + piece.len() as u64 > len {
            return Err(Error::new(
                ErrorKind::Transient,
                format!("long body: more than {len} bytes"),
            ));
        }
        buffer.extend_from_slice(&piece);
        progress(piece.len() as u64);
    }
    if buffer.len() as u64 != len {
        return Err(Error::new(
            ErrorKind::Transient,
            format!("short body: got {} of {len} bytes", buffer.len()),
        ));
    }
    Ok(buffer)
}

/// Plain GET of a whole body of at most `max_bytes`. On a throttling status
/// the error carries the `Retry-After` delay in seconds, if the server sent one.
pub(crate) async fn get(
    client: &reqwest::Client,
    url: &str,
    max_bytes: u64,
    timeout: Duration,
) -> Result<Vec<u8>, (Error, Option<Duration>)> {
    let mut response = client
        .get(url)
        .timeout(timeout)
        .send()
        .await
        .map_err(|e| (request_error(e, "GET failed"), None))?;
    let status = response.status();
    if classify(status.as_u16()).is_some() {
        let retry_after = response
            .headers()
            .get("retry-after")
            .and_then(|v| v.to_str().ok())
            .and_then(|v| v.trim().parse::<u64>().ok())
            .map(Duration::from_secs);
        return Err((status_error(status, url), retry_after));
    }
    let too_large = |what: String| {
        (
            Error::new(
                ErrorKind::TooLarge,
                format!("GET {url}: {what}, over the {max_bytes}-byte limit"),
            ),
            None,
        )
    };
    let declared = response.content_length();
    if let Some(len) = declared.filter(|&len| len > max_bytes) {
        return Err(too_large(format!("Content-Length {len}")));
    }
    let mut buffer = Vec::with_capacity(declared.unwrap_or(0) as usize);
    while let Some(piece) = response
        .chunk()
        .await
        .map_err(|e| (request_error(e, "body interrupted"), None))?
    {
        if buffer.len() as u64 + piece.len() as u64 > max_bytes {
            return Err(too_large("body".to_string()));
        }
        buffer.extend_from_slice(&piece);
    }
    Ok(buffer)
}

/// HEAD `url`; returns `(content-length, etag)`.
pub(crate) async fn head(
    client: &reqwest::Client,
    url: &str,
    timeout: Duration,
) -> Result<(u64, Option<String>), Error> {
    let response = client
        .head(url)
        .timeout(timeout)
        .send()
        .await
        .map_err(|e| request_error(e, "HEAD failed"))?;
    let status = response.status();
    if classify(status.as_u16()).is_some() {
        return Err(status_error(status, url));
    }
    let header = |name: &str| {
        response
            .headers()
            .get(name)
            .and_then(|v| v.to_str().ok())
            .map(str::to_string)
    };
    let len = header("content-length")
        .and_then(|v| v.trim().parse().ok())
        .ok_or_else(|| {
            Error::new(
                ErrorKind::Transient,
                format!("HEAD {url} without Content-Length"),
            )
        })?;
    Ok((len, header("etag")))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn status_classes() {
        assert_eq!(classify(206), None);
        assert_eq!(classify(404), Some(ErrorKind::NotFound));
        assert_eq!(classify(410), Some(ErrorKind::NotFound));
        assert_eq!(classify(403), Some(ErrorKind::Permanent));
        assert_eq!(classify(416), Some(ErrorKind::Permanent));
        assert_eq!(classify(408), Some(ErrorKind::Transient));
        assert_eq!(classify(429), Some(ErrorKind::Throttled));
        assert_eq!(classify(503), Some(ErrorKind::Throttled));
        assert_eq!(classify(500), Some(ErrorKind::Transient));
    }
}
