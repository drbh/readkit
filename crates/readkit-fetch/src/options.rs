//! Tuning knobs.

use std::time::Duration;

/// Connection cap for hosts matching a suffix.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HostLimit {
    /// Matches a host equal to this or ending in `.` + this (ports ignored).
    pub suffix: String,
    /// Most concurrent connections to a matching host (still capped by
    /// [`FetchOptions::max_connections`]).
    pub max_connections: usize,
}

impl HostLimit {
    /// A cap of `max_connections` for hosts matching `suffix`.
    pub fn new(suffix: impl Into<String>, max_connections: usize) -> Self {
        Self {
            suffix: suffix.into(),
            max_connections,
        }
    }

    fn matches(&self, host: &str) -> bool {
        let s = self.suffix.trim_start_matches('.');
        host.eq_ignore_ascii_case(s)
            || (host.len() > s.len()
                && host.as_bytes()[host.len() - s.len() - 1] == b'.'
                && host[host.len() - s.len()..].eq_ignore_ascii_case(s))
    }
}

/// Options for a [`Fetcher`](crate::Fetcher).
///
/// Construct with `FetchOptions { chunk_size: 4 << 20, ..Default::default() }`.
#[derive(Debug, Clone)]
pub struct FetchOptions {
    /// Largest transfer chunk in bytes (default 16 MiB). Smaller objects use
    /// smaller chunks, so that there is at least one per connection the host
    /// allows (and at least 16), with a 1 MiB floor. Every chunk is a request
    /// and pays a round trip, so larger chunks waste less time when there are
    /// enough of them to keep every connection busy.
    pub chunk_size: u64,
    /// Connections a host starts with before AIMD adapts it (default 32,
    /// clamped to the host's cap, so hosts start at their cap: 32 for S3, 8
    /// for others). AIMD halves it after throttling and grows it by 2 every
    /// 5 s while throughput still grows, so a low start leaves most downloads
    /// finished before it has grown.
    pub initial_connections: usize,
    /// Hard cap on connections to any one host (default 48).
    pub max_connections: usize,
    /// Cap for hosts matching no [`host_limits`](Self::host_limits) entry
    /// (default 8).
    pub default_host_connections: usize,
    /// Per-host caps, first match wins. The default allows 32 connections to
    /// `amazonaws.com` (S3); set to empty to treat every host alike.
    pub host_limits: Vec<HostLimit>,
    /// Fail a download after this long without receiving any byte (default
    /// 600 s). Time spent waiting on a slow stream consumer does not count.
    pub no_progress_timeout: Duration,
    /// Re-issue ("hedge") an in-flight chunk older than this once nothing else
    /// is queued (default 3 s), and older than three times the median time
    /// this download's chunks of its size have taken.
    pub hedge_after: Duration,
    /// Lower bound of a chunk request's deadline; the deadline grows with the
    /// chunk size over a quarter of the measured per-connection rate (default
    /// 10 s).
    pub min_chunk_timeout: Duration,
    /// First retry delay of a failing chunk; doubles per attempt (default
    /// 500 ms). A throttled connection also pauses this long.
    pub retry_delay: Duration,
    /// Longest retry delay (default 30 s).
    pub max_retry_delay: Duration,
    /// TCP/TLS connect timeout (default 20 s).
    pub connect_timeout: Duration,
    /// Parallel connections used to repair an object that failed verification
    /// (default 16).
    pub repair_connections: usize,
    /// `User-Agent` header (default `readkit-fetch/<version>`).
    pub user_agent: String,
}

impl Default for FetchOptions {
    fn default() -> Self {
        Self {
            chunk_size: 16 << 20,
            initial_connections: 32,
            max_connections: 48,
            default_host_connections: 8,
            host_limits: vec![HostLimit::new("amazonaws.com", 32)],
            no_progress_timeout: Duration::from_secs(600),
            hedge_after: Duration::from_secs(3),
            min_chunk_timeout: Duration::from_secs(10),
            retry_delay: Duration::from_millis(500),
            max_retry_delay: Duration::from_secs(30),
            connect_timeout: Duration::from_secs(20),
            repair_connections: 16,
            user_agent: concat!("readkit-fetch/", env!("CARGO_PKG_VERSION")).to_string(),
        }
    }
}

impl FetchOptions {
    /// Connection cap for `host` (no port).
    pub fn host_max_connections(&self, host: &str) -> usize {
        self.host_limits
            .iter()
            .find(|l| l.matches(host))
            .map_or(self.default_host_connections, |l| l.max_connections)
            .min(self.max_connections)
            .max(1)
    }
}

/// Options for a small whole-body GET ([`Fetcher::get_with`](crate::Fetcher::get_with)).
///
/// Construct with `GetOptions { max_bytes: 1 << 20, ..Default::default() }`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GetOptions {
    /// Largest body accepted, in bytes (default [`GetOptions::DEFAULT_MAX_BYTES`],
    /// 64 MiB). A larger `Content-Length`, or a body that grows past the
    /// limit, fails with [`ErrorKind::TooLarge`](crate::ErrorKind::TooLarge)
    /// without buffering more than the limit.
    pub max_bytes: u64,
    /// Deadline of one attempt, from connecting to the end of the body
    /// (default 120 s). Exceeding it is [`ErrorKind::Timeout`](crate::ErrorKind::Timeout).
    pub timeout: Duration,
    /// Further attempts after a retryable failure (throttling, 5xx, a dropped
    /// connection, a timeout; default 3). The delays follow
    /// [`FetchOptions::retry_delay`] (doubling per attempt), or the server's
    /// `Retry-After` seconds when longer, never more than
    /// [`FetchOptions::max_retry_delay`].
    pub retries: u32,
}

impl GetOptions {
    /// Default [`max_bytes`](Self::max_bytes): 64 MiB.
    pub const DEFAULT_MAX_BYTES: u64 = 64 << 20;
}

impl Default for GetOptions {
    fn default() -> Self {
        Self {
            max_bytes: Self::DEFAULT_MAX_BYTES,
            timeout: Duration::from_secs(120),
            retries: 3,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_host_classes() {
        let o = FetchOptions::default();
        assert_eq!(
            o.host_max_connections("sra-pub-run-odp.s3.amazonaws.com"),
            32
        );
        assert_eq!(o.host_max_connections("amazonaws.com"), 32);
        assert_eq!(o.host_max_connections("notamazonaws.com"), 8);
        assert_eq!(o.host_max_connections("ftp.sra.ebi.ac.uk"), 8);
    }

    #[test]
    fn host_limits_override() {
        let o = FetchOptions {
            host_limits: vec![HostLimit::new("ebi.ac.uk", 2)],
            max_connections: 16,
            default_host_connections: 64,
            ..Default::default()
        };
        assert_eq!(o.host_max_connections("ftp.sra.ebi.ac.uk"), 2);
        assert_eq!(o.host_max_connections("x.s3.amazonaws.com"), 16);
    }
}
