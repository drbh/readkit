//! Accession resolution (feature `resolve`): run, experiment, sample, study or project
//! accessions → runs → downloadable sources.
//!
//! Metadata comes from the ENA portal API; every run is also probed on
//! NCBI's public S3 bucket for its `.sra` object. Results are typed values;
//! choosing between the `.sra` object and ENA's FASTQ files is the
//! caller's policy.
//!
//! ```no_run
//! use readkit_sra::{FastqGroup, Resolver, SraStatus};
//!
//! let resolver = Resolver::new()?;
//! for run in resolver.resolve(&["PRJNA699558"])? {
//!     if let SraStatus::Available(object) = &run.sra {
//!         println!("{} .sra: {} bytes at {}", run.run, object.size, object.url);
//!     }
//!     if let Some(FastqGroup::Paired(r1, r2)) = run.fastq_group()? {
//!         println!("  FASTQ {} + {}", r1.name, r2.name);
//!     }
//! }
//! # Ok::<(), readkit_sra::Error>(())
//! ```

use std::collections::BTreeMap;
use std::sync::Mutex;
use std::time::Duration;

use readkit_fetch::{Checksum, FetchOptions, Fetcher, GetOptions, Object};

use crate::error::{Error, ErrorKind};
use crate::mates::{Mate, ena_mate};

/// ENA portal API base URL.
pub const ENA_PORTAL: &str = "https://www.ebi.ac.uk/ena/portal/api";
/// NCBI's public S3 bucket for `.sra` objects; a run is at `<base>/<run>/<run>`.
pub const NCBI_S3: &str = "https://sra-pub-run-odp.s3.amazonaws.com/sra";

/// ENA `read_run` fields requested.
pub const ENA_FIELDS: &str = "run_accession,study_accession,sample_accession,experiment_accession,library_layout,instrument_platform,instrument_model,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes,sample_title,library_strategy,library_source,scientific_name";

/// Run accessions per ENA `search` request.
const RUNS_PER_QUERY: usize = 200;
/// Concurrent S3 probes.
const PROBE_THREADS: usize = 8;

/// What an accession names.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum AccessionKind {
    /// `SRR`, `ERR`, `DRR`.
    Run,
    /// `SRX`, `ERX`, `DRX`.
    Experiment,
    /// `SRS`, `ERS`, `DRS`, `SAMN`, `SAMEA`, `SAMD`.
    Sample,
    /// `SRP`, `ERP`, `DRP`.
    Study,
    /// BioProject: `PRJNA`, `PRJEB`, `PRJDB`.
    Project,
}

impl AccessionKind {
    /// Classify an accession, or `None` if it is not an INSDC
    /// run/experiment/sample/study/project accession.
    ///
    /// ```
    /// use readkit_sra::AccessionKind;
    /// assert_eq!(AccessionKind::of("SRR13638690"), Some(AccessionKind::Run));
    /// assert_eq!(AccessionKind::of("PRJNA699558"), Some(AccessionKind::Project));
    /// assert_eq!(AccessionKind::of("SRR1x"), None);
    /// assert_eq!(AccessionKind::of("GSE12345"), None);
    /// ```
    pub fn of(accession: &str) -> Option<Self> {
        const PREFIXES: [(&str, AccessionKind); 17] = [
            ("SRR", AccessionKind::Run),
            ("ERR", AccessionKind::Run),
            ("DRR", AccessionKind::Run),
            ("SRX", AccessionKind::Experiment),
            ("ERX", AccessionKind::Experiment),
            ("DRX", AccessionKind::Experiment),
            ("SRS", AccessionKind::Sample),
            ("ERS", AccessionKind::Sample),
            ("DRS", AccessionKind::Sample),
            ("SAMN", AccessionKind::Sample),
            ("SAMEA", AccessionKind::Sample),
            ("SAMD", AccessionKind::Sample),
            ("SRP", AccessionKind::Study),
            ("ERP", AccessionKind::Study),
            ("DRP", AccessionKind::Study),
            ("PRJNA", AccessionKind::Project),
            ("PRJEB", AccessionKind::Project),
        ];
        let digits = |rest: &str| !rest.is_empty() && rest.bytes().all(|b| b.is_ascii_digit());
        if let Some(rest) = accession.strip_prefix("PRJDB") {
            return digits(rest).then_some(AccessionKind::Project);
        }
        PREFIXES
            .iter()
            .find_map(|(p, k)| accession.strip_prefix(p).filter(|r| digits(r)).map(|_| *k))
    }
}

/// Library layout as ENA records it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum LibraryLayout {
    /// `SINGLE`.
    Single,
    /// `PAIRED`.
    Paired,
}

/// Whether the run's `.sra` object is on NCBI's S3 bucket.
#[derive(Clone, Debug)]
#[non_exhaustive]
pub enum SraStatus {
    /// Present: URL, size and S3 ETag checksum, ready for
    /// [`readkit_fetch::Fetcher::fetch`].
    Available(Object),
    /// Not on the bucket (HTTP 404 or 403).
    Absent,
    /// The probe failed for another reason (network, throttling).
    Unknown(Error),
}

impl SraStatus {
    /// The object, if available.
    pub fn object(&self) -> Option<&Object> {
        match self {
            SraStatus::Available(o) => Some(o),
            _ => None,
        }
    }
}

/// One FASTQ file from ENA.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FastqFile {
    /// The file's basename (`SRR1_1.fastq.gz`).
    pub name: String,
    /// Which mate it holds.
    pub mate: Mate,
    /// HTTPS URL, size and MD5.
    pub object: Object,
}

/// A run's FASTQ files arranged by mate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum FastqGroup<'a> {
    /// One single-end file.
    Single(&'a FastqFile),
    /// R1 and R2.
    Paired(&'a FastqFile, &'a FastqFile),
}

/// A resolved run.
#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct RunRecord {
    /// Run accession.
    pub run: String,
    /// Study accession.
    pub study: String,
    /// Sample accession.
    pub sample: String,
    /// Experiment accession.
    pub experiment: String,
    /// Library layout.
    pub layout: LibraryLayout,
    /// Instrument platform (`ILLUMINA`, ...).
    pub platform: String,
    /// Instrument model.
    pub model: String,
    /// Spots, as ENA reports them.
    pub read_count: Option<u64>,
    /// Bases, as ENA reports them.
    pub base_count: Option<u64>,
    /// The `.sra` object on NCBI S3.
    pub sra: SraStatus,
    /// ENA's FASTQ files, sorted by mate (single, R1, R2, unpaired).
    pub fastq: Vec<FastqFile>,
    /// Every non-empty ENA field, as returned.
    pub fields: BTreeMap<String, String>,
}

impl RunRecord {
    /// The FASTQ files as a complete group for the layout: one
    /// [`Mate::Single`] file for single-end runs, one R1 and one R2 for
    /// paired runs ([`Mate::Unpaired`] files are ignored). `None` when ENA
    /// lists no files; [`ErrorKind::Metadata`] when they do not form a
    /// group.
    pub fn fastq_group(&self) -> Result<Option<FastqGroup<'_>>, Error> {
        if self.fastq.is_empty() {
            return Ok(None);
        }
        let of = |m: Mate| self.fastq.iter().filter(move |f| f.mate == m);
        let one = |m: Mate| {
            let mut it = of(m);
            match (it.next(), it.next()) {
                (Some(f), None) => Some(f),
                _ => None,
            }
        };
        let group = match self.layout {
            LibraryLayout::Single if self.fastq.len() == 1 => {
                one(Mate::Single).map(FastqGroup::Single)
            }
            LibraryLayout::Single => None,
            LibraryLayout::Paired if of(Mate::Single).next().is_none() => one(Mate::R1)
                .zip(one(Mate::R2))
                .map(|(a, b)| FastqGroup::Paired(a, b)),
            LibraryLayout::Paired => None,
        };
        group.map(Some).ok_or_else(|| {
            let mates: Vec<Mate> = self.fastq.iter().map(|f| f.mate).collect();
            Error::new(
                ErrorKind::Metadata,
                format!(
                    "ENA FASTQ files do not form a complete {:?} group: {mates:?}",
                    self.layout
                ),
            )
            .with_input(self.run.clone())
        })
    }
}

/// How the resolver reaches the network. Implement it to resolve against
/// recorded responses, a proxy or a mirror.
pub trait Transport: Send + Sync {
    /// GET `url` and return the body as text. An empty body is fine.
    fn get_text(&self, url: &str) -> Result<String, Error>;

    /// Describe the S3 object at `url` (size and ETag checksum). A missing
    /// object is an error of kind [`ErrorKind::NotFound`].
    fn probe_s3(&self, url: &str) -> Result<Object, Error>;
}

/// The default [`Transport`]: [`readkit_fetch::Fetcher::get_with`] for ENA
/// queries and [`readkit_fetch::Fetcher::probe_s3`] for S3.
///
/// ENA GETs are retried up to 3 times (after 2, 4 and 8 s, or the server's
/// `Retry-After` when longer) on network errors, timeouts, HTTP 429 and 5xx;
/// each attempt has 120 s and a body may be up to 64 MiB. Safe to call from
/// any thread, including inside an async runtime (the call blocks that
/// thread).
pub struct HttpTransport {
    fetcher: Fetcher,
    get: GetOptions,
}

impl HttpTransport {
    /// A transport over its own [`Fetcher`] (a small two-thread I/O runtime).
    pub fn new() -> Result<Self, Error> {
        let fetcher = Fetcher::new(FetchOptions {
            retry_delay: Duration::from_secs(2),
            ..FetchOptions::default()
        })?;
        Ok(Self::with_fetcher(fetcher))
    }

    /// A transport over `fetcher`, sharing its connections and options
    /// (user agent, retry delays). ENA GETs use [`GetOptions::default`].
    pub fn with_fetcher(fetcher: Fetcher) -> Self {
        Self {
            fetcher,
            get: GetOptions::default(),
        }
    }

    /// Use `options` (size limit, deadline, retries) for ENA GETs.
    pub fn with_get_options(mut self, options: GetOptions) -> Self {
        self.get = options;
        self
    }
}

impl Transport for HttpTransport {
    fn get_text(&self, url: &str) -> Result<String, Error> {
        let body = self.fetcher.get_with(url, &self.get)?;
        Ok(match String::from_utf8(body) {
            Ok(text) => text,
            Err(e) => String::from_utf8_lossy(e.as_bytes()).into_owned(),
        })
    }

    fn probe_s3(&self, url: &str) -> Result<Object, Error> {
        Ok(self.fetcher.probe_s3(url)?)
    }
}

impl std::fmt::Debug for HttpTransport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HttpTransport")
            .field("get", &self.get)
            .finish_non_exhaustive()
    }
}

/// Resolves accessions to [`RunRecord`]s.
pub struct Resolver {
    transport: Box<dyn Transport>,
    ena: String,
    s3: String,
}

impl Resolver {
    /// A resolver over [`HttpTransport`] and the public endpoints.
    pub fn new() -> Result<Self, Error> {
        Ok(Self::with_transport(Box::new(HttpTransport::new()?)))
    }

    /// A resolver over `transport` and the public endpoints.
    pub fn with_transport(transport: Box<dyn Transport>) -> Self {
        Self {
            transport,
            ena: ENA_PORTAL.to_string(),
            s3: NCBI_S3.to_string(),
        }
    }

    /// Use other ENA portal and S3 base URLs (no trailing slash).
    pub fn with_endpoints(mut self, ena: impl Into<String>, s3: impl Into<String>) -> Self {
        self.ena = ena.into();
        self.s3 = s3.into();
        self
    }

    /// The ENA `search` URL for a batch of run accessions.
    pub fn search_url(&self, runs: &[&str]) -> String {
        format!(
            "{}/search?result=read_run&includeAccessions={}&fields={ENA_FIELDS}&format=json&limit=0",
            self.ena,
            runs.join(",")
        )
    }

    /// The ENA `filereport` URL for an experiment, sample, study or project.
    pub fn filereport_url(&self, accession: &str) -> String {
        format!(
            "{}/filereport?accession={accession}&result=read_run&fields={ENA_FIELDS}&format=json",
            self.ena
        )
    }

    /// The S3 URL of a run's `.sra` object.
    pub fn s3_url(&self, run: &str) -> String {
        format!("{}/{run}/{run}", self.s3)
    }

    /// Resolve accessions of any supported kind to their runs, sorted by
    /// run accession and without duplicates.
    ///
    /// Runs are looked up in batched ENA `search` requests (200 per
    /// request); every other accession uses one `filereport` request. Every
    /// run is then probed on S3. Errors: [`ErrorKind::Unsupported`] for an
    /// accession of another kind, [`ErrorKind::NotFound`] when ENA knows a
    /// run or group not, [`ErrorKind::Metadata`] for malformed rows,
    /// [`ErrorKind::Network`] when ENA cannot be reached. A failed S3 probe
    /// does not fail the call ([`SraStatus::Unknown`]).
    pub fn resolve(&self, accessions: &[&str]) -> Result<Vec<RunRecord>, Error> {
        let mut runs = Vec::new();
        let mut groups = Vec::new();
        for &a in accessions {
            match AccessionKind::of(a) {
                Some(AccessionKind::Run) => runs.push(a),
                Some(_) => groups.push(a),
                None => {
                    return Err(Error::new(
                        ErrorKind::Unsupported,
                        "not a run, experiment, sample, study or project accession",
                    )
                    .with_input(a));
                }
            }
        }
        let mut rows = Vec::new();
        for chunk in runs.chunks(RUNS_PER_QUERY) {
            let url = self.search_url(chunk);
            rows.extend(parse_ena_rows(&self.transport.get_text(&url)?)?);
        }
        for group in &groups {
            let found = parse_ena_rows(&self.transport.get_text(&self.filereport_url(group))?)?;
            if found.is_empty() {
                return Err(
                    Error::new(ErrorKind::NotFound, "ENA returned no runs").with_input(*group)
                );
            }
            rows.extend(found);
        }
        let mut by_run = BTreeMap::new();
        for row in rows {
            let run = row.get("run_accession").cloned().unwrap_or_default();
            if run.is_empty() {
                return Err(Error::new(
                    ErrorKind::Metadata,
                    "ENA row without run_accession",
                ));
            }
            by_run.insert(run, row);
        }
        for run in &runs {
            if !by_run.contains_key(*run) {
                return Err(
                    Error::new(ErrorKind::NotFound, "ENA has no record of this run")
                        .with_input(*run),
                );
            }
        }
        let names: Vec<&String> = by_run.keys().collect();
        let probes = self.probe_all(&names);
        by_run
            .iter()
            .zip(probes)
            .map(|((run, row), sra)| run_from_row(run, row, sra))
            .collect()
    }

    fn probe_all(&self, runs: &[&String]) -> Vec<SraStatus> {
        let out: Mutex<Vec<Option<SraStatus>>> = Mutex::new(vec![None; runs.len()]);
        let next = std::sync::atomic::AtomicUsize::new(0);
        std::thread::scope(|scope| {
            for _ in 0..PROBE_THREADS.min(runs.len()) {
                scope.spawn(|| {
                    loop {
                        let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        let Some(run) = runs.get(i) else { return };
                        let status = probe_status(self.transport.probe_s3(&self.s3_url(run)));
                        out.lock().unwrap_or_else(|e| e.into_inner())[i] = Some(status);
                    }
                });
            }
        });
        out.into_inner()
            .unwrap_or_else(|e| e.into_inner())
            .into_iter()
            .map(|s| {
                s.unwrap_or_else(|| {
                    SraStatus::Unknown(Error::new(ErrorKind::Internal, "probe thread panicked"))
                })
            })
            .collect()
    }
}

impl std::fmt::Debug for Resolver {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Resolver")
            .field("ena", &self.ena)
            .field("s3", &self.s3)
            .finish()
    }
}

fn probe_status(result: Result<Object, Error>) -> SraStatus {
    match result {
        Ok(o) => SraStatus::Available(o),
        // `From<readkit_fetch::Error>` keeps the HTTP status as the code.
        Err(e) if e.kind() == ErrorKind::NotFound || e.code() == Some(403) => SraStatus::Absent,
        Err(e) => SraStatus::Unknown(e),
    }
}

/// Parse an ENA portal JSON response (an array of objects) into rows of
/// field → value. Scalar values are kept as text; nulls are dropped. An
/// empty or whitespace-only body is no rows.
pub fn parse_ena_rows(text: &str) -> Result<Vec<BTreeMap<String, String>>, Error> {
    if text.trim().is_empty() {
        return Ok(Vec::new());
    }
    let value: serde_json::Value = serde_json::from_str(text)
        .map_err(|e| Error::new(ErrorKind::Metadata, "ENA response is not JSON").with_source(e))?;
    let serde_json::Value::Array(items) = value else {
        return Err(Error::new(
            ErrorKind::Metadata,
            "ENA response is not a JSON array",
        ));
    };
    items
        .into_iter()
        .map(|item| {
            let serde_json::Value::Object(map) = item else {
                return Err(Error::new(
                    ErrorKind::Metadata,
                    "ENA row is not a JSON object",
                ));
            };
            Ok(map
                .into_iter()
                .filter_map(|(k, v)| {
                    let v = match v {
                        serde_json::Value::String(s) => s,
                        serde_json::Value::Null => return None,
                        other => other.to_string(),
                    };
                    Some((k, v))
                })
                .collect())
        })
        .collect()
}

fn split(value: Option<&String>) -> Vec<&str> {
    value
        .map(|v| v.split(';').filter(|s| !s.is_empty()).collect())
        .unwrap_or_default()
}

/// Build a [`RunRecord`] from one ENA row and the run's S3 probe result.
pub fn run_from_row(
    run: &str,
    row: &BTreeMap<String, String>,
    sra: SraStatus,
) -> Result<RunRecord, Error> {
    let get = |k: &str| row.get(k).cloned().unwrap_or_default();
    let meta = |msg: String| Error::new(ErrorKind::Metadata, msg).with_input(run.to_string());
    let layout = match get("library_layout").as_str() {
        "PAIRED" => LibraryLayout::Paired,
        "SINGLE" => LibraryLayout::Single,
        other => return Err(meta(format!("unsupported library layout {other:?}"))),
    };
    let urls = split(row.get("fastq_ftp"));
    let md5s = split(row.get("fastq_md5"));
    let sizes = split(row.get("fastq_bytes"));
    if md5s.len() != urls.len() || sizes.len() != urls.len() {
        return Err(meta(format!(
            "{} FASTQ URLs, {} MD5s, {} sizes",
            urls.len(),
            md5s.len(),
            sizes.len()
        )));
    }
    let mut fastq = Vec::with_capacity(urls.len());
    for ((url, md5), size) in urls.iter().zip(&md5s).zip(&sizes) {
        let name = url.rsplit('/').next().unwrap_or(url).to_string();
        let mate = ena_mate(run, &name, layout == LibraryLayout::Paired)
            .map_err(|e| e.with_input(run.to_string()))?;
        let size: u64 = size
            .parse()
            .map_err(|_| meta(format!("bad fastq_bytes {size:?}")))?;
        if md5.len() != 32 || !md5.bytes().all(|b| b.is_ascii_hexdigit()) {
            return Err(meta(format!("bad fastq_md5 {md5:?}")));
        }
        let rest = url
            .strip_prefix("ftp://")
            .or_else(|| url.strip_prefix("https://"))
            .or_else(|| url.strip_prefix("http://"))
            .unwrap_or(url);
        fastq.push(FastqFile {
            name,
            mate,
            object: Object::new(
                format!("https://{rest}"),
                size,
                Checksum::Md5 {
                    hex: md5.to_ascii_lowercase(),
                },
            ),
        });
    }
    fastq.sort_by_key(|f| f.mate);
    Ok(RunRecord {
        run: run.to_string(),
        study: get("study_accession"),
        sample: get("sample_accession"),
        experiment: get("experiment_accession"),
        layout,
        platform: get("instrument_platform"),
        model: get("instrument_model"),
        read_count: row.get("read_count").and_then(|v| v.parse().ok()),
        base_count: row.get("base_count").and_then(|v| v.parse().ok()),
        sra,
        fastq,
        fields: row
            .iter()
            .filter(|(_, v)| !v.is_empty())
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accession_kinds() {
        use AccessionKind::*;
        for (a, k) in [
            ("SRR13638690", Some(Run)),
            ("ERR1000001", Some(Run)),
            ("DRR2222222", Some(Run)),
            ("SRX9999", Some(Experiment)),
            ("SRS1", Some(Sample)),
            ("SAMN12345", Some(Sample)),
            ("SRP301234", Some(Study)),
            ("ERP1", Some(Study)),
            ("PRJNA699558", Some(Project)),
            ("PRJEB1", Some(Project)),
            ("PRJDB7", Some(Project)),
            ("SRR", None),
            ("SRR1x", None),
            ("srr1", None),
            ("GSE1", None),
            ("PRJNA", None),
        ] {
            assert_eq!(AccessionKind::of(a), k, "{a}");
        }
    }

    #[test]
    fn urls() {
        let r = Resolver::with_transport(Box::new(Never)).with_endpoints("http://ena", "http://s3");
        assert_eq!(r.s3_url("SRR1"), "http://s3/SRR1/SRR1");
        assert!(r.search_url(&["SRR1", "SRR2"]).starts_with(
            "http://ena/search?result=read_run&includeAccessions=SRR1,SRR2&fields=run_accession,"
        ));
        assert!(
            r.filereport_url("PRJNA1")
                .starts_with("http://ena/filereport?accession=PRJNA1&")
        );
    }

    struct Never;
    impl Transport for Never {
        fn get_text(&self, _: &str) -> Result<String, Error> {
            unreachable!()
        }
        fn probe_s3(&self, _: &str) -> Result<Object, Error> {
            unreachable!()
        }
    }

    #[test]
    fn rows_tolerate_numbers_nulls_and_empty_bodies() {
        assert!(parse_ena_rows("").unwrap().is_empty());
        assert!(parse_ena_rows(" \n").unwrap().is_empty());
        assert!(parse_ena_rows("[]").unwrap().is_empty());
        let rows =
            parse_ena_rows(r#"[{"run_accession":"SRR1","read_count":12,"x":null}]"#).unwrap();
        assert_eq!(rows[0]["read_count"], "12");
        assert!(!rows[0].contains_key("x"));
        for bad in ["{", "{}", "[1]"] {
            assert_eq!(
                parse_ena_rows(bad).unwrap_err().kind(),
                ErrorKind::Metadata,
                "{bad}"
            );
        }
    }

    fn row(pairs: &[(&str, &str)]) -> BTreeMap<String, String> {
        pairs
            .iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect()
    }

    #[test]
    fn row_with_mismatched_lists_is_rejected() {
        let r = row(&[
            ("library_layout", "PAIRED"),
            ("fastq_ftp", "a/SRR1_1.fastq.gz;a/SRR1_2.fastq.gz"),
            ("fastq_md5", "0123456789abcdef0123456789abcdef"),
            ("fastq_bytes", "1;2"),
        ]);
        let e = run_from_row("SRR1", &r, SraStatus::Absent).unwrap_err();
        assert_eq!(e.kind(), ErrorKind::Metadata);
    }

    #[test]
    fn unknown_layout_is_rejected() {
        let r = row(&[("library_layout", "")]);
        assert_eq!(
            run_from_row("SRR1", &r, SraStatus::Absent)
                .unwrap_err()
                .kind(),
            ErrorKind::Metadata
        );
    }

    #[test]
    fn probe_outcomes() {
        assert!(matches!(
            probe_status(Err(Error::new(ErrorKind::NotFound, ""))),
            SraStatus::Absent
        ));
        assert!(matches!(
            probe_status(Err(Error::new(ErrorKind::Network, "").with_code(403))),
            SraStatus::Absent
        ));
        assert!(matches!(
            probe_status(Err(Error::new(ErrorKind::Network, "reset"))),
            SraStatus::Unknown(_)
        ));
    }
}
