//! Resolver against recorded ENA and S3 responses (offline), the real
//! `HttpTransport` against a local server (offline), plus one `#[ignore]`d
//! test against the live services. Needs feature `resolve`.

use std::collections::HashMap;
use std::sync::Mutex;

use readkit_fetch::{Checksum, Object};
use readkit_sra::{
    Error, ErrorKind, FastqGroup, HttpTransport, LibraryLayout, Mate, Resolver, SraStatus,
    Transport,
};

const SEARCH: &str = include_str!("fixtures/ena_search_runs.json");
const PROJECT: &str = include_str!("fixtures/ena_filereport_project.json");
const HEADS: &str = include_str!("fixtures/s3_heads.txt");

/// Serves the fixtures and records every URL requested.
struct Recorded {
    seen: Mutex<Vec<String>>,
    heads: HashMap<String, Result<Object, Error>>,
}

impl Recorded {
    fn new() -> Self {
        let mut heads = HashMap::new();
        for line in HEADS
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split_whitespace().collect();
            let url = format!("https://sra-pub-run-odp.s3.amazonaws.com/sra/{0}/{0}", f[0]);
            let result = if f[1] == "200" {
                let size: u64 = f[2].parse().unwrap();
                let part = f[4].parse().unwrap_or(size);
                // What readkit-fetch's probe_s3 builds from the same headers.
                let mut checksum = Checksum::from_s3_etag(f[3], size).unwrap();
                if let Checksum::S3MultipartEtag { part_size, .. } = &mut checksum {
                    *part_size = part;
                }
                Ok(Object::new(url.clone(), size, checksum))
            } else {
                Err(Error::new(ErrorKind::NotFound, "HTTP 404"))
            };
            heads.insert(url, result);
        }
        Self {
            seen: Mutex::new(Vec::new()),
            heads,
        }
    }
}

impl Transport for Recorded {
    fn get_text(&self, url: &str) -> Result<String, Error> {
        self.seen.lock().unwrap().push(url.to_string());
        if url.contains("/search?") {
            Ok(SEARCH.to_string())
        } else if url.contains("accession=PRJNA699558") {
            Ok(PROJECT.to_string())
        } else {
            Ok(String::new())
        }
    }

    fn probe_s3(&self, url: &str) -> Result<Object, Error> {
        self.heads
            .get(url)
            .cloned()
            .unwrap_or_else(|| Err(Error::new(ErrorKind::Network, "connection reset")))
    }
}

#[test]
fn runs_resolve_to_s3_objects_and_ena_files() {
    let resolver = Resolver::with_transport(Box::new(Recorded::new()));
    let runs = resolver.resolve(&["SRR13638690", "SRR12999707"]).unwrap();
    assert_eq!(
        runs.iter().map(|r| r.run.as_str()).collect::<Vec<_>>(),
        ["SRR12999707", "SRR13638690"]
    );

    let yeast = &runs[0];
    assert_eq!(yeast.layout, LibraryLayout::Single);
    assert_eq!(yeast.read_count, Some(89234));
    assert_eq!(yeast.experiment, "SRX9479621");
    assert_eq!(yeast.fields["scientific_name"], "Saccharomyces cerevisiae");
    let sra = yeast.sra.object().unwrap();
    assert_eq!(sra.size, 16_143_261);
    assert!(matches!(
        &sra.checksum,
        Checksum::S3MultipartEtag { parts: 1, .. }
    ));
    let Some(FastqGroup::Single(f)) = yeast.fastq_group().unwrap() else {
        panic!("single-end group")
    };
    assert_eq!(f.mate, Mate::Single);
    assert_eq!(
        f.object.url,
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR129/007/SRR12999707/SRR12999707.fastq.gz"
    );
    assert_eq!(
        f.object.checksum,
        Checksum::Md5 {
            hex: "5d41402abc4b2a76b9719d911017c592".into()
        },
        "MD5 normalized to lowercase"
    );

    let human = &runs[1];
    assert_eq!(human.layout, LibraryLayout::Paired);
    match &human.sra.object().unwrap().checksum {
        Checksum::S3MultipartEtag {
            parts, part_size, ..
        } => assert_eq!((*parts, *part_size), (85, 8_388_608)),
        other => panic!("{other:?}"),
    }
    let Some(FastqGroup::Paired(r1, r2)) = human.fastq_group().unwrap() else {
        panic!("paired group")
    };
    assert_eq!((r1.mate, r2.mate), (Mate::R1, Mate::R2));
    assert_eq!((r1.object.size, r2.object.size), (401_234_567, 412_345_678));
}

#[test]
fn projects_expand_to_runs_with_mates_classified() {
    let transport = Recorded::new();
    let resolver = Resolver::with_transport(Box::new(transport));
    let runs = resolver.resolve(&["PRJNA699558"]).unwrap();
    assert_eq!(runs.len(), 3);

    // SRR13638686: no ENA FASTQ, no S3 object.
    assert!(runs[0].fastq.is_empty());
    assert!(runs[0].fastq_group().unwrap().is_none());
    assert!(matches!(runs[0].sra, SraStatus::Absent));

    // SRR13638687: an orphan file next to the pair is kept as Unpaired and
    // does not break the group.
    let mates: Vec<Mate> = runs[1].fastq.iter().map(|f| f.mate).collect();
    assert_eq!(mates, [Mate::R1, Mate::R2, Mate::Unpaired]);
    assert!(matches!(
        runs[1].fastq_group().unwrap(),
        Some(FastqGroup::Paired(_, _))
    ));

    // SRR13638690: files listed R2 first are sorted; `SRR1…` digits are
    // never read as a mate token.
    let names: Vec<&str> = runs[2].fastq.iter().map(|f| f.name.as_str()).collect();
    assert_eq!(names, ["SRR13638690_1.fastq.gz", "SRR13638690_2.fastq.gz"]);
}

#[test]
fn runs_unknown_to_ena_are_not_found() {
    let resolver = Resolver::with_transport(Box::new(Recorded::new()));
    let e = resolver
        .resolve(&["SRR12999707", "SRR20000000"])
        .unwrap_err();
    assert_eq!(
        (e.kind(), e.input()),
        (ErrorKind::NotFound, Some("SRR20000000"))
    );
}

#[test]
fn search_urls_split_at_200_runs() {
    struct Count(Mutex<Vec<usize>>);
    impl Transport for Count {
        fn get_text(&self, url: &str) -> Result<String, Error> {
            let list = url.split("includeAccessions=").nth(1).unwrap();
            let list = list.split('&').next().unwrap();
            self.0.lock().unwrap().push(list.split(',').count());
            Ok("[]".into())
        }
        fn probe_s3(&self, _: &str) -> Result<Object, Error> {
            unreachable!()
        }
    }
    let counts = std::sync::Arc::new(Count(Mutex::new(Vec::new())));
    struct Shared(std::sync::Arc<Count>);
    impl Transport for Shared {
        fn get_text(&self, url: &str) -> Result<String, Error> {
            self.0.get_text(url)
        }
        fn probe_s3(&self, url: &str) -> Result<Object, Error> {
            self.0.probe_s3(url)
        }
    }
    let resolver = Resolver::with_transport(Box::new(Shared(counts.clone())));
    let many: Vec<String> = (0..450).map(|i| format!("SRR{}", 1_000 + i)).collect();
    let refs: Vec<&str> = many.iter().map(String::as_str).collect();
    let _ = resolver.resolve(&refs);
    assert_eq!(*counts.0.lock().unwrap(), [200, 200, 50]);
}

#[test]
fn bad_accessions_and_empty_groups() {
    let resolver = Resolver::with_transport(Box::new(Recorded::new()));
    let e = resolver.resolve(&["GSE12345"]).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Unsupported);
    let e = resolver.resolve(&["PRJNA1"]).unwrap_err();
    assert_eq!((e.kind(), e.input()), (ErrorKind::NotFound, Some("PRJNA1")));
}

#[test]
fn failed_probes_do_not_fail_resolution() {
    struct NoS3;
    impl Transport for NoS3 {
        fn get_text(&self, _: &str) -> Result<String, Error> {
            Ok(SEARCH.to_string())
        }
        fn probe_s3(&self, _: &str) -> Result<Object, Error> {
            Err(Error::new(ErrorKind::Network, "throttled"))
        }
    }
    let runs = Resolver::with_transport(Box::new(NoS3))
        .resolve(&["SRR12999707"])
        .unwrap();
    assert!(runs.iter().all(|r| matches!(r.sra, SraStatus::Unknown(_))));
}

/// A one-request-per-connection HTTP server: ENA `search` returns the
/// fixture, the yeast run's S3 HEAD succeeds, the human run's is 403, and
/// anything else is 404.
fn serve_ena_and_s3() -> String {
    use std::io::{BufRead, BufReader, Write};
    let listener = std::net::TcpListener::bind("127.0.0.1:0").unwrap();
    let base = format!("http://{}", listener.local_addr().unwrap());
    std::thread::spawn(move || {
        for stream in listener.incoming().flatten() {
            std::thread::spawn(move || {
                let mut reader = BufReader::new(stream.try_clone().unwrap());
                let mut request = String::new();
                reader.read_line(&mut request).unwrap();
                let mut line = String::new();
                while reader.read_line(&mut line).unwrap_or(0) > 2 {
                    line.clear();
                }
                let head = |status: &str, extra: &str, len: usize| {
                    format!(
                        "HTTP/1.1 {status}\r\n{extra}Content-Length: {len}\r\nConnection: close\r\n\r\n"
                    )
                };
                let response = if request.starts_with("GET /ena/search?") {
                    head("200 OK", "Content-Type: application/json\r\n", SEARCH.len()) + SEARCH
                } else if request.starts_with("HEAD /s3/SRR12999707/SRR12999707 ") {
                    head(
                        "200 OK",
                        "ETag: \"1f3870be274f6c49b3e31a0c6728957f\"\r\n",
                        16_143_261,
                    )
                } else if request.starts_with("HEAD /s3/SRR13638690/") {
                    head("403 Forbidden", "", 0)
                } else {
                    head("404 Not Found", "", 0)
                };
                let mut out = stream;
                let _ = out.write_all(response.as_bytes());
            });
        }
    });
    base
}

#[test]
fn http_transport_against_a_local_server() {
    let base = serve_ena_and_s3();
    let fetcher = readkit_fetch::Fetcher::new(readkit_fetch::FetchOptions::default()).unwrap();
    let resolver = Resolver::with_transport(Box::new(HttpTransport::with_fetcher(fetcher.clone())))
        .with_endpoints(format!("{base}/ena"), format!("{base}/s3"));
    let runs = resolver.resolve(&["SRR13638690", "SRR12999707"]).unwrap();
    assert_eq!(runs[0].run, "SRR12999707");
    let sra = runs[0].sra.object().expect("HEAD 200");
    assert_eq!(sra.size, 16_143_261);
    assert!(matches!(
        &sra.checksum,
        Checksum::S3MultipartEtag { parts: 1, .. }
    ));
    assert!(matches!(runs[1].sra, SraStatus::Absent), "HEAD 403");

    // An ENA 404 is a typed NotFound carrying the status; a body over the
    // limit is a Network error whose source is readkit-fetch's TooLarge.
    let transport = HttpTransport::with_fetcher(fetcher);
    let e = transport
        .get_text(&format!("{base}/ena/nothing"))
        .unwrap_err();
    assert_eq!((e.kind(), e.code()), (ErrorKind::NotFound, Some(404)));
    let small = HttpTransport::new()
        .unwrap()
        .with_get_options(readkit_fetch::GetOptions {
            max_bytes: 16,
            ..Default::default()
        });
    let e = small.get_text(&format!("{base}/ena/search?x")).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Network);
    let cause = std::error::Error::source(&e)
        .and_then(|s| s.downcast_ref::<readkit_fetch::Error>())
        .expect("readkit-fetch cause");
    assert_eq!(cause.kind(), readkit_fetch::ErrorKind::TooLarge);
}

/// Live ENA + NCBI S3. Run with `cargo test -p readkit-sra -- --ignored`.
#[test]
#[ignore = "needs the network (ENA portal and NCBI S3)"]
fn live_resolution() {
    let resolver = Resolver::new().unwrap();
    let runs = resolver.resolve(&["SRR12999707", "SRR13638690"]).unwrap();
    assert_eq!(runs.len(), 2);
    let yeast = &runs[0];
    assert_eq!(yeast.sra.object().map(|o| o.size), Some(16_143_261));
    let human = &runs[1];
    assert_eq!(human.layout, LibraryLayout::Paired);
    assert_eq!(human.read_count, Some(4_451_399));
    assert_eq!(human.sra.object().map(|o| o.size), Some(713_058_679));
    assert!(matches!(
        human.fastq_group().unwrap(),
        Some(FastqGroup::Paired(_, _))
    ));
    let project = resolver.resolve(&["PRJNA699558"]).unwrap();
    assert!(project.iter().any(|r| r.run == "SRR13638690"));
}
