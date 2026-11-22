//! `FastqOptions::threads` at full size: inputs spanning many 4 MiB pieces,
//! damage in the middle and near the end, and (with
//! `READKIT_FASTQ_REAL="r1.fastq.gz r2.fastq.gz"`) a real pair. Batches,
//! records and errors must equal the sequential reader's exactly. The
//! generated cases with small pieces are unit tests in `src/par.rs`.

use std::io::{BufRead, BufReader, Cursor};

use readkit_batch::{Batch, Segments};
use readkit_fastq::{FastqInput, FastqOptions, FastqSource};
use readkit_io::BatchSource;

/// Deterministic paired-style text: `records` records of 1..=250 bases.
fn text(records: usize, seed: u64) -> Vec<u8> {
    let mut x = seed;
    let mut next = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    let mut out = Vec::new();
    for i in 0..records {
        let len = 1 + (next() % 250) as usize;
        out.extend_from_slice(format!("@read{i} sample:{}\n", next() % 99).as_bytes());
        out.extend((0..len).map(|_| b"ACGT"[(next() & 3) as usize]));
        out.extend_from_slice(b"\n+\n");
        out.extend((0..len).map(|_| b'!' + (next() % 42) as u8));
        out.push(b'\n');
    }
    out
}

/// Records, batch sizes and the error, read with `threads`.
fn read_all(
    inputs: Vec<Box<dyn BufRead + Send>>,
    opts: FastqOptions,
) -> (Vec<Batch<Segments>>, Option<String>) {
    let inputs = inputs
        .into_iter()
        .enumerate()
        .map(|(i, r)| FastqInput::named(r, format!("r{}", i + 1)))
        .collect();
    let mut src = FastqSource::<Segments>::lockstep(inputs, opts).unwrap();
    let mut batches = Vec::new();
    loop {
        match src.next_batch() {
            Ok(Some(b)) => batches.push(b),
            Ok(None) => return (batches, None),
            Err(e) => return (batches, Some(e.to_string())),
        }
    }
}

fn same(a: &(Vec<Batch<Segments>>, Option<String>), b: &(Vec<Batch<Segments>>, Option<String>)) {
    assert_eq!(a.1, b.1, "errors");
    assert_eq!(a.0.len(), b.0.len(), "batch count");
    for (i, (x, y)) in a.0.iter().zip(&b.0).enumerate() {
        assert_eq!(x.len(), y.len(), "batch {i} size");
        for s in 0..x.segments() {
            let (p, q) = (x.segment(s), y.segment(s));
            assert_eq!(p.seqs().offsets(), q.seqs().offsets(), "batch {i}");
            assert!(p.seqs().data() == q.seqs().data(), "batch {i} bases");
            assert!(p.quals().data() == q.quals().data(), "batch {i} qualities");
            assert_eq!(
                p.names().map(|n| n.data().to_vec()),
                q.names().map(|n| n.data().to_vec()),
                "batch {i} names"
            );
        }
    }
}

fn slices(data: &[Vec<u8>], capacity: usize) -> Vec<Box<dyn BufRead + Send>> {
    data.iter()
        .map(|d| {
            Box::new(BufReader::with_capacity(capacity, Cursor::new(d.clone())))
                as Box<dyn BufRead + Send>
        })
        .collect()
}

/// Damage at `at` (a byte fraction) in input `r`.
fn damaged(mut data: Vec<Vec<u8>>, r: usize, at: f64, what: &str) -> Vec<Vec<u8>> {
    let d = &mut data[r];
    let mut pos = (d.len() as f64 * at) as usize;
    // The start of the record after `pos`.
    while d[pos] != b'@' || d[pos - 1] != b'\n' {
        pos += 1;
    }
    match what {
        "blank" => {
            d.insert(pos, b'\n');
        }
        "quality" => {
            let q = pos + d[pos..].iter().position(|&b| b == b'+').unwrap() + 2;
            d[q] = b' ';
        }
        "short" => {
            d.truncate(pos);
        }
        _ => unreachable!(),
    }
    data
}

#[test]
fn large_inputs_match_the_sequential_reader() {
    // About 18 MB per input: five pieces, a few dozen default batches.
    let clean = vec![text(120_000, 1), text(120_000, 2)];
    let cases = [
        ("clean", clean.clone()),
        (
            "blank line mid-file",
            damaged(clean.clone(), 1, 0.5, "blank"),
        ),
        (
            "bad quality near the end",
            damaged(clean.clone(), 0, 0.97, "quality"),
        ),
        (
            "bad quality mid-file",
            damaged(clean.clone(), 1, 0.4, "quality"),
        ),
        ("one input shorter", damaged(clean.clone(), 1, 0.8, "short")),
    ];
    for (name, data) in cases {
        for capacity in [1 << 20, 4096] {
            for budget in [1 << 20, 8 << 20] {
                let opts = FastqOptions::default().batch_bytes(budget);
                let want = read_all(slices(&data, capacity), opts.clone());
                let got = read_all(slices(&data, capacity), opts.clone().threads(3));
                same(&got, &want);
                eprintln!(
                    "{name}, buffer {capacity}, budget {budget}: {} batches, error {:?}",
                    want.0.len(),
                    want.1
                );
            }
        }
    }
}

#[test]
fn real_pair_matches_the_sequential_reader() {
    let Ok(paths) = std::env::var("READKIT_FASTQ_REAL") else {
        eprintln!("skipping: set READKIT_FASTQ_REAL=\"r1.fastq.gz r2.fastq.gz\"");
        return;
    };
    let paths: Vec<String> = paths.split_whitespace().map(str::to_owned).collect();
    let open = || -> Vec<Box<dyn BufRead + Send>> {
        paths
            .iter()
            .map(|p| {
                let r = readkit_gz::open(p, &Default::default()).unwrap();
                Box::new(BufReader::with_capacity(1 << 20, r)) as Box<dyn BufRead + Send>
            })
            .collect()
    };
    let want = read_all(open(), FastqOptions::default());
    let got = read_all(open(), FastqOptions::default().threads(8));
    same(&got, &want);
    let spots: usize = want.0.iter().map(|b| b.len()).sum();
    eprintln!("{spots} spots in {} batches identical", want.0.len());
}
