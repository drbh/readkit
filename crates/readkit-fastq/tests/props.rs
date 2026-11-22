//! Fuzz-style property tests with a deterministic generator: random record
//! shapes and random corruptions, checked against fasterp's reference
//! reader, and arbitrary bytes, which must never panic.

mod support;

use readkit_batch::{Paired, Segments, Single};
use readkit_fastq::{FastqOptions, FastqSource};
use support::*;

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        // xorshift64*
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
    fn chance(&mut self, percent: usize) -> bool {
        self.below(100) < percent
    }
}

/// A random, mostly valid FASTQ text of `n` records.
fn random_fastq(rng: &mut Rng, n: usize, phred64: bool) -> Vec<u8> {
    let mut out = Vec::new();
    let lo = if phred64 { 64 } else { 33 };
    for i in 0..n {
        if rng.chance(10) {
            let blanks = [&b"\n"[..], b"  \n", b"\t\r\n", b"\x0b\x0c\n", b"\r\n"];
            out.extend_from_slice(blanks[rng.below(blanks.len())]);
        }
        let eol: &[u8] = if rng.chance(20) { b"\r\n" } else { b"\n" };
        let name = format!("r{i}{}", if rng.chance(30) { " desc=1" } else { "" });
        let len = if rng.chance(10) { 0 } else { rng.below(300) };
        out.push(b'@');
        out.extend_from_slice(name.as_bytes());
        out.extend_from_slice(eol);
        for _ in 0..len {
            out.push(b"ACGTNacgtnRYKM"[rng.below(14)]);
        }
        out.extend_from_slice(eol);
        out.push(b'+');
        if rng.chance(30) {
            out.extend_from_slice(name.as_bytes());
        }
        out.extend_from_slice(eol);
        for _ in 0..len {
            out.push(lo + rng.below(127 - lo as usize) as u8);
        }
        if i + 1 < n || len == 0 || rng.chance(50) {
            out.extend_from_slice(eol);
        }
    }
    if rng.chance(10) {
        out.extend_from_slice(b"\n \n");
    }
    out
}

/// Flips, inserts or deletes a few bytes, biased toward FASTQ syntax.
fn corrupt(rng: &mut Rng, data: &mut Vec<u8>) {
    for _ in 0..1 + rng.below(3) {
        if data.is_empty() {
            data.push(b'@');
            continue;
        }
        let at = rng.below(data.len());
        let special = [
            b'\n', b'\r', b'@', b'+', b' ', b'\t', 0, 0xff, b'!', b'~', 0x7f, 31,
        ];
        let byte = special[rng.below(special.len())];
        match rng.below(4) {
            0 => data[at] = byte,
            1 => data.insert(at, byte),
            2 => {
                data.remove(at);
            }
            _ => data.truncate(at),
        }
    }
}

/// The reference's records (spot order) and error.
fn reference_run(
    inputs: &[&[u8]],
    interleaved: bool,
    phred64: bool,
    cap: usize,
) -> (Vec<Triple>, Option<String>) {
    let mut first = reference(inputs[0], cap, phred64);
    let mut second = inputs.get(1).map(|d| reference(d, cap + 1, phred64));
    let mut out = Vec::new();
    let triple = |r: &fasterp::Read| (r.name[1..].to_vec(), r.seq.to_vec(), r.qual.to_vec());
    loop {
        match fasterp::read_pair(&mut first, second.as_mut(), interleaved) {
            Ok(Some((r1, r2))) => {
                out.push(triple(&r1));
                if let Some(r2) = r2 {
                    out.push(triple(&r2));
                }
            }
            Ok(None) => return (out, None),
            Err(e) => return (out, Some(e)),
        }
    }
}

/// Returns the reference error, if the input failed.
fn check_against_reference(
    rng: &mut Rng,
    inputs: &[Vec<u8>],
    interleaved: bool,
    phred64: bool,
) -> Option<String> {
    let refs: Vec<&[u8]> = inputs.iter().map(Vec::as_slice).collect();
    let (want, want_err) = reference_run(&refs, interleaved, phred64, 1 + rng.below(64));
    let small = rng.chance(50);
    let capacity = 1 + rng.below(if small { 16 } else { 4096 });
    let small = rng.chance(50);
    let budget = 1 + rng.below(if small { 64 } else { 1 << 16 });
    let opts = options(phred64)
        .batch_bytes(budget)
        .batch_records(1 + rng.below(1000));
    let (batches, err) = match (refs.len(), interleaved) {
        (1, false) => drain_partial(FastqSource::<Single>::single(
            input(refs[0], capacity),
            opts,
        )),
        (1, true) => drain_partial(
            FastqSource::<Paired>::interleaved(input(refs[0], capacity), 2, opts).unwrap(),
        ),
        _ => drain_partial(FastqSource::<Paired>::paired(
            input(refs[0], capacity),
            input(refs[1], capacity + 7),
            opts,
        )),
    };
    let got = spot_records(&batches);
    let ctx = format!(
        "capacity={capacity} budget={budget} inputs={:?}",
        inputs
            .iter()
            .map(|i| String::from_utf8_lossy(i))
            .collect::<Vec<_>>()
    );
    assert_eq!(err.as_ref().map(fasterp_message), want_err, "{ctx}");
    if want_err.is_none() {
        assert_eq!(got, want, "{ctx}");
    } else {
        assert!(
            got.len() <= want.len() && got[..] == want[..got.len()],
            "{ctx}"
        );
    }
    want_err
}

#[test]
fn random_valid_and_corrupted_inputs_match_fasterp() {
    let mut rng = Rng(0x5eed_1234_abcd_0001);
    let mut failed = 0;
    let mut kinds = std::collections::BTreeSet::new();
    for round in 0..1500 {
        let phred64 = rng.chance(20);
        let n = rng.below(12);
        let corrupt_it = round % 3 != 0;
        let outcome = match rng.below(3) {
            0 => {
                let mut a = random_fastq(&mut rng, n, phred64);
                if corrupt_it {
                    corrupt(&mut rng, &mut a);
                }
                check_against_reference(&mut rng, &[a], false, phred64)
            }
            1 => {
                let odd = usize::from(rng.chance(10));
                let mut a = random_fastq(&mut rng, n * 2 + odd, phred64);
                if corrupt_it {
                    corrupt(&mut rng, &mut a);
                }
                check_against_reference(&mut rng, &[a], true, phred64)
            }
            _ => {
                let m = if rng.chance(10) { n + 1 } else { n };
                let mut a = random_fastq(&mut rng, n, phred64);
                let mut b = random_fastq(&mut rng, m, phred64);
                if corrupt_it {
                    if rng.chance(50) {
                        corrupt(&mut rng, &mut a);
                    } else {
                        corrupt(&mut rng, &mut b);
                    }
                }
                check_against_reference(&mut rng, &[a, b], false, phred64)
            }
        };
        if let Some(e) = outcome {
            failed += 1;
            kinds.insert(
                e.split(" at ")
                    .next()
                    .unwrap()
                    .replace(|c: char| c.is_ascii_digit(), ""),
            );
        }
    }
    // Both outcomes are well represented.
    assert!((300..1200).contains(&failed), "{failed} of 1500 failed");
    // Every kind of input error occurs.
    assert_eq!(kinds.len(), 7, "{kinds:?}");
}

#[test]
fn arbitrary_bytes_never_panic() {
    let mut rng = Rng(0xdead_beef_0000_0042);
    let alphabet = b"@+\n\r \tACGTI!~\x00\xff";
    for _ in 0..4000 {
        let len = rng.below(200);
        let bytes = |rng: &mut Rng| -> Vec<u8> {
            (0..len)
                .map(|_| {
                    if rng.chance(70) {
                        alphabet[rng.below(alphabet.len())]
                    } else {
                        rng.next() as u8
                    }
                })
                .collect()
        };
        let (a, b, c) = (bytes(&mut rng), bytes(&mut rng), bytes(&mut rng));
        let capacity = 1 + rng.below(32);
        let opts = || options(rng_bool(len)).batch_bytes(1 + len % 50);
        let _ = drain_partial(FastqSource::<Single>::single(input(&a, capacity), opts()));
        let _ = drain_partial(FastqSource::<Paired>::paired(
            input(&a, capacity),
            input(&b, capacity),
            opts(),
        ));
        let _ = drain_partial(
            FastqSource::<Segments>::lockstep(
                vec![
                    input(&a, capacity),
                    input(&b, capacity),
                    input(&c, capacity),
                ],
                opts(),
            )
            .unwrap(),
        );
        let _ = drain_partial(
            FastqSource::<Segments>::interleaved(input(&a, capacity), 3, opts()).unwrap(),
        );
    }
}

fn rng_bool(x: usize) -> bool {
    x % 2 == 1
}

#[test]
fn random_valid_records_survive_a_round_trip() {
    let mut rng = Rng(77);
    for _ in 0..300 {
        let n = rng.below(40);
        let text = random_fastq(&mut rng, n, false);
        let capacity = 1 + rng.below(256);
        let budget = 1 + rng.below(4096);
        let source = FastqSource::<Single>::single(
            input(&text, capacity),
            FastqOptions::default().batch_bytes(budget),
        );
        let first = drain(source).unwrap();
        let expected = spot_records(&first);
        let written = write(first, false, readkit_fastq::PlusLine::Bare);
        // Writing is canonical: reading it back yields the same records, and
        // writing again the same bytes.
        let again = drain(FastqSource::<Single>::single(
            input(&written[0], 7),
            FastqOptions::default(),
        ))
        .unwrap();
        assert_eq!(spot_records(&again), expected);
        assert_eq!(write(again, false, readkit_fastq::PlusLine::Bare), written);
    }
}
