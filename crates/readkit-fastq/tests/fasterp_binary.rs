//! Differential test against the fasterp binary with all filtering and
//! trimming disabled (`-A -G -Q -L`): its output files must equal
//! `FastqSource` → `FastqSink` output for the same input, and its error
//! messages must match ours.
//!
//! The binary comes from `FASTERP_BIN` (fasterp 84ba478 or later, which
//! skips blank lines between records). Without it the test prints a note and
//! passes, so the workspace tests do not need fasterp.
//!
//! fastp (and so fasterp) always fails empty reads as too short, even with
//! `-L`, and a pair fails with either mate. The expected output therefore
//! drops spots with an empty segment; everything else passes through
//! unchanged.

mod support;

use std::path::{Path, PathBuf};
use std::process::Command;

use readkit_batch::{BatchBuilder, Paired, Pool, Segments, Single};
use readkit_fastq::{FastqOptions, FastqSource, PlusLine};
use support::*;

fn fasterp() -> Option<PathBuf> {
    let path = PathBuf::from(std::env::var_os("FASTERP_BIN")?);
    path.is_file().then_some(path)
}

fn dir(name: &str) -> PathBuf {
    let d = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("fasterp-diff")
        .join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Runs fasterp in `d` on `inputs` (`in1`, `in2`) and returns its outputs
/// (`out1`, `out2`) or its stderr.
fn run(bin: &Path, d: &Path, inputs: &[&[u8]], extra: &[&str]) -> Result<Vec<Vec<u8>>, String> {
    let mut cmd = Command::new(bin);
    cmd.current_dir(d).args([
        "-A", "-G", "-Q", "-L", "-w", "2", "-j", "r.json", "-h", "r.html",
    ]);
    for (i, data) in inputs.iter().enumerate() {
        let name = format!("in{}.fq", i + 1);
        std::fs::write(d.join(&name), data).unwrap();
        cmd.arg(if i == 0 { "-i" } else { "-I" }).arg(name);
    }
    let outs = if inputs.len() == 2 || extra.contains(&"--interleaved_in") {
        2
    } else {
        1
    };
    for i in 0..outs {
        cmd.arg(if i == 0 { "-o" } else { "-O" })
            .arg(format!("out{}.fq", i + 1));
    }
    cmd.args(extra);
    let output = cmd.output().unwrap();
    if !output.status.success() {
        return Err(String::from_utf8_lossy(&output.stderr).into_owned());
    }
    Ok((0..outs)
        .map(|i| std::fs::read(d.join(format!("out{}.fq", i + 1))).unwrap())
        .collect())
}

/// Our output, without spots that have an empty segment.
fn ours(batches: Vec<readkit_batch::Batch<Segments>>, plus: PlusLine) -> Vec<Vec<u8>> {
    let pool = Pool::new();
    let kept = batches
        .iter()
        .map(|b| {
            let mut out = BatchBuilder::with_segments(&pool, b.segments(), b.len());
            for spot in b.spots().filter(|s| s.iter().all(|r| !r.is_empty())) {
                out.push(&spot.iter().collect::<Vec<_>>()).unwrap();
            }
            out.finish()
        })
        .collect();
    write(kept, true, plus)
}

fn sample(rng: &mut u64, tag: &str, n: usize, name_plus: bool, phred64: bool) -> Vec<u8> {
    let mut next = || {
        *rng ^= *rng << 13;
        *rng ^= *rng >> 7;
        *rng ^= *rng << 17;
        *rng
    };
    let mut out = Vec::new();
    for i in 0..n {
        let len = (next() % 160) as usize;
        let name = format!("{tag}.{i} {i} length={len}");
        let eol: &[u8] = if next() % 7 == 0 { b"\r\n" } else { b"\n" };
        if next() % 11 == 0 {
            out.extend_from_slice(b"\n \t\n");
        }
        out.push(b'@');
        out.extend_from_slice(name.as_bytes());
        out.extend_from_slice(eol);
        for _ in 0..len {
            out.push(b"ACGTNacgtn"[(next() % 10) as usize]);
        }
        out.extend_from_slice(eol);
        out.push(b'+');
        if name_plus {
            out.extend_from_slice(name.as_bytes());
        }
        out.extend_from_slice(eol);
        let lo = if phred64 { 64 } else { 33 };
        for _ in 0..len {
            out.push(lo + (next() % (127 - lo as u64)) as u8);
        }
        out.extend_from_slice(eol);
    }
    out
}

fn plus(name_plus: bool) -> PlusLine {
    if name_plus {
        PlusLine::Name
    } else {
        PlusLine::Bare
    }
}

#[test]
fn fasterp_binary_passes_reads_through_like_readkit() {
    let Some(bin) = fasterp() else {
        eprintln!("fasterp binary not found; set FASTERP_BIN to run this test");
        return;
    };
    let mut rng = 0x1234_5678_9abc_def1u64;
    for (case, name_plus, phred64) in [(0, false, false), (1, true, false), (2, true, true)] {
        let extra: &[&str] = if phred64 { &["-6"] } else { &[] };
        let r1 = sample(&mut rng, "SRR1", 3000, name_plus, phred64);
        let r2 = sample(&mut rng, "SRR2", 3000, name_plus, phred64);
        let opts = || options(phred64).batch_bytes(10_000);

        let d = dir(&format!("single{case}"));
        let want = run(&bin, &d, &[&r1], extra).unwrap();
        let source = FastqSource::<Single>::single(input(&r1, 4096), opts());
        assert_eq!(
            ours(drain(source).unwrap(), plus(name_plus)),
            want,
            "single {case}"
        );

        let d = dir(&format!("paired{case}"));
        let want = run(&bin, &d, &[&r1, &r2], extra).unwrap();
        let source = FastqSource::<Paired>::paired(input(&r1, 4096), input(&r2, 999), opts());
        assert_eq!(
            ours(drain(source).unwrap(), plus(name_plus)),
            want,
            "paired {case}"
        );

        let d = dir(&format!("interleaved{case}"));
        let extra_i: Vec<&str> = extra.iter().copied().chain(["--interleaved_in"]).collect();
        let want = run(&bin, &d, &[&r1], &extra_i).unwrap();
        let source = FastqSource::<Paired>::interleaved(input(&r1, 4096), 2, opts()).unwrap();
        assert_eq!(
            ours(drain(source).unwrap(), plus(name_plus)),
            want,
            "interleaved {case}"
        );
    }
}

#[test]
fn fasterp_binary_reports_the_same_errors() {
    let Some(bin) = fasterp() else {
        eprintln!("fasterp binary not found; set FASTERP_BIN to run this test");
        return;
    };
    let ok = "@a\nAC\n+\nII\n@b\nA\n+\nI\n";
    let cases: [(&str, Option<&str>); 9] = [
        ("@a\nAC\n+\nII\n\n@b\n\nA\n+\nI\n", None),
        ("@a\nAC\n+\nII\nb\nA\n+\nI\n", None),
        ("@a\nAC\n+\nII\n@b\nA\nx\nI\n", None),
        ("@a\nAC\n+\nII\n@b\nA\n+\nII\n", None),
        ("@a\nAC\n+\nII\n@b\nA\n+\n \n", None),
        ("@a\nAC\n+\nII\n@b\nA\n+\n", None),
        ("@a\nAC\n+\nII\n\n   @b\nA\n+\nI\n", None),
        (ok, Some("@a\nAC\n+\nII\n")),
        ("@a\nAC\n+\nII\n", Some(ok)),
    ];
    for (i, (first, second)) in cases.into_iter().enumerate() {
        let d = dir(&format!("error{i}"));
        let inputs: Vec<&[u8]> = [Some(first), second]
            .into_iter()
            .flatten()
            .map(str::as_bytes)
            .collect();
        let stderr = run(&bin, &d, &inputs, &[]).unwrap_err();
        let e = match inputs.len() {
            1 => drain(FastqSource::<Single>::single(
                input(inputs[0], 64),
                FastqOptions::default(),
            )),
            _ => drain(FastqSource::<Paired>::paired(
                input(inputs[0], 64),
                input(inputs[1], 64),
                FastqOptions::default(),
            )),
        }
        .unwrap_err();
        assert!(
            stderr.contains(&fasterp_message(&e)),
            "case {i}: fasterp said {stderr:?}, readkit {e}"
        );
    }
    // Odd interleaved count.
    let d = dir("error-interleaved");
    let stderr = run(&bin, &d, &[b"@a\nA\n+\nI\n"], &["--interleaved_in"]).unwrap_err();
    let e = drain(
        FastqSource::<Paired>::interleaved(input(b"@a\nA\n+\nI\n", 64), 2, FastqOptions::default())
            .unwrap(),
    )
    .unwrap_err();
    assert!(stderr.contains(&fasterp_message(&e)), "{stderr:?}");
}
