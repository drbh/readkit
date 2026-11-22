//! FASTQ → batch → FASTQ is byte-identical to fasterp's output (its
//! vendored reader and `Read::write`) for the same input.

mod support;

use readkit_batch::{Paired, Segments, Single};
use readkit_fastq::{FastqOptions, FastqSink, FastqSource, PlusLine};
use readkit_io::{BatchSink, pump};
use support::*;

/// Records with CRLF, empty reads, descriptions, IUPAC and lowercase bases,
/// blank lines between records and no final newline. `name_plus` repeats
/// the name on the `+` line.
fn sample(tag: &str, n: usize, name_plus: bool) -> String {
    let mut out = String::from("\n");
    for i in 0..n {
        let name = format!("{tag}.{i} {i} length={}", i % 13);
        let seq: String = (0..i % 13)
            .map(|j| b"ACGTNacgtRYKM"[(i + j) % 13] as char)
            .collect();
        let qual: String = (0..i % 13)
            .map(|j| (b'!' + ((i * 7 + j) % 94) as u8) as char)
            .collect();
        let plus = if name_plus { name.as_str() } else { "" };
        let eol = if i % 5 == 0 { "\r\n" } else { "\n" };
        out += &format!("@{name}{eol}{seq}{eol}+{plus}{eol}{qual}");
        // No final newline, unless the last quality line is empty (it would
        // be a truncated record).
        if i + 1 < n || qual.is_empty() {
            out += eol;
        }
        if i % 9 == 4 {
            out += "\n  \t\n";
        }
    }
    out
}

fn style(name_plus: bool) -> PlusLine {
    if name_plus {
        PlusLine::Name
    } else {
        PlusLine::Bare
    }
}

#[test]
fn single_end_matches_fasterp() {
    for name_plus in [false, true] {
        let text = sample("SRR1", 200, name_plus);
        let expected = fasterp_expected(&[text.as_bytes()], false, false);
        assert_eq!(expected.error, None);
        for capacity in [1, 3, 64, 1 << 16] {
            for budget in [1, 500, 1 << 20] {
                let source = FastqSource::<Single>::single(
                    input(text.as_bytes(), capacity),
                    FastqOptions::default().batch_bytes(budget),
                );
                let mut out = Vec::new();
                let sink = FastqSink::new(&mut out).with_plus_line(style(name_plus));
                assert_eq!(pump(source, sink).unwrap(), 200);
                assert_eq!(out, expected.out[0], "capacity={capacity} budget={budget}");
            }
        }
    }
}

#[test]
fn paired_matches_fasterp() {
    for name_plus in [false, true] {
        let (r1, r2) = (
            sample("SRR1", 150, name_plus),
            sample("SRR2", 150, name_plus),
        );
        let expected = fasterp_expected(&[r1.as_bytes(), r2.as_bytes()], false, false);
        assert_eq!(expected.error, None);
        for capacity in [1, 5, 128, 1 << 16] {
            for budget in [1, 700, 1 << 20] {
                let source = FastqSource::<Paired>::paired(
                    input(r1.as_bytes(), capacity),
                    input(r2.as_bytes(), capacity + 3),
                    FastqOptions::default().batch_bytes(budget),
                );
                let mut out = [Vec::new(), Vec::new()];
                let [o1, o2] = &mut out;
                let sink = FastqSink::split(vec![o1, o2]).with_plus_line(style(name_plus));
                pump(source, sink).unwrap();
                assert_eq!(
                    out.to_vec(),
                    expected.out,
                    "capacity={capacity} budget={budget}"
                );
            }
        }
    }
}

#[test]
fn interleaved_matches_fasterp() {
    for name_plus in [false, true] {
        let text = sample("SRR3", 300, name_plus);
        let expected = fasterp_expected(&[text.as_bytes()], true, false);
        assert_eq!(expected.error, None);
        for capacity in [1, 9, 1 << 16] {
            for budget in [1, 300, 1 << 20] {
                let source = FastqSource::<Paired>::interleaved(
                    input(text.as_bytes(), capacity),
                    2,
                    FastqOptions::default().batch_bytes(budget),
                )
                .unwrap();
                let mut out = Vec::new();
                let sink = FastqSink::new(&mut out).with_plus_line(style(name_plus));
                assert_eq!(pump(source, sink).unwrap(), 150);
                assert_eq!(out, expected.out[0]);
            }
        }
    }
}

#[test]
fn four_segments_match_fasterp_per_file() {
    // fasterp reads at most two files; each of the four outputs must equal
    // fasterp's single-end output for that file.
    let files: Vec<String> = ["R1", "R2", "I1", "I2"]
        .iter()
        .map(|t| sample(t, 120, false))
        .collect();
    for capacity in [1, 17, 1 << 16] {
        let inputs = files
            .iter()
            .map(|f| input(f.as_bytes(), capacity))
            .collect();
        let source =
            FastqSource::<Segments>::lockstep(inputs, FastqOptions::default().batch_bytes(999))
                .unwrap();
        let mut sink = FastqSink::split(vec![Vec::new(); 4]);
        let batches = drain(source).unwrap();
        // Interleaved output of the same batches.
        let interleaved = write(
            batches.iter().map(clone_batch).collect(),
            false,
            PlusLine::Bare,
        );
        for b in batches {
            sink.accept(b).unwrap();
        }
        let out = sink.into_inner().unwrap();
        for (o, f) in out.iter().zip(&files) {
            assert_eq!(*o, fasterp_expected(&[f.as_bytes()], false, false).out[0]);
        }
        // 4-way interleaved round trip: read back the interleaved output.
        let source = FastqSource::<Segments>::interleaved(
            input(&interleaved[0], capacity),
            4,
            FastqOptions::default(),
        )
        .unwrap();
        let again = write(drain(source).unwrap(), true, PlusLine::Bare);
        assert_eq!(again, out);
    }
}

#[test]
fn phred64_is_written_as_phred33_like_fasterp() {
    let text = "@a\r\nACGT\r\n+\r\n@@hh\r\n\n@b\nN\n+\n~\n";
    let expected = fasterp_expected(&[text.as_bytes()], false, true);
    assert_eq!(expected.out[0], b"@a\nACGT\n+\n!!II\n@b\nN\n+\n_\n");
    for capacity in [1, 4, 64] {
        let source = FastqSource::<Single>::single(input(text.as_bytes(), capacity), options(true));
        let mut out = Vec::new();
        pump(source, FastqSink::new(&mut out)).unwrap();
        assert_eq!(out, expected.out[0]);
    }
}

#[test]
fn nameless_batches_are_numbered_and_plus_styles_differ() {
    let text = b"@x\nAC\n+x\nII\n@y\nG\n+\nI\n";
    let source =
        FastqSource::<Single>::single(input(text, 64), FastqOptions::default().keep_names(false));
    let out = write(drain(source).unwrap(), false, PlusLine::Name);
    assert_eq!(out[0], b"@1\nAC\n+1\nII\n@2\nG\n+2\nI\n");
    let source = FastqSource::<Single>::single(input(text, 64), FastqOptions::default());
    let out = write(drain(source).unwrap(), false, PlusLine::Name);
    assert_eq!(out[0], b"@x\nAC\n+x\nII\n@y\nG\n+y\nI\n");
    let source = FastqSource::<Single>::single(input(text, 64), FastqOptions::default());
    let out = write(drain(source).unwrap(), false, PlusLine::Bare);
    assert_eq!(out[0], b"@x\nAC\n+\nII\n@y\nG\n+\nI\n");
}

#[test]
fn split_sink_checks_its_writer_count() {
    let source =
        FastqSource::<Single>::single(input(b"@x\nA\n+\nI\n", 64), FastqOptions::default());
    let sink: FastqSink<Single, Vec<u8>> = FastqSink::split(vec![Vec::new(), Vec::new()]);
    let e = pump(source, sink).unwrap_err();
    let cause = e.find_cause::<readkit_fastq::Error>().unwrap();
    assert_eq!(cause.kind(), readkit_fastq::ErrorKind::Options);
}

fn clone_batch(b: &readkit_batch::Batch<Segments>) -> readkit_batch::Batch<Segments> {
    let pool = readkit_batch::Pool::new();
    let mut builder = readkit_batch::BatchBuilder::with_segments(&pool, b.segments(), b.len());
    for spot in b.spots() {
        builder.push(&spot.iter().collect::<Vec<_>>()).unwrap();
    }
    builder.finish()
}

#[test]
fn works_behind_read_ahead() {
    let (r1, r2) = (sample("A", 500, true), sample("B", 500, true));
    let expected = fasterp_expected(&[r1.as_bytes(), r2.as_bytes()], false, false);
    let source = FastqSource::<Paired>::paired(
        input(r1.as_bytes(), 256),
        input(r2.as_bytes(), 256),
        FastqOptions::default().batch_bytes(1000),
    );
    let source = readkit_io::ReadAhead::new(source);
    assert!(source.is_threaded());
    let mut out = [Vec::new(), Vec::new()];
    let [o1, o2] = &mut out;
    pump(
        source,
        FastqSink::split(vec![o1, o2]).with_plus_line(PlusLine::Name),
    )
    .unwrap();
    assert_eq!(out.to_vec(), expected.out);
}
