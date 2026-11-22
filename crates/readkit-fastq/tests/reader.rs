//! fasterp's reader behaviour tests (`compat.rs` and `compat/block.rs`),
//! ported to `FastqSource`, plus the k-segment generalization.

mod support;

use readkit_batch::{Paired, Segments, Single};
use readkit_fastq::{ErrorKind, FastqInput, FastqOptions, FastqSource};
use readkit_io::BatchSource;
use support::*;

fn single(data: &[u8], capacity: usize, opts: FastqOptions) -> FastqSource<Single> {
    FastqSource::single(input(data, capacity), opts)
}

fn first_error(data: &[u8], capacity: usize) -> readkit_fastq::Error {
    drain_partial(single(data, capacity, FastqOptions::default()))
        .1
        .expect("input should fail")
}

#[test]
fn buffered_reader_preserves_lines_at_every_buffer_boundary() {
    let input =
        b"@first description\r\nACGTN\r\n+first\r\n@ABCD\r\n@empty\n\n+\n\n@last\nAC\n+\nHH";
    for capacity in 1..=input.len() + 1 {
        for budget in [1, 31, 1 << 20] {
            let opts = options(true).batch_bytes(budget);
            let batches = drain(single(input, capacity, opts)).unwrap();
            let got = records(&batches, 0);
            let want: Vec<Triple> = vec![
                (
                    b"first description".to_vec(),
                    b"ACGTN".to_vec(),
                    b"!\"#$%".to_vec(),
                ),
                (b"empty".to_vec(), vec![], vec![]),
                (b"last".to_vec(), b"AC".to_vec(), b"))".to_vec()),
            ];
            assert_eq!(got, want, "capacity={capacity} budget={budget}");
        }
    }
}

#[test]
fn reader_skips_blank_lines_between_records_like_fastp() {
    let spaces = " ".repeat(40);
    let input = format!(
        "\n \t\r\n@a\nAC\n+\nII\n\n\n@b\nT\n+\nI\n{spaces}\n\x0b\x0c\r\n@c\nG\n+\nH\n\n\t "
    );
    for capacity in 1..=input.len() + 1 {
        let mut source = single(input.as_bytes(), capacity, FastqOptions::default());
        let mut batches = Vec::new();
        while let Some(b) = source.next_batch().unwrap() {
            batches.push(b.into_segments());
        }
        assert_eq!(
            names(&batches, 0),
            [&b"a"[..], b"b", b"c"],
            "capacity={capacity}"
        );
        assert_eq!(source.spots(), 3);
    }
    for input in [&b""[..], b"\n", b"  \n\t\r\n  "] {
        assert!(
            drain(single(input, 1, FastqOptions::default()))
                .unwrap()
                .is_empty()
        );
    }
}

#[test]
fn reader_rejects_blank_lines_inside_records_at_the_right_record() {
    for (input, kind, error) in [
        (
            &b"\n@a\nA\n+\nI\n\n@b\n\nA\n+\nI\n"[..],
            ErrorKind::Separator,
            "Invalid FASTQ separator at record 2",
        ),
        (
            b"@a\nA\n+\nI\n@b\nA\n\n+\nI\n",
            ErrorKind::Separator,
            "Invalid FASTQ separator at record 2",
        ),
        (
            b"@a\nA\n+\nI\n@b\nA\n+\n\nI\n",
            ErrorKind::LengthMismatch,
            "Sequence/quality length mismatch at record 2",
        ),
        (
            b"@a\nA\n+\nI\n\n  @b\nA\n+\nI\n",
            ErrorKind::Header,
            "Invalid FASTQ header at record 2",
        ),
        (
            b"@a\nA\n+\nI\n\n@b\nA\n+\n",
            ErrorKind::Truncated,
            "Truncated FASTQ record 2",
        ),
    ] {
        for capacity in 1..=input.len() + 1 {
            let e = first_error(input, capacity);
            assert_eq!(e.kind(), kind, "capacity={capacity}, input={input:?}");
            assert_eq!(e.record(), Some(2));
            assert_eq!(fasterp_message(&e), error);
            // And the reference agrees.
            assert_eq!(
                fasterp_expected(&[input], false, false).error.as_deref(),
                Some(error)
            );
        }
    }
}

#[test]
fn every_validation_error_names_its_record() {
    let ok = "@ok\nACGT\n+\nIIII\n";
    for (bad, kind) in [
        (&b"wrong\nA\n+\nI\n"[..], ErrorKind::Header),
        (b"@r\nA\nwrong\nI\n", ErrorKind::Separator),
        (b"@r\nACG\nT\n+\nIII\nI\n", ErrorKind::Separator), // multi-line record
        (b"@r\nAC\n+\nI\n", ErrorKind::LengthMismatch),
        (b"@r\nA\n+\n \n", ErrorKind::Quality),
        (b"@r\nA\n+\n\x7f\n", ErrorKind::Quality),
        (b"@r\nA\n+\n\xff\n", ErrorKind::Quality),
        (b"@r\nA\n+\n", ErrorKind::Truncated),
        (b"@r\nA\n", ErrorKind::Truncated),
        (b"@r\n", ErrorKind::Truncated),
        (b"@r", ErrorKind::Truncated),
    ] {
        for before in 0..3 {
            let data = [ok.repeat(before).as_bytes(), bad].concat();
            for capacity in 1..=data.len() + 1 {
                let e = first_error(&data, capacity);
                assert_eq!(
                    (e.kind(), e.record(), e.segment()),
                    (kind, Some(before as u64 + 1), None),
                    "data={data:?} capacity={capacity}"
                );
                let reference = fasterp_expected(&[&data], false, false);
                assert_eq!(reference.error, Some(fasterp_message(&e)), "data={data:?}");
            }
        }
    }
}

#[test]
fn quality_range_depends_on_phred64() {
    // 33..=126 for Phred+33, 64..=126 for Phred+64 (then converted).
    let data = |q: u8| [&b"@r\nA\n+\n"[..], &[q], b"\n"].concat();
    for q in 0..=255u8 {
        for phred64 in [false, true] {
            let r = drain(single(&data(q), 64, options(phred64)));
            let lo = if phred64 { 64 } else { 33 };
            if (lo..=126).contains(&q) {
                let batches = r.unwrap();
                let want = if phred64 { q - 31 } else { q };
                assert_eq!(batches[0].segment(0).quals().get(0), [want]);
            } else if q == b'\n' || q == b'\r' {
                // An empty quality line after stripping.
                assert_eq!(r.unwrap_err().kind(), ErrorKind::LengthMismatch, "q={q}");
            } else {
                assert_eq!(r.unwrap_err().kind(), ErrorKind::Quality, "q={q}");
            }
        }
    }
}

#[test]
fn crlf_is_stripped_once_per_line() {
    let data = b"@a b\r\nAC\r\n+a b\r\nII\r\n@c\r\r\nG\r\n+\r\nI\r\n";
    for capacity in 1..=data.len() + 1 {
        let batches = drain(single(data, capacity, FastqOptions::default())).unwrap();
        assert_eq!(
            records(&batches, 0),
            vec![
                (b"a b".to_vec(), b"AC".to_vec(), b"II".to_vec()),
                (b"c\r".to_vec(), b"G".to_vec(), b"I".to_vec()),
            ]
        );
    }
}

#[test]
fn names_can_be_dropped() {
    let data = b"@a\nAC\n+\nII\n@b\n\n+\n\n";
    let batches = drain(single(data, 8, FastqOptions::default().keep_names(false))).unwrap();
    assert!(!batches[0].has_names());
    assert_eq!(batches[0].segment(0).seqs().offsets(), &[0, 2, 2]);
}

#[test]
fn errors_name_the_segment_and_input() {
    let ok = b"@a\nA\n+\nI\n@b\nA\n+\nI\n";
    let bad = b"@a\nA\n+\nI\n@b\nA\n+\n!!\n";
    let source = FastqSource::<Paired>::paired(
        named(ok, 64, "r1.fq"),
        named(bad, 64, "r2.fq"),
        FastqOptions::default(),
    );
    let e = drain(source).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::LengthMismatch);
    assert_eq!(
        (e.record(), e.segment(), e.input()),
        (Some(2), Some(1), Some("r2.fq"))
    );
    assert_eq!(
        e.to_string(),
        "sequence/quality length mismatch at record 2, segment 1 in r2.fq"
    );
    let e = first_error(b"@a\nA\n+\nI\n\n x\n", 4);
    assert_eq!(e.to_string(), "invalid FASTQ header at record 2");
}

// ---- compat/block.rs tests ----

#[test]
fn interleaved_blocks_match_owned_records_at_every_buffer_boundary() {
    let data = format!(
        "@first/1\r\nACG\r\n+first\r\nIII\r\n@first/2\n{}\n+\n{}\n@empty/1\n\n+\n\n@empty/2\nT\n+\nH",
        "G".repeat(97),
        "I".repeat(97)
    );
    for phred64 in [false, true] {
        let expected = reference_records(data.as_bytes(), phred64);
        for capacity in 1..=data.len() + 1 {
            for budget in [1, 31, 256, 20000] {
                let source = FastqSource::<Paired>::interleaved(
                    input(data.as_bytes(), capacity),
                    2,
                    options(phred64).batch_bytes(budget),
                )
                .unwrap();
                let batches = drain(source).unwrap();
                assert_eq!(
                    spot_records(&batches),
                    expected,
                    "capacity={capacity} budget={budget} phred64={phred64}"
                );
            }
        }
    }
}

#[test]
fn interleaved_blocks_count_limits_in_pairs_and_reject_dangling_mates() {
    let data = b"@a/1\nAC\n+\nII\n@a/2\nG\n+\nI\n@b/1\nT\n+\nH\n";
    for capacity in 1..=data.len() + 1 {
        let mut source = FastqSource::<Paired>::interleaved(
            input(data, capacity),
            2,
            FastqOptions::default().batch_records(1),
        )
        .unwrap();
        let first = source.next_batch().unwrap().unwrap();
        assert_eq!(first.len(), 1);
        assert_eq!(first.r1().names().unwrap().get(0), b"a/1");
        assert_eq!(first.r2().names().unwrap().get(0), b"a/2");
        let e = source.next_batch().unwrap_err();
        assert_eq!(e.kind(), ErrorKind::Interleaved);
        assert_eq!((e.record(), e.segment()), (Some(4), Some(1)));
    }
}

#[test]
fn blocks_match_owned_reader_at_boundaries_and_asymmetric_lengths() {
    let first = b"@one note\r\nACGTN\r\n+one\r\n@ABCD\r\n@empty\n\n+\n\n@last\nAC\n+\nHH";
    let second = format!(
        "@one\n{}\n+\n{}\n@empty\r\nA\r\n+note\r\nH\r\n@last\nT\n+\nH",
        "G".repeat(4097),
        "I".repeat(4097)
    );
    for phred64 in [false, true] {
        let expected = [
            reference_records(first, phred64),
            reference_records(second.as_bytes(), phred64),
        ];
        for capacity in 1..=first.len() + 1 {
            for budget in [1, 31, 256, 20000] {
                let source = FastqSource::<Paired>::paired(
                    input(first, capacity),
                    input(second.as_bytes(), capacity * 2 + 1),
                    options(phred64).batch_bytes(budget),
                );
                let batches = drain(source).unwrap();
                assert_eq!(
                    [records(&batches, 0), records(&batches, 1)],
                    expected,
                    "capacity={capacity}, budget={budget}, phred64={phred64}"
                );
                if budget == 1 {
                    assert!(batches.iter().all(|b| b.len() == 1));
                }
            }
        }
    }
}

#[test]
fn blocks_stop_at_read_limit_and_reject_unbalanced_mates() {
    let data = b"@one\nAC\n+\nII\n@two\nT\n+\nI\n";
    for capacity in 1..=data.len() + 1 {
        let mut source = FastqSource::<Paired>::paired(
            input(data, capacity),
            input(&data[..13], capacity),
            FastqOptions::default().batch_records(1),
        );
        let batch = source.next_batch().unwrap().unwrap();
        assert_eq!(batch.len(), 1);
        let e = source.next_batch().unwrap_err();
        assert_eq!(e.kind(), ErrorKind::MateCount);
        assert_eq!((e.record(), e.segment()), (Some(2), Some(1)));
        // The other way round.
        let source = FastqSource::<Paired>::paired(
            input(&data[..13], capacity),
            input(data, capacity),
            FastqOptions::default(),
        );
        let e = drain(source).unwrap_err();
        assert_eq!(
            (e.kind(), e.record(), e.segment()),
            (ErrorKind::MateCount, Some(2), Some(0))
        );
    }
}

#[test]
fn blocks_skip_blank_lines_between_records_in_every_mode() {
    let blank = format!("\n \t\r\n{}\n", " ".repeat(40));
    let recs = [
        "@a\nAC\n+\nII\n",
        "@b\nT\n+\nI\n",
        "@c\nG\n+\nH\n",
        "@d\nA\n+\nI",
    ];
    let clean = recs.concat();
    let messy = format!("{blank}{}", recs.join(&blank)) + "\n\n  ";
    let all: Vec<Vec<u8>> = ["a", "b", "c", "d"].map(|n| n.as_bytes().to_vec()).into();
    for capacity in 1..=messy.len() + 1 {
        for budget in [1, 31, 20000] {
            // Blank lines in either mate, or both.
            for (m1, m2) in [(&messy, &clean), (&clean, &messy), (&messy, &messy)] {
                let source = FastqSource::<Paired>::paired(
                    input(m1.as_bytes(), capacity),
                    input(m2.as_bytes(), capacity),
                    FastqOptions::default().batch_bytes(budget),
                );
                let batches = drain(source).unwrap();
                assert_eq!(
                    [names(&batches, 0), names(&batches, 1)],
                    [all.clone(), all.clone()],
                    "capacity={capacity}"
                );
            }
            let source = FastqSource::<Paired>::interleaved(
                input(messy.as_bytes(), capacity),
                2,
                FastqOptions::default().batch_bytes(budget),
            )
            .unwrap();
            let batches = drain(source).unwrap();
            let got: Vec<_> = spot_records(&batches).into_iter().map(|t| t.0).collect();
            assert_eq!(got, all, "capacity={capacity}");
        }
    }
    for data in [&b""[..], b"\n", b"  \n\t\r\n  "] {
        let source =
            FastqSource::<Paired>::paired(input(data, 1), input(data, 1), FastqOptions::default());
        assert!(drain(source).unwrap().is_empty());
    }
}

#[test]
fn blocks_report_blank_lines_inside_records_like_owned_reader() {
    for data in [
        &b"\n@a\nA\n+\nI\n\n@b\n\nA\n+\nI\n@c\nA\n+\nI\n"[..],
        b"@a\nA\n+\nI\n@b\nA\n\n+\nI\n",
        b"@a\nA\n+\nI\n@b\nA\n+\n\nI\n@c\nA\n+\nI\n",
        b"@a\nA\n+\nI\n\n  @b\nA\n+\nI\n",
    ] {
        let expected = fasterp_expected(&[data], false, false).error.unwrap();
        for capacity in 1..=data.len() + 1 {
            for budget in [1, 1024] {
                let e = drain(single(
                    data,
                    capacity,
                    FastqOptions::default().batch_bytes(budget),
                ))
                .unwrap_err();
                assert_eq!(
                    fasterp_message(&e),
                    expected,
                    "capacity={capacity}, data={data:?}"
                );
            }
        }
    }
}

#[test]
fn blocks_reject_the_same_malformed_records_as_owned_reader() {
    for data in [
        &b"wrong\nA\n+\nI\n"[..],
        b"@r\nA\nwrong\nI\n",
        b"@r\nAC\n+\nI\n",
        b"@r\nA\n+\n \n",
        b"@r\nA\n+\n\xff\n",
        b"@r\nA\n+\n",
        b"@r\n",
    ] {
        let expected = reference(data, 7, false).next().unwrap_err();
        for capacity in 1..=data.len() + 1 {
            let e = first_error(data, capacity);
            assert_eq!(
                fasterp_message(&e),
                expected,
                "capacity={capacity}, data={data:?}"
            );
        }
    }
}

// ---- k segments ----

fn spot(names: &[&str], i: usize) -> Vec<String> {
    names
        .iter()
        .enumerate()
        .map(|(s, n)| {
            let len = (i * 7 + s * 3) % 11;
            format!(
                "@{n}{i}\n{}\n+\n{}\n",
                "ACGTN".repeat(5)[..len].to_owned(),
                "I".repeat(len)
            )
        })
        .collect()
}

#[test]
fn four_segments_in_lockstep() {
    let segs = ["R1_", "R2_", "I1_", "I2_"];
    let mut files = vec![String::new(); 4];
    for i in 0..50 {
        for (f, rec) in files.iter_mut().zip(spot(&segs, i)) {
            f.push_str(&rec);
        }
    }
    files[2] = format!("\n\n{}", files[2]); // blank lines in one input
    for capacity in [1, 5, 64, 4096] {
        for budget in [1, 100, 1 << 20] {
            let inputs = files
                .iter()
                .map(|f| input(f.as_bytes(), capacity))
                .collect();
            let source = FastqSource::<Segments>::lockstep(
                inputs,
                FastqOptions::default().batch_bytes(budget),
            )
            .unwrap();
            let batches = drain(source).unwrap();
            assert!(batches.iter().all(|b| b.segments() == 4));
            for (s, f) in files.iter().enumerate() {
                assert_eq!(records(&batches, s), reference_records(f.as_bytes(), false));
            }
        }
    }
    // One input short by one record.
    let mut short = files.clone();
    let cut = short[3].rfind('@').unwrap();
    short[3].truncate(cut);
    for capacity in [1, 64, 4096] {
        let inputs = short
            .iter()
            .enumerate()
            .map(|(s, f)| named(f.as_bytes(), capacity, &format!("seg{s}.fq")))
            .collect();
        let source = FastqSource::<Segments>::lockstep(inputs, FastqOptions::default()).unwrap();
        let e = drain(source).unwrap_err();
        assert_eq!(
            (e.kind(), e.record(), e.segment(), e.input()),
            (ErrorKind::MateCount, Some(50), Some(3), Some("seg3.fq"))
        );
    }
}

#[test]
fn four_way_interleaved() {
    let segs = ["R1_", "R2_", "I1_", "I2_"];
    let text: String = (0..30).flat_map(|i| spot(&segs, i)).collect();
    let expected = reference_records(text.as_bytes(), false);
    for capacity in [1, 7, 64, 4096] {
        for budget in [1, 100, 1 << 20] {
            let source = FastqSource::<Segments>::interleaved(
                input(text.as_bytes(), capacity),
                4,
                FastqOptions::default().batch_bytes(budget),
            )
            .unwrap();
            let batches = drain(source).unwrap();
            assert_eq!(batches.iter().map(|b| b.len()).sum::<usize>(), 30);
            assert_eq!(spot_records(&batches), expected);
        }
        // 30 spots plus two records: not a multiple of 4.
        let extra = text.clone() + &spot(&segs, 99)[..2].concat();
        let source = FastqSource::<Segments>::interleaved(
            input(extra.as_bytes(), capacity),
            4,
            FastqOptions::default(),
        )
        .unwrap();
        let e = drain(source).unwrap_err();
        assert_eq!(
            (e.kind(), e.record(), e.segment()),
            (ErrorKind::Interleaved, Some(123), Some(2))
        );
    }
}

#[test]
fn segment_counts_are_checked() {
    let data = &b"@a\nA\n+\nI\n"[..];
    let e =
        FastqSource::<Paired>::lockstep(vec![input(data, 8)], FastqOptions::default()).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Options);
    let e =
        FastqSource::<Single>::interleaved(input(data, 8), 2, FastqOptions::default()).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Options);
    let e = FastqSource::<Segments>::lockstep(Vec::<FastqInput>::new(), FastqOptions::default())
        .unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Options);
    let e = FastqSource::<Segments>::interleaved(input(data, 8), 0, FastqOptions::default())
        .unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Options);
    // k = 1 through the runtime layout works.
    let s =
        FastqSource::<Segments>::lockstep(vec![input(data, 8)], FastqOptions::default()).unwrap();
    assert_eq!(drain(s).unwrap()[0].len(), 1);
}

#[test]
fn batches_respect_the_byte_budget_and_recycle_buffers() {
    let text: String = (0..1000)
        .map(|i| format!("@r{i}\n{}\n+\n{}\n", "A".repeat(50), "I".repeat(50)))
        .collect();
    let pool = readkit_batch::Pool::new();
    let opts = FastqOptions::default().batch_bytes(4096).pool(&pool);
    let mut source = FastqSource::single(input(text.as_bytes(), 1 << 16), opts);
    let mut total = 0;
    while let Some(b) = source.next_batch().unwrap() {
        let bytes: usize = (0..b.len())
            .map(|i| {
                let r = b.column().read(i);
                r.len() * 2 + r.name().unwrap().len() + 5
            })
            .sum();
        assert!(bytes <= 4096 + 120, "{bytes}");
        total += b.len();
    }
    assert_eq!(total, 1000);
    let stats = pool.stats();
    assert!(stats.reused > stats.fresh, "{stats:?}");
}

/// The byte budget is exact on the slow path too: batch
/// sizes depend on the text and the options only, not on where the input's
/// buffer ends, and a spot that does not fit carries its error into the next
/// batch.
#[test]
fn batches_do_not_depend_on_the_input_buffer() {
    let mut text = Vec::new();
    for i in 0..40 {
        let seq = "ACGT".repeat(1 + i % 7);
        text.extend_from_slice(format!("@r{i}\n{seq}\n+\n{}\n", "I".repeat(seq.len())).as_bytes());
        if i == 20 {
            text.extend_from_slice(b"\n \n");
        }
    }
    let mut bad = text.clone();
    let last_qual = bad.len() - 2;
    bad[last_qual] = b' ';
    for data in [&text, &bad] {
        for budget in [1, 40, 97, 300] {
            let sizes = |capacity: usize| {
                let (batches, error) =
                    drain_partial(single(data, capacity, options(true).batch_bytes(budget)));
                (
                    batches.iter().map(|b| b.len()).collect::<Vec<_>>(),
                    error.map(|e| e.to_string()),
                )
            };
            let want = sizes(1 << 20);
            for capacity in 1..=200 {
                assert_eq!(sizes(capacity), want, "budget={budget} capacity={capacity}");
            }
        }
    }
}
