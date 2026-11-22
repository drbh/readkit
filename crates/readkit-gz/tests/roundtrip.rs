//! Round trips against flate2 through every decoding path, and the errors
//! damaged input must give.
#![cfg(feature = "parallel")]

mod common;

use common::*;
use readkit_gz::{
    DecodeOptions, EncodeOptions, ErrorKind, Fastq, GzReader, GzWriter, ParGzReader, Split,
    Strategy,
};
use std::io::{self, BufRead, Read, Write};

#[test]
fn single_member_large_is_identical_at_every_thread_count() {
    let data = fastq(60_000, 1);
    for (name, gz) in [
        ("flate2-6", flate2_gz(&data, 6)),
        ("libdeflate-6", libdeflate_gz(&data, 6)),
    ] {
        assert!(flate2_decode(&gz).unwrap() == data);
        let sequential = sequential(&gz).unwrap();
        assert!(sequential == data, "{name}: sequential");
        for threads in [1, 2, 4, 8] {
            let got = regions(&gz, 128 << 10, threads).unwrap();
            assert!(got == sequential, "{name}: {threads} threads");
        }
    }
}

#[test]
fn default_regions_split_a_large_member() {
    // Three default regions (6 MiB compressed): `Auto` splits.
    let data = fastq(200_000, 2);
    let gz = flate2_gz(&data, 1);
    assert!(gz.len() >= 6 << 20, "{} bytes", gz.len());
    for threads in [1, 2, 4, 8] {
        let options = DecodeOptions {
            threads: Some(threads),
            ..Default::default()
        };
        let (got, strategy) = parallel_with(&gz, &options).unwrap();
        let expected = if threads == 1 {
            Strategy::Sequential
        } else {
            Strategy::Regions
        };
        assert_eq!(strategy, expected, "{threads} threads");
        assert!(got == data, "{threads} threads");
    }
}

#[test]
fn multi_member_with_empty_members() {
    let a = fastq(3_000, 3);
    let b = fastq(500, 4);
    let gz = [
        flate2_gz(&a, 6),
        flate2_gz(b"", 6),
        libdeflate_gz(&b, 1),
        fq_gz(b""),
        libdeflate_gz(b"", 9),
        fq_gz(&a),
        flate2_gz(b"", 0),
    ]
    .concat();
    let expected = [&a[..], &b[..], &a[..]].concat();
    assert!(flate2_decode(&gz).unwrap() == expected);
    assert!(all_paths(&gz).unwrap() == expected);
    // Only empty members.
    let empty = [flate2_gz(b"", 6), flate2_gz(b"", 1)].concat();
    assert!(all_paths(&empty).unwrap().is_empty());
}

#[test]
fn bgzf_with_eof_block() {
    let data = fastq(20_000, 5);
    let gz = bgzf(&data);
    assert!(flate2_decode(&gz).unwrap() == data);
    for threads in [1, 2, 4, 8] {
        let options = DecodeOptions {
            threads: Some(threads),
            ..Default::default()
        };
        let (got, strategy) = parallel_with(&gz, &options).unwrap();
        assert_eq!(strategy, Strategy::Members);
        assert!(got == data, "{threads} threads");
    }
    assert!(all_paths(&gz).unwrap() == data);
    // The EOF block alone is an empty BGZF file.
    assert!(all_paths(&BGZF_EOF).unwrap().is_empty());
    // A corrupt block in the middle fails with the sequential decoder's error.
    let mut bad = gz.clone();
    let at = bad.len() / 2;
    bad[at] ^= 0x40;
    assert!(all_paths(&bad).is_err());
}

#[test]
fn concatenated_fastq_gz_like_ena() {
    // ENA serves runs submitted in parts as the parts' gzip files
    // concatenated: several large members from different tools, some with a
    // file name in the header.
    let parts = [fastq(15_000, 6), fastq(12_000, 7), fastq(9_000, 8)];
    let gz = [
        named_gz(&parts[0], "lane1_R1.fastq"),
        libdeflate_gz(&parts[1], 6),
        flate2_gz(&parts[2], 9),
    ]
    .concat();
    let expected = parts.concat();
    assert!(flate2_decode(&gz).unwrap() == expected);
    for threads in [1, 2, 4, 8] {
        let got = regions(&gz, 64 << 10, threads).unwrap();
        assert!(got == expected, "{threads} threads");
    }
    let options = DecodeOptions {
        threads: Some(4),
        region_size: 64 << 10,
        ..Default::default()
    };
    let (got, strategy) = parallel_with(&gz, &options).unwrap();
    assert_eq!(strategy, Strategy::Regions);
    assert!(got == expected);
    assert!(all_paths(&gz).unwrap() == expected);
}

#[test]
fn tiny_inputs() {
    let text = fastq(4, 9);
    for len in 0..=40 {
        for (name, gz) in every_compressor(&text[..len]) {
            let got = all_paths(&gz).unwrap_or_else(|e| panic!("{name}, {len} bytes: {e}"));
            assert_eq!(got, &text[..len], "{name}, {len} bytes");
        }
    }
}

#[test]
fn truncated_inputs_are_errors() {
    let data = fastq(2_000, 10);
    for (name, gz) in every_compressor(&data) {
        for cut in [
            0,
            1,
            5,
            9,
            10,
            11,
            100,
            gz.len() / 2,
            gz.len() - 9,
            gz.len() - 8,
            gz.len() - 1,
        ] {
            let e = all_paths(&gz[..cut]).unwrap_err();
            assert_eq!(e.kind(), ErrorKind::Truncated, "{name} cut at {cut}: {e}");
            assert!(flate2_decode(&gz[..cut]).is_err());
        }
    }
    // Empty input is truncated, as for flate2.
    assert_eq!(all_paths(b"").unwrap_err().kind(), ErrorKind::Truncated);
}

#[test]
fn bad_checksums_lengths_and_trailing_bytes_are_errors() {
    let data = fastq(2_000, 11);
    for (name, gz) in every_compressor(&data) {
        let n = gz.len();
        for (what, at) in [
            ("crc", n - 8),
            ("crc", n - 5),
            ("isize", n - 4),
            ("isize", n - 1),
        ] {
            let mut bad = gz.clone();
            bad[at] ^= 1;
            let e = all_paths(&bad).unwrap_err();
            assert_eq!(e.kind(), ErrorKind::Checksum, "{name} {what}");
            assert_eq!(e.offset(), Some(0));
            assert!(flate2_decode(&bad).is_err());
        }
        // A second member with a bad CRC reports where that member starts.
        let mut two = [gz.clone(), gz.clone()].concat();
        let m = two.len();
        two[m - 8] ^= 1;
        let e = all_paths(&two).unwrap_err();
        assert_eq!(e.kind(), ErrorKind::Checksum, "{name} second member");
        assert_eq!(e.offset(), Some(n as u64), "{name} second member");
        // Bytes after the last member must start another member.
        let garbage = [&gz[..], b"not gzip, not at all"].concat();
        assert_eq!(
            all_paths(&garbage).unwrap_err().kind(),
            ErrorKind::Corrupt,
            "{name}"
        );
        assert!(flate2_decode(&garbage).is_err());
        // Bit flips in the deflate data are corrupt data or a checksum
        // mismatch, never different output.
        for at in (20..n - 8).step_by(n / 13) {
            let mut bad = gz.clone();
            bad[at] ^= 0x08;
            if let Ok(out) = all_paths(&bad) {
                assert!(out == flate2_decode(&bad).unwrap(), "{name} flip at {at}");
            }
        }
    }
}

#[test]
fn non_seekable_input_through_a_pipe() {
    let parts = [fastq(5_000, 12), fastq(100, 13)];
    let gz = [
        flate2_gz(&parts[0], 6),
        named_gz(b"", "empty"),
        fq_gz(&parts[1]),
    ]
    .concat();
    let expected = parts.concat();
    // Small, uneven reads with interruptions.
    let mut out = Vec::new();
    GzReader::from_reader(Trickle::new(gz.clone(), 13))
        .read_to_end(&mut out)
        .unwrap();
    assert!(out == expected);
    // A download in progress: chunks arrive from another thread.
    let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<u8>>(1);
    let feeder = {
        let gz = gz.clone();
        std::thread::spawn(move || {
            for chunk in gz.chunks(4093) {
                tx.send(chunk.to_vec()).unwrap();
            }
        })
    };
    struct Channel(std::sync::mpsc::Receiver<Vec<u8>>, io::Cursor<Vec<u8>>);
    impl Read for Channel {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            loop {
                let n = self.1.read(buf)?;
                if n > 0 {
                    return Ok(n);
                }
                match self.0.recv() {
                    Ok(next) => self.1 = io::Cursor::new(next),
                    Err(_) => return Ok(0),
                }
            }
        }
    }
    let reader = Channel(rx, io::Cursor::new(Vec::new()));
    let mut decoded = ParGzReader::from_reader(reader, &DecodeOptions::default()).unwrap();
    assert_eq!(decoded.strategy(), Strategy::Sequential);
    let mut out = Vec::new();
    decoded.read_to_end(&mut out).unwrap();
    assert!(out == expected);
    feeder.join().unwrap();
    // A truncated stream through a pipe.
    let e = GzReader::from_reader(Trickle::new(gz[..gz.len() - 3].to_vec(), 13))
        .read_to_end(&mut Vec::new())
        .unwrap_err();
    assert_eq!(e.kind(), io::ErrorKind::UnexpectedEof);
}

#[cfg(unix)]
#[test]
fn open_reads_a_fifo_sequentially() {
    let dir = std::env::temp_dir().join(format!("readkit-gz-fifo-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let fifo = dir.join("reads.fastq.gz");
    let _ = std::fs::remove_file(&fifo);
    let status = std::process::Command::new("mkfifo")
        .arg(&fifo)
        .status()
        .unwrap();
    assert!(status.success());
    let data = fastq(3_000, 14);
    let gz = flate2_gz(&data, 6);
    let writer = {
        let fifo = fifo.clone();
        std::thread::spawn(move || std::fs::write(fifo, gz).unwrap())
    };
    let mut reader = readkit_gz::open(&fifo, &DecodeOptions::default()).unwrap();
    assert_eq!(reader.strategy(), Strategy::Sequential);
    let mut out = Vec::new();
    reader.read_to_end(&mut out).unwrap();
    writer.join().unwrap();
    assert!(out == data);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn open_a_local_file() {
    let dir = std::env::temp_dir().join(format!("readkit-gz-open-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("reads.fastq.gz");
    let data = fastq(20_000, 15);
    std::fs::write(&path, flate2_gz(&data, 6)).unwrap();
    let options = DecodeOptions {
        threads: Some(4),
        region_size: 128 << 10,
        ..Default::default()
    };
    let reader = readkit_gz::open(&path, &options).unwrap();
    assert_eq!(reader.strategy(), Strategy::Regions);
    // `BufRead` directly: records without another buffer.
    let lines = reader.lines().map(Result::unwrap).count();
    assert_eq!(lines, 4 * 20_000);
    let missing = readkit_gz::open(dir.join("missing.gz"), &options).unwrap_err();
    assert_eq!(missing.kind(), ErrorKind::Io);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn dropping_a_reader_early_stops_decoding() {
    let data = fastq(40_000, 16);
    let gz = flate2_gz(&data, 6);
    for split in [Split::Always, Split::Never] {
        let options = DecodeOptions {
            threads: Some(4),
            split,
            region_size: 64 << 10,
            ..Default::default()
        };
        let mut reader = ParGzReader::new(gz.clone(), gz.len() as u64, &options).unwrap();
        let mut head = [0u8; 1000];
        reader.read_exact(&mut head).unwrap();
        assert_eq!(head, data[..1000]);
        drop(reader); // joins the decoder without reading the rest
    }
    let bgzf = bgzf(&data);
    let mut reader = ParGzReader::new(bgzf, 0, &DecodeOptions::default()).unwrap();
    // A length of 0 means an empty input.
    assert_eq!(
        reader.read(&mut [0; 10]).unwrap_err().kind(),
        io::ErrorKind::UnexpectedEof
    );
}

#[test]
fn errors_carry_their_kind_through_io_error() {
    let gz = flate2_gz(&fastq(100, 17), 6);
    let mut bad = gz.clone();
    let n = bad.len();
    bad[n - 6] ^= 1;
    let e = GzReader::from_reader(&bad[..])
        .read_to_end(&mut Vec::new())
        .unwrap_err();
    assert_eq!(e.kind(), io::ErrorKind::InvalidData);
    let typed = readkit_gz::Error::from_io(&e).unwrap();
    assert_eq!(typed.kind(), ErrorKind::Checksum);
    // The reader keeps returning the error.
    let mut reader = GzReader::from_reader(&bad[..]);
    let _ = reader.read_to_end(&mut Vec::new());
    assert!(reader.read(&mut [0; 4]).is_err());
    // An I/O failure of the input keeps the original error as the source.
    struct Failing;
    impl Read for Failing {
        fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
            Err(io::Error::new(
                io::ErrorKind::ConnectionReset,
                "peer went away",
            ))
        }
    }
    let e = GzReader::from_reader(Failing)
        .read(&mut [0; 4])
        .unwrap_err();
    assert_eq!(e.kind(), io::ErrorKind::ConnectionReset);
    let typed = readkit_gz::Error::from_io(&e).unwrap();
    assert_eq!(typed.kind(), ErrorKind::Io);
    let source = std::error::Error::source(typed).unwrap();
    assert!(source.downcast_ref::<io::Error>().is_some());
    // Bad options.
    let options = DecodeOptions {
        threads: Some(0),
        ..Default::default()
    };
    let e = ParGzReader::new(gz.clone(), gz.len() as u64, &options).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::InvalidInput);
}

#[test]
fn writer_output_does_not_depend_on_threads() {
    let data = fastq(8_000, 18);
    for fastq_mode in [Fastq::Auto, Fastq::Always, Fastq::Never] {
        for level in [0, 1, 4, 6, 12] {
            let encode = |threads: usize| {
                let options = EncodeOptions {
                    level,
                    threads,
                    block_size: 64 << 10,
                    fastq: fastq_mode,
                    in_flight: None,
                };
                let mut writer = GzWriter::new(Vec::new(), &options).unwrap();
                for chunk in data.chunks(10_007) {
                    writer.write_all(chunk).unwrap();
                }
                writer.finish().unwrap()
            };
            let reference = encode(0);
            for threads in [1, 2, 4, 8] {
                assert!(
                    encode(threads) == reference,
                    "{fastq_mode:?} level {level}, {threads} threads"
                );
            }
            assert!(flate2_decode(&reference).unwrap() == data);
            let (out, strategy) = parallel_with(&reference, &DecodeOptions::default()).unwrap();
            assert_eq!(strategy, Strategy::Members);
            assert!(out == data);
        }
    }
    // Binned qualities at level 4 take the FASTQ encoder, and it pays off.
    let binned = |mode| {
        let options = EncodeOptions {
            level: 4,
            fastq: mode,
            ..Default::default()
        };
        let mut writer = GzWriter::new(Vec::new(), &options).unwrap();
        writer.write_all(&data).unwrap();
        writer.finish().unwrap()
    };
    assert!(readkit_gz::binned_qualities(&data));
    assert!(binned(Fastq::Auto) == binned(Fastq::Always));
    assert!(binned(Fastq::Auto) != binned(Fastq::Never));
    assert!(readkit_gz::binned_gzip(&flate2_gz(&data, 6)[..]).unwrap());
}

#[test]
fn dropping_a_writer_finishes_it() {
    let mut out = Vec::new();
    {
        let mut writer = GzWriter::new(&mut out, &EncodeOptions::default()).unwrap();
        writer.write_all(b"@r\nA\n+\nF\n").unwrap();
    }
    assert_eq!(flate2_decode(&out).unwrap(), b"@r\nA\n+\nF\n");
}
