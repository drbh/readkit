//! Members that record their size (BGZF, readkit-gz's `GzWriter`, fasterp)
//! are decoded member-parallel with libdeflate, falling back to the
//! sequential decoder at the first member without a size or the first batch
//! that fails. The result must be exactly flate2's, accepted or rejected.
#![no_main]
#[path = "common.rs"]
mod common;
use libfuzzer_sys::fuzz_target;
use std::io::Write;

fuzz_target!(|data: &[u8]| {
    let mut c = common::Choices::new(data);
    let mut gz = Vec::new();
    for _ in 0..1 + c.below(8) {
        let part = common::payload(&mut c, 40_000);
        match c.byte() % 3 {
            // Sized members, several per call.
            0 => {
                let options = readkit_gz::EncodeOptions {
                    level: u32::from(c.byte() % 13),
                    threads: c.below(3),
                    block_size: 1 + c.below(20_000),
                    ..Default::default()
                };
                let mut w = readkit_gz::GzWriter::new(&mut gz, &options).unwrap();
                w.write_all(&part).unwrap();
                w.finish().unwrap();
            }
            1 => {
                let mut out = Vec::new();
                readkit_gz::fastq_gzip_member(&part, &mut out);
                gz.extend(out);
            }
            // An unsized member: the rest is sequential.
            _ => gz.extend(common::gzip(&mut c, &part)),
        }
    }
    if c.byte() % 2 == 0 {
        for _ in 0..1 + c.below(3) {
            if gz.is_empty() {
                break;
            }
            let at = c.below(gz.len());
            gz[at] ^= 1 << (c.byte() % 8);
        }
    }
    if c.byte() % 8 == 0 {
        let keep = c.below(gz.len() + 1);
        gz.truncate(keep);
    }
    if let Ok(path) = std::env::var("FUZZ_DUMP") {
        std::fs::write(path, &gz).unwrap();
    }
    let options = readkit_gz::DecodeOptions {
        threads: Some(1 + c.below(4)),
        in_flight: Some(1 + c.below(4)),
        buffer_size: 1 + c.below(100_000),
        ..Default::default()
    };
    let ours = common::parallel(&gz, &options);
    let theirs = common::flate2_decode(&gz);
    match (ours, theirs) {
        (Ok(a), Ok(b)) => assert!(a == b, "different bytes than flate2"),
        (Err(_), Err(_)) => {}
        (Ok(_), Err(e)) => panic!("accepted a stream flate2 rejects: {e}"),
        (Err(e), Ok(_)) => panic!("rejected a stream flate2 accepts: {e}"),
    }
});
