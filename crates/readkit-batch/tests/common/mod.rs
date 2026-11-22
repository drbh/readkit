//! Minimal FASTQ helpers for tests: a 4-line record splitter and writer.
//! The real parser and writer belong to readkit-fastq.
#![allow(dead_code)]

use readkit_batch::{Batch, Layout, Read};

/// Splits 4-line FASTQ text into `(name, seq, qual)` without copying.
pub fn records(text: &[u8]) -> Vec<Read<'_>> {
    let mut lines = text.split(|&b| b == b'\n');
    let mut out = Vec::new();
    while let Some(h) = lines.next() {
        if h.is_empty() {
            break; // trailing newline
        }
        let (seq, plus, qual) = (
            lines.next().unwrap(),
            lines.next().unwrap(),
            lines.next().unwrap(),
        );
        assert_eq!(h[0], b'@');
        assert_eq!(plus, b"+");
        out.push(Read::named(&h[1..], seq, qual));
    }
    out
}

/// Writes segment `seg` of a named batch as FASTQ text.
pub fn write_fastq<L: Layout>(batch: &Batch<L>, seg: usize, out: &mut Vec<u8>) {
    let col = batch.segment(seg);
    let names = col.names().expect("batch has names");
    let (seqs, quals) = (col.seqs(), col.quals());
    for i in 0..col.len() {
        out.push(b'@');
        out.extend_from_slice(names.get(i));
        out.push(b'\n');
        out.extend_from_slice(seqs.get(i));
        out.extend_from_slice(b"\n+\n");
        out.extend_from_slice(quals.get(i));
        out.push(b'\n');
    }
}

pub const R1: &[u8] = b"@SRR1.1 1/1\nACGTNacgtRYKM\n+\nIIII#IIIIIIII\n\
@SRR1.2 2/1\n\n+\n\n\
@SRR1.3 3/1\nGATTACA\n+\n!\"#$%&'\n";
pub const R2: &[u8] = b"@SRR1.1 1/2\nTTTT\n+\nABCD\n\
@SRR1.2 2/2\nCCCCCCCCCCCC\n+\nFFFFFFFFFFFF\n\
@SRR1.3 3/2\nN\n+\n~\n";
