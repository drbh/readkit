//! Consumer sketches: a bulk k-mer consumer and a FASTQ writer.

mod common;

use readkit_batch::{BatchBuilder, Binary, Paired, Pool};

fn paired_batch(pool: &Pool) -> readkit_batch::Batch<Paired> {
    let (r1, r2) = (common::records(common::R1), common::records(common::R2));
    let mut b = BatchBuilder::<Paired>::with_pool(pool, r1.len());
    for (a, m) in r1.iter().zip(&r2) {
        b.push_pair(*a, *m).unwrap();
    }
    b.finish()
}

/// k-mers per read computed from offsets alone, and k-mer hashing over one
/// contiguous buffer: no per-record views or allocations.
fn count_kmers(seqs: Binary<'_>, k: usize) -> (usize, u64) {
    let (offsets, values) = (seqs.offsets(), seqs.values());
    let mut n = 0;
    let mut h = 0u64;
    for w in offsets.windows(2) {
        let (a, b) = (w[0] as usize, w[1] as usize);
        if b - a < k {
            continue;
        }
        n += b - a - k + 1;
        for kmer in values[a..b].windows(k) {
            h = h
                .wrapping_mul(31)
                .wrapping_add(kmer.iter().map(|&x| u64::from(x)).sum::<u64>());
        }
    }
    (n, h)
}

/// (d) A kallistors-like consumer: takes the sequence buffer and spans the
/// way `PackedSeqBatch::from_parts(buf, spans)` does, straight from the
/// column.
#[test]
fn bulk_consumer_reads_seq_column_via_offsets() {
    let pool = Pool::new();
    let batch = paired_batch(&pool);
    let seqs = batch.r1().seqs();
    let spans: Vec<std::ops::Range<u32>> = seqs
        .offsets()
        .windows(2)
        .map(|w| w[0] as u32..w[1] as u32)
        .collect();
    assert_eq!(spans, [0..13, 13..13, 13..20]);
    // the buffer is the column itself: no copy
    assert_eq!(seqs.values().as_ptr(), batch.r1().seqs().data().as_ptr());
    let (n, _) = count_kmers(seqs, 5);
    assert_eq!(n, (13 - 5 + 1) + (7 - 5 + 1));
    let (n2, _) = count_kmers(batch.r2().seqs(), 4);
    assert_eq!(n2, 1 + 9);
}

/// (e) FASTQ → batch → FASTQ is byte-identical for valid records.
#[test]
fn fastq_round_trip_is_byte_identical() {
    let pool = Pool::new();
    let batch = paired_batch(&pool);
    let (mut o1, mut o2) = (Vec::new(), Vec::new());
    common::write_fastq(&batch, 0, &mut o1);
    common::write_fastq(&batch, 1, &mut o2);
    assert_eq!(o1, common::R1);
    assert_eq!(o2, common::R2);
}
