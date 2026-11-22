//! The facade's re-exports, used the way a tool would under the default
//! features.

use std::io::Write;

use readkit::fastq::{FastqInput, FastqOptions, FastqSink, FastqSource};
use readkit::gz::{DecodeOptions, EncodeOptions, GzReader, GzWriter, ParGzReader};
use readkit::io::VecSink;
use readkit::{BatchBuilder, BatchSink, BatchSource, Paired, Pool, Read, ReadAhead, Single, pump};

const R1: &[u8] = b"@a/1\nACGT\n+\nIIII\n@b/1\nGG\n+\n##\n";
const R2: &[u8] = b"@a/2\nTTGG\n+\nIIII\n@b/2\nC\n+\n#\n";

#[test]
fn paired_fastq_copies_through_batches() {
    let source = FastqSource::paired(
        FastqInput::new(R1),
        FastqInput::new(R2),
        FastqOptions::default(),
    );
    let mut sink = FastqSink::split(vec![Vec::new(), Vec::new()]);
    let records = pump(ReadAhead::new(source), &mut sink).unwrap();
    assert_eq!(records, 2);
    let out = sink.into_inner().unwrap();
    assert_eq!(out[0], R1);
    assert_eq!(out[1], R2);
}

#[test]
fn gzip_feeds_the_parser() {
    let mut writer = GzWriter::new(Vec::new(), &EncodeOptions::default()).unwrap();
    writer.write_all(R1).unwrap();
    let gz = writer.finish().unwrap();

    // Sequential and parallel decoding give the parser the same text.
    let len = gz.len() as u64;
    let inputs = [
        FastqInput::new(std::io::BufReader::new(GzReader::from_reader(
            std::io::Cursor::new(gz.clone()),
        ))),
        FastqInput::new(ParGzReader::new(gz, len, &DecodeOptions::default()).unwrap()),
    ];
    for input in inputs {
        let mut sink = VecSink::new();
        let source = FastqSource::single(input, FastqOptions::default());
        assert_eq!(pump(source, &mut sink).unwrap(), 2);
        let batch = &sink.batches()[0];
        let seqs = batch.column().seqs();
        assert_eq!(seqs.offsets(), &[0, 4, 6]);
        assert_eq!(seqs.data(), b"ACGTGG");
        assert_eq!(batch.column().names().unwrap().get(1), b"b/1");
    }
}

#[test]
fn built_batches_write_as_fastq() {
    let pool = Pool::new();
    let mut single = BatchBuilder::<Single>::new(&pool);
    single
        .push_read(Read::named(b"a", b"ACGT", b"IIII"))
        .unwrap();
    let mut sink = FastqSink::new(Vec::new());
    sink.accept(single.finish()).unwrap();
    assert_eq!(sink.into_inner().unwrap()[0], b"@a\nACGT\n+\nIIII\n");

    let mut paired = BatchBuilder::<Paired>::new(&pool);
    paired
        .push_pair(
            Read::named(b"a/1", b"AC", b"II"),
            Read::named(b"a/2", b"GT", b"II"),
        )
        .unwrap();
    let batch = paired.finish();
    assert_eq!(batch.r2().seqs().data(), b"GT");

    // A source is anything that yields batches.
    let mut source = readkit::io::iter_source([Ok::<_, std::convert::Infallible>(batch)]);
    assert_eq!(source.next_batch().unwrap().unwrap().len(), 1);
    assert!(
        BatchSource::<Paired>::next_batch(&mut source)
            .unwrap()
            .is_none()
    );
}
