//! Arrow export (feature `arrow`): zero-copy, and buffers return to the pool
//! only after the last Arrow reference drops.
#![cfg(feature = "arrow")]

use readkit_batch::arrow::arrow_array::{Array, BinaryArray, StructArray, UInt16Array};
use readkit_batch::arrow::arrow_schema::DataType;
use readkit_batch::{BatchBuilder, Paired, Pool, Read, Segments};

fn paired(pool: &Pool) -> readkit_batch::Batch<Paired> {
    let mut b = BatchBuilder::<Paired>::with_pool(pool, 3);
    b.push_pair(
        Read::named(b"r0", b"ACGT", b"IIII"),
        Read::named(b"r0", b"GG", b"FF"),
    )
    .unwrap();
    b.push_pair(Read::named(b"r1", b"", b""), Read::named(b"r1", b"T", b"#"))
        .unwrap();
    b.set_tag(7);
    b.finish()
}

fn field<'a>(s: &'a StructArray, name: &str) -> &'a BinaryArray {
    s.column_by_name(name)
        .unwrap()
        .as_any()
        .downcast_ref::<BinaryArray>()
        .unwrap()
}

#[test]
fn export_is_zero_copy() {
    let pool = Pool::new();
    let batch = paired(&pool);
    let (seq_ptr, qual_ptr, off_ptr, name_ptr) = (
        batch.r2().seqs().values().as_ptr(),
        batch.r2().quals().values().as_ptr(),
        batch.r2().seqs().offsets().as_ptr(),
        batch.r1().names().unwrap().values().as_ptr(),
    );
    let tag_ptr = batch.tags().unwrap().as_ptr();
    let rb = batch.into_arrow();
    assert_eq!(rb.num_rows(), 2);
    assert_eq!(rb.num_columns(), 3);
    let schema = rb.schema();
    assert_eq!(schema.field(0).name(), "seg0");
    assert_eq!(schema.field(2).data_type(), &DataType::UInt16);

    let s1 = rb.column(1).as_any().downcast_ref::<StructArray>().unwrap();
    let (seq, qual) = (field(s1, "seq"), field(s1, "qual"));
    assert_eq!(seq.value(0), b"GG");
    assert_eq!(qual.value(1), b"#");
    assert_eq!(seq.values().as_ptr(), seq_ptr);
    assert_eq!(qual.values().as_ptr(), qual_ptr);
    assert_eq!(seq.value_offsets().as_ptr(), off_ptr);
    assert_eq!(qual.value_offsets().as_ptr(), off_ptr); // shared offsets

    let s0 = rb.column(0).as_any().downcast_ref::<StructArray>().unwrap();
    let names = field(s0, "name");
    assert_eq!(names.values().as_ptr(), name_ptr);
    assert_eq!(names.null_count(), 0);
    assert_eq!(field(s0, "seq").value(1), b"");

    let tags = rb.column(2).as_any().downcast_ref::<UInt16Array>().unwrap();
    assert_eq!(tags.values().as_ptr(), tag_ptr);
    assert_eq!(tags.values().as_ref(), &[0, 7]);
}

#[test]
fn buffers_return_after_last_arrow_reference() {
    let pool = Pool::new();
    let batch = paired(&pool);
    let buffers = pool.stats().fresh as usize; // 2 x (offs, seq, qual, name offs, names) + tags
    assert_eq!(buffers, 11);
    let rb = batch.into_arrow();
    assert_eq!(pool.stats().idle, 0);

    let clone = rb.clone();
    let slice = rb.slice(1, 1);
    let seq_col = rb.column(0).clone();
    drop(rb);
    assert_eq!(pool.stats().idle, 0);
    drop(clone);
    assert_eq!(pool.stats().idle, 0);
    assert_eq!(slice.num_rows(), 1);
    drop(slice);
    // only buffers referenced by the seg0 column are still out
    assert_eq!(pool.stats().idle, 11 - 5);
    drop(seq_col);
    assert_eq!(pool.stats().idle, 11);

    // and get reused by the next batch
    let before = pool.stats();
    drop(paired(&pool));
    assert_eq!(pool.stats().fresh, before.fresh);
}

#[test]
fn unnamed_batch_exports_null_names() {
    let pool = Pool::new();
    let mut b = BatchBuilder::<Segments>::with_segments(&pool, 3, 1);
    b.push_with(&[28, 90, 10], None, |_, s, q| {
        s.fill(b'A');
        q.fill(b'I');
    })
    .unwrap();
    let batch = b.finish();
    assert_eq!(batch.arrow_schema(), readkit_batch::arrow::schema(3, false));
    let rb = batch.into_arrow();
    let s = rb.column(2).as_any().downcast_ref::<StructArray>().unwrap();
    assert_eq!(field(s, "name").null_count(), 1);
    assert_eq!(field(s, "seq").value(0).len(), 10);
}
