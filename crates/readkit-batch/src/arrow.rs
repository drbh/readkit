//! Zero-copy export to arrow-rs (feature `arrow`).
//!
//! [`Batch::into_arrow`] wraps each pooled buffer with
//! `Buffer::from_custom_allocation`; the owner is the pooled buffer itself,
//! so it returns to its [`Pool`](crate::Pool) when the last Arrow reference
//! (including clones, slices and C Data Interface exports) drops.
//!
//! Schema: one non-null struct column per segment slot, named `seg0`,
//! `seg1`, ..., with fields `name: Binary` (nullable; all null when the batch
//! has no names), `seq: Binary` and `qual: Binary`; plus `tag: UInt16` when
//! the batch has tags. `seq` and `qual` share one offsets buffer.
//!
//! This module targets arrow-rs 60; the crates are re-exported so callers
//! can match versions.

use std::ptr::NonNull;
use std::sync::Arc;

pub use arrow_array;
pub use arrow_buffer;
pub use arrow_schema;

use arrow_array::{ArrayRef, BinaryArray, RecordBatch, StructArray, UInt16Array};
use arrow_buffer::alloc::Allocation;
use arrow_buffer::{Buffer, OffsetBuffer, ScalarBuffer};
use arrow_schema::{DataType, Field, Fields, Schema};

use crate::Layout;
use crate::batch::Batch;
use crate::pool::{Element, PooledVec};

/// Fields of one segment slot's struct column.
pub fn segment_fields() -> Fields {
    Fields::from(vec![
        Field::new("name", DataType::Binary, true),
        Field::new("seq", DataType::Binary, false),
        Field::new("qual", DataType::Binary, false),
    ])
}

/// Name of segment slot `i`'s column (`seg0`, `seg1`, ...).
pub fn segment_column_name(i: usize) -> String {
    format!("seg{i}")
}

/// The schema of a batch with `segments` slots, with or without tags.
pub fn schema(segments: usize, tags: bool) -> Schema {
    let mut fields: Vec<Field> = (0..segments)
        .map(|i| {
            Field::new(
                segment_column_name(i),
                DataType::Struct(segment_fields()),
                false,
            )
        })
        .collect();
    if tags {
        fields.push(Field::new("tag", DataType::UInt16, false));
    }
    Schema::new(fields)
}

fn buffer<T: Element>(v: PooledVec<T>) -> Buffer {
    let len = std::mem::size_of_val(v.as_slice());
    // A Vec's pointer is never null (dangling and aligned when empty).
    let ptr = NonNull::new(v.as_ptr() as *mut u8).expect("Vec pointer is non-null");
    let owner: Arc<dyn Allocation> = Arc::new(v);
    // SAFETY: `ptr` points to `len` initialized bytes owned by the Vec inside
    // `owner`. Moving the Vec does not move its heap allocation, and nothing
    // can mutate or free it until the last `Arc` to `owner` drops.
    unsafe { Buffer::from_custom_allocation(ptr, len, owner) }
}

fn offsets(v: PooledVec<i32>) -> OffsetBuffer<i32> {
    let n = v.len();
    OffsetBuffer::new(ScalarBuffer::new(buffer(v), 0, n))
}

impl<L: Layout> Batch<L> {
    /// The Arrow schema [`into_arrow`](Self::into_arrow) produces.
    pub fn arrow_schema(&self) -> Schema {
        schema(self.segments(), self.tags().is_some())
    }

    /// Converts into an Arrow `RecordBatch` without copying names,
    /// sequences, qualities, offsets or tags. Buffers return to their pools
    /// when Arrow releases them.
    pub fn into_arrow(self) -> RecordBatch {
        let schema = Arc::new(self.arrow_schema());
        let n = self.len;
        let mut columns: Vec<ArrayRef> = Vec::with_capacity(self.segs.len() + 1);
        for seg in self.segs {
            let names = match seg.names {
                Some((o, v)) => BinaryArray::new(offsets(o), buffer(v), None),
                None => BinaryArray::new_null(n),
            };
            let offs = offsets(seg.offsets);
            let seqs = BinaryArray::new(offs.clone(), buffer(seg.seqs), None);
            let quals = BinaryArray::new(offs, buffer(seg.quals), None);
            columns.push(Arc::new(StructArray::new(
                segment_fields(),
                vec![Arc::new(names), Arc::new(seqs), Arc::new(quals)],
                None,
            )));
        }
        if let Some(t) = self.tags {
            let len = t.len();
            columns.push(Arc::new(UInt16Array::new(
                ScalarBuffer::new(buffer(t), 0, len),
                None,
            )));
        }
        RecordBatch::try_new(schema, columns).expect("columns match the schema by construction")
    }
}
