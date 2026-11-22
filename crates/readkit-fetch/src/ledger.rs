//! Resume ledger: which chunks are on the sink, with their BLAKE3 digests.
//!
//! JSON: `{"version":1,"url":..,"size":..,"checksum":..,"chunk_size":..,
//! "done":"0110..","blake3":[null,"<hex>",..]}`. A ledger whose url, size,
//! checksum, chunk size or chunk count differ from the current download is
//! ignored.

use crate::object::{ChecksumRepr, Layout, Object};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
struct LedgerFile {
    version: u32,
    url: String,
    size: u64,
    checksum: ChecksumRepr,
    chunk_size: u64,
    done: String,
    blake3: Vec<Option<String>>,
}

pub(crate) struct Resume {
    pub done: Vec<bool>,
    pub blake3: Vec<Option<[u8; 32]>>,
}

pub(crate) fn encode(
    object: &Object,
    layout: &Layout,
    done: impl Iterator<Item = bool>,
    blake3: impl Iterator<Item = Option<[u8; 32]>>,
) -> Vec<u8> {
    let file = LedgerFile {
        version: 1,
        url: object.url.clone(),
        size: object.size,
        checksum: ChecksumRepr::from(&object.checksum),
        chunk_size: layout.chunk_bytes,
        done: done.map(|d| if d { '1' } else { '0' }).collect(),
        blake3: blake3.map(|d| d.map(|d| crate::object::hex(&d))).collect(),
    };
    serde_json::to_vec(&file).unwrap_or_default()
}

pub(crate) fn parse(bytes: &[u8], object: &Object, layout: &Layout) -> Option<Resume> {
    let file: LedgerFile = serde_json::from_slice(bytes).ok()?;
    let n = layout.chunks.len();
    if file.version != 1
        || file.url != object.url
        || file.size != object.size
        || file.checksum != ChecksumRepr::from(&object.checksum)
        || file.chunk_size != layout.chunk_bytes
        || file.done.len() != n
    {
        return None;
    }
    let done = file.done.bytes().map(|b| b == b'1').collect();
    let blake3 = (0..n)
        .map(|i| {
            file.blake3
                .get(i)
                .cloned()
                .flatten()
                .and_then(|h| crate::object::unhex32(&h))
        })
        .collect();
    Some(Resume { done, blake3 })
}
