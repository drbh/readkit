//! Remote objects, checksums, chunk planning and verification.

use crate::error::{Error, ErrorKind};
use md5::{Digest as _, Md5};
use serde::{Deserialize, Serialize};

/// A remote object to download: where it is, how big it is, and how to verify
/// it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Object {
    /// HTTP(S) URL; the server must honour `Range` requests.
    pub url: String,
    /// Size in bytes.
    pub size: u64,
    /// Publisher checksum the download is verified against.
    pub checksum: Checksum,
}

impl Object {
    /// An object at `url` of `size` bytes, verified with `checksum`.
    pub fn new(url: impl Into<String>, size: u64, checksum: Checksum) -> Self {
        Self {
            url: url.into(),
            size,
            checksum,
        }
    }
}

/// How a download is verified.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Checksum {
    /// No verification.
    None,
    /// MD5 of the whole object (for example ENA's `fastq_md5`).
    Md5 {
        /// Lowercase or uppercase hex digest.
        hex: String,
    },
    /// An S3 ETag. For multipart uploads (`parts > 1`) it is the MD5 of the
    /// concatenated part MD5s followed by `-<parts>`; otherwise it is the MD5
    /// of the object.
    S3MultipartEtag {
        /// The ETag, without quotes.
        etag: String,
        /// Number of parts, from the ETag's `-N` suffix (1 if absent).
        parts: u32,
        /// Size of every part but the last.
        part_size: u64,
    },
}

impl Checksum {
    /// Build an S3 checksum from an ETag (quotes are stripped), taking the
    /// part count from its `-N` suffix. `part_size` is ignored for
    /// single-part ETags.
    pub fn from_s3_etag(etag: &str, part_size: u64) -> Result<Self, Error> {
        let etag = etag.trim().trim_matches('"').to_string();
        let parts = match etag.rsplit_once('-') {
            Some((_, n)) => n.parse::<u32>().ok().filter(|n| *n > 0).ok_or_else(|| {
                Error::new(
                    ErrorKind::InvalidInput,
                    format!("malformed multipart ETag {etag:?}"),
                )
            })?,
            None => 1,
        };
        Ok(Self::S3MultipartEtag {
            etag,
            parts,
            part_size,
        })
    }

    pub(crate) fn verifies(&self) -> bool {
        !matches!(self, Self::None)
    }
}

/// Serialized form of [`Checksum`] used in the ledger and the verified
/// marker. Externally tagged: `"none"`, `{"md5":{"hex":..}}`,
/// `{"s3_multipart_etag":{"etag":..,"parts":..,"part_size":..}}`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum ChecksumRepr {
    None,
    Md5 {
        hex: String,
    },
    S3MultipartEtag {
        etag: String,
        parts: u32,
        part_size: u64,
    },
}

impl From<&Checksum> for ChecksumRepr {
    fn from(c: &Checksum) -> Self {
        match c {
            Checksum::None => Self::None,
            Checksum::Md5 { hex } => Self::Md5 { hex: hex.clone() },
            Checksum::S3MultipartEtag {
                etag,
                parts,
                part_size,
            } => Self::S3MultipartEtag {
                etag: etag.clone(),
                parts: *parts,
                part_size: *part_size,
            },
        }
    }
}

/// A transfer chunk.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ChunkSpec {
    pub offset: u64,
    pub len: u64,
    pub unit: usize,
}

/// A verification unit: one S3 part, or the whole object.
#[derive(Debug, Clone, Copy)]
pub(crate) struct UnitSpec {
    pub offset: u64,
    pub len: u64,
    pub first_chunk: usize,
    pub chunks: usize,
}

#[derive(Debug, Clone)]
pub(crate) struct Layout {
    pub units: Vec<UnitSpec>,
    pub chunks: Vec<ChunkSpec>,
    pub chunk_bytes: u64,
}

pub(crate) fn hex(bytes: &[u8]) -> String {
    use std::fmt::Write as _;
    let mut s = String::with_capacity(bytes.len() * 2);
    for b in bytes {
        let _ = write!(s, "{b:02x}");
    }
    s
}

pub(crate) fn unhex32(h: &str) -> Option<[u8; 32]> {
    if h.len() != 64 || !h.is_ascii() {
        return None;
    }
    let mut out = [0u8; 32];
    for (i, b) in out.iter_mut().enumerate() {
        *b = u8::from_str_radix(&h[2 * i..2 * i + 2], 16).ok()?;
    }
    Some(out)
}

/// Largest chunk of an object verified by one whole-object MD5.
const WHOLE_MD5_CHUNK: u64 = 4 << 20;

/// Split `object` into verification units and transfer chunks of at most
/// `chunk_bytes`, and small enough that there are at least `min_chunks` of
/// them (with a 1 MiB floor). Chunks never cross units.
pub(crate) fn plan(object: &Object, chunk_bytes: u64, min_chunks: u64) -> Result<Layout, Error> {
    let chunk_bytes = chunk_bytes.max(1);
    let mut units = Vec::new();
    match &object.checksum {
        Checksum::S3MultipartEtag {
            parts, part_size, ..
        } if *parts > 1 => {
            let parts = *parts as u64;
            let last = part_size
                .checked_mul(parts - 1)
                .and_then(|full| object.size.checked_sub(full))
                .filter(|l| *l > 0 && l <= part_size)
                .ok_or_else(|| {
                    Error::new(
                        ErrorKind::InvalidInput,
                        format!(
                            "{}: inconsistent S3 part layout ({parts} parts of {part_size} for {} bytes)",
                            object.url, object.size
                        ),
                    )
                })?;
            for i in 0..parts {
                let len = if i + 1 == parts { last } else { *part_size };
                units.push((i * part_size, len));
            }
        }
        Checksum::S3MultipartEtag { parts: 0, .. } => {
            return Err(Error::new(
                ErrorKind::InvalidInput,
                format!("{}: S3 ETag with zero parts", object.url),
            ));
        }
        _ => units.push((0, object.size)),
    }
    // At least `min_chunks` chunks, so that a small object still uses every
    // connection on hosts that cap per-connection speed. S3 parts are split
    // too: a 64 MiB object of four 16 MiB parts would otherwise use four
    // connections.
    let mut chunk_bytes = chunk_bytes.min(object.size.div_ceil(min_chunks.max(1)).max(1 << 20));
    // One MD5 over the whole object is hashed in order as chunks land. The
    // last round of chunks lands together and is hashed after the transfer
    // ends, so keep that round small.
    if units.len() == 1 && object.checksum.verifies() {
        chunk_bytes = chunk_bytes.min(WHOLE_MD5_CHUNK);
    }
    let mut chunks = Vec::new();
    let mut specs = Vec::new();
    for (u, (offset, len)) in units.into_iter().enumerate() {
        let first = chunks.len();
        let mut at = offset;
        while at < offset + len {
            let l = chunk_bytes.min(offset + len - at);
            chunks.push(ChunkSpec {
                offset: at,
                len: l,
                unit: u,
            });
            at += l;
        }
        if len == 0 {
            chunks.push(ChunkSpec {
                offset,
                len: 0,
                unit: u,
            });
        }
        specs.push(UnitSpec {
            offset,
            len,
            first_chunk: first,
            chunks: chunks.len() - first,
        });
    }
    Ok(Layout {
        units: specs,
        chunks,
        chunk_bytes,
    })
}

/// Check unit digests against the object's checksum. Per-part MD5s are unknown
/// for S3 (only the combined ETag), so units verify collectively.
pub(crate) fn verify(object: &Object, unit_md5: &[[u8; 16]]) -> Result<(), Error> {
    let mismatch = |what: &str, got: &str, want: &str| {
        Err(Error::new(
            ErrorKind::Verification,
            format!("{}: {what} mismatch (got {got}, want {want})", object.url),
        ))
    };
    match &object.checksum {
        Checksum::None => Ok(()),
        Checksum::Md5 { hex: want } => {
            let got = hex(&unit_md5[0]);
            if got.eq_ignore_ascii_case(want) {
                Ok(())
            } else {
                mismatch("MD5", &got, want)
            }
        }
        Checksum::S3MultipartEtag { etag, parts, .. } => {
            let got = if *parts > 1 {
                let mut h = Md5::new();
                for d in unit_md5 {
                    h.update(d);
                }
                format!("{}-{}", hex(&h.finalize()), unit_md5.len())
            } else {
                hex(&unit_md5[0])
            };
            if got.eq_ignore_ascii_case(etag) {
                Ok(())
            } else {
                mismatch("S3 ETag", &got, etag)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn md5_obj(size: u64) -> Object {
        Object::new("http://h/o", size, Checksum::Md5 { hex: String::new() })
    }

    #[test]
    fn small_objects_get_at_least_min_chunks() {
        let l = plan(&md5_obj(64 << 20), 8 << 20, 16).unwrap();
        assert_eq!(l.chunks.len(), 16);
        assert_eq!(l.units.len(), 1);
        let l = plan(&md5_obj(64 << 20), 8 << 20, 32).unwrap();
        assert_eq!(l.chunks.len(), 32);
        // Large objects use `chunk_bytes`, at most 4 MiB for a whole-object MD5.
        let l = plan(&md5_obj(1 << 30), 16 << 20, 32).unwrap();
        assert_eq!(l.chunks.len(), 256);
        let none = Object::new("u", 1 << 30, Checksum::None);
        assert_eq!(plan(&none, 16 << 20, 32).unwrap().chunks.len(), 64);
        // Small files keep a 1 MiB floor.
        let l = plan(&md5_obj(3 << 20), 8 << 20, 16).unwrap();
        assert_eq!(l.chunks.len(), 3);
    }

    #[test]
    fn s3_parts_are_split_to_reach_min_chunks() {
        let part = 16 << 20;
        let o = Object::new(
            "u",
            4 * part,
            Checksum::S3MultipartEtag {
                etag: "x-4".into(),
                parts: 4,
                part_size: part,
            },
        );
        let l = plan(&o, 16 << 20, 32).unwrap();
        assert_eq!(l.units.len(), 4);
        assert_eq!(l.chunks.len(), 32);
        let l = plan(&o, 16 << 20, 4).unwrap();
        assert_eq!(l.chunks.len(), 4);
    }

    #[test]
    fn s3_parts_become_units_and_chunks_do_not_cross_them() {
        let part = 5 << 20;
        let size = 2 * part + 7;
        let o = Object::new(
            "u",
            size,
            Checksum::S3MultipartEtag {
                etag: "x-3".into(),
                parts: 3,
                part_size: part,
            },
        );
        let l = plan(&o, 2 << 20, 1).unwrap();
        assert_eq!(l.units.len(), 3);
        assert_eq!(l.units[2].len, 7);
        for c in &l.chunks {
            let u = l.units[c.unit];
            assert!(c.offset >= u.offset && c.offset + c.len <= u.offset + u.len);
        }
        assert_eq!(l.chunks.iter().map(|c| c.len).sum::<u64>(), size);
    }

    #[test]
    fn inconsistent_s3_layout_is_invalid_input() {
        let o = Object::new(
            "u",
            10,
            Checksum::S3MultipartEtag {
                etag: "x-3".into(),
                parts: 3,
                part_size: 10,
            },
        );
        assert_eq!(plan(&o, 1, 1).unwrap_err().kind(), ErrorKind::InvalidInput);
    }

    #[test]
    fn etag_suffix_gives_part_count() {
        let c = Checksum::from_s3_etag("\"abc-12\"", 8).unwrap();
        assert_eq!(
            c,
            Checksum::S3MultipartEtag {
                etag: "abc-12".into(),
                parts: 12,
                part_size: 8
            }
        );
        let c = Checksum::from_s3_etag("abc", 8).unwrap();
        assert!(matches!(c, Checksum::S3MultipartEtag { parts: 1, .. }));
        assert_eq!(
            Checksum::from_s3_etag("abc-x", 8).unwrap_err().kind(),
            ErrorKind::InvalidInput
        );
    }

    #[test]
    fn empty_object_has_one_empty_chunk() {
        let l = plan(&md5_obj(0), 8 << 20, 16).unwrap();
        assert_eq!(l.chunks.len(), 1);
        assert_eq!(l.chunks[0].len, 0);
    }
}
