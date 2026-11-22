//! The NCBI backend (feature `ngs`): the VDB C API from the SRA Toolkit's
//! `libncbi-ngs`, through a small C shim (`shim/vdb_shim.c`).
//!
//! Each reader owns one VDB manager reference, the run's `SEQUENCE` table
//! and one cursor over `READ_START`, `READ_LEN`, `READ_TYPE` and optionally
//! `NAME`, plus either:
//!
//! - the **fast path**: bases and qualities decoded here from the stored
//!   KDB blobs of the physical `READ` (2na packed), `ALTREAD` (ambiguity
//!   codes) and `QUALITY` (deflated phred) columns ([`crate::vdbfmt`]), or
//! - the **VDB path**: `READ` and `QUALITY` (phred+33 text) from the same
//!   cursor, served as slices into VDB's cursor cache.
//!
//! The fast path is used for unaligned runs with the table schema it was
//! checked against (`NCBI:align:tbl:seq#1.1`, as loaded from FASTQ by
//! `latf-load`) and original qualities, after the first spots decode
//! identically both ways (checked once per file and process). A blob it does
//! not recognise is served through VDB.
//! VDB spends most of its time converting 2na to text through 4na per row and
//! in its bundled zlib; the fast path unpacks with a table and inflates with
//! libdeflate, about 2× faster per thread.

use std::ffi::{CStr, CString, c_int, c_void};
use std::path::Path;

use crate::backend::{OpenRequest, RunInfo, Spot, SpotReader, SraBackend};
use crate::error::{Error, ErrorKind};
use crate::lite::QualityKind;
use crate::rc::rc_kind;
use crate::read_type::ReadType;
use crate::vdbfmt;

mod ffi {
    use std::ffi::{c_char, c_int, c_void};

    #[repr(C)]
    pub struct Spot {
        pub bases: *const c_char,
        pub quals: *const c_char,
        pub read_start: *const i32,
        pub read_len: *const u32,
        pub read_type: *const u8,
        pub name: *const c_char,
        pub bases_len: u32,
        pub quals_len: u32,
        pub nreads: u32,
        pub nstarts: u32,
        pub nlens: u32,
        pub name_len: u32,
    }

    #[repr(C)]
    #[derive(Default)]
    pub struct Info {
        pub first_row: i64,
        pub row_count: u64,
        pub quality: i32,
        pub is_database: i32,
        pub has_names: i32,
    }

    #[repr(C)]
    pub struct Status {
        pub rc: u32,
        pub step: *const c_char,
        pub message: [c_char; 512],
    }

    unsafe extern "C" {
        pub fn rk_open(
            path: *const c_char,
            want_names: c_int,
            want_seq: c_int,
            info: *mut Info,
            st: *mut Status,
        ) -> *mut c_void;
        pub fn rk_typespec(h: *mut c_void, buf: *mut c_char, size: usize, st: *mut Status)
        -> c_int;
        pub fn rk_kcol_open(h: *mut c_void, name: *const c_char, st: *mut Status) -> *mut c_void;
        pub fn rk_kcol_close(kc: *mut c_void);
        pub fn rk_kcol_blob(
            kc: *mut c_void,
            row: i64,
            first: *mut i64,
            count: *mut u32,
            buf: *mut u8,
            cap: usize,
            size: *mut usize,
            st: *mut Status,
        ) -> c_int;
        pub fn rk_read(h: *mut c_void, row: i64, out: *mut Spot, st: *mut Status) -> c_int;
        pub fn rk_close(h: *mut c_void);
    }
}

impl ffi::Status {
    fn new() -> Self {
        Self {
            rc: 0,
            step: std::ptr::null(),
            message: [0; 512],
        }
    }
}

fn status_error(st: &ffi::Status, what: &str, path: &Path) -> Error {
    let step = if st.step.is_null() {
        String::new()
    } else {
        // SAFETY: the shim sets `step` to a static NUL-terminated string.
        unsafe { CStr::from_ptr(st.step) }
            .to_string_lossy()
            .into_owned()
    };
    // SAFETY: the shim always NUL-terminates `message`.
    let msg = unsafe { CStr::from_ptr(st.message.as_ptr()) }.to_string_lossy();
    let (kind, text) = if st.rc == 0 {
        // The shim's own checks (unexpected element size).
        (
            ErrorKind::Unsupported,
            format!("{what}: unexpected column layout for {step}"),
        )
    } else {
        (rc_kind(st.rc), format!("{what} ({step}): {msg}"))
    };
    let mut e = Error::new(kind, text).with_input(path.display().to_string());
    if st.rc != 0 {
        e = e.with_code(st.rc);
    }
    e
}

/// The NCBI VDB backend. Each [`open`](SraBackend::open) makes a new reader.
///
/// [`new`](Self::new) uses the fast path where it applies (module docs);
/// [`vdb_only`](Self::vdb_only) always decodes through VDB, for comparison.
#[derive(Clone, Copy, Debug)]
pub struct NgsBackend {
    fast: bool,
}

impl Default for NgsBackend {
    fn default() -> Self {
        Self::new()
    }
}

impl NgsBackend {
    /// The backend, with the fast path.
    pub fn new() -> Self {
        Self { fast: true }
    }

    /// The backend without the fast path: every spot through VDB.
    pub fn vdb_only() -> Self {
        Self { fast: false }
    }
}

impl SraBackend for NgsBackend {
    fn name(&self) -> &'static str {
        "ncbi-vdb"
    }

    fn open(&self, path: &Path, request: &OpenRequest) -> Result<Box<dyn SpotReader>, Error> {
        let reader = if self.fast {
            NgsReader::open(path, request.names)?
        } else {
            NgsReader::open_vdb(path, request.names)?
        };
        Ok(Box::new(reader))
    }
}

/// A reader over one run. Not shared between threads.
pub struct NgsReader {
    inner: Inner,
}

enum Inner {
    Vdb(VdbReader),
    Fast(Box<FastReader>),
}

impl NgsReader {
    /// Open the run at `path`, with the fast path if it applies. `names`
    /// asks for the `NAME` column.
    pub fn open(path: &Path, names: bool) -> Result<Self, Error> {
        if let Some(fast) = FastReader::open(path, names)? {
            return Ok(Self {
                inner: Inner::Fast(Box::new(fast)),
            });
        }
        Self::open_vdb(path, names)
    }

    /// Open the run at `path`, decoding every spot through VDB.
    pub fn open_vdb(path: &Path, names: bool) -> Result<Self, Error> {
        Ok(Self {
            inner: Inner::Vdb(VdbReader::open(path, names, true)?),
        })
    }

    /// Whether bases and qualities come from the fast path.
    pub fn fast_path(&self) -> bool {
        matches!(self.inner, Inner::Fast(_))
    }

    /// Spots the fast path handed to VDB because their blobs were not in a
    /// form it decodes (0 without the fast path).
    pub fn vdb_fallback_spots(&self) -> u64 {
        match &self.inner {
            Inner::Fast(f) => f.fallback_spots,
            Inner::Vdb(_) => 0,
        }
    }

    /// Why the fast path handed the first such spot to VDB: the column and
    /// what was not recognised.
    pub fn fallback_reason(&self) -> Option<String> {
        let Inner::Fast(f) = &self.inner else {
            return None;
        };
        let cols = [
            ("READ", Some(&f.read)),
            ("QUALITY", Some(&f.quality)),
            ("ALTREAD", f.alt.as_ref()),
        ];
        cols.iter()
            .find_map(|(n, c)| c.and_then(|c| c.reason).map(|r| format!("{n}: {r}")))
            .or_else(|| f.mismatch.map(str::to_string))
    }
}

impl SpotReader for NgsReader {
    fn info(&self) -> &RunInfo {
        match &self.inner {
            Inner::Vdb(r) => &r.info,
            Inner::Fast(r) => &r.small.info,
        }
    }

    fn read_spot(&mut self, spot: u64) -> Result<Spot<'_>, Error> {
        match &mut self.inner {
            Inner::Vdb(r) => r.read_spot(spot),
            Inner::Fast(r) => r.read_spot(spot),
        }
    }
}

impl std::fmt::Debug for NgsReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (path, info) = match &self.inner {
            Inner::Vdb(r) => (&r.path, &r.info),
            Inner::Fast(r) => (&r.small.path, &r.small.info),
        };
        f.debug_struct("NgsReader")
            .field("path", path)
            .field("info", info)
            .field("fast_path", &self.fast_path())
            .finish()
    }
}

/// One shim handle: a VDB cursor over the run.
struct VdbReader {
    handle: *mut c_void,
    info: RunInfo,
    path: std::path::PathBuf,
}

// SAFETY: a VDB cursor may be used from any thread as long as only one
// thread uses it at a time; readers are `!Sync` and need `&mut self`.
unsafe impl Send for VdbReader {}

impl VdbReader {
    /// Open the run at `path`. `names` asks for the `NAME` column; `seq` for
    /// `READ` and `QUALITY` (without them spots have empty bases and
    /// qualities).
    fn open(path: &Path, names: bool, seq: bool) -> Result<Self, Error> {
        let c_path = CString::new(path.as_os_str().as_encoded_bytes()).map_err(|_| {
            Error::new(ErrorKind::Unsupported, "path contains a NUL byte")
                .with_input(path.display().to_string())
        })?;
        let mut info = ffi::Info::default();
        let mut st = ffi::Status::new();
        // SAFETY: valid NUL-terminated path and out-pointers.
        let handle = unsafe {
            ffi::rk_open(
                c_path.as_ptr(),
                c_int::from(names),
                c_int::from(seq),
                &mut info,
                &mut st,
            )
        };
        if handle.is_null() {
            return Err(status_error(&st, "opening the run", path));
        }
        let mut run = RunInfo::new(info.first_row.max(0) as u64, info.row_count);
        run.quality = match info.quality {
            1 => QualityKind::Original,
            2 => QualityKind::Synthetic,
            _ => QualityKind::Unknown,
        };
        run.has_names = info.has_names != 0;
        run.is_database = info.is_database != 0;
        let reader = Self {
            handle,
            info: run,
            path: path.to_path_buf(),
        };
        if info.first_row < 0 {
            return Err(Error::new(
                ErrorKind::Unsupported,
                format!("negative first row id {}", info.first_row),
            )
            .with_input(path.display().to_string()));
        }
        Ok(reader)
    }

    /// The raw spot from the shim; its pointers stay valid until the next
    /// read on this handle.
    fn read_raw(&mut self, spot: u64) -> Result<ffi::Spot, Error> {
        let row = i64::try_from(spot)
            .map_err(|_| Error::new(ErrorKind::NotFound, "spot id out of range").at_spot(spot))?;
        let mut out = std::mem::MaybeUninit::<ffi::Spot>::uninit();
        let mut st = ffi::Status::new();
        // SAFETY: `handle` is live; `out` is written by the shim (memset
        // first) before it returns.
        let rc = unsafe { ffi::rk_read(self.handle, row, out.as_mut_ptr(), &mut st) };
        if rc != 0 {
            return Err(status_error(&st, "reading a spot", &self.path).at_spot(spot));
        }
        // SAFETY: the shim initialized `out` (it zeroes it before filling).
        Ok(unsafe { out.assume_init() })
    }

    fn read_spot(&mut self, spot: u64) -> Result<Spot<'_>, Error> {
        let s = self.read_raw(spot)?;
        // SAFETY: VDB keeps the cell data alive in the cursor cache until the
        // next read on this cursor, which needs `&mut self` again; the
        // returned slices borrow `self`.
        unsafe { Ok(spot_from_raw(&s, None)) }
    }

    fn typespec(&self) -> Result<String, Error> {
        let mut buf = [0 as std::ffi::c_char; 256];
        let mut st = ffi::Status::new();
        // SAFETY: live handle, buffer of the size passed.
        if unsafe { ffi::rk_typespec(self.handle, buf.as_mut_ptr(), buf.len(), &mut st) } != 0 {
            return Err(status_error(&st, "reading the table type", &self.path));
        }
        // SAFETY: VTableTypespec NUL-terminates within the buffer.
        Ok(unsafe { CStr::from_ptr(buf.as_ptr()) }
            .to_string_lossy()
            .into_owned())
    }
}

impl Drop for VdbReader {
    fn drop(&mut self) {
        // SAFETY: `handle` came from `rk_open` and is closed once.
        unsafe { ffi::rk_close(self.handle) }
    }
}

/// A slice from a shim pointer; empty for a null pointer or zero length.
///
/// # Safety
/// `p` must point to `len` valid `T`s that outlive `'a`.
unsafe fn slice<'a, T>(p: *const T, len: u32) -> &'a [T] {
    if p.is_null() || len == 0 {
        &[]
    } else {
        // SAFETY: by the caller's contract.
        unsafe { std::slice::from_raw_parts(p, len as usize) }
    }
}

/// Build a [`Spot`] from shim pointers, with bases and qualities from
/// `seq` when given.
///
/// # Safety
/// The pointers in `s` must stay valid for `'a`.
unsafe fn spot_from_raw<'a>(s: &ffi::Spot, seq: Option<(&'a [u8], &'a [u8])>) -> Spot<'a> {
    // SAFETY: by the caller's contract; `ReadType` is `repr(transparent)`
    // over `u8`.
    unsafe {
        let (bases, quals) = seq.unwrap_or_else(|| {
            (
                slice(s.bases.cast::<u8>(), s.bases_len),
                slice(s.quals.cast::<u8>(), s.quals_len),
            )
        });
        Spot {
            bases,
            quals,
            starts: slice(s.read_start, s.nstarts),
            lens: slice(s.read_len, s.nlens),
            types: slice(s.read_type.cast::<ReadType>(), s.nreads),
            name: (!s.name.is_null()).then(|| slice(s.name.cast::<u8>(), s.name_len)),
        }
    }
}

/// The table type the fast path was checked against.
const FAST_TYPESPEC: &str = "NCBI:align:tbl:seq#1.1";
/// Spots compared with VDB when a run is first opened with the fast path.
const SELF_CHECK_SPOTS: u64 = 64;

/// Self-check verdicts, per file (canonical path, size, modification time),
/// so that a decoder's worker threads do not each repeat it.
type CheckKey = (std::path::PathBuf, u64, Option<std::time::SystemTime>);

fn checked() -> &'static std::sync::Mutex<std::collections::HashMap<CheckKey, bool>> {
    static CHECKED: std::sync::OnceLock<
        std::sync::Mutex<std::collections::HashMap<CheckKey, bool>>,
    > = std::sync::OnceLock::new();
    CHECKED.get_or_init(Default::default)
}

fn check_key(path: &Path) -> Option<CheckKey> {
    let meta = std::fs::metadata(path).ok()?;
    Some((
        std::fs::canonicalize(path).ok()?,
        meta.len(),
        meta.modified().ok(),
    ))
}

/// A physical column read as raw KDB blobs, and its current decoded blob.
struct Column {
    handle: *mut c_void,
    raw: Vec<u8>,
    /// Row ids `[first, end)` of the current blob.
    first: u64,
    end: u64,
    /// `None` if the blob is not in a form the fast path decodes.
    blob: Option<Decoded>,
    /// Why the first blob that was not decoded was not.
    reason: Option<&'static str>,
}

struct Decoded {
    rows: vdbfmt::Rows,
    /// 2na packed bytes, or inflated bytes.
    data: Vec<u8>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Coding {
    /// `INSDC:2na:packed`, stored without an encoding (`.READ`).
    Packed2na,
    /// Bytes under `zip_encoding` (`.QUALITY`, `.ALTREAD`).
    ZipBytes,
}

impl Column {
    /// The physical column `name`, or `None` if the table has none.
    fn open(vdb: &VdbReader, name: &str) -> Option<Self> {
        let c_name = CString::new(name).expect("no NUL");
        let mut st = ffi::Status::new();
        // SAFETY: live handle, NUL-terminated name.
        let handle = unsafe { ffi::rk_kcol_open(vdb.handle, c_name.as_ptr(), &mut st) };
        (!handle.is_null()).then(|| Self {
            handle,
            raw: Vec::new(),
            first: 0,
            end: 0,
            blob: None,
            reason: None,
        })
    }

    /// The decoded blob holding `row` (loaded first if needed) and the
    /// row's offset and length in its data; `None` if that blob cannot be
    /// read or decoded here.
    fn row(
        &mut self,
        row: u64,
        coding: Coding,
        inflater: &mut libdeflater::Decompressor,
    ) -> Option<(&[u8], u64, usize)> {
        if !(self.first..self.end).contains(&row) {
            self.load(row, coding, inflater);
        }
        let b = self.blob.as_ref()?;
        let (off, len) = b.rows[(row - self.first) as usize];
        Some((&b.data, off, len as usize))
    }

    fn load(&mut self, row: u64, coding: Coding, inflater: &mut libdeflater::Decompressor) {
        self.blob = None;
        // A failed read covers just this row: VDB then decides for it.
        (self.first, self.end) = (row, row + 1);
        let Ok(id) = i64::try_from(row) else { return };
        let (mut first, mut count, mut size) = (0i64, 0u32, 0usize);
        loop {
            let mut st = ffi::Status::new();
            // SAFETY: live column handle; `raw` has `raw.len()` writable bytes.
            let rc = unsafe {
                ffi::rk_kcol_blob(
                    self.handle,
                    id,
                    &mut first,
                    &mut count,
                    self.raw.as_mut_ptr(),
                    self.raw.len(),
                    &mut size,
                    &mut st,
                )
            };
            match rc {
                0 => break,
                1 => self.raw.resize(size, 0),
                _ => {
                    self.reason.get_or_insert("the blob could not be read");
                    return;
                }
            }
        }
        let Ok(first) = u64::try_from(first) else {
            return;
        };
        if !(first..first + u64::from(count)).contains(&row) {
            self.reason.get_or_insert("the blob does not hold the row");
            return;
        }
        (self.first, self.end) = (first, first + u64::from(count));
        match decode_blob(&self.raw[..size], count, coding, inflater) {
            Ok(b) => self.blob = Some(b),
            Err(e) => {
                self.reason.get_or_insert(e.0);
            }
        }
    }
}

impl Drop for Column {
    fn drop(&mut self) {
        // SAFETY: `handle` came from `rk_kcol_open` and is closed once.
        unsafe { ffi::rk_kcol_close(self.handle) }
    }
}

fn inflate(
    inflater: &mut libdeflater::Decompressor,
    src: &[u8],
    out: &mut Vec<u8>,
    exact: bool,
) -> Result<(), vdbfmt::FormatError> {
    let n = inflater
        .deflate_decompress(src, out)
        .map_err(|_| vdbfmt::FormatError("deflate stream does not inflate"))?;
    if exact && n != out.len() {
        return Err(vdbfmt::FormatError("deflate stream of an unexpected size"));
    }
    out.truncate(n);
    Ok(())
}

/// Decode one blob of `rows` rows; checks every row against the data.
fn decode_blob(
    blob: &[u8],
    rows: u32,
    coding: Coding,
    inflater: &mut libdeflater::Decompressor,
) -> Result<Decoded, vdbfmt::FormatError> {
    let env = vdbfmt::envelope(blob)?;
    let mut page_inflate = |src: &[u8], bound: usize| {
        let mut out = vec![0u8; bound];
        inflate(inflater, src, &mut out, false)?;
        Ok(out)
    };
    let row_map = vdbfmt::page_map(env.page_map, rows, &mut page_inflate)?;
    let (data, elems) = match coding {
        Coding::Packed2na => {
            if !env.header.is_empty() {
                return Err(vdbfmt::FormatError("encoded 2na column"));
            }
            (
                env.payload.to_vec(),
                vdbfmt::bases_2na(env.payload, env.adjust)?,
            )
        }
        Coding::ZipBytes => {
            let h = vdbfmt::zip_header(env.header)?;
            if env.adjust != 0 {
                return Err(vdbfmt::FormatError("pad bits in a zip payload"));
            }
            let size =
                usize::try_from(h.osize).map_err(|_| vdbfmt::FormatError("blob too large"))?;
            let mut out = vec![0u8; size];
            if size > 0 {
                inflate(inflater, env.payload, &mut out, true)?;
            } else if !env.payload.is_empty() {
                return Err(vdbfmt::FormatError("payload for an empty blob"));
            }
            (out, h.osize)
        }
    };
    if row_map
        .iter()
        .any(|&(off, len)| off + u64::from(len) > elems)
    {
        return Err(vdbfmt::FormatError("page map row beyond the payload"));
    }
    Ok(Decoded {
        rows: row_map,
        data,
    })
}

/// Bases and qualities from KDB blobs; lengths, types and names from a VDB
/// cursor; VDB for anything the blobs do not give.
struct FastReader {
    // Field order is drop order, and it matters: the columns hold KDB
    // references taken from `small`'s table, and closing `small` may free the
    // VDB manager they belong to. Columns must be released first.
    read: Column,
    quality: Column,
    /// Absent when the run stores no ambiguity codes.
    alt: Option<Column>,
    /// The cursor without `READ` and `QUALITY`.
    small: VdbReader,
    inflater: libdeflater::Decompressor,
    bases: Vec<u8>,
    quals: Vec<u8>,
    /// The VDB path, opened on first need.
    vdb: Option<VdbReader>,
    names: bool,
    fallback_spots: u64,
    /// Why the first spot whose decoded columns disagreed was handed to VDB.
    mismatch: Option<&'static str>,
}

// SAFETY: as for `VdbReader`; the column handles belong to this reader only.
unsafe impl Send for FastReader {}

impl FastReader {
    /// A fast reader, or `None` when the run is not one the fast path
    /// applies to (then use VDB).
    fn open(path: &Path, names: bool) -> Result<Option<Self>, Error> {
        let small = VdbReader::open(path, names, false)?;
        if small.info.quality != QualityKind::Original
            || small.typespec().ok().as_deref() != Some(FAST_TYPESPEC)
        {
            return Ok(None);
        }
        // Aligned runs rebuild READ from alignments; their SEQUENCE table has
        // no physical READ but these columns.
        if Column::open(&small, "PRIMARY_ALIGNMENT_ID").is_some()
            || Column::open(&small, "CMP_READ").is_some()
        {
            return Ok(None);
        }
        let (Some(read), Some(quality)) = (
            Column::open(&small, "READ"),
            Column::open(&small, "QUALITY"),
        ) else {
            return Ok(None);
        };
        let alt = Column::open(&small, "ALTREAD");
        let mut fast = Self {
            small,
            read,
            quality,
            alt,
            inflater: libdeflater::Decompressor::new(),
            bases: Vec::new(),
            quals: Vec::new(),
            vdb: None,
            names,
            fallback_spots: 0,
            mismatch: None,
        };
        let key = check_key(path);
        let known = key.as_ref().and_then(|k| {
            checked()
                .lock()
                .unwrap_or_else(|e| e.into_inner())
                .get(k)
                .copied()
        });
        let ok = match known {
            Some(ok) => ok,
            None => {
                let ok = fast.self_check()?;
                if let Some(k) = key {
                    checked()
                        .lock()
                        .unwrap_or_else(|e| e.into_inner())
                        .insert(k, ok);
                }
                ok
            }
        };
        Ok(ok.then_some(fast))
    }

    /// Decode the first spots both ways; `false` if any differs.
    fn self_check(&mut self) -> Result<bool, Error> {
        let info = self.small.info.clone();
        let end = info.end_spot().min(info.first_spot + SELF_CHECK_SPOTS);
        for spot in info.first_spot..end {
            if !self.decode(spot) {
                continue; // served through VDB anyway
            }
            let (bases, quals) = (self.bases.clone(), self.quals.clone());
            let theirs = self.vdb()?.read_spot(spot)?;
            if theirs.bases != bases.as_slice() || theirs.quals != quals.as_slice() {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn vdb(&mut self) -> Result<&mut VdbReader, Error> {
        if self.vdb.is_none() {
            self.vdb = Some(VdbReader::open(&self.small.path, self.names, true)?);
        }
        Ok(self.vdb.as_mut().expect("just opened"))
    }

    /// Decode `spot`'s bases and qualities into `self.bases` / `self.quals`;
    /// `false` if its blobs are not decodable here or disagree.
    fn decode(&mut self, spot: u64) -> bool {
        let inflater = &mut self.inflater;
        let Some((packed, r_off, r_len)) = self.read.row(spot, Coding::Packed2na, inflater) else {
            return false;
        };
        vdbfmt::unpack_2na(packed, r_off, r_len, &mut self.bases);
        if let Some(alt) = self.alt.as_mut() {
            let Some((codes, a_off, a_len)) = alt.row(spot, Coding::ZipBytes, inflater) else {
                return false;
            };
            let row = &codes[a_off as usize..a_off as usize + a_len];
            if let Err(e) = vdbfmt::apply_alt(&mut self.bases, row) {
                self.mismatch.get_or_insert(e.0);
                return false;
            }
        }
        let Some((phred, q_off, q_len)) = self.quality.row(spot, Coding::ZipBytes, inflater) else {
            return false;
        };
        if q_len != r_len {
            self.mismatch
                .get_or_insert("READ and QUALITY rows differ in length");
            return false;
        }
        // phred to phred+33 text (`out_qual_text_phred_33`); a plain loop
        // over equal-length slices, which the compiler vectorizes.
        self.quals.resize(q_len, 0);
        let src = &phred[q_off as usize..q_off as usize + q_len];
        for (d, &q) in self.quals.iter_mut().zip(src) {
            *d = q.wrapping_add(33);
        }
        true
    }

    fn read_spot(&mut self, spot: u64) -> Result<Spot<'_>, Error> {
        if spot < self.small.info.first_spot || spot >= self.small.info.end_spot() {
            // Let VDB report it, as on the VDB path.
            return self.small.read_spot(spot);
        }
        if !self.decode(spot) {
            self.fallback_spots += 1;
            return self.vdb()?.read_spot(spot);
        }
        let raw = self.small.read_raw(spot)?;
        // SAFETY: `raw` points into the small cursor's cache, valid until the
        // next read on it, which needs `&mut self`; bases and qualities
        // borrow `self` too.
        unsafe { Ok(spot_from_raw(&raw, Some((&self.bases, &self.quals)))) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A local run, if present (see tests/fast_path.rs).
    fn run() -> Option<std::path::PathBuf> {
        let dir = std::env::var_os("READKIT_SRA_DIR")
            .map(std::path::PathBuf::from)
            .unwrap_or_else(|| {
                std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../.data/sra")
            });
        ["SRR10971381", "SRR13638690", "SRR12999707"]
            .iter()
            .map(|a| dir.join(a).join(a))
            .find(|p| p.exists())
    }

    /// The raw bytes and row count of the blob holding `row`.
    fn raw_blob(col: &Column, row: u64) -> (Vec<u8>, u32) {
        let mut buf = Vec::new();
        loop {
            let (mut first, mut count, mut n) = (0i64, 0u32, 0usize);
            let mut st = ffi::Status::new();
            // SAFETY: live handle; `buf` has `buf.len()` writable bytes.
            let rc = unsafe {
                ffi::rk_kcol_blob(
                    col.handle,
                    row as i64,
                    &mut first,
                    &mut count,
                    buf.as_mut_ptr(),
                    buf.len(),
                    &mut n,
                    &mut st,
                )
            };
            match rc {
                0 => return (buf[..n].to_vec(), count),
                1 => buf.resize(n, 0),
                _ => panic!("blob read failed"),
            }
        }
    }

    #[test]
    fn corrupt_blobs_are_rejected_not_misread() {
        let Some(path) = run() else {
            eprintln!("skipping: no local run");
            return;
        };
        let small = VdbReader::open(&path, false, false).unwrap();
        let mut inflater = libdeflater::Decompressor::new();
        let read = Column::open(&small, "READ").unwrap();
        let quality = Column::open(&small, "QUALITY").unwrap();
        let (r, rows) = raw_blob(&read, 1);
        let (q, qrows) = raw_blob(&quality, 1);
        assert!(decode_blob(&r, rows, Coding::Packed2na, &mut inflater).is_ok());
        assert!(decode_blob(&q, qrows, Coding::ZipBytes, &mut inflater).is_ok());
        // Truncated 2na payload: the page map points past it.
        assert!(decode_blob(&r[..r.len() - 100], rows, Coding::Packed2na, &mut inflater).is_err());
        // A broken deflate stream, or a truncated one.
        let mut bad = q.clone();
        let mid = bad.len() - bad.len() / 3;
        bad[mid] ^= 0xff;
        bad[mid + 1] ^= 0x5a;
        assert!(decode_blob(&bad, qrows, Coding::ZipBytes, &mut inflater).is_err());
        assert!(decode_blob(&q[..q.len() - 10], qrows, Coding::ZipBytes, &mut inflater).is_err());
        // Wrong row count for the page map, wrong coding for the column.
        assert!(decode_blob(&r, rows + 1, Coding::Packed2na, &mut inflater).is_err());
        assert!(decode_blob(&q, qrows, Coding::Packed2na, &mut inflater).is_err());
        assert!(decode_blob(&r, rows, Coding::ZipBytes, &mut inflater).is_err());
    }

    #[test]
    fn undecodable_blobs_go_through_vdb() {
        let Some(path) = run() else {
            eprintln!("skipping: no local run");
            return;
        };
        let mut fast = FastReader::open(&path, true).unwrap().expect("fast path");
        let mut vdb = VdbReader::open(&path, true, true).unwrap();
        // Make the READ blob for spots 1..2000 look undecodable.
        fast.read.blob = None;
        (fast.read.first, fast.read.end) = (1, 2000);
        for spot in 1..2500u64 {
            let a = fast.read_spot(spot).unwrap();
            let (bases, quals, lens, name) = (
                a.bases.to_vec(),
                a.quals.to_vec(),
                a.lens.to_vec(),
                a.name.map(<[u8]>::to_vec),
            );
            let b = vdb.read_spot(spot).unwrap();
            assert_eq!(
                (bases.as_slice(), quals.as_slice(), lens.as_slice()),
                (b.bases, b.quals, b.lens),
                "spot {spot}"
            );
            assert_eq!(name.as_deref(), b.name, "spot {spot}");
        }
        assert_eq!(fast.fallback_spots, 1999);
    }
}
