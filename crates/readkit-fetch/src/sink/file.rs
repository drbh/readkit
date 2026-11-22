//! File sink: `.partial` file, resume ledger, rename and marker on success.

use super::{Opened, Sink};
use crate::error::{Error, ErrorKind};
use crate::object::{ChecksumRepr, Object};
use serde::{Deserialize, Serialize};
use std::fs::{self, File, OpenOptions};
use std::io;
use std::path::{Path, PathBuf};

/// Downloads to a file, resumably.
///
/// While downloading, bytes go to `<dest>.partial` and resume state to
/// `<dest>.ledger.json` (checkpointed about once per second and on failure),
/// so a later download of the same object to the same `dest` fetches only the
/// missing chunks. On success the partial file is synced and renamed to
/// `dest`, the ledger is removed and, unless disabled, a verified marker is
/// written next to it (see [`FileSink::marker_path`]).
///
/// # Verified marker
///
/// `<dest>.verified.json` is a JSON object:
///
/// ```json
/// {"version":1,"url":"https://...","size":1234,
///  "checksum":{"md5":{"hex":"..."}},"verified_at_unix":1731600000}
/// ```
///
/// `checksum` is `"none"`, `{"md5":{"hex":..}}` or
/// `{"s3_multipart_etag":{"etag":..,"parts":..,"part_size":..}}`. When `dest`
/// exists with the marker's size and the marker's checksum equals the
/// object's, the download is skipped and reported as reused.
///
/// Output: the destination path.
#[derive(Debug)]
pub struct FileSink {
    dest: PathBuf,
    marker: bool,
    open: Option<OpenFile>,
}

#[derive(Debug)]
struct OpenFile {
    file: PosFile,
    partial: PathBuf,
    ledger: PathBuf,
    object: Object,
}

#[derive(Serialize, Deserialize)]
struct Marker {
    version: u32,
    url: String,
    size: u64,
    checksum: ChecksumRepr,
    verified_at_unix: u64,
}

impl FileSink {
    /// Download to `dest`, writing a verified marker.
    pub fn new(dest: impl Into<PathBuf>) -> Self {
        Self {
            dest: dest.into(),
            marker: true,
            open: None,
        }
    }

    /// Whether to write (and trust) the verified marker; on by default.
    /// Without it every download fetches again (resume still works).
    pub fn verified_marker(mut self, enabled: bool) -> Self {
        self.marker = enabled;
        self
    }

    /// The marker path for `dest`: `<dest>.verified.json`.
    pub fn marker_path(dest: &Path) -> PathBuf {
        sibling(dest, ".verified.json")
    }

    fn reusable(&self, object: &Object) -> bool {
        if !self.marker {
            return false;
        }
        let Ok(meta) = fs::metadata(&self.dest) else {
            return false;
        };
        let Ok(bytes) = fs::read(Self::marker_path(&self.dest)) else {
            return false;
        };
        let Ok(m) = serde_json::from_slice::<Marker>(&bytes) else {
            return false;
        };
        meta.is_file()
            && meta.len() == object.size
            && m.size == object.size
            && m.checksum == ChecksumRepr::from(&object.checksum)
    }

    fn opened(&self) -> &OpenFile {
        self.open.as_ref().expect("FileSink used before open")
    }
}

fn sibling(dest: &Path, suffix: &str) -> PathBuf {
    let name = dest.file_name().unwrap_or_default().to_string_lossy();
    dest.with_file_name(format!("{name}{suffix}"))
}

fn write_atomic(path: &Path, bytes: &[u8]) -> io::Result<()> {
    let tmp = sibling(path, ".tmp");
    fs::write(&tmp, bytes)?;
    fs::rename(tmp, path)
}

impl Sink for FileSink {
    type Output = PathBuf;

    fn open(&mut self, object: &Object) -> Result<Opened, Error> {
        if self.reusable(object) {
            return Ok(Opened::Complete);
        }
        if self.dest.file_name().is_none() {
            return Err(Error::new(
                ErrorKind::InvalidInput,
                format!("destination {} has no file name", self.dest.display()),
            ));
        }
        if let Some(dir) = self.dest.parent().filter(|d| !d.as_os_str().is_empty()) {
            fs::create_dir_all(dir)?;
        }
        let partial = sibling(&self.dest, ".partial");
        let ledger = sibling(&self.dest, ".ledger.json");
        let resume = if partial.exists() {
            fs::read(&ledger).ok()
        } else {
            None
        };
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(false)
            .open(&partial)?;
        if resume.is_none() {
            file.set_len(0)?;
        }
        file.set_len(object.size)?;
        self.open = Some(OpenFile {
            file: PosFile::new(file),
            partial,
            ledger,
            object: object.clone(),
        });
        Ok(resume.map_or(Opened::Empty, Opened::Resume))
    }

    fn write_at(&self, offset: u64, data: Vec<u8>) -> Result<(), Error> {
        Ok(self.opened().file.write_all_at(&data, offset)?)
    }

    fn read_at(&self, offset: u64, buf: &mut [u8]) -> Result<(), Error> {
        Ok(self.opened().file.read_exact_at(buf, offset)?)
    }

    fn checkpoint(&self, state: &[u8]) -> Result<(), Error> {
        Ok(write_atomic(&self.opened().ledger, state)?)
    }

    fn commit(self) -> Result<PathBuf, Error> {
        let Some(open) = self.open else {
            // Opened::Complete: nothing to do.
            return Ok(self.dest);
        };
        open.file.sync_all()?;
        drop(open.file);
        fs::rename(&open.partial, &self.dest)?;
        if self.marker {
            let marker = Marker {
                version: 1,
                url: open.object.url.clone(),
                size: open.object.size,
                checksum: ChecksumRepr::from(&open.object.checksum),
                verified_at_unix: std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .map_or(0, |d| d.as_secs()),
            };
            let json = serde_json::to_vec(&marker)
                .map_err(|e| Error::new(ErrorKind::Internal, e.to_string()))?;
            write_atomic(&Self::marker_path(&self.dest), &json)?;
        }
        let _ = fs::remove_file(&open.ledger);
        Ok(self.dest)
    }

    fn abort(self, _error: &Error) {
        // Keep `.partial` and the ledger for a later resume.
    }
}

/// Positional file I/O. On unix this is `pread`/`pwrite` (`FileExt`), safe to
/// use from many threads at once; elsewhere a mutex serializes seek + I/O.
#[derive(Debug)]
struct PosFile {
    #[cfg(unix)]
    file: File,
    #[cfg(not(unix))]
    file: portable::Locked,
}

impl PosFile {
    fn new(file: File) -> Self {
        Self {
            #[cfg(unix)]
            file,
            #[cfg(not(unix))]
            file: portable::Locked::new(file),
        }
    }

    fn write_all_at(&self, buf: &[u8], offset: u64) -> io::Result<()> {
        #[cfg(unix)]
        {
            std::os::unix::fs::FileExt::write_all_at(&self.file, buf, offset)
        }
        #[cfg(not(unix))]
        {
            self.file.write_all_at(buf, offset)
        }
    }

    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> io::Result<()> {
        #[cfg(unix)]
        {
            std::os::unix::fs::FileExt::read_exact_at(&self.file, buf, offset)
        }
        #[cfg(not(unix))]
        {
            self.file.read_exact_at(buf, offset)
        }
    }

    fn sync_all(&self) -> io::Result<()> {
        self.file.sync_all()
    }
}

#[cfg(any(not(unix), test))]
mod portable {
    use std::fs::File;
    use std::io::{self, Read, Seek, SeekFrom, Write};
    use std::sync::Mutex;

    /// Seek + read/write under a lock: a shared cursor must never interleave.
    #[derive(Debug)]
    pub(super) struct Locked(Mutex<File>);

    impl Locked {
        pub(super) fn new(file: File) -> Self {
            Self(Mutex::new(file))
        }

        fn lock(&self) -> std::sync::MutexGuard<'_, File> {
            self.0.lock().unwrap_or_else(|e| e.into_inner())
        }

        pub(super) fn write_all_at(&self, buf: &[u8], offset: u64) -> io::Result<()> {
            let mut f = self.lock();
            f.seek(SeekFrom::Start(offset))?;
            f.write_all(buf)
        }

        pub(super) fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> io::Result<()> {
            let mut f = self.lock();
            f.seek(SeekFrom::Start(offset))?;
            f.read_exact(buf)
        }

        pub(super) fn sync_all(&self) -> io::Result<()> {
            self.lock().sync_all()
        }
    }

    #[cfg(test)]
    #[test]
    fn locked_positional_io_round_trips() {
        let path =
            std::env::temp_dir().join(format!("readkit-fetch-portable-{}", std::process::id()));
        let file = std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(true)
            .open(&path)
            .unwrap();
        file.set_len(8).unwrap();
        let f = Locked::new(file);
        f.write_all_at(b"cd", 2).unwrap();
        f.write_all_at(b"ab", 0).unwrap();
        let mut buf = [0u8; 4];
        f.read_exact_at(&mut buf, 0).unwrap();
        assert_eq!(&buf, b"abcd");
        f.sync_all().unwrap();
        drop(f);
        let _ = std::fs::remove_file(path);
    }
}
