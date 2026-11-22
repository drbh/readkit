//! In-memory sink.

use super::{Opened, Sink};
use crate::error::{Error, ErrorKind};
use crate::object::Object;
use std::sync::Mutex;

/// Downloads into a `Vec<u8>`. Not resumable.
///
/// Output: the object's bytes.
#[derive(Debug, Default)]
pub struct MemorySink {
    buf: Mutex<Vec<u8>>,
}

impl MemorySink {
    /// An empty memory sink; the buffer is allocated when the download starts.
    pub fn new() -> Self {
        Self::default()
    }

    fn lock(&self) -> std::sync::MutexGuard<'_, Vec<u8>> {
        self.buf.lock().unwrap_or_else(|e| e.into_inner())
    }

    fn range(&self, offset: u64, len: usize) -> Result<std::ops::Range<usize>, Error> {
        let start = usize::try_from(offset).ok();
        start
            .and_then(|s| Some(s..s.checked_add(len)?))
            .filter(|r| r.end <= self.lock().len())
            .ok_or_else(|| Error::new(ErrorKind::Io, "write beyond the object's size"))
    }
}

impl Sink for MemorySink {
    type Output = Vec<u8>;

    fn open(&mut self, object: &Object) -> Result<Opened, Error> {
        let size = usize::try_from(object.size)
            .map_err(|_| Error::new(ErrorKind::InvalidInput, "object too large for memory"))?;
        let buf = self.buf.get_mut().unwrap_or_else(|e| e.into_inner());
        buf.clear();
        buf.resize(size, 0);
        Ok(Opened::Empty)
    }

    fn write_at(&self, offset: u64, data: Vec<u8>) -> Result<(), Error> {
        let r = self.range(offset, data.len())?;
        self.lock()[r].copy_from_slice(&data);
        Ok(())
    }

    fn read_at(&self, offset: u64, buf: &mut [u8]) -> Result<(), Error> {
        let r = self.range(offset, buf.len())?;
        buf.copy_from_slice(&self.lock()[r]);
        Ok(())
    }

    fn commit(self) -> Result<Vec<u8>, Error> {
        Ok(self.buf.into_inner().unwrap_or_else(|e| e.into_inner()))
    }

    fn abort(self, _error: &Error) {}
}
