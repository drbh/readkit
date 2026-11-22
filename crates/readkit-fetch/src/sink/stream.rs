//! In-order streaming sink.

use super::{Access, Opened, Sink};
use crate::error::{Error, ErrorKind};
use crate::object::Object;
use std::collections::VecDeque;
use std::io::{self, Read};
use std::sync::{Arc, Condvar, Mutex, MutexGuard};

/// Delivers the object's bytes, in file order, to a [`StreamReader`] while the
/// download is still running.
///
/// Memory is bounded by about twice `buffer` bytes: at most `buffer` bytes
/// are fetched ahead of the reader's position, and at most `buffer` bytes
/// wait in the queue for the reader. The checksum is computed in order; a
/// mismatch is reported as a read error after the last byte (the stream cannot
/// be repaired, since earlier bytes were already handed out). Dropping the
/// reader cancels the download.
///
/// Output: `()`.
#[derive(Debug)]
pub struct StreamSink {
    shared: Arc<Shared>,
    buffer: u64,
}

/// The reading end of a [`StreamSink`]: `impl Read`, yielding bytes in file
/// order, then `Ok(0)` once the download is complete and verified, or an
/// error (an [`io::Error`] wrapping the crate's [`Error`](crate::Error)) if it
/// failed.
#[derive(Debug)]
pub struct StreamReader {
    shared: Arc<Shared>,
    ended: bool,
}

#[derive(Debug)]
struct Shared {
    state: Mutex<State>,
    cv: Condvar,
}

#[derive(Debug, Default)]
struct State {
    queue: VecDeque<Vec<u8>>,
    /// Consumed bytes of `queue[0]`.
    head: usize,
    queued: u64,
    capacity: u64,
    end: Option<Result<(), Error>>,
    reader_gone: bool,
    interrupted: bool,
    next_offset: u64,
}

impl StreamSink {
    /// A sink buffering about `buffer` bytes (at least 1 MiB is sensible: a
    /// chunk larger than the buffer is still accepted, one at a time).
    pub fn new(buffer: usize) -> (StreamSink, StreamReader) {
        let buffer = buffer.max(1) as u64;
        let shared = Arc::new(Shared {
            state: Mutex::new(State {
                capacity: buffer,
                ..State::default()
            }),
            cv: Condvar::new(),
        });
        (
            StreamSink {
                shared: shared.clone(),
                buffer,
            },
            StreamReader {
                shared,
                ended: false,
            },
        )
    }

    fn close(&self, end: Result<(), Error>) {
        let mut s = self.shared.lock();
        s.end.get_or_insert(end);
        self.shared.cv.notify_all();
    }
}

impl Shared {
    fn lock(&self) -> MutexGuard<'_, State> {
        self.state.lock().unwrap_or_else(|e| e.into_inner())
    }
}

impl Sink for StreamSink {
    type Output = ();

    fn open(&mut self, _object: &Object) -> Result<Opened, Error> {
        Ok(Opened::Empty)
    }

    fn access(&self) -> Access {
        Access::Sequential {
            window: self.buffer,
        }
    }

    fn write_at(&self, offset: u64, data: Vec<u8>) -> Result<(), Error> {
        let mut s = self.shared.lock();
        if offset != s.next_offset {
            return Err(Error::new(
                ErrorKind::Internal,
                format!("stream write at {offset}, expected {}", s.next_offset),
            ));
        }
        loop {
            if s.reader_gone {
                return Err(Error::new(ErrorKind::Cancelled, "stream reader dropped"));
            }
            if s.interrupted {
                return Err(Error::cancelled());
            }
            if s.queued == 0 || s.queued + data.len() as u64 <= s.capacity {
                break;
            }
            s = self.shared.cv.wait(s).unwrap_or_else(|e| e.into_inner());
        }
        s.next_offset += data.len() as u64;
        if !data.is_empty() {
            s.queued += data.len() as u64;
            s.queue.push_back(data);
            self.shared.cv.notify_all();
        }
        Ok(())
    }

    fn interrupt(&self) {
        self.shared.lock().interrupted = true;
        self.shared.cv.notify_all();
    }

    fn commit(self) -> Result<(), Error> {
        self.close(Ok(()));
        Ok(())
    }

    fn abort(self, error: &Error) {
        self.close(Err(error.clone()));
    }
}

impl Drop for StreamSink {
    fn drop(&mut self) {
        // A sink dropped without commit or abort must not leave the reader
        // waiting forever.
        self.close(Err(Error::new(
            ErrorKind::Cancelled,
            "stream ended before the download finished",
        )));
    }
}

impl Read for StreamReader {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        let mut s = self.shared.lock();
        loop {
            if let Some(front) = s.queue.front() {
                let avail = &front[s.head..];
                let n = avail.len().min(out.len());
                out[..n].copy_from_slice(&avail[..n]);
                let done = n == avail.len();
                s.head += n;
                s.queued -= n as u64;
                if done {
                    s.queue.pop_front();
                    s.head = 0;
                }
                self.shared.cv.notify_all();
                return Ok(n);
            }
            match &s.end {
                Some(Ok(())) => {
                    self.ended = true;
                    return Ok(0);
                }
                Some(Err(e)) => {
                    self.ended = true;
                    return Err(e.clone().into());
                }
                None => s = self.shared.cv.wait(s).unwrap_or_else(|e| e.into_inner()),
            }
        }
    }
}

impl StreamReader {
    /// Whether a read has returned end of stream or the final error.
    pub fn ended(&self) -> bool {
        self.ended
    }
}

impl Drop for StreamReader {
    fn drop(&mut self) {
        let mut s = self.shared.lock();
        s.reader_gone = true;
        s.queue.clear();
        s.queued = 0;
        self.shared.cv.notify_all();
    }
}
