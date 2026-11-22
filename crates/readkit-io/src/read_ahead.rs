//! The read-ahead helper.

use std::fmt;
use std::marker::PhantomData;

use readkit_batch::{Batch, Layout};

use crate::BatchSource;

/// Wraps a source and reads at most one batch ahead of it on a dedicated
/// thread, so the wrapping stage and its source overlap.
///
/// Tool authors use this inside a stage that consumes a [`BatchSource`];
/// pipeline authors never place it themselves. Memory stays bounded: the
/// thread holds at most one finished batch (handed over through a
/// rendezvous channel) plus the one it is reading.
///
/// - On `wasm32` targets, or with [`ReadAhead::sequential`], it reads on the
///   caller's thread instead (the knob a tool exposes to turn read-ahead off).
/// - Errors are the source's own, unchanged.
/// - After the source returns `None` or an error, `next_batch` returns
///   `Ok(None)`.
/// - Dropping it stops the thread and joins it; the join waits at most for
///   the read in progress. The source is dropped on that thread.
/// - A panic in the source resumes on the caller's thread at the next
///   `next_batch`.
pub struct ReadAhead<L: Layout, S: BatchSource<L>> {
    mode: Mode<L, S>,
}

enum Mode<L: Layout, S: BatchSource<L>> {
    Sequential(Option<S>, PhantomData<fn() -> L>),
    #[cfg(not(target_arch = "wasm32"))]
    Thread(native::Worker<L, S::Error>),
}

impl<L, S> ReadAhead<L, S>
where
    L: Layout,
    S: BatchSource<L> + Send + 'static,
{
    /// Reads ahead on a new thread (sequential on wasm32, or if the OS
    /// refuses to start a thread).
    pub fn new(source: S) -> Self {
        Self::with_read_ahead(source, true)
    }

    /// Reads on the caller's thread; no read-ahead.
    pub fn sequential(source: S) -> Self {
        Self {
            mode: Mode::Sequential(Some(source), PhantomData),
        }
    }

    /// Reads ahead if `enabled`, else sequentially.
    pub fn with_read_ahead(source: S, enabled: bool) -> Self {
        if !enabled {
            return Self::sequential(source);
        }
        #[cfg(not(target_arch = "wasm32"))]
        {
            match native::Worker::spawn(source) {
                Ok(worker) => Self {
                    mode: Mode::Thread(worker),
                },
                Err(source) => Self::sequential(source),
            }
        }
        #[cfg(target_arch = "wasm32")]
        {
            Self::sequential(source)
        }
    }
}

impl<L: Layout, S: BatchSource<L>> ReadAhead<L, S> {
    /// Whether a read-ahead thread is in use.
    pub fn is_threaded(&self) -> bool {
        match &self.mode {
            Mode::Sequential(..) => false,
            #[cfg(not(target_arch = "wasm32"))]
            Mode::Thread(_) => true,
        }
    }
}

impl<L: Layout, S: BatchSource<L>> BatchSource<L> for ReadAhead<L, S> {
    type Error = S::Error;

    fn next_batch(&mut self) -> Result<Option<Batch<L>>, S::Error> {
        match &mut self.mode {
            Mode::Sequential(slot, _) => {
                let Some(source) = slot else {
                    return Ok(None);
                };
                let r = source.next_batch();
                if !matches!(r, Ok(Some(_))) {
                    *slot = None;
                }
                r
            }
            #[cfg(not(target_arch = "wasm32"))]
            Mode::Thread(worker) => worker.next(),
        }
    }
}

impl<L: Layout, S: BatchSource<L>> fmt::Debug for ReadAhead<L, S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ReadAhead")
            .field("threaded", &self.is_threaded())
            .finish()
    }
}

#[cfg(not(target_arch = "wasm32"))]
mod native {
    use std::sync::mpsc::{Receiver, sync_channel};
    use std::sync::{Arc, Mutex};
    use std::thread::JoinHandle;

    use readkit_batch::{Batch, Layout};

    use crate::BatchSource;

    type Item<L, E> = Result<Batch<L>, E>;

    pub(super) struct Worker<L: Layout, E> {
        rx: Option<Receiver<Item<L, E>>>,
        handle: Option<JoinHandle<()>>,
    }

    impl<L: Layout, E: Send + 'static> Worker<L, E> {
        /// Spawns the thread; gives the source back if that fails.
        pub(super) fn spawn<S>(source: S) -> Result<Self, S>
        where
            S: BatchSource<L, Error = E> + Send + 'static,
        {
            // Rendezvous: the thread blocks in `send` with one finished
            // batch until the consumer takes it, so it is at most one ahead.
            let (tx, rx) = sync_channel::<Item<L, E>>(0);
            // The slot lets us recover the source if the spawn fails.
            let slot = Arc::new(Mutex::new(Some(source)));
            let theirs = Arc::clone(&slot);
            let spawned = std::thread::Builder::new()
                .name("readkit-read-ahead".into())
                .spawn(move || {
                    let taken = theirs.lock().unwrap_or_else(|e| e.into_inner()).take();
                    drop(theirs);
                    let Some(mut source) = taken else { return };
                    loop {
                        let item = match source.next_batch() {
                            Ok(Some(b)) => Ok(b),
                            Ok(None) => break,
                            Err(e) => Err(e),
                        };
                        let last = item.is_err();
                        if tx.send(item).is_err() || last {
                            break;
                        }
                    }
                });
            match spawned {
                Ok(handle) => Ok(Self {
                    rx: Some(rx),
                    handle: Some(handle),
                }),
                Err(_) => {
                    let source = slot
                        .lock()
                        .unwrap_or_else(|e| e.into_inner())
                        .take()
                        .expect("source is still in its slot when the spawn fails");
                    Err(source)
                }
            }
        }
    }

    impl<L: Layout, E> Worker<L, E> {
        pub(super) fn next(&mut self) -> Result<Option<Batch<L>>, E> {
            let Some(rx) = &self.rx else {
                return Ok(None);
            };
            match rx.recv() {
                Ok(Ok(b)) => Ok(Some(b)),
                Ok(Err(e)) => {
                    self.stop();
                    Err(e)
                }
                Err(_) => {
                    // Thread finished: end of input, or a panic to resume.
                    self.rx = None;
                    if let Some(h) = self.handle.take()
                        && let Err(panic) = h.join()
                    {
                        std::panic::resume_unwind(panic);
                    }
                    Ok(None)
                }
            }
        }

        fn stop(&mut self) {
            // Dropping the receiver makes a blocked `send` fail, so the
            // thread exits after at most the read in progress.
            self.rx = None;
            if let Some(h) = self.handle.take() {
                let _ = h.join();
            }
        }
    }

    impl<L: Layout, E> Drop for Worker<L, E> {
        fn drop(&mut self) {
            self.stop();
        }
    }
}
