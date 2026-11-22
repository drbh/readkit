//! Buffer pools.
//!
//! A [`Pool`] hands out growable buffers for byte columns, offset arrays and
//! tag arrays. Every buffer remembers its pool; dropping it (usually by
//! dropping the [`Batch`](crate::Batch) that owns it) clears it and returns it
//! with its capacity intact. The pool never blocks: when it has no idle buffer
//! it allocates a fresh one, and when it already holds `max_idle` idle buffers
//! of a kind it frees the returned one instead of keeping it.

use std::fmt;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, MutexGuard};

/// Default number of idle buffers kept per element kind.
pub const DEFAULT_MAX_IDLE: usize = 64;

/// A shared, thread-safe pool of column buffers.
///
/// Cheap to clone (an `Arc`); all clones share the same buffers and counters.
#[derive(Clone)]
pub struct Pool {
    inner: Arc<Inner>,
}

#[doc(hidden)]
pub struct Inner {
    max_idle: usize,
    bytes: Mutex<Vec<Vec<u8>>>,
    offsets: Mutex<Vec<Vec<i32>>>,
    tags: Mutex<Vec<Vec<u16>>>,
    fresh: AtomicU64,
    reused: AtomicU64,
    returned: AtomicU64,
    discarded: AtomicU64,
}

/// Counters describing a pool's activity since it was created.
///
/// `fresh + reused` is the number of buffers handed out. In steady state
/// `fresh` stops growing.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct PoolStats {
    /// Buffers allocated because no idle buffer was available.
    pub fresh: u64,
    /// Buffers handed out from the idle list.
    pub reused: u64,
    /// Buffers returned and kept as idle.
    pub returned: u64,
    /// Buffers returned but freed because the idle cap was reached.
    pub discarded: u64,
    /// Idle buffers currently held (all kinds).
    pub idle: usize,
}

mod sealed {
    pub trait Sealed {}
    impl Sealed for u8 {}
    impl Sealed for i32 {}
    impl Sealed for u16 {}
}

/// Element types a pool stores buffers of (`u8`, `i32`, `u16`). Sealed.
pub trait Element:
    sealed::Sealed + Copy + Default + Send + Sync + std::panic::RefUnwindSafe + 'static
{
    #[doc(hidden)]
    fn idle(inner: &Inner) -> &Mutex<Vec<Vec<Self>>>;
}

impl Element for u8 {
    fn idle(inner: &Inner) -> &Mutex<Vec<Vec<Self>>> {
        &inner.bytes
    }
}
impl Element for i32 {
    fn idle(inner: &Inner) -> &Mutex<Vec<Vec<Self>>> {
        &inner.offsets
    }
}
impl Element for u16 {
    fn idle(inner: &Inner) -> &Mutex<Vec<Vec<Self>>> {
        &inner.tags
    }
}

fn lock<T>(m: &Mutex<T>) -> MutexGuard<'_, T> {
    // A panic while holding the lock cannot leave a Vec<Vec<_>> inconsistent.
    m.lock().unwrap_or_else(|e| e.into_inner())
}

impl Pool {
    /// A pool keeping at most [`DEFAULT_MAX_IDLE`] idle buffers per kind.
    pub fn new() -> Self {
        Self::with_max_idle(DEFAULT_MAX_IDLE)
    }

    /// A pool keeping at most `max_idle` idle buffers of each element kind
    /// (bytes, offsets, tags). `0` disables recycling.
    pub fn with_max_idle(max_idle: usize) -> Self {
        Self {
            inner: Arc::new(Inner {
                max_idle,
                bytes: Mutex::new(Vec::new()),
                offsets: Mutex::new(Vec::new()),
                tags: Mutex::new(Vec::new()),
                fresh: AtomicU64::new(0),
                reused: AtomicU64::new(0),
                returned: AtomicU64::new(0),
                discarded: AtomicU64::new(0),
            }),
        }
    }

    /// The idle cap per element kind.
    pub fn max_idle(&self) -> usize {
        self.inner.max_idle
    }

    /// A snapshot of this pool's counters.
    pub fn stats(&self) -> PoolStats {
        let i = &self.inner;
        PoolStats {
            fresh: i.fresh.load(Ordering::Relaxed),
            reused: i.reused.load(Ordering::Relaxed),
            returned: i.returned.load(Ordering::Relaxed),
            discarded: i.discarded.load(Ordering::Relaxed),
            idle: lock(&i.bytes).len() + lock(&i.offsets).len() + lock(&i.tags).len(),
        }
    }

    /// Takes an empty buffer with room for at least `capacity` elements.
    /// Never blocks: allocates when no idle buffer is available.
    pub fn take<T: Element>(&self, capacity: usize) -> PooledVec<T> {
        let recycled = lock(T::idle(&self.inner)).pop();
        let vec = match recycled {
            Some(mut v) => {
                self.inner.reused.fetch_add(1, Ordering::Relaxed);
                v.reserve(capacity);
                v
            }
            None => {
                self.inner.fresh.fetch_add(1, Ordering::Relaxed);
                Vec::with_capacity(capacity)
            }
        };
        PooledVec {
            vec,
            pool: Arc::clone(&self.inner),
        }
    }

    /// Whether two handles refer to the same pool.
    pub fn same_pool(&self, other: &Pool) -> bool {
        Arc::ptr_eq(&self.inner, &other.inner)
    }
}

impl Default for Pool {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Debug for Pool {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Pool")
            .field("max_idle", &self.inner.max_idle)
            .field("stats", &self.stats())
            .finish()
    }
}

/// A `Vec` borrowed from a [`Pool`]; returned to it (cleared, capacity kept)
/// on drop.
pub struct PooledVec<T: Element> {
    vec: Vec<T>,
    pool: Arc<Inner>,
}

impl<T: Element> PooledVec<T> {
    /// The buffer contents.
    pub fn as_slice(&self) -> &[T] {
        &self.vec
    }

    /// Mutable access to the underlying `Vec` (for producers filling columns).
    pub fn vec_mut(&mut self) -> &mut Vec<T> {
        &mut self.vec
    }

    /// Number of elements.
    pub fn len(&self) -> usize {
        self.vec.len()
    }

    /// Whether the buffer is empty.
    pub fn is_empty(&self) -> bool {
        self.vec.is_empty()
    }
}

impl<T: Element> std::ops::Deref for PooledVec<T> {
    type Target = [T];
    fn deref(&self) -> &[T] {
        &self.vec
    }
}

impl<T: Element> Drop for PooledVec<T> {
    fn drop(&mut self) {
        let mut vec = std::mem::take(&mut self.vec);
        vec.clear();
        let mut idle = lock(T::idle(&self.pool));
        if idle.len() < self.pool.max_idle {
            idle.push(vec);
            drop(idle);
            self.pool.returned.fetch_add(1, Ordering::Relaxed);
        } else {
            drop(idle);
            self.pool.discarded.fetch_add(1, Ordering::Relaxed);
        }
    }
}

impl<T: Element + fmt::Debug> fmt::Debug for PooledVec<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PooledVec")
            .field("len", &self.vec.len())
            .field("capacity", &self.vec.capacity())
            .finish()
    }
}
