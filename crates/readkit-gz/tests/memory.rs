//! Region decoding stays within its documented memory bound
//! ([`DecodeOptions`]), counted by the allocator, so the check is the same
//! on every platform. Its own test binary: the counting allocator sees the
//! whole process, so nothing else may run alongside.
#![cfg(feature = "parallel")]

mod common;

use readkit_gz::{DecodeOptions, ParGzReader, Split};
use std::alloc::{GlobalAlloc, Layout, System};
use std::io::Read;
use std::sync::atomic::{AtomicUsize, Ordering};

/// The system allocator, counting the bytes allocated, their peak, and the
/// total of large allocations.
struct Counting;

static LIVE: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);
static LARGE: AtomicUsize = AtomicUsize::new(0);

/// Allocations this large are buffers, not bookkeeping. Most allocators
/// serve them from separate mappings, and keep them resident when freed.
const LARGE_SIZE: usize = 256 << 10;

fn grew(by: usize) {
    let now = LIVE.fetch_add(by, Ordering::Relaxed) + by;
    PEAK.fetch_max(now, Ordering::Relaxed);
    if by >= LARGE_SIZE {
        LARGE.fetch_add(by, Ordering::Relaxed);
    }
}

// SAFETY: every call forwards to `System` with the same arguments.
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let p = unsafe { System.alloc(layout) };
        if !p.is_null() {
            grew(layout.size());
        }
        p
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        let p = unsafe { System.alloc_zeroed(layout) };
        if !p.is_null() {
            grew(layout.size());
        }
        p
    }
    unsafe fn dealloc(&self, p: *mut u8, layout: Layout) {
        unsafe { System.dealloc(p, layout) };
        LIVE.fetch_sub(layout.size(), Ordering::Relaxed);
    }
    unsafe fn realloc(&self, p: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        let q = unsafe { System.realloc(p, layout, size) };
        if !q.is_null() {
            // Counted as both blocks at once, as a copying realloc holds them.
            grew(size);
            LIVE.fetch_sub(layout.size(), Ordering::Relaxed);
        }
        q
    }
}

#[global_allocator]
static ALLOCATOR: Counting = Counting;

#[test]
fn region_decoding_stays_within_its_memory_bound() {
    const MIB: usize = 1 << 20;
    let data = common::fastq(60_000, 11);
    let gz = common::flate2_gz(&data, 6);
    let region = 256 << 10;
    let mut chunk = vec![0u8; 64 << 10];
    for threads in [1, 4] {
        let in_flight = threads + 2;
        let options = DecodeOptions {
            threads: Some(threads),
            split: Split::Always,
            region_size: region,
            ..Default::default()
        };
        let input = gz.clone();
        let before = LIVE.load(Ordering::Relaxed);
        PEAK.store(before, Ordering::Relaxed);
        LARGE.store(0, Ordering::Relaxed);
        let mut reader = ParGzReader::new(input, gz.len() as u64, &options).unwrap();
        let mut at = 0;
        loop {
            let n = reader.read(&mut chunk).unwrap();
            if n == 0 {
                break;
            }
            assert!(chunk[..n] == data[at..at + n], "output at {at}");
            at += n;
        }
        assert_eq!(at, data.len());
        drop(reader);
        let peak = PEAK.load(Ordering::Relaxed) - before;
        let large = LARGE.load(Ordering::Relaxed);
        // The bound documented on `DecodeOptions`: each region in flight up
        // to three times its decoded size, three more decoded regions with
        // the reader, and a compressed buffer per thread. Region buffers are
        // reserved for five times `region_size` until a region has been
        // decoded, and 1/8 more than needed.
        let ratio = data.len().div_ceil(gz.len());
        let decoded = (region * ratio).max(5 * region) * 9 / 8;
        let bound = (3 * in_flight + 3) * decoded + (threads + 1) * (region + MIB) + MIB;
        assert!(
            peak <= bound,
            "{threads} threads: peak {peak} bytes over the bound {bound}"
        );
        // Buffers are reused, not allocated per region: all large
        // allocations together fit in the bound too, for any input length.
        // Allocating them per region came to 3.5 times the output here.
        assert!(
            large <= bound,
            "{threads} threads: {large} bytes of large allocations, over {bound}"
        );
        eprintln!("{threads} threads: peak {peak}, large {large}, bound {bound} bytes");
    }
}
