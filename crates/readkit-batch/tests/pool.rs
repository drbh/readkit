//! Pool behaviour.

use readkit_batch::{Batch, BatchBuilder, Pool, Read, Single};

fn make(pool: &Pool, n: usize) -> Batch<Single> {
    let mut b = BatchBuilder::<Single>::with_pool(pool, n);
    for _ in 0..n {
        b.push_read(Read::named(b"name", b"ACGTACGT", b"IIIIIIII"))
            .unwrap();
    }
    b.finish()
}

#[test]
fn dropping_a_batch_returns_its_buffers() {
    let pool = Pool::new();
    let batch = make(&pool, 10);
    let s = pool.stats();
    assert_eq!((s.fresh, s.idle), (5, 0)); // seq/qual/offsets + name offsets/bytes
    drop(batch);
    let s = pool.stats();
    assert_eq!((s.returned, s.idle), (5, 5));
    let _batch = make(&pool, 10);
    let s = pool.stats();
    assert_eq!((s.fresh, s.reused, s.idle), (5, 5, 0));
}

#[test]
fn steady_state_allocates_nothing_new() {
    let pool = Pool::new();
    for _ in 0..3 {
        drop(make(&pool, 100));
    }
    let warm = pool.stats();
    for _ in 0..100 {
        drop(make(&pool, 100));
    }
    let s = pool.stats();
    assert_eq!(s.fresh, warm.fresh);
    assert_eq!(s.reused - warm.reused, 500);
}

#[test]
fn never_blocks_when_consumers_hold_batches() {
    let pool = Pool::with_max_idle(2);
    // a consumer (say Python holding Arrow arrays) keeps everything
    let held: Vec<_> = (0..50).map(|_| make(&pool, 4)).collect();
    let s = pool.stats();
    assert_eq!(s.fresh, 250);
    assert_eq!(held.iter().map(Batch::len).sum::<usize>(), 200);
}

#[test]
fn idle_cap_is_respected() {
    let pool = Pool::with_max_idle(3);
    let held: Vec<_> = (0..10).map(|_| make(&pool, 4)).collect();
    drop(held);
    let s = pool.stats();
    // three idle buffers per element kind: u8 and i32 (no tags here)
    assert_eq!(s.idle, 6);
    assert_eq!(s.returned + s.discarded, 50);
    assert_eq!(s.discarded, 44);
}

#[test]
fn pool_is_shared_across_threads() {
    let pool = Pool::new();
    let handles: Vec<_> = (0..4)
        .map(|_| {
            let pool = pool.clone();
            std::thread::spawn(move || make(&pool, 8))
        })
        .collect();
    for h in handles {
        drop(h.join().unwrap());
    }
    // every buffer ever allocated is back (threads may have reused some)
    let s = pool.stats();
    assert_eq!(s.idle as u64, s.fresh);
    assert_eq!(s.fresh + s.reused, 20);
    assert!(pool.same_pool(&pool.clone()));
}
