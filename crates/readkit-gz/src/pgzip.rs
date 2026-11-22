//! Parallel decompression of a single gzip member, the approach of pugz and
//! rapidgzip. A deflate stream can only be decoded from its start, so:
//!
//! 1. The compressed file is cut into regions. A worker finds a dynamic
//!    Huffman block that starts in its region by trying every bit offset and
//!    keeping the first whose header is valid and whose whole block decodes.
//! 2. It decodes from there, without the 32 KiB of history before the region:
//!    a back-reference into that unknown window becomes a placeholder symbol.
//!    Once 32 KiB in a row have no placeholder no later copy can reach one, so
//!    from there it decodes straight to bytes.
//! 3. It stops at the first block boundary at or after the next region's
//!    start, which is where the next worker should have started.
//! 4. A coordinator takes regions in order, fills in placeholders from the
//!    previous region's last 32 KiB, checks the CRC and length, and sends the
//!    bytes on. A region whose start does not match the previous region's end
//!    (no block found, or a false start) is decoded again from the known end.
//!
//! A wrong guess costs time, never output: only regions that continue exactly
//! from verified output are used, and the trailer CRC covers everything.
//!
//! Members after the first are decoded the same way while they stay large
//! (see [`decompress`]); otherwise the rest of the input is decoded by the
//! sequential [`GzReader`].
use crate::error::Error;
use crate::read_at::{Cursor, ReadAt};
use crate::sequential::GzReader;
use std::io::BufReader;

const WINDOW: usize = 32768;
/// Compressed bytes per region.
const REGION: usize = 2 << 20;
/// Bytes read past a region for the block that crosses its end.
const MARGIN: usize = 1 << 20;

// bits

struct Bits<'a> {
    data: &'a [u8],
    /// Next byte to load.
    pos: usize,
    buf: u64,
    n: u32,
}

impl<'a> Bits<'a> {
    fn at(data: &'a [u8], bit: usize) -> Self {
        let mut bits = Self {
            data,
            pos: bit / 8,
            buf: 0,
            n: 0,
        };
        bits.refill();
        bits.consume((bit % 8) as u32);
        bits
    }

    /// Ensure at least 56 bits. Past the end it reads zeros; `overrun` tells.
    #[inline(always)]
    fn refill(&mut self) {
        if self.pos + 8 <= self.data.len() {
            let v = u64::from_le_bytes(self.data[self.pos..self.pos + 8].try_into().unwrap());
            // Bytes beyond those counted are loaded again at the same bit
            // positions next time, so or-ing them in early is harmless.
            self.buf |= v << self.n;
            let take = (63 - self.n) / 8;
            self.pos += take as usize;
            self.n += take * 8;
        } else {
            while self.n <= 56 {
                let b = self.data.get(self.pos).copied().unwrap_or(0);
                self.buf |= u64::from(b) << self.n;
                self.pos += 1;
                self.n += 8;
            }
        }
    }

    #[inline(always)]
    fn consume(&mut self, k: u32) {
        self.buf >>= k;
        self.n -= k;
    }

    #[inline(always)]
    fn take(&mut self, k: u32) -> u32 {
        if self.n < k {
            self.refill();
        }
        let v = (self.buf & ((1u64 << k) - 1)) as u32;
        self.consume(k);
        v
    }

    fn bit_pos(&self) -> usize {
        self.pos * 8 - self.n as usize
    }

    fn overrun(&self) -> bool {
        self.bit_pos() > self.data.len() * 8
    }
}

// huffman

const PRIMARY: u32 = 10;

/// Two-level decode table. A primary entry is `symbol << 8 | length`, or
/// `offset << 8 | 0x80 | sub_bits` pointing at a subtable for longer codes,
/// whose entries are `symbol << 8 | (length - PRIMARY)`. Length 0 is invalid.
struct Huffman {
    table: Vec<u32>,
}

impl Huffman {
    /// Canonical code for `lens`, or `None` if it is over-subscribed, or
    /// incomplete with more than a single one-bit code (zlib's rule).
    fn new(lens: &[u8]) -> Option<Self> {
        let mut count = [0u16; 16];
        for &l in lens {
            count[usize::from(l)] += 1;
        }
        count[0] = 0;
        let mut left = 1i32;
        let mut max = 0;
        for (len, &c) in count.iter().enumerate().skip(1) {
            left = (left << 1) - i32::from(c);
            if left < 0 {
                return None;
            }
            if c > 0 {
                max = len;
            }
        }
        if left > 0 && max > 1 {
            return None;
        }
        let mut next = [0u32; 16];
        let mut code = 0u32;
        for len in 1..16 {
            code = (code + u32::from(count[len - 1])) << 1;
            next[len] = code;
        }
        let primary = 1usize << PRIMARY;
        let mut table = vec![0u32; primary];
        let reverse = |code: u32, len: usize| code.reverse_bits() >> (32 - len);
        let mut long = Vec::new();
        for (sym, &l) in lens.iter().enumerate() {
            let len = usize::from(l);
            if len == 0 {
                continue;
            }
            let rev = reverse(next[len], len) as usize;
            next[len] += 1;
            if len <= PRIMARY as usize {
                for j in (rev..primary).step_by(1 << len) {
                    table[j] = (sym as u32) << 8 | len as u32;
                }
            } else {
                long.push((sym, len, rev));
            }
        }
        if !long.is_empty() {
            let mut sub_len = vec![0usize; primary];
            for &(_, len, rev) in &long {
                let p = rev & (primary - 1);
                sub_len[p] = sub_len[p].max(len - PRIMARY as usize);
            }
            for (p, &bits) in sub_len.iter().enumerate() {
                if bits > 0 {
                    table[p] = (table.len() as u32) << 8 | 0x80 | bits as u32;
                    table.resize(table.len() + (1 << bits), 0);
                }
            }
            for (sym, len, rev) in long {
                let p = rev & (primary - 1);
                let bits = sub_len[p];
                let offset = (table[p] >> 8) as usize;
                let extra = len - PRIMARY as usize;
                for j in ((rev >> PRIMARY)..(1 << bits)).step_by(1 << extra) {
                    table[offset + j] = (sym as u32) << 8 | extra as u32;
                }
            }
        }
        Some(Self { table })
    }

    /// Needs at least 15 bits in `bits`.
    #[inline(always)]
    fn decode(&self, bits: &mut Bits) -> Option<u32> {
        let e = self.table[(bits.buf & ((1 << PRIMARY) - 1)) as usize];
        if e & 0x80 == 0 {
            let len = e & 15;
            if len == 0 {
                return None;
            }
            bits.consume(len);
            Some(e >> 8)
        } else {
            let sub = ((bits.buf >> PRIMARY) & ((1 << (e & 15)) - 1)) as usize;
            let e2 = self.table[(e >> 8) as usize + sub];
            let len = e2 & 15;
            if len == 0 {
                return None;
            }
            bits.consume(PRIMARY + len);
            Some(e2 >> 8)
        }
    }
}

const LEN_BASE: [u16; 29] = [
    3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131,
    163, 195, 227, 258,
];
const LEN_EXTRA: [u8; 29] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0,
];
const DIST_BASE: [u16; 30] = [
    1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537,
    2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,
];
const DIST_EXTRA: [u8; 30] = [
    0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13,
    13,
];
const CODE_ORDER: [usize; 19] = [
    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
];

fn fixed_tables() -> (Huffman, Huffman) {
    let mut lens = [0u8; 288];
    lens[..144].fill(8);
    lens[144..256].fill(9);
    lens[256..280].fill(7);
    lens[280..].fill(8);
    (
        Huffman::new(&lens).unwrap(),
        Huffman::new(&[5u8; 32]).unwrap(),
    )
}

/// The litlen and distance codes of a dynamic block, after its 3 header bits.
fn dynamic_tables(bits: &mut Bits) -> Option<(Huffman, Huffman)> {
    let hlit = bits.take(5) as usize + 257;
    let hdist = bits.take(5) as usize + 1;
    let hclen = bits.take(4) as usize + 4;
    if hlit > 286 || hdist > 30 {
        return None;
    }
    let mut code_lens = [0u8; 19];
    for &i in &CODE_ORDER[..hclen] {
        code_lens[i] = bits.take(3) as u8;
    }
    let codes = Huffman::new(&code_lens)?;
    let mut lens = [0u8; 286 + 30];
    let total = hlit + hdist;
    let mut i = 0;
    while i < total {
        bits.refill();
        let sym = codes.decode(bits)?;
        let (value, repeat) = match sym {
            0..=15 => (sym as u8, 1),
            16 => (*lens[..i].last()?, 3 + bits.take(2) as usize),
            17 => (0, 3 + bits.take(3) as usize),
            _ => (0, 11 + bits.take(7) as usize),
        };
        if i + repeat > total {
            return None;
        }
        lens[i..i + repeat].fill(value);
        i += repeat;
    }
    if lens[256] == 0 {
        return None;
    }
    Some((
        Huffman::new(&lens[..hlit])?,
        Huffman::new(&lens[hlit..total])?,
    ))
}

// inflate

/// Output symbols: bytes, or (in `u16`) placeholders for the unknown window.
trait Symbol: Copy {
    fn byte(b: u8) -> Self;
    /// The symbol for window byte `index` (0 is 32 KiB back), if representable.
    fn window(index: usize) -> Option<Self>;
}

impl Symbol for u8 {
    fn byte(b: u8) -> Self {
        b
    }
    fn window(_: usize) -> Option<Self> {
        None
    }
}

impl Symbol for u16 {
    fn byte(b: u8) -> Self {
        u16::from(b)
    }
    fn window(index: usize) -> Option<Self> {
        Some(256 + index as u16)
    }
}

#[derive(Debug)]
enum Inflate {
    /// Invalid deflate data.
    Bad,
    /// Ran past the end of the bytes given.
    More,
}

struct Inflater<'a> {
    bits: Bits<'a>,
    fixed: Option<(Huffman, Huffman)>,
}

impl<'a> Inflater<'a> {
    /// Decode one block into `out[*n..]`, advancing `*n`; returns whether it
    /// was the final block. `out` is used as a buffer with spare length past
    /// `*n` and grows as needed; callers truncate it to `*n` at the end.
    /// A copy reaching before `out[0]` becomes window placeholders if `T`
    /// has them, and is an error otherwise. `last_window` is kept at or past
    /// the end of the last placeholder written.
    fn block<T: Symbol>(
        &mut self,
        out: &mut Vec<T>,
        n: &mut usize,
        last_window: &mut usize,
    ) -> Result<bool, Inflate> {
        let bits = &mut self.bits;
        let last = bits.take(1) == 1;
        let dynamic;
        let (litlen, dist) = match bits.take(2) {
            0 => {
                let skip = bits.n % 8;
                bits.consume(skip);
                let len = bits.take(16) as usize;
                if len as u32 != !bits.take(16) & 0xffff {
                    return Err(Inflate::Bad);
                }
                grow(out, *n + len);
                for slot in &mut out[*n..*n + len] {
                    *slot = T::byte(bits.take(8) as u8);
                }
                *n += len;
                return if bits.overrun() {
                    Err(Inflate::More)
                } else {
                    Ok(last)
                };
            }
            1 => {
                let fixed = self.fixed.get_or_insert_with(fixed_tables);
                (&fixed.0, &fixed.1)
            }
            2 => {
                dynamic = dynamic_tables(bits).ok_or(Inflate::Bad)?;
                (&dynamic.0, &dynamic.1)
            }
            _ => return Err(Inflate::Bad),
        };
        loop {
            // Room for three literals and the longest copy, which may write up
            // to a chunk past its end.
            if *n + 3 + 258 + CHUNK > out.len() {
                grow(out, *n + 3 + 258 + CHUNK);
            }
            // A refill leaves at least 56 bits: three literal codes of up to
            // 15 bits each, or one length, distance and their extra bits.
            bits.refill();
            let mut sym = litlen.decode(bits).ok_or(Inflate::Bad)?;
            if sym < 256 {
                out[*n] = T::byte(sym as u8);
                *n += 1;
                sym = litlen.decode(bits).ok_or(Inflate::Bad)?;
                if sym < 256 {
                    out[*n] = T::byte(sym as u8);
                    *n += 1;
                    sym = litlen.decode(bits).ok_or(Inflate::Bad)?;
                    if sym < 256 {
                        out[*n] = T::byte(sym as u8);
                        *n += 1;
                        // Past the end the bits read as zeros, which can
                        // decode as literals forever: stop there as copies do.
                        if bits.overrun() {
                            return Err(Inflate::More);
                        }
                        continue;
                    }
                }
                if bits.n < 48 {
                    bits.refill();
                }
            }
            if sym == 256 {
                break;
            }
            let s = (sym - 257) as usize;
            if s >= 29 {
                return Err(Inflate::Bad);
            }
            let len = usize::from(LEN_BASE[s]) + bits.take(u32::from(LEN_EXTRA[s])) as usize;
            let d = dist.decode(bits).ok_or(Inflate::Bad)? as usize;
            if d >= 30 {
                return Err(Inflate::Bad);
            }
            let distance = usize::from(DIST_BASE[d]) + bits.take(u32::from(DIST_EXTRA[d])) as usize;
            let p = *n;
            // Copying from before the last placeholder may copy placeholders.
            if p.saturating_sub(distance) < *last_window {
                *last_window = p + len;
            }
            if distance <= p {
                let from = p - distance;
                if distance >= CHUNK {
                    // Copy whole chunks, which compile to vector loads and
                    // stores, and may write past `len` into spare buffer.
                    // Chunks never overlap the bytes they read.
                    let mut i = 0;
                    while i < len {
                        let chunk: [T; CHUNK] = out[from + i..from + i + CHUNK].try_into().unwrap();
                        out[p + i..p + i + CHUNK].copy_from_slice(&chunk);
                        i += CHUNK;
                    }
                } else if distance >= len {
                    out.copy_within(from..from + len, p);
                } else if distance == 1 {
                    let v = out[from];
                    out[p..p + len].fill(v);
                } else {
                    for i in 0..len {
                        out[p + i] = out[from + i];
                    }
                }
            } else {
                for i in 0..len {
                    let at = (p + i) as isize - distance as isize;
                    out[p + i] = if at >= 0 {
                        out[at as usize]
                    } else {
                        T::window((WINDOW as isize + at) as usize).ok_or(Inflate::Bad)?
                    };
                }
                *last_window = p + len;
            }
            *n += len;
            if bits.overrun() {
                return Err(Inflate::More);
            }
        }
        if bits.overrun() {
            return Err(Inflate::More);
        }
        Ok(last)
    }
}

/// Elements per copy step for back-references at least this far back.
const CHUNK: usize = 16;

/// Make `out` at least `need` long. Capacity doubles so growth is amortized,
/// but only a little past `need` is initialized: zero-filling the whole
/// doubled buffer would make all of it resident and write it twice.
fn grow<T: Symbol>(out: &mut Vec<T>, need: usize) {
    if out.len() < need {
        let len = need + (1 << 18);
        if out.capacity() < len {
            out.reserve(len.max(out.capacity() * 2) - out.len());
        }
        out.resize(len, T::byte(0));
    }
}

/// Buffers kept for reuse between regions, so that the decoder does not
/// allocate and free several megabytes per region. Freed large blocks stay
/// resident in most allocators (macOS kept 0.5 to 1 GB of them), and fresh
/// ones fault in new pages. A buffer is allocated only when no idle one is
/// available, so the buffers that exist never outnumber the most that were
/// in use at once; each channel's capacity caps the idle ones.
///
/// Reused buffers keep their length and old contents: the decoder tracks how
/// much it has written and grows only past that, so pages are not zeroed
/// again for every region.
struct Pool {
    /// Compressed input, one per decoding thread.
    input: (
        crossbeam_channel::Sender<Vec<u8>>,
        crossbeam_channel::Receiver<Vec<u8>>,
    ),
    /// Output that may hold window placeholders.
    words: (
        crossbeam_channel::Sender<Vec<u16>>,
        crossbeam_channel::Receiver<Vec<u16>>,
    ),
    /// Decoded output: buffers the reader is done with come back here.
    output: (
        crossbeam_channel::Sender<Vec<u8>>,
        crossbeam_channel::Receiver<Vec<u8>>,
    ),
    /// Capacity for new output and placeholder buffers.
    output_size: Estimate,
    words_size: Estimate,
}

/// How large a buffer must be: the largest size seen so far, or a guess
/// before the first.
struct Estimate {
    largest: std::sync::atomic::AtomicUsize,
    guess: usize,
}

impl Estimate {
    fn new(guess: usize) -> Self {
        Self {
            largest: std::sync::atomic::AtomicUsize::new(0),
            guess,
        }
    }
    fn get(&self) -> usize {
        match self.largest.load(std::sync::atomic::Ordering::Relaxed) {
            0 => self.guess,
            n => n,
        }
    }
    fn saw(&self, n: usize) {
        self.largest
            .fetch_max(n, std::sync::atomic::Ordering::Relaxed);
    }
}

/// A buffer from `idle` with capacity for at least `need` elements. A new
/// one gets 1/8 more, so that later regions rarely outgrow it, and it is
/// reserved exactly rather than grown by doubling, which copies and frees
/// every smaller block on the way (freed blocks stay resident).
fn reuse<T>(idle: &crossbeam_channel::Receiver<Vec<T>>, need: usize) -> Vec<T> {
    let v = idle.try_recv().unwrap_or_default();
    if v.capacity() >= need {
        return v;
    }
    // Not `reserve`: reallocating would copy the old contents.
    drop(v);
    Vec::with_capacity(need + need / 8)
}

impl Pool {
    fn new(plan: Plan, output: Recycle) -> Self {
        Self {
            input: crossbeam_channel::bounded(plan.threads + 2),
            words: crossbeam_channel::bounded(plan.threads + 2),
            output,
            // About a FASTQ region's output (ratio 4 to 5); for regions with
            // less output the rest is only address space, never touched.
            output_size: Estimate::new(plan.region.saturating_mul(5)),
            // Placeholder output is often the whole region: copies of copies
            // keep reaching the unknown window (FASTQ names and qualities).
            words_size: Estimate::new(plan.region.saturating_mul(5)),
        }
    }
    fn input(&self) -> Vec<u8> {
        self.input.1.try_recv().unwrap_or_default()
    }
    /// A placeholder buffer with room for a region's output.
    fn words(&self) -> Vec<u16> {
        reuse(&self.words.1, self.words_size.get())
    }
    /// An output buffer with room for a region's output.
    fn output(&self) -> Vec<u8> {
        reuse(&self.output.1, self.output_size.get())
    }
    fn give_input(&self, v: Vec<u8>) {
        if v.capacity() > 0 {
            let _ = self.input.0.try_send(v);
        }
    }
    fn give_words(&self, v: Vec<u16>) {
        if v.capacity() > 0 {
            let _ = self.words.0.try_send(v);
        }
    }
    fn give_output(&self, v: Vec<u8>) {
        if v.capacity() > 0 {
            let _ = self.output.0.try_send(v);
        }
    }
}

/// Output buffers handed back for reuse: the reader sends on one end, the
/// decoder takes from the other.
pub(crate) type Recycle = (
    crossbeam_channel::Sender<Vec<u8>>,
    crossbeam_channel::Receiver<Vec<u8>>,
);

/// Output buffers that can be in use at once with `ahead` regions in
/// flight: one per region in flight, the output queue's [`OUTPUT_DEPTH`],
/// and the reader's own, plus one to spare. The reader's recycling channel
/// holds this many, so a returned buffer is never freed.
pub(crate) fn output_buffers(ahead: usize) -> usize {
    ahead + OUTPUT_DEPTH + 2
}

/// Regions queued between the region decoder and the reader.
pub(crate) const OUTPUT_DEPTH: usize = 2;

/// A region decoded from a block start up to a block boundary.
struct Decoded {
    start: usize,
    end: usize,
    /// The member ended here (the final block was decoded).
    last: bool,
    /// Output that may hold window placeholders.
    head: Vec<u16>,
    /// The region's whole output: `head.len()` bytes for the head, filled in
    /// by [`resolve`], then the plain bytes decoded after it (the tail).
    /// Empty until needed if there is no tail.
    out: Vec<u8>,
    /// CRC32 of the tail.
    tail_crc: u32,
}

impl Decoded {
    fn tail(&self) -> &[u8] {
        self.out.get(self.head.len()..).unwrap_or_default()
    }
    fn len(&self) -> usize {
        self.head.len() + self.tail().len()
    }
}

/// Decode from bit `start` of `data` until a block boundary at or after bit
/// `stop`, or the final block. `first` means the stream starts at `start`, so
/// there is no window and copies may not reach before it.
///
/// The tail is decoded straight into the output buffer after room for the
/// head, so the region is delivered as one buffer without a copy.
fn decode(
    data: &[u8],
    start: usize,
    stop: usize,
    first: bool,
    pool: &Pool,
) -> Result<Decoded, Inflate> {
    let mut inflater = Inflater {
        bits: Bits::at(data, start),
        fixed: None,
    };
    let mut head = Vec::new();
    // Taken when the tail starts: a region that is all head needs its
    // output buffer only once it can be resolved.
    let mut out = Vec::new();
    if first {
        out = pool.output();
    } else {
        head = pool.words();
    }
    let (mut head_n, mut out_n) = (0, 0);
    let mut in_tail = first;
    let mut last_window = 0;
    let mut last = false;
    while inflater.bits.bit_pos() < stop || inflater.bits.bit_pos() == start {
        let result = if in_tail {
            inflater.block(&mut out, &mut out_n, &mut last_window)
        } else {
            inflater.block(&mut head, &mut head_n, &mut last_window)
        };
        last = match result {
            Ok(last) => last,
            Err(e) => {
                pool.give_words(head);
                pool.give_output(out);
                // Past the end the bits read as zeros, which can form an
                // invalid code before the overrun is noticed: that is
                // missing input.
                return Err(match e {
                    Inflate::Bad if inflater.bits.overrun() => Inflate::More,
                    e => e,
                });
            }
        };
        if !in_tail && head_n - last_window.min(head_n) >= WINDOW {
            // No later copy can reach a placeholder: continue in bytes after
            // the head, whose last 32 KiB are plain bytes and the history.
            out = pool.output();
            grow(&mut out, head_n);
            for (o, &s) in out[head_n - WINDOW..head_n]
                .iter_mut()
                .zip(&head[head_n - WINDOW..head_n])
            {
                *o = s as u8;
            }
            out_n = head_n;
            in_tail = true;
        }
        if last {
            break;
        }
    }
    head.truncate(head_n);
    out.truncate(out_n);
    pool.words_size.saw(head_n);
    pool.output_size.saw(out_n.max(head_n));
    let tail_crc = crc32fast::hash(out.get(head_n..).unwrap_or_default());
    Ok(Decoded {
        start,
        end: inflater.bits.bit_pos(),
        last,
        head,
        out,
        tail_crc,
    })
}

/// The first bit at or after `from` and before `to` where a non-final
/// dynamic block starts whose header is valid and whose whole block decodes.
fn find_block(data: &[u8], from: usize, to: usize, scratch: &mut Vec<u16>) -> Option<usize> {
    let to = to.min(data.len().saturating_sub(8) * 8);
    for bit in from..to {
        let byte = bit / 8;
        let w = u32::from_le_bytes(data[byte..byte + 4].try_into().unwrap()) >> (bit % 8);
        // BFINAL = 0, BTYPE = 2, HLIT <= 29, HDIST <= 29.
        if w & 7 != 4 || (w >> 3) & 31 > 29 || (w >> 8) & 31 > 29 {
            continue;
        }
        // The code-length code must be complete or a single one-bit code.
        let hclen = ((w >> 13) & 15) as usize + 4;
        let mut bits = Bits::at(data, bit + 17);
        let mut lens = [0u8; 19];
        for &i in &CODE_ORDER[..hclen] {
            lens[i] = bits.take(3) as u8;
        }
        if Huffman::new(&lens).is_none() {
            continue;
        }
        let mut inflater = Inflater {
            bits: Bits::at(data, bit),
            fixed: None,
        };
        let (mut n, mut last_window) = (0, 0);
        if inflater.block(scratch, &mut n, &mut last_window).is_ok() {
            return Some(bit);
        }
    }
    None
}

// driver

/// Default compressed bytes per region.
pub(crate) const DEFAULT_REGION: usize = REGION;

/// Decoded output, in order; an empty buffer marks the end.
pub(crate) type Output = crossbeam_channel::Sender<Result<Vec<u8>, Error>>;

/// How to split.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Plan {
    /// Compressed bytes per region (at least 1).
    pub(crate) region: usize,
    /// Decoding workers (at least 1).
    pub(crate) threads: usize,
    /// Regions in flight, from queued for decoding to delivered (at
    /// least 1).
    pub(crate) ahead: usize,
    /// Bytes per output buffer for members decoded sequentially.
    pub(crate) buffer: usize,
}

/// Whether a member starting `remaining` bytes before the end of the input
/// spans enough regions to decode in parallel.
pub(crate) fn worth_it(remaining: u64, region: usize) -> bool {
    remaining >= 3 * region as u64
}

fn read_range_into(
    src: &dyn ReadAt,
    from: usize,
    to: usize,
    buf: &mut Vec<u8>,
) -> Result<(), Error> {
    buf.resize(to - from, 0);
    src.read_exact_at(buf, from as u64).map_err(Error::input)
}

fn read_range(src: &dyn ReadAt, from: usize, to: usize) -> Result<Vec<u8>, Error> {
    let mut buf = Vec::new();
    read_range_into(src, from, to, &mut buf)?;
    Ok(buf)
}

/// Bit offset of the deflate stream of the member at byte `start`.
fn deflate_start(src: &dyn ReadAt, len: usize, start: usize) -> Result<usize, Error> {
    let mut want = 64 << 10;
    loop {
        let to = (start + want).min(len);
        match crate::header::parse(&read_range(src, start, to)?)? {
            Some(n) => return Ok((start + n) * 8),
            None if to == len => return Err(Error::truncated()),
            None => want *= 4,
        }
    }
}

/// Region `k` of the member whose deflate stream starts at bit `deflate`:
/// decode from its first block, or from `start` when given, reading more
/// past the region if the block crossing its end needs it. `Ok(None)` means
/// no block start was found in the region.
fn region(
    src: &dyn ReadAt,
    len: usize,
    size: usize,
    k: usize,
    deflate: usize,
    start: Option<usize>,
    pool: &Pool,
) -> Result<Option<Decoded>, Error> {
    let begin = deflate / 8 + k * size;
    let end = (begin + size).min(len);
    let from = start.map_or(begin, |s| s / 8);
    let mut margin = MARGIN;
    let mut data = pool.input();
    let result = loop {
        let read_to = (end + margin).min(len);
        if let Err(e) = read_range_into(src, from, read_to, &mut data) {
            break Err(e);
        }
        let base = from * 8;
        let local_start = match start {
            Some(s) => s - base,
            None if k == 0 => deflate - base,
            None => {
                let mut scratch = pool.words();
                let found = find_block(&data, 0, (end - from) * 8, &mut scratch);
                pool.give_words(scratch);
                match found {
                    Some(bit) => bit,
                    None => break Ok(None),
                }
            }
        };
        let first = k == 0 && start.is_none_or(|s| s == deflate);
        match decode(&data, local_start, (end - from) * 8, first, pool) {
            Ok(mut decoded) => {
                decoded.start += base;
                decoded.end += base;
                break Ok(Some(decoded));
            }
            Err(Inflate::More) if read_to < len => margin *= 4,
            Err(Inflate::More) => break Err(Error::truncated()),
            Err(Inflate::Bad) => break Err(Error::corrupt("invalid deflate data")),
        }
    };
    pool.give_input(data);
    result
}

/// A resolved region: its dispatch number, the head's CRC32 (or the error
/// resolving it), the head's length, and the region.
type Resolved = (usize, Result<u32, Error>, usize, Decoded);

enum Job {
    Decode(usize),
    /// Fill a decoded region's placeholders from the 32 KiB before it.
    Resolve(usize, Decoded, std::sync::Arc<Vec<u8>>),
}

/// Write the bytes for `head` into `out` (as long as `head`), filling
/// placeholders from `window`, the bytes before it.
fn resolve(head: &[u16], window: &[u8], out: &mut [u8]) -> Result<(), Error> {
    let mut lut = [0u8; 256 + WINDOW];
    for (i, b) in lut[..256].iter_mut().enumerate() {
        *b = i as u8;
    }
    let known = window.len().min(WINDOW);
    lut[256 + WINDOW - known..].copy_from_slice(&window[window.len() - known..]);
    // Near the start of the stream the window is short: a placeholder for a
    // byte before the stream is invalid data.
    if known < WINDOW
        && head
            .iter()
            .any(|&s| s >= 256 && usize::from(s - 256) < WINDOW - known)
    {
        return Err(Error::corrupt("invalid deflate distance"));
    }
    for (o, &s) in out.iter_mut().zip(head) {
        *o = lut[usize::from(s)];
    }
    Ok(())
}

/// Decompress `src` (`len` bytes), sending output in order on `sender`,
/// then an empty buffer at the end, or the first error. The reader hands
/// buffers it is done with back on `recycle`.
///
/// The first member is decoded in parallel regions. Each later member is too
/// while the one before it spanned at least three regions and enough input
/// remains (concatenated files, as ENA serves them); from the first member
/// that does not, the rest is decoded sequentially.
pub(crate) fn decompress(
    src: &dyn ReadAt,
    len: usize,
    plan: Plan,
    sender: &Output,
    recycle: Recycle,
) {
    let pool = Pool::new(plan, recycle);
    let result = (|| -> Result<bool, Error> {
        let mut start = 0;
        loop {
            let Some(end) =
                member(src, len, start, plan, sender, &pool).map_err(|e| e.at(start as u64))?
            else {
                return Ok(false);
            };
            if end == len {
                return Ok(true);
            }
            if end - start >= 3 * plan.region && worth_it((len - end) as u64, plan.region) {
                start = end;
                continue;
            }
            return sequential_rest(src, len, end, plan.buffer, sender, &pool.output.1);
        }
    })();
    match result {
        Ok(true) => {
            let _ = sender.send(Ok(Vec::new()));
        }
        Ok(false) => {}
        Err(e) => {
            let _ = sender.send(Err(e));
        }
    }
}

/// Decode `src[start..]` with [`GzReader`] into buffers of `buffer` bytes,
/// reusing those that come back on `reused`; `false` if the reader went
/// away.
pub(crate) fn sequential_rest(
    src: &dyn ReadAt,
    len: usize,
    start: usize,
    buffer: usize,
    sender: &Output,
    reused: &crossbeam_channel::Receiver<Vec<u8>>,
) -> Result<bool, Error> {
    let cursor = Cursor {
        src,
        pos: start as u64,
        end: len as u64,
    };
    let mut decoder = GzReader::new(BufReader::with_capacity(
        crate::sequential::INPUT_BUFFER,
        cursor,
    ))
    .with_base(start as u64);
    loop {
        let mut chunk = reused.try_recv().unwrap_or_default();
        chunk.resize(buffer, 0);
        let mut filled = 0;
        while filled < chunk.len() {
            match decoder.read_typed(&mut chunk[filled..])? {
                0 => break,
                n => filled += n,
            }
        }
        if filled == 0 {
            return Ok(true);
        }
        chunk.truncate(filled);
        if sender.send(Ok(chunk)).is_err() {
            return Ok(false);
        }
    }
}

/// Decode the member at byte `start` in parallel regions. Returns the
/// offset just past its trailer, or `None` if the reader went away.
fn member(
    src: &dyn ReadAt,
    len: usize,
    start: usize,
    plan: Plan,
    sender: &Output,
    pool: &Pool,
) -> Result<Option<usize>, Error> {
    let Plan {
        region: size,
        threads,
        ahead,
        ..
    } = plan;
    let deflate = deflate_start(src, len, start)?;
    let regions = (len - deflate / 8).div_ceil(size).max(1);
    // Stops the workers early if the reader goes away or the member ends.
    let stopped = std::sync::atomic::AtomicBool::new(false);
    let stopped = &stopped;
    // Set once the member's final block is decoded: skip further regions.
    let ended = std::sync::atomic::AtomicBool::new(false);
    let ended = &ended;
    std::thread::scope(|scope| -> Result<Option<usize>, Error> {
        let (jobs, queue) = crossbeam_channel::bounded::<Job>(2 * ahead);
        // Resolve jobs go first: each holds a region's placeholder output,
        // the largest buffer, until it is done.
        let (resolves, urgent) = crossbeam_channel::unbounded::<Job>();
        let (decoded_tx, decoded_rx) = crossbeam_channel::unbounded();
        let (resolved_tx, resolved_rx) = crossbeam_channel::unbounded();
        for _ in 0..threads {
            let (queue, urgent, decoded_tx, resolved_tx) = (
                queue.clone(),
                urgent.clone(),
                decoded_tx.clone(),
                resolved_tx.clone(),
            );
            scope.spawn(move || {
                let next = || match urgent.try_recv() {
                    Ok(job) => Ok(job),
                    Err(_) => {
                        let mut select = crossbeam_channel::Select::new();
                        select.recv(&urgent);
                        select.recv(&queue);
                        let op = select.select();
                        match op.index() {
                            0 => op.recv(&urgent),
                            _ => op.recv(&queue),
                        }
                    }
                };
                // Both channels close together, at the end of the member.
                while let Ok(job) = next() {
                    if stopped.load(std::sync::atomic::Ordering::Relaxed) {
                        break;
                    }
                    match job {
                        // A failure here may only be a false block start, so
                        // the coordinator decides; it re-decodes from known output.
                        Job::Decode(_) if ended.load(std::sync::atomic::Ordering::Relaxed) => {}
                        Job::Decode(k) => {
                            let _ = decoded_tx
                                .send((k, region(src, len, size, k, deflate, None, pool)));
                        }
                        Job::Resolve(k, mut decoded, window) => {
                            let n = decoded.head.len();
                            if decoded.out.len() < n {
                                decoded.out = pool.output();
                                decoded.out.resize(n, 0);
                                decoded.out.truncate(n);
                            }
                            let head = &mut decoded.out[..n];
                            let crc = resolve(&decoded.head, &window, head)
                                .map(|()| crc32fast::hash(head));
                            // Keep only the length; the buffer is free now.
                            let head = std::mem::take(&mut decoded.head);
                            let _ = resolved_tx.send((k, crc, n, decoded));
                            pool.give_words(head);
                        }
                    }
                }
            });
        }
        drop((queue, decoded_tx, resolved_tx));
        let stop = || stopped.store(true, std::sync::atomic::Ordering::Relaxed);
        let workers_stopped = || Error::internal("decompression workers stopped");
        let mut queued = 0;
        let mut decoded_pending = std::collections::BTreeMap::new();
        let mut resolved_pending = std::collections::BTreeMap::new();
        let mut window: Vec<u8> = Vec::with_capacity(2 * WINDOW);
        let mut crc = crc32fast::Hasher::new();
        let mut total = 0u64;
        let mut pos = deflate;
        // Resolve jobs are numbered in dispatch order: skipped regions have
        // none; `next_send` is the next of them to send.
        let (mut k, mut dispatched, mut next_send) = (0, 0, 0);
        let mut finished = None;
        // Take a resolved region and send what is now in order; `false` if
        // the reader went away.
        let deliver = |item: Resolved,
                       next_send: &mut usize,
                       resolved_pending: &mut std::collections::BTreeMap<_, _>,
                       crc: &mut crc32fast::Hasher|
         -> Result<bool, Error> {
            let (i, head_crc, head_len, decoded) = item;
            resolved_pending.insert(i, (head_crc, head_len, decoded));
            while let Some((head_crc, head_len, decoded)) = resolved_pending.remove(next_send) {
                let head_crc: u32 = head_crc?;
                let decoded: Decoded = decoded;
                crc.combine(&crc32fast::Hasher::new_with_initial_len(
                    head_crc,
                    head_len as u64,
                ));
                crc.combine(&crc32fast::Hasher::new_with_initial_len(
                    decoded.tail_crc,
                    (decoded.out.len() - head_len) as u64,
                ));
                if decoded.out.is_empty() {
                    pool.give_output(decoded.out);
                } else if sender.send(Ok(decoded.out)).is_err() {
                    return Ok(false);
                }
                *next_send += 1;
            }
            Ok(true)
        };
        macro_rules! deliver {
            ($item:expr) => {
                if !deliver($item, &mut next_send, &mut resolved_pending, &mut crc)? {
                    stop();
                    return Ok(None);
                }
            };
        }
        // At most `ahead` regions are in flight, from queued to sent: queued
        // from `k` on, and taken but not yet sent. Regions decoded out of
        // order wait in `decoded_pending` within that bound.
        macro_rules! top_up {
            () => {
                while queued < regions && queued - k + dispatched - next_send < ahead {
                    jobs.send(Job::Decode(queued)).ok();
                    queued += 1;
                }
            };
        }
        while k < regions && finished.is_none() {
            top_up!();
            if queued == k {
                // Every slot is taken by regions being resolved.
                deliver!(resolved_rx.recv().map_err(|_| workers_stopped())?);
                continue;
            }
            let result = loop {
                if let Some(r) = decoded_pending.remove(&k) {
                    break r;
                }
                crossbeam_channel::select! {
                    recv(decoded_rx) -> r => {
                        let (i, r) = r.map_err(|_| workers_stopped())?;
                        decoded_pending.insert(i, r);
                    }
                    recv(resolved_rx) -> r => {
                        deliver!(r.map_err(|_| workers_stopped())?);
                        top_up!();
                    }
                }
            };
            let region_end = (deflate / 8 + (k + 1) * size).min(len) * 8;
            if pos >= region_end {
                // An earlier region's last block covered this one.
                if let Ok(Some(d)) = result {
                    pool.give_words(d.head);
                    pool.give_output(d.out);
                }
                k += 1;
                continue;
            }
            let decoded = match result {
                Ok(Some(d)) if d.start == pos => d,
                other => {
                    // A false start: its buffers are reused for the redo.
                    if let Ok(Some(d)) = other {
                        pool.give_words(d.head);
                        pool.give_output(d.out);
                    }
                    region(src, len, size, k, deflate, Some(pos), pool)?
                        .expect("a region decoded from a known start always has a block")
                }
            };
            // The next region's window is this region's last 32 KiB. Only that
            // much of the head is resolved here; workers resolve the rest.
            let need = WINDOW.saturating_sub(decoded.tail().len());
            let head_end = &decoded.head[decoded.head.len().saturating_sub(need)..];
            let mut head_bytes = vec![0; head_end.len()];
            resolve(head_end, &window, &mut head_bytes)?;
            let this_window = std::sync::Arc::new(window.clone());
            for part in [&head_bytes[..], decoded.tail()] {
                window.extend_from_slice(&part[part.len().saturating_sub(WINDOW)..]);
                if window.len() > WINDOW {
                    window.drain(..window.len() - WINDOW);
                }
            }
            total += decoded.len() as u64;
            pos = decoded.end;
            if decoded.last {
                finished = Some(pos);
                // Regions queued past the member's end are not needed.
                ended.store(true, std::sync::atomic::Ordering::Relaxed);
            }
            resolves
                .send(Job::Resolve(dispatched, decoded, this_window))
                .ok();
            dispatched += 1;
            k += 1;
            while let Ok(item) = resolved_rx.try_recv() {
                deliver!(item);
            }
        }
        while next_send < dispatched {
            deliver!(resolved_rx.recv().map_err(|_| workers_stopped())?);
        }
        stop();
        drop((jobs, resolves));
        let end = finished.ok_or_else(Error::truncated)?;
        let trailer_at = end.div_ceil(8);
        if trailer_at + 8 > len {
            return Err(Error::truncated());
        }
        let trailer = read_range(src, trailer_at, trailer_at + 8)?;
        if u32::from_le_bytes(trailer[..4].try_into().unwrap()) != crc.finalize()
            || u32::from_le_bytes(trailer[4..].try_into().unwrap()) != total as u32
        {
            return Err(Error::checksum(
                "gzip member CRC32 or length does not match its trailer",
            ));
        }
        Ok(Some(trailer_at + 8))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn fastq(n: usize, seed: u64) -> Vec<u8> {
        let mut state = seed;
        let mut next = |m: u64| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) % m
        };
        let mut out = Vec::new();
        for i in 0..n {
            let len = 50 + next(100) as usize;
            writeln!(out, "@read.{i} length={len}").unwrap();
            out.extend((0..len).map(|_| b"ACGTN"[next(if i % 50 == 0 { 5 } else { 4 }) as usize]));
            out.extend_from_slice(b"\n+\n");
            out.extend((0..len).map(|_| b"#,:FFFF"[next(7) as usize]));
            out.push(b'\n');
        }
        out
    }

    /// Decode `gz` by regions of `region` bytes on `threads` workers.
    fn pgzip(gz: &[u8], region: usize, threads: usize) -> Result<Vec<u8>, Error> {
        let (sender, receiver) = crossbeam_channel::unbounded();
        let plan = Plan {
            region,
            threads,
            ahead: threads + 2,
            buffer: 1 << 20,
        };
        decompress(&gz, gz.len(), plan, &sender, crossbeam_channel::bounded(1));
        drop(sender);
        let mut out = Vec::new();
        for part in receiver {
            let part = part?;
            if part.is_empty() {
                return Ok(out);
            }
            out.extend(part);
        }
        Err(Error::internal("no end marker"))
    }

    fn zlib_gzip(data: &[u8], level: u32) -> Vec<u8> {
        let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(level));
        e.write_all(data).unwrap();
        e.finish().unwrap()
    }

    #[test]
    fn regions_reassemble_exactly() {
        // Many regions at every level, including stored (0) and fast (1).
        let data = fastq(20_000, 3);
        for level in [0, 1, 6, 9] {
            let gz = zlib_gzip(&data, level);
            let region = 64 << 10;
            assert!(gz.len() > 8 * region, "level {level} spans regions");
            assert_eq!(pgzip(&gz, region, 4).unwrap(), data, "level {level}");
        }
        let mut c = libdeflater::Compressor::new(libdeflater::CompressionLvl::new(4).unwrap());
        let mut gz = vec![0; c.gzip_compress_bound(data.len())];
        let n = c.gzip_compress(&data, &mut gz).unwrap();
        gz.truncate(n);
        assert_eq!(pgzip(&gz, 64 << 10, 3).unwrap(), data, "libdeflate");
    }

    #[test]
    fn damaged_stream_fails_instead_of_decoding_zeros() {
        // Found by fuzzing: a 40-byte damaged member whose block runs past
        // the end of the input. The bits read there are zeros, which decoded
        // as literals forever, until the decoder ran out of memory.
        let gz = [
            0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0x03, 0x08, 0, 0x46, 0x50, 0x04, 0, 0x28, 0, 0,
            0, 0x05, 0xc1, 0x81, 0x08, 0, 0, 0, 0, 0xa0, 0xfd, 0xa9, 0x2f, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        assert!(pgzip(&gz, 64, 1).is_err());
    }

    #[test]
    fn later_members_and_errors() {
        let data = fastq(20_000, 5);
        let region = 64 << 10;
        let mut gz = zlib_gzip(&data, 1);
        gz.extend(zlib_gzip(b"second member\n", 6));
        assert_eq!(
            pgzip(&gz, region, 4).unwrap(),
            [&data[..], b"second member\n"].concat()
        );
        // Two large members: the second is split too.
        let mut two = zlib_gzip(&data, 1);
        two.extend(zlib_gzip(&data, 6));
        assert_eq!(
            pgzip(&two, region, 4).unwrap(),
            [&data[..], &data[..]].concat()
        );
        let gz = zlib_gzip(&data, 1);
        let mut corrupt = gz.clone();
        let at = corrupt.len() / 2;
        corrupt[at] ^= 0x10;
        assert!(pgzip(&corrupt, region, 4).is_err(), "corrupt");
        let e = pgzip(&gz[..gz.len() - 3], region, 4).unwrap_err();
        assert_eq!(e.kind(), crate::ErrorKind::Truncated, "truncated trailer");
        let e = pgzip(&gz[..gz.len() / 2], region, 4).unwrap_err();
        assert_eq!(e.kind(), crate::ErrorKind::Truncated, "truncated stream");
        let mut bad_crc = gz.clone();
        let n = bad_crc.len();
        bad_crc[n - 8] ^= 1;
        let e = pgzip(&bad_crc, region, 4).unwrap_err();
        assert_eq!(e.kind(), crate::ErrorKind::Checksum, "bad crc");
        assert_eq!(e.offset(), Some(0));
        let mut bad_second = two.clone();
        let n = bad_second.len();
        bad_second[n - 1] ^= 1;
        let e = pgzip(&bad_second, region, 4).unwrap_err();
        assert_eq!(
            e.kind(),
            crate::ErrorKind::Checksum,
            "bad isize in member 2"
        );
        assert_eq!(e.offset(), Some(gz.len() as u64));
    }
}
