# fuzzing

fuzz targets for the gzip reader and writer. each one compares against flate2 or libdeflate, so a target fails on wrong bytes and not only on crashes

- `fqdeflate_roundtrip` the FASTQ-aware encoder's output decodes back to its input
- `pgzip_differential` region decoding (`Split::Always`) gives flate2's bytes
- `pgzip_corrupt` damaged streams error or give flate2's bytes, and never panic or hang
- `gzip_members` `GzReader` accepts exactly the streams flate2 accepts
- `parallel_members` `ParGzReader` accepts exactly what flate2 accepts, on sized and unsized members

inputs are shaped like FASTQ by `fuzz_targets/common.rs` so the targets reach the FASTQ-specific paths

this directory is its own workspace. run a target with nightly rust and https://github.com/rust-fuzz/cargo-fuzz

```bash
cargo install cargo-fuzz
cd crates/readkit-gz/fuzz
cargo +nightly fuzz run -O pgzip_differential -- -max_total_time=600 -max_len=65536 -timeout=20 \
  -rss_limit_mb=8192 -malloc_limit_mb=2048
```

note: the rss limit is loose since the decoder targets start threads on every input, and on macOS AddressSanitizer's per-thread bookkeeping grows by about a gigabyte a minute (the decoder itself does not leak)

a failing input is saved under `artifacts/<target>/`. replay it with `FUZZ_DUMP` to see the gzip stream it produced

```bash
FUZZ_DUMP=/tmp/stream.gz cargo +nightly fuzz run -O pgzip_corrupt artifacts/pgzip_corrupt/<file>
```
