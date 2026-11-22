# readkit

zero-copy batches of sequencing reads, and the stages that produce and consume them

stages in one program pass reads as owned, columnar batches. no FASTQ text in between and no re-parsing. batches use arrow's layout, so pyarrow, polars or duckdb can read them without copying

## install

```toml
[dependencies]
readkit = { git = "https://github.com/drbh/readkit.git" }
```

## usage

```rust
use readkit::fastq::{FastqInput, FastqOptions, FastqSink, FastqSource};

let r1 = b"@a/1\nACGT\n+\nIIII\n".as_slice();
let r2 = b"@a/2\nTTGG\n+\nIIII\n".as_slice();
let source = FastqSource::paired(FastqInput::new(r1), FastqInput::new(r2), FastqOptions::default());
let mut sink = FastqSink::split(vec![Vec::new(), Vec::new()]);
let records = readkit::pump(source, &mut sink)?;
```

anything that is a `BatchSource` can be chained, and each stage reads at most one batch ahead

`fastq`, `gz` and `parallel` are on by default. `fetch`, `sra`, `resolve`, `ngs` and `arrow` are opt in

the core is pure rust and builds for wasm32 with `--no-default-features --features fastq,gz`
