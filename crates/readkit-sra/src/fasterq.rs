//! The `fasterq-dump` subprocess backend; see [`FasterqDumpBackend`].

use std::ffi::OsString;
use std::io::{BufReader, Read as _};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, ExitStatus, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};
use std::thread::JoinHandle;

use readkit_batch::{Batch, Single};
use readkit_fastq::{FastqInput, FastqOptions, FastqSource};
use readkit_io::BatchSource;

use crate::backend::{OpenRequest, RunInfo, Spot, SpotReader, SraBackend};
use crate::error::{Error, ErrorKind};
use crate::read_type::ReadType;

/// How much of a tool's stderr is kept for error messages.
const STDERR_KEEP: usize = 16 << 10;
/// Lines of stderr quoted in an error message.
const STDERR_LINES: usize = 6;

/// The `fasterq-dump` subprocess backend: the SRA Toolkit's own dumper, its
/// FASTQ output parsed with `readkit-fastq` and regrouped into spots.
///
/// - **Command.** `fasterq-dump --split-spot --stdout --seq-defline
///   '@$si $sn' --qual-defline '+' -t <private temp dir> [-e N]
///   [--include-technical] [extra args] <path>`. `--split-spot` prints every
///   read of a spot as its own record, consecutively, so a spot's reads
///   stay together whatever the layout. The custom defline carries exactly
///   what the backend needs: the spot id (`$si`) and the stored name (`$sn`,
///   the spot id when the run stores none). `$sn` is only requested when
///   names are.
/// - **Spot count.** FASTQ does not say how many spots a run has, and the
///   [`Decoder`](crate::Decoder) needs it at open, so the backend first
///   runs `vdb-dump --id_range <path>` (first row and row count of the
///   `SEQUENCE` table). `fasterq-dump` itself starts on the first spot read.
/// - **Grouping.** Consecutive records with the same spot id form one
///   spot; spot ids must increase. A spot id that `fasterq-dump` skips (all
///   its reads empty or filtered) is served with as many empty reads as the
///   run's first printed spot has, which is what the VDB backend serves for
///   an empty single-end read. Records beyond the run's last spot are
///   [`ErrorKind::Corrupt`].
/// - **Not served:** read types (FASTQ has none; every read is reported
///   biological) and the run's quality declaration (`QualityKind::Unknown`,
///   so the decoder samples qualities for SRA Lite).
/// - **Process lifetime.** The child runs in its own process group, because
///   the SRA Toolkit's `fasterq-dump` is a driver that forks the real
///   dumper. Dropping the reader, or any error, kills the whole group with
///   `SIGKILL`, reaps the child and removes the private temp dir.
///
/// ```no_run
/// use std::sync::Arc;
/// use readkit_io::BatchSource;
/// use readkit_sra::{DecodeOptions, FasterqDumpBackend, SingleDecoder};
///
/// let backend = FasterqDumpBackend::new().threads(4).temp_dir("/scratch");
/// let mut d = SingleDecoder::open_with(Arc::new(backend), "SRR12999707", DecodeOptions::default())?;
/// while let Some(batch) = d.next_batch()? {
///     let _ = batch.column().seqs();
/// }
/// # Ok::<(), readkit_sra::Error>(())
/// ```
#[derive(Clone, Debug)]
pub struct FasterqDumpBackend {
    program: PathBuf,
    vdb_dump: Option<PathBuf>,
    threads: Option<usize>,
    temp_dir: Option<PathBuf>,
    include_technical: bool,
    args: Vec<OsString>,
}

impl Default for FasterqDumpBackend {
    fn default() -> Self {
        Self::new()
    }
}

impl FasterqDumpBackend {
    /// `fasterq-dump` and `vdb-dump` from `PATH`, the tool's default thread
    /// count, temp files under [`std::env::temp_dir`], biological reads only.
    pub fn new() -> Self {
        Self {
            program: PathBuf::from("fasterq-dump"),
            vdb_dump: None,
            threads: None,
            temp_dir: None,
            include_technical: false,
            args: Vec::new(),
        }
    }

    /// The `fasterq-dump` executable: a bare name is looked up on `PATH`.
    pub fn program(mut self, program: impl Into<PathBuf>) -> Self {
        self.program = program.into();
        self
    }

    /// The `vdb-dump` executable used to count spots. Default: `vdb-dump`
    /// next to [`program`](Self::program) when that is a path and the file
    /// exists there, else `vdb-dump` on `PATH`.
    pub fn vdb_dump(mut self, program: impl Into<PathBuf>) -> Self {
        self.vdb_dump = Some(program.into());
        self
    }

    /// `fasterq-dump -e <n>` (the tool's default is 6).
    pub fn threads(mut self, n: usize) -> Self {
        self.threads = Some(n.max(1));
        self
    }

    /// The directory under which each reader creates its private temp dir
    /// (`readkit-fasterq-<pid>-<n>`, passed as `-t` and used as the child's
    /// working directory, removed when the reader is dropped). It must
    /// exist. Default [`std::env::temp_dir`].
    pub fn temp_dir(mut self, dir: impl Into<PathBuf>) -> Self {
        self.temp_dir = Some(dir.into());
        self
    }

    /// `--include-technical`: print technical reads too. FASTQ carries no
    /// read type, so they are reported biological like every other read,
    /// and every [`SegmentSelection`](crate::SegmentSelection) keeps them.
    /// Empty technical reads (the common case) are still not printed.
    pub fn include_technical(mut self, on: bool) -> Self {
        self.include_technical = on;
        self
    }

    /// Append one argument before the run path (for example `-m` and a
    /// memory limit). Arguments that change the record format or order
    /// (`--fasta`, `--split-3`, `--seq-defline`, ...) break regrouping.
    pub fn arg(mut self, arg: impl Into<OsString>) -> Self {
        self.args.push(arg.into());
        self
    }

    /// Append several arguments; see [`arg`](Self::arg).
    pub fn args<I, S>(mut self, args: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<OsString>,
    {
        self.args.extend(args.into_iter().map(Into::into));
        self
    }

    /// Whether the `fasterq-dump` and `vdb-dump` executables exist (on
    /// `PATH` for bare names). Does not run them.
    pub fn is_available(&self) -> bool {
        find_program(&self.program).is_some() && find_program(&self.vdb_dump_program()).is_some()
    }

    fn vdb_dump_program(&self) -> PathBuf {
        if let Some(p) = &self.vdb_dump {
            return p.clone();
        }
        if self.program.components().count() > 1
            && let Some(dir) = self.program.parent()
        {
            let sibling = dir.join("vdb-dump");
            if sibling.is_file() {
                return sibling;
            }
        }
        PathBuf::from("vdb-dump")
    }

    /// Open `path` as a concrete reader (the trait's
    /// [`open`](SraBackend::open) boxes this).
    pub fn open_reader(
        &self,
        path: &Path,
        request: &OpenRequest,
    ) -> Result<FasterqDumpReader, Error> {
        let input = path.display().to_string();
        std::fs::metadata(path).map_err(|e| Error::from(e).with_input(input.clone()))?;
        // `fasterq-dump` runs in its temp dir, where a relative path does not
        // exist and a bare one is looked up as an accession.
        let path =
            &std::path::absolute(path).map_err(|e| Error::from(e).with_input(input.clone()))?;
        let (first, spots) = self.id_range(path)?;
        let mut info = RunInfo::new(first, spots);
        info.has_names = request.names;
        Ok(FasterqDumpReader {
            backend: self.clone(),
            path: path.to_path_buf(),
            input,
            names: request.names,
            info,
            stream: None,
            proc: None,
            temp: None,
            cur: SpotBuf::default(),
            cur_id: None,
            ahead: SpotBuf::default(),
            ahead_full: false,
            last_id: None,
            template: None,
            failed: None,
        })
    }

    /// `vdb-dump --id_range`: (first row, row count).
    fn id_range(&self, path: &Path) -> Result<(u64, u64), Error> {
        let program = self.vdb_dump_program();
        let input = path.display().to_string();
        let out = Command::new(&program)
            .arg("--id_range")
            .arg(path)
            .stdin(Stdio::null())
            .output()
            .map_err(|e| spawn_error(&program, e).with_input(input.clone()))?;
        let stderr = String::from_utf8_lossy(&out.stderr);
        if !out.status.success() {
            return Err(tool_error(&program, out.status, &stderr).with_input(input));
        }
        let stdout = String::from_utf8_lossy(&out.stdout);
        parse_id_range(&stdout).ok_or_else(|| {
            Error::new(
                ErrorKind::Backend,
                format!(
                    "unexpected `{} --id_range` output: {:?}",
                    program.display(),
                    stdout.trim()
                ),
            )
            .with_input(input)
        })
    }

    /// The full `fasterq-dump` command for `path`, temp dir `tmp`.
    fn command(&self, path: &Path, tmp: &Path, names: bool) -> Command {
        let mut cmd = Command::new(&self.program);
        cmd.arg("--split-spot")
            .arg("--stdout")
            .arg("--seq-defline")
            .arg(if names { "@$si $sn" } else { "@$si" })
            .arg("--qual-defline")
            .arg("+")
            .arg("-t")
            .arg(tmp);
        if let Some(n) = self.threads {
            cmd.arg("-e").arg(n.to_string());
        }
        if self.include_technical {
            cmd.arg("--include-technical");
        }
        cmd.args(&self.args)
            .arg(path)
            .current_dir(tmp)
            .stdin(Stdio::null())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped());
        #[cfg(unix)]
        {
            use std::os::unix::process::CommandExt as _;
            cmd.process_group(0);
        }
        cmd
    }
}

impl SraBackend for FasterqDumpBackend {
    fn name(&self) -> &'static str {
        "fasterq-dump"
    }

    fn open(&self, path: &Path, request: &OpenRequest) -> Result<Box<dyn SpotReader>, Error> {
        Ok(Box::new(self.open_reader(path, request)?))
    }
}

/// `program` on its own path, or found on `PATH` for a bare name.
fn find_program(program: &Path) -> Option<PathBuf> {
    if program.components().count() > 1 {
        return program.is_file().then(|| program.to_path_buf());
    }
    let path = std::env::var_os("PATH")?;
    std::env::split_paths(&path)
        .map(|d| d.join(program))
        .find(|p| p.is_file())
}

/// Parse `id-range: first-row = 1, row-count = 89,234`.
fn parse_id_range(text: &str) -> Option<(u64, u64)> {
    let number = |key: &str| -> Option<u64> {
        let rest = &text[text.find(key)? + key.len()..];
        let rest = rest.trim_start().strip_prefix('=')?.trim_start();
        let digits: String = rest
            .chars()
            .take_while(|c| c.is_ascii_digit() || *c == ',')
            .filter(char::is_ascii_digit)
            .collect();
        digits.parse().ok()
    };
    Some((number("first-row")?, number("row-count")?))
}

fn spawn_error(program: &Path, e: std::io::Error) -> Error {
    if e.kind() == std::io::ErrorKind::NotFound {
        Error::new(
            ErrorKind::Unsupported,
            format!(
                "`{}` not found: install the SRA Toolkit or set FasterqDumpBackend::program",
                program.display()
            ),
        )
        .with_source(e)
    } else {
        Error::new(
            ErrorKind::Io,
            format!("starting `{}`: {e}", program.display()),
        )
        .with_source(e)
    }
}

/// The error for a tool that exited unsuccessfully.
fn tool_error(program: &Path, status: ExitStatus, stderr: &str) -> Error {
    let tail = stderr_tail(stderr);
    let mut msg = format!("`{}` failed ({status})", program.display());
    if !tail.is_empty() {
        msg.push_str(": ");
        msg.push_str(&tail);
    }
    Error::new(classify_stderr(stderr), msg)
}

/// The last few non-empty stderr lines, joined with ` | `.
fn stderr_tail(stderr: &str) -> String {
    let lines: Vec<&str> = stderr
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.chars().all(|c| c == '='))
        .collect();
    lines[lines.len().saturating_sub(STDERR_LINES)..].join(" | ")
}

/// Classify a failed SRA Toolkit run by its stderr.
///
/// This is the only place the backend reads tool text, and it reads two
/// things:
///
/// - the state field of the last `RC(module,target,context,object,state)`
///   that the tools print for VDB status codes (the same field
///   [`rc_kind`](crate::rc) maps for the VDB backend);
/// - two messages the tools print without an `RC(...)` when the file is not
///   a VDB object they can dump: vdb-dump's "cannot be opened as
///   vdb-database or vdb-table" and fasterq-dump's "is missing the
///   QUALITY-column".
///
/// Anything else is [`ErrorKind::Backend`].
pub(crate) fn classify_stderr(stderr: &str) -> ErrorKind {
    if let Some(at) = stderr.rfind("RC(") {
        let inner = &stderr[at + 3..];
        if let Some(end) = inner.find(')')
            && let Some(state) = inner[..end].split(',').nth(4)
        {
            return rc_state_kind(state.trim());
        }
    }
    if stderr.contains("cannot be opened as vdb-database or vdb-table")
        || stderr.contains("is missing the QUALITY-column")
    {
        return ErrorKind::Unsupported;
    }
    ErrorKind::Backend
}

/// A VDB `RCState` name to an error kind.
fn rc_state_kind(state: &str) -> ErrorKind {
    match state {
        "rcNotFound" => ErrorKind::NotFound,
        "rcCorrupt" | "rcIncorrect" | "rcInconsistent" | "rcIncomplete" | "rcTooShort"
        | "rcInvalid" => ErrorKind::Corrupt,
        "rcUnsupported" | "rcBadVersion" | "rcWrongType" => ErrorKind::Unsupported,
        "rcUnauthorized" | "rcNoPerm" | "rcExhausted" => ErrorKind::Io,
        _ => ErrorKind::Backend,
    }
}

/// The reader's private temp dir, removed on drop.
#[derive(Debug)]
struct TempDir(PathBuf);

impl TempDir {
    fn create(root: &Path) -> Result<Self, Error> {
        static SEQ: AtomicU64 = AtomicU64::new(0);
        let seq = SEQ.fetch_add(1, Ordering::Relaxed);
        let dir = root.join(format!("readkit-fasterq-{}-{seq}", std::process::id()));
        std::fs::create_dir(&dir).map_err(|e| {
            Error::new(
                ErrorKind::Io,
                format!("creating temp dir {}: {e}", dir.display()),
            )
            .with_source(e)
        })?;
        Ok(Self(dir))
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A running `fasterq-dump` and its stderr collector.
struct Proc {
    child: Child,
    program: PathBuf,
    stderr: Option<JoinHandle<Vec<u8>>>,
    stderr_text: String,
    status: Option<ExitStatus>,
}

impl Proc {
    fn join_stderr(&mut self) {
        if let Some(h) = self.stderr.take() {
            let bytes = h.join().unwrap_or_default();
            self.stderr_text = String::from_utf8_lossy(&bytes).into_owned();
        }
    }

    /// Wait for a child whose stdout has ended; an error if it failed.
    fn finish(&mut self) -> Result<(), Error> {
        if self.status.is_none() {
            let status = self.child.wait().map_err(Error::from)?;
            self.status = Some(status);
            self.join_stderr();
        }
        match self.status {
            Some(s) if !s.success() => Err(tool_error(&self.program, s, &self.stderr_text)),
            _ => Ok(()),
        }
    }

    /// Kill the process group and reap the child. Idempotent.
    fn kill(&mut self) {
        if self.status.is_some() {
            return;
        }
        // The child has not been reaped, so its pid (= process group id)
        // cannot have been reused.
        #[cfg(unix)]
        if let Ok(pgid) = libc::pid_t::try_from(self.child.id()) {
            // SAFETY: plain syscall; a negative pid names the process group.
            unsafe {
                libc::kill(-pgid, libc::SIGKILL);
            }
        }
        let _ = self.child.kill();
        self.status = self.child.wait().ok();
        self.join_stderr();
    }

    fn pid(&self) -> u32 {
        self.child.id()
    }
}

/// The FASTQ stream from the child's stdout.
struct Stream {
    source: FastqSource<Single>,
    batch: Option<Batch<Single>>,
    pos: usize,
    eof: bool,
}

/// One regrouped spot.
#[derive(Default)]
struct SpotBuf {
    id: u64,
    bases: Vec<u8>,
    quals: Vec<u8>,
    starts: Vec<i32>,
    lens: Vec<u32>,
    types: Vec<ReadType>,
    name: Vec<u8>,
}

impl SpotBuf {
    fn reset(&mut self, id: u64) {
        self.id = id;
        self.bases.clear();
        self.quals.clear();
        self.starts.clear();
        self.lens.clear();
        self.types.clear();
        self.name.clear();
    }

    fn reads(&self) -> usize {
        self.lens.len()
    }
}

/// A reader over one `fasterq-dump` run. Sequential only
/// ([`random_access`](SpotReader::random_access) is false).
///
/// `fasterq-dump` is started on the first [`read_spot`](SpotReader::read_spot),
/// so opening (and [`run_info`](crate::run_info)) costs only the
/// `vdb-dump --id_range` call. Dropping the reader kills the child's
/// process group, reaps it and removes the temp dir.
pub struct FasterqDumpReader {
    backend: FasterqDumpBackend,
    path: PathBuf,
    input: String,
    names: bool,
    info: RunInfo,
    // Drop order: `proc` is killed in `Drop::drop`, then `stream` (the
    // stdout pipe) and `temp` go.
    stream: Option<Stream>,
    proc: Option<Proc>,
    temp: Option<TempDir>,
    cur: SpotBuf,
    cur_id: Option<u64>,
    ahead: SpotBuf,
    ahead_full: bool,
    last_id: Option<u64>,
    /// Read count of the first printed spot; empty spots get this many.
    template: Option<usize>,
    failed: Option<Error>,
}

impl FasterqDumpReader {
    /// The child's process id, once it has been started.
    pub fn child_pid(&self) -> Option<u32> {
        self.proc.as_ref().map(Proc::pid)
    }

    fn corrupt(&self, msg: impl Into<String>) -> Error {
        Error::new(ErrorKind::Corrupt, msg).with_input(self.input.clone())
    }

    fn start(&mut self) -> Result<(), Error> {
        let root = self
            .backend
            .temp_dir
            .clone()
            .unwrap_or_else(std::env::temp_dir);
        let temp = TempDir::create(&root)?;
        let mut cmd = self.backend.command(&self.path, &temp.0, self.names);
        self.temp = Some(temp);
        let program = self.backend.program.clone();
        let mut child = cmd.spawn().map_err(|e| spawn_error(&program, e))?;
        let stdout = child.stdout.take();
        let stderr = child.stderr.take();
        let stderr = stderr.and_then(|mut err| {
            std::thread::Builder::new()
                .name("readkit-fasterq-stderr".into())
                .spawn(move || {
                    let mut tail = Vec::new();
                    let mut buf = [0u8; 4096];
                    while let Ok(n) = err.read(&mut buf) {
                        if n == 0 {
                            break;
                        }
                        tail.extend_from_slice(&buf[..n]);
                        if tail.len() > 2 * STDERR_KEEP {
                            tail.drain(..tail.len() - STDERR_KEEP);
                        }
                    }
                    tail
                })
                .ok()
        });
        self.proc = Some(Proc {
            child,
            program: program.clone(),
            stderr,
            stderr_text: String::new(),
            status: None,
        });
        let Some(stdout) = stdout else {
            return Err(Error::new(ErrorKind::Internal, "child stdout not piped"));
        };
        let input = FastqInput::named(
            BufReader::with_capacity(1 << 20, stdout),
            format!("{} stdout", program.display()),
        );
        self.stream = Some(Stream {
            source: FastqSource::single(input, FastqOptions::default()),
            batch: None,
            pos: 0,
            eof: false,
        });
        Ok(())
    }

    /// The child's exit status at the end of its output: an error if it
    /// failed.
    fn finish_child(&mut self) -> Result<(), Error> {
        match &mut self.proc {
            Some(p) => p.finish().map_err(|e| e.with_input(self.input.clone())),
            None => Ok(()),
        }
    }

    /// A FASTQ parse error: prefer the child's own failure if it has one.
    fn parse_failure(&mut self, e: readkit_fastq::Error) -> Error {
        let truncated = matches!(
            e.kind(),
            readkit_fastq::ErrorKind::Truncated | readkit_fastq::ErrorKind::Io
        );
        if truncated {
            // The output ended early; the child is exiting or has exited.
            if let Err(child) = self.finish_child() {
                return child;
            }
        }
        let kind = if e.kind() == readkit_fastq::ErrorKind::Io {
            ErrorKind::Io
        } else {
            ErrorKind::Corrupt
        };
        Error::new(kind, format!("parsing fasterq-dump output: {e}"))
            .with_source(e)
            .with_input(self.input.clone())
    }

    /// The spot id of the next unconsumed record, or `None` at the end of a
    /// successful run.
    fn peek(&mut self) -> Result<Option<u64>, Error> {
        loop {
            let stream = self.stream.as_mut().expect("started");
            if let Some(b) = &stream.batch
                && stream.pos < b.len()
            {
                let col = b.column();
                let head = col.names().map(|n| n.get(stream.pos)).unwrap_or(b"");
                let id = head.split(|&c| c == b' ').next().unwrap_or(b"");
                if let Some(id) = std::str::from_utf8(id).ok().and_then(|s| s.parse().ok()) {
                    return Ok(Some(id));
                }
                let msg = format!(
                    "fasterq-dump record header {:?} does not start with a spot id",
                    String::from_utf8_lossy(head)
                );
                return Err(self.corrupt(msg));
            }
            if stream.eof {
                return Ok(None);
            }
            stream.batch = None;
            match stream.source.next_batch() {
                Ok(Some(b)) => {
                    stream.batch = Some(b);
                    stream.pos = 0;
                }
                Ok(None) => {
                    stream.eof = true;
                    self.finish_child()?;
                    return Ok(None);
                }
                Err(e) => return Err(self.parse_failure(e)),
            }
        }
    }

    /// Group the next spot's records into `ahead`. False at the end.
    fn fill_ahead(&mut self) -> Result<bool, Error> {
        let Some(id) = self.peek()? else {
            return Ok(false);
        };
        if let Some(last) = self.last_id
            && id <= last
        {
            return Err(self.corrupt(format!("fasterq-dump printed spot {id} after spot {last}")));
        }
        if id < self.info.first_spot || id >= self.info.end_spot() {
            return Err(self.corrupt(format!(
                "fasterq-dump printed spot {id}, outside the run's spots {}..{}",
                self.info.first_spot,
                self.info.end_spot()
            )));
        }
        self.ahead.reset(id);
        self.last_id = Some(id);
        loop {
            {
                let stream = self.stream.as_mut().expect("started");
                let b = stream.batch.as_ref().expect("peeked");
                let col = b.column();
                let i = stream.pos;
                let (seq, qual) = (col.seqs().get(i), col.quals().get(i));
                let a = &mut self.ahead;
                let start = i32::try_from(a.bases.len()).map_err(|_| {
                    Error::new(ErrorKind::Corrupt, "spot longer than 2^31 bases").at_spot(id)
                })?;
                if a.reads() == 0 && self.names {
                    let head = col.names().map(|n| n.get(i)).unwrap_or(b"");
                    if let Some(sp) = head.iter().position(|&c| c == b' ') {
                        a.name.extend_from_slice(&head[sp + 1..]);
                    }
                }
                a.bases.extend_from_slice(seq);
                a.quals.extend_from_slice(qual);
                a.starts.push(start);
                a.lens.push(seq.len() as u32);
                a.types.push(ReadType::BIOLOGICAL_READ);
                stream.pos += 1;
            }
            match self.peek()? {
                Some(next) if next == id => {}
                // Checked now rather than at the next spot, so that records
                // after the run's last spot are not silently dropped.
                Some(next) if next < id => {
                    return Err(
                        self.corrupt(format!("fasterq-dump printed spot {next} after spot {id}"))
                    );
                }
                Some(next) if next >= self.info.end_spot() => {
                    return Err(self.corrupt(format!(
                        "fasterq-dump printed spot {next}, outside the run's spots {}..{}",
                        self.info.first_spot,
                        self.info.end_spot()
                    )));
                }
                _ => break,
            }
        }
        self.ahead_full = true;
        Ok(true)
    }

    fn advance(&mut self, spot: u64) -> Result<(), Error> {
        if self.stream.is_none() {
            self.start()?;
        }
        // Skip printed spots before `spot` (a spot range not at the start).
        loop {
            if !self.ahead_full && !self.fill_ahead()? {
                break;
            }
            if self.ahead.id >= spot {
                break;
            }
            self.ahead_full = false;
        }
        if self.ahead_full && self.ahead.id == spot {
            std::mem::swap(&mut self.cur, &mut self.ahead);
            self.ahead_full = false;
            self.template.get_or_insert(self.cur.reads());
        } else {
            // Not printed: every read empty (or filtered out).
            let n = match self.template {
                Some(n) => n,
                None if self.ahead_full => self.ahead.reads(),
                None => 0,
            };
            self.cur.reset(spot);
            self.cur.starts.resize(n, 0);
            self.cur.lens.resize(n, 0);
            self.cur.types.resize(n, ReadType::BIOLOGICAL_READ);
            if self.names {
                self.cur.name.extend_from_slice(spot.to_string().as_bytes());
            }
        }
        self.cur_id = Some(spot);
        Ok(())
    }

    fn fail(&mut self, e: Error) -> Error {
        if let Some(p) = &mut self.proc {
            p.kill();
        }
        self.failed = Some(e.clone());
        e
    }
}

impl SpotReader for FasterqDumpReader {
    fn info(&self) -> &RunInfo {
        &self.info
    }

    fn random_access(&self) -> bool {
        false
    }

    fn read_spot(&mut self, spot: u64) -> Result<Spot<'_>, Error> {
        if let Some(e) = &self.failed {
            return Err(e.clone());
        }
        if spot < self.info.first_spot || spot >= self.info.end_spot() {
            return Err(Error::new(ErrorKind::NotFound, "spot id out of range")
                .with_input(self.input.clone())
                .at_spot(spot));
        }
        match self.cur_id {
            Some(c) if c == spot => {}
            Some(c) if spot < c => {
                return Err(Error::new(
                    ErrorKind::Unsupported,
                    format!("fasterq-dump is sequential: spot {spot} asked for after spot {c}"),
                )
                .with_input(self.input.clone())
                .at_spot(spot));
            }
            _ => {
                if let Err(e) = self.advance(spot) {
                    let e = if e.spot().is_none() {
                        e.at_spot(spot)
                    } else {
                        e
                    };
                    return Err(self.fail(e));
                }
            }
        }
        let c = &self.cur;
        Ok(Spot {
            bases: &c.bases,
            quals: &c.quals,
            starts: &c.starts,
            lens: &c.lens,
            types: &c.types,
            name: self.names.then_some(c.name.as_slice()),
        })
    }
}

impl Drop for FasterqDumpReader {
    fn drop(&mut self) {
        if let Some(p) = &mut self.proc {
            p.kill();
        }
    }
}

impl std::fmt::Debug for FasterqDumpReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FasterqDumpReader")
            .field("path", &self.path)
            .field("info", &self.info)
            .field("child_pid", &self.child_pid())
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn id_range_parses() {
        assert_eq!(
            parse_id_range("id-range: first-row = 1, row-count = 89,234\n"),
            Some((1, 89_234))
        );
        assert_eq!(
            parse_id_range("id-range: first-row = 7, row-count = 3"),
            Some((7, 3))
        );
        assert_eq!(parse_id_range("nothing"), None);
    }

    #[test]
    fn stderr_classification() {
        let nf = "2026-11-15T01:34:45 fasterq-dump.3.2.1 err: cannot find 'x' -> \
                  RC(rcVFS,rcQuery,rcResolving,rcName,rcNotFound)";
        assert_eq!(classify_stderr(nf), ErrorKind::NotFound);
        let bad = "err: file invalid while opening directory\n\
                   err: insp_seq_data().VCursorOpen( 'SEQUENCE' ) -> \
                   RC(rcFS,rcDirectory,rcOpening,rcFile,rcInvalid)\n";
        assert_eq!(classify_stderr(bad), ErrorKind::Corrupt);
        // The last RC wins.
        let two = format!("{nf}\n{bad}");
        assert_eq!(classify_stderr(&two), ErrorKind::Corrupt);
        assert_eq!(
            classify_stderr("RC(rcVDB,rcMgr,rcOpening,rcSchema,rcBadVersion)"),
            ErrorKind::Unsupported
        );
        assert_eq!(
            classify_stderr("RC(rcFS,rcFile,rcWriting,rcStorage,rcExhausted)"),
            ErrorKind::Io
        );
        assert_eq!(
            classify_stderr(
                "vdb-dump int: item not found while constructing - the path 'x' cannot be opened as vdb-database or vdb-table"
            ),
            ErrorKind::Unsupported
        );
        assert_eq!(
            classify_stderr("fasterq-dump err: the input data is missing the QUALITY-column"),
            ErrorKind::Unsupported
        );
        assert_eq!(classify_stderr("segmentation fault"), ErrorKind::Backend);
        assert_eq!(classify_stderr("RC(broken"), ErrorKind::Backend);
    }

    #[test]
    fn stderr_tail_keeps_last_lines() {
        let text = (0..20).map(|i| format!("line {i}\n")).collect::<String>() + "\n=====\n";
        let tail = stderr_tail(&text);
        assert!(tail.starts_with("line 14 | "), "{tail}");
        assert!(tail.ends_with("line 19"), "{tail}");
    }
}
