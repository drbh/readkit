//! SRA Lite detection.
//!
//! SRA Lite runs keep bases but replace qualities with synthetic values.
//! NCBI documents them as one value for reads that pass the quality filter
//! (Q30, `?` in phred+33) and another for rejected reads (Q3, `$`). Decoded
//! qualities of a Lite run are therefore not real base qualities and must
//! not be fed to quality trimming as if they were.
//!
//! Two signals, strongest first:
//! 1. The run's own declaration: VDB reports whether a table can deliver
//!    original qualities (`VTableGetQualityCapability`). Backends map it to
//!    [`QualityKind`] in [`RunInfo`](crate::RunInfo).
//! 2. A histogram heuristic over the first [`LITE_SAMPLE`] reads
//!    ([`QualityProbe`]), used when the backend cannot tell
//!    ([`QualityKind::Unknown`]), for example a FASTQ stream.

/// What the qualities of a run are.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum QualityKind {
    /// Original (full) qualities as submitted.
    Original,
    /// Synthetic qualities only: an SRA Lite run.
    Synthetic,
    /// The backend cannot tell; the decoder samples qualities instead.
    Unknown,
}

/// What a decoder does with an SRA Lite run.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum LitePolicy {
    /// Fail with [`ErrorKind::Lite`](crate::ErrorKind::Lite): at open time
    /// when the run declares synthetic qualities, or after
    /// [`LITE_SAMPLE`] reads when sampling detects them. The default, so
    /// synthetic qualities are never mistaken for real ones.
    #[default]
    Reject,
    /// Decode anyway. [`Decoder::quality`](crate::Decoder::quality) still
    /// reports what the run is.
    Allow,
}

/// Reads sampled by [`QualityProbe`] before it decides.
pub const LITE_SAMPLE: u64 = 100_000;

/// Quality-histogram heuristic for SRA Lite.
///
/// Counts quality bytes of up to [`LITE_SAMPLE`] reads. Qualities look
/// synthetic when at most two values each make up more than 0.1 % of the
/// sampled bytes.
///
/// ```
/// use readkit_sra::QualityProbe;
/// let mut p = QualityProbe::new();
/// for _ in 0..10 { p.observe(b"????????"); }
/// assert!(p.looks_synthetic());
/// p.observe(b"ABCDEFGHIJ");
/// assert_eq!(p.distinct_values(), 11);
/// ```
#[derive(Clone)]
pub struct QualityProbe {
    hist: [u64; 256],
    sampled: u64,
}

impl QualityProbe {
    /// An empty probe.
    pub fn new() -> Self {
        Self {
            hist: [0; 256],
            sampled: 0,
        }
    }

    /// Count one read's quality bytes, until [`LITE_SAMPLE`] reads were seen.
    pub fn observe(&mut self, qual: &[u8]) {
        if self.sampled >= LITE_SAMPLE {
            return;
        }
        self.sampled += 1;
        for &b in qual {
            self.hist[b as usize] += 1;
        }
    }

    /// Reads counted so far.
    pub fn sampled(&self) -> u64 {
        self.sampled
    }

    /// Whether the sample is complete.
    pub fn is_full(&self) -> bool {
        self.sampled >= LITE_SAMPLE
    }

    /// Distinct quality bytes seen.
    pub fn distinct_values(&self) -> usize {
        self.hist.iter().filter(|c| **c > 0).count()
    }

    /// Whether the qualities seen so far look synthetic (see type docs).
    /// `false` when nothing was sampled.
    pub fn looks_synthetic(&self) -> bool {
        let total: u64 = self.hist.iter().sum();
        let significant = self.hist.iter().filter(|c| **c * 1000 > total).count();
        total > 0 && significant <= 2
    }
}

impl Default for QualityProbe {
    fn default() -> Self {
        Self::new()
    }
}

impl std::fmt::Debug for QualityProbe {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("QualityProbe")
            .field("sampled", &self.sampled)
            .field("distinct_values", &self.distinct_values())
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lite_pattern_is_synthetic() {
        // Q30 for passing reads, Q3 for rejected ones.
        let mut p = QualityProbe::new();
        for i in 0..1000 {
            p.observe(if i % 10 == 0 {
                b"$$$$$$$$"
            } else {
                b"????????"
            });
        }
        assert!(p.looks_synthetic());
        assert_eq!(p.distinct_values(), 2);
    }

    #[test]
    fn real_qualities_are_not() {
        let mut p = QualityProbe::new();
        for i in 0..1000u32 {
            let q: Vec<u8> = (0..50).map(|j| b'#' + ((i + j) % 40) as u8).collect();
            p.observe(&q);
        }
        assert!(!p.looks_synthetic());
    }

    #[test]
    fn rare_outliers_do_not_hide_lite() {
        // A handful of odd bytes (< 0.1 %) do not count as significant.
        let mut p = QualityProbe::new();
        for _ in 0..10_000 {
            p.observe(b"??????????");
        }
        p.observe(b"ABC");
        assert!(p.looks_synthetic());
        assert_eq!(p.distinct_values(), 4);
    }

    #[test]
    fn empty_and_cap() {
        let mut p = QualityProbe::new();
        assert!(!p.looks_synthetic());
        for _ in 0..LITE_SAMPLE + 10 {
            p.observe(b"?");
        }
        assert_eq!(p.sampled(), LITE_SAMPLE);
        assert!(p.is_full());
        p.observe(b"ABCDEFGHIJKLMNOP");
        assert_eq!(
            p.distinct_values(),
            1,
            "nothing counted after the sample is full"
        );
    }
}
