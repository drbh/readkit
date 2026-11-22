//! Mate classification of FASTQ file names.
//!
//! A mate token is `R1`, `R2`, `1` or `2`, bounded on both sides by `_`,
//! `.`, `-` or the start or end of the name without its FASTQ extension.
//! Tokens are whole delimiter-separated words, so the digits inside an
//! accession never count: RNAA's `contains("R1")` matched every `SRR1…`
//! accession and labelled `SRR1…_2.fastq.gz` as R1.

use crate::error::{Error, ErrorKind};

/// Which mate a FASTQ file holds.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum Mate {
    /// Single-end reads (no mate token).
    Single,
    /// First mate.
    R1,
    /// Second mate.
    R2,
    /// Reads without a mate in a paired run (ENA's `<run>.fastq.gz` next to
    /// `<run>_1` and `<run>_2`).
    Unpaired,
}

const EXTENSIONS: [&str; 6] = [
    ".fastq.gz",
    ".fq.gz",
    ".fastq.bz2",
    ".fq.bz2",
    ".fastq",
    ".fq",
];

fn stem(basename: &str) -> &str {
    EXTENSIONS
        .iter()
        .find_map(|ext| basename.strip_suffix(ext))
        .unwrap_or(basename)
}

/// Classify a FASTQ file's basename by the mate grammar.
///
/// No token gives [`Mate::Single`]. Conflicting tokens (`x_R1_R2`) and
/// glued tokens (`R1R2`) are [`ErrorKind::Metadata`] errors.
///
/// ```
/// use readkit_sra::{Mate, classify_mate};
/// assert_eq!(classify_mate("SRR13638686_2.fastq.gz")?, Mate::R2);
/// assert_eq!(classify_mate("lane1_R1_001.fastq.gz")?, Mate::R1);
/// assert_eq!(classify_mate("SRR1.fastq.gz")?, Mate::Single); // not R1
/// # Ok::<(), readkit_sra::Error>(())
/// ```
pub fn classify_mate(basename: &str) -> Result<Mate, Error> {
    let mut found = None;
    for token in stem(basename).split(['_', '.', '-']) {
        let mate = match token {
            "R1" | "1" => Mate::R1,
            "R2" | "2" => Mate::R2,
            _ if token.starts_with("R1R2") || token.starts_with("R2R1") => {
                return Err(Error::new(
                    ErrorKind::Metadata,
                    format!("ambiguous mate token `{token}` in {basename}"),
                ));
            }
            _ => continue,
        };
        match found {
            Some(previous) if previous != mate => {
                return Err(Error::new(
                    ErrorKind::Metadata,
                    format!("conflicting mate tokens in {basename}"),
                ));
            }
            _ => found = Some(mate),
        }
    }
    Ok(found.unwrap_or(Mate::Single))
}

/// Classify one of ENA's `fastq_ftp` files for `run`: `<run>_1` / `<run>_2`
/// are R1 / R2, `<run>.fastq.gz` is single-end in a single-end run and
/// [`Mate::Unpaired`] in a paired one. Other names use [`classify_mate`].
#[cfg(feature = "resolve")]
pub(crate) fn ena_mate(run: &str, basename: &str, paired: bool) -> Result<Mate, Error> {
    let stem = basename.strip_suffix(".fastq.gz").unwrap_or("");
    if let Some(rest) = stem.strip_prefix(run) {
        match rest {
            "_1" => return Ok(Mate::R1),
            "_2" => return Ok(Mate::R2),
            "" => return Ok(if paired { Mate::Unpaired } else { Mate::Single }),
            _ => {}
        }
    }
    classify_mate(basename)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grammar_vectors() {
        assert_eq!(classify_mate("SRR13638686_2.fastq.gz").unwrap(), Mate::R2);
        assert_eq!(classify_mate("SRR1_1.fastq.gz").unwrap(), Mate::R1);
        assert_eq!(classify_mate("lane1_R2_001.fastq.gz").unwrap(), Mate::R2);
        assert_eq!(classify_mate("sample.R1.fq.gz").unwrap(), Mate::R1);
        assert_eq!(classify_mate("sample-2.fq").unwrap(), Mate::R2);
        assert_eq!(classify_mate("ERR123456.fastq.gz").unwrap(), Mate::Single);
        assert_eq!(classify_mate("reads.fastq.bz2").unwrap(), Mate::Single);
        assert_eq!(classify_mate("x_R1_1.fastq.gz").unwrap(), Mate::R1);
    }

    #[test]
    fn accession_digits_are_never_mates() {
        // The RNAA bug: `contains("R1")` matched these.
        for run in [
            "SRR13638686",
            "ERR1000001",
            "DRR2222222",
            "SRR1",
            "SRR2",
            "SRR12",
        ] {
            assert_eq!(
                classify_mate(&format!("{run}.fastq.gz")).unwrap(),
                Mate::Single
            );
            assert_eq!(
                classify_mate(&format!("{run}_1.fastq.gz")).unwrap(),
                Mate::R1
            );
            assert_eq!(
                classify_mate(&format!("{run}_2.fastq.gz")).unwrap(),
                Mate::R2
            );
        }
        assert_eq!(classify_mate("R1sample.fastq.gz").unwrap(), Mate::Single);
    }

    #[test]
    fn conflicts_are_errors() {
        for name in ["x_R1_R2.fastq.gz", "R1R2_x.fastq.gz", "a_1_2.fq", "R2R1.fq"] {
            let e = classify_mate(name).unwrap_err();
            assert_eq!(e.kind(), ErrorKind::Metadata, "{name}");
        }
    }

    #[cfg(feature = "resolve")]
    #[test]
    fn ena_names() {
        assert_eq!(ena_mate("SRR1", "SRR1_1.fastq.gz", true).unwrap(), Mate::R1);
        assert_eq!(ena_mate("SRR1", "SRR1_2.fastq.gz", true).unwrap(), Mate::R2);
        assert_eq!(
            ena_mate("SRR1", "SRR1.fastq.gz", true).unwrap(),
            Mate::Unpaired
        );
        assert_eq!(
            ena_mate("SRR1", "SRR1.fastq.gz", false).unwrap(),
            Mate::Single
        );
        // A run whose accession is a prefix of another's.
        assert_eq!(
            ena_mate("SRR1", "SRR12_2.fastq.gz", true).unwrap(),
            Mate::R2
        );
        assert_eq!(
            ena_mate("SRR1", "SRR12.fastq.gz", true).unwrap(),
            Mate::Single
        );
    }
}
