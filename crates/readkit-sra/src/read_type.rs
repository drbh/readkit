//! Per-read type bits (INSDC `xread_type`).

use std::fmt;

/// A read's type within its spot, as stored in VDB's `READ_TYPE` column
/// (`INSDC:SRA:xread_type`).
///
/// Bit 0 is biological (set) vs technical (clear). Bits 1 and 2 are the
/// orientation flags *forward* and *reverse*, which may be set on technical
/// reads too: `0b010` is a forward **technical** read. Only bit 0 decides
/// biological vs technical. (ncbi-vdb-sys treats any non-zero value as
/// biological and so labels such reads biological.)
///
/// ```
/// use readkit_sra::{Orientation, ReadType};
/// let barcode = ReadType::from_bits(0b010);
/// assert!(barcode.is_technical());
/// assert_eq!(barcode.orientation(), Orientation::Forward);
/// assert!(ReadType::from_bits(0b101).is_biological());
/// ```
#[derive(Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(transparent)]
pub struct ReadType(u8);

/// Orientation flags of a [`ReadType`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Orientation {
    /// Neither flag set.
    Unspecified,
    /// `SRA_READ_TYPE_FORWARD` (bit 1).
    Forward,
    /// `SRA_READ_TYPE_REVERSE` (bit 2).
    Reverse,
    /// Both flags set (invalid per INSDC, reported as stored).
    Both,
}

impl ReadType {
    /// `SRA_READ_TYPE_BIOLOGICAL`.
    pub const BIOLOGICAL: u8 = 1;
    /// `SRA_READ_TYPE_FORWARD`.
    pub const FORWARD: u8 = 2;
    /// `SRA_READ_TYPE_REVERSE`.
    pub const REVERSE: u8 = 4;

    /// A technical read with no orientation.
    pub const TECHNICAL_READ: ReadType = ReadType(0);
    /// A biological read with no orientation.
    pub const BIOLOGICAL_READ: ReadType = ReadType(1);

    /// Wrap raw `xread_type` bits.
    pub const fn from_bits(bits: u8) -> Self {
        Self(bits)
    }

    /// The raw bits.
    pub const fn bits(self) -> u8 {
        self.0
    }

    /// Whether bit 0 (biological) is set.
    pub const fn is_biological(self) -> bool {
        self.0 & Self::BIOLOGICAL != 0
    }

    /// Whether bit 0 is clear (barcodes, UMIs, adapters, linkers).
    pub const fn is_technical(self) -> bool {
        !self.is_biological()
    }

    /// The orientation flags.
    pub const fn orientation(self) -> Orientation {
        match (self.0 & Self::FORWARD != 0, self.0 & Self::REVERSE != 0) {
            (false, false) => Orientation::Unspecified,
            (true, false) => Orientation::Forward,
            (false, true) => Orientation::Reverse,
            (true, true) => Orientation::Both,
        }
    }
}

impl fmt::Debug for ReadType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let kind = if self.is_biological() {
            "biological"
        } else {
            "technical"
        };
        match self.orientation() {
            Orientation::Unspecified => write!(f, "ReadType({kind})"),
            o => write!(f, "ReadType({kind}, {o:?})"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_bit_pattern() {
        // (bits, biological, orientation)
        let table = [
            (0, false, Orientation::Unspecified),
            (1, true, Orientation::Unspecified),
            (2, false, Orientation::Forward),
            (3, true, Orientation::Forward),
            (4, false, Orientation::Reverse),
            (5, true, Orientation::Reverse),
            (6, false, Orientation::Both),
            (7, true, Orientation::Both),
        ];
        for (bits, bio, orient) in table {
            let t = ReadType::from_bits(bits);
            assert_eq!(t.is_biological(), bio, "bits {bits}");
            assert_eq!(t.is_technical(), !bio, "bits {bits}");
            assert_eq!(t.orientation(), orient, "bits {bits}");
            assert_eq!(t.bits(), bits);
        }
    }

    #[test]
    fn oriented_technical_reads_are_not_biological() {
        // The ncbi-vdb-sys bug: `bits != 0` would call these biological.
        for bits in [ReadType::FORWARD, ReadType::REVERSE] {
            assert!(ReadType::from_bits(bits).is_technical());
        }
    }

    #[test]
    fn debug_names() {
        assert_eq!(
            format!("{:?}", ReadType::from_bits(2)),
            "ReadType(technical, Forward)"
        );
        assert_eq!(
            format!("{:?}", ReadType::BIOLOGICAL_READ),
            "ReadType(biological)"
        );
    }
}
