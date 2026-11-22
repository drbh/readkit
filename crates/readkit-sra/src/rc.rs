//! VDB status codes (`rc_t`) mapped to [`ErrorKind`]s.

#![cfg_attr(not(feature = "ngs"), allow(dead_code))]

use crate::error::ErrorKind;

// VDB `RCState` values (ncbi-vdb `interfaces/kfc/rc.h`).
const RC_UNSUPPORTED: u32 = 3;
const RC_BAD_VERSION: u32 = 8;
const RC_CORRUPT: u32 = 11;
const RC_INCORRECT: u32 = 12;
const RC_INCONSISTENT: u32 = 13;
const RC_INCOMPLETE: u32 = 15;
const RC_NOT_FOUND: u32 = 24;
const RC_UNAUTHORIZED: u32 = 29;
const RC_NO_PERM: u32 = 32;
const RC_TOO_SHORT: u32 = 34;
const RC_WRONG_TYPE: u32 = 48;

/// Map a VDB `rc_t` to an error kind by its state field (low 6 bits).
pub(crate) fn rc_kind(rc: u32) -> ErrorKind {
    match rc & 0x3F {
        RC_NOT_FOUND => ErrorKind::NotFound,
        RC_CORRUPT | RC_INCORRECT | RC_INCONSISTENT | RC_INCOMPLETE | RC_TOO_SHORT => {
            ErrorKind::Corrupt
        }
        RC_UNSUPPORTED | RC_BAD_VERSION | RC_WRONG_TYPE => ErrorKind::Unsupported,
        RC_UNAUTHORIZED | RC_NO_PERM => ErrorKind::Io,
        _ => ErrorKind::Backend,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rc_states_map_to_kinds() {
        // RC(rcVFS,rcMgr,rcOpening,rcDirectory,rcNotFound), seen for a
        // missing file.
        assert_eq!(rc_kind(2615477144), ErrorKind::NotFound);
        assert_eq!(rc_kind(RC_CORRUPT), ErrorKind::Corrupt);
        assert_eq!(rc_kind((5 << 27) | RC_UNSUPPORTED), ErrorKind::Unsupported);
        assert_eq!(rc_kind(RC_NO_PERM), ErrorKind::Io);
        assert_eq!(rc_kind(46), ErrorKind::Backend);
    }
}
