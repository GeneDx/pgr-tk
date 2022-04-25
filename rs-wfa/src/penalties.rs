use libc::c_int;

use crate::bindings::affine_penalties_t;

#[repr(C)]
pub struct AffinePenalties {
    pub match_: c_int,
    pub mismatch: c_int,
    pub gap_opening: c_int,
    pub gap_extension: c_int,
}

impl AffinePenalties {
    pub fn as_ptr(&mut self) -> *mut affine_penalties_t {
        let ptr = self as *mut AffinePenalties;
        ptr.cast()
    }
}
