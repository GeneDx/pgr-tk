use std::os::raw::{c_char, c_int};

use crate::{
    bindings::*, mm_allocator::MMAllocator, penalties::AffinePenalties,
};

#[derive(Debug, Clone, Copy)]
pub enum WavefrontError {
    InputLengthError,
}

impl std::fmt::Display for WavefrontError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WavefrontError::InputLengthError => write!(
                f,
                "Pattern and text must be no longer than \
                 the lengths used to allocate the wavefront"
            ),
        }
    }
}

impl std::error::Error for WavefrontError {}

/// Safe wrapper over an `affine_wavefronts_t` instance allocated by
/// a libwfa `mm_allocator`.
pub struct AffineWavefronts<'a> {
    ptr: *mut affine_wavefronts_t,
    // This allocator ref is mainly kept to force the wavefronts to be
    // dropped before the allocator is freed
    allocator: &'a MMAllocator,
    pattern_len: usize,
    text_len: usize,
}

impl<'a> AffineWavefronts<'a> {
    /// Construct a new set of complete wavefronts
    pub fn new_complete(
        pattern_len: usize,
        text_len: usize,
        penalties: &mut AffinePenalties,
        alloc: &'a MMAllocator,
    ) -> Self {
        assert!(pattern_len > 0 && text_len > 0);
        let stats_ptr = std::ptr::null_mut() as *mut wavefronts_stats_t;
        let ptr = unsafe {
            affine_wavefronts_new_complete(
                pattern_len as c_int,
                text_len as c_int,
                penalties.as_ptr(),
                stats_ptr,
                alloc.alloc_ptr(),
            )
        };
        AffineWavefronts {
            ptr,
            allocator: alloc,
            pattern_len,
            text_len,
        }
    }

    /// Construct a new set of reduced wavefronts
    pub fn new_reduced(
        pattern_len: usize,
        text_len: usize,
        penalties: &mut AffinePenalties,
        min_wavefront_len: i32,
        min_dist_threshold: i32,
        alloc: &'a MMAllocator,
    ) -> Self {
        assert!(pattern_len > 0 && text_len > 0);
        let stats_ptr = std::ptr::null_mut() as *mut wavefronts_stats_t;
        let ptr = unsafe {
            affine_wavefronts_new_reduced(
                pattern_len as c_int,
                text_len as c_int,
                penalties.as_ptr(),
                min_wavefront_len as c_int,
                min_dist_threshold as c_int,
                stats_ptr,
                alloc.alloc_ptr(),
            )
        };

        Self {
            ptr,
            allocator: alloc,
            pattern_len,
            text_len,
        }
    }

    /// Clear the wavefronts
    pub fn clear(&mut self) {
        unsafe {
            affine_wavefronts_clear(self.ptr);
        }
    }

    /// Align the given pattern and text string. Callers need to make
    /// sure the byteslices have the correct length compared to the
    /// lengths used to construct thing wavefronts object.
    ///
    /// Does *not* check that `pattern` and `text` are nul-terminated
    /// CStrings, since the C function used takes the string lengths
    /// as arguments.
    pub fn align(
        &mut self,
        pattern: &[u8],
        text: &[u8],
    ) -> Result<(), WavefrontError> {
        if pattern.len() > self.pattern_len || text.len() > self.text_len {
            return Err(WavefrontError::InputLengthError);
        }
        unsafe {
            affine_wavefronts_align(
                self.ptr,
                pattern.as_ptr() as *const c_char,
                pattern.len() as c_int,
                text.as_ptr() as *const c_char,
                text.len() as c_int,
            );
        }
        Ok(())
    }

    fn edit_cigar(&self) -> &edit_cigar_t {
        unsafe {
            let wf_ref = self.ptr.as_ref().unwrap();
            &wf_ref.edit_cigar
        }
    }

    /// Returns the cigar string for the wavefront alignment as a
    /// vector of bytes. Note that each operation is repeated however
    /// many times it applies, i.e. instead of "3M1X" you get "MMMX".
    pub fn cigar_bytes_raw(&self) -> Vec<u8> {
        let slice = unsafe { self.cigar_slice() };
        slice.into()
    }

    pub fn cigar_bytes(&self) -> Vec<u8> {
        let slice = unsafe { self.cigar_slice() };
        if slice.is_empty() {
            Vec::new()
        } else {
            compress_cigar(slice).unwrap()
        }
    }

    /// Returns a slice to the cigar string for the wavefront
    /// alignment. Unsafe as the slice is pointing to the
    /// `edit_cigar_t` managed by libwfa.
    pub unsafe fn cigar_slice(&self) -> &[u8] {
        let cigar = self.edit_cigar();
        let ops_ptr = cigar.operations as *mut u8;
        let start = ops_ptr.offset(cigar.begin_offset as isize);
        let len = (cigar.end_offset - cigar.begin_offset) as usize;
        std::slice::from_raw_parts(start, len)
    }

    /// Returns the alignment score
    pub fn edit_cigar_score(
        &mut self,
        penalties: &mut AffinePenalties,
    ) -> isize {
        let penalties = penalties as *mut AffinePenalties;
        let penalties_ptr: *mut affine_penalties_t = penalties.cast();
        let score = unsafe {
            let wf_ref = self.ptr.as_mut().unwrap();
            let cigar = &mut wf_ref.edit_cigar as *mut edit_cigar_t;
            edit_cigar_score_gap_affine(cigar, penalties_ptr)
        };

        score as isize
    }

    /// Prints the alignment using the C library pretty printer. For
    /// now it only prints to stderr.
    pub fn print_cigar(&mut self, pattern: &[u8], text: &[u8]) {
        unsafe {
            let wf_ref = self.ptr.as_mut().unwrap();
            let cg_mut = &mut wf_ref.edit_cigar as *mut edit_cigar_t;
            edit_cigar_print_pretty(
                stderr,
                pattern.as_ptr() as *const c_char,
                pattern.len() as c_int,
                text.as_ptr() as *const c_char,
                text.len() as c_int,
                cg_mut,
                self.allocator.alloc_ptr(),
            );
        }
    }
}

/// "Compresses" a cigar string produced by WFA (like "MMMXXII") into
/// a regular cigar string with operation counts (e.g. "3M2X2I")
fn compress_cigar(cigar: &[u8]) -> Option<Vec<u8>> {
    let mut result = Vec::new();

    let mut iter = cigar.iter();
    let mut last_op = iter.next().copied()?;
    let mut last_count = 1;
    for &op in iter {
        if op == last_op {
            last_count += 1;
        } else {
            let op_char = char::from(last_op);
            let string = format!("{}{}", last_count, op_char);
            result.extend(string.as_bytes());
            last_op = op;
            last_count = 1;
        }
    }

    let op_char = char::from(last_op);
    let string = format!("{}{}", last_count, op_char);
    result.extend(string.as_bytes());

    Some(result)
}

impl Drop for AffineWavefronts<'_> {
    fn drop(&mut self) {
        unsafe { affine_wavefronts_delete(self.ptr) }
    }
}
