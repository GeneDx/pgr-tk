use crate::bindings::*;

#[repr(C)]
pub struct MMAllocator {
    allocator: *mut mm_allocator_t,
}

// TODO need to track the items that have actually been allocated, too
impl Drop for MMAllocator {
    fn drop(&mut self) {
        unsafe {
            mm_allocator_delete(self.allocator);
        }
    }
}

impl MMAllocator {
    pub fn new(buffer_size: u64) -> Self {
        let allocator = unsafe { mm_allocator_new(buffer_size) };
        Self { allocator }
    }

    pub fn alloc_ptr(&self) -> *mut mm_allocator_t {
        self.allocator
    }
}
