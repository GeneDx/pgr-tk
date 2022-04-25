libwfa
========

Rust bindings for the [wavefront
algorithm](https://github.com/smarco/WFA) for pairwise sequence
alignment.

## Usage

This crate will handle compiling the C library, and statically link it.

Just add `libwfa` to your Cargo dependencies:

```toml
[dependencies]
libwfa = "0.1"
```

## Dependencies

As a binding, llvm and libclang are required on Unix systems. These can be installed by package managers, for example on Ubuntu with:

```bash
sudo apt install llvm
sudo apt install libclang-dev
```

## Example

At this stage, usage maps closely to the C library.

This is equivalent to the basic example from the [WFA readme](https://github.com/smarco/WFA):

```rust
use libwfa::{affine_wavefront::*, bindings::*, mm_allocator::*, penalties::*};

fn main() {
    let alloc = MMAllocator::new(BUFFER_SIZE_8M as u64);

    let pattern = String::from("TCTTTACTCGCGCGTTGGAGAAATACAATAGT");
    let text = String::from("TCTATACTGCGCGTTTGGAGAAATAAAATAGT");

    let mut penalties = AffinePenalties {
        match_: 0,
        mismatch: 4,
        gap_opening: 6,
        gap_extension: 2,
    };

    let pat_len = pattern.as_bytes().len();
    let text_len = text.as_bytes().len();

    let mut wavefronts = AffineWavefronts::new_complete(
        pat_len,
        text_len,
        &mut penalties,
        &alloc,
    );

    wavefronts
        .align(pattern.as_bytes(), text.as_bytes())
        .unwrap();

    let score = wavefronts.edit_cigar_score(&mut penalties);

    println!("score: {}", score);
    wavefronts.print_cigar(pattern.as_bytes(), text.as_bytes());

    // The cigar can also be extracted as a byte vector
    let cigar = wavefronts.cigar_bytes_raw();
    let cg_str = std::str::from_utf8(&cigar).unwrap();
    println!("cigar: {}", cg_str);

    // Or as a prettier byte vector

    let cigar = wavefronts.cigar_bytes();
    let cg_str = std::str::from_utf8(&cigar).unwrap();
    println!("cigar: {}", cg_str);

}
```

See the tests for more examples.

## Build from source

Make sure to clone with the WFA submodule:

```bash
git clone --recursive https://github.com/chfi/wfa-rs
cd wfa-rs
cargo build
```
