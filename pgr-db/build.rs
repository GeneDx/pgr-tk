extern crate bindgen;

use std::{
    env,
    fs::{read_dir, remove_dir_all},
    path::PathBuf,
    process::Command,
};

fn build_agc() -> Option<()> {
    let mut wfa_dir = read_dir(&"../agc").ok()?;
    if !wfa_dir.any(|f| f.unwrap().file_name() == "makefile") {
        return None;
    }

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());

    let agc_path = out_path.join("agc");

    let _ = remove_dir_all(agc_path.as_path());

    // copy the WFA dir to OUT_PATH and build it there... clunky, but
    // don't want to pull in the entire 100MB WFA repo, since git2
    // doesn't seem to support shallow clones, and build scripts
    // should only modify things inside OUT_PATH. since the WFA folder
    // is just a couple MB, this is fine for now.
    let _cp_wfa = Command::new("cp")
        .arg("-r")
        .arg("../agc")
        .arg(&out_path)
        .output()
        .unwrap();

    let output = Command::new("make")
        .arg("-f")
        .arg("makefile.release")
        .arg("clean")
        .arg("libagc")
        .current_dir(&agc_path)
        .output()
        .unwrap();
    if output.status.success() {
        Some(())
    } else {
        panic!("make error: {}", String::from_utf8_lossy(&output.stderr));
    }
}
fn main() {
    if build_agc().is_none() {
        panic!("Error building AGC C library");
    }
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    let agc_path = out_path.join("agc");

    // shared library.
    println!("cargo:rustc-link-lib=agc");
    println!("cargo:rustc-link-search={}", agc_path.display());
    println!("cargo:rustc-link-lib=zstd");
    println!("cargo:rustc-link-search={}/libs", agc_path.display());
    println!("cargo:rustc-link-lib=stdc++");
    println!("cargo:rustc-link-search=/usr/lib/gcc/x86_64-linux-gnu/9/");

    // Tell cargo to invalidate the built crate whenever the wrapper changes
    println!("cargo:rerun-if-changed=wrapper.h");

    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header("wrapper.h")
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
