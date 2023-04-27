extern crate bindgen;
use std::env::consts::{ARCH, OS};

#[cfg(debug_assertions)]
const BUILD_TYPE: &str = "debug";
#[cfg(not(debug_assertions))]
const BUILD_TYPE: &str = "release";

#[cfg(feature = "with_agc")]
use std::fs::{read_dir, remove_dir_all};

#[cfg(feature = "with_agc")]
use std::path::PathBuf;

use std::{env, process::Command};

#[cfg(feature = "with_agc")]
fn build_agc() -> Option<()> {
    let mut agc_dir = read_dir("../agc").ok()?;
    if !agc_dir.any(|f| f.unwrap().file_name() == "makefile") {
        return None;
    }

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());

    let agc_path = out_path.join("agc");

    let _ = remove_dir_all(agc_path.as_path());

    // copy the AGC dir to OUT_PATH and build it there... clunky, but
    // don't want to pull in the entire 100MB WFA repo, since git2
    // doesn't seem to support shallow clones, and build scripts
    // should only modify things inside OUT_PATH. since the WFA folder
    // is just a couple MB, this is fine for now.
    let _cp_agc = Command::new("cp")
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

// fn wfa() {
//     let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
//     let _cp_agc = Command::new("cp")
//         .arg("../WFA2-lib/lib/libwfa.a")
//         .arg(&out_path)
//         .output()
//         .unwrap();
//     // The directory of the WFA libraries, added to the search path.
//     println!("cargo:rustc-link-search={}", out_path.display());
//     // Link the `wfa-lib` library.
//     println!("cargo:rustc-link-lib=wfa");
//     // Also link `omp`.
//     println!("cargo:rustc-link-lib=omp5");
// }

fn main() {
    //wfa();

    #[cfg(feature = "with_agc")]
    if build_agc().is_none() {
        panic!("Error building AGC C library");
    } else {
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
    // from https://vallentin.dev/2019/06/06/versioning
    let branch_name = get_branch_name();
    if branch_name != *"bioconda" {
        let version_string = format!(
            "{} {} ({}:{}{}, {} build, {} [{}] [{}])",
            env!("CARGO_PKG_NAME"),
            env!("CARGO_PKG_VERSION"),
            get_branch_name(),
            get_commit_hash(),
            if is_working_tree_clean() { "" } else { "+" },
            BUILD_TYPE,
            OS,
            ARCH,
            get_rustc_version()
        );

        println!("cargo:rustc-env=VERSION_STRING={}", version_string);
    } else {
        let version_string = format!(
            "{} {} (bioconda {} build ({}:{}{}), {} [{}] [{}])",
            env!("CARGO_PKG_NAME"),
            env!("CARGO_PKG_VERSION"),
            BUILD_TYPE,
            get_branch_name(),
            get_commit_hash(),
            if is_working_tree_clean() { "" } else { "+" },
            OS,
            ARCH,
            get_rustc_version()
        );
        println!("cargo:rustc-env=VERSION_STRING={}", version_string);
    }
}

fn get_rustc_version() -> String {
    let output = Command::new("rustc")
        .arg("--version")
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .unwrap();

    assert!(output.status.success());

    String::from_utf8_lossy(&output.stdout)
        .trim_end()
        .to_string()
}

fn get_commit_hash() -> String {
    let output = Command::new("git")
        .arg("log")
        .arg("-1")
        .arg("--pretty=format:%h") // Abbreviated commit hash
        // .arg("--pretty=format:%H") // Full commit hash
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .unwrap();

    // assert!(output.status.success());
    if output.status.success() {
        String::from_utf8_lossy(&output.stdout).to_string()
    } else {
        String::from("bioconda")
    }
}

fn get_branch_name() -> String {
    let output = Command::new("git")
        .arg("rev-parse")
        .arg("--abbrev-ref")
        .arg("HEAD")
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .unwrap();

    //assert!(output.status.success());
    if output.status.success() {
        String::from_utf8_lossy(&output.stdout)
            .trim_end()
            .to_string()
    } else {
        String::from("bioconda")
    }
}

fn is_working_tree_clean() -> bool {
    let status = Command::new("git")
        .arg("diff")
        .arg("--quiet")
        .arg("--exit-code")
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .status()
        .unwrap();

    if status.success() {
        status.code().unwrap() == 0
    } else {
        true
    }
}
