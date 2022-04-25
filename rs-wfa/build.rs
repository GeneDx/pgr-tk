use std::{
    env,
    fs::{read_dir, remove_dir_all},
    path::PathBuf,
    process::Command,
};

fn build_wfa() -> Option<()> {
    let mut wfa_dir = read_dir(&"./WFA").ok()?;
    if !wfa_dir.any(|f| f.unwrap().file_name() == "Makefile") {
        return None;
    }

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());

    let wfa_path = out_path.join("WFA");

    let _ = remove_dir_all(wfa_path.as_path());

    // copy the WFA dir to OUT_PATH and build it there... clunky, but
    // don't want to pull in the entire 100MB WFA repo, since git2
    // doesn't seem to support shallow clones, and build scripts
    // should only modify things inside OUT_PATH. since the WFA folder
    // is just a couple MB, this is fine for now.
    let _cp_wfa = Command::new("cp")
        .arg("-r")
        .arg("./WFA")
        .arg(&out_path)
        .output()
        .unwrap();

    let output = Command::new("make")
        .arg("clean")
        .arg("all")
        .current_dir(&wfa_path)
        .output()
        .unwrap();
    if output.status.success() {
        Some(())
    } else {
        panic!("make error: {}", String::from_utf8_lossy(&output.stderr));
    }
}

fn main() {
    if build_wfa().is_none() {
        panic!("Error building WFA C library");
    }

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());

    let wfa_path = out_path.join("WFA").join("build");

    println!("cargo:rustc-link-lib=wfa");
    println!("cargo:rustc-link-search={}", wfa_path.display());

    let bindings = bindgen::Builder::default()
        .clang_arg("-IWFA")
        .header("wrapper.h")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        .allowlist_var("stderr")
        .allowlist_var("stdout")
        .allowlist_function("mm_.*")
        .allowlist_type("mm_.*")
        .allowlist_function("affine_.*")
        .allowlist_type("affine_.*")
        .allowlist_function("edit_.*")
        .allowlist_type("edit_.*")
        // .whitelist_function("backtrace_.*")
        // .whitelist_type("backtrace_.*")
        .allowlist_function("alignment_.*")
        .allowlist_type("alignment_.*")
        // .whitelist_function("wavefront_.*")
        // .whitelist_type("wavefront_.*")
        // .whitelist_function("swg_.*")
        // .whitelist_var("METRIC_FACTOR_*")
        // .whitelist_var("NUM_LINES_*")
        .allowlist_var("BUFFER_SIZE_.*")
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());

    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
