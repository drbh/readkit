//! Builds the VDB shim only with the `ngs` feature.
//!
//! Links the SRA Toolkit's `libncbi-ngs` (which exports the VDB C API) from
//! `$NGS_PREFIX/lib`, default `/opt/homebrew`. The shim declares the few VDB
//! prototypes it uses itself, because the toolkit installs no VDB headers.

fn main() {
    println!("cargo:rerun-if-changed=shim/vdb_shim.c");
    println!("cargo:rerun-if-env-changed=NGS_PREFIX");
    #[cfg(feature = "ngs")]
    build_shim();
}

#[cfg(feature = "ngs")]
fn build_shim() {
    let prefix = std::env::var("NGS_PREFIX").unwrap_or_else(|_| "/opt/homebrew".to_string());
    cc::Build::new()
        .file("shim/vdb_shim.c")
        .flag_if_supported("-std=c99")
        .warnings(true)
        .compile("readkit_vdb_shim");
    println!("cargo:rustc-link-search=native={prefix}/lib");
    println!("cargo:rustc-link-lib=dylib=ncbi-ngs");
    // The library is a dylib outside the default search path; record it so
    // test and example binaries find it without DYLD_LIBRARY_PATH.
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("macos")
        || std::env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("linux")
    {
        println!("cargo:rustc-link-arg=-Wl,-rpath,{prefix}/lib");
    }
}
