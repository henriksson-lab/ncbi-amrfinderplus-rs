fn main() {
    // Phase 1: Compile C++ sources for FFI
    // This will be populated as we add FFI bindings
    // For now, just ensure the build script exists

    // Tell cargo to re-run if C++ sources change
    println!("cargo:rerun-if-changed=../amr/common.cpp");
    println!("cargo:rerun-if-changed=../amr/common.hpp");
    println!("cargo:rerun-if-changed=../amr/seq.cpp");
    println!("cargo:rerun-if-changed=../amr/seq.hpp");
    println!("cargo:rerun-if-changed=../amr/alignment.cpp");
    println!("cargo:rerun-if-changed=../amr/alignment.hpp");
    println!("cargo:rerun-if-changed=../amr/graph.cpp");
    println!("cargo:rerun-if-changed=../amr/graph.hpp");
    println!("cargo:rerun-if-changed=../amr/gff.cpp");
    println!("cargo:rerun-if-changed=../amr/gff.hpp");
    println!("cargo:rerun-if-changed=../amr/tsv.cpp");
    println!("cargo:rerun-if-changed=../amr/tsv.hpp");
    println!("cargo:rerun-if-changed=../amr/amr_report.cpp");
    println!("cargo:rerun-if-changed=../amr/amrfinder.cpp");
    println!("cargo:rerun-if-changed=../amr/ffi_shim.cpp");
}
