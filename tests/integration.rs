// Integration tests comparing Rust output against C++ golden files

use std::path::PathBuf;
use std::process::Command;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

fn cpp_binary_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amr")
}

fn rust_binary() -> PathBuf {
    // Built by cargo test
    PathBuf::from(env!("CARGO_BIN_EXE_amrfinder"))
}

fn cpp_binary(name: &str) -> PathBuf {
    cpp_binary_dir().join(name)
}

// --- fasta_check tests ---

#[test]
fn test_fasta_check_protein_golden() {
    let cpp_bin = cpp_binary("fasta_check");
    let input = test_data_dir().join("test_prot.fa");
    if !cpp_bin.exists() || !input.exists() {
        return;
    }

    let cpp_out = Command::new(&cpp_bin)
        .arg(&input)
        .arg("-aa")
        .output()
        .expect("C++ fasta_check failed");

    let rust_out = Command::new(rust_binary())
        .args(["check-fasta", input.to_str().unwrap(), "--aa"])
        .output()
        .expect("Rust fasta_check failed");

    assert_eq!(
        String::from_utf8_lossy(&cpp_out.stdout),
        String::from_utf8_lossy(&rust_out.stdout),
        "fasta_check protein output differs"
    );
}

#[test]
fn test_fasta_check_dna_golden() {
    let cpp_bin = cpp_binary("fasta_check");
    let input = test_data_dir().join("test_dna.fa");
    if !cpp_bin.exists() || !input.exists() {
        return;
    }

    let cpp_out = Command::new(&cpp_bin)
        .arg(&input)
        .arg("-ambig")
        .output()
        .expect("C++ fasta_check failed");

    let rust_out = Command::new(rust_binary())
        .args(["check-fasta", input.to_str().unwrap(), "--ambig"])
        .output()
        .expect("Rust fasta_check failed");

    assert_eq!(
        String::from_utf8_lossy(&cpp_out.stdout),
        String::from_utf8_lossy(&rust_out.stdout),
        "fasta_check DNA output differs"
    );
}

// --- Rust pipeline test ---

#[test]
fn test_rust_pipeline_protein_matches_expected() {
    let db_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amrfinder_db/2026-03-24.1");
    let input = test_data_dir().join("test_prot.fa");
    let gff = test_data_dir().join("test_prot.gff");
    let expected = test_data_dir().join("test_prot.expected");

    if !db_dir.exists() || !input.exists() || !expected.exists() {
        return;
    }

    // Check if the database is indexed
    if !db_dir.join("AMRProt.fa.phr").exists() {
        return;
    }

    let output = Command::new(rust_binary())
        .args([
            "run",
            "-p", input.to_str().unwrap(),
            "-g", gff.to_str().unwrap(),
            "-O", "Escherichia",
            "--plus",
            "--print_node",
            "-d", db_dir.to_str().unwrap(),
            "--threads", "6",
        ])
        .output()
        .expect("Rust pipeline failed");

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        eprintln!("Rust pipeline stderr: {}", stderr);
        panic!("Rust pipeline failed: {}", stderr);
    }

    let actual = String::from_utf8_lossy(&output.stdout);
    let expected_content = std::fs::read_to_string(&expected)
        .expect("Failed to read expected file");

    assert_eq!(
        actual, expected_content,
        "Rust pipeline protein output doesn't match expected"
    );
}
