use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amr")
}

fn cpp_binary(name: &str) -> PathBuf {
    test_data_dir().join(name)
}

#[test]
fn fasta_check_protein_accepts_golden() {
    let path = test_data_dir().join("test_prot.fa");
    if !path.exists() {
        return;
    }
    let result = amrfinder::fasta_check::body(&path, true, false, false, 0, false, None, None);
    assert!(result.is_ok(), "fasta_check failed: {:?}", result.err());
}

#[test]
fn fasta_check_dna_accepts_golden() {
    let path = test_data_dir().join("test_dna.fa");
    if !path.exists() {
        return;
    }
    let result = amrfinder::fasta_check::body(&path, false, false, true, 0, false, None, None);
    assert!(result.is_ok(), "fasta_check failed: {:?}", result.err());
}

#[test]
fn fasta_check_protein_matches_cpp() {
    let path = test_data_dir().join("test_prot.fa");
    let cpp_bin = cpp_binary("fasta_check");
    if !path.exists() || !cpp_bin.exists() {
        return;
    }

    let cpp_output = Command::new(&cpp_bin)
        .arg(&path)
        .arg("-aa")
        .output()
        .expect("failed to run C++ fasta_check");
    let cpp_stdout = String::from_utf8_lossy(&cpp_output.stdout);

    let (num_seqs, max_len, total_len) =
        amrfinder::fasta_check::body(&path, true, false, false, 0, false, None, None).unwrap();
    let rust_stdout = format!("{}\n{}\n{}\n", num_seqs, max_len, total_len);

    assert_eq!(
        cpp_stdout, rust_stdout,
        "C++ and Rust fasta_check output differ for protein"
    );
}

#[test]
fn fasta_check_dna_matches_cpp() {
    let path = test_data_dir().join("test_dna.fa");
    let cpp_bin = cpp_binary("fasta_check");
    if !path.exists() || !cpp_bin.exists() {
        return;
    }

    let cpp_output = Command::new(&cpp_bin)
        .arg(&path)
        .arg("-ambig")
        .output()
        .expect("failed to run C++ fasta_check");
    let cpp_stdout = String::from_utf8_lossy(&cpp_output.stdout);

    let (num_seqs, max_len, total_len) =
        amrfinder::fasta_check::body(&path, false, false, true, 0, false, None, None).unwrap();
    let rust_stdout = format!("{}\n{}\n{}\n", num_seqs, max_len, total_len);

    assert_eq!(
        cpp_stdout, rust_stdout,
        "C++ and Rust fasta_check output differ for DNA"
    );
}

#[test]
fn fasta_check_error_line_numbers_are_one_based() {
    let mut fasta = tempfile::NamedTempFile::new().unwrap();
    writeln!(fasta, ">seq1").unwrap();
    writeln!(fasta, "ACGT").unwrap();
    writeln!(fasta, "Z").unwrap();

    let err = amrfinder::fasta_check::body(fasta.path(), false, false, true, 0, false, None, None)
        .unwrap_err();

    assert!(err.to_string().contains("line 3"));
}

#[test]
fn fasta_check_output_preserves_non_utf8_header_comment_bytes() {
    let mut fasta = tempfile::NamedTempFile::new().unwrap();
    fasta.write_all(b">seq1 \xFFcomment\nACGT\n").unwrap();
    let out = tempfile::NamedTempFile::new().unwrap();

    amrfinder::fasta_check::body(
        fasta.path(),
        false,
        false,
        true,
        0,
        false,
        None,
        Some(out.path()),
    )
    .unwrap();

    let output = std::fs::read(out.path()).unwrap();
    assert!(output.starts_with(b">seq1 \xFFcomment\n"));
}

#[test]
fn fasta2parts_preserves_non_utf8_sequence_bytes() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("input.fa");
    let out_dir = dir.path().join("parts");
    std::fs::create_dir(&out_dir).unwrap();
    let mut file = std::fs::File::create(&fasta).unwrap();
    file.write_all(b">seq1\nAC\xFFT\t \n>seq2\nGG\n").unwrap();

    amrfinder::fasta2parts::body(&fasta, 2, &out_dir).unwrap();

    let first = std::fs::read(out_dir.join("1")).unwrap();
    assert!(
        first.windows(4).any(|window| window == b"AC\xFFT"),
        "non-UTF-8 sequence bytes should be preserved"
    );
    assert!(
        !first.windows(2).any(|window| window == b"\t\n"),
        "C-space trailing bytes should be trimmed"
    );
}
