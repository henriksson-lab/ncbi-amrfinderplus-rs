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

fn db_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amrfinder_db/2026-03-24.1")
}

fn require_pipeline_fixture(paths: &[PathBuf]) -> bool {
    let db = db_dir();
    db.exists() && db.join("AMRProt.fa.phr").exists() && paths.iter().all(|path| path.exists())
}

fn run_amrfinder(args: &[&str]) -> String {
    let output = Command::new(rust_binary())
        .args(args)
        .output()
        .expect("Rust pipeline failed");

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        panic!("Rust pipeline failed: {stderr}");
    }

    String::from_utf8(output.stdout).expect("pipeline output should be UTF-8")
}

fn run_cpp_amrfinder(args: &[&str]) -> String {
    let output = Command::new(cpp_binary("amrfinder"))
        .args(args)
        .output()
        .expect("C++ pipeline failed");

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        panic!("C++ pipeline failed: {stderr}");
    }

    String::from_utf8(output.stdout).expect("C++ pipeline output should be UTF-8")
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

#[test]
fn test_rust_pipeline_protein_current_regression() {
    let db_dir = db_dir();
    let input = test_data_dir().join("test_prot.fa");
    let gff = test_data_dir().join("test_prot.gff");

    if !db_dir.exists() || !input.exists() || !gff.exists() {
        return;
    }
    if !db_dir.join("AMRProt.fa.phr").exists() {
        return;
    }

    let output = Command::new(rust_binary())
        .args([
            "run",
            "-p",
            input.to_str().unwrap(),
            "-g",
            gff.to_str().unwrap(),
            "-O",
            "Escherichia",
            "--plus",
            "--print_node",
            "-d",
            db_dir.to_str().unwrap(),
            "--threads",
            "6",
        ])
        .output()
        .expect("Rust pipeline failed");

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        panic!("Rust pipeline failed: {stderr}");
    }

    let actual = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = actual.lines().collect();
    assert_eq!(
        lines.len(),
        15,
        "protein pipeline should emit header plus 14 rows"
    );
    assert!(actual.contains(
        "aph3pp-Ib_partial_5p_neg\tcontig09\t1\t675\t-\taph(3'')-Ib\taminoglycoside O-phosphotransferase APH(3'')-Ib\tcore\tAMR\tAMR\tAMINOGLYCOSIDE\tSTREPTOMYCIN\tPARTIAL_CONTIG_ENDP\t225\t267\t81.27\t100.00"
    ));
    assert!(actual.contains(
        "pmrB_C84R\tcontig14\t1093\t2181\t+\tpmrB_C84R\tEscherichia colistin resistant PmrB\tcore\tAMR\tPOINT\tCOLISTIN\tCOLISTIN\tPOINTP"
    ));
    assert!(actual.contains(
        "nfsA_R15C_K141STOP\tcontig16\t1\t423\t+\tnfsA_R15C\tEscherichia nitrofurantoin resistant NfsA\tcore\tAMR\tPOINT\tNITROFURAN\tNITROFURANTOIN\tPOINTP"
    ));
    assert!(actual.contains(
        "blaOXA-436_partial\tcontig03\t101\t802\t+\tblaOXA\tOXA-48 family class D beta-lactamase"
    ));
}

#[test]
fn test_rust_pipeline_nucleotide_matches_cpp_expected() {
    let db = db_dir();
    let input = test_data_dir().join("test_dna.fa");
    let expected = test_data_dir().join("test_dna.expected");
    if !require_pipeline_fixture(&[input.clone(), expected.clone()]) {
        return;
    }

    let mutation_all = tempfile::NamedTempFile::new().expect("mutation_all tempfile");
    let actual = run_amrfinder(&[
        "run",
        "-n",
        input.to_str().unwrap(),
        "-O",
        "Escherichia",
        "--plus",
        "--print_node",
        "-d",
        db.to_str().unwrap(),
        "--threads",
        "6",
        "--mutation_all",
        mutation_all.path().to_str().unwrap(),
    ]);
    let expected = std::fs::read_to_string(expected).expect("expected fixture should be readable");

    assert_eq!(
        actual, expected,
        "nucleotide pipeline output should match C++ expected fixture"
    );

    let mutation_all =
        std::fs::read_to_string(mutation_all.path()).expect("mutation_all should be written");
    assert!(
        mutation_all.contains("\tblaTEMp_C141C\tEscherichia blaTEM promoter region [WILDTYPE]\t")
            && mutation_all
                .contains("\tampC_T-14TGT\tEscherichia cephalosporin resistant ampC promoter\t"),
        "mutation_all should include wildtype and mutant DNA mutation rows"
    );
    assert!(
        mutation_all
            .lines()
            .next()
            .is_some_and(|header| header.ends_with("\tHierarchy node")),
        "mutation_all should preserve --print_node"
    );
}

#[test]
fn test_rust_pipeline_mutation_all_matches_cpp() {
    let db = db_dir();
    let cpp = cpp_binary("amrfinder");
    let input = test_data_dir().join("test_dna.fa");
    if !cpp.exists() || !require_pipeline_fixture(std::slice::from_ref(&input)) {
        return;
    }

    let cpp_mutation_all = tempfile::NamedTempFile::new().expect("C++ mutation_all tempfile");
    let rust_mutation_all = tempfile::NamedTempFile::new().expect("Rust mutation_all tempfile");

    let common = [
        "-n",
        input.to_str().unwrap(),
        "-O",
        "Escherichia",
        "--print_node",
        "-d",
        db.to_str().unwrap(),
        "--threads",
        "6",
    ];

    let mut cpp_args = common.to_vec();
    cpp_args.extend(["--mutation_all", cpp_mutation_all.path().to_str().unwrap()]);
    let _ = run_cpp_amrfinder(&cpp_args);

    let mut rust_args = vec!["run"];
    rust_args.extend(common);
    rust_args.extend(["--mutation_all", rust_mutation_all.path().to_str().unwrap()]);
    let _ = run_amrfinder(&rust_args);

    let cpp =
        std::fs::read_to_string(cpp_mutation_all.path()).expect("C++ mutation_all should exist");
    let rust =
        std::fs::read_to_string(rust_mutation_all.path()).expect("Rust mutation_all should exist");
    assert_eq!(
        rust, cpp,
        "pipeline --mutation_all side output should match C++ byte-for-byte"
    );
}

#[test]
fn test_rust_pipeline_protein_and_nucleotide_matches_cpp_expected() {
    let db = db_dir();
    let protein = test_data_dir().join("test_prot.fa");
    let nucleotide = test_data_dir().join("test_dna.fa");
    let gff = test_data_dir().join("test_prot.gff");
    let expected = test_data_dir().join("test_both.expected");
    if !require_pipeline_fixture(&[
        protein.clone(),
        nucleotide.clone(),
        gff.clone(),
        expected.clone(),
    ]) {
        return;
    }

    let actual = run_amrfinder(&[
        "run",
        "-p",
        protein.to_str().unwrap(),
        "-n",
        nucleotide.to_str().unwrap(),
        "-g",
        gff.to_str().unwrap(),
        "-O",
        "Escherichia",
        "--plus",
        "--print_node",
        "-d",
        db.to_str().unwrap(),
        "--threads",
        "6",
    ]);
    let expected = std::fs::read_to_string(expected).expect("expected fixture should be readable");

    assert_eq!(
        actual, expected,
        "combined protein+nucleotide pipeline output should match C++ expected fixture"
    );
}

#[test]
fn test_rust_pipeline_report_common_keeps_organism_suppressed_rows() {
    let db = db_dir();
    let input = test_data_dir().join("test_prot.fa");
    let gff = test_data_dir().join("test_prot.gff");
    if !require_pipeline_fixture(&[input.clone(), gff.clone()]) {
        return;
    }

    let actual = run_amrfinder(&[
        "run",
        "-p",
        input.to_str().unwrap(),
        "-g",
        gff.to_str().unwrap(),
        "-O",
        "Escherichia",
        "--plus",
        "--print_node",
        "--report_common",
        "-d",
        db.to_str().unwrap(),
        "--threads",
        "6",
    ]);

    assert!(
        actual.contains("\narsR-suppressed-in-escherichia\t"),
        "--report_common should retain rows listed in AMRProt-suppress.tsv"
    );
}

#[test]
fn test_rust_pipeline_campylobacter_matches_cpp() {
    let db = db_dir();
    let cpp = cpp_binary("amrfinder");
    let input = test_data_dir().join("test_prot.fa");
    if !cpp.exists() || !require_pipeline_fixture(std::slice::from_ref(&input)) {
        return;
    }

    let common = [
        "-p",
        input.to_str().unwrap(),
        "-O",
        "Campylobacter",
        "--plus",
        "--print_node",
        "-d",
        db.to_str().unwrap(),
        "--threads",
        "6",
    ];

    let cpp = run_cpp_amrfinder(&common);
    let mut rust_args = vec!["run"];
    rust_args.extend(common);
    let rust = run_amrfinder(&rust_args);

    assert_eq!(
        rust, cpp,
        "Campylobacter protein pipeline output should match C++ byte-for-byte"
    );
}

#[test]
fn test_rust_pipeline_does_not_synthesize_stx_operon_without_plus() {
    let db = db_dir();
    let input = test_data_dir().join("test_dna.fa");
    if !require_pipeline_fixture(std::slice::from_ref(&input)) {
        return;
    }

    let actual = run_amrfinder(&[
        "run",
        "-n",
        input.to_str().unwrap(),
        "-O",
        "Escherichia",
        "--print_node",
        "-d",
        db.to_str().unwrap(),
        "--threads",
        "6",
    ]);

    assert!(
        !actual.contains("\tstx2a_operon\t"),
        "stx2a_operon should only be synthesized for plus reports"
    );
}
