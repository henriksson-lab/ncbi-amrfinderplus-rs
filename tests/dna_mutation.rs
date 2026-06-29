use std::path::PathBuf;
use std::process::Command;

use amrfinder::dna_mutation::body;

#[test]
fn dna_mutation_matches_cpp() {
    let test_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    let blastn_file = test_dir.join("blastn");
    let db = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amrfinder_db/2026-03-24.1");
    let mutation_table = db.join("AMR_DNA-Escherichia.tsv");
    let expected_file = test_dir.join("dna_mutation_expected.tsv");

    if !blastn_file.exists() || !mutation_table.exists() || !expected_file.exists() {
        return;
    }

    let mut output = Vec::new();
    body(
        &blastn_file,
        &mutation_table,
        "Escherichia",
        None,
        "",
        true,
        &mut output,
    )
    .unwrap();

    let rust_output = String::from_utf8(output).unwrap();
    let cpp_output = std::fs::read_to_string(&expected_file).unwrap();

    assert_eq!(
        rust_output, cpp_output,
        "DNA mutation output must be byte-identical to the C++ golden output"
    );
}

#[test]
fn dna_mutation_all_matches_cpp_byte_for_byte() {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("amr/dna_mutation");
    let blastn_file = manifest.join("tests/golden/blastn");
    let mutation_table = manifest.join("amrfinder_db/2026-03-24.1/AMR_DNA-Escherichia.tsv");

    if !cpp.exists() || !blastn_file.exists() || !mutation_table.exists() {
        return;
    }

    let cpp_mut_all = tempfile::NamedTempFile::new().unwrap();
    let cpp_out = Command::new(cpp)
        .arg(&blastn_file)
        .arg(&mutation_table)
        .arg("Escherichia")
        .arg("-mutation_all")
        .arg(cpp_mut_all.path())
        .arg("-print_node")
        .output()
        .expect("C++ dna_mutation failed");
    assert!(
        cpp_out.status.success(),
        "C++ dna_mutation failed: {}",
        String::from_utf8_lossy(&cpp_out.stderr)
    );

    let mut rust_output = Vec::new();
    let rust_mut_all = tempfile::NamedTempFile::new().unwrap();
    body(
        &blastn_file,
        &mutation_table,
        "Escherichia",
        Some(rust_mut_all.path()),
        "",
        true,
        &mut rust_output,
    )
    .unwrap();

    let cpp_mut_all = std::fs::read_to_string(cpp_mut_all.path()).unwrap();
    let rust_mut_all = std::fs::read_to_string(rust_mut_all.path()).unwrap();

    assert_eq!(
        rust_mut_all, cpp_mut_all,
        "DNA mutation_all output must be byte-identical to C++"
    );
}

#[test]
fn ampc_insertion_replaces_generic_unknown_in_mutation_all() {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let blastn_file = manifest.join("tests/golden/blastn");
    let mutation_table = manifest.join("amrfinder_db/2026-03-24.1/AMR_DNA-Escherichia.tsv");

    if !blastn_file.exists() || !mutation_table.exists() {
        return;
    }

    let blastn = std::fs::read_to_string(&blastn_file).unwrap();
    let amp_c_line = blastn
        .lines()
        .find(|line| line.contains("\tcontig17\t"))
        .expect("golden BLASTN should contain the ampC insertion case");
    let focused_blastn = tempfile::NamedTempFile::new().unwrap();
    std::fs::write(focused_blastn.path(), format!("{amp_c_line}\n")).unwrap();

    let mut report = Vec::new();
    let mutation_all = tempfile::NamedTempFile::new().unwrap();
    body(
        focused_blastn.path(),
        &mutation_table,
        "Escherichia",
        Some(mutation_all.path()),
        "",
        true,
        &mut report,
    )
    .unwrap();

    let mutation_all = std::fs::read_to_string(mutation_all.path()).unwrap();
    assert!(
        mutation_all
            .contains("\tampC_T-14TGT\tEscherichia cephalosporin resistant ampC promoter\t"),
        "known ampC promoter insertion should be reported"
    );
    assert!(
        !mutation_all.contains("\tampC_C38CGT\t"),
        "known insertion should suppress the generic anchored unknown insertion artifact"
    );

    let lines: Vec<&str> = mutation_all.lines().collect();
    let c11 = lines
        .iter()
        .position(|line| line.contains("\tampC_C-11C\t"))
        .expect("wildtype ampC_C-11C row should be present");
    let insertion = lines
        .iter()
        .position(|line| line.contains("\tampC_T-14TGT\t"))
        .expect("ampC insertion row should be present");
    assert!(
        c11 < insertion,
        "C++ sorts empty wildtype SeqChanges before the real insertion change"
    );
}
