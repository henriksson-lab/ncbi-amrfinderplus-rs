use std::io::Write;

use amrfinder::mutate::{body, main};

#[test]
fn mutates_protein_and_keeps_original_when_requested() {
    let mut fasta = tempfile::NamedTempFile::new().unwrap();
    let mut mutations = tempfile::NamedTempFile::new().unwrap();
    writeln!(fasta, ">seq1 comment").unwrap();
    writeln!(fasta, "MST").unwrap();
    writeln!(mutations, "seq1 2 gene_S2L gene_S2L").unwrap();

    let mut out = Vec::new();
    body(fasta.path(), mutations.path(), true, true, &mut out).unwrap();

    assert_eq!(
        String::from_utf8(out).unwrap(),
        ">seq1 comment\nMST\n>seq1 comment:2:gene_S2L\nMLT\n"
    );
}

#[test]
fn mutates_dna_and_lowercases_result() {
    let mut fasta = tempfile::NamedTempFile::new().unwrap();
    let mut mutations = tempfile::NamedTempFile::new().unwrap();
    writeln!(fasta, ">dna1").unwrap();
    writeln!(fasta, "ACGT").unwrap();
    writeln!(mutations, "dna1 2 gene_C2T gene_C2T").unwrap();

    let mut out = Vec::new();
    body(fasta.path(), mutations.path(), false, false, &mut out).unwrap();

    assert_eq!(String::from_utf8(out).unwrap(), ">dna1:2:gene_C2T\natgt\n");
}

#[test]
fn main_parses_flags_and_writes_to_output() {
    let mut fasta = tempfile::NamedTempFile::new().unwrap();
    let mut mutations = tempfile::NamedTempFile::new().unwrap();
    writeln!(fasta, ">seq1").unwrap();
    writeln!(fasta, "MST").unwrap();
    writeln!(mutations, "seq1 2 gene_S2L gene_S2L").unwrap();

    let mut out = Vec::new();
    assert_eq!(
        main(
            vec![
                "mutate".into(),
                fasta.path().into(),
                mutations.path().into(),
                "-aa".into(),
            ],
            &mut out,
        )
        .unwrap(),
        0
    );

    assert_eq!(String::from_utf8(out).unwrap(), ">seq1:2:gene_S2L\nMLT\n");
}

#[test]
fn wraps_saved_sequences_at_fasta_line_length() {
    let mut fasta = tempfile::NamedTempFile::new().unwrap();
    let mut mutations = tempfile::NamedTempFile::new().unwrap();
    let seq = "A".repeat(81);
    writeln!(fasta, ">seq1").unwrap();
    writeln!(fasta, "{seq}").unwrap();
    writeln!(mutations, "seq1 81 gene_A81T gene_A81T").unwrap();

    let mut out = Vec::new();
    body(fasta.path(), mutations.path(), true, true, &mut out).unwrap();

    assert_eq!(
        String::from_utf8(out).unwrap(),
        format!(
            ">seq1\n{}\nA\n>seq1:81:gene_A81T\n{}\nT\n",
            "A".repeat(80),
            "A".repeat(80)
        )
    );
}

#[test]
fn rejects_bad_sequence_character_with_sequence_context() {
    let mut fasta = tempfile::NamedTempFile::new().unwrap();
    let mut mutations = tempfile::NamedTempFile::new().unwrap();
    writeln!(fasta, ">seq1 comment").unwrap();
    writeln!(fasta, "ACGZ").unwrap();
    writeln!(mutations, "seq1 2 gene_C2T gene_C2T").unwrap();

    let mut out = Vec::new();
    let err = body(fasta.path(), mutations.path(), false, false, &mut out).unwrap_err();

    assert!(err.to_string().contains("seq1 comment"));
    assert!(err.to_string().contains("Bad sequence character"));
}

#[test]
fn mutation_apply_errors_include_sequence_name() {
    let mut fasta = tempfile::NamedTempFile::new().unwrap();
    let mut mutations = tempfile::NamedTempFile::new().unwrap();
    writeln!(fasta, ">seq1 comment").unwrap();
    writeln!(fasta, "MST").unwrap();
    writeln!(mutations, "seq1 9 gene_S9L gene_S9L").unwrap();

    let mut out = Vec::new();
    let err = body(fasta.path(), mutations.path(), true, false, &mut out).unwrap_err();

    assert!(err.to_string().contains("seq1 comment"));
}

#[test]
fn blank_line_ends_sequence_and_progresses_to_next_header() {
    let mut fasta = tempfile::NamedTempFile::new().unwrap();
    let mut mutations = tempfile::NamedTempFile::new().unwrap();
    writeln!(fasta, ">seq1").unwrap();
    writeln!(fasta, "AC").unwrap();
    writeln!(fasta).unwrap();
    writeln!(fasta, ">seq2").unwrap();
    writeln!(fasta, "GT").unwrap();
    writeln!(mutations, "seq1 2 gene_C2T gene_C2T").unwrap();
    writeln!(mutations, "seq2 1 gene_G1A gene_G1A").unwrap();

    let mut out = Vec::new();
    body(fasta.path(), mutations.path(), false, false, &mut out).unwrap();

    assert_eq!(
        String::from_utf8(out).unwrap(),
        ">seq1:2:gene_C2T\nat\n>seq2:1:gene_G1A\nat\n"
    );
}

#[test]
fn rejects_non_header_after_blank_separator() {
    let mut fasta = tempfile::NamedTempFile::new().unwrap();
    let mut mutations = tempfile::NamedTempFile::new().unwrap();
    writeln!(fasta, ">seq1").unwrap();
    writeln!(fasta, "AC").unwrap();
    writeln!(fasta).unwrap();
    writeln!(fasta, "not-a-header").unwrap();
    writeln!(mutations, "seq1 2 gene_C2T gene_C2T").unwrap();

    let mut out = Vec::new();
    let err = body(fasta.path(), mutations.path(), false, false, &mut out).unwrap_err();

    assert!(err.to_string().contains("Error in Multifasta"));
}

#[test]
fn fasta_name_tabs_are_replaced_with_spaces() {
    let mut fasta = tempfile::NamedTempFile::new().unwrap();
    let mut mutations = tempfile::NamedTempFile::new().unwrap();
    writeln!(fasta, ">seq1\tcomment").unwrap();
    writeln!(fasta, "MST").unwrap();
    writeln!(mutations, "seq1 2 gene_S2L gene_S2L").unwrap();

    let mut out = Vec::new();
    body(fasta.path(), mutations.path(), true, true, &mut out).unwrap();

    assert_eq!(
        String::from_utf8(out).unwrap(),
        ">seq1 comment\nMST\n>seq1 comment:2:gene_S2L\nMLT\n"
    );
}

#[test]
fn sparse_gaps_are_removed_before_sequence_qc() {
    let mut fasta = tempfile::NamedTempFile::new().unwrap();
    let mut mutations = tempfile::NamedTempFile::new().unwrap();
    writeln!(fasta, ">seq1").unwrap();
    writeln!(fasta, "A-C").unwrap();
    writeln!(mutations, "seq1 2 gene_C2T gene_C2T").unwrap();

    let mut out = Vec::new();
    body(fasta.path(), mutations.path(), false, true, &mut out).unwrap();

    assert_eq!(
        String::from_utf8(out).unwrap(),
        ">seq1\nac\n>seq1:2:gene_C2T\nat\n"
    );
}
