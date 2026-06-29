use std::io::Write;

use amrfinder::fasta_extract::body;

#[test]
fn dna_header_uses_ascii_range_separator() {
    let mut fasta = tempfile::NamedTempFile::new().unwrap();
    let mut target = tempfile::NamedTempFile::new().unwrap();
    writeln!(fasta, ">seq1").unwrap();
    writeln!(fasta, "ACGTACGT").unwrap();
    writeln!(target, "seq1 2 5 + gene product").unwrap();

    let mut out = Vec::new();
    body(fasta.path(), target.path(), false, false, &mut out).unwrap();

    let out = String::from_utf8(out).unwrap();
    assert!(out.starts_with(">seq1:2-5 strand:+ gene product\n"));
}

#[test]
fn dna_stop_is_clamped_before_slicing() {
    let mut fasta = tempfile::NamedTempFile::new().unwrap();
    let mut target = tempfile::NamedTempFile::new().unwrap();
    writeln!(fasta, ">seq1").unwrap();
    writeln!(fasta, "ACGT").unwrap();
    writeln!(target, "seq1 2 10 + gene product").unwrap();

    let mut out = Vec::new();
    body(fasta.path(), target.path(), false, false, &mut out).unwrap();

    let out = String::from_utf8(out).unwrap();
    assert_eq!(out, ">seq1:2-4 strand:+ gene product\nCGT\n");
}

#[test]
fn target_product_name_preserves_internal_whitespace() {
    let mut fasta = tempfile::NamedTempFile::new().unwrap();
    let mut target = tempfile::NamedTempFile::new().unwrap();
    writeln!(fasta, ">seq1").unwrap();
    writeln!(fasta, "ACGT").unwrap();
    writeln!(target, "seq1 1 4 + gene   product\twith  spacing").unwrap();

    let mut out = Vec::new();
    body(fasta.path(), target.path(), false, false, &mut out).unwrap();

    let out = String::from_utf8(out).unwrap();
    assert!(out.starts_with(">seq1:1-4 strand:+ gene product\twith  spacing\n"));
}

#[test]
fn blank_target_line_is_rejected() {
    let mut fasta = tempfile::NamedTempFile::new().unwrap();
    let mut target = tempfile::NamedTempFile::new().unwrap();
    writeln!(fasta, ">seq1").unwrap();
    writeln!(fasta, "ACGT").unwrap();
    writeln!(target).unwrap();

    let mut out = Vec::new();
    let err = body(fasta.path(), target.path(), false, false, &mut out).unwrap_err();

    assert!(err.to_string().contains("Empty target line"));
}

#[test]
fn bad_reverse_complement_nucleotide_returns_error() {
    let mut fasta = tempfile::NamedTempFile::new().unwrap();
    let mut target = tempfile::NamedTempFile::new().unwrap();
    writeln!(fasta, ">seq1").unwrap();
    writeln!(fasta, "ACGZ").unwrap();
    writeln!(target, "seq1 1 4 - gene product").unwrap();

    let mut out = Vec::new();
    let err = body(fasta.path(), target.path(), false, false, &mut out).unwrap_err();

    assert!(err.to_string().contains("Bad nucleotide Z"));
}
