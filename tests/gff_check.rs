use std::io::Write;
use std::path::Path;

use amrfinder::gff::GffType;
use amrfinder::gff_check::{body, main, NO_FILE};

#[test]
fn no_file_exits_without_inputs() {
    body(
        NO_FILE,
        GffType::Genbank,
        Path::new(""),
        Path::new(""),
        Path::new(""),
        Path::new(""),
        false,
        false,
        &mut Vec::new(),
    )
    .unwrap();
}

#[test]
fn genbank_protein_match_uses_locus_tag() {
    let mut gff = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        gff,
        "ctg\tRefSeq\tCDS\t1\t3\t.\t+\t0\tID=cds1;locus_tag=LT1;gene=g;product=p"
    )
    .unwrap();
    let mut prot = tempfile::NamedTempFile::new().unwrap();
    writeln!(prot, ">prot1 [locus_tag=LT1]").unwrap();
    writeln!(prot, "MAA").unwrap();
    let prot_match = tempfile::NamedTempFile::new().unwrap();

    body(
        gff.path().to_str().unwrap(),
        GffType::Genbank,
        prot.path(),
        Path::new(""),
        prot_match.path(),
        Path::new(""),
        false,
        false,
        &mut Vec::new(),
    )
    .unwrap();

    assert_eq!(
        std::fs::read_to_string(prot_match.path()).unwrap(),
        "prot1\tLT1\n"
    );
}

#[test]
fn verbose_reports_protein_count() {
    let mut gff = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        gff,
        "ctg\tRefSeq\tCDS\t1\t3\t.\t+\t0\tID=cds1;locus_tag=LT1;gene=g;product=p"
    )
    .unwrap();
    let mut prot = tempfile::NamedTempFile::new().unwrap();
    writeln!(prot, ">prot1 [locus_tag=LT1]").unwrap();
    writeln!(prot, "MAA").unwrap();
    let prot_match = tempfile::NamedTempFile::new().unwrap();

    let mut out = Vec::new();
    body(
        gff.path().to_str().unwrap(),
        GffType::Genbank,
        prot.path(),
        Path::new(""),
        prot_match.path(),
        Path::new(""),
        false,
        true,
        &mut out,
    )
    .unwrap();

    assert_eq!(String::from_utf8(out).unwrap(), "# Proteins in GFF: 1\n");
}

#[test]
fn main_parses_gfftype_and_no_file() {
    assert_eq!(
        main(vec![
            "gff_check".into(),
            NO_FILE.into(),
            "-gfftype".into(),
            "pgap".into(),
            "-lcl".into(),
        ])
        .unwrap(),
        0
    );
}
