use std::fs;

use amrfinder::disruption2genesymbol::{body, main, SymbolRaw};
use amrfinder::seq::DisruptionType;
use tempfile::tempdir;

#[test]
fn symbol_raw_new_parses_genesymbol_raw_and_preserves_rest() {
    let symbol = SymbolRaw::new("contig1 protA fs_4_5_30_90_0_STOP extra fields").unwrap();

    assert_eq!(symbol.contig, "contig1");
    assert_eq!(symbol.prot, "protA");
    assert_eq!(symbol.dtype, DisruptionType::Frameshift);
    assert_eq!(symbol.qstart, 4);
    assert_eq!(symbol.qend, 5);
    assert_eq!(symbol.sstart, 30);
    assert_eq!(symbol.send, 90);
    assert_eq!(symbol.strand, -1);
    assert!(symbol.stop);
    assert_eq!(symbol.rest, " extra fields");
}

#[test]
fn save_text_formats_deletion_symbol_and_raw_disruption() {
    let mut symbol = SymbolRaw::new("contig1 protA del_1_3_0_0_1 tail").unwrap();
    symbol.ref_ = "ABC".to_string();

    let mut out = Vec::new();
    symbol.save_text(&mut out, false).unwrap();

    assert_eq!(
        String::from_utf8(out).unwrap(),
        "contig1\tprotA\tABC2del\tdel_1_3_0_0_1\t tail\n"
    );
}

#[test]
fn contig2aa_translates_plus_and_minus_strands() {
    let plus = SymbolRaw::new("contig1 protA ins_0_0_0_6_1").unwrap();
    assert_eq!(plus.contig2aa("ATGTAA", 0, 11).unwrap(), 'M');
    assert_eq!(plus.contig2aa("ATGTAA", 1, 11).unwrap(), '*');
    assert_eq!(plus.contig2aa("ATGTAA", 2, 11).unwrap(), '?');

    let minus = SymbolRaw::new("contig1 protA ins_0_0_0_6_0").unwrap();
    assert_eq!(minus.contig2aa("TTACAT", 0, 11).unwrap(), 'M');
    assert_eq!(minus.contig2aa("TTACAT", 1, 11).unwrap(), '*');
    assert_eq!(minus.contig2aa("TTACAT", 2, 11).unwrap(), '?');
}

#[test]
fn body_reads_fasta_inputs_and_inserts_gene_symbol() {
    let dir = tempdir().unwrap();
    let nucl = dir.path().join("nucl.fa");
    let prot = dir.path().join("prot.fa");
    let tab = dir.path().join("tab.tsv");

    fs::write(&nucl, ">contig1 description\nATGTAA\n").unwrap();
    fs::write(&prot, ">protA description\nMK\n").unwrap();
    fs::write(&tab, "contig1 protA fs_0_1_0_6_1_STOP more\n").unwrap();

    let mut out = Vec::new();
    body(&nucl, &prot, &tab, 11, 0, &mut out).unwrap();

    assert_eq!(
        String::from_utf8(out).unwrap(),
        "contig1\tprotA\tM1MfsTer1\tfs_0_1_0_6_1_STOP\t more\n"
    );
}

#[test]
fn main_parses_gencode_and_prot_id_position() {
    let dir = tempdir().unwrap();
    let nucl = dir.path().join("nucl.fa");
    let prot = dir.path().join("prot.fa");
    let tab = dir.path().join("tab.tsv");

    fs::write(&nucl, ">contig1 description\nATGTAA\n").unwrap();
    fs::write(&prot, ">acc|protA|gene description\nMK\n").unwrap();
    fs::write(&tab, "contig1 protA fs_0_1_0_6_1_STOP more\n").unwrap();

    let mut out = Vec::new();
    assert_eq!(
        main(
            vec![
                "disruption2genesymbol".into(),
                nucl.into(),
                prot.into(),
                tab.into(),
                "-gencode".into(),
                "11".into(),
                "-prot_id_pos".into(),
                "2".into(),
            ],
            &mut out,
        )
        .unwrap(),
        0
    );

    assert_eq!(
        String::from_utf8(out).unwrap(),
        "contig1\tprotA\tM1MfsTer1\tfs_0_1_0_6_1_STOP\t more\n"
    );
}

#[test]
fn body_normalizes_protein_reference_to_uppercase() {
    let dir = tempdir().unwrap();
    let nucl = dir.path().join("nucl.fa");
    let prot = dir.path().join("prot.fa");
    let tab = dir.path().join("tab.tsv");

    fs::write(&nucl, ">contig1\nATGTAA\n").unwrap();
    fs::write(&prot, ">protA\nmk\n").unwrap();
    fs::write(&tab, "contig1 protA fs_0_1_0_6_1_STOP more\n").unwrap();

    let mut out = Vec::new();
    body(&nucl, &prot, &tab, 11, 0, &mut out).unwrap();

    assert!(String::from_utf8(out).unwrap().contains("M1MfsTer1"));
}

#[test]
fn body_rejects_bad_nucleotide_character() {
    let dir = tempdir().unwrap();
    let nucl = dir.path().join("nucl.fa");
    let prot = dir.path().join("prot.fa");
    let tab = dir.path().join("tab.tsv");

    fs::write(&nucl, ">contig1\nATZ\n").unwrap();
    fs::write(&prot, ">protA\nMK\n").unwrap();
    fs::write(&tab, "contig1 protA fs_0_1_0_3_1 more\n").unwrap();

    let mut out = Vec::new();
    let err = body(&nucl, &prot, &tab, 11, 0, &mut out).unwrap_err();

    assert!(err.to_string().contains("Bad sequence character"));
}

#[test]
fn body_rejects_sequence_text_after_blank_line() {
    let dir = tempdir().unwrap();
    let nucl = dir.path().join("nucl.fa");
    let prot = dir.path().join("prot.fa");
    let tab = dir.path().join("tab.tsv");

    fs::write(&nucl, ">contig1\nATG\n\nTAA\n").unwrap();
    fs::write(&prot, ">protA\nMK\n").unwrap();
    fs::write(&tab, "contig1 protA fs_0_1_0_3_1 more\n").unwrap();

    let mut out = Vec::new();
    let err = body(&nucl, &prot, &tab, 11, 0, &mut out).unwrap_err();

    assert!(err.to_string().contains("Error in Multifasta, line 5"));
}

#[test]
fn body_reports_protein_multifasta_line_number_after_blank_line() {
    let dir = tempdir().unwrap();
    let nucl = dir.path().join("nucl.fa");
    let prot = dir.path().join("prot.fa");
    let tab = dir.path().join("tab.tsv");

    fs::write(&nucl, ">contig1\nATG\n").unwrap();
    fs::write(&prot, ">protA\nMK\n\nnot-a-header\n").unwrap();
    fs::write(&tab, "contig1 protA fs_0_1_0_3_1 more\n").unwrap();

    let mut out = Vec::new();
    let err = body(&nucl, &prot, &tab, 11, 0, &mut out).unwrap_err();

    assert!(err.to_string().contains("Error in Multifasta, line 5"));
}

#[test]
fn body_does_not_trim_leading_sequence_whitespace() {
    let dir = tempdir().unwrap();
    let nucl = dir.path().join("nucl.fa");
    let prot = dir.path().join("prot.fa");
    let tab = dir.path().join("tab.tsv");

    fs::write(&nucl, ">contig1\n ATG\n").unwrap();
    fs::write(&prot, ">protA\nMK\n").unwrap();
    fs::write(&tab, "contig1 protA fs_0_1_0_3_1 more\n").unwrap();

    let mut out = Vec::new();
    let err = body(&nucl, &prot, &tab, 11, 0, &mut out).unwrap_err();

    assert!(err.to_string().contains("Bad sequence character"));
}

#[test]
fn body_normalizes_fasta_header_tabs_like_seq_constructor() {
    let dir = tempdir().unwrap();
    let nucl = dir.path().join("nucl.fa");
    let prot = dir.path().join("prot.fa");
    let tab = dir.path().join("tab.tsv");

    fs::write(&nucl, ">contig1\tdescription\nATGTAA\n").unwrap();
    fs::write(&prot, ">protA\tdescription\nMK\n").unwrap();
    fs::write(&tab, "contig1 protA fs_0_1_0_6_1_STOP more\n").unwrap();

    let mut out = Vec::new();
    body(&nucl, &prot, &tab, 11, 0, &mut out).unwrap();

    assert_eq!(
        String::from_utf8(out).unwrap(),
        "contig1\tprotA\tM1MfsTer1\tfs_0_1_0_6_1_STOP\t more\n"
    );
}
