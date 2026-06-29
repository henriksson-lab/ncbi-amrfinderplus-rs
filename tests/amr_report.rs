use std::collections::{BTreeMap, HashMap};
use std::io::Write;
use std::path::PathBuf;

use amrfinder::amr_report::{
    load_blast_results, load_hmm_results, run_amr_report, AmrReportConfig, Batch, BlastAlignment,
    BlastRule, Fam, HmmAlignment, Susceptible,
};

fn db_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amrfinder_db/2026-03-24.1")
}

fn empty_batch() -> Batch {
    Batch {
        fam_map: HashMap::new(),
        hmm2fam: HashMap::new(),
        reportable_min: 0,
        ident_min: -1.0,
        coverage_min: 0.5,
        report_all_equal: false,
        print_node_raw: false,
        name: String::new(),
        cds_exist: false,
        suppress_prots: Vec::new(),
        alien_prots: Vec::new(),
        dna_lens: HashMap::new(),
        blast_als: Vec::new(),
        hmm_als: Vec::new(),
        hmm_exist: false,
        domains: HashMap::new(),
        accession2mutations: HashMap::new(),
        accession2susceptible: HashMap::new(),
        target2blast_als: BTreeMap::new(),
        target2good_blast_als: BTreeMap::new(),
        target2hmm_als: BTreeMap::new(),
        target2good_hmm_als: BTreeMap::new(),
    }
}

#[test]
fn test_blast_rule() {
    let br = BlastRule::new(0.9, 0.8);
    assert!(!br.empty());
    assert_eq!(br.ident, 0.9);
    assert_eq!(br.target_coverage, 0.0);
    assert_eq!(br.ref_coverage, 0.8);

    let empty = BlastRule::default();
    assert!(empty.empty());

    let cxx_empty = BlastRule::new(0.0, 90.0);
    assert!(cxx_empty.empty());
}

#[test]
fn save_text_methods_match_cpp_field_order() {
    let br = BlastRule::new(90.0, 80.0);
    let mut out = Vec::new();
    br.save_text(&mut out).unwrap();
    assert_eq!(String::from_utf8(out).unwrap(), "0.9 0.8");

    let fam = Fam {
        id: "fam".to_string(),
        genesymbol: "gene".to_string(),
        family_name: "Family".to_string(),
        reportable: 2,
        hmm: "HMM".to_string(),
        tc1: 10.0,
        tc2: 5.0,
        complete_br: br.clone(),
        partial_br: br.clone(),
        type_: "AMR".to_string(),
        subtype: "AMR".to_string(),
        class: "CLASS".to_string(),
        subclass: "SUBCLASS".to_string(),
        parent_id: String::new(),
    };
    let mut out = Vec::new();
    fam.save_text(&mut out).unwrap();
    assert_eq!(String::from_utf8(out).unwrap(), "HMM 10 5 Family 2");

    let susceptible = Susceptible {
        genesymbol: "gene".to_string(),
        cutoff: 0.8,
        class: "CLASS".to_string(),
        subclass: "SUBCLASS".to_string(),
        name: "Susceptible Name".to_string(),
    };
    let mut out = Vec::new();
    susceptible.save_text(&mut out).unwrap();
    assert_eq!(
        String::from_utf8(out).unwrap(),
        "gene\t0.8\tCLASS\tSUBCLASS\tSusceptible Name\n"
    );

    let hmm_al = HmmAlignment {
        sseqid: "target".to_string(),
        score1: 42.0,
        score2: 21.0,
        fam_id: "fam".to_string(),
        domain: None,
        blast_al_idx: None,
    };
    let fam_map = HashMap::from([("fam".to_string(), fam)]);
    let mut out = Vec::new();
    hmm_al.save_text(&mut out, &fam_map).unwrap();
    assert_eq!(String::from_utf8(out).unwrap(), "target 42 21 HMM");
}

#[test]
fn test_fam_loading() {
    let fam_path = db_dir().join("fam.tsv");
    if !fam_path.exists() {
        return;
    }

    let batch = Batch::from_fam_file(&fam_path, 0).unwrap();
    assert!(!batch.fam_map.is_empty(), "FAM map should not be empty");
    assert!(!batch.hmm2fam.is_empty(), "HMM map should not be empty");
    assert!(
        batch.fam_map.contains_key("blaTEM"),
        "Should have blaTEM family"
    );
}

#[test]
fn fam_loader_normalizes_null_family_name() {
    let mut fam_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        fam_file,
        "fam\t\tgene\t-\t0\t0\t90\t0\t90\t90\t0\t50\t2\tAMR\tSTYPE\tCLASS\tSUBCLASS\tNULL"
    )
    .unwrap();

    let batch = Batch::from_fam_file(fam_file.path(), 0).unwrap();

    assert_eq!(batch.fam_map["fam"].family_name, "");
}

#[test]
fn fam_loader_matches_cpp_rule_defaults_and_parent_pass() {
    let mut fam_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        fam_file,
        "root\t\troot_gene\tHMM_ROOT\t40\t20\t0\t0\t0\t0\t0\t0\t2\tAMR\tSTYPE\tCLASS\tSUBCLASS\tRoot"
    )
    .unwrap();
    writeln!(
        fam_file,
        "child\troot\tchild_gene\tHMM_CHILD\t50\t25\t95\t10\t70\t80\t20\t30\t1\tAMR\tSTYPE\tCLASS\tSUBCLASS\tChild"
    )
    .unwrap();
    writeln!(
        fam_file,
        "onepct\troot\tone_gene\tHMM_ONE\t50\t25\t1\t1\t1\t1\t1\t1\t1\tAMR\tSTYPE\tCLASS\tSUBCLASS\tOne"
    )
    .unwrap();

    let batch = Batch::from_fam_file(fam_file.path(), 0).unwrap();

    let root = &batch.fam_map["root"];
    assert_eq!(root.complete_br.ident, 0.9);
    assert_eq!(root.complete_br.ref_coverage, 0.9);
    assert_eq!(root.partial_br.ident, 0.9);
    assert_eq!(root.partial_br.ref_coverage, 0.5);

    let child = &batch.fam_map["child"];
    assert_eq!(child.parent_id, "root");
    assert_eq!(child.complete_br.ident, 0.95);
    assert_eq!(child.complete_br.target_coverage, 0.10);
    assert_eq!(child.complete_br.ref_coverage, 0.9);
    assert_eq!(child.partial_br.ident, 0.80);
    assert_eq!(child.partial_br.target_coverage, 0.20);
    assert_eq!(child.partial_br.ref_coverage, 0.5);
    assert_eq!(batch.hmm2fam["HMM_CHILD"], "child");

    let onepct = &batch.fam_map["onepct"];
    assert_eq!(onepct.complete_br.ident, 0.01);
    assert_eq!(onepct.complete_br.target_coverage, 0.01);
    assert_eq!(onepct.partial_br.ident, 0.01);
    assert_eq!(onepct.partial_br.target_coverage, 0.01);
}

#[test]
fn fam_loader_rejects_cpp_constructor_errors() {
    let mut duplicate = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        duplicate,
        "fam\t\tgene\t-\t0\t0\t0\t0\t0\t0\t0\t0\t2\tAMR\tSTYPE\tCLASS\tSUBCLASS\tFamily"
    )
    .unwrap();
    writeln!(
        duplicate,
        "fam\t\tgene\t-\t0\t0\t0\t0\t0\t0\t0\t0\t2\tAMR\tSTYPE\tCLASS\tSUBCLASS\tFamily"
    )
    .unwrap();
    let err = match Batch::from_fam_file(duplicate.path(), 0) {
        Ok(_) => panic!("duplicate family should fail"),
        Err(err) => err,
    };
    assert!(err.to_string().contains("Family fam is duplicated"));

    let mut missing_parent = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        missing_parent,
        "child\tmissing\tgene\t-\t0\t0\t0\t0\t0\t0\t0\t0\t2\tAMR\tSTYPE\tCLASS\tSUBCLASS\tFamily"
    )
    .unwrap();
    let err = match Batch::from_fam_file(missing_parent.path(), 0) {
        Ok(_) => panic!("missing parent should fail"),
        Err(err) => err,
    };
    assert!(err.to_string().contains("parentFamId \"missing\""));

    let mut malformed = tempfile::NamedTempFile::new().unwrap();
    writeln!(malformed, "fam\ttoo_short").unwrap();
    let err = match Batch::from_fam_file(malformed.path(), 0) {
        Ok(_) => panic!("malformed row should fail"),
        Err(err) => err,
    };
    assert!(err.to_string().contains("expected at least 13"));
}

#[test]
fn test_mutation_loading() {
    let db = db_dir();
    let fam_path = db.join("fam.tsv");
    let mut_path = db.join("AMRProt-mutation.tsv");
    if !fam_path.exists() || !mut_path.exists() {
        return;
    }

    let mut batch = Batch::from_fam_file(&fam_path, 0).unwrap();
    batch.load_mutations(&mut_path, "Escherichia").unwrap();
    assert!(
        !batch.accession2mutations.is_empty(),
        "Should have mutations for Escherichia"
    );
}

#[test]
fn mutation_and_susceptible_loaders_normalize_organism_underscores() {
    let mut mutation_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        mutation_file,
        "Escherichia_coli WP_1 1 gene_A1T gene_A1T CLASS SUBCLASS Name"
    )
    .unwrap();
    let mut susceptible_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        susceptible_file,
        "Escherichia_coli gene WP_2 80 CLASS SUBCLASS Susceptible_Name"
    )
    .unwrap();

    let mut batch = empty_batch();
    batch
        .load_mutations(mutation_file.path(), "Escherichia coli")
        .unwrap();
    batch
        .load_susceptible(susceptible_file.path(), "Escherichia coli")
        .unwrap();

    assert!(batch.accession2mutations.contains_key("WP_1"));
    let susceptible = &batch.accession2susceptible["WP_2"];
    assert_eq!(susceptible.cutoff, 0.8);
    assert_eq!(susceptible.name, "Susceptible Name");
}

#[test]
fn susceptible_loader_divides_cutoff_percentage_like_cpp() {
    let mut susceptible_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        susceptible_file,
        "Escherichia gene WP_2 75 CLASS SUBCLASS Susceptible_Name"
    )
    .unwrap();

    let mut batch = empty_batch();
    batch
        .load_susceptible(susceptible_file.path(), "Escherichia")
        .unwrap();

    assert_eq!(batch.accession2susceptible["WP_2"].cutoff, 0.75);
}

#[test]
fn susceptible_loader_divides_fractional_cutoff_like_cpp() {
    let mut susceptible_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        susceptible_file,
        "Escherichia gene WP_2 0.75 CLASS SUBCLASS Susceptible_Name"
    )
    .unwrap();

    let mut batch = empty_batch();
    batch
        .load_susceptible(susceptible_file.path(), "Escherichia")
        .unwrap();

    assert_eq!(batch.accession2susceptible["WP_2"].cutoff, 0.0075);
}

#[test]
fn mutation_loader_rejects_duplicate_mutations_like_cpp() {
    let mut mutation_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        mutation_file,
        "Escherichia_coli WP_1 1 gene_A1T gene_A1T CLASS SUBCLASS Name"
    )
    .unwrap();
    writeln!(
        mutation_file,
        "Escherichia_coli WP_1 1 gene_A1T gene_A1T CLASS SUBCLASS Name"
    )
    .unwrap();

    let mut batch = empty_batch();
    let err = batch
        .load_mutations(mutation_file.path(), "Escherichia coli")
        .unwrap_err();
    assert!(err.to_string().contains("Duplicate mutations for WP_1"));
}

#[test]
fn susceptible_loader_rejects_duplicate_accessions_like_cpp() {
    let mut susceptible_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        susceptible_file,
        "Escherichia_coli gene WP_2 80 CLASS SUBCLASS Susceptible_Name"
    )
    .unwrap();
    writeln!(
        susceptible_file,
        "Escherichia_coli gene WP_2 70 CLASS SUBCLASS Susceptible_Name"
    )
    .unwrap();

    let mut batch = empty_batch();
    let err = batch
        .load_susceptible(susceptible_file.path(), "Escherichia coli")
        .unwrap_err();
    assert!(err.to_string().contains("Duplicate protein accession WP_2"));
}

#[test]
fn suppress_loader_reads_accession_tokens_without_organism_column() {
    let mut suppress_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(suppress_file, "WP_1\tWP_2").unwrap();
    writeln!(suppress_file, "WP_3").unwrap();

    let mut batch = empty_batch();
    batch
        .load_suppress(suppress_file.path(), "Escherichia")
        .unwrap();

    assert_eq!(batch.suppress_prots, vec!["WP_1", "WP_2", "WP_3"]);
}

#[test]
fn suppress_loader_reads_all_whitespace_tokens_like_cpp() {
    let mut suppress_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(suppress_file, "#taxgroup\tprotein_accession").unwrap();
    writeln!(suppress_file, "Escherichia\tBAE77793.1").unwrap();
    writeln!(suppress_file, "Vibrio_cholerae\tABQ18953.1").unwrap();

    let mut batch = empty_batch();
    batch
        .load_suppress(suppress_file.path(), "Escherichia")
        .unwrap();

    assert_eq!(
        batch.suppress_prots,
        vec!["ABQ18953.1", "BAE77793.1", "Escherichia", "Vibrio_cholerae"]
    );
}

#[test]
fn test_blast_alignment_parsing() {
    let line = "WP_061158039.1|1|1|blaTEM-156|blaTEM|hydrolase|2|BETA-LACTAM|BETA-LACTAM|class_A_beta-lactamase_TEM-156\tblaTEM-156\t1\t286\t287\t1\t286\t286\tMSIQH\tMSIQH";
    let default_br = BlastRule::default();
    let al = BlastAlignment::from_blast_line(line, true, true, &default_br, &default_br);
    assert!(al.is_ok(), "Failed to parse BLAST line: {:?}", al.err());
    let al = al.unwrap();
    assert_eq!(al.hsp.sseqid, "blaTEM-156");
    assert_eq!(al.ref_accession, "WP_061158039.1");
    assert_eq!(al.fam_id, "blaTEM-156");
    assert_eq!(al.gene, "blaTEM");
}

#[test]
fn blast_alignment_normalizes_header_class_and_subclass_like_cpp() {
    let line = "WP_1|1|1|fam|gene|AMR|2|SUB_CLASS|CLASS_NAME|Product\tprot\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*";
    let default_br = BlastRule::default();

    let al = BlastAlignment::from_blast_line(line, true, true, &default_br, &default_br)
        .expect("BLAST line should parse");

    assert_eq!(al.subclass, "SUB CLASS");
    assert_eq!(al.class, "CLASS NAME");
}

#[test]
fn blast_alignment_splits_qseqid_from_right_like_cpp() {
    let line = "WP|EXTRA|1|1|fam|gene|AMR|2|SUB_CLASS|CLASS_NAME|Product\tprot\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*";
    let default_br = BlastRule::default();

    let al = BlastAlignment::from_blast_line(line, true, true, &default_br, &default_br)
        .expect("BLAST line should parse");

    assert_eq!(al.ref_accession, "WP|EXTRA");
    assert_eq!(al.part, 1);
    assert_eq!(al.parts, 1);
    assert_eq!(al.fam_id, "fam");
    assert_eq!(al.gene, "gene");
    assert_eq!(al.product, "Product");
}

#[test]
fn blast_alignment_rejects_malformed_qseqid_like_cpp() {
    let line =
        "WP_1|1|fam|gene|AMR|2|SUB_CLASS|CLASS_NAME|Product\tprot\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*";
    let default_br = BlastRule::default();

    let err = BlastAlignment::from_blast_line(line, true, true, &default_br, &default_br)
        .expect_err("malformed qseqid metadata should fail");

    assert!(err.to_string().contains("Bad AMRFinder database"));
}

#[test]
fn blast_alignment_rejects_nonnumeric_qseqid_fields_like_cpp() {
    let line = "WP_1|part|1|fam|gene|AMR|2|SUB_CLASS|CLASS_NAME|Product\tprot\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*";
    let default_br = BlastRule::default();

    let err = BlastAlignment::from_blast_line(line, true, true, &default_br, &default_br)
        .expect_err("nonnumeric qseqid metadata should fail");

    assert!(err.to_string().contains("Bad AMRFinder database"));
}

#[test]
fn blast_alignment_parses_declarative_mutation_ref_accession() {
    let line = "WP_002851214.1:85:rplV_R85RARAR|1|1|rplV|rplV|mutation|2|||50S_ribosomal_protein_L22\t50S_L22\t1\t5\t6\t1\t5\t5\tRARAR\tRARAR";
    let default_br = BlastRule::default();

    let al = BlastAlignment::from_blast_line(line, true, true, &default_br, &default_br)
        .expect("declarative mutation BLAST line should parse");

    assert_eq!(al.ref_accession, "WP_002851214.1");
    assert_eq!(al.ref_mutation.pos_real, 84);
    assert_eq!(al.ref_mutation.gene_mutation_std, "rplV_R85RARAR");
    assert_eq!(al.ref_mutation.reference, "R");
    assert_eq!(al.ref_mutation.allele, "RARAR");
    assert_eq!(al.ref_mutation.gene, "rplV");
}

#[test]
fn blast_alignment_rejects_bad_declarative_mutation_ref_accession_like_cpp() {
    let default_br = BlastRule::default();
    let missing_pos = "WP_002851214.1:rplV_R85RARAR|1|1|rplV|rplV|mutation|2|||50S_ribosomal_protein_L22\t50S_L22\t1\t5\t6\t1\t5\t5\tRARAR\tRARAR";
    let non_mutation = "WP_002851214.1:85:rplV_R85RARAR|1|1|rplV|rplV|AMR|2|||50S_ribosomal_protein_L22\t50S_L22\t1\t5\t6\t1\t5\t5\tRARAR\tRARAR";

    let err = BlastAlignment::from_blast_line(missing_pos, true, true, &default_br, &default_br)
        .expect_err("missing mutation position should fail");
    assert!(err.to_string().contains("Bad AMRFinder database"));

    let err = BlastAlignment::from_blast_line(non_mutation, true, true, &default_br, &default_br)
        .expect_err("declarative mutation syntax on non-mutation protein should fail");
    assert!(err.to_string().contains("Bad AMRFinder database"));
}

#[test]
fn blast_loader_rejects_missing_cpp_hierarchy_like_get_fam() {
    let mut fam_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        fam_file,
        "known\t\tknown\t-\t0\t0\t70\t0\t90\t70\t0\t50\t2\tAMR\tAMR\tCLASS\tSUBCLASS\tKnown"
    )
    .unwrap();
    let mut blast_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        blast_file,
        "WP_UNKNOWN|1|1|missing_fam|missing_gene|AMR|2|SUBCLASS|CLASS|Product\tprot0\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
    )
    .unwrap();

    let mut batch = Batch::from_fam_file(fam_file.path(), 0).unwrap();
    let err = load_blast_results(blast_file.path(), &mut batch, true, false)
        .expect_err("C++ getFam() rejects missing famId and gene hierarchy");

    assert!(
        err.to_string()
            .contains("Cannot find hierarchy for: missing_fam (genesymbol: missing_gene)"),
        "{err}"
    );
}

#[test]
fn hmm_loader_rejects_unknown_hmm_accession_like_cpp() {
    let mut fam_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        fam_file,
        "fam0\t\tgene0\tHMM_NAME\t100\t50\t90\t0\t90\t90\t0\t50\t2\tAMR\tAMR\tCLASS\tSUBCLASS\tFamily 0"
    )
    .unwrap();
    let mut hmmsearch_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        hmmsearch_file,
        "target - HMM_NAME UNKNOWN_ACC 1e-20 120.0 0.0 1e-20 110.0 0.0 1.0 1 0 0 1 1 1 1 description"
    )
    .unwrap();

    let mut batch = Batch::from_fam_file(fam_file.path(), 0).unwrap();
    let err = load_hmm_results(hmmsearch_file.path(), &mut batch)
        .expect_err("C++ HmmAlignment constructor rejects unknown HMM accession");

    assert!(err.to_string().contains("No family for HMM UNKNOWN_ACC"));
}

#[test]
fn run_amr_report_defaults_to_reportable_entries_only() {
    let mut fam_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        fam_file,
        "fam0\t\tgene0\t-\t0\t0\t90\t0\t90\t90\t0\t50\t0\tAMR\tAMR\tCLASS\tSUBCLASS\tFamily 0"
    )
    .unwrap();

    let mut blast_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        blast_file,
        "WP_0|1|1|fam0|gene0|AMR|0|SUBCLASS|CLASS|Product\tprot0\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
    )
    .unwrap();

    let config = AmrReportConfig {
        fam_file: fam_file.path(),
        blastp_file: Some(blast_file.path()),
        blastx_file: None,
        dna_len_file: None,
        hmmsearch_file: None,
        hmmdom_file: None,
        gff_file: None,
        gff_type: "genbank",
        organism: "",
        mutation_file: None,
        susceptible_file: None,
        suppress_file: None,
        coverage_min: 0.5,
        ident_min: -1.0,
        print_node: false,
        print_node_raw: false,
        mutation_all: None,
        name: "",
        non_reportable: false,
        report_core_only: false,
        report_all_equal: false,
        cds_exist: false,
        nosame: false,
        noblast: false,
        nohmm: false,
        retain_blasts: false,
        skip_hmm_check: false,
    };
    let mut output = Vec::new();

    run_amr_report(&config, &mut output).unwrap();

    let output = String::from_utf8(output).unwrap();
    assert_eq!(output.lines().count(), 1, "{output}");
}

#[test]
fn run_amr_report_writes_mutation_all_after_main_report_like_cpp() {
    struct FailingWriter;

    impl Write for FailingWriter {
        fn write(&mut self, _buf: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::new(
                std::io::ErrorKind::BrokenPipe,
                "forced report failure",
            ))
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    let mut fam_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        fam_file,
        "fam0\t\tgene0\t-\t0\t0\t90\t0\t90\t90\t0\t50\t0\tAMR\tAMR\tCLASS\tSUBCLASS\tFamily 0"
    )
    .unwrap();
    let temp_dir = tempfile::tempdir().unwrap();
    let mutation_all_path = temp_dir.path().join("mutation_all.tsv");

    let config = AmrReportConfig {
        fam_file: fam_file.path(),
        blastp_file: None,
        blastx_file: None,
        dna_len_file: None,
        hmmsearch_file: None,
        hmmdom_file: None,
        gff_file: None,
        gff_type: "genbank",
        organism: "",
        mutation_file: None,
        susceptible_file: None,
        suppress_file: None,
        coverage_min: 0.5,
        ident_min: -1.0,
        print_node: false,
        print_node_raw: false,
        mutation_all: Some(&mutation_all_path),
        name: "",
        non_reportable: false,
        report_core_only: false,
        report_all_equal: false,
        cds_exist: false,
        nosame: false,
        noblast: false,
        nohmm: false,
        retain_blasts: false,
        skip_hmm_check: false,
    };

    let err = run_amr_report(&config, &mut FailingWriter).unwrap_err();
    assert!(
        err.to_string().contains("forced report failure"),
        "unexpected error: {err:#}"
    );
    assert!(
        !mutation_all_path.exists(),
        "C++ writes mutation_all only after the main report succeeds"
    );
}

#[test]
fn run_amr_report_non_reportable_matches_cpp_reportable_min_zero() {
    let mut fam_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        fam_file,
        "fam0\t\tgene0\t-\t0\t0\t90\t0\t90\t90\t0\t50\t0\tAMR\tAMR\tCLASS\tSUBCLASS\tFamily 0"
    )
    .unwrap();

    let mut blast_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        blast_file,
        "WP_0|1|1|fam0|gene0|AMR|0|SUBCLASS|CLASS|Product\tprot0\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
    )
    .unwrap();

    let config = AmrReportConfig {
        fam_file: fam_file.path(),
        blastp_file: Some(blast_file.path()),
        blastx_file: None,
        dna_len_file: None,
        hmmsearch_file: None,
        hmmdom_file: None,
        gff_file: None,
        gff_type: "genbank",
        organism: "",
        mutation_file: None,
        susceptible_file: None,
        suppress_file: None,
        coverage_min: 0.5,
        ident_min: -1.0,
        print_node: false,
        print_node_raw: false,
        mutation_all: None,
        name: "",
        non_reportable: true,
        report_core_only: false,
        report_all_equal: false,
        cds_exist: false,
        nosame: false,
        noblast: false,
        nohmm: false,
        retain_blasts: false,
        skip_hmm_check: false,
    };
    let mut output = Vec::new();

    run_amr_report(&config, &mut output).unwrap();

    let output = String::from_utf8(output).unwrap();
    assert!(output.contains("\nprot0\t"), "{output}");
}

#[test]
fn run_amr_report_noblast_skips_blast_inputs_like_cpp_body() {
    let mut fam_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        fam_file,
        "fam0\t\tgene0\t-\t0\t0\t90\t0\t90\t90\t0\t50\t2\tAMR\tAMR\tCLASS\tSUBCLASS\tFamily 0"
    )
    .unwrap();

    let mut blast_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        blast_file,
        "WP_0|1|1|fam0|gene0|AMR|2|SUBCLASS|CLASS|Product\tprot0\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
    )
    .unwrap();

    let config = AmrReportConfig {
        fam_file: fam_file.path(),
        blastp_file: Some(blast_file.path()),
        blastx_file: None,
        dna_len_file: None,
        hmmsearch_file: None,
        hmmdom_file: None,
        gff_file: None,
        gff_type: "genbank",
        organism: "",
        mutation_file: None,
        susceptible_file: None,
        suppress_file: None,
        coverage_min: 0.5,
        ident_min: -1.0,
        print_node: false,
        print_node_raw: false,
        mutation_all: None,
        name: "",
        non_reportable: false,
        report_core_only: false,
        report_all_equal: false,
        cds_exist: false,
        nosame: false,
        noblast: true,
        nohmm: false,
        retain_blasts: false,
        skip_hmm_check: false,
    };
    let mut output = Vec::new();

    run_amr_report(&config, &mut output).unwrap();

    let output = String::from_utf8(output).unwrap();
    assert_eq!(output.lines().count(), 1, "{output}");
    assert!(!output.contains("\nprot0\t"), "{output}");
}

#[test]
fn run_amr_report_blastx_forces_cds_columns_like_cpp_body() {
    let mut fam_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        fam_file,
        "blaTEM-156\t\tblaTEM\t-\t0\t0\t90\t0\t90\t90\t0\t50\t2\tAMR\tAMR\tBETA-LACTAM\tBETA-LACTAM\tclass A beta-lactamase TEM-156"
    )
    .unwrap();

    let mut blastx_file = tempfile::NamedTempFile::new().unwrap();
    let qseq = format!("M{}", "A".repeat(285));
    writeln!(
        blastx_file,
        "WP_061158039.1|1|1|blaTEM-156|blaTEM|beta-lactamase|2|BETA-LACTAM|BETA-LACTAM|class_A_beta-lactamase_TEM-156\tcontig_blastx\t1\t286\t287\t101\t958\t1200\t{}\t{}",
        qseq, qseq
    )
    .unwrap();

    let config = AmrReportConfig {
        fam_file: fam_file.path(),
        blastp_file: None,
        blastx_file: Some(blastx_file.path()),
        dna_len_file: None,
        hmmsearch_file: None,
        hmmdom_file: None,
        gff_file: None,
        gff_type: "genbank",
        organism: "",
        mutation_file: None,
        susceptible_file: None,
        suppress_file: None,
        coverage_min: 0.5,
        ident_min: -1.0,
        print_node: false,
        print_node_raw: false,
        mutation_all: None,
        name: "",
        non_reportable: false,
        report_core_only: false,
        report_all_equal: false,
        cds_exist: false,
        nosame: false,
        noblast: false,
        nohmm: false,
        retain_blasts: false,
        skip_hmm_check: false,
    };
    let mut output = Vec::new();

    run_amr_report(&config, &mut output).unwrap();

    let output = String::from_utf8(output).unwrap();
    let header = output.lines().next().unwrap();
    assert!(
        header.starts_with("Protein id\tContig id\tStart\tStop\tStrand\t"),
        "{output}"
    );
}

#[test]
fn run_amr_report_validates_thresholds_like_cpp_body() {
    let fam_file = tempfile::NamedTempFile::new().unwrap();
    let make_config = |ident_min, coverage_min| AmrReportConfig {
        fam_file: fam_file.path(),
        blastp_file: None,
        blastx_file: None,
        dna_len_file: None,
        hmmsearch_file: None,
        hmmdom_file: None,
        gff_file: None,
        gff_type: "genbank",
        organism: "",
        mutation_file: None,
        susceptible_file: None,
        suppress_file: None,
        coverage_min,
        ident_min,
        print_node: false,
        print_node_raw: false,
        mutation_all: None,
        name: "",
        non_reportable: false,
        report_core_only: false,
        report_all_equal: false,
        cds_exist: false,
        nosame: false,
        noblast: false,
        nohmm: false,
        retain_blasts: false,
        skip_hmm_check: false,
    };

    let invalid_ident = make_config(-2.0, 0.5);
    let err = run_amr_report(&invalid_ident, &mut Vec::new()).unwrap_err();
    assert_eq!(err.to_string(), "-ident_min must be -1 or between 0 and 1");

    let invalid_coverage = make_config(-1.0, -1.0);
    let err = run_amr_report(&invalid_coverage, &mut Vec::new()).unwrap_err();
    assert_eq!(
        err.to_string(),
        "-coverage_min must be -1 or between 0 and 1"
    );

    let too_high_coverage = make_config(-1.0, 0.95);
    let err = run_amr_report(&too_high_coverage, &mut Vec::new()).unwrap_err();
    assert_eq!(
        err.to_string(),
        "-coverage_min must be less than 0.9 - threshod for complete matches"
    );
}

#[test]
fn run_amr_report_requires_gff_for_combined_blast_inputs_like_cpp_body() {
    let fam_file = tempfile::NamedTempFile::new().unwrap();
    let blastp_file = tempfile::NamedTempFile::new().unwrap();
    let blastx_file = tempfile::NamedTempFile::new().unwrap();
    let config = AmrReportConfig {
        fam_file: fam_file.path(),
        blastp_file: Some(blastp_file.path()),
        blastx_file: Some(blastx_file.path()),
        dna_len_file: None,
        hmmsearch_file: None,
        hmmdom_file: None,
        gff_file: None,
        gff_type: "genbank",
        organism: "",
        mutation_file: None,
        susceptible_file: None,
        suppress_file: None,
        coverage_min: 0.5,
        ident_min: -1.0,
        print_node: false,
        print_node_raw: false,
        mutation_all: None,
        name: "",
        non_reportable: false,
        report_core_only: false,
        report_all_equal: false,
        cds_exist: false,
        nosame: false,
        noblast: false,
        nohmm: false,
        retain_blasts: false,
        skip_hmm_check: false,
    };

    let err = run_amr_report(&config, &mut Vec::new()).unwrap_err();
    assert_eq!(
        err.to_string(),
        "If BLASTP and BLASTX files are present then a GFF file must be present"
    );
}

#[test]
fn run_amr_report_requires_hmmsearch_and_hmmdom_as_pair_like_cpp_body() {
    let fam_file = tempfile::NamedTempFile::new().unwrap();
    let hmmsearch_file = tempfile::NamedTempFile::new().unwrap();
    let config = AmrReportConfig {
        fam_file: fam_file.path(),
        blastp_file: None,
        blastx_file: None,
        dna_len_file: None,
        hmmsearch_file: Some(hmmsearch_file.path()),
        hmmdom_file: None,
        gff_file: None,
        gff_type: "genbank",
        organism: "",
        mutation_file: None,
        susceptible_file: None,
        suppress_file: None,
        coverage_min: 0.5,
        ident_min: -1.0,
        print_node: false,
        print_node_raw: false,
        mutation_all: None,
        name: "",
        non_reportable: false,
        report_core_only: false,
        report_all_equal: false,
        cds_exist: false,
        nosame: false,
        noblast: false,
        nohmm: false,
        retain_blasts: false,
        skip_hmm_check: false,
    };

    let err = run_amr_report(&config, &mut Vec::new()).unwrap_err();
    assert_eq!(
        err.to_string(),
        "hmmsearch and hmmdom must be both present or both absent"
    );
}

#[test]
fn run_amr_report_requires_mutation_table_for_organism_like_cpp_body() {
    let fam_file = tempfile::NamedTempFile::new().unwrap();
    let config = AmrReportConfig {
        fam_file: fam_file.path(),
        blastp_file: None,
        blastx_file: None,
        dna_len_file: None,
        hmmsearch_file: None,
        hmmdom_file: None,
        gff_file: None,
        gff_type: "genbank",
        organism: "Escherichia",
        mutation_file: None,
        susceptible_file: None,
        suppress_file: None,
        coverage_min: 0.5,
        ident_min: -1.0,
        print_node: false,
        print_node_raw: false,
        mutation_all: None,
        name: "",
        non_reportable: false,
        report_core_only: false,
        report_all_equal: false,
        cds_exist: false,
        nosame: false,
        noblast: false,
        nohmm: false,
        retain_blasts: false,
        skip_hmm_check: false,
    };

    let err = run_amr_report(&config, &mut Vec::new()).unwrap_err();
    assert_eq!(err.to_string(), "mutation_tab is empty");
}

#[test]
fn fam_loader_trims_lines_and_rejects_bad_numeric_fields_like_cpp() {
    let mut fam_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(fam_file, "   # skipped after C++ trim").unwrap();
    writeln!(
        fam_file,
        "fam0\t\tgene0\t-\t0\t0\t90\t0\t90\t90\t0\t50\t2\tAMR\tAMR\tCLASS\tSUBCLASS\tFamily 0   "
    )
    .unwrap();

    let batch = Batch::from_fam_file(fam_file.path(), 0).unwrap();
    assert_eq!(batch.fam_map["fam0"].family_name, "Family 0");

    let mut bad_fam_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        bad_fam_file,
        "fam0\t\tgene0\t-\tbad\t0\t90\t0\t90\t90\t0\t50\t2\tAMR\tAMR\tCLASS\tSUBCLASS\tFamily 0"
    )
    .unwrap();

    let err = Batch::from_fam_file(bad_fam_file.path(), 0)
        .err()
        .expect("bad FAM numeric field should fail");
    assert!(
        err.to_string().contains(&format!(
            "Cannot read {}, line 1",
            bad_fam_file.path().display()
        )),
        "{err}"
    );
}

#[test]
fn run_amr_report_nosame_filters_same_reference_blast_hits_like_cpp() {
    let mut fam_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        fam_file,
        "fam0\t\tgene0\t-\t0\t0\t90\t0\t90\t90\t0\t50\t1\tAMR\tAMR\tCLASS\tSUBCLASS\tFamily 0"
    )
    .unwrap();

    let mut blast_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        blast_file,
        "WP_0|1|1|fam0|gene0|AMR|1|SUBCLASS|CLASS|Product\tWP_0\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
    )
    .unwrap();
    writeln!(
        blast_file,
        "WP_0|1|1|fam0|gene0|AMR|1|SUBCLASS|CLASS|Product\tprot_keep\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
    )
    .unwrap();

    let config = AmrReportConfig {
        fam_file: fam_file.path(),
        blastp_file: Some(blast_file.path()),
        blastx_file: None,
        dna_len_file: None,
        hmmsearch_file: None,
        hmmdom_file: None,
        gff_file: None,
        gff_type: "genbank",
        organism: "",
        mutation_file: None,
        susceptible_file: None,
        suppress_file: None,
        coverage_min: 0.5,
        ident_min: -1.0,
        print_node: false,
        print_node_raw: false,
        mutation_all: None,
        name: "",
        non_reportable: false,
        report_core_only: false,
        report_all_equal: false,
        cds_exist: false,
        nosame: true,
        noblast: false,
        nohmm: false,
        retain_blasts: false,
        skip_hmm_check: false,
    };

    let mut output = Vec::new();
    run_amr_report(&config, &mut output).unwrap();
    let output = String::from_utf8(output).unwrap();

    assert!(output.contains("\nprot_keep\t"), "{output}");
    assert!(!output.contains("\nWP_0\t"), "{output}");
}

#[test]
fn run_amr_report_user_identity_overrides_curated_family_identity() {
    let mut fam_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        fam_file,
        "fam0\t\tgene0\t-\t0\t0\t90\t0\t90\t90\t0\t50\t2\tAMR\tAMR\tCLASS\tSUBCLASS\tFamily 0"
    )
    .unwrap();

    let mut blast_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        blast_file,
        "WP_0|1|1|fam0|gene0|AMR|2|SUBCLASS|CLASS|Product\tprot0\t1\t5\t5\t1\t5\t5\tAAAA*\tAATA*"
    )
    .unwrap();

    let config = AmrReportConfig {
        fam_file: fam_file.path(),
        blastp_file: Some(blast_file.path()),
        blastx_file: None,
        dna_len_file: None,
        hmmsearch_file: None,
        hmmdom_file: None,
        gff_file: None,
        gff_type: "genbank",
        organism: "",
        mutation_file: None,
        susceptible_file: None,
        suppress_file: None,
        coverage_min: 0.5,
        ident_min: 0.7,
        print_node: false,
        print_node_raw: false,
        mutation_all: None,
        name: "",
        non_reportable: false,
        report_core_only: false,
        report_all_equal: false,
        cds_exist: false,
        nosame: false,
        noblast: false,
        nohmm: false,
        retain_blasts: false,
        skip_hmm_check: false,
    };

    let mut output = Vec::new();
    run_amr_report(&config, &mut output).unwrap();
    let output = String::from_utf8(output).unwrap();
    assert!(output.contains("\nprot0\t"), "{output}");
}

#[test]
fn run_amr_report_print_node_uses_gene_for_non_exact_allele_like_cpp() {
    let mut fam_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        fam_file,
        "fam\t\tgene\t-\t0\t0\t70\t0\t90\t70\t0\t50\t2\tAMR\tAMR\tCLASS\tSUBCLASS\tFamily"
    )
    .unwrap();
    writeln!(
        fam_file,
        "fam-allele\tfam\tallele_gene\t-\t0\t0\t70\t0\t90\t70\t0\t50\t2\tAMR\tAMR\tCLASS\tSUBCLASS\tAllele Family"
    )
    .unwrap();

    let mut blast_file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        blast_file,
        "WP_0|1|1|fam-allele|fam|AMR|2|SUBCLASS|CLASS|Product\tprot0\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAT*"
    )
    .unwrap();

    let mut config = AmrReportConfig {
        fam_file: fam_file.path(),
        blastp_file: Some(blast_file.path()),
        blastx_file: None,
        dna_len_file: None,
        hmmsearch_file: None,
        hmmdom_file: None,
        gff_file: None,
        gff_type: "genbank",
        organism: "",
        mutation_file: None,
        susceptible_file: None,
        suppress_file: None,
        coverage_min: 0.5,
        ident_min: 0.7,
        print_node: true,
        print_node_raw: false,
        mutation_all: None,
        name: "",
        non_reportable: false,
        report_core_only: false,
        report_all_equal: false,
        cds_exist: false,
        nosame: false,
        noblast: false,
        nohmm: false,
        retain_blasts: false,
        skip_hmm_check: false,
    };
    let mut output = Vec::new();

    run_amr_report(&config, &mut output).unwrap();

    let output = String::from_utf8(output).unwrap();
    let row = output.lines().nth(1).expect("expected one report row");
    let node = row.split('\t').next_back().unwrap();
    assert_eq!(node, "fam", "{output}");

    config.print_node_raw = true;
    let mut raw_output = Vec::new();
    run_amr_report(&config, &mut raw_output).unwrap();

    let raw_output = String::from_utf8(raw_output).unwrap();
    let raw_row = raw_output.lines().nth(1).expect("expected one report row");
    let raw_node = raw_row.split('\t').next_back().unwrap();
    assert_eq!(raw_node, "fam-allele", "{raw_output}");
}
