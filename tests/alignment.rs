use amrfinder::alignment::{Alignment, AmrMutation, SeqChange};
use amrfinder::seq::{Disruption, Hsp, Interval, NO_INDEX};

#[test]
fn amr_mutation_parse_simple() {
    let m = AmrMutation::new(
        100,
        "gyrA_S83L",
        "gyrA_S83L",
        "QUINOLONE",
        "FLUOROQUINOLONE",
        "name",
    );
    assert_eq!(m.pos_real, 99);
    assert_eq!(m.gene, "gyrA");
    assert_eq!(m.reference, "S");
    assert_eq!(m.allele, "L");
    assert_eq!(m.pos_std, 82);
    assert_eq!(m.get_stop(), 100);
    assert_eq!(m.wildtype(), "gyrA_S83S");
}

#[test]
fn amr_mutation_parse_negative_pos() {
    let m = AmrMutation::new(
        40,
        "ampC_T-14TGT",
        "ampC_T-14TGT",
        "BETA-LACTAM",
        "CEPH",
        "name",
    );
    assert_eq!(m.gene, "ampC");
    assert_eq!(m.reference, "T");
    assert_eq!(m.allele, "TGT");
    assert_eq!(m.pos_std, -15);
    assert_eq!(m.pos_real, 39);
}

#[test]
fn amr_mutation_parse_stop_codon() {
    let m = AmrMutation::new(
        141,
        "nfsA_K141Ter",
        "nfsA_K141Ter",
        "NITRO",
        "NITRO",
        "name",
    );
    assert_eq!(m.gene, "nfsA");
    assert_eq!(m.reference, "K");
    assert_eq!(m.allele, "*");
    assert_eq!(m.pos_std, 140);
    assert_eq!(m.pos_real, 140);
}

#[test]
fn amr_mutation_parse_deletion() {
    let m = AmrMutation::new(
        5,
        "pmrB_RPISLR6del",
        "pmrB_RPISLR6del",
        "COLISTIN",
        "COLISTIN",
        "name",
    );
    assert_eq!(m.gene, "pmrB");
    assert_eq!(m.reference, "RPISLR");
    assert_eq!(m.allele, "");
    assert_eq!(m.pos_std, 5);
    assert_eq!(m.frameshift, NO_INDEX);
    assert_eq!(m.frameshift_insertion, 0);
}

#[test]
fn amr_mutation_parse_frameshift_insertion() {
    let m = AmrMutation::new(
        252,
        "cirA_Y253CfsTer5ins1",
        "cirA_Y253CfsTer5ins1",
        "BETA-LACTAM",
        "CEFIDEROCOL",
        "name",
    );
    assert_eq!(m.gene, "cirA");
    assert_eq!(m.reference, "Y");
    assert_eq!(m.allele, "C");
    assert_eq!(m.pos_std, 252);
    assert_eq!(m.frameshift, 5);
    assert_eq!(m.frameshift_insertion, 1);
}

#[test]
fn amr_mutation_parse_frameshift_deletion() {
    let m = AmrMutation::new(
        632,
        "cirA_K633LfsTer8del2",
        "cirA_K633LfsTer8del2",
        "BETA-LACTAM",
        "CEFIDEROCOL",
        "name",
    );
    assert_eq!(m.gene, "cirA");
    assert_eq!(m.reference, "K");
    assert_eq!(m.allele, "L");
    assert_eq!(m.pos_std, 632);
    assert_eq!(m.frameshift, 8);
    assert_eq!(m.frameshift_insertion, -2);
}

#[test]
fn amr_mutation_parse_promoter() {
    let m = AmrMutation::new(162, "blaTEMp_G162T", "blaTEMp_G162T", "BL", "BL", "name");
    assert_eq!(m.gene, "blaTEMp");
    assert_eq!(m.reference, "G");
    assert_eq!(m.allele, "T");
    assert_eq!(m.pos_std, 161);
}

#[test]
fn amr_mutation_empty() {
    let m = AmrMutation::default();
    assert!(m.empty());
}

#[test]
fn amr_mutation_apply_uses_zero_based_internal_pos() {
    let m = AmrMutation::new(3, "gene_C3T", "gene_C3T", "X", "X", "name");
    let mut seq = "AACAA".to_string();
    m.apply(&mut seq).unwrap();
    assert_eq!(seq, "AATAA");
}

#[test]
fn seq_change_get_mutation_str_matches_original_format() {
    let replacement = SeqChange {
        start_ref: 82,
        reference: "S".to_string(),
        allele: "L".to_string(),
        ..Default::default()
    };
    assert_eq!(replacement.get_mutation_str(), "S83L");

    let deletion = SeqChange {
        start_ref: 5,
        reference: "RP".to_string(),
        allele: String::new(),
        ..Default::default()
    };
    assert_eq!(deletion.get_mutation_str(), "RP6DEL");

    let stop = SeqChange {
        start_ref: 140,
        reference: "K".to_string(),
        allele: "*".to_string(),
        ..Default::default()
    };
    assert_eq!(stop.get_mutation_str(), "K141Ter");
}

#[test]
fn seq_change_orders_by_start_target_and_prefers_original_better_rules() {
    let mut left = SeqChange {
        start_target: 10,
        len: 1,
        neighborhood_mismatch: 0.01,
        reference: "A".to_string(),
        allele: "T".to_string(),
        mutations: vec![0],
        ..Default::default()
    };
    let right = SeqChange {
        start_target: 20,
        len: 1,
        neighborhood_mismatch: 0.02,
        reference: "A".to_string(),
        allele: "G".to_string(),
        ..Default::default()
    };

    assert!(left < right);
    assert!(left.better(&right, 0.9, 1.0));
    left.mutations.clear();
    let right = SeqChange {
        mutations: vec![0],
        ..right
    };
    assert!(!left.better(&right, 1.0, 0.9));
}

#[test]
fn disruption_genesymbol_raw_uses_original_interval_encoding() {
    let disr = Disruption {
        prev_hsp_idx: Some(0),
        next_hsp_idx: Some(1),
        prev_start: 0,
        next_stop: 0,
        intron: false,
        q_interval: Interval::new(4, 7, 0),
        s_interval: Interval::new(30, 39, 1),
    };

    let seq_change = SeqChange {
        disr: Some(disr),
        ..Default::default()
    };

    assert_eq!(seq_change.get_mutation_str(), "del_4_7_30_39_1");
}

#[test]
fn seq_change_finish_sets_sequence_and_forward_coordinates() {
    let hsp = Hsp {
        q_int: Interval::new(10, 14, 0),
        s_int: Interval::new(100, 104, 1),
        qlen: 20,
        slen: 200,
        qseq: "ACGT".to_string(),
        sseq: "ATGT".to_string(),
        ..Default::default()
    };
    let mut change = SeqChange {
        start: 1,
        len: 1,
        ..Default::default()
    };

    assert!(change.finish(&hsp, &hsp.qseq, 0));
    assert_eq!(change.reference, "C");
    assert_eq!(change.allele, "T");
    assert_eq!(change.start_ref, 11);
    assert_eq!(change.stop_ref, 12);
    assert_eq!(change.start_target, 101);
}

#[test]
fn seq_change_finish_uses_reverse_strand_target_coordinate() {
    let hsp = Hsp {
        q_int: Interval::new(10, 14, 0),
        s_int: Interval::new(100, 104, -1),
        qlen: 20,
        slen: 200,
        qseq: "ACGT".to_string(),
        sseq: "ATGT".to_string(),
        ..Default::default()
    };
    let mut change = SeqChange {
        start: 1,
        len: 1,
        ..Default::default()
    };

    assert!(change.finish(&hsp, &hsp.qseq, 0));
    assert_eq!(change.start_target, 102);
}

#[test]
fn seq_change_finish_computes_original_neighborhood_mismatch() {
    let hsp = Hsp {
        q_int: Interval::new(1, 6, 0),
        s_int: Interval::new(1, 6, 1),
        qlen: 7,
        slen: 7,
        qseq: "AACGT".to_string(),
        sseq: "ATGGT".to_string(),
        ..Default::default()
    };
    let mut change = SeqChange {
        start: 2,
        len: 1,
        ..Default::default()
    };

    assert!(!change.finish(&hsp, &hsp.qseq, 1));
    assert_eq!(change.reference, "C");
    assert_eq!(change.allele, "G");
    assert!(change.neighborhood_mismatch > 0.04);
}

#[test]
fn seq_change_matches_mutation_checks_reference_and_allele() {
    let change = SeqChange {
        start_ref: 10,
        stop_ref: 11,
        reference: "C".to_string(),
        allele: "T".to_string(),
        len: 1,
        ..Default::default()
    };
    let mutation = AmrMutation::new(11, "gene_C11T", "gene_C11T", "X", "X", "name");
    assert!(change.matches_mutation(&mutation).unwrap());

    let different_allele = AmrMutation::new(11, "gene_C11G", "gene_C11G", "X", "X", "name");
    assert!(!change.matches_mutation(&different_allele).unwrap());

    let wrong_reference = AmrMutation::new(11, "gene_A11T", "gene_A11T", "X", "X", "name");
    let err = change.matches_mutation(&wrong_reference).unwrap_err();
    assert!(err
        .to_string()
        .contains("Reference sequence has 'C', but mutation is: gene_A11T"));
}

#[test]
fn seq_change_save_text_matches_original_field_order() {
    let mutation = AmrMutation::new(11, "gene_C11T", "gene_C11T", "X", "X", "species_name");
    let change = SeqChange {
        start: 2,
        len: 1,
        reference: "C".to_string(),
        allele: "T".to_string(),
        start_ref: 10,
        stop_ref: 11,
        start_target: 100,
        neighborhood_mismatch: 0.25,
        mutations: vec![0],
        ..Default::default()
    };
    let mut out = Vec::new();

    change.save_text(&mut out, &[mutation]).unwrap();

    assert_eq!(
        String::from_utf8(out).unwrap(),
        "3 1 \"C\" -> \"T\" 11..11 101 0.25 11 gene_C11T 0 species name\n"
    );
}

#[test]
fn seq_change_qc_accepts_finished_replacement_change() {
    let hsp = Hsp {
        q_prot: true,
        q_int: Interval::new(10, 14, 0),
        s_int: Interval::new(100, 104, 1),
        qlen: 20,
        slen: 200,
        qseq: "ACGT".to_string(),
        sseq: "ATGT".to_string(),
        ..Default::default()
    };
    let mutation = AmrMutation::new(12, "gene_C12T", "gene_C12T", "X", "X", "name");
    let change = SeqChange {
        start: 1,
        len: 1,
        reference: "C".to_string(),
        allele: "T".to_string(),
        start_ref: 11,
        stop_ref: 12,
        start_target: 101,
        mutations: vec![0],
        ..Default::default()
    };

    change.qc(&hsp, &[mutation]);
}

#[test]
fn alignment_ref_mutation2ref_seq_pos_detects_declarative_mutation() {
    let mutation = AmrMutation::new(12, "gene_T12G", "gene_T12G", "X", "X", "name");
    let alignment = Alignment {
        hsp: Hsp {
            qseqid: "ref".to_string(),
            sseqid: "target".to_string(),
            q_int: Interval::new(10, 14, 0),
            s_int: Interval::new(100, 104, 1),
            qlen: 20,
            slen: 200,
            qseq: "AGGA".to_string(),
            sseq: "AGGA".to_string(),
            ..Default::default()
        },
        ref_mutation: mutation,
        seq_changes: Vec::new(),
    };

    assert_eq!(alignment.ref_mutation2ref_seq_pos(), 1);
}

#[test]
fn alignment_ref_mutation2ref_seq_pos_rejects_changed_flank() {
    let mutation = AmrMutation::new(12, "gene_T12G", "gene_T12G", "X", "X", "name");
    let alignment = Alignment {
        hsp: Hsp {
            q_int: Interval::new(10, 14, 0),
            s_int: Interval::new(100, 104, 1),
            qlen: 20,
            slen: 200,
            qseq: "AGGA".to_string(),
            sseq: "AGTA".to_string(),
            ..Default::default()
        },
        ref_mutation: mutation,
        seq_changes: Vec::new(),
    };

    assert_eq!(alignment.ref_mutation2ref_seq_pos(), NO_INDEX);
}

#[test]
fn alignment_set_seq_changes_records_observed_matching_mutation() {
    let mutations = vec![AmrMutation::new(
        12,
        "gene_C12T",
        "gene_C12T",
        "X",
        "X",
        "name",
    )];
    let mut alignment = Alignment {
        hsp: Hsp::from_blast_line(
            "ref\ttarget\t11\t14\t20\t101\t104\t200\tACGT\tATGT",
            true,
            true,
            true,
            false,
            false,
        )
        .unwrap(),
        ref_mutation: AmrMutation::default(),
        seq_changes: Vec::new(),
    };

    alignment.set_seq_changes(&mutations, 0).unwrap();

    assert_eq!(alignment.seq_changes.len(), 1);
    assert_eq!(alignment.seq_changes[0].reference, "C");
    assert_eq!(alignment.seq_changes[0].allele, "T");
    assert_eq!(alignment.seq_changes[0].mutations, vec![0]);
    assert!(alignment.has_mutation());
    alignment.qc(&mutations);

    let mut out = Vec::new();
    alignment.save_text(&mut out, &mutations).unwrap();
    let text = String::from_utf8(out).unwrap();
    assert!(text.starts_with("Hsp: merged=0 ref(20) 11-14 target(200) 101-104"));
    assert!(text.contains(" #seqChanges:1\n"));
    assert!(text.contains("\"C\" -> \"T\""));
}

#[test]
fn alignment_set_seq_changes_records_wildtype_mutation() {
    let mutations = vec![AmrMutation::new(
        12,
        "gene_C12T",
        "gene_C12T",
        "X",
        "X",
        "name",
    )];
    let mut alignment = Alignment {
        hsp: Hsp {
            q_int: Interval::new(10, 14, 0),
            s_int: Interval::new(100, 104, 1),
            qlen: 20,
            slen: 200,
            qseq: "ACGT".to_string(),
            sseq: "ACGT".to_string(),
            ..Default::default()
        },
        ref_mutation: AmrMutation::default(),
        seq_changes: Vec::new(),
    };

    alignment.set_seq_changes(&mutations, 0).unwrap();

    assert_eq!(alignment.seq_changes.len(), 1);
    assert!(alignment.seq_changes[0].empty());
    assert_eq!(alignment.seq_changes[0].mutations, vec![0]);
    assert!(!alignment.has_mutation());
}

#[test]
fn alignment_set_seq_changes_handles_declarative_mutation() {
    let mutation = AmrMutation::new(12, "gene_T12G", "gene_T12G", "X", "X", "name");
    let mut alignment = Alignment {
        hsp: Hsp {
            q_int: Interval::new(10, 14, 0),
            s_int: Interval::new(100, 104, 1),
            qlen: 20,
            slen: 200,
            qseq: "AGGA".to_string(),
            sseq: "AGGA".to_string(),
            ..Default::default()
        },
        ref_mutation: mutation.clone(),
        seq_changes: Vec::new(),
    };

    alignment.set_seq_changes(&[mutation], 0).unwrap();

    assert_eq!(alignment.seq_changes.len(), 1);
    assert_eq!(alignment.seq_changes[0].start, 1);
    assert_eq!(alignment.seq_changes[0].reference, "T");
    assert_eq!(alignment.seq_changes[0].allele, "G");
    assert_eq!(alignment.seq_changes[0].mutations, vec![0]);
}
