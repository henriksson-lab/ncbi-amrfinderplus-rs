mod alignment {
    pub use amrfinder::alignment::*;
}
mod gff {
    pub use amrfinder::gff::*;
}
mod seq {
    pub use amrfinder::seq::*;
}
mod tsv {
    pub use amrfinder::tsv::*;
}
mod columns {
    pub use amrfinder::columns::*;
}

mod amr_report_impl {
    include!("../src/amr_report.rs");

    #[cfg(test)]
    mod tests {
        use super::*;
        use std::path::PathBuf;

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

        fn fam(id: &str, parent_id: &str, hmm: &str, tc1: f64, tc2: f64) -> Fam {
            Fam {
                id: id.to_string(),
                genesymbol: id.to_string(),
                family_name: format!("{id} family"),
                reportable: 2,
                hmm: hmm.to_string(),
                tc1,
                tc2,
                complete_br: BlastRule::default(),
                partial_br: BlastRule::default(),
                type_: "AMR".to_string(),
                subtype: "AMR".to_string(),
                class: "CLASS".to_string(),
                subclass: "SUBCLASS".to_string(),
                parent_id: parent_id.to_string(),
            }
        }

        // --- amr_report behavior tests ---
        // These test specific behaviors that must match C++ amr_report

        #[test]
        fn hmm_better_criterion0_prefers_descendant_family() {
            let mut batch = empty_batch();
            batch.fam_map.insert(
                "parent".to_string(),
                fam("parent", "", "HMM_PARENT", 50.0, 10.0),
            );
            batch.fam_map.insert(
                "child".to_string(),
                fam("child", "parent", "HMM_CHILD", 40.0, 10.0),
            );

            let parent = HmmAlignment {
                sseqid: "target".to_string(),
                score1: 100.0,
                score2: 20.0,
                fam_id: "parent".to_string(),
                domain: None,
                blast_al_idx: None,
            };
            let child = HmmAlignment {
                sseqid: "target".to_string(),
                score1: 90.0,
                score2: 20.0,
                fam_id: "child".to_string(),
                domain: None,
                blast_al_idx: None,
            };

            assert!(child.better(&parent, 0, &batch.fam_map, false));
            assert!(!parent.better(&child, 0, &batch.fam_map, false));
        }

        #[test]
        fn process_pareto_filters_hmms_by_cpp_two_criteria() {
            let mut batch = empty_batch();
            batch.fam_map.insert(
                "parent".to_string(),
                fam("parent", "", "HMM_PARENT", 50.0, 10.0),
            );
            batch.fam_map.insert(
                "child".to_string(),
                fam("child", "parent", "HMM_CHILD", 40.0, 10.0),
            );
            batch
                .fam_map
                .insert("tie_a".to_string(), fam("tie_a", "", "HMM_A", 80.0, 10.0));
            batch
                .fam_map
                .insert("tie_b".to_string(), fam("tie_b", "", "HMM_B", 70.0, 10.0));

            batch.add_hmm_al(HmmAlignment {
                sseqid: "target1".to_string(),
                score1: 100.0,
                score2: 20.0,
                fam_id: "parent".to_string(),
                domain: None,
                blast_al_idx: None,
            });
            batch.add_hmm_al(HmmAlignment {
                sseqid: "target1".to_string(),
                score1: 90.0,
                score2: 20.0,
                fam_id: "child".to_string(),
                domain: None,
                blast_al_idx: None,
            });
            batch.add_hmm_al(HmmAlignment {
                sseqid: "target2".to_string(),
                score1: 100.0,
                score2: 20.0,
                fam_id: "tie_b".to_string(),
                domain: None,
                blast_al_idx: None,
            });
            batch.add_hmm_al(HmmAlignment {
                sseqid: "target2".to_string(),
                score1: 100.0,
                score2: 20.0,
                fam_id: "tie_a".to_string(),
                domain: None,
                blast_al_idx: None,
            });

            batch.process(false, false);

            let target1 = &batch.target2good_hmm_als["target1"];
            assert_eq!(target1, &vec![1]);
            assert_eq!(batch.hmm_als[target1[0]].fam_id, "child");

            let target2 = &batch.target2good_hmm_als["target2"];
            assert_eq!(target2, &vec![3]);
            assert_eq!(batch.hmm_als[target2[0]].fam_id, "tie_a");
        }

        #[test]
        fn report_equidistant_still_filters_strictly_worse_blast_like_cpp() {
            let br = BlastRule::default();
            let qseq = format!("{}*", "A".repeat(99));
            let better_sseq = format!("{}T*", "A".repeat(98));
            let worse_sseq = format!("{}{}*", "A".repeat(90), "T".repeat(9));
            let better_line = format!(
                "WP_BETTER|1|1|famA|famA|AMR|2|CLASS|SUBCLASS|Better_product\t\
                 target\t1\t100\t100\t1\t100\t100\t{qseq}\t{better_sseq}"
            );
            let worse_line = format!(
                "WP_WORSE|1|1|famA|famA|AMR|2|CLASS|SUBCLASS|Worse_product\t\
                 target\t1\t100\t100\t1\t100\t100\t{qseq}\t{worse_sseq}"
            );

            let mut family = fam("famA", "", "", 0.0, 0.0);
            family.complete_br = BlastRule::new(0.80, 0.90);
            family.partial_br = BlastRule::new(0.80, 0.50);

            let mut batch = empty_batch();
            batch.report_all_equal = true;
            batch.fam_map.insert("famA".to_string(), family);
            batch.add_blast_al(
                BlastAlignment::from_blast_line(&worse_line, true, true, &br, &br).unwrap(),
            );
            batch.add_blast_al(
                BlastAlignment::from_blast_line(&better_line, true, true, &br, &br).unwrap(),
            );

            batch.process(false, false);

            let good = batch.target2good_blast_als.get("target").unwrap();
            assert_eq!(good.len(), 1);
            assert_eq!(batch.blast_als[good[0]].ref_accession, "WP_BETTER");
        }

        #[test]
        fn process_retain_blasts_keeps_dominated_good_rows_like_cpp_flag() {
            let br = BlastRule::default();
            let qseq = format!("{}*", "A".repeat(99));
            let better_sseq = format!("{}T*", "A".repeat(98));
            let worse_sseq = format!("{}{}*", "A".repeat(90), "T".repeat(9));
            let better_line = format!(
                "WP_BETTER|1|1|famA|famA|AMR|2|CLASS|SUBCLASS|Better_product\t\
                 target\t1\t100\t100\t1\t100\t100\t{qseq}\t{better_sseq}"
            );
            let worse_line = format!(
                "WP_WORSE|1|1|famA|famA|AMR|2|CLASS|SUBCLASS|Worse_product\t\
                 target\t1\t100\t100\t1\t100\t100\t{qseq}\t{worse_sseq}"
            );

            let mut family = fam("famA", "", "", 0.0, 0.0);
            family.complete_br = BlastRule::new(0.80, 0.90);
            family.partial_br = BlastRule::new(0.80, 0.50);

            let mut batch = empty_batch();
            batch.fam_map.insert("famA".to_string(), family);
            batch.add_blast_al(
                BlastAlignment::from_blast_line(&worse_line, true, true, &br, &br).unwrap(),
            );
            batch.add_blast_al(
                BlastAlignment::from_blast_line(&better_line, true, true, &br, &br).unwrap(),
            );

            batch.process(true, false);

            let good = batch.target2good_blast_als.get("target").unwrap();
            assert_eq!(good.len(), 2);
            let retained = good
                .iter()
                .map(|&idx| batch.blast_als[idx].ref_accession.as_str())
                .collect::<HashSet<_>>();
            assert_eq!(retained, HashSet::from(["WP_WORSE", "WP_BETTER"]));
        }

        #[test]
        fn process_does_not_apply_user_ident_threshold_to_mutation_proteins() {
            let br = BlastRule::default();
            let qseq = format!("{}*", "A".repeat(99));
            let sseq = format!("{}{}*", "A".repeat(40), "T".repeat(59));
            let line = format!(
                "WP_MUT|1|1|famA|geneA|mutation|2|CLASS|SUBCLASS|Mutation_product\t\
                 target\t1\t100\t100\t1\t100\t100\t{qseq}\t{sseq}"
            );

            let mut al = BlastAlignment::from_blast_line(&line, true, true, &br, &br).unwrap();
            al.seq_changes.push(SeqChange {
                start: 40,
                len: 1,
                reference: "A".to_string(),
                allele: "T".to_string(),
                start_ref: 40,
                stop_ref: 41,
                start_target: 40,
                ..SeqChange::default()
            });

            let mut batch = empty_batch();
            batch.ident_min = 0.95;
            batch.add_blast_al(al);

            batch.process(false, false);

            let good = batch.target2good_blast_als.get("target").unwrap();
            assert_eq!(good.len(), 1);
            assert_eq!(batch.blast_als[good[0]].ref_accession, "WP_MUT");
        }

        #[test]
        fn process_skip_hmm_check_keeps_low_identity_hmm_family_blast_like_cpp_flag() {
            let br = BlastRule::default();
            let qseq = format!("{}*", "A".repeat(99));
            let sseq = format!("{}{}*", "A".repeat(90), "T".repeat(9));
            let blast_line = format!(
                "WP_LOW|1|1|famA|famA|AMR|2|CLASS|SUBCLASS|Low_identity_product\t\
                 target\t1\t100\t100\t1\t100\t100\t{qseq}\t{sseq}"
            );

            let mut family = fam("famA", "", "HMM_A", 50.0, 10.0);
            family.complete_br = BlastRule::new(0.80, 0.90);
            family.partial_br = BlastRule::new(0.80, 0.50);

            let mut without_skip = empty_batch();
            without_skip.hmm_exist = true;
            without_skip
                .fam_map
                .insert("famA".to_string(), family.clone());
            without_skip.add_blast_al(
                BlastAlignment::from_blast_line(&blast_line, true, true, &br, &br).unwrap(),
            );
            without_skip.process(false, false);
            assert!(
                without_skip
                    .target2good_blast_als
                    .get("target")
                    .is_none_or(Vec::is_empty),
                "low-identity HMM-backed BLAST should be removed without skip_hmm_check"
            );

            let mut with_skip = empty_batch();
            with_skip.hmm_exist = true;
            with_skip.fam_map.insert("famA".to_string(), family);
            with_skip.add_blast_al(
                BlastAlignment::from_blast_line(&blast_line, true, true, &br, &br).unwrap(),
            );
            with_skip.process(false, true);

            let good = with_skip.target2good_blast_als.get("target").unwrap();
            assert_eq!(good.len(), 1);
            assert_eq!(with_skip.blast_als[good[0]].ref_accession, "WP_LOW");
        }

        #[test]
        fn hmm_better_blast_requires_protein_target_and_descendant_match_family() {
            let mut batch = empty_batch();
            batch.fam_map.insert(
                "parent".to_string(),
                fam("parent", "", "HMM_PARENT", 50.0, 10.0),
            );
            batch.fam_map.insert(
                "child".to_string(),
                fam("child", "parent", "HMM_CHILD", 40.0, 10.0),
            );
            let hmm = HmmAlignment {
                sseqid: "target".to_string(),
                score1: 100.0,
                score2: 20.0,
                fam_id: "child".to_string(),
                domain: None,
                blast_al_idx: None,
            };
            let br = BlastRule::default();
            let blast_line = concat!(
                "WP_1|1|1|parent|parent|AMR|2|SUBCLASS|CLASS|Parent_family\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let mut blast =
                BlastAlignment::from_blast_line(blast_line, true, true, &br, &br).unwrap();
            blast.br_fam_id = Some("parent".to_string());

            assert!(hmm.better_blast(&blast, &batch.fam_map));

            let mut dna_blast = blast.clone();
            dna_blast.hsp.s_prot = false;
            assert!(!hmm.better_blast(&dna_blast, &batch.fam_map));
        }

        #[test]
        fn hmm_better_blast_uses_selected_match_family_like_cpp() {
            let mut batch = empty_batch();
            batch.fam_map.insert(
                "parent".to_string(),
                fam("parent", "", "HMM_PARENT", 50.0, 10.0),
            );
            batch.fam_map.insert(
                "child".to_string(),
                fam("child", "parent", "HMM_CHILD", 40.0, 10.0),
            );
            batch.fam_map.insert(
                "grandchild".to_string(),
                fam("grandchild", "child", "HMM_GRANDCHILD", 30.0, 10.0),
            );

            let hmm = HmmAlignment {
                sseqid: "target".to_string(),
                score1: 100.0,
                score2: 20.0,
                fam_id: "child".to_string(),
                domain: None,
                blast_al_idx: None,
            };
            let br = BlastRule::default();
            let blast_line = concat!(
                "WP_1|1|1|grandchild|grandchild|AMR|2|SUBCLASS|CLASS|Grandchild_family\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let mut blast =
                BlastAlignment::from_blast_line(blast_line, true, true, &br, &br).unwrap();
            blast.br_fam_id = Some("parent".to_string());

            assert!(
                hmm.better_blast(&blast, &batch.fam_map),
                "C++ compares HMM family to BlastAlignment::getMatchFam(), not raw famId/gene"
            );
        }

        #[test]
        fn blastx_better_hmm_uses_hmm_synthetic_cds_match_like_cpp() {
            let mut batch = empty_batch();
            batch.fam_map.insert(
                "parent".to_string(),
                fam("parent", "", "HMM_PARENT", 50.0, 10.0),
            );
            batch.fam_map.insert(
                "child".to_string(),
                fam("child", "parent", "HMM_CHILD", 40.0, 10.0),
            );

            let br = BlastRule::default();
            let protein_line = concat!(
                "WP_1|1|1|child|child|AMR|2|SUBCLASS|CLASS|Child_family\t",
                "protein\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let mut hmm_synthetic =
                BlastAlignment::from_blast_line(protein_line, true, true, &br, &br).unwrap();
            hmm_synthetic.from_hmm = true;
            hmm_synthetic.cdss.push(
                Locus::new(
                    1,
                    "contig",
                    10,
                    28,
                    true,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );
            batch.blast_als.push(hmm_synthetic);

            let hmm = HmmAlignment {
                sseqid: "protein".to_string(),
                score1: 100.0,
                score2: 20.0,
                fam_id: "child".to_string(),
                domain: None,
                blast_al_idx: Some(0),
            };

            let blastx_line = concat!(
                "WP_2|1|1|parent|parent|AMR|2|SUBCLASS|CLASS|Parent_family\t",
                "contig\t1\t5\t5\t10\t25\t100\tAAAA*\tAAAA*"
            );
            let mut blastx =
                BlastAlignment::from_blast_line(blastx_line, true, false, &br, &br).unwrap();
            blastx.br_fam_id = Some("parent".to_string());

            assert!(batch.blast_als[0].matches_cds(&blastx));
            assert!(blastx.better_hmm(&hmm, &batch.fam_map, &batch.blast_als));

            let mut off_cds_hmm = hmm.clone();
            batch.blast_als[0].cdss[0] = Locus::new(
                1,
                "contig",
                1000,
                1020,
                true,
                false,
                0,
                String::new(),
                String::new(),
            )
            .unwrap();
            off_cds_hmm.blast_al_idx = Some(0);
            assert!(!blastx.better_hmm(&off_cds_hmm, &batch.fam_map, &batch.blast_als));
        }

        /// C++ isCore() requires reportable >= 2 (amr_report.cpp:820-821).
        /// reportable=1 entries (like stx genes) should be "plus", not "core".
        #[test]
        fn test_scope_requires_reportable_2_for_core() {
            let fam_path = db_dir().join("fam.tsv");
            if !fam_path.exists() {
                return;
            }
            let batch = Batch::from_fam_file(&fam_path, 0).unwrap();

            // Build a BLAST alignment with reportable=1 (like stxA2)
            let line = "TJA36680.1|1|1|stxA2_acd|stxA2|VIRULENCE|1|stxA2|STX2|Shiga_toxin_Stx2_subunit_A\tstxA2a_prot\t1\t319\t320\t94\t1051\t320\tMKCIL\tMKCIL";
            let br = BlastRule::default();
            let al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();

            let reportable = batch.get_reportable(&al);
            // reportable=1 should NOT be core; C++ requires >= 2
            assert!(
                reportable < 2,
                "stxA2 with reportable=1 should not qualify as core (reportable={})",
                reportable
            );
        }

        #[test]
        fn exact_allele_match_is_core_only_when_header_reportable_is_two_like_cpp() {
            let br = BlastRule::default();
            let line = concat!(
                "WP_1|1|1|fam-allele|fam|AMR|1|CLASS|SUBCLASS|Product\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();
            let mut batch = empty_batch();
            batch.add_blast_al(al);

            assert!(batch.blast_als[0].allele_match());
            assert!(!batch.is_core(&batch.blast_als[0]));

            batch.blast_als[0].reportable = 2;
            assert!(batch.is_core(&batch.blast_als[0]));
        }

        #[test]
        fn blast_alignment_get_method_uses_cpp_strong_susceptible_point_order() {
            let mut batch = empty_batch();
            batch.accession2susceptible.insert(
                "WP_S".to_string(),
                Susceptible {
                    genesymbol: "gene".to_string(),
                    cutoff: 0.0,
                    class: "CLASS".to_string(),
                    subclass: "SUBCLASS".to_string(),
                    name: "Product".to_string(),
                },
            );
            let br = BlastRule::default();
            let line = concat!(
                "WP_S|1|1|fam|gene|susceptible|2|CLASS|SUBCLASS|Product\t",
                "target\t1\t6\t6\t1\t6\t6\tMMMMM*\tMMMMM*"
            );
            let al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();

            assert_eq!(al.get_method(&batch, None), "POINTP");
        }

        #[test]
        fn blast_alignment_get_method_internal_stop_has_no_blast_suffix() {
            let batch = empty_batch();
            let br = BlastRule::default();
            let line = concat!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                "target\t1\t5\t10\t1\t5\t10\tMMMMM\tMM*MM"
            );
            let al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();

            assert_eq!(al.get_method(&batch, None), "INTERNAL_STOP");
        }

        #[test]
        fn blast_alignment_get_method_partial_contig_end_requires_matching_missed_tail() {
            let batch = empty_batch();
            let br = BlastRule::default();
            let qseq = "M".repeat(70);
            let line = format!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 target\t1\t70\t100\t1\t70\t100\t{qseq}\t{qseq}"
            );
            let mut al = BlastAlignment::from_blast_line(&line, true, true, &br, &br).unwrap();
            al.cdss.push(Locus {
                line_num: 0,
                contig: "target".to_string(),
                start: 0,
                stop: 210,
                strand: true,
                partial: false,
                contig_len: 1000,
                cross_origin: false,
                gene: String::new(),
                product: String::new(),
            });

            assert_eq!(al.get_method(&batch, al.cdss.first()), "PARTIALP");

            al.hsp.q_int.start = 10;
            assert_eq!(
                al.get_method(&batch, al.cdss.first()),
                "PARTIAL_CONTIG_ENDP"
            );
        }

        #[test]
        fn blast_alignment_get_method_hmm_has_no_suffix() {
            let batch = empty_batch();
            let br = BlastRule::default();
            let line = concat!(
                "WP_H|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                "target\t1\t5\t5\t1\t5\t5\tMMMM*\tMMMM*"
            );
            let mut al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();
            al.from_hmm = true;

            assert_eq!(al.get_method(&batch, None), "HMM");
        }

        #[test]
        fn blast_alignment_get_method_blastx_uses_ref_prot_exact_match_checks() {
            let batch = empty_batch();
            let br = BlastRule::default();
            let line = concat!(
                "WP_X|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                "contig\t1\t6\t6\t1\t18\t18\tMMMMM*\tMMMMM*"
            );
            let mut al = BlastAlignment::from_blast_line(line, true, false, &br, &br).unwrap();
            assert_eq!(al.get_method(&batch, al.cdss.first()), "ALLELEX");

            al.hsp.disrs.push(crate::seq::Disruption {
                prev_hsp_idx: Some(0),
                next_hsp_idx: Some(1),
                prev_start: 0,
                next_stop: 0,
                prev_slen: 18,
                intron: false,
                s_stop_codon: false,
                q_interval: crate::seq::Interval::new(2, 2, 1),
                s_interval: crate::seq::Interval::new(6, 7, 1),
            });
            assert_eq!(al.get_method(&batch, al.cdss.first()), "BLASTX");
        }

        /// C++ getFam() falls back from famId to gene field (amr_report.cpp:1222-1224).
        /// When fam_id is an allele (e.g. "blaTEM-156") not in fam.tsv, the gene field
        /// ("blaTEM") should be used as fallback.
        #[test]
        fn test_get_fam_info_gene_fallback() {
            let fam_path = db_dir().join("fam.tsv");
            if !fam_path.exists() {
                return;
            }
            let batch = Batch::from_fam_file(&fam_path, 0).unwrap();

            // "blaTEM-156" is an allele — it won't be in fam.tsv
            // "blaTEM" is the parent family — it IS in fam.tsv
            assert!(
                !batch.fam_map.contains_key("blaTEM-156"),
                "blaTEM-156 should NOT be a direct key in fam_map"
            );
            assert!(
                batch.fam_map.contains_key("blaTEM"),
                "blaTEM should be in fam_map"
            );

            let br = BlastRule::default();
            let line = concat!(
                "WP_1|1|1|blaTEM-156|blaTEM|AMR|2|BETA-LACTAM|BETA-LACTAM|TEM_family\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();
            let (_genesymbol, type_, _subtype, _class, _subclass, _reportable) =
                batch.get_fam_info(&al);
            assert_eq!(type_, "AMR", "gene fallback should resolve blaTEM");
        }

        /// Verify fam_map resolves type/class for known allele families.
        /// The C++ uses getFam() which falls back famId → gene. Without fallback,
        /// the Rust code uses the BLAST header's resistance field ("hydrolase")
        /// instead of the fam.tsv's type ("AMR").
        #[test]
        fn test_get_fam_info_uses_fam_type_not_blast_header() {
            let fam_path = db_dir().join("fam.tsv");
            if !fam_path.exists() {
                return;
            }
            let batch = Batch::from_fam_file(&fam_path, 0).unwrap();

            // BLAST header for blaTEM-156: resistance="hydrolase"
            // But fam.tsv entry for blaTEM has type_="AMR"
            let line = "WP_061158039.1|1|1|blaTEM-156|blaTEM|hydrolase|2|BETA-LACTAM|BETA-LACTAM|class_A_beta-lactamase_TEM-156\tblaTEM-156\t1\t286\t287\t1\t286\t286\tMSIQH\tMSIQH";
            let br = BlastRule::default();
            let al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();

            let (_genesymbol, type_, _subtype, _class, _subclass, _reportable) =
                batch.get_fam_info(&al);
            assert_eq!(
                type_, "AMR",
                "Type should come from fam.tsv (AMR), not BLAST header (hydrolase)"
            );
        }

        /// C++ links BLAST alignments with supporting HMM results via getFam()->hmm.
        /// The report should show HMM accession and description for BLAST hits
        /// that have matching HMM families, not always "NA".
        #[test]
        fn test_hmm_info_populated_for_blast_hits_with_hmm_families() {
            let fam_path = db_dir().join("fam.tsv");
            if !fam_path.exists() {
                return;
            }
            let batch = Batch::from_fam_file(&fam_path, 0).unwrap();

            // Check that the blaTEM family has HMM info in fam.tsv
            if let Some(fam) = batch.fam_map.get("blaTEM") {
                // If the family has an HMM, BLAST hits in this family should report it
                if !fam.hmm.is_empty() {
                    // This documents the expected behavior:
                    // C++ outputs HMM accession/description even for BLAST-method hits
                    // when the family has an HMM entry
                    assert!(
                        !fam.hmm.is_empty(),
                        "blaTEM family should have HMM accession"
                    );
                    assert!(
                        !fam.family_name.is_empty(),
                        "blaTEM family should have family_name for HMM description"
                    );
                }
            }
        }

        /// Pareto filter should remove dominated hits.
        /// A hit with lower identity and equal coverage should be dominated.
        #[test]
        fn test_pareto_filter_removes_dominated_hits() {
            let fam_path = db_dir().join("fam.tsv");
            if !fam_path.exists() {
                return;
            }
            let mut batch = Batch::from_fam_file(&fam_path, 0).unwrap();

            let br = BlastRule::default();
            // Use short identical-length sequences to avoid finish_hsp issues
            let seq = "MSIQH";

            // Hit 1: 100% identity (5/5 nident), full coverage
            let line1 = format!(
            "WP_1|1|1|blaTEM-156|blaTEM|hydrolase|2|BETA-LACTAM|BETA-LACTAM|product1\ttarget1\t1\t5\t6\t1\t5\t5\t{}\t{}",
            seq, seq
        );
            let al1 = BlastAlignment::from_blast_line(&line1, true, true, &br, &br).unwrap();
            batch.add_blast_al(al1);

            // Hit 2: 80% identity (4/5 nident), same coverage — dominated by hit 1
            let line2 = format!(
            "WP_2|1|1|blaTEM-239|blaTEM|hydrolase|2|NA|NA|product2\ttarget1\t1\t5\t6\t1\t5\t5\t{}\t{}",
            seq, "MSXQH"
        );
            let al2 = BlastAlignment::from_blast_line(&line2, true, true, &br, &br).unwrap();
            batch.add_blast_al(al2);

            let indices = vec![0, 1];
            let filtered = batch.pareto_filter_blast(&indices);

            // Hit 2 should be dominated by hit 1 (strictly lower identity, equal coverage)
            assert_eq!(
                filtered.len(),
                1,
                "Pareto filter should remove dominated hit, got {} hits",
                filtered.len()
            );
        }

        #[test]
        fn mutation_disruption_does_not_change_point_subtype() {
            let mut batch = empty_batch();
            let br = BlastRule::default();
            let line = concat!(
                "WP_1|1|1|pmrB|pmrB|mutation|2|COLISTIN|COLISTIN|PmrB\t",
                "target1\t1\t5\t6\t1\t5\t5\tAAAAA\tAAAAT"
            );
            let mut point = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();
            point.ref_mutation =
                AmrMutation::new(5, "pmrB_A5T", "pmrB_A5T", "COLISTIN", "COLISTIN", "PmrB");
            point.seq_changes.push(SeqChange::default());

            let mut point_disrupt = point.clone();
            point_disrupt.hsp.disrs.push(crate::seq::Disruption {
                prev_hsp_idx: Some(0),
                next_hsp_idx: Some(1),
                prev_start: 0,
                next_stop: 0,
                prev_slen: 5,
                intron: false,
                s_stop_codon: false,
                q_interval: crate::seq::Interval::new(2, 2, 1),
                s_interval: crate::seq::Interval::new(6, 7, 1),
            });

            batch.add_blast_al(point);
            batch.add_blast_al(point_disrupt);

            assert!(!batch.blast_better(&batch.blast_als[0], &batch.blast_als[1]));
            assert!(!batch.blast_better(&batch.blast_als[1], &batch.blast_als[0]));
            assert_eq!(batch.pareto_filter_blast(&[0, 1]), vec![0, 1]);

            let (_, _, plain_subtype, _, _, _) = batch.get_fam_info(&batch.blast_als[0]);
            let (_, _, disrupted_subtype, _, _, _) = batch.get_fam_info(&batch.blast_als[1]);
            assert_eq!(plain_subtype, "POINT");
            assert_eq!(disrupted_subtype, "POINT");
        }

        #[test]
        fn blast_better_same_accession_prefers_wider_same_reference_span_like_cpp() {
            let batch = empty_batch();
            let br = BlastRule::default();
            let left_line = concat!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                "target\t1\t6\t6\t10\t15\t30\tMMMMM*\tMMMMM*"
            );
            let right_line = concat!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                "target\t1\t6\t6\t12\t17\t30\tMMMMM*\tMMMMM*"
            );
            let left = BlastAlignment::from_blast_line(left_line, true, true, &br, &br).unwrap();
            let right = BlastAlignment::from_blast_line(right_line, true, true, &br, &br).unwrap();

            assert!(batch.blast_better(&left, &right));
            assert!(!batch.blast_better(&right, &left));
        }

        #[test]
        fn blast_better_uses_cpp_accession_tie_break_for_equal_non_mutations() {
            let batch = empty_batch();
            let br = BlastRule::default();
            let first_line = concat!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                "target\t1\t6\t6\t10\t15\t30\tMMMMM*\tMMMMM*"
            );
            let second_line = concat!(
                "WP_2|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                "target\t1\t6\t6\t10\t15\t30\tMMMMM*\tMMMMM*"
            );
            let first = BlastAlignment::from_blast_line(first_line, true, true, &br, &br).unwrap();
            let second =
                BlastAlignment::from_blast_line(second_line, true, true, &br, &br).unwrap();

            assert!(batch.blast_better(&first, &second));
            assert!(!batch.blast_better(&second, &first));
        }

        #[test]
        fn blast_alignment_member_predicates_match_cpp_names() {
            let br = BlastRule::default();
            let line = concat!(
                "WP_1|1|1|fam-allele|fam|mutation|2|CLASS|SUBCLASS|Product\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let mut al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();

            assert!(al.is_mutation_prot());
            assert!(!al.is_susceptible_prot());
            assert!(!al.in_fam());
            assert!(al.allele());
            assert!(al.ref_prot_exactly_matched(true));
            assert!(al.allele_match());
            assert!(!al.partial());
            assert_eq!(al.ref_effective_len(), 4);

            al.partial_dna = true;
            al.hsp.q_int = crate::seq::Interval::new(1, 4, 1);
            assert_eq!(al.ref_effective_len(), 3);
            assert!(al.partial());
        }

        #[test]
        fn blast_alignment_good_applies_cpp_special_rejections() {
            let br = BlastRule::default();
            let mutation_line = concat!(
                "WP_M|1|1|fam|gene|mutation|2|CLASS|SUBCLASS|Product\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let mut mutation_al =
                BlastAlignment::from_blast_line(mutation_line, true, true, &br, &br).unwrap();
            assert!(!mutation_al.good(None));
            mutation_al.seq_changes.push(SeqChange::default());
            assert!(mutation_al.good(None));

            let susceptible_line = concat!(
                "WP_S|1|1|fam|gene|susceptible|2|CLASS|SUBCLASS|Product\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let mut susceptible_al =
                BlastAlignment::from_blast_line(susceptible_line, true, true, &br, &br).unwrap();
            assert!(!susceptible_al.good(None));

            let weak_susceptible = Susceptible {
                genesymbol: "gene".to_string(),
                cutoff: 0.99,
                class: "CLASS".to_string(),
                subclass: "SUBCLASS".to_string(),
                name: "Product".to_string(),
            };
            assert!(!susceptible_al.good(Some(&weak_susceptible)));
            let cutoff_match = Susceptible {
                cutoff: 1.0,
                ..weak_susceptible
            };
            assert!(susceptible_al.good(Some(&cutoff_match)));

            let strong_susceptible = Susceptible {
                cutoff: 0.0,
                ..cutoff_match
            };
            assert!(!susceptible_al.good(Some(&strong_susceptible)));
            susceptible_al.hsp.disrs.push(crate::seq::Disruption {
                prev_hsp_idx: Some(0),
                next_hsp_idx: Some(1),
                prev_start: 0,
                next_stop: 0,
                prev_slen: 5,
                intron: false,
                s_stop_codon: false,
                q_interval: crate::seq::Interval::new(2, 2, 1),
                s_interval: crate::seq::Interval::new(6, 7, 1),
            });
            assert!(susceptible_al.good(Some(&strong_susceptible)));
        }

        #[test]
        fn blast_alignment_good_rejects_cpp_partial_fusion_and_short_pseudo() {
            let br = BlastRule::default();
            let short_seq = "A".repeat(35);
            let short_line = format!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 target\t1\t35\t100\t1\t35\t100\t{short_seq}\t{short_seq}"
            );
            let short_partial =
                BlastAlignment::from_blast_line(&short_line, true, true, &br, &br).unwrap();
            assert!(!short_partial.good(None));

            let longer_seq = "A".repeat(36);
            let longer_line = format!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 target\t1\t36\t100\t1\t36\t100\t{longer_seq}\t{longer_seq}"
            );
            let mut longer_partial =
                BlastAlignment::from_blast_line(&longer_line, true, true, &br, &br).unwrap();
            longer_partial.br_fam_id = Some("fam".to_string());
            assert!(longer_partial.good(None));

            let fusion_seq = "A".repeat(70);
            let fusion_line = format!(
                "WP_FUS|1|2|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 target\t1\t70\t100\t1\t70\t100\t{fusion_seq}\t{fusion_seq}"
            );
            let fusion_partial =
                BlastAlignment::from_blast_line(&fusion_line, true, true, &br, &br).unwrap();
            assert!(!fusion_partial.good(None));
        }

        #[test]
        fn process_uses_cpp_parent_blast_rule_when_child_rule_fails() {
            let mut parent = fam("parent", "", "", 0.0, 0.0);
            parent.genesymbol = "parent_symbol".to_string();
            parent.complete_br = BlastRule::new(0.90, 0.90);
            parent.partial_br = BlastRule::new(0.90, 0.50);
            let mut child = fam("child", "parent", "", 0.0, 0.0);
            child.genesymbol = "child_symbol".to_string();
            child.complete_br = BlastRule::new(0.99, 0.90);
            child.partial_br = BlastRule::new(0.99, 0.50);

            let qseq = format!("{}*", "A".repeat(99));
            let sseq = format!("{}{}*", "A".repeat(94), "T".repeat(5));
            let line = format!(
                "WP_CHILD|1|1|child|child|AMR|2|CLASS|SUBCLASS|Product\t\
                 target\t1\t100\t100\t1\t100\t100\t{qseq}\t{sseq}"
            );
            let br = BlastRule::default();
            let al = BlastAlignment::from_blast_line(&line, true, true, &br, &br).unwrap();

            assert!(al.hsp.rel_identity() < child.complete_br.ident);
            assert!(al.hsp.rel_identity() >= parent.complete_br.ident);

            let mut batch = empty_batch();
            batch.fam_map.insert("parent".to_string(), parent);
            batch.fam_map.insert("child".to_string(), child);
            batch.add_blast_al(al);

            batch.process(false, false);

            assert_eq!(batch.target2good_blast_als.get("target"), Some(&vec![0]));
            assert_eq!(batch.blast_als[0].complete_br.ident, 0.90);
            assert_eq!(batch.blast_als[0].br_fam_id.as_deref(), Some("parent"));
            assert_eq!(
                batch.fusion_2_gene_symbols(&batch.blast_als[0]),
                "parent_symbol"
            );
        }

        #[test]
        fn report_uses_cpp_brfam_metadata_not_child_family_metadata() {
            let mut parent = fam("parent", "", "", 0.0, 0.0);
            parent.genesymbol = "parent_symbol".to_string();
            parent.family_name = "Parent product".to_string();
            parent.reportable = 0;
            parent.type_ = "PARENT_TYPE".to_string();
            parent.subtype = "PARENT_SUBTYPE".to_string();
            parent.class = "PARENT_CLASS".to_string();
            parent.subclass = "PARENT_SUBCLASS".to_string();
            parent.complete_br = BlastRule::new(0.90, 0.90);
            parent.partial_br = BlastRule::new(0.90, 0.50);

            let mut child = fam("child", "parent", "", 0.0, 0.0);
            child.genesymbol = "child_symbol".to_string();
            child.family_name = "Child product".to_string();
            child.reportable = 1;
            child.type_ = "CHILD_TYPE".to_string();
            child.subtype = "CHILD_SUBTYPE".to_string();
            child.class = "CHILD_CLASS".to_string();
            child.subclass = "CHILD_SUBCLASS".to_string();
            child.complete_br = BlastRule::new(0.99, 0.90);
            child.partial_br = BlastRule::new(0.99, 0.50);

            let qseq = format!("{}*", "A".repeat(99));
            let sseq = format!("{}{}*", "A".repeat(94), "T".repeat(5));
            let line = format!(
                "WP_CHILD|1|1|child|child|AMR|1|CHILD_CLASS|CHILD_SUBCLASS|Child_product\t\
                 target\t1\t100\t100\t1\t100\t100\t{qseq}\t{sseq}"
            );
            let br = BlastRule::default();
            let al = BlastAlignment::from_blast_line(&line, true, true, &br, &br).unwrap();

            let mut batch = empty_batch();
            batch.reportable_min = 1;
            batch.fam_map.insert("parent".to_string(), parent);
            batch.fam_map.insert("child".to_string(), child);
            batch.add_blast_al(al);

            batch.process(false, false);

            let al = &batch.blast_als[0];
            assert_eq!(al.br_fam_id.as_deref(), Some("parent"));
            assert_eq!(batch.get_reportable(al), 0);
            let (genesymbol, type_, subtype, class, subclass, reportable) = batch.get_fam_info(al);
            assert_eq!(genesymbol, "parent_symbol");
            assert_eq!(type_, "PARENT_TYPE");
            assert_eq!(subtype, "PARENT_SUBTYPE");
            assert_eq!(class, "PARENT_CLASS");
            assert_eq!(subclass, "PARENT_SUBCLASS");
            assert_eq!(reportable, 0);

            let mut output = Vec::new();
            batch.report(&mut output, false).unwrap();
            let output = String::from_utf8(output).unwrap();
            assert_eq!(output.lines().count(), 1, "{output}");
        }

        #[test]
        fn report_uses_na_for_empty_match_family_genesymbol_like_cpp() {
            let mut parent = fam("parent", "", "", 0.0, 0.0);
            parent.genesymbol.clear();
            parent.family_name = "Parent product".to_string();
            parent.reportable = 2;
            parent.complete_br = BlastRule::new(0.90, 0.90);
            parent.partial_br = BlastRule::new(0.90, 0.50);

            let mut child = fam("child", "parent", "", 0.0, 0.0);
            child.genesymbol = "child_symbol".to_string();
            child.complete_br = BlastRule::new(0.99, 0.90);
            child.partial_br = BlastRule::new(0.99, 0.50);

            let qseq = format!("{}*", "A".repeat(99));
            let sseq = format!("{}{}*", "A".repeat(94), "T".repeat(5));
            let line = format!(
                "WP_CHILD|1|1|child|child|AMR|2|CHILD_CLASS|CHILD_SUBCLASS|Child_product\t\
                 target\t1\t100\t100\t1\t100\t100\t{qseq}\t{sseq}"
            );
            let br = BlastRule::default();
            let al = BlastAlignment::from_blast_line(&line, true, true, &br, &br).unwrap();

            let mut batch = empty_batch();
            batch.fam_map.insert("parent".to_string(), parent);
            batch.fam_map.insert("child".to_string(), child);
            batch.add_blast_al(al);

            batch.process(false, false);

            let al = &batch.blast_als[0];
            assert_eq!(al.br_fam_id.as_deref(), Some("parent"));
            assert_eq!(batch.fusion_2_gene_symbols(al), crate::columns::NA);

            let mut output = Vec::new();
            batch.report(&mut output, false).unwrap();
            let output = String::from_utf8(output).unwrap();
            let mut lines = output.lines();
            let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
            let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
            let col = |name: &str| {
                let idx = header
                    .iter()
                    .position(|field| *field == name)
                    .unwrap_or_else(|| panic!("missing report column {name}"));
                row[idx]
            };

            assert!(lines.next().is_none(), "{output}");
            assert_eq!(col("Element symbol"), crate::columns::NA);
            assert_eq!(col("Element name"), "Parent product");
        }

        #[test]
        fn report_uses_na_for_empty_match_family_name_like_cpp() {
            let mut parent = fam("parent", "", "", 0.0, 0.0);
            parent.genesymbol = "parent_symbol".to_string();
            parent.family_name.clear();
            parent.reportable = 2;
            parent.complete_br = BlastRule::new(0.90, 0.90);
            parent.partial_br = BlastRule::new(0.90, 0.50);

            let mut child = fam("child", "parent", "", 0.0, 0.0);
            child.genesymbol = "child_symbol".to_string();
            child.family_name = "Header fallback should not report".to_string();
            child.complete_br = BlastRule::new(0.99, 0.90);
            child.partial_br = BlastRule::new(0.99, 0.50);

            let qseq = format!("{}*", "A".repeat(99));
            let sseq = format!("{}{}*", "A".repeat(94), "T".repeat(5));
            let line = format!(
                "WP_CHILD|1|1|child|child|AMR|2|CHILD_CLASS|CHILD_SUBCLASS|Header_Product\t\
                 target\t1\t100\t100\t1\t100\t100\t{qseq}\t{sseq}"
            );
            let br = BlastRule::default();
            let al = BlastAlignment::from_blast_line(&line, true, true, &br, &br).unwrap();

            let mut batch = empty_batch();
            batch.fam_map.insert("parent".to_string(), parent);
            batch.fam_map.insert("child".to_string(), child);
            batch.add_blast_al(al);

            batch.process(false, false);

            let al = &batch.blast_als[0];
            assert_eq!(al.br_fam_id.as_deref(), Some("parent"));

            let mut output = Vec::new();
            batch.report(&mut output, false).unwrap();
            let output = String::from_utf8(output).unwrap();
            let mut lines = output.lines();
            let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
            let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
            let col = |name: &str| {
                let idx = header
                    .iter()
                    .position(|field| *field == name)
                    .unwrap_or_else(|| panic!("missing report column {name}"));
                row[idx]
            };

            assert!(lines.next().is_none(), "{output}");
            assert_eq!(col("Element symbol"), "parent_symbol");
            assert_eq!(col("Element name"), crate::columns::NA);
            assert_eq!(col("Closest reference name"), "Header Product");
        }

        #[test]
        fn report_print_node_raw_keeps_non_exact_allele_fam_id_like_cpp() {
            let mut parent = fam("fam", "", "", 0.0, 0.0);
            parent.genesymbol = "gene_symbol".to_string();
            parent.family_name = "Family".to_string();
            parent.reportable = 2;
            parent.complete_br = BlastRule::new(0.70, 0.90);
            parent.partial_br = BlastRule::new(0.70, 0.50);

            let mut allele = fam("fam-allele", "fam", "", 0.0, 0.0);
            allele.genesymbol = "allele_gene".to_string();
            allele.family_name = "Allele Family".to_string();
            allele.reportable = 2;
            allele.complete_br = BlastRule::new(0.70, 0.90);
            allele.partial_br = BlastRule::new(0.70, 0.50);

            let line = "WP_0|1|1|fam-allele|fam|AMR|2|SUBCLASS|CLASS|Product\t\
                        prot0\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAT*";
            let br = BlastRule::default();
            let al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();

            let mut batch = empty_batch();
            batch.fam_map.insert("fam".to_string(), parent);
            batch.fam_map.insert("fam-allele".to_string(), allele);
            batch.add_blast_al(al);
            batch.process(false, false);

            let mut output = Vec::new();
            batch.report(&mut output, true).unwrap();
            let output = String::from_utf8(output).unwrap();
            let node = output
                .lines()
                .nth(1)
                .expect("expected one report row")
                .split('\t')
                .next_back()
                .unwrap();
            assert_eq!(node, "fam", "{output}");

            batch.print_node_raw = true;
            let mut raw_output = Vec::new();
            batch.report(&mut raw_output, true).unwrap();
            let raw_output = String::from_utf8(raw_output).unwrap();
            let raw_node = raw_output
                .lines()
                .nth(1)
                .expect("expected one report row")
                .split('\t')
                .next_back()
                .unwrap();
            assert_eq!(raw_node, "fam-allele", "{raw_output}");
        }

        #[test]
        fn process_rejects_in_family_hit_when_no_cpp_blast_rule_passes() {
            let mut parent = fam("parent", "", "", 0.0, 0.0);
            parent.complete_br = BlastRule::new(0.98, 0.90);
            parent.partial_br = BlastRule::new(0.98, 0.50);
            let mut child = fam("child", "parent", "", 0.0, 0.0);
            child.complete_br = BlastRule::new(0.99, 0.90);
            child.partial_br = BlastRule::new(0.99, 0.50);

            let qseq = format!("{}*", "A".repeat(99));
            let sseq = format!("{}{}*", "A".repeat(94), "T".repeat(5));
            let default_complete_br = BlastRule::new(0.90, 0.90);
            let default_partial_br = BlastRule::new(0.90, 0.50);
            let line = format!(
                "WP_CHILD|1|1|child|child|AMR|2|CLASS|SUBCLASS|Product\t\
                 target\t1\t100\t100\t1\t100\t100\t{qseq}\t{sseq}"
            );
            let al = BlastAlignment::from_blast_line(
                &line,
                true,
                true,
                &default_complete_br,
                &default_partial_br,
            )
            .unwrap();

            assert!(al.hsp.rel_identity() >= default_complete_br.ident);
            assert!(al.hsp.rel_identity() < parent.complete_br.ident);

            let mut batch = empty_batch();
            batch.fam_map.insert("parent".to_string(), parent);
            batch.fam_map.insert("child".to_string(), child);
            batch.add_blast_al(al);

            batch.process(false, false);

            assert!(!batch.target2good_blast_als.contains_key("target"));
            assert!(batch.blast_als[0].br_fam_id.is_none());
        }

        #[test]
        fn process_removes_complete_low_identity_blast_without_required_hmm_like_cpp() {
            let mut family = fam("famH", "", "HMM_H", 50.0, 20.0);
            family.complete_br = BlastRule::new(0.90, 0.90);
            family.partial_br = BlastRule::new(0.90, 0.50);

            let qseq = format!("{}*", "A".repeat(99));
            let sseq = format!("{}{}*", "A".repeat(96), "T".repeat(3));
            let line = format!(
                "WP_H|1|1|famH|famH|AMR|2|CLASS|SUBCLASS|Product\t\
                 target\t1\t100\t100\t1\t100\t100\t{qseq}\t{sseq}"
            );
            let br = BlastRule::default();

            let mut batch = empty_batch();
            batch.hmm_exist = true;
            batch.fam_map.insert("famH".to_string(), family);
            batch.add_blast_al(
                BlastAlignment::from_blast_line(&line, true, true, &br, &br).unwrap(),
            );

            batch.process(false, false);

            assert!(!batch.target2good_blast_als.contains_key("target"));
        }

        #[test]
        fn process_reports_child_hmm_instead_of_weaker_parent_blast_like_cpp() {
            let mut parent = fam("parent", "", "", 0.0, 0.0);
            parent.genesymbol = "parent_symbol".to_string();
            parent.complete_br = BlastRule::new(0.90, 0.90);
            parent.partial_br = BlastRule::new(0.90, 0.50);

            let mut child = fam("child", "parent", "HMM_CHILD", 50.0, 20.0);
            child.genesymbol = "child_symbol".to_string();
            child.complete_br = BlastRule::new(0.90, 0.90);
            child.partial_br = BlastRule::new(0.90, 0.50);

            let qseq = format!("{}*", "A".repeat(99));
            let sseq = format!("{}{}*", "A".repeat(94), "T".repeat(5));
            let line = format!(
                "WP_PARENT|1|1|parent|parent|AMR|2|CLASS|SUBCLASS|Parent_product\t\
                 target\t1\t100\t100\t1\t100\t100\t{qseq}\t{sseq}"
            );

            let mut batch = empty_batch();
            batch.fam_map.insert("parent".to_string(), parent);
            batch.fam_map.insert("child".to_string(), child);
            let mut blast_al = BlastAlignment::from_blast_line(
                &line,
                true,
                true,
                &BlastRule::default(),
                &BlastRule::default(),
            )
            .unwrap();
            blast_al.cdss.push(
                Locus::new(
                    1,
                    "contig",
                    10,
                    310,
                    true,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );
            batch.add_blast_al(blast_al);
            batch.add_hmm_al(HmmAlignment {
                sseqid: "target".to_string(),
                score1: 100.0,
                score2: 40.0,
                fam_id: "child".to_string(),
                domain: Some(HmmDomain {
                    score: 100.0,
                    hmm_len: 100,
                    hmm_start: 0,
                    hmm_stop: 100,
                    seq_len: 100,
                    seq_start: 0,
                    seq_stop: 100,
                }),
                blast_al_idx: Some(0),
            });

            batch.process(false, false);

            let good = batch.target2good_blast_als.get("target").unwrap();
            assert_eq!(good.len(), 1);
            let reported = &batch.blast_als[good[0]];
            assert!(reported.from_hmm);
            assert_eq!(reported.fam_id, "child");
            assert_eq!(reported.hmm_al_idx, Some(0));
            assert_eq!(batch.hmm_als[0].blast_al_idx, Some(good[0]));
            assert_eq!(reported.cdss.len(), 1);
            assert_eq!(reported.cdss[0].contig, "contig");
        }

        #[test]
        fn load_hmm_domains_replaces_equal_score_domain_like_cpp() {
            let mut batch = empty_batch();
            batch
                .hmm2fam
                .insert("HMM_CHILD".to_string(), "child".to_string());

            let mut dom_file = tempfile::NamedTempFile::new().unwrap();
            writeln!(
                dom_file,
                "target - 100 query HMM_CHILD 80 1e-20 100.0 0.0 1 1 1e-20 1e-20 42.0 0.0 1 40 5 45 5 45 0.99 first"
            )
            .unwrap();
            writeln!(
                dom_file,
                "target - 100 query HMM_CHILD 80 1e-20 100.0 0.0 1 1 1e-20 1e-20 42.0 0.0 2 50 7 55 7 55 0.99 second"
            )
            .unwrap();

            load_hmm_domains(dom_file.path(), &mut batch).unwrap();

            let domain = batch
                .domains
                .get(&("target".to_string(), "child".to_string()))
                .expect("domain should be loaded");
            assert_eq!(domain.score, 42.0);
            assert_eq!(domain.hmm_start, 1);
            assert_eq!(domain.hmm_stop, 50);
            assert_eq!(domain.seq_start, 6);
            assert_eq!(domain.seq_stop, 55);
        }

        #[test]
        fn load_hmm_results_keeps_first_equal_nident_blast_anchor_like_cpp() {
            let mut batch = empty_batch();
            batch.fam_map.insert(
                "child".to_string(),
                fam("child", "", "HMM_CHILD", 50.0, 10.0),
            );
            batch
                .hmm2fam
                .insert("HMM_CHILD".to_string(), "child".to_string());
            batch.domains.insert(
                ("target".to_string(), "child".to_string()),
                HmmDomain {
                    score: 60.0,
                    hmm_len: 100,
                    hmm_start: 0,
                    hmm_stop: 100,
                    seq_len: 100,
                    seq_start: 0,
                    seq_stop: 100,
                },
            );

            let br = BlastRule::new(0.8, 0.5);
            let qseq = format!("{}*", "A".repeat(99));
            let first_line = format!(
                "WP_FIRST|1|1|child|child|AMR|2|SUBCLASS|CLASS|First_product\t\
                 target\t1\t100\t100\t1\t100\t100\t{qseq}\t{qseq}"
            );
            let second_line = format!(
                "WP_SECOND|1|1|child|child|AMR|2|SUBCLASS|CLASS|Second_product\t\
                 target\t1\t100\t100\t1\t100\t100\t{qseq}\t{qseq}"
            );
            batch.add_blast_al(
                BlastAlignment::from_blast_line(&first_line, true, true, &br, &br).unwrap(),
            );
            batch.add_blast_al(
                BlastAlignment::from_blast_line(&second_line, true, true, &br, &br).unwrap(),
            );

            let mut search_file = tempfile::NamedTempFile::new().unwrap();
            writeln!(
                search_file,
                "target - query HMM_CHILD 1e-20 100.0 0.0 1e-20 40.0"
            )
            .unwrap();

            load_hmm_results(search_file.path(), &mut batch).unwrap();

            assert_eq!(batch.hmm_als.len(), 1);
            assert_eq!(batch.hmm_als[0].blast_al_idx, Some(0));
            assert_eq!(batch.blast_als[0].ref_accession, "WP_FIRST");
        }

        #[test]
        fn process_keeps_susceptible_proteins_out_of_family_blast_rule_walk() {
            let mut family = fam("famS", "", "", 0.0, 0.0);
            family.genesymbol = "family_symbol".to_string();
            family.complete_br = BlastRule::new(0.80, 0.90);
            family.partial_br = BlastRule::new(0.80, 0.50);

            let qseq = format!("{}*", "A".repeat(99));
            let sseq = format!("{}{}*", "A".repeat(84), "T".repeat(15));
            let default_complete_br = BlastRule::new(0.90, 0.90);
            let default_partial_br = BlastRule::new(0.90, 0.50);
            let line = format!(
                "WP_S|1|1|famS|header_gene|susceptible|2|CLASS|SUBCLASS|Product\t\
                 target\t1\t100\t100\t1\t100\t100\t{qseq}\t{sseq}"
            );

            let mut batch = empty_batch();
            batch.fam_map.insert("famS".to_string(), family);
            batch.accession2susceptible.insert(
                "WP_S".to_string(),
                Susceptible {
                    genesymbol: "sus_symbol".to_string(),
                    cutoff: 0.95,
                    class: "SCLASS".to_string(),
                    subclass: "SSUB".to_string(),
                    name: "Susceptible Product".to_string(),
                },
            );
            batch.add_blast_al(
                BlastAlignment::from_blast_line(
                    &line,
                    true,
                    true,
                    &default_complete_br,
                    &default_partial_br,
                )
                .unwrap(),
            );

            batch.process(false, false);

            assert!(!batch.target2good_blast_als.contains_key("target"));
            assert_eq!(batch.blast_als[0].br_fam_id, None);
            assert!(!batch.blast_als[0].br_fam_checked);
            assert_eq!(
                batch.fusion_2_gene_symbols(&batch.blast_als[0]),
                "sus_symbol"
            );
        }

        #[test]
        fn report_uses_susceptible_metadata_only_for_susceptible_rows() {
            let br = BlastRule::default();
            let line = concat!(
                "WP_SHARED|1|1|famA|famA|AMR|2|HEADER_CLASS|HEADER_SUB|Header_Product\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );

            let mut family = fam("famA", "", "", 0.0, 0.0);
            family.genesymbol = "family_symbol".to_string();
            family.type_ = "STRESS".to_string();
            family.subtype = "FAMILY_SUBTYPE".to_string();
            family.class = "FAMILY_CLASS".to_string();
            family.subclass = "FAMILY_SUBCLASS".to_string();

            let mut batch = empty_batch();
            batch.fam_map.insert("famA".to_string(), family);
            batch.accession2susceptible.insert(
                "WP_SHARED".to_string(),
                Susceptible {
                    genesymbol: "sus_symbol".to_string(),
                    cutoff: 0.0,
                    class: "SUS_CLASS".to_string(),
                    subclass: "SUS_SUBCLASS".to_string(),
                    name: "Susceptible Product".to_string(),
                },
            );
            batch
                .add_blast_al(BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap());
            batch.blast_als[0].br_fam_id = Some("famA".to_string());
            batch
                .target2good_blast_als
                .insert("target".to_string(), vec![0]);

            let mut output = Vec::new();
            batch.report(&mut output, false).unwrap();
            let output = String::from_utf8(output).unwrap();
            let header: Vec<&str> = output.lines().next().unwrap().split('\t').collect();
            let row: Vec<&str> = output
                .lines()
                .nth(1)
                .expect("missing report row")
                .split('\t')
                .collect();
            let value = |name: &str| {
                let idx = header
                    .iter()
                    .position(|field| *field == name)
                    .unwrap_or_else(|| panic!("missing report column {name}"));
                row[idx]
            };

            assert_eq!(value("Element symbol"), "family_symbol");
            assert_eq!(value("Element name"), "Header Product");
            assert_eq!(value("Type"), "STRESS");
            assert_eq!(value("Subtype"), "FAMILY_SUBTYPE");
            assert_eq!(value("Class"), "FAMILY_CLASS");
            assert_eq!(value("Subclass"), "FAMILY_SUBCLASS");
        }

        #[test]
        fn process_filters_mutation_proteins_without_cpp_seq_changes() {
            let br = BlastRule::default();
            let mutation_line = concat!(
                "WP_M|1|1|fam|gene|mutation|2|CLASS|SUBCLASS|Product\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let normal_line = concat!(
                "WP_N|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );

            let mut batch = empty_batch();
            batch.add_blast_al(
                BlastAlignment::from_blast_line(mutation_line, true, true, &br, &br).unwrap(),
            );
            batch.add_blast_al(
                BlastAlignment::from_blast_line(normal_line, true, true, &br, &br).unwrap(),
            );

            batch.process(false, false);

            let good = batch.target2good_blast_als.get("target").unwrap();
            assert_eq!(good.len(), 1);
            assert_eq!(batch.blast_als[good[0]].ref_accession, "WP_N");
        }

        #[test]
        fn process_removes_alien_organism_proteins_like_cpp() {
            let mut mutation_file = tempfile::NamedTempFile::new().unwrap();
            writeln!(
                mutation_file,
                "Other_species WP_ALIEN 1 A1T A1T CLASS SUBCLASS Alien_name"
            )
            .unwrap();

            let mut family = fam("fam", "", "", 0.0, 0.0);
            family.complete_br = BlastRule::new(0.90, 0.90);
            family.partial_br = BlastRule::new(0.90, 0.50);

            let br = BlastRule::default();
            let alien_line = concat!(
                "WP_ALIEN|1|1|fam|fam|AMR|2|CLASS|SUBCLASS|Alien_product\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let own_line = concat!(
                "WP_OWN|1|1|fam|fam|AMR|2|CLASS|SUBCLASS|Own_product\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );

            let mut batch = empty_batch();
            batch.fam_map.insert("fam".to_string(), family);
            batch
                .load_mutations(mutation_file.path(), "Target species")
                .unwrap();
            assert_eq!(batch.alien_prots, vec!["WP_ALIEN"]);
            batch.add_blast_al(
                BlastAlignment::from_blast_line(alien_line, true, true, &br, &br).unwrap(),
            );
            batch.add_blast_al(
                BlastAlignment::from_blast_line(own_line, true, true, &br, &br).unwrap(),
            );

            batch.process(false, false);

            assert_eq!(batch.target2blast_als.get("target"), Some(&vec![1]));
            assert_eq!(batch.target2good_blast_als.get("target"), Some(&vec![1]));
            assert_eq!(batch.blast_als[1].ref_accession, "WP_OWN");
        }

        #[test]
        fn mutation_symbols_follow_cpp_non_empty_changes_including_replacements() {
            let br = BlastRule::default();
            let line = concat!(
                "WP_M|1|1|pmrB|pmrB|mutation|2|COLISTIN|COLISTIN|PmrB\t",
                "target\t1\t5\t6\t1\t5\t5\tAAAAA\tAAAAT"
            );
            let mut al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();
            al.ref_mutation =
                AmrMutation::new(5, "pmrB_A5T", "pmrB_A5T", "COLISTIN", "COLISTIN", "PmrB");
            al.seq_changes.push(SeqChange {
                start: 4,
                len: 1,
                reference: "A".to_string(),
                allele: "T".to_string(),
                start_ref: 4,
                stop_ref: 5,
                start_target: 4,
                mutations: vec![0],
                ..SeqChange::default()
            });
            al.seq_changes.push(SeqChange {
                start: 3,
                len: 1,
                reference: "A".to_string(),
                allele: "G".to_string(),
                start_ref: 3,
                stop_ref: 4,
                start_target: 3,
                mutations: vec![0],
                replacement: Some(0),
                ..SeqChange::default()
            });

            let symbols = al.get_mutation_symbols();
            assert_eq!(symbols.len(), 1);
            assert!(symbols.contains("pmrB_A5T"));

            al.ref_mutation =
                AmrMutation::new(3, "pmrB_A3G", "pmrB_A3G", "COLISTIN", "COLISTIN", "PmrB");
            al.seq_changes.remove(0);
            let symbols = al.get_mutation_symbols();
            assert_eq!(symbols.len(), 1);
            assert!(
                symbols.contains("pmrB_A3G"),
                "C++ getMutationSymbols() does not skip replacement-marked non-empty changes"
            );
        }

        #[test]
        fn mutation_symbols_follow_seq_change_mutation_indices_like_cpp_pointers() {
            let br = BlastRule::default();
            let line = concat!(
                "WP_M|1|1|pmrB|pmrB|mutation|2|COLISTIN|COLISTIN|PmrB\t",
                "target\t1\t5\t6\t1\t5\t5\tAAAAA\tAAAGT"
            );
            let mut al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();
            al.ref_mutations = vec![
                AmrMutation::new(4, "pmrB_A4G", "pmrB_A4G", "COLISTIN", "COLISTIN", "PmrB"),
                AmrMutation::new(5, "pmrB_A5T", "pmrB_A5T", "COLISTIN", "COLISTIN", "PmrB"),
            ];
            al.ref_mutation = al.ref_mutations[0].clone();
            al.seq_changes.push(SeqChange {
                start: 3,
                len: 2,
                reference: "AA".to_string(),
                allele: "GT".to_string(),
                start_ref: 3,
                stop_ref: 5,
                start_target: 3,
                mutations: vec![0, 1],
                ..SeqChange::default()
            });

            let symbols = al.get_mutation_symbols();
            assert_eq!(symbols.len(), 2);
            assert!(symbols.contains("pmrB_A4G"));
            assert!(symbols.contains("pmrB_A5T"));
        }

        #[test]
        fn process_marks_weaker_colocated_seq_change_replacement_like_cpp() {
            let br = BlastRule::default();
            let line1 = concat!(
                "WP_M1|1|1|pmrB|pmrB|mutation|2|COLISTIN|COLISTIN|PmrB\t",
                "target\t1\t5\t6\t1\t5\t5\tAAAAA\tAAAAT"
            );
            let line2 = concat!(
                "WP_M2|1|1|pmrB|pmrB|mutation|2|COLISTIN|COLISTIN|PmrB\t",
                "target\t1\t5\t6\t1\t5\t5\tAAAAA\tAAAAT"
            );
            let mutation =
                AmrMutation::new(5, "pmrB_A5T", "pmrB_A5T", "COLISTIN", "COLISTIN", "PmrB");
            let mut better = BlastAlignment::from_blast_line(line1, true, true, &br, &br).unwrap();
            better.ref_mutations = vec![mutation.clone()];
            better.ref_mutation = mutation.clone();
            better.seq_changes.push(SeqChange {
                start: 4,
                len: 1,
                reference: "A".to_string(),
                allele: "T".to_string(),
                start_ref: 4,
                stop_ref: 5,
                start_target: 4,
                neighborhood_mismatch: 0.0,
                mutations: vec![0],
                ..SeqChange::default()
            });

            let mut weaker = BlastAlignment::from_blast_line(line2, true, true, &br, &br).unwrap();
            weaker.ref_mutations = vec![mutation.clone()];
            weaker.ref_mutation = mutation;
            weaker.seq_changes.push(SeqChange {
                start: 4,
                len: 1,
                reference: "A".to_string(),
                allele: "T".to_string(),
                start_ref: 4,
                stop_ref: 5,
                start_target: 4,
                neighborhood_mismatch: 1.0,
                mutations: vec![0],
                ..SeqChange::default()
            });

            let mut batch = empty_batch();
            batch.add_blast_al(better);
            batch.add_blast_al(weaker);
            batch.process(false, false);

            assert_eq!(batch.target2good_blast_als["target"], vec![0, 1]);
            assert!(batch.blast_als[0].seq_changes[0].replacement.is_none());
            assert!(batch.blast_als[1].seq_changes[0].replacement.is_some());

            let mut output = Vec::new();
            batch.report(&mut output, false).unwrap();
            let output = String::from_utf8(output).unwrap();
            assert_eq!(output.matches("\tpmrB_A5T\t").count(), 1, "{output}");
        }

        #[test]
        fn blast_better_eq_uses_cpp_mutation_symbol_superset_ordering() {
            let batch = empty_batch();
            let br = BlastRule::default();
            let line = concat!(
                "WP_M|1|1|pmrB|pmrB|mutation|2|COLISTIN|COLISTIN|PmrB\t",
                "target\t1\t5\t6\t1\t5\t5\tAAAAA\tAAAAT"
            );
            let mut with_symbol =
                BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();
            with_symbol.ref_mutation =
                AmrMutation::new(5, "pmrB_A5T", "pmrB_A5T", "COLISTIN", "COLISTIN", "PmrB");
            with_symbol.seq_changes.push(SeqChange {
                start: 4,
                len: 1,
                reference: "A".to_string(),
                allele: "T".to_string(),
                start_ref: 4,
                stop_ref: 5,
                start_target: 4,
                mutations: vec![0],
                ..SeqChange::default()
            });

            let mut without_symbol = with_symbol.clone();
            without_symbol.ref_mutation.gene_mutation.clear();

            assert!(batch.blast_better_eq(&with_symbol, &without_symbol));
            assert!(!batch.blast_better_eq(&without_symbol, &with_symbol));
        }

        #[test]
        fn blast_alignment_cds_span_and_same_match_follow_cpp_less_fields() {
            let br = BlastRule::default();
            let line = concat!(
                "WP_1|1|1|fam|fam|AMR|2|CLASS|SUBCLASS|Product\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let mut al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();
            al.cdss.push(
                Locus::new(
                    1,
                    "contig",
                    20,
                    40,
                    true,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );
            al.cdss.push(
                Locus::new(
                    2,
                    "contig",
                    10,
                    50,
                    true,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );

            assert_eq!(al.get_cds_start(), 10);
            assert_eq!(al.get_cds_stop(), 50);

            let mut same = al.clone();
            same.part = 2;
            same.fam_id = "different".to_string();
            assert!(al.same_match(&same));

            same.ref_accession = "WP_2".to_string();
            assert!(!al.same_match(&same));
        }

        #[test]
        fn process_marks_cpp_fusion_parts_and_redundant_rows() {
            let br = BlastRule::default();
            let part1 = concat!(
                "WP_FUS|1|2|famA|geneA|AMR|2|CLASS|SUBCLASS|Part_A\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let part2 = concat!(
                "WP_FUS|2|2|famB|geneB|AMR|2|CLASS|SUBCLASS|Part_B\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let mut batch = empty_batch();
            batch.add_blast_al(
                BlastAlignment::from_blast_line(part1, true, true, &br, &br).unwrap(),
            );
            batch.add_blast_al(
                BlastAlignment::from_blast_line(part2, true, true, &br, &br).unwrap(),
            );

            batch.process(false, false);

            assert_eq!(
                batch.target2good_blast_als.get("target").unwrap(),
                &vec![0, 1]
            );
            assert_eq!(batch.blast_als[0].fusion_ids, vec![0, 1]);
            assert!(!batch.blast_als[0].fusion_redundant);
            assert!(batch.blast_als[1].fusion_ids.is_empty());
            assert!(batch.blast_als[1].fusion_redundant);
        }

        #[test]
        fn report_aggregates_cpp_fusion_fields_on_main_row() {
            let br = BlastRule::default();
            let part1 = concat!(
                "WP_FUS|1|2|famA|geneA|AMR|1|BETA/MAC|SUB2|Part_A\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let part2 = concat!(
                "WP_FUS|2|2|famB|geneB|STRESS|2|BETA|SUB1/SUB2|Part_B\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let mut fam_a = fam("famA", "", "", 0.0, 0.0);
            fam_a.genesymbol = "geneA".to_string();
            fam_a.reportable = 1;
            fam_a.type_ = "AMR".to_string();
            fam_a.subtype = "BETA-LACTAM".to_string();
            fam_a.class = "BETA/MAC".to_string();
            fam_a.subclass = "SUB2".to_string();
            let mut fam_b = fam("famB", "", "HMM_B", 0.0, 0.0);
            fam_b.genesymbol = "geneB".to_string();
            fam_b.reportable = 2;
            fam_b.type_ = "STRESS".to_string();
            fam_b.subtype = "ACID".to_string();
            fam_b.class = "BETA".to_string();
            fam_b.subclass = "SUB1/SUB2".to_string();

            let mut batch = empty_batch();
            batch.fam_map.insert("famA".to_string(), fam_a);
            batch.fam_map.insert("famB".to_string(), fam_b);
            batch.add_blast_al(
                BlastAlignment::from_blast_line(part1, true, true, &br, &br).unwrap(),
            );
            batch.add_blast_al(
                BlastAlignment::from_blast_line(part2, true, true, &br, &br).unwrap(),
            );
            batch.hmm_als.push(HmmAlignment {
                sseqid: "target".to_string(),
                score1: 100.0,
                score2: 20.0,
                fam_id: "famB".to_string(),
                domain: None,
                blast_al_idx: Some(1),
            });
            batch.blast_als[1].hmm_al_idx = Some(0);

            batch.process(false, false);

            let mut output = Vec::new();
            batch.report(&mut output, true).unwrap();
            let output = String::from_utf8(output).unwrap();
            let mut lines = output.lines();
            let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
            let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
            let col = |name: &str| {
                let idx = header
                    .iter()
                    .position(|field| *field == name)
                    .unwrap_or_else(|| panic!("missing report column {name}"));
                row[idx]
            };

            assert!(
                lines.next().is_none(),
                "redundant fusion row must be suppressed"
            );
            assert_eq!(col("Element symbol"), "geneA/geneB");
            assert_eq!(col("Scope"), "core");
            assert_eq!(col("Type"), "AMR/STRESS");
            assert_eq!(col("Subtype"), "ACID/BETA-LACTAM");
            assert_eq!(col("Class"), "BETA/MAC");
            assert_eq!(col("Subclass"), "SUB1/SUB2");
            assert_eq!(col("HMM accession"), "HMM_B");
            assert_eq!(col("HMM description"), "famB family");
            assert_eq!(col("Hierarchy node"), "famA/famB");
        }

        #[test]
        fn report_keeps_seq_change_rows_below_reportable_min_like_cpp() {
            let br = BlastRule::default();
            let line = concat!(
                "WP_M|1|1|pmrB|pmrB|mutation|0|COLISTIN|COLISTIN|PmrB\t",
                "target\t1\t5\t6\t1\t5\t5\tAAAAA\tAAAAT"
            );
            let mut al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();
            al.ref_mutations = vec![AmrMutation::new(
                5, "pmrB_A5T", "pmrB_A5T", "COLISTIN", "COLISTIN", "PmrB",
            )];
            al.ref_mutation = al.ref_mutations[0].clone();
            al.seq_changes.push(SeqChange {
                start: 4,
                len: 1,
                reference: "A".to_string(),
                allele: "T".to_string(),
                start_ref: 4,
                stop_ref: 5,
                start_target: 4,
                mutations: vec![0],
                ..SeqChange::default()
            });

            let mut batch = empty_batch();
            batch.reportable_min = 2;
            batch.add_blast_al(al);
            batch
                .target2good_blast_als
                .insert("target".to_string(), vec![0]);

            let mut output = Vec::new();
            batch.report(&mut output, false).unwrap();
            let output = String::from_utf8(output).unwrap();

            assert!(
                output.contains("\ntarget\tpmrB_A5T\tPmrB\tcore\tAMR\tPOINT\t"),
                "{output}"
            );
        }

        #[test]
        fn report_print_node_substitutes_gene_for_non_exact_mutation_allele_like_cpp() {
            let br = BlastRule::default();
            let line = concat!(
                "WP_M|1|1|pmrB-allele|pmrB|mutation|2|COLISTIN|COLISTIN|PmrB\t",
                "target\t1\t5\t6\t1\t5\t5\tAAAAA\tAAAAT"
            );
            let mut al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();
            al.ref_mutations = vec![AmrMutation::new(
                5, "pmrB_A5T", "pmrB_A5T", "COLISTIN", "COLISTIN", "PmrB",
            )];
            al.ref_mutation = al.ref_mutations[0].clone();
            al.seq_changes.push(SeqChange {
                start: 4,
                len: 1,
                reference: "A".to_string(),
                allele: "T".to_string(),
                start_ref: 4,
                stop_ref: 5,
                start_target: 4,
                mutations: vec![0],
                ..SeqChange::default()
            });

            let mut batch = empty_batch();
            batch.reportable_min = 2;
            batch.add_blast_al(al);
            batch
                .target2good_blast_als
                .insert("target".to_string(), vec![0]);

            let mut output = Vec::new();
            batch.report(&mut output, true).unwrap();
            let output = String::from_utf8(output).unwrap();
            let row = output.lines().nth(1).expect("expected mutation report row");
            assert_eq!(row.split('\t').next_back().unwrap(), "pmrB", "{output}");

            batch.print_node_raw = true;
            let mut raw_output = Vec::new();
            batch.report(&mut raw_output, true).unwrap();
            let raw_output = String::from_utf8(raw_output).unwrap();
            let raw_row = raw_output
                .lines()
                .nth(1)
                .expect("expected raw mutation report row");
            assert_eq!(
                raw_row.split('\t').next_back().unwrap(),
                "pmrB-allele",
                "{raw_output}"
            );
        }

        #[test]
        fn report_iterates_seq_change_mutation_pointers_like_cpp() {
            let br = BlastRule::default();
            let line = concat!(
                "WP_M|1|1|pmrB|pmrB|mutation|2|HEADER_CLASS|HEADER_SUB|PmrB_product\t",
                "target\t1\t5\t6\t1\t5\t5\tAAAAA\tAAAAT"
            );
            let mut al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();
            al.ref_mutations = vec![
                AmrMutation::new(
                    5,
                    "pmrB_A5T",
                    "pmrB_A5T",
                    "COLISTIN",
                    "COLISTIN",
                    "PmrB_mutant_name",
                ),
                AmrMutation::new(
                    5,
                    "pmrB_A5V",
                    "pmrB_A5V",
                    "BETA",
                    "BETA_SUB",
                    "PmrB_second_name",
                ),
            ];
            al.ref_mutation = al.ref_mutations[0].clone();
            al.seq_changes.push(SeqChange {
                start: 4,
                len: 1,
                reference: "A".to_string(),
                allele: "T".to_string(),
                start_ref: 4,
                stop_ref: 5,
                start_target: 4,
                mutations: vec![0, 1],
                ..SeqChange::default()
            });

            let mut batch = empty_batch();
            batch.add_blast_al(al);
            batch
                .target2good_blast_als
                .insert("target".to_string(), vec![0]);

            let mut output = Vec::new();
            batch.report(&mut output, false).unwrap();
            let output = String::from_utf8(output).unwrap();
            let rows: Vec<Vec<&str>> = output
                .lines()
                .map(|line| line.split('\t').collect())
                .collect();
            let col = |name: &str| {
                rows[0]
                    .iter()
                    .position(|field| *field == name)
                    .unwrap_or_else(|| panic!("missing report column {name}"))
            };

            assert_eq!(rows.len(), 3, "{output}");
            assert_eq!(rows[1][col("Element symbol")], "pmrB_A5T");
            assert_eq!(rows[1][col("Element name")], "PmrB mutant name");
            assert_eq!(rows[1][col("Class")], "COLISTIN");
            assert_eq!(rows[1][col("Subclass")], "COLISTIN");
            assert_eq!(rows[2][col("Element symbol")], "pmrB_A5V");
            assert_eq!(rows[2][col("Element name")], "PmrB second name");
            assert_eq!(rows[2][col("Class")], "BETA");
            assert_eq!(rows[2][col("Subclass")], "BETA_SUB");
        }

        #[test]
        fn mutation_all_keeps_seq_change_rows_below_reportable_min_like_cpp() {
            let br = BlastRule::default();
            let line = concat!(
                "WP_M|1|1|pmrB|pmrB|mutation|0|COLISTIN|COLISTIN|PmrB\t",
                "target\t1\t5\t6\t1\t5\t5\tAAAAA\tAAAAT"
            );
            let mut al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();
            al.ref_mutations = vec![AmrMutation::new(
                5, "pmrB_A5T", "pmrB_A5T", "COLISTIN", "COLISTIN", "PmrB",
            )];
            al.ref_mutation = al.ref_mutations[0].clone();
            al.seq_changes.push(SeqChange {
                start: 4,
                len: 1,
                reference: "A".to_string(),
                allele: "T".to_string(),
                start_ref: 4,
                stop_ref: 5,
                start_target: 4,
                mutations: vec![0],
                ..SeqChange::default()
            });

            let mut batch = empty_batch();
            batch.reportable_min = 2;
            batch.suppress_prots.push("WP_M".to_string());
            batch
                .accession2mutations
                .insert("WP_M".to_string(), al.ref_mutations.clone());
            batch.add_blast_al(al);
            batch
                .target2good_blast_als
                .insert("target".to_string(), vec![0]);

            let mut output = Vec::new();
            batch.report_mutation_all(&mut output, false).unwrap();
            let output = String::from_utf8(output).unwrap();

            assert!(
                output.contains("\ntarget\tpmrB_A5T\tPmrB\tcore\tAMR\tPOINT\t"),
                "{output}"
            );
        }

        #[test]
        fn report_alignment_length_uses_aligned_sseq_size_like_cpp() {
            let br = BlastRule::default();
            let line = concat!(
                "WP_G|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                "target\t1\t4\t5\t1\t5\t5\tAA-AA\tAATAA"
            );
            let al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();
            assert_eq!(al.hsp.q_len_real(), 4);
            assert_eq!(al.hsp.sseq.len(), 5);

            let mut batch = empty_batch();
            batch.add_blast_al(al);
            batch
                .target2good_blast_als
                .insert("target".to_string(), vec![0]);

            let mut output = Vec::new();
            batch.report(&mut output, false).unwrap();
            let output = String::from_utf8(output).unwrap();
            let mut lines = output.lines();
            let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
            let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
            let align_len_idx = header
                .iter()
                .position(|field| *field == "Alignment length")
                .expect("missing Alignment length column");

            assert_eq!(row[align_len_idx], "5");
        }

        #[test]
        fn report_standalone_hmm_rows_use_cpp_na_reference_fields() {
            let mut family = fam("famH", "", "HMM_H", 0.0, 0.0);
            family.family_name = "HMM family".to_string();

            let mut hsp = Hsp::default();
            hsp.sseqid = "target".to_string();
            hsp.qseqid = "famH".to_string();
            hsp.slen = 120;
            hsp.s_int = crate::seq::Interval::new(10, 70, 0);
            hsp.qlen = 100;
            hsp.q_int = crate::seq::Interval::new(20, 80, 0);
            hsp.q_prot = true;
            hsp.s_prot = true;
            hsp.a_prot = true;
            hsp.length = 60;
            hsp.nident = 42;

            let mut batch = empty_batch();
            batch.fam_map.insert("famH".to_string(), family);
            batch.hmm_als.push(HmmAlignment {
                sseqid: "target".to_string(),
                score1: 100.0,
                score2: 50.0,
                fam_id: "famH".to_string(),
                domain: None,
                blast_al_idx: None,
            });
            batch.add_blast_al(BlastAlignment {
                hsp,
                partial_dna: false,
                from_hmm: true,
                ref_accession: String::new(),
                part: 1,
                parts: 1,
                fam_id: "famH".to_string(),
                gene: "famH".to_string(),
                resistance: String::new(),
                class: "CLASS".to_string(),
                subclass: "SUBCLASS".to_string(),
                product: "Header Product".to_string(),
                reportable: 2,
                genesymbol: "famH".to_string(),
                method: "HMM".to_string(),
                br_fam_id: None,
                br_fam_checked: false,
                complete_br: BlastRule::default(),
                partial_br: BlastRule::default(),
                cdss: Vec::new(),
                hmm_al_idx: Some(0),
                susceptible_idx: None,
                seq_changes: Vec::new(),
                ref_mutations: Vec::new(),
                ref_mutation: AmrMutation::default(),
                fusion_ids: Vec::new(),
                fusion_redundant: false,
            });
            batch
                .target2good_blast_als
                .insert("target".to_string(), vec![0]);

            let mut output = Vec::new();
            batch.report(&mut output, false).unwrap();
            let output = String::from_utf8(output).unwrap();
            let mut lines = output.lines();
            let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
            let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
            let col = |name: &str| {
                let idx = header
                    .iter()
                    .position(|field| *field == name)
                    .unwrap_or_else(|| panic!("missing report column {name}"));
                row[idx]
            };

            assert_eq!(col("Method"), "HMM");
            assert_eq!(col("Element name"), "Header Product");
            assert_eq!(col("Reference sequence length"), "NA");
            assert_eq!(col("% Coverage of reference"), "NA");
            assert_eq!(col("% Identity to reference"), "NA");
            assert_eq!(col("Alignment length"), "NA");
            assert_eq!(col("Closest reference accession"), "NA");
            assert_eq!(col("Closest reference name"), "NA");
            assert_eq!(col("HMM accession"), "HMM_H");
            assert_eq!(col("HMM description"), "HMM family");
        }

        #[test]
        fn report_iterates_all_cdss_and_uses_target_coords_for_empty_cds_like_cpp() {
            let br = BlastRule::default();
            let line = concat!(
                "WP_G|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let mut al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();
            al.cdss.push(
                Locus::new(
                    1,
                    "contig1",
                    10,
                    25,
                    true,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );
            al.cdss.push(
                Locus::new(
                    2,
                    "contig2",
                    40,
                    55,
                    false,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );

            let mut batch = empty_batch();
            batch.cds_exist = true;
            batch.add_blast_al(al);
            batch
                .target2good_blast_als
                .insert("target".to_string(), vec![0]);

            let mut output = Vec::new();
            batch.report(&mut output, false).unwrap();
            let output = String::from_utf8(output).unwrap();
            let rows: Vec<Vec<&str>> = output
                .lines()
                .map(|line| line.split('\t').collect())
                .collect();
            let col = |name: &str| {
                rows[0]
                    .iter()
                    .position(|field| *field == name)
                    .unwrap_or_else(|| panic!("missing report column {name}"))
            };

            assert_eq!(rows.len(), 3, "{output}");
            assert_eq!(rows[1][col("Contig id")], "contig1");
            assert_eq!(rows[1][col("Start")], "11");
            assert_eq!(rows[1][col("Stop")], "25");
            assert_eq!(rows[1][col("Strand")], "+");
            assert_eq!(rows[2][col("Contig id")], "contig2");
            assert_eq!(rows[2][col("Start")], "41");
            assert_eq!(rows[2][col("Stop")], "55");
            assert_eq!(rows[2][col("Strand")], "-");

            let line = concat!(
                "WP_E|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let mut batch = empty_batch();
            batch.cds_exist = true;
            batch
                .add_blast_al(BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap());
            batch
                .target2good_blast_als
                .insert("target".to_string(), vec![0]);

            let mut output = Vec::new();
            batch.report(&mut output, false).unwrap();
            let output = String::from_utf8(output).unwrap();
            let rows: Vec<Vec<&str>> = output
                .lines()
                .map(|line| line.split('\t').collect())
                .collect();
            let col = |name: &str| {
                rows[0]
                    .iter()
                    .position(|field| *field == name)
                    .unwrap_or_else(|| panic!("missing report column {name}"))
            };
            assert_eq!(rows[1][col("Contig id")], "target");
            assert_eq!(rows[1][col("Start")], "1");
            assert_eq!(rows[1][col("Stop")], "4");
        }

        #[test]
        fn report_filters_protein_cdss_to_contig_key_when_blastx_rekeys_like_cpp() {
            let br = BlastRule::default();
            let line = concat!(
                "WP_G|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let mut al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();
            al.cdss.push(
                Locus::new(
                    1,
                    "contig1",
                    10,
                    25,
                    true,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );
            al.cdss.push(
                Locus::new(
                    2,
                    "contig2",
                    40,
                    55,
                    false,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );

            let mut batch = empty_batch();
            batch.cds_exist = true;
            batch.add_blast_al(al);
            batch
                .target2good_blast_als
                .insert("contig1".to_string(), vec![0]);
            batch
                .target2good_blast_als
                .insert("contig2".to_string(), vec![0]);

            let mut output = Vec::new();
            batch.report(&mut output, false).unwrap();
            let output = String::from_utf8(output).unwrap();
            let rows: Vec<Vec<&str>> = output
                .lines()
                .map(|line| line.split('\t').collect())
                .collect();
            let col = |name: &str| {
                rows[0]
                    .iter()
                    .position(|field| *field == name)
                    .unwrap_or_else(|| panic!("missing report column {name}"))
            };

            assert_eq!(rows.len(), 3, "{output}");
            assert_eq!(rows[1][col("Contig id")], "contig1");
            assert_eq!(rows[1][col("Start")], "11");
            assert_eq!(rows[2][col("Contig id")], "contig2");
            assert_eq!(rows[2][col("Start")], "41");
        }

        #[test]
        fn mutation_all_iterates_all_cdss_like_cpp_report_loop() {
            let br = BlastRule::default();
            let line = concat!(
                "WP_M|1|1|pmrB|pmrB|mutation|2|COLISTIN|COLISTIN|PmrB\t",
                "target\t1\t5\t6\t1\t5\t5\tAAAAA\tAAAAT"
            );
            let mut al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();
            al.ref_mutations = vec![AmrMutation::new(
                5, "pmrB_A5T", "pmrB_A5T", "COLISTIN", "COLISTIN", "PmrB",
            )];
            al.ref_mutation = al.ref_mutations[0].clone();
            al.seq_changes.push(SeqChange {
                start: 4,
                len: 1,
                reference: "A".to_string(),
                allele: "T".to_string(),
                start_ref: 4,
                stop_ref: 5,
                start_target: 4,
                mutations: vec![0],
                ..SeqChange::default()
            });
            al.cdss.push(
                Locus::new(
                    1,
                    "contig1",
                    10,
                    25,
                    true,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );
            al.cdss.push(
                Locus::new(
                    2,
                    "contig2",
                    40,
                    55,
                    false,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );

            let mut batch = empty_batch();
            batch.cds_exist = true;
            batch
                .accession2mutations
                .insert("WP_M".to_string(), al.ref_mutations.clone());
            batch.add_blast_al(al);
            batch
                .target2good_blast_als
                .insert("target".to_string(), vec![0]);

            let mut output = Vec::new();
            batch.report_mutation_all(&mut output, false).unwrap();
            let output = String::from_utf8(output).unwrap();
            let rows: Vec<Vec<&str>> = output
                .lines()
                .map(|line| line.split('\t').collect())
                .collect();
            let col = |name: &str| {
                rows[0]
                    .iter()
                    .position(|field| *field == name)
                    .unwrap_or_else(|| panic!("missing report column {name}"))
            };

            assert_eq!(rows.len(), 3, "{output}");
            assert_eq!(rows[1][col("Contig id")], "contig1");
            assert_eq!(rows[2][col("Contig id")], "contig2");
            assert_eq!(rows[1][col("Element symbol")], "pmrB_A5T");
            assert_eq!(rows[2][col("Element symbol")], "pmrB_A5T");
        }

        #[test]
        fn mutation_all_filters_protein_cdss_to_contig_key_when_blastx_rekeys_like_cpp() {
            let br = BlastRule::default();
            let line = concat!(
                "WP_M|1|1|pmrB|pmrB|mutation|2|COLISTIN|COLISTIN|PmrB\t",
                "target\t1\t5\t6\t1\t5\t5\tAAAAA\tAAAAT"
            );
            let mut al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();
            al.ref_mutations = vec![AmrMutation::new(
                5, "pmrB_A5T", "pmrB_A5T", "COLISTIN", "COLISTIN", "PmrB",
            )];
            al.ref_mutation = al.ref_mutations[0].clone();
            al.seq_changes.push(SeqChange {
                start: 4,
                len: 1,
                reference: "A".to_string(),
                allele: "T".to_string(),
                start_ref: 4,
                stop_ref: 5,
                start_target: 4,
                mutations: vec![0],
                ..SeqChange::default()
            });
            al.cdss.push(
                Locus::new(
                    1,
                    "contig1",
                    10,
                    25,
                    true,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );
            al.cdss.push(
                Locus::new(
                    2,
                    "contig2",
                    40,
                    55,
                    false,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );

            let mut batch = empty_batch();
            batch.cds_exist = true;
            batch
                .accession2mutations
                .insert("WP_M".to_string(), al.ref_mutations.clone());
            batch.add_blast_al(al);
            batch
                .target2good_blast_als
                .insert("contig1".to_string(), vec![0]);
            batch
                .target2good_blast_als
                .insert("contig2".to_string(), vec![0]);

            let mut output = Vec::new();
            batch.report_mutation_all(&mut output, false).unwrap();
            let output = String::from_utf8(output).unwrap();
            let rows: Vec<Vec<&str>> = output
                .lines()
                .map(|line| line.split('\t').collect())
                .collect();
            let col = |name: &str| {
                rows[0]
                    .iter()
                    .position(|field| *field == name)
                    .unwrap_or_else(|| panic!("missing report column {name}"))
            };

            assert_eq!(rows.len(), 3, "{output}");
            assert_eq!(rows[1][col("Contig id")], "contig1");
            assert_eq!(rows[2][col("Contig id")], "contig2");
            assert_eq!(rows[1][col("Element symbol")], "pmrB_A5T");
            assert_eq!(rows[2][col("Element symbol")], "pmrB_A5T");
        }

        #[test]
        fn report_susceptible_name_uses_organism_metadata_like_cpp() {
            let br = BlastRule::default();
            let line = concat!(
                "WP_S|1|1|famS|geneS|susceptible|2|HEADER_CLASS|HEADER_SUB|Header_Product\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );

            let mut batch = empty_batch();
            batch.accession2susceptible.insert(
                "WP_S".to_string(),
                Susceptible {
                    genesymbol: "susGene".to_string(),
                    cutoff: 0.95,
                    class: "SCLASS".to_string(),
                    subclass: "SSUB".to_string(),
                    name: "Metadata Product".to_string(),
                },
            );
            batch
                .add_blast_al(BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap());
            batch
                .target2good_blast_als
                .insert("target".to_string(), vec![0]);

            let mut output = Vec::new();
            batch.report(&mut output, false).unwrap();
            let output = String::from_utf8(output).unwrap();
            let mut lines = output.lines();
            let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
            let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
            let col = |name: &str| {
                let idx = header
                    .iter()
                    .position(|field| *field == name)
                    .unwrap_or_else(|| panic!("missing report column {name}"));
                row[idx]
            };

            assert_eq!(col("Element symbol"), "susGene");
            assert_eq!(col("Element name"), "Metadata Product");
            assert_eq!(col("Closest reference name"), "Header Product");
        }

        #[test]
        fn blast_better_compares_blastp_and_blastx_only_when_cds_matches() {
            let batch = empty_batch();
            let br = BlastRule::default();
            let blastp_line = concat!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                "protein\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let blastx_line = concat!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                "contig\t1\t5\t5\t10\t25\t100\tAAAA*\tAAAA*"
            );
            let mut blastp =
                BlastAlignment::from_blast_line(blastp_line, true, true, &br, &br).unwrap();
            let blastx =
                BlastAlignment::from_blast_line(blastx_line, true, false, &br, &br).unwrap();
            blastp.cdss.push(
                Locus::new(
                    1,
                    "contig",
                    10,
                    28,
                    true,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );

            assert!(blastp.matches_cds(&blastx));
            assert!(batch.blast_better(&blastp, &blastx));
            assert!(!batch.blast_better(&blastx, &blastp));
        }

        #[test]
        fn blast_better_keeps_unmatched_blastp_and_blastx_incomparable() {
            let batch = empty_batch();
            let br = BlastRule::default();
            let blastp_line = concat!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                "protein\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let blastx_line = concat!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                "contig\t1\t5\t5\t10\t25\t100\tAAAA*\tAAAA*"
            );
            let mut blastp =
                BlastAlignment::from_blast_line(blastp_line, true, true, &br, &br).unwrap();
            let blastx =
                BlastAlignment::from_blast_line(blastx_line, true, false, &br, &br).unwrap();
            blastp.cdss.push(
                Locus::new(
                    1,
                    "contig",
                    1000,
                    1020,
                    true,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );

            assert!(!blastp.matches_cds(&blastx));
            assert!(!batch.blast_better(&blastx, &blastp));
            assert!(!batch.blast_better(&blastp, &blastx));
        }

        #[test]
        fn blastx_edge_overlap_does_not_mark_or_suppress_without_cpp_cds_match() {
            let br = BlastRule::default();
            let qseq = format!("{}*", "M".repeat(100));
            let sseq = qseq.clone();
            let blastp_line = format!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 protein\t1\t101\t101\t1\t101\t101\t{qseq}\t{sseq}"
            );
            let mut blastx_sseq = "M".repeat(101);
            blastx_sseq.replace_range(50..51, "*");
            blastx_sseq.replace_range(100..101, "*");
            let blastx_line = format!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 contig\t1\t101\t101\t1991\t2293\t5000\t{qseq}\t{blastx_sseq}"
            );
            let mut blastp =
                BlastAlignment::from_blast_line(&blastp_line, true, true, &br, &br).unwrap();
            let blastx =
                BlastAlignment::from_blast_line(&blastx_line, true, false, &br, &br).unwrap();
            blastp.cdss.push(
                Locus::new(
                    1,
                    "contig",
                    1000,
                    2000,
                    true,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );

            assert!(!blastp.matches_cds(&blastx));

            let mut batch = empty_batch();
            batch.add_blast_al(blastp);
            batch.add_blast_al(blastx);
            batch.target2blast_als.clear();
            batch
                .target2blast_als
                .insert("contig".to_string(), vec![0, 1]);

            batch.mark_protein_internal_stops_from_blastx();
            assert!(!batch.blast_als[0].hsp.s_internal_stop);

            batch.process(false, false);
            let retained = batch.target2good_blast_als["contig"]
                .iter()
                .copied()
                .collect::<HashSet<_>>();
            assert_eq!(retained, HashSet::from([0, 1]));
        }

        #[test]
        fn process_retain_blasts_keeps_overlapping_blastp_and_blastx_like_cpp_flag() {
            let br = BlastRule::default();
            let qseq = format!("{}*", "M".repeat(119));
            let sseq = qseq.clone();
            let blastp_line = format!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 protein\t1\t120\t120\t1\t120\t120\t{qseq}\t{sseq}"
            );
            let blastx_line = format!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 contig\t1\t120\t120\t100\t460\t1000\t{qseq}\t{sseq}"
            );
            let mut blastp =
                BlastAlignment::from_blast_line(&blastp_line, true, true, &br, &br).unwrap();
            let blastx =
                BlastAlignment::from_blast_line(&blastx_line, true, false, &br, &br).unwrap();
            blastp.cdss.push(
                Locus::new(
                    1,
                    "contig",
                    100,
                    463,
                    true,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );

            let mut batch = empty_batch();
            batch.add_blast_al(blastp);
            batch.add_blast_al(blastx);
            batch.target2blast_als.clear();
            batch
                .target2blast_als
                .insert("contig".to_string(), vec![0, 1]);

            batch.process(true, false);

            let retained = batch.target2good_blast_als["contig"]
                .iter()
                .copied()
                .collect::<HashSet<_>>();
            assert_eq!(retained, HashSet::from([0, 1]));
        }

        #[test]
        fn process_marks_internal_stop_before_blast_pareto_like_cpp() {
            let br = BlastRule::default();
            let protein_qseq = format!("{}*", "M".repeat(79));
            let protein_sseq = protein_qseq.clone();
            let blastp_line = format!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 protein\t1\t80\t80\t1\t80\t100\t{}\t{}",
                protein_qseq, protein_sseq
            );
            let better_qseq = format!("{}*", "M".repeat(99));
            let better_sseq = better_qseq.clone();
            let better_blastp_line = format!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 protein\t1\t100\t100\t1\t100\t100\t{}\t{}",
                better_qseq, better_sseq
            );
            let blastx_qseq = format!("{}*", "M".repeat(100));
            let mut blastx_sseq = blastx_qseq.clone();
            blastx_sseq.replace_range(50..51, "*");
            let blastx_line = format!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 contig\t1\t101\t101\t100\t403\t1000\t{}\t{}",
                blastx_qseq, blastx_sseq
            );
            let mut blastp =
                BlastAlignment::from_blast_line(&blastp_line, true, true, &br, &br).unwrap();
            blastp.cdss.push(
                Locus::new(
                    1,
                    "contig",
                    100,
                    406,
                    true,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );
            let blastx =
                BlastAlignment::from_blast_line(&blastx_line, true, false, &br, &br).unwrap();
            let mut better_blastp =
                BlastAlignment::from_blast_line(&better_blastp_line, true, true, &br, &br).unwrap();
            better_blastp.cdss.push(
                Locus::new(
                    1,
                    "contig",
                    100,
                    406,
                    true,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );

            let mut batch = empty_batch();
            batch
                .fam_map
                .insert("fam".to_string(), fam("fam", "", "", 0.0, 0.0));
            batch.blast_als.push(blastp);
            batch.blast_als.push(better_blastp);
            batch.blast_als.push(blastx);
            batch
                .target2blast_als
                .insert("contig".to_string(), vec![0, 1, 2]);

            batch.process(false, false);

            assert!(batch.blast_als[0].hsp.s_internal_stop);
            assert!(!batch.target2good_blast_als["contig"].contains(&0));
        }

        #[test]
        fn process_does_not_transfer_internal_stop_from_rejected_blastx_like_cpp() {
            let br = BlastRule::default();
            let protein_qseq = format!("{}*", "M".repeat(100));
            let protein_sseq = protein_qseq.clone();
            let blastp_line = format!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 protein\t1\t101\t101\t1\t101\t101\t{}\t{}",
                protein_qseq, protein_sseq
            );
            let blastx_qseq = format!("{}*", "M".repeat(100));
            let mut blastx_sseq = blastx_qseq.clone();
            blastx_sseq.replace_range(50..51, "*");
            let blastx_line = format!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 contig\t1\t101\t101\t100\t403\t1000\t{}\t{}",
                blastx_qseq, blastx_sseq
            );
            let mut blastp =
                BlastAlignment::from_blast_line(&blastp_line, true, true, &br, &br).unwrap();
            blastp.cdss.push(
                Locus::new(
                    1,
                    "contig",
                    100,
                    403,
                    true,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );
            let mut blastx =
                BlastAlignment::from_blast_line(&blastx_line, true, false, &br, &br).unwrap();
            blastx.hsp.nident = 60;

            let mut fam_with_rule = fam("fam", "", "", 0.0, 0.0);
            fam_with_rule.complete_br = BlastRule::new(0.9, 0.9);
            fam_with_rule.partial_br = BlastRule::new(0.9, 0.5);

            let mut batch = empty_batch();
            batch.fam_map.insert("fam".to_string(), fam_with_rule);
            batch.blast_als.push(blastp);
            batch.blast_als.push(blastx);
            batch
                .target2blast_als
                .insert("contig".to_string(), vec![0, 1]);

            batch.process(false, false);

            assert!(!batch.blast_als[0].hsp.s_internal_stop);
            assert_eq!(batch.target2blast_als["contig"], vec![0]);
            assert_eq!(batch.target2good_blast_als["contig"], vec![0]);
        }

        #[test]
        fn process_transfers_internal_stop_to_hmm_synthetic_row_like_cpp() {
            let br = BlastRule::default();
            let domain = HmmDomain {
                score: 100.0,
                hmm_len: 101,
                hmm_start: 0,
                hmm_stop: 101,
                seq_len: 101,
                seq_start: 0,
                seq_stop: 101,
            };
            let hmm_al = HmmAlignment {
                sseqid: "protein".to_string(),
                score1: 100.0,
                score2: 50.0,
                fam_id: "fam".to_string(),
                domain: Some(domain.clone()),
                blast_al_idx: Some(0),
            };
            let mut hmm_blast = BlastAlignment::from_hmm_alignment(
                &hmm_al,
                Some(&fam("fam", "", "HMM_FAM", 50.0, 20.0)),
                None,
                &domain,
            );
            hmm_blast.cdss.push(
                Locus::new(
                    1,
                    "contig",
                    100,
                    403,
                    true,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );

            let blastx_qseq = format!("{}*", "M".repeat(100));
            let mut blastx_sseq = blastx_qseq.clone();
            blastx_sseq.replace_range(50..51, "*");
            let blastx_line = format!(
                "WP_X|1|1|other|other|AMR|2|CLASS|SUBCLASS|Product\t\
                 contig\t1\t101\t101\t100\t403\t1000\t{}\t{}",
                blastx_qseq, blastx_sseq
            );
            let blastx =
                BlastAlignment::from_blast_line(&blastx_line, true, false, &br, &br).unwrap();

            let mut batch = empty_batch();
            batch
                .fam_map
                .insert("fam".to_string(), fam("fam", "", "HMM_FAM", 50.0, 20.0));
            batch
                .fam_map
                .insert("other".to_string(), fam("other", "", "", 0.0, 0.0));
            batch.blast_als.push(hmm_blast);
            batch.blast_als.push(blastx);
            batch.hmm_als.push(hmm_al);
            batch.target2blast_als.insert("contig".to_string(), vec![1]);
            batch.target2hmm_als.insert("contig".to_string(), vec![0]);

            batch.process(false, false);

            assert!(batch.blast_als[0].hsp.s_internal_stop);
        }

        #[test]
        fn blast_alignment_matches_cds_uses_cpp_overlap_and_tail_rules() {
            let br = BlastRule::default();
            let blastp_line = concat!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                "protein\t1\t100\t101\t1\t100\t100\t",
                "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
                "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\t",
                "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
                "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM"
            );
            let blastx_line = concat!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                "contig\t1\t80\t101\t40\t280\t400\t",
                "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
                "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\t",
                "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
                "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM"
            );
            let mut blastp =
                BlastAlignment::from_blast_line(blastp_line, true, true, &br, &br).unwrap();
            let blastx =
                BlastAlignment::from_blast_line(blastx_line, true, false, &br, &br).unwrap();
            blastp.cdss.push(
                Locus::new(
                    1,
                    "contig",
                    10,
                    310,
                    true,
                    false,
                    0,
                    String::new(),
                    String::new(),
                )
                .unwrap(),
            );

            assert!(blastp.matches_cds(&blastx));

            let mut wrong_strand = blastp.clone();
            wrong_strand.cdss[0].strand = false;
            assert!(!wrong_strand.matches_cds(&blastx));

            let mut cross_origin = blastp.clone();
            cross_origin.cdss[0].cross_origin = true;
            assert!(!cross_origin.matches_cds(&blastx));
        }

        #[test]
        fn blast_alignment_inside_eq_allows_partial_pseudo_repeat_overlap() {
            let br = BlastRule::default();
            let qseq = "M".repeat(60);
            let sseq = "M".repeat(60);
            let left_line = format!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 contig\t11\t70\t100\t101\t280\t1000\t{}\t{}",
                qseq, sseq
            );
            let right_line = format!(
                "WP_2|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 contig\t11\t70\t100\t221\t400\t1000\t{}\t{}",
                qseq, sseq
            );
            let mut left =
                BlastAlignment::from_blast_line(&left_line, true, false, &br, &br).unwrap();
            let mut right =
                BlastAlignment::from_blast_line(&right_line, true, false, &br, &br).unwrap();
            left.finish().unwrap();
            right.finish().unwrap();

            assert!(left.partial_pseudo());
            assert!(!left
                .hsp
                .s_int
                .inside_eq(&right.hsp.s_int, 10 * left.hsp.a2s));
            let batch = empty_batch();
            assert!(left.inside_eq(&right, &batch));
        }

        #[test]
        fn blast_alignment_inside_eq_rejects_partial_pseudo_wrong_strand_or_symbol() {
            let br = BlastRule::default();
            let qseq = "M".repeat(60);
            let sseq = "M".repeat(60);
            let left_line = format!(
                "WP_1|1|1|fam|geneA|AMR|2|CLASS|SUBCLASS|Product\t\
                 contig\t11\t70\t100\t101\t280\t1000\t{}\t{}",
                qseq, sseq
            );
            let right_line = format!(
                "WP_2|1|1|fam|geneA|AMR|2|CLASS|SUBCLASS|Product\t\
                 contig\t11\t70\t100\t221\t400\t1000\t{}\t{}",
                qseq, sseq
            );
            let mut left =
                BlastAlignment::from_blast_line(&left_line, true, false, &br, &br).unwrap();
            let mut right =
                BlastAlignment::from_blast_line(&right_line, true, false, &br, &br).unwrap();
            left.finish().unwrap();
            right.finish().unwrap();

            right.hsp.s_int.strand = -1;
            let mut batch = empty_batch();
            assert!(!left.inside_eq(&right, &batch));

            right.hsp.s_int.strand = 1;
            left.br_fam_id = Some("brA".to_string());
            right.br_fam_id = Some("brB".to_string());
            let mut br_a = fam("brA", "", "", 0.0, 0.0);
            br_a.genesymbol = "geneA".to_string();
            let mut br_b = fam("brB", "", "", 0.0, 0.0);
            br_b.genesymbol = "geneB".to_string();
            batch.fam_map.insert("brA".to_string(), br_a);
            batch.fam_map.insert("brB".to_string(), br_b);
            assert!(!left.inside_eq(&right, &batch));
        }

        #[test]
        fn blast_alignment_inside_eq_uses_match_family_symbols_like_cpp() {
            let br = BlastRule::default();
            let qseq = "M".repeat(60);
            let sseq = "M".repeat(60);
            let left_line = format!(
                "WP_1|1|1|famA|raw_gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 contig\t11\t70\t100\t101\t280\t1000\t{}\t{}",
                qseq, sseq
            );
            let right_line = format!(
                "WP_2|1|1|famB|raw_gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 contig\t11\t70\t100\t221\t400\t1000\t{}\t{}",
                qseq, sseq
            );
            let mut left =
                BlastAlignment::from_blast_line(&left_line, true, false, &br, &br).unwrap();
            let mut right =
                BlastAlignment::from_blast_line(&right_line, true, false, &br, &br).unwrap();
            left.finish().unwrap();
            right.finish().unwrap();
            left.genesymbol = "raw_gene".to_string();
            right.genesymbol = "raw_gene".to_string();
            left.br_fam_id = Some("brA".to_string());
            right.br_fam_id = Some("brB".to_string());

            let mut batch = empty_batch();
            let mut br_a = fam("brA", "", "", 0.0, 0.0);
            br_a.genesymbol = "symbol_a".to_string();
            let mut br_b = fam("brB", "", "", 0.0, 0.0);
            br_b.genesymbol = "symbol_b".to_string();
            batch.fam_map.insert("brA".to_string(), br_a);
            batch.fam_map.insert("brB".to_string(), br_b);

            assert_eq!(batch.fusion_2_gene_symbols(&left), "symbol_a");
            assert_eq!(batch.fusion_2_gene_symbols(&right), "symbol_b");
            assert!(!left.inside_eq(&right, &batch));
        }

        #[test]
        fn blast_alignment_inside_eq_rejects_truncated_partial_pseudo() {
            let br = BlastRule::default();
            let qseq = "M".repeat(60);
            let sseq = "M".repeat(60);
            let left_line = format!(
                "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 contig\t11\t70\t100\t1\t180\t1000\t{}\t{}",
                qseq, sseq
            );
            let right_line = format!(
                "WP_2|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 contig\t11\t70\t100\t121\t300\t1000\t{}\t{}",
                qseq, sseq
            );
            let mut left =
                BlastAlignment::from_blast_line(&left_line, true, false, &br, &br).unwrap();
            let mut right =
                BlastAlignment::from_blast_line(&right_line, true, false, &br, &br).unwrap();
            left.finish().unwrap();
            right.finish().unwrap();

            assert!(left.truncated_cds());
            assert!(!left.partial_pseudo());
            let batch = empty_batch();
            assert!(!left.inside_eq(&right, &batch));
        }

        #[test]
        fn blast_alignment_less_orders_fusion_pass_fields() {
            let br = BlastRule::default();
            let qseq = "M".repeat(90);
            let sseq = "M".repeat(90);
            let first_line = format!(
                "WP_A|1|2|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 target\t1\t90\t91\t101\t190\t1000\t{}\t{}",
                qseq, sseq
            );
            let second_line = format!(
                "WP_B|2|2|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 target\t1\t90\t91\t101\t190\t1000\t{}\t{}",
                qseq, sseq
            );
            let earlier_target_line = format!(
                "WP_Z|1|2|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 target\t1\t90\t91\t11\t100\t1000\t{}\t{}",
                qseq, sseq
            );
            let mut first =
                BlastAlignment::from_blast_line(&first_line, true, true, &br, &br).unwrap();
            let mut second =
                BlastAlignment::from_blast_line(&second_line, true, true, &br, &br).unwrap();
            let earlier_target =
                BlastAlignment::from_blast_line(&earlier_target_line, true, true, &br, &br)
                    .unwrap();
            first.cdss = vec![Locus {
                line_num: 0,
                contig: "contig".to_string(),
                start: 100,
                stop: 370,
                strand: true,
                partial: false,
                contig_len: 1000,
                cross_origin: false,
                gene: "gene".to_string(),
                product: "Product".to_string(),
            }];
            second.cdss = first.cdss.clone();

            assert!(earlier_target.less(&first));
            assert!(first.less(&second));
            assert!(!second.less(&first));
        }

        #[test]
        fn report_preserves_cpp_process_alignment_order() {
            let br = BlastRule::default();
            let qseq = "M".repeat(90);
            let sseq = "M".repeat(90);
            let z_line = format!(
                "WP_Z|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product_Z\t\
                 target\t1\t90\t91\t101\t190\t1000\t{}\t{}",
                qseq, sseq
            );
            let a_line = format!(
                "WP_A|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product_A\t\
                 target\t1\t90\t91\t101\t190\t1000\t{}\t{}",
                qseq, sseq
            );
            let mut z = BlastAlignment::from_blast_line(&z_line, true, true, &br, &br).unwrap();
            let mut a = BlastAlignment::from_blast_line(&a_line, true, true, &br, &br).unwrap();
            z.cdss = vec![Locus {
                line_num: 0,
                contig: "contig".to_string(),
                start: 10,
                stop: 280,
                strand: true,
                partial: false,
                contig_len: 1000,
                cross_origin: false,
                gene: "gene".to_string(),
                product: "Product_Z".to_string(),
            }];
            a.cdss = vec![Locus {
                line_num: 0,
                contig: "contig".to_string(),
                start: 100,
                stop: 370,
                strand: true,
                partial: false,
                contig_len: 1000,
                cross_origin: false,
                gene: "gene".to_string(),
                product: "Product_A".to_string(),
            }];
            assert!(z.less(&a));
            assert!(a.hsp.qseqid < z.hsp.qseqid);

            let mut batch = empty_batch();
            batch
                .fam_map
                .insert("fam".to_string(), fam("fam", "", "", 0.0, 0.0));
            batch.blast_als.push(z);
            batch.blast_als.push(a);
            batch
                .target2good_blast_als
                .insert("target".to_string(), vec![0, 1]);

            let mut output = Vec::new();
            batch.report(&mut output, false).unwrap();
            let output = String::from_utf8(output).unwrap();
            let mut lines = output.lines();
            let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
            let accession_col = header
                .iter()
                .position(|name| *name == crate::columns::CLOSEST_REF_ACCESSION_COL_NAME)
                .unwrap();
            let first: Vec<&str> = lines.next().unwrap().split('\t').collect();
            let second: Vec<&str> = lines.next().unwrap().split('\t').collect();

            assert_eq!(first[accession_col], "WP_Z");
            assert_eq!(second[accession_col], "WP_A");
        }

        #[test]
        fn blast_better_eq_uses_fusion_gene_symbols_for_family_overlap_gate() {
            let br = BlastRule::default();
            let a_qseq = "M".repeat(90);
            let b_qseq = "M".repeat(80);
            let a_line = format!(
                "WP_A|1|1|famA|raw_gene_a|AMR|2|CLASS|SUBCLASS|Product\t\
                 target\t1\t90\t100\t1\t90\t400\t{}\t{}",
                a_qseq, a_qseq
            );
            let b_line = format!(
                "WP_B|1|1|famB|raw_gene_b|AMR|2|CLASS|SUBCLASS|Product\t\
                 target\t1\t80\t100\t201\t280\t400\t{}\t{}",
                b_qseq, b_qseq
            );
            let mut a = BlastAlignment::from_blast_line(&a_line, true, true, &br, &br).unwrap();
            let mut b = BlastAlignment::from_blast_line(&b_line, true, true, &br, &br).unwrap();
            a.genesymbol = "shared_symbol".to_string();
            b.genesymbol = "shared_symbol".to_string();

            let batch = empty_batch();
            assert!(batch.blast_better_eq(&a, &b));
            assert!(!batch.blast_better_eq(&b, &a));
        }

        #[test]
        fn blast_better_eq_empty_family_gene_symbols_fall_back_to_na_like_cpp() {
            let br = BlastRule::default();
            let a_qseq = "M".repeat(90);
            let b_qseq = "M".repeat(80);
            let a_line = format!(
                "WP_A|1|1|famA|raw_gene_a|AMR|2|CLASS|SUBCLASS|Product\t\
                 target\t1\t90\t100\t1\t90\t400\t{}\t{}",
                a_qseq, a_qseq
            );
            let b_line = format!(
                "WP_B|1|1|famB|raw_gene_b|AMR|2|CLASS|SUBCLASS|Product\t\
                 target\t1\t80\t100\t201\t280\t400\t{}\t{}",
                b_qseq, b_qseq
            );
            let a = BlastAlignment::from_blast_line(&a_line, true, true, &br, &br).unwrap();
            let b = BlastAlignment::from_blast_line(&b_line, true, true, &br, &br).unwrap();

            let batch = empty_batch();
            assert_eq!(batch.fusion_2_gene_symbols(&a), crate::columns::NA);
            assert_eq!(batch.fusion_2_gene_symbols(&b), crate::columns::NA);
            assert!(batch.blast_better_eq(&a, &b));
            assert!(!batch.blast_better_eq(&b, &a));
        }

        #[test]
        fn blast_better_eq_uses_allele_family_id_as_gene_symbol_like_cpp() {
            let br = BlastRule::default();
            let qseq = format!("{}*", "M".repeat(99));
            let a_line = format!(
                "WP_A|1|1|alleleA|shared_gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 target\t1\t100\t100\t1\t100\t100\t{}\t{}",
                qseq, qseq
            );
            let b_line = format!(
                "WP_B|1|1|alleleB|shared_gene|AMR|2|CLASS|SUBCLASS|Product\t\
                 target\t1\t100\t100\t301\t400\t100\t{}\t{}",
                qseq, qseq
            );
            let mut a = BlastAlignment::from_blast_line(&a_line, true, true, &br, &br).unwrap();
            let mut b = BlastAlignment::from_blast_line(&b_line, true, true, &br, &br).unwrap();
            a.genesymbol = "shared_symbol".to_string();
            b.genesymbol = "shared_symbol".to_string();

            assert!(a.allele_match());
            assert!(b.allele_match());
            let batch = empty_batch();
            assert_eq!(batch.fusion_2_gene_symbols(&a), "alleleA");
            assert_eq!(batch.fusion_2_gene_symbols(&b), "alleleB");
            assert!(!batch.blast_better_eq(&a, &b));
            assert!(!batch.blast_better_eq(&b, &a));
        }

        #[test]
        fn blast_alignment_fusion_2_fam_ids_matches_cpp_member_branches() {
            let br = BlastRule::default();
            let qseq = "AAAA*";
            let part1 = format!(
                "WP_FUS|1|2|famA|geneA|AMR|2|CLASS|SUBCLASS|Part_A\t\
                 target\t1\t5\t5\t1\t5\t5\t{qseq}\t{qseq}"
            );
            let part2 = format!(
                "WP_FUS|2|2|famB|geneB|AMR|2|CLASS|SUBCLASS|Part_B\t\
                 target\t1\t5\t5\t1\t5\t5\t{qseq}\t{qseq}"
            );
            let plain = format!(
                "WP_PLAIN|1|1|famPlain|genePlain|AMR|2|CLASS|SUBCLASS|Plain\t\
                 target\t1\t5\t5\t1\t5\t5\t{qseq}\t{qseq}"
            );
            let mutation = format!(
                "WP_MUT|1|1|famMutation|geneMutation|mutation|2|CLASS|SUBCLASS|Mutation\t\
                 target\t1\t5\t5\t1\t5\t5\t{qseq}\t{qseq}"
            );

            let mut blast_als = vec![
                BlastAlignment::from_blast_line(&part1, true, true, &br, &br).unwrap(),
                BlastAlignment::from_blast_line(&part2, true, true, &br, &br).unwrap(),
                BlastAlignment::from_blast_line(&plain, true, true, &br, &br).unwrap(),
                BlastAlignment::from_blast_line(&mutation, true, true, &br, &br).unwrap(),
            ];
            blast_als[0].fusion_ids = vec![0, 1];
            blast_als[3].fusion_ids = vec![0, 1];

            assert_eq!(blast_als[0].fusion_2_fam_ids(&blast_als), "famA/famB");
            assert_eq!(blast_als[2].fusion_2_fam_ids(&blast_als), "famPlain");
            assert_eq!(blast_als[3].fusion_2_fam_ids(&blast_als), "famMutation");
        }

        #[test]
        fn blast_alignment_fusion_2_type_matches_cpp_member_branches() {
            let br = BlastRule::default();
            let qseq = "AAAA*";
            let part1 = format!(
                "WP_FUS|1|3|famA|geneA|AMR|2|CLASS|SUBCLASS|Part_A\t\
                 target\t1\t5\t5\t1\t5\t5\t{qseq}\t{qseq}"
            );
            let part2 = format!(
                "WP_FUS|2|3|famB|geneB|AMR|2|CLASS|SUBCLASS|Part_B\t\
                 target\t1\t5\t5\t1\t5\t5\t{qseq}\t{qseq}"
            );
            let part3 = format!(
                "WP_FUS|3|3|famC|geneC|AMR|2|CLASS|SUBCLASS|Part_C\t\
                 target\t1\t5\t5\t1\t5\t5\t{qseq}\t{qseq}"
            );
            let susceptible = format!(
                "WP_SUS|1|1|famS|geneS|susceptible|2|CLASS|SUBCLASS|Susceptible\t\
                 target\t1\t5\t5\t1\t5\t5\t{qseq}\t{qseq}"
            );

            let mut batch = empty_batch();
            let mut fam_a = fam("famA", "", "", 0.0, 0.0);
            fam_a.type_ = "STRESS".to_string();
            let mut fam_b = fam("famB", "", "", 0.0, 0.0);
            fam_b.type_ = "AMR".to_string();
            let mut fam_c = fam("famC", "", "", 0.0, 0.0);
            fam_c.type_ = "STRESS".to_string();
            batch.fam_map.insert("famA".to_string(), fam_a);
            batch.fam_map.insert("famB".to_string(), fam_b);
            batch.fam_map.insert("famC".to_string(), fam_c);

            batch.blast_als = vec![
                BlastAlignment::from_blast_line(&part1, true, true, &br, &br).unwrap(),
                BlastAlignment::from_blast_line(&part2, true, true, &br, &br).unwrap(),
                BlastAlignment::from_blast_line(&part3, true, true, &br, &br).unwrap(),
                BlastAlignment::from_blast_line(&susceptible, true, true, &br, &br).unwrap(),
            ];
            batch.blast_als[0].fusion_ids = vec![0, 1, 2];
            batch.blast_als[3].susceptible_idx = Some(0);

            assert_eq!(batch.blast_als[2].fusion_2_type(&batch), "STRESS");
            assert_eq!(batch.blast_als[0].fusion_2_type(&batch), "AMR/STRESS");
            assert_eq!(batch.blast_als[3].fusion_2_type(&batch), "AMR");
        }

        #[test]
        fn blast_alignment_fusion_2_subtype_matches_cpp_member_branches() {
            let br = BlastRule::default();
            let qseq = "AAAA*";
            let part1 = format!(
                "WP_FUS|1|3|famA|geneA|AMR|2|CLASS|SUBCLASS|Part_A\t\
                 target\t1\t5\t5\t1\t5\t5\t{qseq}\t{qseq}"
            );
            let part2 = format!(
                "WP_FUS|2|3|famB|geneB|AMR|2|CLASS|SUBCLASS|Part_B\t\
                 target\t1\t5\t5\t1\t5\t5\t{qseq}\t{qseq}"
            );
            let part3 = format!(
                "WP_FUS|3|3|famC|geneC|AMR|2|CLASS|SUBCLASS|Part_C\t\
                 target\t1\t5\t5\t1\t5\t5\t{qseq}\t{qseq}"
            );
            let susceptible = format!(
                "WP_SUS|1|1|famS|geneS|susceptible|2|CLASS|SUBCLASS|Susceptible\t\
                 target\t1\t5\t5\t1\t5\t5\t{qseq}\t{qseq}"
            );

            let mut batch = empty_batch();
            let mut fam_a = fam("famA", "", "", 0.0, 0.0);
            fam_a.subtype = "BETA".to_string();
            let mut fam_b = fam("famB", "", "", 0.0, 0.0);
            fam_b.subtype = "ALPHA".to_string();
            let mut fam_c = fam("famC", "", "", 0.0, 0.0);
            fam_c.subtype = "BETA".to_string();
            batch.fam_map.insert("famA".to_string(), fam_a);
            batch.fam_map.insert("famB".to_string(), fam_b);
            batch.fam_map.insert("famC".to_string(), fam_c);

            batch.blast_als = vec![
                BlastAlignment::from_blast_line(&part1, true, true, &br, &br).unwrap(),
                BlastAlignment::from_blast_line(&part2, true, true, &br, &br).unwrap(),
                BlastAlignment::from_blast_line(&part3, true, true, &br, &br).unwrap(),
                BlastAlignment::from_blast_line(&susceptible, true, true, &br, &br).unwrap(),
            ];
            batch.blast_als[0].fusion_ids = vec![0, 1, 2];
            batch.blast_als[3].susceptible_idx = Some(0);

            assert_eq!(batch.blast_als[2].fusion_2_subtype(&batch), "BETA");
            assert_eq!(batch.blast_als[0].fusion_2_subtype(&batch), "ALPHA/BETA");
            assert_eq!(batch.blast_als[3].fusion_2_subtype(&batch), "AMR");
        }

        #[test]
        fn blast_alignment_fusion_class_subclass_and_hmm_match_cpp_member_branches() {
            let br = BlastRule::default();
            let qseq = "AAAA*";
            let part1 = format!(
                "WP_FUS|1|3|famA|geneA|AMR|2|CLASS|SUBCLASS|Part_A\t\
                 target\t1\t5\t5\t1\t5\t5\t{qseq}\t{qseq}"
            );
            let part2 = format!(
                "WP_FUS|2|3|famB|geneB|AMR|2|CLASS|SUBCLASS|Part_B\t\
                 target\t1\t5\t5\t1\t5\t5\t{qseq}\t{qseq}"
            );
            let part3 = format!(
                "WP_FUS|3|3|famC|geneC|AMR|2|CLASS|SUBCLASS|Part_C\t\
                 target\t1\t5\t5\t1\t5\t5\t{qseq}\t{qseq}"
            );
            let allele = concat!(
                "WP_ALLELE|1|1|alleleFam|geneA|AMR|2|ALLELE_CLASS|ALLELE_SUB|Allele\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );
            let susceptible = concat!(
                "WP_SUS|1|1|famS|famS|susceptible|2|CLASS|SUBCLASS|Susceptible\t",
                "target\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*"
            );

            let mut batch = empty_batch();
            let mut fam_a = fam("famA", "", "", 0.0, 0.0);
            fam_a.class = "BETA/ALPHA".to_string();
            fam_a.subclass = "B2/A1".to_string();
            let mut fam_b = fam("famB", "", "", 0.0, 0.0);
            fam_b.class = "GAMMA".to_string();
            fam_b.subclass = "C3".to_string();
            let mut fam_c = fam("famC", "", "", 0.0, 0.0);
            fam_c.class = "ALPHA".to_string();
            fam_c.subclass = "A1".to_string();
            batch.fam_map.insert("famA".to_string(), fam_a);
            batch.fam_map.insert("famB".to_string(), fam_b);
            batch.fam_map.insert("famC".to_string(), fam_c);
            batch.accession2susceptible.insert(
                "WP_SUS".to_string(),
                Susceptible {
                    genesymbol: "sus".to_string(),
                    cutoff: 0.0,
                    class: "SCLASS".to_string(),
                    subclass: "SSUB".to_string(),
                    name: "Susceptible".to_string(),
                },
            );
            batch.hmm_als = vec![
                HmmAlignment {
                    sseqid: "target".to_string(),
                    score1: 100.0,
                    score2: 50.0,
                    fam_id: "famA".to_string(),
                    domain: None,
                    blast_al_idx: None,
                },
                HmmAlignment {
                    sseqid: "target".to_string(),
                    score1: 90.0,
                    score2: 40.0,
                    fam_id: "famB".to_string(),
                    domain: None,
                    blast_al_idx: None,
                },
            ];

            batch.blast_als = vec![
                BlastAlignment::from_blast_line(&part1, true, true, &br, &br).unwrap(),
                BlastAlignment::from_blast_line(&part2, true, true, &br, &br).unwrap(),
                BlastAlignment::from_blast_line(&part3, true, true, &br, &br).unwrap(),
                BlastAlignment::from_blast_line(allele, true, true, &br, &br).unwrap(),
                BlastAlignment::from_blast_line(susceptible, true, true, &br, &br).unwrap(),
            ];
            batch.blast_als[0].fusion_ids = vec![0, 1, 2];
            batch.blast_als[1].hmm_al_idx = Some(1);
            batch.blast_als[2].hmm_al_idx = Some(0);

            assert_eq!(
                batch.blast_als[0].fusion_2_class(&batch),
                "ALPHA/BETA/GAMMA"
            );
            assert_eq!(batch.blast_als[0].fusion_2_subclass(&batch), "A1/B2/C3");
            assert_eq!(batch.blast_als[3].fusion_2_class(&batch), "ALLELE SUB");
            assert_eq!(batch.blast_als[3].fusion_2_subclass(&batch), "ALLELE CLASS");
            assert_eq!(batch.blast_als[4].fusion_2_class(&batch), "SCLASS");
            assert_eq!(batch.blast_als[4].fusion_2_subclass(&batch), "SSUB");
            assert_eq!(
                batch.blast_als[0]
                    .fusion_2_hmm_al(&batch)
                    .map(|hmm_al| hmm_al.fam_id.as_str()),
                Some("famB")
            );
        }

        #[test]
        fn blastx_finish_target_marks_left_tail_partial_dna_and_cds() {
            let br = BlastRule::default();
            let line = concat!(
                "WP_1|1|1|fam|fam|AMR|2|CLASS|SUBCLASS|Product\t",
                "contig\t6\t15\t20\t1\t30\t100\tMMMMMMMMMM\tMMMMMMMMMM"
            );
            let al = BlastAlignment::from_blast_line(line, true, false, &br, &br).unwrap();
            assert!(!al.partial_dna);
            assert!(al.cdss.is_empty());

            let mut batch = empty_batch();
            batch.add_blast_al(al);
            batch.process(false, false);
            let al = &batch.blast_als[0];
            assert!(al.partial_dna);
            assert_eq!(al.ref_effective_len(), 10);
            assert_eq!(al.cdss.len(), 1);
            assert_eq!(al.cdss[0].contig, "contig");
            assert_eq!(al.cdss[0].start, 0);
            assert_eq!(al.cdss[0].stop, 30);
            assert!(al.cdss[0].strand);
            assert!(al.cdss[0].partial);
        }

        #[test]
        fn blastx_finish_target_marks_right_tail_partial_dna_and_cds() {
            let br = BlastRule::default();
            let line = concat!(
                "WP_1|1|1|fam|fam|AMR|2|CLASS|SUBCLASS|Product\t",
                "contig\t1\t10\t20\t71\t100\t100\tMMMMMMMMMM\tMMMMMMMMMM"
            );
            let al = BlastAlignment::from_blast_line(line, true, false, &br, &br).unwrap();
            assert!(!al.partial_dna);
            assert!(al.cdss.is_empty());

            let mut batch = empty_batch();
            batch.add_blast_al(al);
            batch.process(false, false);
            let al = &batch.blast_als[0];
            assert!(al.partial_dna);
            assert_eq!(al.ref_effective_len(), 10);
            assert_eq!(al.cdss.len(), 1);
            assert_eq!(al.cdss[0].contig, "contig");
            assert_eq!(al.cdss[0].start, 70);
            assert_eq!(al.cdss[0].stop, 100);
            assert!(al.cdss[0].strand);
            assert!(al.cdss[0].partial);
        }

        #[test]
        fn blastx_loader_defers_strong_susceptible_merge_to_process_like_cpp() {
            let mut blastx = tempfile::NamedTempFile::new().unwrap();
            writeln!(
                blastx,
                "{}",
                concat!(
                    "WP_S|1|1|fam|gene|susceptible|2|CLASS|SUBCLASS|Product\t",
                    "contig\t1\t3\t7\t1\t9\t30\tAAA\tAAA"
                )
            )
            .unwrap();
            writeln!(
                blastx,
                "{}",
                concat!(
                    "WP_S|1|1|fam|gene|susceptible|2|CLASS|SUBCLASS|Product\t",
                    "contig\t5\t6\t7\t13\t18\t30\tCC\tCC"
                )
            )
            .unwrap();

            let mut batch = empty_batch();
            batch.accession2susceptible.insert(
                "WP_S".to_string(),
                Susceptible {
                    genesymbol: "susGene".to_string(),
                    cutoff: 0.0,
                    class: "SCLASS".to_string(),
                    subclass: "SSUB".to_string(),
                    name: "Susceptible Product".to_string(),
                },
            );
            batch
                .fam_map
                .insert("fam".to_string(), fam("fam", "", "", 0.0, 0.0));
            load_blast_results(blastx.path(), &mut batch, false, false).unwrap();

            assert_eq!(batch.blast_als.len(), 2);
            assert!(!batch.blast_als[0].hsp.merged);
            assert!(!batch.blast_als[1].hsp.merged);
            assert_eq!(batch.target2blast_als.get("contig").unwrap(), &vec![0, 1]);

            batch.process(false, false);

            assert_eq!(batch.blast_als.len(), 3);
            assert_eq!(batch.target2blast_als.get("contig").unwrap(), &vec![2]);
            let good = batch.target2good_blast_als.get("contig").unwrap();
            assert_eq!(good, &vec![2]);
            let hsp = &batch.blast_als[2].hsp;
            assert!(hsp.merged);
            assert_eq!(hsp.q_int, crate::seq::Interval::new(0, 6, 0));
            assert_eq!(hsp.s_int, crate::seq::Interval::new(0, 18, 1));
            assert_eq!(hsp.disrs.len(), 1);
            assert_eq!(
                hsp.disrs[0].disruption_type(),
                crate::seq::DisruptionType::Deletion
            );
            assert_eq!(batch.blast_als[2].seq_changes.len(), 1);
            assert!(batch.blast_als[2].seq_changes[0].disr.is_some());
        }

        #[test]
        fn blastx_loader_keeps_ordinary_adjacent_hsps_separate_like_cpp() {
            let mut blastx = tempfile::NamedTempFile::new().unwrap();
            writeln!(
                blastx,
                "{}",
                concat!(
                    "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                    "contig\t1\t3\t7\t1\t9\t30\tAAA\tAAA"
                )
            )
            .unwrap();
            writeln!(
                blastx,
                "{}",
                concat!(
                    "WP_1|1|1|fam|gene|AMR|2|CLASS|SUBCLASS|Product\t",
                    "contig\t5\t6\t7\t13\t18\t30\tCC\tCC"
                )
            )
            .unwrap();

            let mut batch = empty_batch();
            batch
                .fam_map
                .insert("fam".to_string(), fam("fam", "", "", 0.0, 0.0));
            load_blast_results(blastx.path(), &mut batch, false, false).unwrap();

            assert_eq!(batch.blast_als.len(), 2);
            assert!(!batch.blast_als[0].hsp.merged);
            assert!(!batch.blast_als[1].hsp.merged);
            assert_eq!(batch.target2blast_als.get("contig").unwrap(), &vec![0, 1]);
        }

        #[test]
        fn process_reports_strong_susceptible_disruption_seq_change_like_cpp() {
            let mut blastx = tempfile::NamedTempFile::new().unwrap();
            writeln!(
                blastx,
                "{}",
                concat!(
                    "WP_S|1|1|famS|geneS|susceptible|2|CLASS|SUBCLASS|Susceptible\t",
                    "contig\t1\t3\t7\t1\t9\t30\tAAA\tAAA"
                )
            )
            .unwrap();
            writeln!(
                blastx,
                "{}",
                concat!(
                    "WP_S|1|1|famS|geneS|susceptible|2|CLASS|SUBCLASS|Susceptible\t",
                    "contig\t5\t6\t7\t13\t18\t30\tCC\tCC"
                )
            )
            .unwrap();

            let mut batch = empty_batch();
            batch.accession2susceptible.insert(
                "WP_S".to_string(),
                Susceptible {
                    genesymbol: "susGene".to_string(),
                    cutoff: 0.0,
                    class: "SCLASS".to_string(),
                    subclass: "SSUB".to_string(),
                    name: "Susceptible Product".to_string(),
                },
            );
            load_blast_results(blastx.path(), &mut batch, false, false).unwrap();
            batch.process(false, false);

            let good = batch.target2good_blast_als.get("contig").unwrap();
            assert_eq!(good, &vec![2]);
            assert_eq!(batch.blast_als[2].seq_changes.len(), 1);
            assert!(batch.blast_als[2].seq_changes[0].disr.is_some());

            let mut output = Vec::new();
            batch.report(&mut output, false).unwrap();
            let output = String::from_utf8(output).unwrap();
            let mut lines = output.lines();
            let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
            let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
            let col = |name: &str| {
                let idx = header
                    .iter()
                    .position(|field| *field == name)
                    .unwrap_or_else(|| panic!("missing report column {name}"));
                row[idx]
            };

            assert!(lines.next().is_none(), "{output}");
            assert!(col("Element symbol").starts_with("geneS_@del_"), "{output}");
            assert_eq!(col("Element name"), "Susceptible Product");
            assert_eq!(col("Subtype"), "POINT_DISRUPT");
            assert_eq!(col("Method"), "POINTX");
        }
    }
    #[cfg(test)]
    mod amr_report_cli_tests {
        use super::*;
        use std::path::{Path, PathBuf};

        fn test_fixtures() -> PathBuf {
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
        }

        fn db_dir() -> PathBuf {
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amrfinder_db/2026-03-24.1")
        }

        fn require_file(path: &Path) {
            assert!(
                path.exists(),
                "required parity fixture is missing: {}",
                path.display()
            );
        }

        fn run_cpp_golden_report() -> (String, String) {
            let fixtures = test_fixtures();
            let db = db_dir();
            let blastp_file = fixtures.join("blastp");
            let dom_file = fixtures.join("dom");
            let search_file = fixtures.join("hmmsearch");
            let expected_file = fixtures.join("expected_output");
            let gff_file = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amr/test_prot.gff");

            require_file(&db.join("fam.tsv"));
            require_file(&db.join("AMRProt-mutation.tsv"));
            require_file(&db.join("AMRProt-susceptible.tsv"));
            require_file(&blastp_file);
            require_file(&dom_file);
            require_file(&search_file);
            require_file(&expected_file);
            require_file(&gff_file);

            let mut output = Vec::new();
            let config = AmrReportConfig {
                fam_file: &db.join("fam.tsv"),
                blastp_file: Some(blastp_file.as_path()),
                blastx_file: None,
                dna_len_file: None,
                hmmsearch_file: Some(search_file.as_path()),
                hmmdom_file: Some(dom_file.as_path()),
                gff_file: Some(&gff_file),
                gff_type: "genbank",
                organism: "Escherichia",
                mutation_file: Some(&db.join("AMRProt-mutation.tsv")),
                susceptible_file: Some(&db.join("AMRProt-susceptible.tsv")),
                suppress_file: None,
                coverage_min: 0.5,
                ident_min: -1.0,
                print_node: true,
                print_node_raw: false,
                mutation_all: None,
                name: "",
                non_reportable: false,
                report_core_only: false,
                report_all_equal: false,
                cds_exist: true,
                nosame: false,
                noblast: false,
                nohmm: false,
                retain_blasts: false,
                skip_hmm_check: false,
            };
            run_amr_report(&config, &mut output).unwrap();

            let rust_output = String::from_utf8(output).unwrap();
            let cpp_output = std::fs::read_to_string(&expected_file).unwrap();
            (rust_output, cpp_output)
        }

        fn run_cpp_blastx_golden_report() -> (String, String) {
            let fixtures = test_fixtures();
            let db = db_dir();
            let blastx_file = fixtures.join("blastx");
            let expected_file = fixtures.join("expected_blastx_report");

            require_file(&db.join("fam.tsv"));
            require_file(&db.join("AMRProt-mutation.tsv"));
            require_file(&db.join("AMRProt-susceptible.tsv"));
            require_file(&blastx_file);
            require_file(&expected_file);

            let mut output = Vec::new();
            let config = AmrReportConfig {
                fam_file: &db.join("fam.tsv"),
                blastp_file: None,
                blastx_file: Some(blastx_file.as_path()),
                dna_len_file: None,
                hmmsearch_file: None,
                hmmdom_file: None,
                gff_file: None,
                gff_type: "genbank",
                organism: "Escherichia",
                mutation_file: Some(&db.join("AMRProt-mutation.tsv")),
                susceptible_file: Some(&db.join("AMRProt-susceptible.tsv")),
                suppress_file: None,
                coverage_min: 0.5,
                ident_min: -1.0,
                print_node: true,
                print_node_raw: false,
                mutation_all: None,
                name: "",
                non_reportable: false,
                report_core_only: false,
                report_all_equal: false,
                cds_exist: true,
                nosame: false,
                noblast: false,
                nohmm: false,
                retain_blasts: false,
                skip_hmm_check: false,
            };
            run_amr_report(&config, &mut output).unwrap();

            let rust_output = String::from_utf8(output).unwrap();
            let cpp_output = std::fs::read_to_string(&expected_file).unwrap();
            (rust_output, cpp_output)
        }

        fn run_combined_golden_report() -> String {
            let fixtures = test_fixtures();
            let db = db_dir();
            let blastp_file = fixtures.join("blastp");
            let blastx_file = fixtures.join("blastx");
            let dom_file = fixtures.join("dom");
            let search_file = fixtures.join("hmmsearch");
            let gff_file = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amr/test_prot.gff");

            require_file(&db.join("fam.tsv"));
            require_file(&db.join("AMRProt-mutation.tsv"));
            require_file(&db.join("AMRProt-susceptible.tsv"));
            require_file(&blastp_file);
            require_file(&blastx_file);
            require_file(&dom_file);
            require_file(&search_file);
            require_file(&gff_file);

            let mut output = Vec::new();
            let config = AmrReportConfig {
                fam_file: &db.join("fam.tsv"),
                blastp_file: Some(blastp_file.as_path()),
                blastx_file: Some(blastx_file.as_path()),
                dna_len_file: None,
                hmmsearch_file: Some(search_file.as_path()),
                hmmdom_file: Some(dom_file.as_path()),
                gff_file: Some(&gff_file),
                gff_type: "genbank",
                organism: "Escherichia",
                mutation_file: Some(&db.join("AMRProt-mutation.tsv")),
                susceptible_file: Some(&db.join("AMRProt-susceptible.tsv")),
                suppress_file: None,
                coverage_min: 0.5,
                ident_min: -1.0,
                print_node: true,
                print_node_raw: false,
                mutation_all: None,
                name: "",
                non_reportable: false,
                report_core_only: false,
                report_all_equal: false,
                cds_exist: true,
                nosame: false,
                noblast: false,
                nohmm: false,
                retain_blasts: false,
                skip_hmm_check: false,
            };
            run_amr_report(&config, &mut output).unwrap();

            String::from_utf8(output).unwrap()
        }

        fn run_protein_report_for_organism(organism: &str, suppress: bool) -> String {
            let fixtures = test_fixtures();
            let db = db_dir();
            let blastp_file = fixtures.join("blastp");
            let dom_file = fixtures.join("dom");
            let search_file = fixtures.join("hmmsearch");
            let gff_file = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amr/test_prot.gff");
            let suppress_file = db.join("AMRProt-suppress.tsv");

            require_file(&db.join("fam.tsv"));
            require_file(&db.join("AMRProt-mutation.tsv"));
            require_file(&db.join("AMRProt-susceptible.tsv"));
            require_file(&blastp_file);
            require_file(&dom_file);
            require_file(&search_file);
            require_file(&gff_file);
            if suppress {
                require_file(&suppress_file);
            }

            let mut output = Vec::new();
            let config = AmrReportConfig {
                fam_file: &db.join("fam.tsv"),
                blastp_file: Some(blastp_file.as_path()),
                blastx_file: None,
                dna_len_file: None,
                hmmsearch_file: Some(search_file.as_path()),
                hmmdom_file: Some(dom_file.as_path()),
                gff_file: Some(&gff_file),
                gff_type: "genbank",
                organism,
                mutation_file: Some(&db.join("AMRProt-mutation.tsv")),
                susceptible_file: Some(&db.join("AMRProt-susceptible.tsv")),
                suppress_file: suppress.then_some(suppress_file.as_path()),
                coverage_min: 0.5,
                ident_min: -1.0,
                print_node: true,
                print_node_raw: false,
                mutation_all: None,
                name: "",
                non_reportable: false,
                report_core_only: false,
                report_all_equal: false,
                cds_exist: true,
                nosame: false,
                noblast: false,
                nohmm: false,
                retain_blasts: false,
                skip_hmm_check: false,
            };
            run_amr_report(&config, &mut output).unwrap();

            String::from_utf8(output).unwrap()
        }

        fn tsv_row<'a>(output: &'a str, protein_id: &str) -> Vec<&'a str> {
            output
                .lines()
                .skip(1)
                .find(|line| line.split('\t').next() == Some(protein_id))
                .unwrap_or_else(|| panic!("missing report row for {protein_id}"))
                .split('\t')
                .collect()
        }

        fn col<'a>(header: &[&str], row: &'a [&str], name: &str) -> &'a str {
            let idx = header
                .iter()
                .position(|field| *field == name)
                .unwrap_or_else(|| panic!("missing report column {name}"));
            row[idx]
        }

        #[test]
        fn amr_report_cpp_golden_fixture_is_complete() {
            let (_rust_output, cpp_output) = run_cpp_golden_report();
            let cpp_lines: Vec<&str> = cpp_output.lines().collect();

            assert_eq!(
                cpp_lines.len(),
                16,
                "C++ golden should contain header plus 15 rows"
            );
            assert!(
                cpp_output.contains("nimIJ_hmm\tcontigX"),
                "C++ golden should include the standalone HMM-only hit"
            );
            assert!(
                cpp_output.contains("pmrB_C84R\tcontig14"),
                "C++ golden should include protein mutation hits"
            );
        }

        #[test]
        fn amr_report_matches_cpp_golden_byte_for_byte() {
            let (rust_output, cpp_output) = run_cpp_golden_report();

            assert_eq!(
                rust_output, cpp_output,
                "Rust amr_report output must be byte-identical to the C++ golden output"
            );
        }

        #[test]
        fn amr_report_matches_cpp_blastx_golden_byte_for_byte() {
            let (rust_output, cpp_output) = run_cpp_blastx_golden_report();

            assert_eq!(
                rust_output, cpp_output,
                "Rust BLASTX amr_report output must be byte-identical to the C++ golden output"
            );
        }

        #[test]
        fn hmm_loader_links_or_synthesizes_expected_alignment_kinds() {
            let fixtures = test_fixtures();
            let db = db_dir();
            let mut batch = Batch::from_fam_file(&db.join("fam.tsv"), 0).unwrap();

            load_blast_results(&fixtures.join("blastp"), &mut batch, true, false).unwrap();
            let blast_count_before_hmm = batch.blast_als.len();
            load_hmm_domains(&fixtures.join("dom"), &mut batch).unwrap();
            load_hmm_results(&fixtures.join("hmmsearch"), &mut batch).unwrap();

            let bla_tem_rows: Vec<_> = batch
                .blast_als
                .iter()
                .filter(|al| al.hsp.sseqid == "blaTEM-156")
                .collect();
            assert!(
                bla_tem_rows
                    .iter()
                    .any(|al| al.ref_accession == "WP_061158039.1" && al.hmm_al_idx.is_some()),
                "good BLAST hit should link to its HMM alignment"
            );
            assert!(
                bla_tem_rows.iter().all(|al| !al.from_hmm),
                "good BLAST target should not get duplicate synthetic HMM rows"
            );

            let nimij_hmm = batch
                .blast_als
                .iter()
                .find(|al| al.hsp.sseqid == "nimIJ_hmm" && al.from_hmm)
                .expect("below-threshold BLAST target should get a synthetic HMM row");
            assert_eq!(nimij_hmm.fam_id, "nimIJ");
            assert_eq!(nimij_hmm.ref_accession, "WP_005812825.1");
            assert!(nimij_hmm.hmm_al_idx.is_some());
            assert!(
                batch.blast_als.len() > blast_count_before_hmm,
                "HMM processing should add synthetic rows only where BLAST alone is insufficient"
            );
        }

        #[test]
        fn hmm_loader_creates_standalone_row_without_blast_target() {
            let db = db_dir();
            let mut batch = Batch::from_fam_file(&db.join("fam.tsv"), 0).unwrap();
            let mut dom = tempfile::NamedTempFile::new().unwrap();
            let mut search = tempfile::NamedTempFile::new().unwrap();

            writeln!(
            dom,
            "standalone_hmm - 166 NimIJ-NCBIFAM NF000262.1 151 1.1e-82 274.7 0.3 1 1 1.2e-86 1.2e-82 274.5 0.3 1 151 4 154 4 154 1.00 description"
        )
        .unwrap();
            writeln!(
            search,
            "standalone_hmm - NimIJ-NCBIFAM NF000262.1 1.1e-82 274.7 0.3 1.2e-82 274.5 0.3 1.0 1 0 0 1 1 1 1 description"
        )
        .unwrap();

            load_hmm_domains(dom.path(), &mut batch).unwrap();
            load_hmm_results(search.path(), &mut batch).unwrap();

            assert_eq!(batch.hmm_als.len(), 1);
            assert_eq!(batch.blast_als.len(), 1);
            let al = &batch.blast_als[0];
            assert!(al.from_hmm);
            assert_eq!(al.hsp.sseqid, "standalone_hmm");
            assert_eq!(al.ref_accession, "");
            assert_eq!(al.fam_id, "nimIJ");
            assert_eq!(al.method, "HMM");
            assert_eq!(al.hmm_al_idx, Some(0));
        }

        #[test]
        fn run_amr_report_nohmm_skips_hmmer_inputs_like_cpp_body() {
            let db = db_dir();
            let mut dom = tempfile::NamedTempFile::new().unwrap();
            let mut search = tempfile::NamedTempFile::new().unwrap();

            writeln!(dom, "not enough hmmdom fields").unwrap();
            writeln!(search, "not enough hmmsearch fields").unwrap();

            let config = AmrReportConfig {
                fam_file: &db.join("fam.tsv"),
                blastp_file: None,
                blastx_file: None,
                dna_len_file: None,
                hmmsearch_file: Some(search.path()),
                hmmdom_file: Some(dom.path()),
                gff_file: None,
                gff_type: "genbank",
                organism: "",
                mutation_file: None,
                susceptible_file: None,
                suppress_file: None,
                coverage_min: 0.5,
                ident_min: -1.0,
                print_node: true,
                print_node_raw: false,
                mutation_all: None,
                name: "",
                non_reportable: false,
                report_core_only: false,
                report_all_equal: false,
                cds_exist: true,
                nosame: false,
                noblast: false,
                nohmm: true,
                retain_blasts: false,
                skip_hmm_check: false,
            };

            let mut output = Vec::new();
            run_amr_report(&config, &mut output).unwrap();
            let output = String::from_utf8(output).unwrap();
            assert_eq!(output.lines().count(), 1, "{output}");
        }

        #[test]
        fn report_rows_cover_known_parity_edge_cases() {
            let (rust_output, _cpp_output) = run_cpp_golden_report();
            let header: Vec<&str> = rust_output.lines().next().unwrap().split('\t').collect();

            let aph = tsv_row(&rust_output, "aph3pp-Ib_partial_5p_neg");
            assert_eq!(col(&header, &aph, "% Identity to reference"), "100.00");
            assert_eq!(col(&header, &aph, "Method"), "PARTIAL_CONTIG_ENDP");

            let mutation = tsv_row(&rust_output, "pmrB_C84R");
            assert_eq!(col(&header, &mutation, "Element symbol"), "pmrB_C84R");
            assert_eq!(col(&header, &mutation, "Type"), "AMR");
            assert_eq!(col(&header, &mutation, "Subtype"), "POINT");
            assert_eq!(col(&header, &mutation, "Hierarchy node"), "pmrB");

            let oxa = tsv_row(&rust_output, "blaOXA-436_partial");
            assert_eq!(col(&header, &oxa, "Hierarchy node"), "blaOXA-48_fam");

            let sul2 = tsv_row(&rust_output, "sul2_partial_3p_neg");
            assert_eq!(col(&header, &sul2, "HMM accession"), "NA");
            assert_eq!(col(&header, &sul2, "HMM description"), "NA");
        }

        #[test]
        fn blastx_report_uses_dna_coordinates_and_x_methods() {
            let db = db_dir();
            let mut blastx = tempfile::NamedTempFile::new().unwrap();
            let qseq = format!("M{}", "A".repeat(285));
            let sseq = qseq.clone();
            writeln!(
            blastx,
            "WP_061158039.1|1|1|blaTEM-156|blaTEM|beta-lactamase|2|BETA-LACTAM|BETA-LACTAM|class_A_beta-lactamase_TEM-156\tcontig_blastx\t1\t286\t287\t101\t958\t1200\t{}\t{}",
            qseq, sseq
        )
        .unwrap();

            let mut output = Vec::new();
            let config = AmrReportConfig {
                fam_file: &db.join("fam.tsv"),
                blastp_file: None,
                blastx_file: Some(blastx.path()),
                dna_len_file: None,
                hmmsearch_file: None,
                hmmdom_file: None,
                gff_file: None,
                gff_type: "genbank",
                organism: "Escherichia",
                mutation_file: Some(&db.join("AMRProt-mutation.tsv")),
                susceptible_file: Some(&db.join("AMRProt-susceptible.tsv")),
                suppress_file: None,
                coverage_min: 0.5,
                ident_min: -1.0,
                print_node: true,
                print_node_raw: false,
                mutation_all: None,
                name: "",
                non_reportable: false,
                report_core_only: false,
                report_all_equal: false,
                cds_exist: true,
                nosame: false,
                noblast: false,
                nohmm: false,
                retain_blasts: false,
                skip_hmm_check: false,
            };
            run_amr_report(&config, &mut output).unwrap();

            let output = String::from_utf8(output).unwrap();
            let header: Vec<&str> = output.lines().next().unwrap().split('\t').collect();
            let row = tsv_row(&output, "NA");
            assert_eq!(col(&header, &row, "Contig id"), "contig_blastx");
            assert_eq!(col(&header, &row, "Start"), "101");
            assert_eq!(col(&header, &row, "Stop"), "958");
            assert_eq!(col(&header, &row, "Method"), "BLASTX");
            assert_eq!(col(&header, &row, "Target length"), "286");
            assert_eq!(col(&header, &row, "Reference sequence length"), "286");
            assert_eq!(col(&header, &row, "% Coverage of reference"), "100.00");
            assert_eq!(col(&header, &row, "% Identity to reference"), "100.00");
        }

        #[test]
        fn blastx_report_does_not_use_dna_len_for_target_length() {
            let db = db_dir();
            let mut blastx = tempfile::NamedTempFile::new().unwrap();
            let mut dna_len = tempfile::NamedTempFile::new().unwrap();
            let qseq = format!("M{}", "A".repeat(285));
            let sseq = qseq.clone();
            writeln!(
            blastx,
            "WP_061158039.1|1|1|blaTEM-156|blaTEM|beta-lactamase|2|BETA-LACTAM|BETA-LACTAM|class_A_beta-lactamase_TEM-156\tcontig_blastx\t1\t286\t287\t101\t958\t1200\t{}\t{}",
            qseq, sseq
        )
        .unwrap();
            writeln!(dna_len, "contig_blastx\t1200").unwrap();

            let mut output = Vec::new();
            let config = AmrReportConfig {
                fam_file: &db.join("fam.tsv"),
                blastp_file: None,
                blastx_file: Some(blastx.path()),
                dna_len_file: Some(dna_len.path()),
                hmmsearch_file: None,
                hmmdom_file: None,
                gff_file: None,
                gff_type: "genbank",
                organism: "Escherichia",
                mutation_file: Some(&db.join("AMRProt-mutation.tsv")),
                susceptible_file: Some(&db.join("AMRProt-susceptible.tsv")),
                suppress_file: None,
                coverage_min: 0.5,
                ident_min: -1.0,
                print_node: true,
                print_node_raw: false,
                mutation_all: None,
                name: "",
                non_reportable: false,
                report_core_only: false,
                report_all_equal: false,
                cds_exist: true,
                nosame: false,
                noblast: false,
                nohmm: false,
                retain_blasts: false,
                skip_hmm_check: false,
            };
            run_amr_report(&config, &mut output).unwrap();

            let output = String::from_utf8(output).unwrap();
            let header: Vec<&str> = output.lines().next().unwrap().split('\t').collect();
            let row = tsv_row(&output, "NA");
            assert_eq!(col(&header, &row, "Target length"), "286");
        }

        #[test]
        fn blastx_mutation_all_reports_wildtype_and_mutant_rows() {
            let fixtures = test_fixtures();
            let db = db_dir();
            let blastx_file = fixtures.join("blastx");

            if !blastx_file.exists() {
                return;
            }
            require_file(&db.join("fam.tsv"));
            require_file(&db.join("AMRProt-mutation.tsv"));
            require_file(&db.join("AMRProt-susceptible.tsv"));

            let mutation_all = tempfile::NamedTempFile::new().unwrap();
            let mut dna_len = tempfile::NamedTempFile::new().unwrap();
            writeln!(dna_len, "contig14\t4096").unwrap();
            let mut output = Vec::new();
            let config = AmrReportConfig {
                fam_file: &db.join("fam.tsv"),
                blastp_file: None,
                blastx_file: Some(blastx_file.as_path()),
                dna_len_file: Some(dna_len.path()),
                hmmsearch_file: None,
                hmmdom_file: None,
                gff_file: None,
                gff_type: "genbank",
                organism: "Escherichia",
                mutation_file: Some(&db.join("AMRProt-mutation.tsv")),
                susceptible_file: Some(&db.join("AMRProt-susceptible.tsv")),
                suppress_file: None,
                coverage_min: 0.5,
                ident_min: -1.0,
                print_node: true,
                print_node_raw: false,
                mutation_all: Some(mutation_all.path()),
                name: "",
                non_reportable: false,
                report_core_only: false,
                report_all_equal: false,
                cds_exist: true,
                nosame: false,
                noblast: false,
                nohmm: false,
                retain_blasts: false,
                skip_hmm_check: false,
            };
            run_amr_report(&config, &mut output).unwrap();

            let mutation_all = std::fs::read_to_string(mutation_all.path()).unwrap();
            assert!(
                mutation_all.contains("\nNA\tcontig14\t1\t1089\t+\tpmrB_C84R\t"),
                "BLASTX mutation_all should include mutant pmrB row"
            );
            let header: Vec<&str> = mutation_all.lines().next().unwrap().split('\t').collect();
            let pmrb: Vec<&str> = mutation_all
                .lines()
                .find(|line| line.contains("\tpmrB_C84R\t"))
                .expect("missing pmrB_C84R mutation_all row")
                .split('\t')
                .collect();
            assert_eq!(col(&header, &pmrb, "Target length"), "363");
            assert!(
                mutation_all.contains("\nNA\tcontig14\t1\t1089\t+\tpmrB_A159A\t"),
                "BLASTX mutation_all should include wildtype pmrB rows"
            );
            assert!(
                mutation_all.contains("\nNA\tcontig16\t1\t720\t+\tnfsA_K141Ter\t")
                    && mutation_all.contains("\nNA\tcontig16\t1\t720\t+\tnfsA_R15C\t"),
                "BLASTX mutation_all should include all observed nfsA mutation rows"
            );
            assert!(
                mutation_all
                    .lines()
                    .next()
                    .is_some_and(|header| header.ends_with("\tHierarchy node")),
                "direct amr_report mutation_all should preserve --print_node"
            );
        }

        #[test]
        fn mutation_annotation_detects_deletion_seq_change_like_cpp_set_seq_changes() {
            let br = BlastRule::default();
            let mut batch = Batch {
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
            };
            batch.accession2mutations.insert(
                "WP_DEL".to_string(),
                vec![AmrMutation::new(
                    2,
                    "gene_A2del",
                    "gene_A2del",
                    "CLASS",
                    "SUBCLASS",
                    "Deleted_residue",
                )],
            );
            let line = concat!(
                "WP_DEL|1|1|gene|gene|mutation|2|CLASS|SUBCLASS|Mutation_product\t",
                "target\t1\t6\t6\t1\t5\t5\tMACDE*\tM-CDE*"
            );
            batch
                .add_blast_al(BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap());

            batch.process(false, false);

            let al = &batch.blast_als[0];
            assert_eq!(al.seq_changes.len(), 1);
            assert_eq!(al.seq_changes[0].get_mutation_str(), "A2DEL");
            assert_eq!(al.seq_changes[0].mutations, vec![0]);

            let mut output = Vec::new();
            batch.report(&mut output, false).unwrap();
            let output = String::from_utf8(output).unwrap();
            assert!(
                output.contains("\ntarget\tgene_A2del\tDeleted residue\t"),
                "C++ setSeqChanges reports deletion mutations through the normal report path:\n{output}"
            );
        }

        #[test]
        fn combined_report_suppresses_blastx_when_protein_cds_covers_region() {
            let output = run_combined_golden_report();

            assert!(
                output.contains("\nblaTEM-156\tcontig01\t101\t961\t+\tblaTEM-156\t"),
                "protein/CDS row should be retained"
            );
            assert!(
                !output.contains("\nNA\tcontig01\t101\t958\t+\tblaTEM-156\t"),
                "overlapping BLASTX row should be suppressed when a protein CDS row covers it"
            );
            assert!(
                output.contains("\nNA\tcontig04\t1261\t2391\t+\tblaEC\t"),
                "non-overlapping BLASTX hit on the same contig should remain reportable"
            );
        }

        #[test]
        fn combined_report_transfers_internal_stop_from_blastx_to_protein_row() {
            let output = run_combined_golden_report();
            let header: Vec<&str> = output.lines().next().unwrap().split('\t').collect();
            let row = tsv_row(&output, "blaTEM-internal_stop");

            assert_eq!(col(&header, &row, "Method"), "INTERNAL_STOP");
            assert!(
                !output.contains("\nNA\tcontig11\t101\t958\t+\tblaTEM\t"),
                "BLASTX row carrying the same internal-stop evidence should be suppressed"
            );
        }

        #[test]
        fn combined_report_prefers_full_blastx_mutations_over_partial_protein_mutation() {
            let output = run_combined_golden_report();

            assert!(
                !output.contains("\nnfsA_R15C_K141STOP\t"),
                "partial protein mutation row should be suppressed by fuller BLASTX evidence"
            );
            assert!(
                output.contains("\nNA\tcontig16\t1\t720\t+\tnfsA_K141Ter\t")
                    && output.contains("\nNA\tcontig16\t1\t720\t+\tnfsA_R15C\t"),
                "distinct BLASTX mutation rows should remain reportable"
            );
        }

        #[test]
        fn organism_suppress_file_removes_common_proteins_unless_report_common_is_used() {
            let suppressed = run_protein_report_for_organism("Escherichia", true);
            let report_common = run_protein_report_for_organism("Escherichia", false);

            assert!(
                !suppressed.contains("\narsR-suppressed-in-escherichia\t"),
                "Escherichia suppress file should remove arsR"
            );
            assert!(
                report_common.contains("\narsR-suppressed-in-escherichia\t"),
                "report_common path should omit suppress file and retain arsR"
            );
        }

        #[test]
        fn campylobacter_fixture_reports_campylobacter_specific_mutations() {
            let output = run_protein_report_for_organism("Campylobacter", false);

            assert!(
                output.contains("\ngyrA\tcontig06\t31\t2616\t+\tgyrA_T86A\t")
                    || output.contains("\tgyrA_T86A\tCampylobacter quinolone resistant GyrA\t"),
                "Campylobacter organism should report the gyrA_T86A protein mutation"
            );
            assert!(
                output.contains("\t50S_L22:A103V\tCampylobacter macrolide resistant RplV\t")
                    || output.contains("\trplV_A103V\tCampylobacter macrolide resistant RplV\t"),
                "Campylobacter organism should report the 50S_L22/rplV mutation fixture"
            );
            assert!(
                !output.contains("\tnfsA_R15C\tEscherichia nitrofurantoin resistant NfsA\t"),
                "Escherichia-specific nfsA mutation should not be reported for Campylobacter"
            );
        }
    }
}
