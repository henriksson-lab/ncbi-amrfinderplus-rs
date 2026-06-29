mod cli {
    include!("../src/main.rs");

    #[cfg(test)]
    mod tests {
        use super::*;
        use std::process::Command;

        #[test]
        fn update_cli_accepts_index_tool_directories() {
            let cli = Cli::try_parse_from([
                "amrfinder",
                "update",
                "--force",
                "--database",
                "/db",
                "--blast_bin",
                "/opt/blast/bin",
                "--hmmer_bin",
                "/opt/hmmer/bin",
                "-q",
                "--debug",
            ])
            .unwrap();

            match cli.command {
                Commands::Update {
                    force,
                    database,
                    blast_bin,
                    hmmer_bin,
                    quiet,
                    debug,
                    top_level_update,
                    ..
                } => {
                    assert!(force);
                    assert_eq!(database, Some(PathBuf::from("/db")));
                    assert_eq!(blast_bin, Some(PathBuf::from("/opt/blast/bin")));
                    assert_eq!(hmmer_bin, Some(PathBuf::from("/opt/hmmer/bin")));
                    assert!(quiet);
                    assert!(debug);
                    assert!(!top_level_update);
                }
                _ => panic!("expected update command"),
            }
        }

        #[test]
        fn update_debug_forwarding_preserves_original_debug_flag() {
            let source = include_str!("../src/main.rs");
            assert!(source.contains("argv.push(OsString::from(\"--debug\"));"));
            assert!(!source.contains("if debug {\n                argv.push(OsString::from(\"--qc\"));\n            }\n            amrfinder::amrfinder_update::main(argv)"));
        }

        #[test]
        fn top_level_short_version_does_not_route_to_run_parser() {
            let output = Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .arg("-v")
                .output()
                .unwrap();

            assert!(output.status.success());
            assert_eq!(String::from_utf8(output.stderr).unwrap(), "");
            assert!(String::from_utf8(output.stdout)
                .unwrap()
                .contains(env!("CARGO_PKG_VERSION")));
        }

        #[test]
        fn upstream_update_flags_route_to_default_update_command() {
            let cli = Cli::parse_from(args_with_default_run_subcommand_from([
                "amrfinder",
                "-update",
            ]));
            match cli.command {
                Commands::Update {
                    force,
                    database,
                    top_level_update,
                    ..
                } => {
                    assert!(!force);
                    assert_eq!(database, None);
                    assert!(top_level_update);
                }
                _ => panic!("expected update command"),
            }

            let cli = Cli::parse_from(args_with_default_run_subcommand_from([
                "amrfinder",
                "-force_update",
            ]));
            match cli.command {
                Commands::Update {
                    force,
                    database,
                    top_level_update,
                    ..
                } => {
                    assert!(force);
                    assert_eq!(database, None);
                    assert!(top_level_update);
                }
                _ => panic!("expected update command"),
            }
        }

        #[test]
        fn upstream_update_flags_route_independently_of_position() {
            let cli = Cli::parse_from(args_with_default_run_subcommand_from([
                "amrfinder",
                "-p",
                "input.fa",
                "-force_update",
            ]));
            match cli.command {
                Commands::Update {
                    force,
                    protein,
                    top_level_update,
                    ..
                } => {
                    assert!(force);
                    assert_eq!(protein, Some(PathBuf::from("input.fa")));
                    assert!(top_level_update);
                }
                _ => panic!("expected update command"),
            }
        }

        #[test]
        fn top_level_update_normalizes_original_long_option_spellings() {
            let cli = Cli::parse_from(args_with_default_run_subcommand_from([
                "amrfinder",
                "-protein=input.fa",
                "-database=/db",
                "-blast_bin=/blast",
                "-hmmer_bin",
                "/hmmer",
                "-qc",
                "-u",
            ]));
            match cli.command {
                Commands::Update {
                    protein,
                    database,
                    blast_bin,
                    hmmer_bin,
                    qc,
                    top_level_update,
                    ..
                } => {
                    assert_eq!(protein, Some(PathBuf::from("input.fa")));
                    assert_eq!(database, Some(PathBuf::from("/db")));
                    assert_eq!(blast_bin, Some(PathBuf::from("/blast")));
                    assert_eq!(hmmer_bin, Some(PathBuf::from("/hmmer")));
                    assert!(qc);
                    assert!(top_level_update);
                }
                _ => panic!("expected update command"),
            }
        }

        #[test]
        fn top_level_update_equal_form_sequence_input_gets_upstream_error() {
            let output = Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .args(["-protein=input.fa", "-u"])
                .output()
                .unwrap();

            assert!(!output.status.success());
            let stderr = String::from_utf8(output.stderr).unwrap();
            assert!(stderr.contains(
                "AMRFinder -u/--update option cannot be run with -n/--nucleotide or -p/--protein options"
            ));
            assert!(!stderr.contains("unexpected argument"));
        }

        #[test]
        fn top_level_update_equal_form_database_gets_upstream_error() {
            let output = Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .args(["-database=/db", "-u"])
                .output()
                .unwrap();

            assert!(!output.status.success());
            let stderr = String::from_utf8(output.stderr).unwrap();
            assert!(stderr.contains(
                "AMRFinder update option (-u/--update) only operates on the default database directory. The -d/--database option is not permitted"
            ));
            assert!(!stderr.contains("unexpected argument"));
        }

        #[test]
        fn upstream_run_accepts_default_application_flags() {
            let cli = Cli::parse_from(args_with_default_run_subcommand_from([
                "amrfinder",
                "-p",
                "prot.fa",
                "-q",
                "-debug",
                "-qc",
                "-verbose=2",
            ]));
            match cli.command {
                Commands::Run {
                    protein,
                    quiet,
                    debug,
                    qc,
                    verbose,
                    ..
                } => {
                    assert_eq!(protein, Some(PathBuf::from("prot.fa")));
                    assert!(quiet);
                    assert!(debug);
                    assert!(qc);
                    assert_eq!(verbose, 2);
                }
                _ => panic!("expected run command"),
            }
        }

        #[test]
        fn upstream_run_accepts_pgap_and_gpipe_org_flags() {
            let cli = Cli::parse_from(args_with_default_run_subcommand_from([
                "amrfinder",
                "-p",
                "prot.fa",
                "-pgap",
                "-gpipe_org",
            ]));
            match cli.command {
                Commands::Run {
                    protein,
                    pgap,
                    gpipe_org,
                    ..
                } => {
                    assert_eq!(protein, Some(PathBuf::from("prot.fa")));
                    assert!(pgap);
                    assert!(gpipe_org);
                }
                _ => panic!("expected run command"),
            }
        }

        #[test]
        fn upstream_run_accepts_original_full_key_spellings() {
            let cli = Cli::parse_from(args_with_default_run_subcommand_from([
                "amrfinder",
                "-protein=prot.fa",
                "-nucleotide",
                "nucl.fa",
                "-gff=annot.gff",
                "-dir",
                "/inputs",
                "-database",
                "/db",
                "-organism=Escherichia",
                "-ident_min",
                "0.9",
                "-coverage_min=0.8",
                "-threads=2",
                "-plus",
                "-report_common",
                "-report_all_equal",
                "-print_node",
                "-mutation_all=mut.tsv",
                "-annotation_format",
                "pgap",
                "-translation_table=4",
                "-name",
                "sample1",
                "-blast_bin=/blast",
                "-hmmer_bin",
                "/hmmer",
                "-parm",
                "-nosame -noblast -skip_hmm_check -bed",
                "-output=out.tsv",
                "-protein_output",
                "prot.out.fa",
                "-nucleotide_output=nucl.out.fa",
                "-nucleotide_flank5_output",
                "flank.out.fa",
                "-nucleotide_flank5_size=25",
                "-log=run.log",
                "-database_version",
                "-list_organisms",
            ]));

            match cli.command {
                Commands::Run {
                    protein,
                    nucleotide,
                    gff,
                    dir,
                    database,
                    organism,
                    ident_min,
                    coverage_min,
                    threads,
                    plus,
                    report_common,
                    report_all_equal,
                    print_node,
                    mutation_all,
                    annotation_format,
                    translation_table,
                    name,
                    blast_bin,
                    hmmer_bin,
                    parm,
                    output,
                    protein_output,
                    nucleotide_output,
                    nucleotide_flank5_output,
                    nucleotide_flank5_size,
                    log,
                    database_version,
                    list_organisms,
                    ..
                } => {
                    assert_eq!(protein, Some(PathBuf::from("prot.fa")));
                    assert_eq!(nucleotide, Some(PathBuf::from("nucl.fa")));
                    assert_eq!(gff, Some(PathBuf::from("annot.gff")));
                    assert_eq!(dir, Some(PathBuf::from("/inputs")));
                    assert_eq!(database, Some(PathBuf::from("/db")));
                    assert_eq!(organism, "Escherichia");
                    assert_eq!(ident_min, 0.9);
                    assert_eq!(coverage_min, 0.8);
                    assert_eq!(threads, 2);
                    assert!(plus);
                    assert!(report_common);
                    assert!(report_all_equal);
                    assert!(print_node);
                    assert_eq!(mutation_all, Some(PathBuf::from("mut.tsv")));
                    assert_eq!(annotation_format, "pgap");
                    assert_eq!(translation_table, 4);
                    assert_eq!(name, "sample1");
                    assert_eq!(blast_bin, "/blast");
                    assert_eq!(hmmer_bin, "/hmmer");
                    assert_eq!(parm, "-nosame -noblast -skip_hmm_check -bed");
                    assert_eq!(output, Some(PathBuf::from("out.tsv")));
                    assert_eq!(protein_output, Some(PathBuf::from("prot.out.fa")));
                    assert_eq!(nucleotide_output, Some(PathBuf::from("nucl.out.fa")));
                    assert_eq!(
                        nucleotide_flank5_output,
                        Some(PathBuf::from("flank.out.fa"))
                    );
                    assert_eq!(nucleotide_flank5_size, 25);
                    assert_eq!(log, Some(PathBuf::from("run.log")));
                    assert!(database_version);
                    assert!(list_organisms);
                }
                _ => panic!("expected run command"),
            }
        }

        #[test]
        fn upstream_run_checks_database_before_pgap_conflict_like_cpp() {
            let output = Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .args(["-p", "prot.fa", "-pgap", "-a", "standard"])
                .output()
                .unwrap();

            assert!(!output.status.success());
            let stderr = String::from_utf8(output.stderr).unwrap();
            assert!(stderr.contains("No valid AMRFinder database is found."));
            assert!(!stderr.contains("--pgap conflicts with GFF type \"standard\""));
        }

        #[test]
        fn upstream_update_flags_reject_database_override_like_cpp() {
            let output = Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .args(["-U", "-d", "/db"])
                .output()
                .unwrap();

            assert!(!output.status.success());
            let stderr = String::from_utf8(output.stderr).unwrap();
            assert!(stderr.contains(
                "AMRFinder update option (-u/--update) only operates on the default database directory. The -d/--database option is not permitted"
            ));
        }

        #[test]
        fn upstream_update_flags_reject_sequence_inputs_like_cpp() {
            let output = Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .args(["-u", "-p", "input.fa"])
                .output()
                .unwrap();

            assert!(!output.status.success());
            let stderr = String::from_utf8(output.stderr).unwrap();
            assert!(stderr.contains(
                "AMRFinder -u/--update option cannot be run with -n/--nucleotide or -p/--protein options"
            ));
        }

        #[test]
        fn late_update_flag_rejects_sequence_inputs_like_cpp() {
            let output = Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .args(["-p", "input.fa", "-u"])
                .output()
                .unwrap();

            assert!(!output.status.success());
            let stderr = String::from_utf8(output.stderr).unwrap();
            assert!(stderr.contains(
                "AMRFinder -u/--update option cannot be run with -n/--nucleotide or -p/--protein options"
            ));
        }

        #[test]
        fn explicit_rust_update_subcommand_still_accepts_database() {
            let cli = Cli::parse_from(["amrfinder", "update", "-d", "/db"]);
            match cli.command {
                Commands::Update {
                    database,
                    top_level_update,
                    ..
                } => {
                    assert_eq!(database, Some(PathBuf::from("/db")));
                    assert!(!top_level_update);
                }
                _ => panic!("expected update command"),
            }
        }

        #[test]
        fn extract_accepts_original_aa_flag_spelling() {
            let cli = Cli::parse_from(args_with_default_run_subcommand_from([
                "amrfinder",
                "extract",
                "input.fa",
                "targets.tsv",
                "-aa",
                "-verbose=1",
            ]));
            match cli.command {
                Commands::Extract {
                    fasta,
                    target,
                    aa,
                    verbose,
                    ..
                } => {
                    assert_eq!(fasta, PathBuf::from("input.fa"));
                    assert_eq!(target, PathBuf::from("targets.tsv"));
                    assert!(aa);
                    assert_eq!(verbose, 1);
                }
                _ => panic!("expected extract command"),
            }
        }

        #[test]
        fn executable_name_routes_to_fasta_check_flow_with_original_flags() {
            let cli = Cli::parse_from(args_with_default_run_subcommand_from([
                "/opt/amr/fasta_check",
                "input.fa",
                "-aa",
                "-stop_codon",
                "-ambig_max=7",
                "-out",
                "clean.fa",
            ]));

            match cli.command {
                Commands::CheckFasta {
                    input,
                    aa,
                    stop_codon,
                    ambig_max,
                    out,
                    ..
                } => {
                    assert_eq!(input, PathBuf::from("input.fa"));
                    assert!(aa);
                    assert!(stop_codon);
                    assert_eq!(ambig_max, 7);
                    assert_eq!(out, Some(PathBuf::from("clean.fa")));
                }
                _ => panic!("expected check-fasta command"),
            }
        }

        #[test]
        fn executable_name_routes_to_fasta_extract_flow_with_original_flags() {
            let cli = Cli::parse_from(args_with_default_run_subcommand_from([
                "/opt/amr/fasta_extract",
                "input.fa",
                "targets.tsv",
                "-aa",
                "-verbose=1",
            ]));

            match cli.command {
                Commands::Extract {
                    fasta,
                    target,
                    aa,
                    verbose,
                    ..
                } => {
                    assert_eq!(fasta, PathBuf::from("input.fa"));
                    assert_eq!(target, PathBuf::from("targets.tsv"));
                    assert!(aa);
                    assert_eq!(verbose, 1);
                }
                _ => panic!("expected extract command"),
            }
        }

        #[test]
        fn executable_name_routes_to_fasta2parts_flow() {
            let cli = Cli::parse_from(args_with_default_run_subcommand_from([
                "/opt/amr/fasta2parts",
                "input.fa",
                "3",
                "parts",
            ]));

            match cli.command {
                Commands::SplitFasta {
                    input,
                    parts_max,
                    dir,
                    ..
                } => {
                    assert_eq!(input, PathBuf::from("input.fa"));
                    assert_eq!(parts_max, 3);
                    assert_eq!(dir, PathBuf::from("parts"));
                }
                _ => panic!("expected split-fasta command"),
            }
        }

        #[test]
        fn executable_name_routes_to_mutate_flow_with_original_flags() {
            let cli = Cli::parse_from(args_with_default_run_subcommand_from([
                "/opt/amr/mutate",
                "input.fa",
                "mutations.tsv",
                "-aa",
                "-orig",
            ]));

            match cli.command {
                Commands::Mutate {
                    input,
                    mutations,
                    aa,
                    orig,
                    ..
                } => {
                    assert_eq!(input, PathBuf::from("input.fa"));
                    assert_eq!(mutations, PathBuf::from("mutations.tsv"));
                    assert!(aa);
                    assert!(orig);
                }
                _ => panic!("expected mutate command"),
            }
        }

        #[test]
        fn executable_name_routes_to_gff_check_flow_with_original_flags() {
            let cli = Cli::parse_from(args_with_default_run_subcommand_from([
                "/opt/amr/gff_check",
                "annot.gff",
                "-gfftype",
                "pgap",
                "-prot",
                "prot.fa",
                "-dna=nucl.fa",
                "-lcl",
                "-gff_prot_match=prot.match",
                "-gff_dna_match",
                "dna.match",
                "-verbose=1",
            ]));

            match cli.command {
                Commands::GffCheck {
                    gff,
                    gfftype,
                    prot,
                    dna,
                    lcl,
                    gff_prot_match,
                    gff_dna_match,
                    verbose,
                    ..
                } => {
                    assert_eq!(gff, PathBuf::from("annot.gff"));
                    assert_eq!(gfftype, "pgap");
                    assert_eq!(prot, Some(PathBuf::from("prot.fa")));
                    assert_eq!(dna, Some(PathBuf::from("nucl.fa")));
                    assert!(lcl);
                    assert_eq!(gff_prot_match, Some(PathBuf::from("prot.match")));
                    assert_eq!(gff_dna_match, Some(PathBuf::from("dna.match")));
                    assert_eq!(verbose, 1);
                }
                _ => panic!("expected gff-check command"),
            }
        }

        #[test]
        fn executable_name_routes_to_dna_mutation_flow_with_original_flags() {
            let cli = Cli::parse_from(args_with_default_run_subcommand_from([
                "/opt/amr/dna_mutation",
                "blastn.tsv",
                "mutation.tsv",
                "Escherichia_coli",
                "-mutation_all=all.tsv",
                "-name",
                "sample1",
                "-print_node",
            ]));

            match cli.command {
                Commands::DnaMutation {
                    blastn,
                    mutation,
                    organism,
                    mutation_all,
                    name,
                    print_node,
                    ..
                } => {
                    assert_eq!(blastn, PathBuf::from("blastn.tsv"));
                    assert_eq!(mutation, PathBuf::from("mutation.tsv"));
                    assert_eq!(organism, "Escherichia_coli");
                    assert_eq!(mutation_all, Some(PathBuf::from("all.tsv")));
                    assert_eq!(name, "sample1");
                    assert!(print_node);
                }
                _ => panic!("expected dna-mutation command"),
            }
        }

        #[test]
        fn executable_name_routes_to_disruption2genesymbol_flow_with_original_flags() {
            let cli = Cli::parse_from(args_with_default_run_subcommand_from([
                "/opt/amr/disruption2genesymbol",
                "nucl.fa",
                "prot.fa",
                "disruptions.tsv",
                "-gencode=4",
                "-prot_id_pos",
                "2",
            ]));

            match cli.command {
                Commands::Disruption2GeneSymbol {
                    nucl,
                    prot,
                    tab,
                    gencode,
                    prot_id_pos,
                    ..
                } => {
                    assert_eq!(nucl, PathBuf::from("nucl.fa"));
                    assert_eq!(prot, PathBuf::from("prot.fa"));
                    assert_eq!(tab, PathBuf::from("disruptions.tsv"));
                    assert_eq!(gencode, "4");
                    assert_eq!(prot_id_pos, "2");
                }
                _ => panic!("expected disruption2-gene-symbol command"),
            }
        }

        #[test]
        fn helper_executable_routes_accept_default_qc_and_log_flags() {
            let cli = Cli::parse_from(args_with_default_run_subcommand_from([
                "/opt/amr/fasta_check",
                "input.fa",
                "-qc",
                "-log=check.log",
            ]));
            match cli.command {
                Commands::CheckFasta { qc, log, .. } => {
                    assert!(qc);
                    assert_eq!(log, Some(PathBuf::from("check.log")));
                }
                _ => panic!("expected check-fasta command"),
            }

            let cli = Cli::parse_from(args_with_default_run_subcommand_from([
                "/opt/amr/dna_mutation",
                "blastn.tsv",
                "mutation.tsv",
                "Escherichia_coli",
                "-qc",
                "-log",
                "dna.log",
            ]));
            match cli.command {
                Commands::DnaMutation { qc, log, .. } => {
                    assert!(qc);
                    assert_eq!(log, Some(PathBuf::from("dna.log")));
                }
                _ => panic!("expected dna-mutation command"),
            }
        }

        #[test]
        fn executable_name_routes_to_amrfinder_update_flow_with_original_flags() {
            let cli = Cli::parse_from(args_with_default_run_subcommand_from([
                "/opt/amr/amrfinder_update",
                "-database=/db",
                "-blast_bin",
                "/blast",
                "-hmmer_bin=/hmmer",
                "-force_update",
                "-debug",
                "-q",
            ]));

            match cli.command {
                Commands::AmrFinderUpdate {
                    database,
                    blast_bin,
                    hmmer_bin,
                    force_update,
                    quiet,
                    debug,
                } => {
                    assert_eq!(database, Some(PathBuf::from("/db")));
                    assert_eq!(blast_bin, Some(PathBuf::from("/blast")));
                    assert_eq!(hmmer_bin, Some(PathBuf::from("/hmmer")));
                    assert!(force_update);
                    assert!(quiet);
                    assert!(debug);
                }
                _ => panic!("expected amr-finder-update command"),
            }
        }

        #[test]
        fn executable_name_routes_to_amrfinder_index_flow_with_original_flags() {
            let cli = Cli::parse_from(args_with_default_run_subcommand_from([
                "/opt/amr/amrfinder_index",
                "/db",
                "-blast_bin=/blast",
                "-hmmer_bin",
                "/hmmer",
                "-debug",
                "-q",
                "-log=index.log",
            ]));

            match cli.command {
                Commands::AmrFinderIndex {
                    database,
                    blast_bin,
                    hmmer_bin,
                    quiet,
                    debug,
                    log,
                } => {
                    assert_eq!(database, PathBuf::from("/db"));
                    assert_eq!(blast_bin, Some(PathBuf::from("/blast")));
                    assert_eq!(hmmer_bin, Some(PathBuf::from("/hmmer")));
                    assert!(quiet);
                    assert!(debug);
                    assert_eq!(log, Some(PathBuf::from("index.log")));
                }
                _ => panic!("expected amr-finder-index command"),
            }
        }

        #[test]
        fn database_version_rejects_sequence_inputs_like_upstream() {
            let db = tempfile::tempdir().unwrap();
            std::fs::write(db.path().join("version.txt"), "2026-03-24.1\n").unwrap();
            std::fs::write(db.path().join("database_format_version.txt"), "0.1.0\n").unwrap();
            std::fs::write(db.path().join("AMRProt.fa.phr"), "").unwrap();

            let output = Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .arg("-V")
                .arg("-p")
                .arg("input.fa")
                .arg("-d")
                .arg(db.path())
                .output()
                .unwrap();

            assert!(!output.status.success());
            let stdout = String::from_utf8(output.stdout).unwrap();
            assert!(stdout.contains("Database version: 2026-03-24.1"));
            let stderr = String::from_utf8(output.stderr).unwrap();
            assert!(stderr.contains("No processing of input.fa is done"));
        }

        #[test]
        fn database_version_reports_directory_and_version_like_upstream() {
            let db = tempfile::tempdir().unwrap();
            std::fs::write(db.path().join("version.txt"), "2026-03-24.1\n").unwrap();
            std::fs::write(db.path().join("database_format_version.txt"), "0.1.0\n").unwrap();
            std::fs::write(db.path().join("AMRProt.fa.phr"), "").unwrap();

            let output = Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .arg("-V")
                .arg("-d")
                .arg(db.path())
                .output()
                .unwrap();

            assert!(output.status.success());
            let stdout = String::from_utf8(output.stdout).unwrap();
            assert!(stdout.contains("Database directory: '"));
            assert!(stdout.contains("Database version: 2026-03-24.1\n"));
        }

        #[test]
        fn aggregate_run_accepts_upstream_default_log_key_before_database_version_exit() {
            let db = tempfile::tempdir().unwrap();
            std::fs::write(db.path().join("version.txt"), "2026-03-24.1\n").unwrap();
            std::fs::write(db.path().join("database_format_version.txt"), "0.1.0\n").unwrap();
            std::fs::write(db.path().join("AMRProt.fa.phr"), "").unwrap();
            let log = db.path().join("amrfinder.log");

            let output = Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .arg("-V")
                .arg("-d")
                .arg(db.path())
                .arg("--log")
                .arg(&log)
                .output()
                .unwrap();

            assert!(output.status.success());
            assert_eq!(String::from_utf8(output.stderr).unwrap(), "");
            assert!(String::from_utf8(output.stdout)
                .unwrap()
                .contains("Database version: 2026-03-24.1\n"));
            assert!(std::fs::read_to_string(log).unwrap().contains("$ "));
        }

        #[test]
        fn aggregate_debug_flag_routes_to_translated_integrity_flag() {
            let db = tempfile::tempdir().unwrap();
            std::fs::write(db.path().join("version.txt"), "2026-03-24.1\n").unwrap();
            std::fs::write(db.path().join("database_format_version.txt"), "0.1.0\n").unwrap();
            std::fs::write(db.path().join("AMRProt.fa.phr"), "").unwrap();

            let output = Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .arg("-debug")
                .arg("-qc")
                .arg("-V")
                .arg("-d")
                .arg(db.path())
                .output()
                .unwrap();

            assert!(output.status.success());
            assert_eq!(String::from_utf8(output.stderr).unwrap(), "");
            assert!(String::from_utf8(output.stdout)
                .unwrap()
                .contains("Database version: 2026-03-24.1\n"));
        }

        #[test]
        fn database_version_requires_blast_database_marker_like_upstream() {
            let db = tempfile::tempdir().unwrap();
            std::fs::write(db.path().join("version.txt"), "2026-03-24.1\n").unwrap();
            std::fs::write(db.path().join("database_format_version.txt"), "0.1.0\n").unwrap();

            let output = Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .arg("-V")
                .arg("-d")
                .arg(db.path())
                .output()
                .unwrap();

            assert!(!output.status.success());
            let stderr = String::from_utf8(output.stderr).unwrap();
            assert!(stderr.contains("The BLAST database for AMRProt.fa"));
        }

        #[test]
        fn database_version_checks_minimum_versions_like_upstream() {
            let db = tempfile::tempdir().unwrap();
            std::fs::write(db.path().join("version.txt"), "2025-09-22.1\n").unwrap();
            std::fs::write(db.path().join("database_format_version.txt"), "0.1.0\n").unwrap();
            std::fs::write(db.path().join("AMRProt.fa.phr"), "").unwrap();

            let output = Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .arg("-V")
                .arg("-d")
                .arg(db.path())
                .output()
                .unwrap();

            assert!(!output.status.success());
            let stderr = String::from_utf8(output.stderr).unwrap();
            assert!(stderr.contains("Software requires database version at least 2025-09-22.2"));
        }

        #[test]
        fn list_organisms_uses_translated_db2organisms_sources() {
            let db = tempfile::tempdir().unwrap();
            std::fs::write(db.path().join("version.txt"), "2026-03-24.1\n").unwrap();
            std::fs::write(db.path().join("database_format_version.txt"), "0.1.0\n").unwrap();
            std::fs::write(db.path().join("AMRProt.fa.phr"), "").unwrap();
            std::fs::write(
                db.path().join("taxgroup.tsv"),
                "taxgroup\tgpipe_org\nSalmonella\tSalmonella enterica\nEscherichia\tEscherichia coli\n",
            )
            .unwrap();
            std::fs::write(
                db.path().join("AMRProt-mutation.tsv"),
                "taxgroup\taccession\tgenesymbol\nKlebsiella\tWP_1\tompK36\nEscherichia\tWP_2\tgyrA\n",
            )
            .unwrap();
            std::fs::write(
                db.path().join("AMRProt-susceptible.tsv"),
                "taxgroup\taccession\tgenesymbol\nAcinetobacter\tWP_3\tcarO\nSalmonella\tWP_4\tompC\n",
            )
            .unwrap();

            let output = Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .arg("-l")
                .arg("-d")
                .arg(db.path())
                .output()
                .unwrap();

            assert!(output.status.success());
            assert_eq!(
                String::from_utf8(output.stdout).unwrap(),
                "\nAvailable --organism options: Acinetobacter, Escherichia, Klebsiella, Salmonella\n"
            );
        }

        #[test]
        fn list_organisms_validates_database_like_upstream() {
            let db = tempfile::tempdir().unwrap();
            std::fs::write(db.path().join("version.txt"), "2026-03-24.1\n").unwrap();
            std::fs::write(db.path().join("database_format_version.txt"), "0.1.0\n").unwrap();
            std::fs::write(db.path().join("taxgroup.tsv"), "taxgroup\nEscherichia\n").unwrap();
            std::fs::write(
                db.path().join("AMRProt-mutation.tsv"),
                "taxgroup\taccession\tgenesymbol\n",
            )
            .unwrap();
            std::fs::write(
                db.path().join("AMRProt-susceptible.tsv"),
                "taxgroup\taccession\tgenesymbol\n",
            )
            .unwrap();

            let output = Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .arg("-l")
                .arg("-d")
                .arg(db.path())
                .output()
                .unwrap();

            assert!(!output.status.success());
            let stderr = String::from_utf8(output.stderr).unwrap();
            assert!(stderr.contains("The BLAST database for AMRProt.fa"));
        }

        fn args_with_default_run_subcommand_from<const N: usize>(args: [&str; N]) -> Vec<OsString> {
            let mut args: Vec<OsString> = args.into_iter().map(OsString::from).collect();
            if let Some(command) = args
                .first()
                .and_then(|arg| std::path::Path::new(arg).file_stem())
                .and_then(|arg| arg.to_str())
                .and_then(|program| match program {
                    "fasta_check" => Some("check-fasta"),
                    "fasta_extract" => Some("extract"),
                    "fasta2parts" => Some("split-fasta"),
                    "mutate" => Some("mutate"),
                    "gff_check" => Some("gff-check"),
                    "dna_mutation" => Some("dna-mutation"),
                    "disruption2genesymbol" => Some("disruption2-gene-symbol"),
                    "amrfinder_update" => Some("amr-finder-update"),
                    "amrfinder_index" => Some("amr-finder-index"),
                    _ => None,
                })
            {
                args.insert(1, OsString::from(command));
                normalize_upstream_single_dash_flags(&mut args);
                return args;
            }
            let Some(first_arg) = args.get(1) else {
                return args;
            };
            let first_arg = first_arg.to_str();
            let known_subcommand = first_arg.is_some_and(|arg| {
                matches!(
                    arg,
                    "run"
                        | "update"
                        | "check-fasta"
                        | "extract"
                        | "split-fasta"
                        | "mutate"
                        | "gff-check"
                        | "dna-mutation"
                        | "disruption2-gene-symbol"
                        | "amr-finder-update"
                        | "amr-finder-index"
                )
            });
            if !known_subcommand {
                let mut update_positions = Vec::new();
                let mut force_update = false;
                let mut skip_value = false;
                for (pos, arg) in args.iter().enumerate().skip(1) {
                    let Some(value) = arg.to_str() else {
                        skip_value = false;
                        continue;
                    };
                    if skip_value {
                        skip_value = false;
                        continue;
                    }
                    match value {
                        "-u" | "--update" | "-update" => update_positions.push(pos),
                        "-U" | "--force_update" | "-force_update" => {
                            force_update = true;
                            update_positions.push(pos);
                        }
                        "-p"
                        | "--protein"
                        | "-protein"
                        | "-n"
                        | "--nucleotide"
                        | "-nucleotide"
                        | "-g"
                        | "--gff"
                        | "-gff"
                        | "-d"
                        | "--database"
                        | "-database"
                        | "--dir"
                        | "-dir"
                        | "-O"
                        | "--organism"
                        | "-organism"
                        | "-i"
                        | "--ident_min"
                        | "-ident_min"
                        | "-c"
                        | "--coverage_min"
                        | "-coverage_min"
                        | "--threads"
                        | "-threads"
                        | "--mutation_all"
                        | "-mutation_all"
                        | "-a"
                        | "--annotation_format"
                        | "-annotation_format"
                        | "-t"
                        | "--translation_table"
                        | "-translation_table"
                        | "--name"
                        | "-name"
                        | "--blast_bin"
                        | "-blast_bin"
                        | "--hmmer_bin"
                        | "-hmmer_bin"
                        | "--parm"
                        | "-parm"
                        | "-o"
                        | "--output"
                        | "-output"
                        | "--protein_output"
                        | "-protein_output"
                        | "--nucleotide_output"
                        | "-nucleotide_output"
                        | "--nucleotide_flank5_output"
                        | "-nucleotide_flank5_output"
                        | "--nucleotide_flank5_size"
                        | "-nucleotide_flank5_size"
                        | "--verbose"
                        | "-verbose"
                        | "--log"
                        | "-log" => skip_value = true,
                        _ => {}
                    }
                }
                if !update_positions.is_empty() {
                    for pos in update_positions.into_iter().rev() {
                        args.remove(pos);
                    }
                    args.insert(1, OsString::from("update"));
                    let mut next = 2;
                    if force_update {
                        args.insert(next, OsString::from("--force_update"));
                        next += 1;
                    }
                    args.insert(next, OsString::from("--__top_level_update"));
                    normalize_upstream_single_dash_flags(&mut args);
                    return args;
                }
            }
            match first_arg {
                Some("-v") => {
                    args[1] = OsString::from("--version");
                }
                Some(arg)
                    if arg.starts_with('-')
                        && arg != "-h"
                        && arg != "--help"
                        && arg != "--version" =>
                {
                    args.insert(1, OsString::from("run"));
                }
                _ => {}
            }
            normalize_upstream_single_dash_flags(&mut args);
            args
        }
    }
}
