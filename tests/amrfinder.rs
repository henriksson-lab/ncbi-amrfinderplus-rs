mod seq {
    pub use amrfinder::seq::*;
}
mod common {
    pub use amrfinder::common::*;
}
mod tsv {
    pub use amrfinder::tsv::*;
}
mod fasta2parts {
    pub use amrfinder::fasta2parts::*;
}
mod fasta_extract {
    pub use amrfinder::fasta_extract::*;
}
mod fasta_check {
    pub use amrfinder::fasta_check::*;
}
mod gff {
    pub use amrfinder::gff::*;
}
mod gff_check {
    pub use amrfinder::gff_check::*;
}
mod amr_report {
    pub use amrfinder::amr_report::*;
}
mod amrfinder_update {
    pub use amrfinder::amrfinder_update::*;
}
mod dna_mutation {
    pub use amrfinder::dna_mutation::*;
}
mod disruption2genesymbol {
    pub use amrfinder::disruption2genesymbol::*;
}
mod columns {
    pub use amrfinder::columns::*;
}

mod amrfinder_impl {
    include!("../src/amrfinder.rs");

    #[cfg(test)]
    mod tests {
        use super::*;
        use std::sync::Mutex;

        const AMRFINDER_VERSION: &str = include_str!("../amr/version.txt");
        static ENV_LOCK: Mutex<()> = Mutex::new(());

        #[test]
        fn test_pipeline_config_default() {
            let config = PipelineConfig::default();
            assert_eq!(config.threads, 4);
            assert_eq!(config.coverage_min, 0.5);
            assert!(!config.plus);
        }

        #[test]
        fn this_application_defaults_threads_to_original_amrfinder_value() {
            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from("/opt/amr/amrfinder"),
                    std::ffi::OsString::from("-p"),
                    std::ffi::OsString::from("prot.fa"),
                ])
                .unwrap();

            assert_eq!(run.key_args["threads"], "4");
        }

        #[test]
        fn this_application_parses_original_amrfinder_keys() {
            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from("/opt/amr/amrfinder"),
                    std::ffi::OsString::from("-p"),
                    std::ffi::OsString::from("prot.fa"),
                    std::ffi::OsString::from("-n"),
                    std::ffi::OsString::from("dna.fa"),
                    std::ffi::OsString::from("-g"),
                    std::ffi::OsString::from("annot.gff"),
                    std::ffi::OsString::from("-d"),
                    std::ffi::OsString::from("/db"),
                    std::ffi::OsString::from("-O"),
                    std::ffi::OsString::from("Escherichia coli"),
                    std::ffi::OsString::from("-i"),
                    std::ffi::OsString::from("0.9"),
                    std::ffi::OsString::from("-c"),
                    std::ffi::OsString::from("0.8"),
                    std::ffi::OsString::from("--threads"),
                    std::ffi::OsString::from("3"),
                    std::ffi::OsString::from("--plus"),
                    std::ffi::OsString::from("--report_common"),
                    std::ffi::OsString::from("--print_node"),
                    std::ffi::OsString::from("--pgap"),
                    std::ffi::OsString::from("--gpipe_org"),
                    std::ffi::OsString::from("--mutation_all"),
                    std::ffi::OsString::from("mut.tsv"),
                    std::ffi::OsString::from("-a"),
                    std::ffi::OsString::from("pgap"),
                    std::ffi::OsString::from("-t"),
                    std::ffi::OsString::from("4"),
                    std::ffi::OsString::from("--blast_bin"),
                    std::ffi::OsString::from("/blast"),
                    std::ffi::OsString::from("--hmmer_bin"),
                    std::ffi::OsString::from("/hmmer"),
                    std::ffi::OsString::from("--debug"),
                    std::ffi::OsString::from("-qc"),
                    std::ffi::OsString::from("--parm"),
                    std::ffi::OsString::from("-nosame -noblast -skip_hmm_check -bed"),
                    std::ffi::OsString::from("-o"),
                    std::ffi::OsString::from("out.tsv"),
                ])
                .unwrap();

            assert_eq!(run.key_args["protein"], "prot.fa");
            assert_eq!(run.key_args["nucleotide"], "dna.fa");
            assert_eq!(run.key_args["gff"], "annot.gff");
            assert_eq!(run.key_args["database"], "/db");
            assert_eq!(run.key_args["organism"], "Escherichia coli");
            assert_eq!(run.key_args["ident_min"], "0.9");
            assert_eq!(run.key_args["coverage_min"], "0.8");
            assert_eq!(run.key_args["threads"], "3");
            assert_eq!(run.key_args["plus"], "true");
            assert_eq!(run.key_args["report_common"], "true");
            assert_eq!(run.key_args["print_node"], "true");
            assert_eq!(run.key_args["pgap"], "true");
            assert_eq!(run.key_args["gpipe_org"], "true");
            assert_eq!(run.key_args["mutation_all"], "mut.tsv");
            assert_eq!(run.key_args["annotation_format"], "pgap");
            assert_eq!(run.key_args["translation_table"], "4");
            assert_eq!(run.key_args["blast_bin"], "/blast");
            assert_eq!(run.key_args["hmmer_bin"], "/hmmer");
            assert_eq!(run.key_args["debug"], "true");
            assert_eq!(run.key_args["qc"], "true");
            assert_eq!(
                run.key_args["parm"],
                "-nosame -noblast -skip_hmm_check -bed"
            );
            assert_eq!(run.key_args["output"], "out.tsv");
        }

        #[test]
        fn update_delegation_preserves_original_debug_flag() {
            let source = include_str!("../src/amrfinder.rs");
            assert!(source.contains("argv.push(OsString::from(\"--debug\"));"));
            assert!(source.contains("} else if run.key_args[\"qc\"] == \"true\" {"));
        }

        #[test]
        fn update_delegation_passes_parent_of_latest_database_like_upstream() {
            let source = include_str!("../src/amrfinder.rs");
            assert!(source.contains(".file_name()"));
            assert!(source.contains("|name| name == \"latest\""));
            assert!(source.contains(".parent()"));
            assert!(source.contains("amrfinder_update"));
        }

        #[test]
        fn update_delegation_uses_original_short_database_key() {
            let source = include_str!("../src/amrfinder.rs");
            let update_branch = source
                .find("if update {")
                .expect("translated update branch should be present");
            let update_call = source[update_branch..]
                .find("crate::amrfinder_update::main(argv)?;")
                .expect("translated update delegation should be present");
            let branch = &source[update_branch..update_branch + update_call];

            assert!(branch.contains("argv.push(OsString::from(\"-d\"));"));
            assert!(!branch.contains("argv.push(OsString::from(\"-database\"));"));
        }

        #[test]
        fn update_delegation_forwards_tool_directories_like_upstream() {
            let source = include_str!("../src/amrfinder.rs");
            assert!(source.contains("argv.push(OsString::from(\"-blast_bin\"));"));
            assert!(source.contains("argv.push(OsString::from(&blast_bin));"));
            assert!(source.contains("argv.push(OsString::from(\"-hmmer_bin\"));"));
            assert!(source.contains("argv.push(OsString::from(&run.key_args[\"hmmer_bin\"]));"));
        }

        #[test]
        fn update_delegation_resolves_blast_bin_env_before_updater_like_upstream() {
            let source = include_str!("../src/amrfinder.rs");
            let blast_env = source
                .find("std::env::var_os(\"BLAST_BIN\")")
                .expect("translated shell body should consult BLAST_BIN");
            let update_branch = source
                .find("if update {")
                .expect("translated update branch should be present");
            assert!(blast_env < update_branch);

            let update_call = source[update_branch..]
                .find("crate::amrfinder_update::main(argv)?;")
                .expect("translated update delegation should be present");
            let branch = &source[update_branch..update_branch + update_call];
            assert!(branch.contains("if !blast_bin.is_empty()"));
            assert!(branch.contains("argv.push(OsString::from(&blast_bin));"));
            assert!(!branch.contains("argv.push(OsString::from(&run.key_args[\"blast_bin\"]));"));
        }

        #[test]
        fn update_delegation_forwards_original_force_update_flag_spelling() {
            let source = include_str!("../src/amrfinder.rs");
            let update_branch = source
                .find("if update {")
                .expect("translated update branch should be present");
            let update_call = source[update_branch..]
                .find("crate::amrfinder_update::main(argv)?;")
                .expect("translated update delegation should be present");
            let branch = &source[update_branch..update_branch + update_call];

            assert!(branch.contains("argv.push(OsString::from(\"--force_update\"));"));
            assert!(!branch.contains("argv.push(OsString::from(\"-force_update\"));"));
        }

        #[test]
        fn update_delegation_continues_to_database_validation_like_upstream() {
            let source = include_str!("../src/amrfinder.rs");
            let update_branch = source
                .find("if update {")
                .expect("translated update branch should be present");
            let database_version_branch = source[update_branch..]
                .find("if run.key_args[\"database_version\"] == \"true\"")
                .expect("database-version branch should follow update delegation");
            let update_block = &source[update_branch..update_branch + database_version_branch];
            assert!(update_block.contains("crate::amrfinder_update::main(argv)?;"));
            assert!(!update_block.contains("return Ok(());"));
            let database_validation = source[update_branch + database_version_branch..]
                .find("let _ = database_version(&database)?;")
                .expect("translated update path should validate the database after update");
            let no_input_return = source[update_branch + database_version_branch..]
                .find("if protein.is_none() && nucleotide.is_none()")
                .expect("translated update path should return from the no-input branch");
            assert!(database_validation < no_input_return);
        }

        #[test]
        fn update_delegation_skips_updater_for_non_latest_database_like_upstream() {
            let source = include_str!("../src/amrfinder.rs");
            let update_database_branch = source
                .find("let update_database = PathBuf::from(&run.key_args[\"database\"]);")
                .expect("translated update database branch should be present");
            let database_version_branch = source[update_database_branch..]
                .find("if run.key_args[\"database_version\"] == \"true\"")
                .expect("database-version branch should follow update delegation");
            let branch =
                &source[update_database_branch..update_database_branch + database_version_branch];
            let latest = branch
                .find("|name| name == \"latest\"")
                .expect("latest directory gate should be present");
            let updater = branch
                .find("crate::amrfinder_update::main(argv)?;")
                .expect("updater should be called inside latest branch");
            let warning = branch
                .find("Updating database directory works only for databases with the default data directory format")
                .expect("non-latest branch should warn like upstream");
            assert!(latest < updater);
            assert!(updater < warning);
        }

        #[test]
        fn update_delegation_ignores_amrfinder_db_like_upstream() {
            let source = include_str!("../src/amrfinder.rs");
            let update_branch = source
                .find("if update {")
                .expect("translated update branch should be present");
            let update_call = source[update_branch..]
                .find("crate::amrfinder_update::main(argv)?;")
                .expect("translated update delegation should be present");
            let branch = &source[update_branch..update_branch + update_call];
            assert!(branch.contains("std::env::var_os(\"AMRFINDER_DB\").is_some()"));
            assert!(branch.contains(
                "AMRFinder auto-update only downloads to the default database directory"
            ));
            assert!(branch.contains("database = PathBuf::from(&run.key_args[\"database\"]);"));
        }

        #[test]
        fn blast_thread_caps_match_upstream_integer_division_without_rounding_up() {
            let source = include_str!("../src/amrfinder.rs");
            assert!(source.contains("std::cmp::min(n_prot, prot_len_total / 10000),"));
            assert!(source.contains("std::cmp::min(n_dna, dna_len_total / 10002)"));
            assert!(!source.contains("prot_len_total / 10000 + 1"));
            assert!(!source.contains("dna_len_total / 10002 + 1"));
        }

        #[test]
        fn blast_num_threads_is_omitted_for_zero_or_one_thread_like_upstream() {
            let source = include_str!("../src/amrfinder.rs");
            let blastp_threads = source
                .find("let blastp_threads = std::cmp::min(")
                .expect("translated blastp thread limit should be present");
            let blastp_run = source[blastp_threads..]
                .find("let blastp_output = Command::new(&blastp_prog)")
                .expect("translated blastp command should be present");
            let blastp_block = &source[blastp_threads..blastp_threads + blastp_run];
            assert!(blastp_block.contains("if blastp_threads > 1"));
            assert!(
                blastp_block.contains("\"-num_threads\".to_string(), blastp_threads.to_string()")
            );

            let blast_threads = source
                .find("let blast_threads =")
                .expect("translated blastx thread limit should be present");
            let blastx_run = source[blast_threads..]
                .find("let blast_status = Command::new(&blast_prog)")
                .expect("translated blastx command should be present");
            let blastx_block = &source[blast_threads..blast_threads + blastx_run];
            assert!(blastx_block.contains("if blast_threads > 1"));
            assert!(
                blastx_block.contains("\"-num_threads\".to_string(), blast_threads.to_string()")
            );
        }

        #[test]
        fn chunked_searches_use_upstream_threads_minus_one_worker_limit() {
            let source = include_str!("../src/amrfinder.rs");
            let hmm_fn = source
                .find("fn run_hmmsearch(")
                .expect("translated hmmsearch runner should be present");
            let tblastn_fn = source
                .find("fn run_tblastn_search(")
                .expect("translated tblastn runner should be present");
            let parse_report_fn = source
                .find("fn parse_report_blastx_args(")
                .expect("tblastn runner should end before report argument parser");
            let hmm_block = &source[hmm_fn..tblastn_fn];
            let tblastn_block = &source[tblastn_fn..parse_report_fn];

            assert!(hmm_block.contains("let worker_limit = threads.saturating_sub(1).max(1);"));
            assert!(hmm_block.contains("Vec::with_capacity(worker_limit)"));
            assert!(hmm_block.contains("if handles.len() == worker_limit"));
            assert!(tblastn_block.contains("let worker_limit = threads.saturating_sub(1).max(1);"));
            assert!(tblastn_block.contains("Vec::with_capacity(worker_limit)"));
            assert!(tblastn_block.contains("if handles.len() == worker_limit"));
            assert!(!hmm_block.contains("let mut handles = Vec::with_capacity(chunks.len())"));
            assert!(!tblastn_block.contains("let mut handles = Vec::with_capacity(chunks.len())"));
        }

        #[test]
        fn protein_search_runs_blastp_before_hmmsearch_like_upstream() {
            let source = include_str!("../src/amrfinder.rs");
            let protein_block = source
                .find("// Process protein input")
                .expect("protein processing block should be present");
            let nucleotide_block = source[protein_block..]
                .find("// Process nucleotide input")
                .expect("nucleotide processing block should follow protein processing");
            let block = &source[protein_block..protein_block + nucleotide_block];

            let blastp_status = block
                .find("let blastp_output = Command::new(&blastp_prog)")
                .expect("blastp command should be present");
            let blastp_error = block
                .find("bail!(\"blastp failed: {}\", stderr);")
                .expect("blastp failure check should be present");
            let hmm_run = block
                .find("run_hmmsearch(")
                .expect("hmmsearch invocation should be present");

            assert!(blastp_status < blastp_error);
            assert!(blastp_error < hmm_run);
            assert!(!block.contains("hmm_handle"));
        }

        #[test]
        fn mutation_all_without_organism_warns_before_database_failure_like_upstream() {
            let output = std::process::Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .arg("--mutation_all")
                .arg("mut.tsv")
                .arg("-d")
                .arg("/missing/amrfinder/db/latest")
                .output()
                .unwrap();

            assert!(!output.status.success());
            let stderr = String::from_utf8(output.stderr).unwrap();
            let warning = stderr
                .find("WARNING: --mutation_all option used without -O/--organism option. No point mutations will be screened")
                .expect("upstream mutation_all warning should be printed");
            let database_error = stderr
                .find("No valid AMRFinder database is found.")
                .expect("missing database should still fail before search execution");
            assert!(warning < database_error);
        }

        #[test]
        fn output_preflight_precedes_database_validation_like_upstream() {
            let tmp = tempfile::tempdir().unwrap();
            let output = tmp.path().join("missing_parent").join("out.tsv");

            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from("amrfinder"),
                    std::ffi::OsString::from("-o"),
                    std::ffi::OsString::from(&output),
                    std::ffi::OsString::from("-d"),
                    std::ffi::OsString::from("/missing/amrfinder/db/latest"),
                    std::ffi::OsString::from("-p"),
                    std::ffi::OsString::from("prot.fa"),
                ])
                .unwrap();

            let err = app.shell_body(&run, &mut Vec::new()).unwrap_err();

            assert_eq!(
                err.to_string(),
                format!("Cannot create output file '{}'", output.display())
            );
            assert!(!output.exists());
        }

        #[test]
        fn shell_body_uses_blast_bin_environment_fallback_like_upstream() {
            let source = include_str!("../src/amrfinder.rs");
            assert!(source.contains("let mut blast_bin = run.key_args[\"blast_bin\"].clone();"));
            assert!(source.contains("std::env::var_os(\"BLAST_BIN\")"));
            assert!(!source.contains("blast_bin: run.key_args[\"blast_bin\"].clone(),"));
            assert!(source.contains("blast_bin,"));
        }

        #[test]
        fn this_application_accepts_original_dir_key() {
            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from("amrfinder"),
                    std::ffi::OsString::from("-dir"),
                    std::ffi::OsString::from("/inputs"),
                    std::ffi::OsString::from("-protein"),
                    std::ffi::OsString::from("prot.fa"),
                ])
                .unwrap();

            assert_eq!(run.key_args["dir"], "/inputs");
            assert_eq!(run.key_args["protein"], "prot.fa");
        }

        #[test]
        fn this_application_shell_body_handles_database_version_without_pipeline() {
            let db = tempfile::tempdir().unwrap();
            fs::write(db.path().join("version.txt"), "2026-03-24.1\n").unwrap();
            fs::write(db.path().join("database_format_version.txt"), "0.1.0\n").unwrap();
            fs::write(db.path().join("AMRProt.fa.phr"), "").unwrap();

            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from("amrfinder"),
                    std::ffi::OsString::from("-V"),
                    std::ffi::OsString::from("-d"),
                    std::ffi::OsString::from(db.path()),
                ])
                .unwrap();
            let mut out = Vec::new();

            app.shell_body(&run, &mut out).unwrap();

            let out = String::from_utf8(out).unwrap();
            assert!(out.contains("Software directory: '"));
            assert!(out.contains(&format!(
                "Software version: {}\n",
                AMRFINDER_VERSION.trim_end()
            )));
            assert!(out.contains("Database directory: '"));
            assert!(out.contains("Database version: 2026-03-24.1\n"));
        }

        #[test]
        fn database_version_uses_amrfinder_db_when_database_arg_empty() {
            let _guard = ENV_LOCK.lock().unwrap();
            let db = tempfile::tempdir().unwrap();
            fs::write(db.path().join("version.txt"), "2026-03-24.1\n").unwrap();
            fs::write(db.path().join("database_format_version.txt"), "0.1.0\n").unwrap();
            fs::write(db.path().join("AMRProt.fa.phr"), "").unwrap();

            let previous = std::env::var_os("AMRFINDER_DB");
            std::env::set_var("AMRFINDER_DB", db.path());

            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from("amrfinder"),
                    std::ffi::OsString::from("-V"),
                ])
                .unwrap();
            let mut out = Vec::new();
            let result = app.shell_body(&run, &mut out);

            if let Some(previous) = previous {
                std::env::set_var("AMRFINDER_DB", previous);
            } else {
                std::env::remove_var("AMRFINDER_DB");
            }

            result.unwrap();
            let out = String::from_utf8(out).unwrap();
            assert!(out.contains("Software directory: '"));
            assert!(out.contains(&format!(
                "Software version: {}\n",
                AMRFINDER_VERSION.trim_end()
            )));
            assert!(out.contains(&format!("Database directory: '{}'", db.path().display())));
            assert!(out.contains("Database version: 2026-03-24.1\n"));
        }

        #[test]
        fn database_version_with_protein_reports_db_then_no_processing_like_upstream() {
            let db = tempfile::tempdir().unwrap();
            fs::write(db.path().join("version.txt"), "2026-03-24.1\n").unwrap();
            fs::write(db.path().join("database_format_version.txt"), "0.1.0\n").unwrap();
            fs::write(db.path().join("AMRProt.fa.phr"), "").unwrap();

            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from("amrfinder"),
                    std::ffi::OsString::from("-V"),
                    std::ffi::OsString::from("-d"),
                    std::ffi::OsString::from(db.path()),
                    std::ffi::OsString::from("-p"),
                    std::ffi::OsString::from("prot.fa"),
                ])
                .unwrap();
            let mut out = Vec::new();

            let err = app.shell_body(&run, &mut out).unwrap_err();

            assert_eq!(err.to_string(), "No processing of prot.fa is done");
            let out = String::from_utf8(out).unwrap();
            assert!(out.contains("Software directory: '"));
            assert!(out.contains(&format!(
                "Software version: {}\n",
                AMRFINDER_VERSION.trim_end()
            )));
            assert!(out.contains("Database directory: '"));
            assert!(out.contains("Database version: 2026-03-24.1\n"));
        }

        #[test]
        fn shell_body_prepends_dir_to_sequence_inputs_like_upstream() {
            let db = tempfile::tempdir().unwrap();
            fs::write(db.path().join("version.txt"), "2026-03-24.1\n").unwrap();
            fs::write(db.path().join("database_format_version.txt"), "0.1.0\n").unwrap();
            fs::write(db.path().join("AMRProt.fa.phr"), "").unwrap();
            let inputs = tempfile::tempdir().unwrap();

            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from("amrfinder"),
                    std::ffi::OsString::from("-V"),
                    std::ffi::OsString::from("-database"),
                    std::ffi::OsString::from(db.path()),
                    std::ffi::OsString::from("-dir"),
                    std::ffi::OsString::from(inputs.path()),
                    std::ffi::OsString::from("-protein"),
                    std::ffi::OsString::from("prot.fa"),
                ])
                .unwrap();
            let mut out = Vec::new();

            let err = app.shell_body(&run, &mut out).unwrap_err();

            assert_eq!(
                err.to_string(),
                format!(
                    "No processing of {} is done",
                    inputs.path().join("prot.fa").display()
                )
            );
        }

        #[test]
        fn database_version_missing_database_precedes_no_processing_error() {
            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from("amrfinder"),
                    std::ffi::OsString::from("-V"),
                    std::ffi::OsString::from("-d"),
                    std::ffi::OsString::from("/missing/amrfinder/db/latest"),
                    std::ffi::OsString::from("-p"),
                    std::ffi::OsString::from("prot.fa"),
                ])
                .unwrap();
            let mut out = Vec::new();

            let err = app.shell_body(&run, &mut out).unwrap_err();

            assert!(err.to_string().contains(
                "No valid AMRFinder database is found.\nThis directory (or symbolic link to directory) is not found: /missing/amrfinder/db/latest"
            ));
            let out = String::from_utf8(out).unwrap();
            assert!(out.contains("Software directory: '"));
            assert!(out.contains(&format!(
                "Software version: {}\n",
                AMRFINDER_VERSION.trim_end()
            )));
            assert!(!out.contains("Database directory: '"));
        }

        #[test]
        fn database_version_rejects_name_tab_before_database_work_like_upstream() {
            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from("amrfinder"),
                    std::ffi::OsString::from("-V"),
                    std::ffi::OsString::from("-d"),
                    std::ffi::OsString::from("/missing/amrfinder/db/latest"),
                    std::ffi::OsString::from("--name"),
                    std::ffi::OsString::from("bad\tname"),
                ])
                .unwrap();
            let mut out = Vec::new();

            let err = app.shell_body(&run, &mut out).unwrap_err();

            assert_eq!(err.to_string(), "NAME cannot contain a tab character");
            let out = String::from_utf8(out).unwrap();
            assert!(out.contains("Software directory: '"));
            assert!(out.contains(&format!(
                "Software version: {}\n",
                AMRFINDER_VERSION.trim_end()
            )));
            assert!(!out.contains("Database directory: '"));
        }

        #[test]
        fn database_version_rejects_search_thresholds_before_database_work_like_upstream() {
            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from("amrfinder"),
                    std::ffi::OsString::from("-V"),
                    std::ffi::OsString::from("-d"),
                    std::ffi::OsString::from("/missing/amrfinder/db/latest"),
                    std::ffi::OsString::from("-i"),
                    std::ffi::OsString::from("1.5"),
                ])
                .unwrap();
            let err = app.shell_body(&run, &mut Vec::new()).unwrap_err();

            assert_eq!(err.to_string(), "ident_min must be between 0 and 1");

            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from("amrfinder"),
                    std::ffi::OsString::from("-V"),
                    std::ffi::OsString::from("-d"),
                    std::ffi::OsString::from("/missing/amrfinder/db/latest"),
                    std::ffi::OsString::from("-c"),
                    std::ffi::OsString::from("-0.1"),
                ])
                .unwrap();
            let err = app.shell_body(&run, &mut Vec::new()).unwrap_err();

            assert_eq!(err.to_string(), "coverage_min must be between 0 and 1");
        }

        #[test]
        fn database_version_rejects_annotation_format_before_database_work_like_upstream() {
            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from("amrfinder"),
                    std::ffi::OsString::from("-V"),
                    std::ffi::OsString::from("-d"),
                    std::ffi::OsString::from("/missing/amrfinder/db/latest"),
                    std::ffi::OsString::from("-a"),
                    std::ffi::OsString::from("not_a_gff_type"),
                ])
                .unwrap();
            let mut out = Vec::new();

            let err = app.shell_body(&run, &mut out).unwrap_err();

            assert_eq!(err.to_string(), "Unknown GFF type: \"not_a_gff_type\"");
            assert!(String::from_utf8(out).unwrap().is_empty());
        }

        #[test]
        fn database_version_missing_marker_prints_directory_before_error_like_upstream() {
            let db = tempfile::tempdir().unwrap();
            fs::write(db.path().join("version.txt"), "2026-03-24.1\n").unwrap();
            fs::write(db.path().join("database_format_version.txt"), "0.1.0\n").unwrap();

            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from("amrfinder"),
                    std::ffi::OsString::from("-V"),
                    std::ffi::OsString::from("-d"),
                    std::ffi::OsString::from(db.path()),
                ])
                .unwrap();
            let mut out = Vec::new();

            let err = app.shell_body(&run, &mut out).unwrap_err();

            assert!(err
                .to_string()
                .contains("The BLAST database for AMRProt.fa"));
            let out = String::from_utf8(out).unwrap();
            assert!(out.contains("Database directory: '"));
            assert!(!out.contains("Database version: "));
        }

        #[test]
        fn database_version_minimum_data_error_prints_version_before_error_like_upstream() {
            let db = tempfile::tempdir().unwrap();
            fs::write(db.path().join("version.txt"), "2025-09-22.1\n").unwrap();
            fs::write(db.path().join("database_format_version.txt"), "0.1.0\n").unwrap();
            fs::write(db.path().join("AMRProt.fa.phr"), "").unwrap();

            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from("amrfinder"),
                    std::ffi::OsString::from("-V"),
                    std::ffi::OsString::from("-d"),
                    std::ffi::OsString::from(db.path()),
                ])
                .unwrap();
            let mut out = Vec::new();

            let err = app.shell_body(&run, &mut out).unwrap_err();

            assert!(err
                .to_string()
                .contains("Software requires database version at least 2025-09-22.2"));
            let out = String::from_utf8(out).unwrap();
            assert!(out.contains("Database directory: '"));
            assert!(out.contains("Database version: 2025-09-22.1\n"));
        }

        #[test]
        fn this_application_update_rejects_sequence_inputs_before_network() {
            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from("amrfinder"),
                    std::ffi::OsString::from("-u"),
                    std::ffi::OsString::from("-p"),
                    std::ffi::OsString::from("prot.fa"),
                ])
                .unwrap();
            let err = app.shell_body(&run, &mut Vec::new()).unwrap_err();

            assert_eq!(
                err.to_string(),
                "AMRFinder -u/--update option cannot be run with -n/--nucleotide or -p/--protein options"
            );
        }

        #[test]
        fn this_application_update_rejects_database_override_before_updater() {
            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from("amrfinder"),
                    std::ffi::OsString::from("-U"),
                    std::ffi::OsString::from("-d"),
                    std::ffi::OsString::from("/db"),
                ])
                .unwrap();
            let err = app.shell_body(&run, &mut Vec::new()).unwrap_err();

            assert_eq!(
                err.to_string(),
                "AMRFinder update option (-u/--update) only operates on the default database directory. The -d/--database option is not permitted"
            );
        }

        #[test]
        fn blastx_report_args_include_dna_len_path() {
            let parsed = parse_report_blastx_args("-blastx /tmp/blastx -dna_len /tmp/len");
            assert_eq!(parsed.blastx_path, Some(PathBuf::from("/tmp/blastx")));
            assert_eq!(parsed.dna_len_path, Some(PathBuf::from("/tmp/len")));
        }

        #[test]
        fn rust_amr_report_handoff_applies_parm_nosame_like_upstream() {
            let tmp = tempfile::tempdir().unwrap();
            let db = tmp.path().join("db");
            fs::create_dir(&db).unwrap();
            fs::write(
                db.join("fam.tsv"),
                "fam0\t\tgene0\t-\t0\t0\t90\t0\t90\t90\t0\t50\t1\tAMR\tAMR\tCLASS\tSUBCLASS\tFamily 0\n",
            )
            .unwrap();
            fs::write(db.join("AMRProt-mutation.tsv"), "# no mutation rows\n").unwrap();
            fs::write(
                db.join("AMRProt-susceptible.tsv"),
                "# no susceptible rows\n",
            )
            .unwrap();
            let blastp = tmp.path().join("blastp");
            fs::write(
                &blastp,
                concat!(
                    "WP_0|1|1|fam0|gene0|AMR|1|SUBCLASS|CLASS|Product\tWP_0\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*\n",
                    "WP_0|1|1|fam0|gene0|AMR|1|SUBCLASS|CLASS|Product\tprot_keep\t1\t5\t5\t1\t5\t5\tAAAA*\tAAAA*\n",
                ),
            )
            .unwrap();
            let config = PipelineConfig {
                database: db.clone(),
                parm: "-nosame".to_string(),
                plus: true,
                ..PipelineConfig::default()
            };

            let (output, _fasta_source) = run_rust_amr_report(
                &config,
                tmp.path(),
                db.to_str().unwrap(),
                "",
                &format!("-blastp {}", blastp.display()),
                "",
                None,
                false,
            )
            .unwrap();

            assert!(output.contains("\nprot_keep\t"), "{output}");
            assert!(!output.contains("\nWP_0\t"), "{output}");
        }

        #[test]
        fn rust_amr_report_handoff_parses_retain_and_skip_hmm_flags_like_upstream_parm() {
            let source = include_str!("../src/amrfinder.rs");
            let handoff = source
                .find("fn run_rust_amr_report(")
                .expect("translated report handoff should be present");
            let report_config = source[handoff..]
                .find("let report_config = AmrReportConfig")
                .expect("report handoff should construct AmrReportConfig");
            let block = &source[handoff..handoff + report_config];

            assert!(block.contains("let mut retain_blasts = false;"));
            assert!(block.contains("let mut skip_hmm_check = false;"));
            assert!(block.contains("let mut print_node = config.print_node;"));
            assert!(block.contains("let mut print_node_raw = false;"));
            assert!(block.contains("\"-retain_blasts\" => retain_blasts = true,"));
            assert!(block.contains("\"-skip_hmm_check\" => skip_hmm_check = true,"));
            assert!(block.contains("\"-print_node\" => print_node = true,"));
            assert!(block.contains("\"-print_node_raw\" => print_node_raw = true,"));

            let config_block = &source[handoff + report_config..];
            assert!(config_block.contains("retain_blasts,"));
            assert!(config_block.contains("skip_hmm_check,"));
            assert!(config_block.contains("print_node,"));
            assert!(config_block.contains("print_node_raw,"));
        }

        #[test]
        fn rust_amr_report_handoff_forwards_print_node_raw_from_parm_like_upstream() {
            let tmp = tempfile::tempdir().unwrap();
            let db = tmp.path().join("db");
            fs::create_dir(&db).unwrap();
            fs::write(
                db.join("fam.tsv"),
                "fam0\t\tgene0\t-\t0\t0\t90\t0\t90\t90\t0\t50\t1\tAMR\tAMR\tCLASS\tSUBCLASS\tFamily 0\n",
            )
            .unwrap();
            fs::write(db.join("AMRProt-mutation.tsv"), "# no mutation rows\n").unwrap();
            fs::write(
                db.join("AMRProt-susceptible.tsv"),
                "# no susceptible rows\n",
            )
            .unwrap();

            let config = PipelineConfig {
                database: db.clone(),
                parm: "-print_node_raw".to_string(),
                plus: true,
                ..PipelineConfig::default()
            };

            let err = run_rust_amr_report(
                &config,
                tmp.path(),
                db.to_str().unwrap(),
                "",
                "",
                "",
                None,
                false,
            )
            .unwrap_err();

            assert_eq!(err.to_string(), "print_node_raw requires print_node");
        }

        #[test]
        fn susceptible_taxgroup_detection_is_inlined_like_upstream() {
            let source = include_str!("../src/amrfinder.rs");
            assert!(!source.contains("fn susceptible_taxgroup_exists("));

            let susceptible_table = source
                .find("let tab = TextTable::from_file(&format!(\"{}/AMRProt-susceptible.tsv\", db))?;")
                .expect("susceptible table should be loaded inline like upstream");
            let slow_search = source[susceptible_table..]
                .find("run_tblastn_search(")
                .expect("susceptible slow tblastn should follow table probe");
            let block = &source[susceptible_table..susceptible_table + slow_search];

            assert!(block.contains("tab.qc()?;"));
            assert!(block.contains("let taxgroup_col = tab.col2num(\"taxgroup\")?;"));
            assert!(block.contains("let mut found = false;"));
            assert!(
                block.contains("ensure!(!row[taxgroup_col].is_empty(), \"QC assertion failed\");")
            );
            assert!(block.contains("if row[taxgroup_col] == organism1"));
            assert!(!block.contains(".exists()"));
        }

        #[test]
        fn amr_report_handoff_filters_suppress_table_by_organism_like_upstream() {
            let source = include_str!("../src/amrfinder.rs");
            let suppress_block = source
                .find("let suppress_source = PathBuf::from(format!(\"{}/AMRProt-suppress.tsv\", db));")
                .expect("translated suppress_prot creation should be present");
            let report_config = source[suppress_block..]
                .find("let report_config = AmrReportConfig")
                .expect("suppress handoff should precede amr_report config");
            let block = &source[suppress_block..suppress_block + report_config];

            assert!(block.contains("let suppress_filtered = tmp.join(\"suppress_prot\");"));
            assert!(block.contains("if line.starts_with('#')"));
            assert!(block.contains("if org == organism1"));
            assert!(block.contains("writeln!(suppress_out, \"{accver}\")?;"));
            assert!(block.contains("Some(suppress_filtered)"));
            assert!(!block.contains("Some(suppress_source)"));
        }

        #[test]
        fn db2organisms_reads_mutation_sources_sorts_and_deduplicates() {
            let tmp = tempfile::tempdir().unwrap();
            let db = tmp.path();
            write_organism_db(db);

            let organisms = db2organisms(db).unwrap();

            assert_eq!(
                organisms,
                vec![
                    "Acinetobacter".to_string(),
                    "Escherichia".to_string(),
                    "Escherichia_coli".to_string(),
                    "Klebsiella".to_string(),
                    "Salmonella".to_string(),
                ]
            );
        }

        fn write_organism_db(db: &Path) {
            write_database_version_files(db);
            fs::write(
                db.join("taxgroup.tsv"),
                concat!(
                    "taxgroup\tgpipe_org\n",
                    "Salmonella\tSalmonella enterica\n",
                    "Escherichia\tEscherichia coli\n",
                    "Escherichia_coli\tEscherichia coli str. K-12\n",
                ),
            )
            .unwrap();
            fs::write(
                db.join("AMRProt-mutation.tsv"),
                concat!(
                    "taxgroup\taccession\tgenesymbol\n",
                    "Klebsiella\tWP_1\tompK36\n",
                    "Escherichia\tWP_2\tgyrA\n",
                ),
            )
            .unwrap();
            fs::write(
                db.join("AMRProt-susceptible.tsv"),
                concat!(
                    "taxgroup\taccession\tgenesymbol\n",
                    "Acinetobacter\tWP_3\tcarO\n",
                    "Salmonella\tWP_4\tompC\n",
                    "\tWP_5\tblank_taxgroup_is_not_reported\n",
                ),
            )
            .unwrap();
        }

        fn write_gpipe_taxgroup_db(db: &Path) {
            write_database_version_files(db);
            fs::write(
                db.join("taxgroup.tsv"),
                concat!(
                    "taxgroup\tgpipe_taxgroup\n",
                    "Salmonella\tSalmonella enterica\n",
                    "Escherichia\tEscherichia coli, Escherichia coli str. K-12\n",
                ),
            )
            .unwrap();
            fs::write(
                db.join("AMRProt-mutation.tsv"),
                "taxgroup\taccession\tgenesymbol\nEscherichia\tWP_2\tgyrA\n",
            )
            .unwrap();
            fs::write(
                db.join("AMRProt-susceptible.tsv"),
                "taxgroup\taccession\tgenesymbol\nSalmonella\tWP_4\tompC\n",
            )
            .unwrap();
        }

        fn write_database_version_files(db: &Path) {
            fs::write(db.join("version.txt"), "2026-03-24.1\n").unwrap();
            fs::write(db.join("database_format_version.txt"), "0.1.0\n").unwrap();
            fs::write(db.join("AMRProt.fa.phr"), "").unwrap();
        }

        fn test_stxtyper_path() -> PathBuf {
            let stx_dir = std::env::current_exe()
                .unwrap()
                .parent()
                .unwrap()
                .join("stx");
            fs::create_dir_all(&stx_dir).unwrap();
            stx_dir.join("stxtyper")
        }

        #[test]
        fn normal_run_prints_database_metadata_before_search_mode_error_like_upstream() {
            let tmp = tempfile::tempdir().unwrap();
            write_database_version_files(tmp.path());

            let output = std::process::Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .arg("-d")
                .arg(tmp.path())
                .output()
                .unwrap();

            assert!(!output.status.success());
            let stderr = String::from_utf8(output.stderr).unwrap();
            let directory = stderr
                .find("Database directory: ")
                .expect("normal execution should print database directory");
            let version = stderr
                .find("Database version: 2026-03-24.1")
                .expect("normal execution should print database version");
            let search_mode = stderr
                .find("Parameter --protein or --nucleotide must be present")
                .expect("missing input should still fail before search execution");

            assert!(directory < version);
            assert!(version < search_mode);
        }

        #[test]
        fn list_organisms_prints_database_metadata_like_upstream() {
            let tmp = tempfile::tempdir().unwrap();
            write_organism_db(tmp.path());

            let output = std::process::Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .arg("-l")
                .arg("-d")
                .arg(tmp.path())
                .output()
                .unwrap();

            assert!(output.status.success());
            let stderr = String::from_utf8(output.stderr).unwrap();
            assert!(stderr.contains("Database directory: "));
            assert!(stderr.contains("Database version: 2026-03-24.1"));
            let stdout = String::from_utf8(output.stdout).unwrap();
            assert!(stdout.contains(
                "Available --organism options: Acinetobacter, Escherichia, Escherichia_coli, Klebsiella, Salmonella"
            ));
        }

        #[test]
        fn post_report_raw_disruptions_are_replaced_in_element_symbols() {
            let tmp = tempfile::tempdir().unwrap();
            let db = tmp.path().join("db");
            fs::create_dir(&db).unwrap();
            let dna = tmp.path().join("dna.fa");
            fs::write(&dna, ">contig1\nATGTAA\n").unwrap();
            fs::write(db.join("AMRProt-susceptible.fa"), ">WP_1\nM\n").unwrap();
            let output = concat!(
                "Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\tSubtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference\t% Identity to reference\tAlignment length\tClosest reference accession\tClosest reference name\tHMM accession\tHMM description\tHierarchy node\n",
                "NA\tcontig1\t1\t6\t+\tcirA_@fs_0_1_0_6_1_STOP\tKlebsiella pneumoniae cefiderocol resistant CirA\tcore\tAMR\tPOINT_DISRUPT\tBETA-LACTAM\tCEFIDEROCOL\tPOINTX\t1\t1\t100.00\t99.85\t1\tWP_1\tcatecholate siderophore receptor CirA\tNA\tNA\tcirA\n",
                "NA\tcontig\t1\t9\t+\tblaTEM\tTEM family class A beta-lactamase\tcore\tAMR\tAMR\tBETA-LACTAM\tBETA-LACTAM\tBLAST\t286\t286\t100.00\t100.00\t286\tWP_1\tTEM\tNA\tNA\tblaTEM\n",
            )
            .to_string();
            let mut table = text_table_from_tsv(&output).unwrap();

            amr_tab_disruptions(&mut table, db.to_str().unwrap(), Some(&dna), 11).unwrap();
            let output = text_table_to_tsv(&table).unwrap();

            assert!(output
                .contains("\tcirA_M1MfsTer1\tKlebsiella pneumoniae cefiderocol resistant CirA\t"));
            assert!(output.contains("\tblaTEM\tTEM family class A beta-lactamase\t"));
        }

        #[test]
        fn disruption_postprocessing_requires_original_reference_column() {
            let tmp = tempfile::tempdir().unwrap();
            let db = tmp.path().join("db");
            fs::create_dir(&db).unwrap();
            let dna = tmp.path().join("dna.fa");
            fs::write(&dna, ">contig1\nATGTAA\n").unwrap();
            fs::write(db.join("AMRProt-susceptible.fa"), ">WP_1\nM\n").unwrap();
            let mut table = text_table_from_tsv(concat!(
                "Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\n",
                "NA\tcontig1\t1\t6\t+\tcirA_@fs_0_1_0_6_1_STOP\tCirA\n",
            ))
            .unwrap();

            let err = amr_tab_disruptions(&mut table, db.to_str().unwrap(), Some(&dna), 11)
                .unwrap_err()
                .to_string();

            assert!(err.contains(crate::columns::CLOSEST_REF_ACCESSION_COL_NAME));
        }

        #[test]
        fn prepare_fasta_extract_protein_targets_drop_na_sort_and_dedup() {
            let table = text_table_from_tsv(&minimal_report_tsv(&[
                "prot2\tNA\t0\t0\t+\tblaTEM\tTEM family class A beta-lactamase\tcore\tAMR\tAMR",
                "NA\tcontig\t1\t9\t+\tgyrA\tDNA gyrase\tcore\tAMR\tPOINT",
                "prot1\tNA\t0\t0\t+\taac\taminoglycoside acetyltransferase\tcore\tAMR\tAMR",
                "prot1\tNA\t0\t0\t+\taac\taminoglycoside acetyltransferase\tcore\tAMR\tAMR",
            ]))
            .unwrap();

            let targets = prepare_fasta_extract(
                &table,
                &[
                    crate::columns::PROT_COL_NAME,
                    crate::columns::GENESYMBOL_COL_NAME,
                    crate::columns::ELEM_NAME_COL_NAME,
                ],
                false,
            )
            .unwrap();
            let tsv = text_table_to_tsv(&targets).unwrap();

            assert_eq!(
                tsv,
                concat!(
                    "prot1\taac\taminoglycoside acetyltransferase\n",
                    "prot2\tblaTEM\tTEM family class A beta-lactamase\n",
                )
            );
        }

        #[test]
        fn prepare_fasta_extract_nucleotide_targets_can_keep_header_for_flank_edit() {
            let table = text_table_from_tsv(&minimal_report_tsv(&[
                "NA\tcontig2\t20\t30\t-\tgyrA\tDNA gyrase\tcore\tAMR\tPOINT",
                "prot1\tNA\t0\t0\t+\tblaTEM\tTEM family class A beta-lactamase\tcore\tAMR\tAMR",
                "NA\tcontig1\t1\t9\t+\taac\taminoglycoside acetyltransferase\tcore\tAMR\tAMR",
            ]))
            .unwrap();

            let targets = prepare_fasta_extract(
                &table,
                &[
                    crate::columns::CONTIG_COL_NAME,
                    crate::columns::START_COL_NAME,
                    crate::columns::STOP_COL_NAME,
                    crate::columns::STRAND_COL_NAME,
                    crate::columns::GENESYMBOL_COL_NAME,
                    crate::columns::ELEM_NAME_COL_NAME,
                ],
                true,
            )
            .unwrap();
            let tsv = text_table_to_tsv(&targets).unwrap();

            assert_eq!(
                tsv,
                concat!(
                    "Contig id\tStart\tStop\tStrand\tElement symbol\tElement name\n",
                    "contig1\t1\t9\t+\taac\taminoglycoside acetyltransferase\n",
                    "contig2\t20\t30\t-\tgyrA\tDNA gyrase\n",
                )
            );
        }

        #[test]
        fn report_fasta_outputs_use_translated_fasta_extract_targets() {
            let tmp = tempfile::tempdir().unwrap();
            let protein = tmp.path().join("prot.fa");
            let nucleotide = tmp.path().join("dna.fa");
            let protein_output = tmp.path().join("prot_out.fa");
            let nucleotide_output = tmp.path().join("dna_out.fa");
            let flank_output = tmp.path().join("flank_out.fa");
            fs::write(&protein, ">prot1 description\nMKK\n>prot2\nAAA\n").unwrap();
            fs::write(&nucleotide, ">contig1\nAACGTTA\n>contig2\nTTTTTTTT\n").unwrap();
            let config = PipelineConfig {
                protein: Some(protein),
                nucleotide: Some(nucleotide),
                protein_output: Some(protein_output.clone()),
                nucleotide_output: Some(nucleotide_output.clone()),
                nucleotide_flank5_output: Some(flank_output.clone()),
                nucleotide_flank5_size: 2,
                ..PipelineConfig::default()
            };
            let report = minimal_report_tsv(&[
                "prot1\tNA\t0\t0\t+\tblaTEM\tTEM family class A beta-lactamase\tcore\tAMR\tAMR",
                "NA\tcontig1\t3\t5\t+\tgyrA\tDNA gyrase\tcore\tAMR\tPOINT",
            ]);

            write_report_fasta_outputs(
                &report,
                &config,
                tmp.path(),
                config.protein.as_deref(),
                config.nucleotide.as_deref(),
            )
            .unwrap();

            assert_eq!(
                fs::read_to_string(protein_output).unwrap(),
                ">prot1 blaTEM TEM family class A beta-lactamase\nMKK\n"
            );
            assert_eq!(
                fs::read_to_string(nucleotide_output).unwrap(),
                ">contig1:3-5 strand:+ gyrA DNA gyrase\nCGT\n"
            );
            assert_eq!(
                fs::read_to_string(flank_output).unwrap(),
                ">contig1:1-5 strand:+ gyrA DNA gyrase\nAACGT\n"
            );
        }

        #[test]
        fn run_pipeline_rejects_output_fasta_options_without_required_inputs() {
            let tmp = tempfile::tempdir().unwrap();
            write_database_version_files(tmp.path());
            let err = run_pipeline(&PipelineConfig {
                nucleotide: Some(PathBuf::from("dna.fa")),
                protein_output: Some(PathBuf::from("prot_out.fa")),
                database: tmp.path().to_path_buf(),
                ..PipelineConfig::default()
            })
            .unwrap_err();
            assert!(err
                .to_string()
                .contains("Parameter --protein must be present for --protein_output"));

            let err = run_pipeline(&PipelineConfig {
                protein: Some(PathBuf::from("prot.fa")),
                nucleotide_flank5_size: 10,
                database: tmp.path().to_path_buf(),
                ..PipelineConfig::default()
            })
            .unwrap_err();
            assert!(err.to_string().contains(
                "Parameter --nucleotide_flank5_output must be present for --nucleotide_flank5_size"
            ));
        }

        #[test]
        fn run_pipeline_validates_database_before_search_mode_options_like_upstream() {
            let err = run_pipeline(&PipelineConfig {
                nucleotide: Some(PathBuf::from("dna.fa")),
                protein_output: Some(PathBuf::from("prot_out.fa")),
                database: PathBuf::from("/missing/amrfinder/db/latest"),
                ..PipelineConfig::default()
            })
            .unwrap_err();

            assert!(err.to_string().contains(
                "No valid AMRFinder database is found.\nThis directory (or symbolic link to directory) is not found: /missing/amrfinder/db/latest"
            ));

            let err = run_pipeline(&PipelineConfig {
                database: PathBuf::from("/missing/amrfinder/db/latest"),
                ..PipelineConfig::default()
            })
            .unwrap_err();

            assert!(err.to_string().contains(
                "No valid AMRFinder database is found.\nThis directory (or symbolic link to directory) is not found: /missing/amrfinder/db/latest"
            ));
            assert!(!err
                .to_string()
                .contains("Parameter --protein or --nucleotide must be present"));

            let err = run_pipeline(&PipelineConfig {
                protein: Some(PathBuf::from("prot.fa")),
                pgap: true,
                annotation_format: "standard".to_string(),
                database: PathBuf::from("/missing/amrfinder/db/latest"),
                ..PipelineConfig::default()
            })
            .unwrap_err();

            assert!(err.to_string().contains(
                "No valid AMRFinder database is found.\nThis directory (or symbolic link to directory) is not found: /missing/amrfinder/db/latest"
            ));
            assert!(!err.to_string().contains("--pgap conflicts"));
        }

        #[test]
        fn run_pipeline_validates_organism_before_pgap_conflict_like_upstream() {
            let tmp = tempfile::tempdir().unwrap();
            write_organism_db(tmp.path());

            let err = run_pipeline(&PipelineConfig {
                protein: Some(PathBuf::from("prot.fa")),
                database: tmp.path().to_path_buf(),
                organism: "Missing".to_string(),
                pgap: true,
                annotation_format: "standard".to_string(),
                ..PipelineConfig::default()
            })
            .unwrap_err();

            assert!(err.to_string().contains("Possible organisms:"));
            assert!(err.to_string().contains("Escherichia"));
            assert!(!err.to_string().contains("--pgap conflicts"));
        }

        #[test]
        fn pipeline_fasta_check_uses_translated_body_for_protein_and_dna() {
            let tmp = tempfile::tempdir().unwrap();
            let prot = tmp.path().join("prot.fa");
            let dna = tmp.path().join("dna.fa");
            let len = tmp.path().join("dna.len");
            let cleaned = tmp.path().join("clean.fa");
            fs::write(&prot, ">p1\nMK*\n>p2\nAA\n").unwrap();
            fs::write(&dna, ">c1\nACG-TN\n").unwrap();

            assert_eq!(
                fasta_check(&prot, true, 20, None, Some(&cleaned)).unwrap(),
                (2, 3, 5)
            );
            assert_eq!(fs::read_to_string(&cleaned).unwrap(), ">p1\nMK*\n>p2\nAA\n");

            assert_eq!(
                fasta_check(&dna, false, 0, Some(&len), None).unwrap(),
                (1, 6, 6)
            );
            assert_eq!(fs::read_to_string(&len).unwrap(), "c1\t6\n");
        }

        #[test]
        fn protein_fasta_check_recovery_cleans_fixable_inputs() {
            let tmp = tempfile::tempdir().unwrap();
            let hyphen = tmp.path().join("hyphen.fa");
            let hyphen_clean = tmp.path().join("hyphen.clean.fa");
            fs::write(&hyphen, ">p1\nM-K\n").unwrap();

            let err = fasta_check(&hyphen, true, 20, None, None).unwrap_err();
            assert!(err.to_string().contains("Hyphen in the sequence"));
            assert_eq!(
                fasta_check(&hyphen, true, 20, None, Some(&hyphen_clean)).unwrap(),
                (1, 2, 2)
            );
            assert_eq!(fs::read_to_string(&hyphen_clean).unwrap(), ">p1\nMK\n");

            let ambig = tmp.path().join("ambig.fa");
            let ambig_clean = tmp.path().join("ambig.clean.fa");
            fs::write(&ambig, format!(">bad\n{}\n>good\nMKQ\n", "X".repeat(21))).unwrap();

            let err = fasta_check(&ambig, true, 20, None, None).unwrap_err();
            assert!(err.to_string().contains("Too many ambiguities"));
            assert_eq!(
                fasta_check(&ambig, true, 20, None, Some(&ambig_clean)).unwrap(),
                (2, 3, 3)
            );
            assert_eq!(fs::read_to_string(&ambig_clean).unwrap(), ">good\nMKQ\n");
        }

        #[test]
        fn run_pipeline_rejects_upstream_option_validation_before_database_work() {
            let err = run_pipeline(&PipelineConfig {
                name: "bad\tname".to_string(),
                ..PipelineConfig::default()
            })
            .unwrap_err();
            assert!(err
                .to_string()
                .contains("NAME cannot contain a tab character"));

            let err = run_pipeline(&PipelineConfig {
                protein: Some(PathBuf::from("prot.fa")),
                ident_min: 1.5,
                ..PipelineConfig::default()
            })
            .unwrap_err();
            assert!(err
                .to_string()
                .contains("ident_min must be between 0 and 1"));

            let err = run_pipeline(&PipelineConfig {
                nucleotide: Some(PathBuf::from("dna.fa")),
                coverage_min: -0.1,
                ..PipelineConfig::default()
            })
            .unwrap_err();
            assert!(err
                .to_string()
                .contains("coverage_min must be between 0 and 1"));

            let err = run_pipeline(&PipelineConfig {
                protein: Some(PathBuf::from("prot.fa")),
                report_common: true,
                ..PipelineConfig::default()
            })
            .unwrap_err();
            assert!(err
                .to_string()
                .contains("--report_common requires --organism"));

            let err = run_pipeline(&PipelineConfig {
                protein: Some(PathBuf::from("prot.fa")),
                organism: "Escherichia".to_string(),
                report_common: true,
                ..PipelineConfig::default()
            })
            .unwrap_err();
            assert!(err.to_string().contains("--report_common requires --plus"));
        }

        #[test]
        fn run_pipeline_rejects_upstream_search_mode_option_combinations() {
            let tmp = tempfile::tempdir().unwrap();
            write_database_version_files(tmp.path());

            let err = run_pipeline(&PipelineConfig {
                database: tmp.path().to_path_buf(),
                ..PipelineConfig::default()
            })
            .unwrap_err();
            assert!(err
                .to_string()
                .contains("Parameter --protein or --nucleotide must be present"));

            let err = run_pipeline(&PipelineConfig {
                nucleotide: Some(PathBuf::from("dna.fa")),
                gff: Some(PathBuf::from("ann.gff")),
                database: tmp.path().to_path_buf(),
                ..PipelineConfig::default()
            })
            .unwrap_err();
            assert!(err.to_string().contains("Parameter --gff is redundant"));

            let err = run_pipeline(&PipelineConfig {
                protein: Some(PathBuf::from("prot.fa")),
                nucleotide: Some(PathBuf::from("dna.fa")),
                database: tmp.path().to_path_buf(),
                ..PipelineConfig::default()
            })
            .unwrap_err();
            assert!(err.to_string().contains(
                "If parameters --protein and --nucleotide are present then parameter --gff must be present"
            ));
        }

        #[test]
        fn run_pipeline_validates_gff_against_protein_before_blast_like_upstream() {
            let tmp = tempfile::tempdir().unwrap();
            let db = tmp.path().join("db");
            fs::create_dir(&db).unwrap();
            write_database_version_files(&db);
            let protein = tmp.path().join("prot.fa");
            let gff = tmp.path().join("annot.gff");
            fs::write(&protein, ">prot1\nMAA\n").unwrap();
            fs::write(
                &gff,
                "ctg\tRefSeq\tCDS\t1\t3\t.\t+\t0\tName=other;gene=g;product=p\n",
            )
            .unwrap();

            let err = run_pipeline(&PipelineConfig {
                protein: Some(protein),
                gff: Some(gff),
                database: db,
                annotation_format: "standard".to_string(),
                ..PipelineConfig::default()
            })
            .unwrap_err();

            assert!(err.to_string().starts_with("GFF file mismatch.\n"));
            assert!(err
                .to_string()
                .contains("Protein FASTA id \"prot1\" is not in the GFF file"));
            assert!(!err.to_string().contains("blastp failed"));
        }

        #[test]
        fn run_pipeline_missing_database_uses_upstream_download_message() {
            let err = run_pipeline(&PipelineConfig {
                protein: Some(PathBuf::from("prot.fa")),
                database: PathBuf::from("/missing/amrfinder/db/latest"),
                ..PipelineConfig::default()
            })
            .unwrap_err();

            assert!(err.to_string().contains(
                "No valid AMRFinder database is found.\nThis directory (or symbolic link to directory) is not found: /missing/amrfinder/db/latest\nTo download the latest version to the default directory run: amrfinder -u"
            ));
        }

        #[test]
        fn run_pipeline_requires_amrprot_blast_database_marker_before_inputs_like_upstream() {
            let tmp = tempfile::tempdir().unwrap();
            fs::write(tmp.path().join("version.txt"), "2026-03-24.1\n").unwrap();
            fs::write(tmp.path().join("database_format_version.txt"), "0.1.0\n").unwrap();

            let err = run_pipeline(&PipelineConfig {
                protein: Some(PathBuf::from("prot.fa")),
                database: tmp.path().to_path_buf(),
                ..PipelineConfig::default()
            })
            .unwrap_err();

            assert!(err.to_string().contains(&format!(
                "The BLAST database for AMRProt.fa ({}) was not found.",
                tmp.path().join("AMRProt.fa.phr").display()
            )));
            assert!(!err
                .to_string()
                .contains("Protein FASTA file is empty or missing"));
        }

        #[test]
        fn run_pipeline_allows_empty_protein_file_like_upstream() {
            let tmp = tempfile::tempdir().unwrap();
            let db = tmp.path().join("db");
            fs::create_dir(&db).unwrap();
            write_database_version_files(&db);
            fs::write(
                db.join("fam.tsv"),
                "fam0\t\tgene0\t-\t0\t0\t90\t0\t90\t90\t0\t50\t1\tAMR\tAMR\tCLASS\tSUBCLASS\tFamily 0\n",
            )
            .unwrap();
            fs::write(db.join("AMRProt-mutation.tsv"), "# no mutation rows\n").unwrap();
            fs::write(
                db.join("AMRProt-susceptible.tsv"),
                "# no susceptible rows\n",
            )
            .unwrap();
            let protein = tmp.path().join("empty.fa");
            fs::write(&protein, "").unwrap();

            let output = run_pipeline(&PipelineConfig {
                protein: Some(protein),
                database: db,
                plus: true,
                ..PipelineConfig::default()
            })
            .unwrap();

            assert!(output.starts_with("Protein id\tElement symbol\tElement name\t"));
            assert_eq!(output.matches('\n').count(), 1);
        }

        #[test]
        fn run_pipeline_allows_empty_nucleotide_file_like_upstream() {
            let tmp = tempfile::tempdir().unwrap();
            let db = tmp.path().join("db");
            fs::create_dir(&db).unwrap();
            write_database_version_files(&db);
            fs::write(
                db.join("fam.tsv"),
                "fam0\t\tgene0\t-\t0\t0\t90\t0\t90\t90\t0\t50\t1\tAMR\tAMR\tCLASS\tSUBCLASS\tFamily 0\n",
            )
            .unwrap();
            fs::write(db.join("AMRProt-mutation.tsv"), "# no mutation rows\n").unwrap();
            fs::write(
                db.join("AMRProt-susceptible.tsv"),
                "# no susceptible rows\n",
            )
            .unwrap();
            fs::write(db.join("AMRProt-suppress.tsv"), "# no suppress rows\n").unwrap();
            let nucleotide = tmp.path().join("empty.fa");
            fs::write(&nucleotide, "").unwrap();

            let output = run_pipeline(&PipelineConfig {
                nucleotide: Some(nucleotide),
                database: db,
                plus: true,
                ..PipelineConfig::default()
            })
            .unwrap();

            assert!(output.starts_with("Protein id\tContig id\tStart\tStop\tStrand\t"));
            assert_eq!(output.matches('\n').count(), 1);
        }

        #[test]
        fn normal_run_prints_search_mode_and_include_hints_like_upstream() {
            let tmp = tempfile::tempdir().unwrap();
            let db = tmp.path().join("db");
            fs::create_dir(&db).unwrap();
            write_database_version_files(&db);
            fs::write(
                db.join("fam.tsv"),
                "fam0\t\tgene0\t-\t0\t0\t90\t0\t90\t90\t0\t50\t1\tAMR\tAMR\tCLASS\tSUBCLASS\tFamily 0\n",
            )
            .unwrap();
            fs::write(db.join("AMRProt-mutation.tsv"), "# no mutation rows\n").unwrap();
            fs::write(
                db.join("AMRProt-susceptible.tsv"),
                "# no susceptible rows\n",
            )
            .unwrap();
            fs::write(db.join("AMRProt-suppress.tsv"), "# no suppress rows\n").unwrap();
            let nucleotide = tmp.path().join("empty.fa");
            fs::write(&nucleotide, "").unwrap();

            let output = std::process::Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .arg("-n")
                .arg(&nucleotide)
                .arg("-d")
                .arg(&db)
                .arg("--plus")
                .output()
                .unwrap();

            assert!(output.status.success());
            let stderr = String::from_utf8(output.stderr).unwrap();
            assert!(stderr.contains("AMRFinder translated nucleotide search"));
            assert!(stderr.contains(
                "  - include -O/--organism option to add mutation searches and suppress common proteins"
            ));
        }

        #[test]
        fn run_pipeline_uncompresses_gzipped_nucleotide_input_like_upstream() {
            let tmp = tempfile::tempdir().unwrap();
            let db = tmp.path().join("db");
            fs::create_dir(&db).unwrap();
            write_database_version_files(&db);
            fs::write(
                db.join("fam.tsv"),
                "fam0\t\tgene0\t-\t0\t0\t90\t0\t90\t90\t0\t50\t1\tAMR\tAMR\tCLASS\tSUBCLASS\tFamily 0\n",
            )
            .unwrap();
            fs::write(db.join("AMRProt-mutation.tsv"), "# no mutation rows\n").unwrap();
            fs::write(
                db.join("AMRProt-susceptible.tsv"),
                "# no susceptible rows\n",
            )
            .unwrap();
            fs::write(db.join("AMRProt-suppress.tsv"), "# no suppress rows\n").unwrap();
            let empty = tmp.path().join("empty.fa");
            let gz = tmp.path().join("empty.fa.gz");
            fs::write(&empty, "").unwrap();
            let gzip = std::process::Command::new("gzip")
                .arg("-c")
                .arg(&empty)
                .output()
                .unwrap();
            assert!(gzip.status.success());
            fs::write(&gz, gzip.stdout).unwrap();

            let output = run_pipeline(&PipelineConfig {
                nucleotide: Some(gz),
                database: db,
                plus: true,
                ..PipelineConfig::default()
            })
            .unwrap();

            assert!(output.starts_with("Protein id\tContig id\tStart\tStop\tStrand\t"));
            assert_eq!(output.matches('\n').count(), 1);
        }

        #[test]
        fn empty_gzipped_nucleotide_warning_reports_original_path_like_upstream() {
            let tmp = tempfile::tempdir().unwrap();
            let db = tmp.path().join("db");
            fs::create_dir(&db).unwrap();
            write_database_version_files(&db);
            fs::write(
                db.join("fam.tsv"),
                "fam0\t\tgene0\t-\t0\t0\t90\t0\t90\t90\t0\t50\t1\tAMR\tAMR\tCLASS\tSUBCLASS\tFamily 0\n",
            )
            .unwrap();
            fs::write(db.join("AMRProt-mutation.tsv"), "# no mutation rows\n").unwrap();
            fs::write(
                db.join("AMRProt-susceptible.tsv"),
                "# no susceptible rows\n",
            )
            .unwrap();
            fs::write(db.join("AMRProt-suppress.tsv"), "# no suppress rows\n").unwrap();
            let empty = tmp.path().join("empty.fa");
            let gz = tmp.path().join("empty.fa.gz");
            fs::write(&empty, "").unwrap();
            let gzip = std::process::Command::new("gzip")
                .arg("-c")
                .arg(&empty)
                .output()
                .unwrap();
            assert!(gzip.status.success());
            fs::write(&gz, gzip.stdout).unwrap();

            let output = std::process::Command::new(env!("CARGO_BIN_EXE_amrfinder"))
                .arg("-n")
                .arg(&gz)
                .arg("-d")
                .arg(&db)
                .arg("--plus")
                .output()
                .unwrap();

            assert!(output.status.success());
            let stderr = String::from_utf8(output.stderr).unwrap();
            assert!(stderr.contains(&format!("WARNING: Empty file: {}", gz.display())));
            assert!(!stderr.contains("dna_flat"));
        }

        #[test]
        fn escherichia_plus_without_dna_database_does_not_require_stxtyper_like_upstream() {
            let tmp = tempfile::tempdir().unwrap();
            let db = tmp.path().join("db");
            fs::create_dir(&db).unwrap();
            write_database_version_files(&db);
            fs::write(
                db.join("fam.tsv"),
                "fam0\t\tgene0\t-\t0\t0\t90\t0\t90\t90\t0\t50\t1\tAMR\tAMR\tCLASS\tSUBCLASS\tFamily 0\n",
            )
            .unwrap();
            fs::write(
                db.join("taxgroup.tsv"),
                "taxgroup\tgpipe_taxgroup\nEscherichia\tEscherichia coli\n",
            )
            .unwrap();
            fs::write(db.join("AMRProt-mutation.tsv"), "# no mutation rows\n").unwrap();
            fs::write(
                db.join("AMRProt-susceptible.tsv"),
                "# no susceptible rows\n",
            )
            .unwrap();
            fs::write(db.join("AMRProt-suppress.tsv"), "# no suppress rows\n").unwrap();
            let nucleotide = tmp.path().join("empty.fa");
            fs::write(&nucleotide, "").unwrap();

            let output = run_pipeline(&PipelineConfig {
                nucleotide: Some(nucleotide),
                database: db,
                organism: "Escherichia".to_string(),
                plus: true,
                ..PipelineConfig::default()
            })
            .unwrap();

            assert!(output.starts_with("Protein id\tContig id\tStart\tStop\tStrand\t"));
            assert_eq!(output.matches('\n').count(), 1);
        }

        #[cfg(unix)]
        #[test]
        fn stxtyper_version_is_checked_before_empty_nucleotide_report_like_upstream() {
            use std::os::unix::fs::PermissionsExt;

            let _guard = ENV_LOCK.lock().unwrap();
            let tmp = tempfile::tempdir().unwrap();
            let db = tmp.path().join("db");
            fs::create_dir(&db).unwrap();
            write_database_version_files(&db);
            fs::write(
                db.join("taxgroup.tsv"),
                "taxgroup\tgpipe_taxgroup\nEscherichia\tEscherichia coli\n",
            )
            .unwrap();
            fs::write(db.join("AMRProt-mutation.tsv"), "# no mutation rows\n").unwrap();
            fs::write(
                db.join("AMRProt-susceptible.tsv"),
                "# no susceptible rows\n",
            )
            .unwrap();
            fs::write(db.join("AMR_DNA-Escherichia.fa"), ">ref\nATGC\n").unwrap();
            fs::write(db.join("AMR_DNA-Escherichia.fa.ndb"), "").unwrap();
            let stxtyper = test_stxtyper_path();
            fs::write(
                &stxtyper,
                "#!/bin/sh\nif [ \"$1\" = \"-v\" ]; then echo 0.0.0; exit 0; fi\nexit 1\n",
            )
            .unwrap();
            let mut perms = fs::metadata(&stxtyper).unwrap().permissions();
            perms.set_mode(0o755);
            fs::set_permissions(&stxtyper, perms).unwrap();
            let nucleotide = tmp.path().join("empty.fa");
            fs::write(&nucleotide, "").unwrap();

            let err = run_pipeline(&PipelineConfig {
                nucleotide: Some(nucleotide),
                database: db,
                organism: "Escherichia".to_string(),
                plus: true,
                ..PipelineConfig::default()
            })
            .unwrap_err();

            assert_eq!(
                err.to_string(),
                "AMRFinder invokes StxTyper version 0.0.0. Expected StxTyper version is 1.0.45"
            );
        }

        #[test]
        fn empty_nucleotide_branch_creates_empty_blastn_when_dna_database_exists_like_upstream() {
            let source = include_str!("../src/amrfinder.rs");
            let nucleotide_empty_branch = source
                .find("if nucleotide_empty {")
                .expect("translated empty nucleotide branch should be present");
            let non_empty_branch = source[nucleotide_empty_branch..]
                .find("} else {")
                .expect("empty nucleotide branch should have a non-empty counterpart");
            let branch =
                &source[nucleotide_empty_branch..nucleotide_empty_branch + non_empty_branch];

            assert!(branch.contains("File::create(&blastx_out)?;"));
            assert!(branch.contains("File::create(&dna_len_path)?;"));
            assert!(branch.contains("if blastn"));
            assert!(branch.contains("File::create(tmp.join(\"blastn\"))?;"));
        }

        #[test]
        fn susceptible_tblastn_handoff_requires_organism_like_upstream() {
            let source = include_str!("../src/amrfinder.rs");
            let susceptible_probe = source
                .find("let tab = TextTable::from_file(&format!(\"{}/AMRProt-susceptible.tsv\", db))?;")
                .expect("translated susceptible search probe should be present");
            let branch = &source[..susceptible_probe];
            let guard = branch
                .rfind("!organism1.is_empty()")
                .expect("susceptible search should be guarded by a non-empty organism");
            let blastx_choice = branch
                .rfind("let use_tblastn = dna_len_max > 100_000;")
                .expect("susceptible search should remain in the non-empty DNA branch");

            assert!(blastx_choice < guard);
        }

        #[test]
        fn run_pipeline_validates_annotation_format_before_database_work() {
            let err = run_pipeline(&PipelineConfig {
                protein: Some(PathBuf::from("prot.fa")),
                annotation_format: "not_a_gff_type".to_string(),
                database: PathBuf::from("/missing/amrfinder/db/latest"),
                ..PipelineConfig::default()
            })
            .unwrap_err();

            assert!(err
                .to_string()
                .contains("Unknown GFF type: \"not_a_gff_type\""));
        }

        #[test]
        fn run_pipeline_translates_pgap_flag_gff_type_slice() {
            let tmp = tempfile::tempdir().unwrap();
            write_database_version_files(tmp.path());

            let err = run_pipeline(&PipelineConfig {
                protein: Some(PathBuf::from("prot.fa")),
                pgap: true,
                annotation_format: "standard".to_string(),
                database: tmp.path().to_path_buf(),
                ..PipelineConfig::default()
            })
            .unwrap_err();
            assert!(err
                .to_string()
                .contains("--pgap conflicts with GFF type \"standard\""));

            let err = run_pipeline(&PipelineConfig {
                protein: Some(PathBuf::from("prot.fa")),
                pgap: true,
                database: PathBuf::from("/missing/amrfinder/db/latest"),
                ..PipelineConfig::default()
            })
            .unwrap_err();
            assert!(err.to_string().contains(
                "No valid AMRFinder database is found.\nThis directory (or symbolic link to directory) is not found: /missing/amrfinder/db/latest"
            ));
        }

        #[test]
        fn pgap_lcl_nucleotide_headers_are_handed_to_report_like_upstream() {
            let source = include_str!("../src/amrfinder.rs");
            let lcl_probe = source
                .find("line.starts_with(\">lcl|\")")
                .expect("PGAP lcl nucleotide header probe should be present");
            let report_handoff = source[lcl_probe..]
                .find("amr_report_blastp.push_str(\" -lcl\")")
                .expect("PGAP lcl state should be handed to the report stage");
            let report_parser = source[lcl_probe + report_handoff..]
                .find("\"-lcl\"")
                .expect("report handoff parser should accept the lcl marker");
            let normalizer = source[lcl_probe + report_handoff + report_parser..]
                .find("normalized.push_str(\"lcl|\")")
                .expect("report GFF normalization should prefix PGAP contigs with lcl");

            assert!(report_handoff > 0);
            assert!(report_parser > 0);
            assert!(normalizer > 0);
        }

        #[test]
        fn run_pipeline_validates_organism_against_translated_db_sources_before_blast() {
            let tmp = tempfile::tempdir().unwrap();
            write_organism_db(tmp.path());

            let err = run_pipeline(&PipelineConfig {
                protein: Some(PathBuf::from("prot.fa")),
                database: tmp.path().to_path_buf(),
                organism: "MissingTaxgroup".to_string(),
                ..PipelineConfig::default()
            })
            .unwrap_err();

            assert!(err.to_string().contains(
                "Possible organisms: Acinetobacter, Escherichia, Escherichia_coli, Klebsiella, Salmonella"
            ));
        }

        #[test]
        fn run_pipeline_normalizes_organism_spaces_like_original_organism1() {
            let tmp = tempfile::tempdir().unwrap();
            write_organism_db(tmp.path());

            let err = run_pipeline(&PipelineConfig {
                protein: Some(PathBuf::from("prot.fa")),
                database: tmp.path().to_path_buf(),
                organism: "Escherichia coli".to_string(),
                ..PipelineConfig::default()
            })
            .unwrap_err();

            assert!(err
                .to_string()
                .contains("Protein FASTA file is empty or missing: prot.fa"));
        }

        #[test]
        fn run_pipeline_maps_gpipe_taxgroup_before_organism_validation() {
            let tmp = tempfile::tempdir().unwrap();
            write_gpipe_taxgroup_db(tmp.path());

            let err = run_pipeline(&PipelineConfig {
                protein: Some(PathBuf::from("prot.fa")),
                database: tmp.path().to_path_buf(),
                organism: "Escherichia coli str. K-12".to_string(),
                gpipe_org: true,
                ..PipelineConfig::default()
            })
            .unwrap_err();

            assert!(err
                .to_string()
                .contains("Protein FASTA file is empty or missing: prot.fa"));
        }

        #[test]
        fn run_pipeline_missing_gpipe_taxgroup_clears_organism_like_upstream() {
            let tmp = tempfile::tempdir().unwrap();
            write_gpipe_taxgroup_db(tmp.path());

            let err = run_pipeline(&PipelineConfig {
                protein: Some(PathBuf::from("prot.fa")),
                database: tmp.path().to_path_buf(),
                organism: "Unknown GPipe".to_string(),
                gpipe_org: true,
                ..PipelineConfig::default()
            })
            .unwrap_err();

            assert!(err
                .to_string()
                .contains("Protein FASTA file is empty or missing: prot.fa"));
            assert!(!err.to_string().contains("Possible organisms:"));
        }

        #[test]
        fn run_pipeline_checks_dna_mutation_blast_database_markers_like_upstream() {
            let tmp = tempfile::tempdir().unwrap();
            write_organism_db(tmp.path());
            let dna = tmp.path().join("query.fa");
            fs::write(&dna, ">contig\nATGC\n").unwrap();
            fs::write(tmp.path().join("AMR_DNA-Escherichia.fa"), ">ref\nATGC\n").unwrap();

            let err = run_pipeline(&PipelineConfig {
                nucleotide: Some(dna),
                database: tmp.path().to_path_buf(),
                organism: "Escherichia".to_string(),
                ..PipelineConfig::default()
            })
            .unwrap_err();

            assert!(err.to_string().contains(
                "The BLAST database for AMR_DNA-Escherichia.fa was not found.\nUse amrfinder -u or amrfinder --force_update to download and prepare database for AMRFinderPlus"
            ));
        }

        #[test]
        fn amr_report_handoff_runs_dna_mutation_only_when_blastn_gate_is_true() {
            let tmp = tempfile::tempdir().unwrap();
            let db = tmp.path().join("db");
            fs::create_dir(&db).unwrap();
            fs::write(
                db.join("fam.tsv"),
                "fam0\t\tgene0\t-\t0\t0\t90\t0\t90\t90\t0\t50\t1\tAMR\tAMR\tCLASS\tSUBCLASS\tFamily 0\n",
            )
            .unwrap();
            fs::write(db.join("AMRProt-mutation.tsv"), "# no mutation rows\n").unwrap();
            fs::write(
                db.join("AMRProt-susceptible.tsv"),
                "# no susceptible rows\n",
            )
            .unwrap();
            fs::write(db.join("AMRProt-suppress.tsv"), "# no suppress rows\n").unwrap();
            fs::write(tmp.path().join("blastn"), "").unwrap();
            fs::write(db.join("AMR_DNA-Escherichia.tsv"), "mutation_id\n").unwrap();
            let config = PipelineConfig {
                nucleotide: Some(PathBuf::from("dna.fa")),
                database: db.clone(),
                organism: "Escherichia".to_string(),
                plus: true,
                ..PipelineConfig::default()
            };

            let (output, _fasta_source) = run_rust_amr_report(
                &config,
                tmp.path(),
                db.to_str().unwrap(),
                "Escherichia",
                "",
                "",
                None,
                false,
            )
            .unwrap();

            assert!(output.starts_with("Protein id\tContig id\tStart\tStop\tStrand\t"));
            assert_eq!(output.matches('\n').count(), 1);
        }

        fn minimal_report_tsv(rows: &[&str]) -> String {
            let mut tsv = concat!(
                "Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\t",
                "Scope\tType\tSubtype\n"
            )
            .to_string();
            for row in rows {
                tsv.push_str(row);
                tsv.push('\n');
            }
            tsv
        }

        #[test]
        fn nucleotide_report_deredundify_prefers_point_over_point_disrupt_plus_start() {
            let mut config = PipelineConfig::default();
            config.nucleotide = Some(PathBuf::from("dna.fa"));
            let tsv = minimal_report_tsv(&[
                "NA\tcontig\t10\t20\t+\tgyrA\tname\tcore\tAMR\tPOINT_DISRUPT",
                "NA\tcontig\t10\t18\t+\tgyrA\tname\tcore\tAMR\tPOINT",
            ]);

            let sorted = sort_tsv_output(&tsv, &config, true).unwrap();

            assert!(sorted.contains("\tPOINT\n"));
            assert!(!sorted.contains("\tPOINT_DISRUPT\n"));
        }

        #[test]
        fn nucleotide_report_deredundify_prefers_point_over_point_disrupt_minus_stop() {
            let mut config = PipelineConfig::default();
            config.nucleotide = Some(PathBuf::from("dna.fa"));
            let tsv = minimal_report_tsv(&[
                "NA\tcontig\t10\t30\t-\tpmrB\tname\tcore\tAMR\tPOINT_DISRUPT",
                "NA\tcontig\t12\t30\t-\tpmrB\tname\tcore\tAMR\tPOINT",
            ]);

            let sorted = sort_tsv_output(&tsv, &config, true).unwrap();

            assert!(sorted.contains("\tPOINT\n"));
            assert!(!sorted.contains("\tPOINT_DISRUPT\n"));
        }

        #[test]
        fn mutation_all_sort_does_not_deredundify_point_disrupt_rows() {
            let mut config = PipelineConfig::default();
            config.nucleotide = Some(PathBuf::from("dna.fa"));
            let tsv = minimal_report_tsv(&[
                "NA\tcontig\t10\t20\t+\tgyrA\tname\tcore\tAMR\tPOINT_DISRUPT",
                "NA\tcontig\t10\t18\t+\tgyrA\tname\tcore\tAMR\tPOINT",
            ]);

            let sorted = sort_tsv_output(&tsv, &config, false).unwrap();

            assert!(sorted.contains("\tPOINT\n"));
            assert!(sorted.contains("\tPOINT_DISRUPT\n"));
        }

        #[test]
        fn mutation_all_postprocessing_replaces_disruptions_before_sort_like_upstream() {
            let source = include_str!("../src/amrfinder.rs");
            let mutation_all_block = source
                .find("if let Some(path) = &config.mutation_all")
                .expect("translated mutation_all post-processing should be present");
            let block = &source[mutation_all_block..];
            let disruption = block
                .find("amr_tab_disruptions(")
                .expect("mutation_all post-processing should replace disruptions");
            let sort = block
                .find("sort_tsv_output(&mutation_all, config, false)")
                .expect("mutation_all post-processing should sort after replacement");

            assert!(disruption < sort);
        }

        #[test]
        fn report_fasta_outputs_use_raw_tmp_amr_order_like_upstream() {
            let source = include_str!("../src/amrfinder.rs");
            let report_fn = source
                .find("fn run_rust_amr_report(")
                .expect("run_rust_amr_report should remain present");
            let report_block = &source[report_fn..];
            let stxtyper = report_block
                .find("append_stxtyper_rows_if_available(config, tmp, db, &mut raw_result, nucleotide_flat)?")
                .expect("StxTyper rows should be appended to tmp/amr equivalent");
            let fasta_source = report_block
                .find("let fasta_source = raw_result.clone();")
                .expect("raw tmp/amr equivalent should be preserved for FASTA extraction");
            let disruption = report_block
                .find("amr_tab_disruptions(&mut amr_tab, db, nucleotide_flat, config.translation_table)?")
                .expect("main report disruption post-processing should remain present");

            assert!(stxtyper < fasta_source);
            assert!(fasta_source < disruption);
            assert!(source.contains("write_report_fasta_outputs(\n        &fasta_source,"));
        }

        #[test]
        fn final_report_sort_uses_text_table_numeric_tie_breakers_like_upstream() {
            let config = PipelineConfig::default();
            let tsv = concat!(
                "Protein id\tElement symbol\tElement name\tScope\tType\tSubtype\tTarget length\n",
                "prot\tgene\tname\tcore\tAMR\tPARTIAL\t100\n",
                "prot\tgene\tname\tcore\tAMR\tPARTIAL\t2\n",
            );

            let sorted = sort_tsv_output(tsv, &config, true).unwrap();
            let rows: Vec<&str> = sorted.lines().skip(1).collect();

            assert_eq!(rows[0].split('\t').next_back(), Some("2"));
            assert_eq!(rows[1].split('\t').next_back(), Some("100"));
        }

        #[test]
        fn stxtyper_rows_are_appended_without_headers_or_duplicates() {
            let mut output = concat!(
                "Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\tSubtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference\t% Identity to reference\tAlignment length\tClosest reference accession\tClosest reference name\tHMM accession\tHMM description\tHierarchy node\n",
                "NA\tcontig18\t279\t1235\t+\tstxA2\tShiga toxin Stx2 subunit A\tplus\tVIRULENCE\tVIRULENCE\tSTX2\tstxA2\tEXACTX\t319\t319\t100.00\t100.00\t319\tTJA36680.1\tShiga toxin Stx2 subunit A\tNA\tNA\tstxA2_acd\n",
            )
            .to_string();
            let stxtyper_output = concat!(
                "\n",
                "Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\tSubtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference\t% Identity to reference\tAlignment length\tClosest reference accession\tClosest reference name\tHMM accession\tHMM description\tHierarchy node\n",
                "NA\tcontig18\t279\t1516\t+\tstx2a_operon\tstx2a operon\tplus\tVIRULENCE\tSTX_TYPE\tSTX2\tSTX2A\tCOMPLETE\t1238\tNA\tNA\t100.00\t408\tTJA36680.1,AAM90978.1\tShiga toxin stx2a\tNA\tNA\tstxA2_acd::stxB2a\n",
                "NA\tcontig18\t279\t1516\t+\tstx2a_operon\tstx2a operon\tplus\tVIRULENCE\tSTX_TYPE\tSTX2\tSTX2A\tCOMPLETE\t1238\tNA\tNA\t100.00\t408\tTJA36680.1,AAM90978.1\tShiga toxin stx2a\tNA\tNA\tstxA2_acd::stxB2a\n",
            );

            append_report_rows(&mut output, stxtyper_output);
            append_report_rows(&mut output, stxtyper_output);

            assert_eq!(
                output.matches("\tstx2a_operon\t").count(),
                1,
                "stxtyper appends should be idempotent"
            );
            assert_eq!(
                output.matches("Protein id\tContig id").count(),
                1,
                "stxtyper headers should not be appended"
            );
        }

        #[test]
        #[cfg(unix)]
        fn stxtyper_runs_even_when_operon_row_already_exists_like_upstream() {
            use std::os::unix::fs::PermissionsExt;

            let _guard = ENV_LOCK.lock().unwrap();
            let tmp = tempfile::tempdir().unwrap();
            let db = tmp.path().join("db");
            let stxtyper_log = tmp.path().join("stxtyper-args");
            let stxtyper = test_stxtyper_path();
            fs::write(
                &stxtyper,
                format!(
                    concat!(
                        "#!/bin/sh\n",
                        "printf '%s\\n' \"$(pwd -P)\" \"$@\" > '{}'\n",
                        "case \" $* \" in *' --nucleotide '*|*' --output '*|*' --quiet '*) exit 9;; esac\n",
                        "out=\n",
                        "name_seen=0\n",
                        "while [ \"$#\" -gt 0 ]; do\n",
                        "  if [ \"$1\" = \"-o\" ]; then out=\"$2\"; shift; fi\n",
                        "  if [ \"$1\" = \"--name\" ]; then name_seen=1; shift; fi\n",
                        "  shift\n",
                        "done\n",
                        "if [ \"$name_seen\" != 1 ]; then exit 8; fi\n",
                        "printf '%s\\n' 'Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\tSubtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference\t% Identity to reference\tAlignment length\tClosest reference accession\tClosest reference name\tHMM accession\tHMM description\tHierarchy node' > \"$out\"\n",
                        "printf '%s\\n' 'NA\tcontig18\t1600\t1700\t+\tstx2a_context\tstx2a context\tplus\tVIRULENCE\tSTX_TYPE\tSTX2\tSTX2A\tCOMPLETE\t100\tNA\tNA\t100.00\t100\tCTX1\tStxTyper context\tNA\tNA\tstx_context' >> \"$out\"\n",
                    ),
                    stxtyper_log.display(),
                ),
            )
            .unwrap();
            let mut perms = fs::metadata(&stxtyper).unwrap().permissions();
            perms.set_mode(0o755);
            fs::set_permissions(&stxtyper, perms).unwrap();

            let nucleotide = tmp.path().join("query.fa");
            fs::write(&nucleotide, ">contig18\nACGT\n").unwrap();
            let mut output = concat!(
                "Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\tSubtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference\t% Identity to reference\tAlignment length\tClosest reference accession\tClosest reference name\tHMM accession\tHMM description\tHierarchy node\n",
                "NA\tcontig18\t279\t1516\t+\tstx2a_operon\tstx2a operon\tplus\tVIRULENCE\tSTX_TYPE\tSTX2\tSTX2A\tCOMPLETE\t1238\tNA\tNA\t100.00\t408\tTJA36680.1,AAM90978.1\tShiga toxin stx2a\tNA\tNA\tstxA2_acd::stxB2a\n",
            )
            .to_string();
            let config = PipelineConfig {
                nucleotide: Some(nucleotide),
                organism: "Escherichia".to_string(),
                plus: true,
                ..PipelineConfig::default()
            };

            append_stxtyper_rows_if_available(
                &config,
                tmp.path(),
                db.to_str().unwrap(),
                &mut output,
                None,
            )
            .unwrap();

            assert!(output.contains("\tstx2a_operon\t"));
            assert!(output.contains("\tstx2a_context\t"));
            let log = fs::read_to_string(stxtyper_log).unwrap();
            let logged: Vec<&str> = log.lines().collect();
            assert_eq!(
                logged.first().copied(),
                Some(
                    std::env::current_dir()
                        .unwrap()
                        .canonicalize()
                        .unwrap()
                        .to_str()
                        .unwrap()
                )
            );
            assert!(logged
                .windows(2)
                .any(|pair| pair == ["-n", config.nucleotide.as_ref().unwrap().to_str().unwrap()]));
            assert!(logged
                .windows(2)
                .any(|pair| pair == ["-o", tmp.path().join("stxtyper.tsv").to_str().unwrap()]));
            assert!(logged.windows(2).any(|pair| pair == ["--name", ""]));
            assert!(logged.contains(&"-q"));
        }

        #[test]
        fn stxtyper_unavailable_is_an_error() {
            let _guard = ENV_LOCK.lock().unwrap();
            let tmp = tempfile::tempdir().unwrap();
            let db = tmp.path().join("db");
            fs::create_dir(&db).unwrap();
            let _ = fs::remove_file(test_stxtyper_path());
            let nucleotide = tmp.path().join("query.fa");
            fs::write(&nucleotide, ">contig18\nACGT\n").unwrap();
            let mut output = concat!(
                "Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\tSubtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference\t% Identity to reference\tAlignment length\tClosest reference accession\tClosest reference name\tHMM accession\tHMM description\tHierarchy node\n",
                "NA\tcontig18\t279\t1235\t+\tstxA2\tShiga toxin Stx2 subunit A\tplus\tVIRULENCE\tVIRULENCE\tSTX2\tstxA2\tEXACTX\t319\t319\t100.00\t100.00\t319\tTJA36680.1\tShiga toxin Stx2 subunit A\tNA\tNA\tstxA2_acd\n",
                "NA\tcontig18\t1250\t1516\t+\tstxB2\tShiga toxin Stx2a subunit B\tplus\tVIRULENCE\tVIRULENCE\tSTX2\tstxB2a\tEXACTX\t89\t89\t100.00\t100.00\t89\tAAM90978.1\tShiga toxin Stx2a subunit B\tNA\tNA\tstxB2a\n"
            )
            .to_string();

            let original = output.clone();
            let config = PipelineConfig {
                nucleotide: Some(nucleotide),
                organism: "Escherichia".to_string(),
                plus: true,
                ..PipelineConfig::default()
            };

            let err = append_stxtyper_rows_if_available(
                &config,
                tmp.path(),
                db.to_str().unwrap(),
                &mut output,
                None,
            )
            .unwrap_err();

            assert!(err.to_string().contains("stxtyper is required"));
            assert_eq!(output, original);
            assert!(!output.contains("\tstx2a_operon\t"));
        }

        #[test]
        #[cfg(unix)]
        fn non_executable_stxtyper_is_an_error() {
            use std::os::unix::fs::PermissionsExt;

            let _guard = ENV_LOCK.lock().unwrap();
            let tmp = tempfile::tempdir().unwrap();
            let db = tmp.path().join("db");
            let stxtyper = test_stxtyper_path();
            fs::write(&stxtyper, "#!/bin/sh\nexit 0\n").unwrap();
            let mut perms = fs::metadata(&stxtyper).unwrap().permissions();
            perms.set_mode(0o644);
            fs::set_permissions(&stxtyper, perms).unwrap();
            let nucleotide = tmp.path().join("query.fa");
            fs::write(&nucleotide, ">contig18\nACGT\n").unwrap();
            let mut output = concat!(
                "Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\tSubtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference\t% Identity to reference\tAlignment length\tClosest reference accession\tClosest reference name\tHMM accession\tHMM description\tHierarchy node\n",
                "NA\tcontig18\t279\t1235\t+\tstxA2\tShiga toxin Stx2 subunit A\tplus\tVIRULENCE\tVIRULENCE\tSTX2\tstxA2\tEXACTX\t319\t319\t100.00\t100.00\t319\tTJA36680.1\tShiga toxin Stx2 subunit A\tNA\tNA\tstxA2_acd\n",
                "NA\tcontig18\t1250\t1516\t+\tstxB2\tShiga toxin Stx2a subunit B\tplus\tVIRULENCE\tVIRULENCE\tSTX2\tstxB2a\tEXACTX\t89\t89\t100.00\t100.00\t89\tAAM90978.1\tShiga toxin Stx2a subunit B\tNA\tNA\tstxB2a\n"
            )
            .to_string();

            let original = output.clone();
            let config = PipelineConfig {
                nucleotide: Some(nucleotide),
                organism: "Escherichia".to_string(),
                plus: true,
                ..PipelineConfig::default()
            };

            let err = append_stxtyper_rows_if_available(
                &config,
                tmp.path(),
                db.to_str().unwrap(),
                &mut output,
                None,
            )
            .unwrap_err();

            assert!(err.to_string().contains("stxtyper is required"));
            assert_eq!(output, original);
            assert!(!output.contains("\tstx2a_operon\t"));
        }

        #[test]
        #[cfg(unix)]
        fn database_local_stxtyper_is_ignored_like_upstream() {
            use std::os::unix::fs::PermissionsExt;

            let _guard = ENV_LOCK.lock().unwrap();
            let tmp = tempfile::tempdir().unwrap();
            let db = tmp.path().join("db");
            let db_stx_dir = db.join("stx");
            fs::create_dir_all(&db_stx_dir).unwrap();
            let db_stxtyper = db_stx_dir.join("stxtyper");
            fs::write(&db_stxtyper, "#!/bin/sh\nexit 0\n").unwrap();
            let mut perms = fs::metadata(&db_stxtyper).unwrap().permissions();
            perms.set_mode(0o755);
            fs::set_permissions(&db_stxtyper, perms).unwrap();
            let _ = fs::remove_file(test_stxtyper_path());

            assert!(find_stxtyper(&db).is_none());
        }
    }
}
