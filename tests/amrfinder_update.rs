mod curl_easy {
    pub use amrfinder::curl_easy::*;
}

mod common {
    pub use amrfinder::common::*;
}

mod amrfinder_update {
    include!("../src/amrfinder_update.rs");

    #[cfg(test)]
    mod tests {
        use super::*;
        use std::io::{Read, Write};
        use std::net::{TcpListener, TcpStream};
        use std::thread;

        #[test]
        fn parses_greatest_minor_from_directory_listing() {
            let listing = r#"
<a href="3.10/">3.10/</a>
<a href="3.9/">3.9/</a>
<a href="latest/">latest/</a>
<a href="3.11/">3.11/</a>
"#;

            assert_eq!(parse_latest_minor(listing).unwrap(), "3.11");
        }

        #[test]
        fn directory_parsers_return_canonical_version_strings_like_cpp_objects() {
            assert_eq!(
                parse_latest_minor(
                    r#"
<a href="03.011/">03.011/</a>
<a href="3.10/">3.10/</a>
"#
                )
                .unwrap(),
                "3.11"
            );

            assert_eq!(
                parse_latest_data_version(
                    r#"
<a href="2024-7-2.1/">2024-7-2.1/</a>
<a href="2024-07-01.9/">2024-07-01.9/</a>
"#
                )
                .unwrap(),
                "2024-07-02.1"
            );
        }

        #[test]
        fn directory_parsers_ignore_version_like_non_anchor_text() {
            let listing = r#"
            latest prose mentions 99.99 but it is not a directory
<a href="3.11/">3.11/</a>
"#;
            assert_eq!(parse_latest_minor(listing).unwrap(), "3.11");

            let listing = r#"
            latest prose mentions 2099-01-01.9 but it is not a directory
<a href="2024-07-22.1/">2024-07-22.1/</a>
"#;
            assert_eq!(parse_latest_data_version(listing).unwrap(), "2024-07-22.1");
        }

        #[test]
        fn directory_parsers_trim_listing_lines_before_left_anchor_check_like_cpp() {
            let listing = r#"
  <a href="3.12/">3.12/</a>
prefix <a href="3.13/">3.13</a>
<a href="3.10/">3.10/</a><a href="3.11/">3.11/</a>
<a href="3.9/">3.9/</a>
"#;

            assert_eq!(parse_latest_minor(listing).unwrap(), "3.12");

            let listing = r#"
  <a href="2024-07-22.2/">2024-07-22.2/</a>
prefix <a href="2099-01-01.1/">2099-01-01.1/</a>
<a href="2024-07-22.1/">2024-07-22.1/</a>
"#;

            assert_eq!(parse_latest_data_version(listing).unwrap(), "2024-07-22.2");
        }

        #[test]
        fn directory_parsers_ignore_left_anchor_without_directory_body_like_cpp() {
            let listing = r#"
<a href="3.99">3.99</a>
<a href="3.12/">3.12/</a>
"#;

            assert_eq!(parse_latest_minor(listing).unwrap(), "3.12");

            let listing = r#"
<a href="2099-01-01.1">2099-01-01.1</a>
<a href="2024-07-22.2/">2024-07-22.2/</a>
"#;

            assert_eq!(parse_latest_data_version(listing).unwrap(), "2024-07-22.2");
        }

        #[test]
        fn directory_parsers_ignore_malformed_left_anchor_before_body_like_cpp() {
            let listing = r#"
<a href="/<broken">3.99</a>
<a href="3.12/">3.12/</a>
"#;

            assert_eq!(parse_latest_minor(listing).unwrap(), "3.12");

            let listing = r#"
<a href="/<broken">2099-01-01.1</a>
<a href="2024-07-22.2/">2024-07-22.2/</a>
"#;

            assert_eq!(parse_latest_data_version(listing).unwrap(), "2024-07-22.2");
        }

        #[test]
        fn normalizes_first_trimmed_version_line_like_cpp_string_vector() {
            assert_eq!(
                normalize_version("  2024-07-22.1\ncomment\n").unwrap(),
                "2024-07-22.1"
            );
            assert_eq!(normalize_version("\n2024-07-22.1\n").unwrap(), "");
        }

        #[test]
        fn parses_greatest_data_version_from_directory_listing() {
            assert_eq!(
                parse_latest_data_version(
                    r#"
<a href="2024-07-22.1/">2024-07-22.1/</a>
<a href="2024-07-22.2/">2024-07-22.2/</a>
<a href="2023-12-01.1/">2023-12-01.1/</a>
"#
                )
                .unwrap(),
                "2024-07-22.2"
            );
        }

        #[test]
        fn builds_expected_db_file_url() {
            assert_eq!(
                db_file_url(
                    "https://example.test/db/",
                    "3.11",
                    "2024-07-22.1",
                    "AMR.LIB"
                ),
                "https://example.test/db/3.11/2024-07-22.1/AMR.LIB"
            );
        }

        #[test]
        fn fetch_amr_file_downloads_file_url_through_translated_curl() {
            let temp = tempfile::tempdir().unwrap();
            let source = temp.path().join("source.tar.gz");
            let dest_dir = temp.path().join("out");
            let dest = dest_dir.join("source.tar.gz");
            fs::write(&source, b"archive").unwrap();
            fs::create_dir(&dest_dir).unwrap();
            let mut curl = Curl::new().unwrap();

            fetch_amr_file(
                &mut curl,
                &format!("file://{}/", source.parent().unwrap().display()),
                &dir_name_path(&dest_dir),
                source.file_name().unwrap().to_str().unwrap(),
            )
            .unwrap();

            assert_eq!(fs::read(dest).unwrap(), b"archive");
        }

        #[test]
        fn fetch_amr_file_concatenates_local_dir_and_file_name_like_cpp() {
            let temp = tempfile::tempdir().unwrap();
            let source_dir = temp.path().join("source");
            let dest_dir = temp.path().join("out");
            fs::create_dir(&source_dir).unwrap();
            fs::create_dir(&dest_dir).unwrap();
            fs::write(source_dir.join("artifact"), b"archive").unwrap();
            let mut curl = Curl::new().unwrap();

            fetch_amr_file(
                &mut curl,
                &format!("file://{}/", source_dir.display()),
                &dir_name_path(&dest_dir),
                "/artifact",
            )
            .unwrap();

            assert_eq!(fs::read(dest_dir.join("artifact")).unwrap(), b"archive");
        }

        #[test]
        #[should_panic(expected = "urlDir must be a directory name")]
        fn fetch_amr_file_requires_directory_url_like_cpp_assert() {
            let temp = tempfile::tempdir().unwrap();
            let mut curl = Curl::new().unwrap();

            let _ = fetch_amr_file(&mut curl, "file:///tmp/amr-db", temp.path(), "AMR.LIB");
        }

        #[test]
        #[should_panic(expected = "localDir must be a directory name")]
        fn fetch_amr_file_requires_directory_local_dir_like_cpp_assert() {
            let temp = tempfile::tempdir().unwrap();
            let mut curl = Curl::new().unwrap();

            let _ = fetch_amr_file(&mut curl, "file:///tmp/", temp.path(), "AMR.LIB");
        }

        #[test]
        fn file_directory_urls_are_not_rewritten_to_index_html() {
            assert_eq!(
                latest_minor_index_url("file:///tmp/amr-db"),
                "file:///tmp/amr-db/"
            );
            assert_eq!(
                latest_data_version_index_url("file:///tmp/amr-db", "3.11"),
                "file:///tmp/amr-db/3.11/"
            );
        }

        #[test]
        fn dir_name_path_normalizes_like_cpp_dir_before_adding_slash() {
            assert_eq!(dir_name_path(Path::new("db/../db")), PathBuf::from("db/"));
            assert_eq!(dir_name_path(Path::new("missing/..")), PathBuf::from("./"));
            assert_eq!(
                dir_name_path(Path::new("/tmp//amr/./old/../db")),
                PathBuf::from("/tmp/amr/db/")
            );
            assert_eq!(dir_name_path(Path::new("/")), PathBuf::from("/"));
        }

        #[test]
        fn rust_only_compatibility_api_is_not_exposed() {
            let source = include_str!("../src/amrfinder_update.rs");

            assert!(!source.contains("pub fn download_file"));
            assert!(!source.contains("pub fn get_latest_version"));
            assert!(!source.contains("pub fn check_update"));
            assert!(!source.contains("pub fn update("));
            assert!(!source.contains("pub fn update_with_index_options"));
            assert!(!source.contains("pub fn shell_body_with_index_runner"));
            assert!(!source.contains("fn shell_body_with_index_runner"));
            assert!(!source.contains("pub fn read_url"));
            assert!(!source.contains("pub fn get_latest_minor()"));
            assert!(!source.contains("pub fn get_latest_data_version()"));
            assert!(!source.contains("ureq"));
            assert!(!source.contains("use std::process::Command"));
            assert!(!source.contains("Command::new"));
        }

        #[test]
        fn updater_lookup_functions_keep_original_curl_argument_shape() {
            let source = include_str!("../src/amrfinder_update.rs");

            assert!(source.contains("pub fn get_latest_minor(curl: &mut Curl)"));
            assert!(source.contains("pub fn get_latest_data_version(curl: &mut Curl, minor: &str)"));
        }

        #[test]
        fn rust_only_update_option_constructors_are_not_present() {
            let source = include_str!("../src/amrfinder_update.rs");

            assert!(!source.contains("impl UpdateOptions"));
            assert!(!source.contains("impl Default for AmrFinderIndexOptions"));
        }

        #[test]
        fn current_software_minor_uses_upstream_svn_revision_source() {
            assert_eq!(
                current_software_minor(),
                minor_from_software_version(include_str!("../amr/version.txt").trim())
            );
            assert_ne!(
                current_software_minor(),
                minor_from_software_version(env!("CARGO_PKG_VERSION"))
            );
        }

        #[test]
        fn shell_body_skips_when_current_without_fetching() {
            let source = tempfile::tempdir().unwrap();
            fs::write(
                source.path().join("index.html"),
                r#"<a href="3.11/">3.11/</a>"#,
            )
            .unwrap();
            fs::create_dir(source.path().join("3.11")).unwrap();
            fs::write(
                source.path().join("3.11").join("index.html"),
                r#"<a href="2024-07-22.1/">2024-07-22.1/</a>"#,
            )
            .unwrap();
            let source_version = source.path().join("3.11").join("2024-07-22.1");
            fs::create_dir(&source_version).unwrap();
            fs::write(source_version.join(VERSION_FILE), "2024-07-22.1\n").unwrap();
            let server = HttpFixture::serve(source.path());

            let database = tempfile::tempdir().unwrap();
            fs::create_dir(database.path().join("2024-07-22.1")).unwrap();
            fs::write(
                database.path().join("2024-07-22.1").join(VERSION_FILE),
                "2024-07-22.1\n",
            )
            .unwrap();
            create_latest_link(database.path(), "2024-07-22.1").unwrap();

            let outcome = shell_body(UpdateOptions {
                database_dir: database.path().to_path_buf(),
                force: false,
                source_base: server.url(),
                current_minor: "3.11".to_string(),
                amrfinder_index: None,
            })
            .unwrap();

            assert_eq!(outcome.status, UpdateStatus::AlreadyCurrent);
            assert!(outcome.fetched_files.is_empty());
            assert!(outcome.warnings.is_empty());
            assert_latest_points_to(database.path(), "2024-07-22.1");
        }

        #[test]
        fn shell_body_compares_first_version_file_line_like_cpp() {
            let source = tempfile::tempdir().unwrap();
            fs::write(
                source.path().join("index.html"),
                r#"<a href="3.11/">3.11/</a>"#,
            )
            .unwrap();
            fs::create_dir(source.path().join("3.11")).unwrap();
            fs::write(
                source.path().join("3.11").join("index.html"),
                r#"<a href="2024-07-22.1/">2024-07-22.1/</a>"#,
            )
            .unwrap();
            let source_version = source.path().join("3.11").join("2024-07-22.1");
            fs::create_dir(&source_version).unwrap();
            fs::write(source_version.join(VERSION_FILE), "\nremote-build\n").unwrap();
            let server = HttpFixture::serve(source.path());

            let database = tempfile::tempdir().unwrap();
            fs::create_dir(database.path().join("2024-07-22.1")).unwrap();
            fs::write(
                database.path().join("2024-07-22.1").join(VERSION_FILE),
                "\nlocal-build\n",
            )
            .unwrap();

            let outcome = shell_body(UpdateOptions {
                database_dir: database.path().to_path_buf(),
                force: false,
                source_base: server.url(),
                current_minor: "3.11".to_string(),
                amrfinder_index: None,
            })
            .unwrap();

            assert_eq!(outcome.status, UpdateStatus::AlreadyCurrent);
            assert!(outcome.fetched_files.is_empty());
            assert_latest_points_to(database.path(), "2024-07-22.1");
        }

        #[test]
        fn shell_body_fetches_when_existing_target_version_file_differs_from_remote() {
            let source = tempfile::tempdir().unwrap();
            fs::write(
                source.path().join("index.html"),
                r#"<a href="3.11/">3.11/</a>"#,
            )
            .unwrap();
            fs::create_dir(source.path().join("3.11")).unwrap();
            fs::write(
                source.path().join("3.11").join("index.html"),
                r#"<a href="2024-07-22.1/">2024-07-22.1/</a>"#,
            )
            .unwrap();
            let source_version = source.path().join("3.11").join("2024-07-22.1");
            fs::create_dir(&source_version).unwrap();
            for file_name in CORE_DB_FILES {
                fs::write(
                    source_version.join(file_name),
                    format!("remote {file_name}\n"),
                )
                .unwrap();
            }
            fs::write(source_version.join(VERSION_FILE), "remote-build\n").unwrap();
            fs::write(
                source_version.join("taxgroup.tsv"),
                "#taxgroup\tgpipe\tmutations\n",
            )
            .unwrap();
            fs::write(source_version.join(CHANGES_FILE), b"changes\n").unwrap();
            let server = HttpFixture::serve(source.path());

            let database = tempfile::tempdir().unwrap();
            fs::create_dir(database.path().join("2024-07-22.1")).unwrap();
            fs::write(
                database.path().join("2024-07-22.1").join(VERSION_FILE),
                "local-build\n",
            )
            .unwrap();

            let outcome = shell_body(UpdateOptions {
                database_dir: database.path().to_path_buf(),
                force: false,
                source_base: server.url(),
                current_minor: "3.11".to_string(),
                amrfinder_index: None,
            })
            .unwrap();

            assert_eq!(outcome.status, UpdateStatus::Updated);
            assert_eq!(
                fs::read_to_string(database.path().join("2024-07-22.1").join(VERSION_FILE))
                    .unwrap(),
                "remote-build\n"
            );
        }

        #[test]
        fn shell_body_errors_when_existing_target_version_file_is_missing_like_cpp() {
            let source = tempfile::tempdir().unwrap();
            fs::write(
                source.path().join("index.html"),
                r#"<a href="3.11/">3.11/</a>"#,
            )
            .unwrap();
            fs::create_dir(source.path().join("3.11")).unwrap();
            fs::write(
                source.path().join("3.11").join("index.html"),
                r#"<a href="2024-07-22.1/">2024-07-22.1/</a>"#,
            )
            .unwrap();
            let source_version = source.path().join("3.11").join("2024-07-22.1");
            fs::create_dir(&source_version).unwrap();
            fs::write(source_version.join(VERSION_FILE), "remote-build\n").unwrap();
            let server = HttpFixture::serve(source.path());

            let database = tempfile::tempdir().unwrap();
            let version_dir = database.path().join("2024-07-22.1");
            fs::create_dir(&version_dir).unwrap();

            let err = shell_body(UpdateOptions {
                database_dir: database.path().to_path_buf(),
                force: false,
                source_base: server.url(),
                current_minor: "3.11".to_string(),
                amrfinder_index: None,
            })
            .unwrap_err();

            assert!(err.to_string().contains("Reading file"));
            assert!(format!("{err:#}").contains("version.txt"));
            assert!(!version_dir.join("AMR.LIB").exists());
        }

        #[test]
        fn planner_does_not_fallback_when_current_minor_listing_cannot_be_read() {
            let source = tempfile::tempdir().unwrap();
            fs::write(
                source.path().join("index.html"),
                "<a href=\"3.11/\">3.11/</a>\n<a href=\"3.12/\">3.12/</a>\n",
            )
            .unwrap();
            fs::create_dir(source.path().join("3.12")).unwrap();
            fs::write(
                source.path().join("3.12").join("index.html"),
                r#"<a href="2024-08-01.1/">2024-08-01.1/</a>"#,
            )
            .unwrap();
            let server = HttpFixture::serve(source.path());

            let mut curl = Curl::new().unwrap();
            let err = plan_update(&mut curl, &server.url(), "3.11").unwrap_err();

            assert!(err.to_string().contains("CURL: Cannot read"));
        }

        #[test]
        fn planner_reports_cpp_error_when_published_minor_listing_is_empty() {
            let source = tempfile::tempdir().unwrap();
            fs::write(source.path().join("index.html"), "").unwrap();
            let server = HttpFixture::serve(source.path());

            let mut curl = Curl::new().unwrap();
            let err = plan_update(&mut curl, &server.url(), "3.11").unwrap_err();

            assert_eq!(
                err.to_string(),
                "Cannot get the software minor version of the latest published database version"
            );
        }

        #[test]
        fn planner_reports_cpp_error_when_published_data_listing_is_empty() {
            let source = tempfile::tempdir().unwrap();
            fs::write(
                source.path().join("index.html"),
                r#"<a href="3.12/">3.12/</a>"#,
            )
            .unwrap();
            fs::create_dir(source.path().join("3.12")).unwrap();
            fs::write(source.path().join("3.12").join("index.html"), "").unwrap();
            let server = HttpFixture::serve(source.path());

            let mut curl = Curl::new().unwrap();
            let err = plan_update(&mut curl, &server.url(), "3.11").unwrap_err();

            assert_eq!(
                err.to_string(),
                "Cannot get the latest published database version for the software minor version 3.12"
            );
        }

        #[test]
        fn shell_body_fetches_files_into_version_dir_and_updates_latest() {
            let source = tempfile::tempdir().unwrap();
            fs::write(
                source.path().join("index.html"),
                r#"<a href="3.11/">3.11/</a>"#,
            )
            .unwrap();
            fs::create_dir(source.path().join("3.11")).unwrap();
            fs::write(
                source.path().join("3.11").join("index.html"),
                r#"<a href="2024-07-22.1/">2024-07-22.1/</a>"#,
            )
            .unwrap();
            let source_version = source.path().join("3.11").join("2024-07-22.1");
            fs::create_dir(&source_version).unwrap();
            for file_name in CORE_DB_FILES {
                fs::write(source_version.join(file_name), format!("{file_name}\n")).unwrap();
            }
            fs::write(source_version.join(VERSION_FILE), "2024-07-22.1\n").unwrap();
            fs::write(
                source_version.join("taxgroup.tsv"),
                "#taxgroup\tgpipe\tmutations\nEscherichia\t-\t2\nSalmonella\t-\t0\n",
            )
            .unwrap();
            fs::write(source_version.join("AMR_DNA-Escherichia.fa"), b">dna\n").unwrap();
            fs::write(source_version.join("AMR_DNA-Escherichia.tsv"), b"gene\n").unwrap();
            fs::write(source_version.join(CHANGES_FILE), b"changes\n").unwrap();
            let server = HttpFixture::serve(source.path());

            let database = tempfile::tempdir().unwrap();

            let outcome = shell_body(UpdateOptions {
                database_dir: database.path().to_path_buf(),
                force: false,
                source_base: server.url(),
                current_minor: "3.11".to_string(),
                amrfinder_index: None,
            })
            .unwrap();

            assert_eq!(outcome.status, UpdateStatus::Updated);
            assert_eq!(outcome.latest_minor, "3.11");
            assert_eq!(outcome.latest_data_version, "2024-07-22.1");
            assert!(outcome
                .fetched_files
                .contains(&database.path().join("2024-07-22.1").join("AMR.LIB")));
            assert_eq!(
                fs::read_to_string(database.path().join("2024-07-22.1").join(VERSION_FILE))
                    .unwrap(),
                "2024-07-22.1\n"
            );
            assert_eq!(
                fs::read(
                    database
                        .path()
                        .join("2024-07-22.1")
                        .join("AMR_DNA-Escherichia.fa")
                )
                .unwrap(),
                b">dna\n"
            );
            assert!(!database
                .path()
                .join("2024-07-22.1")
                .join("AMR_DNA-Salmonella.fa")
                .exists());
            assert_latest_points_to(database.path(), "2024-07-22.1");
        }

        #[test]
        fn planner_prefers_current_minor_and_warns_when_newer_minor_exists() {
            let source = tempfile::tempdir().unwrap();
            fs::write(
                source.path().join("index.html"),
                "<a href=\"3.11/\">3.11/</a>\n<a href=\"3.12/\">3.12/</a>\n",
            )
            .unwrap();
            fs::create_dir(source.path().join("3.11")).unwrap();
            fs::write(
                source.path().join("3.11").join("index.html"),
                r#"<a href="2024-07-22.1/">2024-07-22.1/</a>"#,
            )
            .unwrap();
            fs::create_dir(source.path().join("3.12")).unwrap();
            fs::write(
                source.path().join("3.12").join("index.html"),
                r#"<a href="2024-08-01.1/">2024-08-01.1/</a>"#,
            )
            .unwrap();
            let server = HttpFixture::serve(source.path());

            let mut curl = Curl::new().unwrap();
            let planned = plan_update(&mut curl, &server.url(), "3.11").unwrap();

            assert_eq!(planned.minor, "3.11");
            assert_eq!(planned.data_version, "2024-07-22.1");
            assert_eq!(
                planned.warnings,
                [concat!(
                    "A newer version of the database exists (2024-08-01.1), but it requires ",
                    "a newer version of the software (3.12) to install.\n",
                    "See https://github.com/ncbi/amr/wiki/Upgrading for more information."
                )]
            );
        }

        #[test]
        #[should_panic(expected = "current data version must be older than published data version")]
        fn planner_asserts_when_current_minor_database_is_newer_than_published_like_cpp() {
            let source = tempfile::tempdir().unwrap();
            fs::write(
                source.path().join("index.html"),
                "<a href=\"3.11/\">3.11/</a>\n<a href=\"3.12/\">3.12/</a>\n",
            )
            .unwrap();
            fs::create_dir(source.path().join("3.11")).unwrap();
            fs::write(
                source.path().join("3.11").join("index.html"),
                r#"<a href="2024-08-01.1/">2024-08-01.1/</a>"#,
            )
            .unwrap();
            fs::create_dir(source.path().join("3.12")).unwrap();
            fs::write(
                source.path().join("3.12").join("index.html"),
                r#"<a href="2024-07-22.1/">2024-07-22.1/</a>"#,
            )
            .unwrap();
            let server = HttpFixture::serve(source.path());

            let mut curl = Curl::new().unwrap();
            let _ = plan_update(&mut curl, &server.url(), "3.11");
        }

        #[test]
        fn planner_falls_back_to_published_minor_when_current_minor_has_no_database() {
            let source = tempfile::tempdir().unwrap();
            fs::write(
                source.path().join("index.html"),
                "<a href=\"3.11/\">3.11/</a>\n<a href=\"3.12/\">3.12/</a>\n",
            )
            .unwrap();
            fs::create_dir(source.path().join("3.11")).unwrap();
            fs::write(source.path().join("3.11").join("index.html"), "").unwrap();
            fs::create_dir(source.path().join("3.12")).unwrap();
            fs::write(
                source.path().join("3.12").join("index.html"),
                r#"<a href="2024-08-01.1/">2024-08-01.1/</a>"#,
            )
            .unwrap();
            let server = HttpFixture::serve(source.path());

            let mut curl = Curl::new().unwrap();
            let planned = plan_update(&mut curl, &server.url(), "3.11").unwrap();

            assert_eq!(planned.minor, "3.12");
            assert_eq!(planned.data_version, "2024-08-01.1");
            assert_eq!(
                planned.warnings,
                [concat!(
                    "Cannot get the latest published database version for the current software minor version 3.11.\n",
                    "The latest published database version 2024-08-01.1 for the latest published software minor version 3.12 ",
                    "will be used instead"
                )]
            );
        }

        #[test]
        fn shell_body_runs_amrfinder_index_command_after_download() {
            let source = tempfile::tempdir().unwrap();
            fs::write(
                source.path().join("index.html"),
                r#"<a href="3.11/">3.11/</a>"#,
            )
            .unwrap();
            fs::create_dir(source.path().join("3.11")).unwrap();
            fs::write(
                source.path().join("3.11").join("index.html"),
                r#"<a href="2024-07-22.1/">2024-07-22.1/</a>"#,
            )
            .unwrap();
            let source_version = source.path().join("3.11").join("2024-07-22.1");
            fs::create_dir(&source_version).unwrap();
            for file_name in CORE_DB_FILES {
                fs::write(source_version.join(file_name), format!("{file_name}\n")).unwrap();
            }
            fs::write(source_version.join(VERSION_FILE), "2024-07-22.1\n").unwrap();
            fs::write(
                source_version.join("taxgroup.tsv"),
                "#taxgroup\tgpipe\tmutations\n",
            )
            .unwrap();
            fs::write(source_version.join(CHANGES_FILE), b"changes\n").unwrap();
            let server = HttpFixture::serve(source.path());

            let database = tempfile::tempdir().unwrap();
            let script_dir = tempfile::tempdir().unwrap();
            let script = script_dir.path().join("custom_amrfinder_index");
            let args_file = script_dir.path().join("args.txt");
            fs::write(
                &script,
                format!(
                    "#!/bin/sh\nfor arg in \"$@\"; do printf '%s\\n' \"$arg\"; done > {}\n",
                    crate::common::shell_quote(&args_file.to_string_lossy())
                ),
            )
            .unwrap();
            #[cfg(unix)]
            {
                use std::os::unix::fs::PermissionsExt;
                let mut permissions = fs::metadata(&script).unwrap().permissions();
                permissions.set_mode(0o755);
                fs::set_permissions(&script, permissions).unwrap();
            }

            let outcome = shell_body(UpdateOptions {
                database_dir: database.path().to_path_buf(),
                force: false,
                source_base: server.url(),
                current_minor: "3.11".to_string(),
                amrfinder_index: Some(AmrFinderIndexOptions {
                    program: script,
                    blast_bin: Some(PathBuf::from("/blast")),
                    hmmer_bin: Some(PathBuf::from("/hmmer")),
                    quiet: true,
                    debug: true,
                }),
            })
            .unwrap();

            assert_eq!(outcome.status, UpdateStatus::Updated);
            let args = fs::read_to_string(args_file).unwrap();
            assert_eq!(
                args.lines().collect::<Vec<_>>(),
                vec![
                    format!("{}/", database.path().join("2024-07-22.1").display()),
                    "-blast_bin".to_string(),
                    "/blast/".to_string(),
                    "-hmmer_bin".to_string(),
                    "/hmmer/".to_string(),
                    "-q".to_string(),
                    "--debug".to_string(),
                ]
            );
        }

        #[test]
        fn this_application_shell_body_uses_amrfinder_index_from_exec_dir() {
            let source = tempfile::tempdir().unwrap();
            fs::write(
                source.path().join("index.html"),
                r#"<a href="3.11/">3.11/</a>"#,
            )
            .unwrap();
            fs::create_dir(source.path().join("3.11")).unwrap();
            fs::write(
                source.path().join("3.11").join("index.html"),
                r#"<a href="2024-07-22.1/">2024-07-22.1/</a>"#,
            )
            .unwrap();
            let source_version = source.path().join("3.11").join("2024-07-22.1");
            fs::create_dir(&source_version).unwrap();
            for file_name in CORE_DB_FILES {
                fs::write(source_version.join(file_name), format!("{file_name}\n")).unwrap();
            }
            fs::write(source_version.join(VERSION_FILE), "2024-07-22.1\n").unwrap();
            fs::write(
                source_version.join("taxgroup.tsv"),
                "#taxgroup\tgpipe\tmutations\n",
            )
            .unwrap();
            fs::write(source_version.join(CHANGES_FILE), b"changes\n").unwrap();
            let server = HttpFixture::serve(source.path());

            let database = tempfile::tempdir().unwrap();
            let exe_dir = tempfile::tempdir().unwrap();
            let amrfinder_update = exe_dir.path().join("amrfinder_update");
            fs::write(&amrfinder_update, b"").unwrap();
            let amrfinder_index = exe_dir.path().join("amrfinder_index");
            let args_file = exe_dir.path().join("args.txt");
            fs::write(
                &amrfinder_index,
                format!(
                    "#!/bin/sh\nfor arg in \"$@\"; do printf '%s\\n' \"$arg\"; done > {}\n",
                    crate::common::shell_quote(&args_file.to_string_lossy())
                ),
            )
            .unwrap();
            #[cfg(unix)]
            {
                use std::os::unix::fs::PermissionsExt;
                let mut permissions = fs::metadata(&amrfinder_index).unwrap().permissions();
                permissions.set_mode(0o755);
                fs::set_permissions(&amrfinder_index, permissions).unwrap();
            }

            let mut app = ThisApplication::new();
            app.cur_minor = "3.11".to_string();
            app.source_base = server.url();
            app.amrfinder_index = amrfinder_update
                .canonicalize()
                .unwrap()
                .parent()
                .unwrap()
                .join("amrfinder_index");
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from(&amrfinder_update),
                    std::ffi::OsString::from("-d"),
                    std::ffi::OsString::from(database.path()),
                    std::ffi::OsString::from("--blast_bin"),
                    std::ffi::OsString::from("/blast"),
                    std::ffi::OsString::from("--hmmer_bin"),
                    std::ffi::OsString::from("/hmmer"),
                    std::ffi::OsString::from("-q"),
                    std::ffi::OsString::from("-debug"),
                ])
                .unwrap();

            let outcome = app.shell_body(&run).unwrap();

            assert_eq!(outcome.status, UpdateStatus::Updated);
            assert_eq!(
                fs::read_to_string(args_file)
                    .unwrap()
                    .lines()
                    .collect::<Vec<_>>(),
                vec![
                    format!("{}/", database.path().join("2024-07-22.1").display()),
                    "-blast_bin".to_string(),
                    "/blast/".to_string(),
                    "-hmmer_bin".to_string(),
                    "/hmmer/".to_string(),
                    "-q".to_string(),
                    "--debug".to_string(),
                ]
            );
        }

        #[test]
        fn this_application_parses_original_update_keys_and_database_default() {
            let temp = tempfile::tempdir().unwrap();
            let exe_dir = temp.path().join("bin");
            fs::create_dir(&exe_dir).unwrap();
            let exe = exe_dir.join("amrfinder_update");
            fs::write(&exe, b"").unwrap();

            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from(&exe),
                    std::ffi::OsString::from("-d"),
                    std::ffi::OsString::from("custom-db"),
                    std::ffi::OsString::from("--blast_bin"),
                    std::ffi::OsString::from("/blast"),
                    std::ffi::OsString::from("--hmmer_bin"),
                    std::ffi::OsString::from("/hmmer"),
                    std::ffi::OsString::from("--force_update"),
                    std::ffi::OsString::from("-q"),
                    std::ffi::OsString::from("--qc"),
                ])
                .unwrap();

            assert!(run.positional_args.is_empty());
            assert_eq!(run.key_args["database"], "custom-db");
            assert_eq!(run.key_args["blast_bin"], "/blast");
            assert_eq!(run.key_args["hmmer_bin"], "/hmmer");
            assert_eq!(run.key_args["force_update"], "true");
            assert_eq!(run.key_args["quiet"], "true");
            assert_eq!(run.key_args["qc"], "true");

            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![
                    std::ffi::OsString::from(&exe),
                    std::ffi::OsString::from("-debug"),
                ])
                .unwrap();
            assert_eq!(run.key_args["debug"], "true");

            let mut app = ThisApplication::new();
            let run = app
                .application
                .run(vec![std::ffi::OsString::from(&exe)])
                .unwrap();
            assert_eq!(
                run.key_args["database"],
                format!("{}/data", exe_dir.display())
            );
        }

        #[test]
        fn standalone_main_resolves_bare_program_name_through_path_like_cpp() {
            let source = include_str!("../src/amrfinder_update.rs");

            assert!(source.contains("which(&program.to_string_lossy())"));
            assert!(source.contains("found.canonicalize().unwrap_or(found)"));
            assert!(source.contains("argv[0] = OsString::from(&program_path);"));
            assert!(source.contains("app.amrfinder_index = exec_dir.join(\"amrfinder_index\");"));
        }

        #[test]
        fn this_application_description_matches_upstream_source_and_requirement() {
            let app = ThisApplication::new();
            assert!(app.application.description.contains(
                "Update the database for AMRFinder from https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/"
            ));
            assert!(app.application.description.contains(
                "Requirement: the database directory contains subdirectories named by database versions."
            ));
        }

        #[test]
        fn build_amrfinder_index_command_includes_tool_dirs_and_redirects_stdout() {
            let invocation = AmrFinderIndexInvocation {
                version_dir: PathBuf::from("/db/2024-07-22.1"),
                program: PathBuf::from("/tool dir/custom_amrfinder_index"),
                blast_bin: Some(PathBuf::from("/opt/blast/bin")),
                hmmer_bin: Some(PathBuf::from("/opt/hmmer/bin")),
                quiet: true,
                debug: true,
            };

            assert_eq!(
                build_amrfinder_index_command(
                    &invocation,
                    Path::new("/tmp/amrfinder update/amrfinder_index.err")
                ),
                "'/tool dir/custom_amrfinder_index' '/db/2024-07-22.1/' -blast_bin '/opt/blast/bin/' -hmmer_bin '/opt/hmmer/bin/' -q --debug > '/tmp/amrfinder update/amrfinder_index.err'"
            );
        }

        #[test]
        fn build_amrfinder_index_command_omits_absent_tool_dirs() {
            let invocation = AmrFinderIndexInvocation {
                version_dir: PathBuf::from("/db/2024-07-22.1"),
                program: PathBuf::from("amrfinder_index"),
                blast_bin: None,
                hmmer_bin: None,
                quiet: false,
                debug: false,
            };

            assert_eq!(
                build_amrfinder_index_command(&invocation, Path::new("/tmp/amrfinder_index.err")),
                "'amrfinder_index' '/db/2024-07-22.1/' > '/tmp/amrfinder_index.err'"
            );
        }

        #[test]
        fn dna_point_mutation_taxgroups_use_third_taxgroup_column() {
            let temp = tempfile::tempdir().unwrap();
            let taxgroup = temp.path().join("taxgroup.tsv");
            fs::write(
                &taxgroup,
                "# comment\nEscherichia\t-\t2\nSalmonella\t-\t0\nKlebsiella\t-\t1\n",
            )
            .unwrap();

            assert_eq!(
                dna_point_mutation_taxgroups(&taxgroup).unwrap(),
                vec!["Escherichia".to_string(), "Klebsiella".to_string()]
            );
        }

        #[test]
        fn dna_point_mutation_taxgroups_reject_negative_counts() {
            let temp = tempfile::tempdir().unwrap();
            let taxgroup = temp.path().join("taxgroup.tsv");
            fs::write(&taxgroup, "Escherichia\t-\t-1\n").unwrap();

            assert_eq!(
                dna_point_mutation_taxgroups(&taxgroup)
                    .unwrap_err()
                    .to_string(),
                format!("Bad {}", taxgroup.display())
            );
        }

        #[test]
        fn dna_point_mutation_taxgroups_parse_integer_prefix_like_cpp_streams() {
            let temp = tempfile::tempdir().unwrap();
            let taxgroup = temp.path().join("taxgroup.tsv");
            fs::write(&taxgroup, "Escherichia\t-\t2extra\nSalmonella\t-\t0zero\n").unwrap();

            assert_eq!(
                dna_point_mutation_taxgroups(&taxgroup).unwrap(),
                vec!["Escherichia".to_string()]
            );
        }

        #[test]
        fn dna_point_mutation_taxgroups_reject_blank_lines_like_cpp() {
            let temp = tempfile::tempdir().unwrap();
            let taxgroup = temp.path().join("taxgroup.tsv");
            fs::write(&taxgroup, "# comment\n\nEscherichia\t-\t2\n").unwrap();

            assert_eq!(
                dna_point_mutation_taxgroups(&taxgroup)
                    .unwrap_err()
                    .to_string(),
                format!("Bad {}", taxgroup.display())
            );
        }

        #[test]
        fn dna_point_mutation_taxgroups_only_skip_hash_in_first_column() {
            let temp = tempfile::tempdir().unwrap();
            let taxgroup = temp.path().join("taxgroup.tsv");
            fs::write(&taxgroup, "  # indented comment\nEscherichia\t-\t2\n").unwrap();

            assert_eq!(
                dna_point_mutation_taxgroups(&taxgroup)
                    .unwrap_err()
                    .to_string(),
                format!("Bad {}", taxgroup.display())
            );
        }

        #[test]
        fn create_latest_link_replaces_existing_empty_directory_like_cpp_remove() {
            let temp = tempfile::tempdir().unwrap();
            fs::create_dir(temp.path().join("2024-07-22.1")).unwrap();
            fs::create_dir(temp.path().join(LATEST_DIR)).unwrap();

            create_latest_link(temp.path(), "2024-07-22.1").unwrap();

            assert_latest_points_to(temp.path(), "2024-07-22.1");
        }

        #[test]
        fn create_latest_link_does_not_replace_existing_nonempty_directory_like_cpp_remove() {
            let temp = tempfile::tempdir().unwrap();
            fs::create_dir(temp.path().join("2024-07-22.1")).unwrap();
            fs::create_dir(temp.path().join(LATEST_DIR)).unwrap();
            fs::write(temp.path().join(LATEST_DIR).join("keep"), "").unwrap();

            let err = create_latest_link(temp.path(), "2024-07-22.1").unwrap_err();

            assert!(err.to_string().contains("latest"));
            assert!(temp.path().join(LATEST_DIR).is_dir());
        }

        #[test]
        fn create_latest_link_rejects_missing_relative_target_like_cpp_set_symlink() {
            let temp = tempfile::tempdir().unwrap();

            let err = create_latest_link(temp.path(), "2024-07-22.1").unwrap_err();

            assert!(err.to_string().contains("Cannot make a relative symlink"));
            assert!(err.to_string().contains("does not exist"));
            assert!(!temp.path().join(LATEST_DIR).exists());
        }

        fn assert_latest_points_to(database_dir: &Path, data_version: &str) {
            let latest = database_dir.join(LATEST_DIR);
            #[cfg(unix)]
            assert_eq!(fs::read_link(latest).unwrap(), PathBuf::from(data_version));

            #[cfg(not(unix))]
            assert_eq!(fs::read_to_string(latest).unwrap(), data_version);
        }

        struct HttpFixture {
            base_url: String,
        }

        impl HttpFixture {
            fn serve(root: &Path) -> Self {
                let listener = TcpListener::bind("127.0.0.1:0").unwrap();
                let base_url = format!("http://{}", listener.local_addr().unwrap());
                let root = root.to_path_buf();
                thread::spawn(move || {
                    for stream in listener.incoming().flatten() {
                        handle_http_request(stream, &root);
                    }
                });
                Self { base_url }
            }

            fn url(&self) -> String {
                self.base_url.clone()
            }
        }

        fn handle_http_request(mut stream: TcpStream, root: &Path) {
            let mut request = [0_u8; 2048];
            let Ok(n_read) = stream.read(&mut request) else {
                return;
            };
            let request = String::from_utf8_lossy(&request[..n_read]);
            let Some(target) = request
                .lines()
                .next()
                .and_then(|line| line.split_whitespace().nth(1))
            else {
                return;
            };

            let rel = target.trim_start_matches('/');
            let path = if rel.is_empty() || rel.ends_with('/') {
                root.join(rel).join("index.html")
            } else {
                root.join(rel)
            };
            let Ok(body) = fs::read(path) else {
                return;
            };

            let header = format!(
                "HTTP/1.0 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                body.len()
            );
            let _ = stream.write_all(header.as_bytes());
            let _ = stream.write_all(&body);
        }
    }
}
