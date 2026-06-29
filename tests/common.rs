use std::ffi::{OsStr, OsString};
use std::path::PathBuf;
use std::process::Command;

use amrfinder::common::{
    exec, get_command_line, get_threads_max, make_temp_dir, profile_enabled, progress_enabled,
    qc_on, quiet_enabled, sigpipe_enabled, Application, DataVersion, Key, ShellApplication,
    SoftwareVersion, DEFAULT_VERSION,
};

fn os_args(args: &[&str]) -> Vec<OsString> {
    args.iter().map(OsString::from).collect()
}

fn app(positionals: usize, key_args: Vec<Key>) -> Application {
    Application {
        description: "common test",
        usage: "common test",
        positionals,
        needs_arg: false,
        threads_used: false,
        key_args,
    }
}

#[test]
fn software_version_parses_orders_and_saves_text_like_upstream() {
    let mut out = Vec::new();
    let version = SoftwareVersion::from_text(" 4.2.6 trailing", false).unwrap();
    version.save_text(&mut out).unwrap();

    assert_eq!(version, SoftwareVersion::new(4, 2, 6));
    assert_eq!(String::from_utf8(out).unwrap(), "4.2.6");
    assert_eq!(version.get_minor(), "4.2");
    assert!(SoftwareVersion::new(4, 2, 5) < version);
    assert_eq!(
        SoftwareVersion::from_text("\n4.2 ignored", true).unwrap(),
        SoftwareVersion::new(4, 2, 0)
    );
}

#[test]
fn data_version_parses_orders_and_saves_text_like_upstream() {
    let mut out = Vec::new();
    let version = DataVersion::from_text(" 2026-03-24.1 trailing").unwrap();
    version.save_text(&mut out).unwrap();

    assert_eq!(version, DataVersion::new(2026, 3, 24, 1));
    assert_eq!(String::from_utf8(out).unwrap(), "2026-03-24.1");
    assert!(DataVersion::new(2025, 9, 22, 2) < version);
}

#[cfg(unix)]
fn make_executable(path: &std::path::Path) {
    use std::os::unix::fs::PermissionsExt;

    let mut permissions = std::fs::metadata(path).unwrap().permissions();
    permissions.set_mode(0o755);
    std::fs::set_permissions(path, permissions).unwrap();
}

#[cfg(unix)]
const SIGPIPE_NUM: i32 = 13;

#[cfg(unix)]
unsafe extern "C" {
    fn raise(sig: i32) -> i32;
}

#[cfg(unix)]
fn spawn_sigpipe_child(enable_sigpipe: bool) -> std::process::ExitStatus {
    Command::new(std::env::current_exe().unwrap())
        .arg("--exact")
        .arg("sigpipe_child_helper")
        .arg("--nocapture")
        .env("AMRFINDER_COMMON_SIGPIPE_CHILD", "1")
        .env(
            "AMRFINDER_COMMON_SIGPIPE_ENABLE",
            if enable_sigpipe { "1" } else { "0" },
        )
        .status()
        .unwrap()
}

fn spawn_shell_drop_child(quiet: bool) -> std::process::Output {
    Command::new(std::env::current_exe().unwrap())
        .arg("--exact")
        .arg("shell_drop_child_helper")
        .arg("--nocapture")
        .env("AMRFINDER_COMMON_SHELL_DROP_CHILD", "1")
        .env(
            "AMRFINDER_COMMON_SHELL_DROP_QUIET",
            if quiet { "1" } else { "0" },
        )
        .output()
        .unwrap()
}

#[cfg(unix)]
#[test]
fn sigpipe_child_helper() {
    if std::env::var_os("AMRFINDER_COMMON_SIGPIPE_CHILD").is_none() {
        return;
    }

    let mut app = app(0, Vec::new());
    let mut args = vec![OsString::from("prog")];
    if std::env::var_os("AMRFINDER_COMMON_SIGPIPE_ENABLE").as_deref() == Some(OsStr::new("1")) {
        args.push(OsString::from("-sigpipe"));
    }
    app.run(args).unwrap();

    unsafe {
        raise(SIGPIPE_NUM);
    }
    unreachable!("SIGPIPE handler should exit the process");
}

#[test]
fn shell_drop_child_helper() {
    if std::env::var_os("AMRFINDER_COMMON_SHELL_DROP_CHILD").is_none() {
        return;
    }

    let quiet =
        std::env::var_os("AMRFINDER_COMMON_SHELL_DROP_QUIET").as_deref() == Some(OsStr::new("1"));
    let mut common_app = app(
        0,
        vec![Key {
            name: "quiet",
            flag: true,
            default_value: "false",
            acronym: None,
        }],
    );
    let mut args = vec![OsString::from("shell_drop_prog")];
    if quiet {
        args.push(OsString::from("-quiet"));
    }
    common_app.run(args).unwrap();

    let mut shell = ShellApplication {
        application: app(0, Vec::new()),
        use_tmp: false,
        tmp: None,
        stderr_quiet: false,
        start_time: None,
    };
    shell.init_var().unwrap();
}

#[test]
fn application_rejects_duplicate_key_names_like_cpp_add_key() {
    let mut app = Application {
        description: "duplicate key test",
        usage: "duplicate key test",
        positionals: 0,
        needs_arg: false,
        threads_used: false,
        key_args: vec![
            Key {
                name: "same",
                flag: false,
                default_value: "",
                acronym: None,
            },
            Key {
                name: "same",
                flag: true,
                default_value: "false",
                acronym: None,
            },
        ],
    };

    assert!(matches!(app.run(os_args(&["prog"])), Err(1)));
}

#[test]
fn application_rejects_duplicate_key_acronyms_like_cpp_add_key() {
    let mut app = Application {
        description: "duplicate acronym test",
        usage: "duplicate acronym test",
        positionals: 0,
        needs_arg: false,
        threads_used: false,
        key_args: vec![
            Key {
                name: "quiet",
                flag: true,
                default_value: "false",
                acronym: Some('q'),
            },
            Key {
                name: "query",
                flag: false,
                default_value: "",
                acronym: Some('q'),
            },
        ],
    };

    assert!(matches!(app.run(os_args(&["prog"])), Err(1)));
}

#[test]
fn application_rejects_reserved_key_names_like_cpp_arg_qc() {
    let mut app = Application {
        description: "reserved key test",
        usage: "reserved key test",
        positionals: 0,
        needs_arg: false,
        threads_used: false,
        key_args: vec![Key {
            name: "help",
            flag: false,
            default_value: "",
            acronym: None,
        }],
    };

    assert!(matches!(app.run(os_args(&["prog"])), Err(1)));
}

#[test]
fn application_rejects_bad_key_names_like_cpp_arg_qc() {
    let mut app = Application {
        description: "bad key test",
        usage: "bad key test",
        positionals: 0,
        needs_arg: false,
        threads_used: false,
        key_args: vec![Key {
            name: "-bad",
            flag: false,
            default_value: "",
            acronym: None,
        }],
    };

    assert!(matches!(app.run(os_args(&["prog"])), Err(1)));
}

#[test]
fn application_reports_zero_seed_as_cpp_run_does() {
    let mut app = app(0, Vec::new());

    assert!(matches!(app.run(os_args(&["prog", "-seed", "0"])), Err(1)));
}

#[test]
fn application_adds_default_args_once() {
    let mut app = app(0, Vec::new());

    let first = app.run(os_args(&["prog"])).unwrap();
    let second = app.run(os_args(&["prog"])).unwrap();

    assert_eq!(first.key_args["seed"], "1");
    assert_eq!(first.key_args["verbose"], "0");
    assert_eq!(first.key_args["qc"], "false");
    assert_eq!(second.key_args["seed"], "1");
    assert_eq!(
        app.key_args.iter().filter(|key| key.name == "seed").count(),
        1
    );
}

#[test]
fn application_substitutes_base_in_default_values_like_cpp_init_environment() {
    let mut app = app(
        0,
        vec![Key {
            name: "db",
            flag: false,
            default_value: "$BASE/data",
            acronym: None,
        }],
    );

    let run = app.run(os_args(&["/bin/sh"])).unwrap();
    let expected_base = std::fs::canonicalize("/bin/sh")
        .unwrap()
        .parent()
        .unwrap()
        .join("data");

    assert_eq!(PathBuf::from(&run.key_args["db"]), expected_base);
}

#[test]
fn application_log_arg_opens_and_appends_header_on_run() {
    let dir = tempfile::tempdir().unwrap();
    let log = dir.path().join("common.log");
    let mut app = app(0, Vec::new());

    let run = app
        .run(vec![
            OsString::from("prog"),
            OsString::from("-log"),
            log.as_os_str().to_os_string(),
        ])
        .unwrap();

    assert_eq!(PathBuf::from(&run.key_args["log"]), log);
    let text = std::fs::read_to_string(log).unwrap();
    assert!(text.contains("$ prog -log "));
}

#[test]
fn application_accepts_gnu_double_dash_key_value_like_cpp_run() {
    let dir = tempfile::tempdir().unwrap();
    let log = dir.path().join("gnu-common.log");
    let mut app = app(0, Vec::new());

    let run = app
        .run(vec![
            OsString::from("prog"),
            OsString::from("--log"),
            log.as_os_str().to_os_string(),
        ])
        .unwrap();

    assert_eq!(PathBuf::from(&run.key_args["log"]), log);
    assert!(std::fs::read_to_string(log)
        .unwrap()
        .contains("$ prog --log "));
}

#[test]
fn application_accepts_gnu_double_dash_key_equals_value_like_cpp_run() {
    let dir = tempfile::tempdir().unwrap();
    let log = dir.path().join("gnu-equals-common.log");
    let mut app = app(0, Vec::new());

    let run = app
        .run(vec![
            OsString::from("prog"),
            OsString::from(format!("--log={}", log.display())),
        ])
        .unwrap();

    assert_eq!(PathBuf::from(&run.key_args["log"]), log);
    assert!(std::fs::read_to_string(log)
        .unwrap()
        .contains("$ prog --log "));
}

#[test]
fn application_accepts_gnu_double_dash_help_and_version_keywords() {
    let mut help_app = app(0, Vec::new());
    assert!(matches!(help_app.run(os_args(&["prog", "--help"])), Err(0)));

    let mut version_app = app(0, Vec::new());
    assert!(matches!(
        version_app.run(os_args(&["prog", "--version"])),
        Err(0)
    ));
}

#[test]
fn application_accepts_key_acronyms_like_cpp_char2arg() {
    let mut app = app(
        0,
        vec![
            Key {
                name: "quiet",
                flag: true,
                default_value: "false",
                acronym: Some('q'),
            },
            Key {
                name: "database",
                flag: false,
                default_value: "",
                acronym: Some('d'),
            },
        ],
    );

    let run = app.run(os_args(&["prog", "-q", "-d", "db"])).unwrap();

    assert_eq!(run.key_args["quiet"], "true");
    assert_eq!(run.key_args["database"], "db");
}

#[test]
fn command_line_quotes_shell_sensitive_arguments() {
    let mut app = app(1, Vec::new());

    app.run(os_args(&["prog", "two words"])).unwrap();

    assert_eq!(get_command_line(), "prog 'two words'");
}

#[test]
fn version_default_matches_cpp_application_default() {
    assert_eq!(DEFAULT_VERSION, "0.0.0");
}

#[test]
fn make_temp_dir_uses_program_name_prefix_and_directory_is_removable() {
    let mut app = app(0, Vec::new());
    app.run(os_args(&["common_prog"])).unwrap();

    let tmp = make_temp_dir().unwrap();
    let name = tmp.file_name().unwrap().to_string_lossy();
    assert!(name.starts_with("common_prog."));
    assert!(tmp.is_dir());
    std::fs::remove_dir(tmp).unwrap();
}

#[test]
fn shell_application_drop_removes_temp_directory_without_log() {
    let tmp = make_temp_dir().unwrap();
    let tmp_copy = tmp.clone();
    {
        let _shell = ShellApplication {
            application: app(0, Vec::new()),
            use_tmp: true,
            tmp: Some(tmp),
            stderr_quiet: false,
            start_time: None,
        };
    }

    assert!(!tmp_copy.exists());
}

#[test]
fn shell_application_drop_reports_elapsed_time_like_cpp_destructor() {
    let output = spawn_shell_drop_child(false);

    assert!(output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("Running: shell_drop_prog"));
    assert!(stderr.contains("shell_drop_prog took "));
    assert!(stderr.contains(" seconds to complete"));
}

#[test]
fn shell_application_drop_honors_quiet_stderr_like_cpp_destructor() {
    let output = spawn_shell_drop_child(true);

    assert!(output.status.success());
    assert!(String::from_utf8(output.stderr).unwrap().is_empty());
}

#[test]
fn exec_reports_nonzero_status_like_cpp_exec() {
    let err = exec("exit 7", "").unwrap_err();

    assert!(err.contains("exit 7"));
    assert!(err.contains("status = "));
}

#[test]
fn application_updates_common_runtime_flags_like_cpp_run() {
    let mut app = app(
        0,
        vec![Key {
            name: "threads",
            flag: false,
            default_value: "1",
            acronym: None,
        }],
    );

    app.run(os_args(&[
        "prog",
        "-qc",
        "-noprogress",
        "-sigpipe",
        "-threads",
        "3",
    ]))
    .unwrap();

    assert!(qc_on());
    assert!(!progress_enabled());
    assert!(sigpipe_enabled());
    assert!(!profile_enabled());
    assert_eq!(get_threads_max(), 3);
}

#[test]
fn application_threads_used_adds_default_threads_key_like_cpp() {
    let mut app = Application {
        description: "threads used test",
        usage: "threads used test",
        positionals: 0,
        needs_arg: false,
        threads_used: true,
        key_args: Vec::new(),
    };

    let run = app.run(os_args(&["prog", "-threads", "4"])).unwrap();

    assert_eq!(run.key_args["threads"], "4");
    assert_eq!(get_threads_max(), 4);
    assert_eq!(
        app.key_args
            .iter()
            .filter(|key| key.name == "threads")
            .count(),
        1
    );
}

#[test]
fn application_needs_arg_prints_usage_on_empty_invocation_like_cpp() {
    let mut app = Application {
        description: "needs arg test",
        usage: "needs arg usage",
        positionals: 0,
        needs_arg: true,
        threads_used: false,
        key_args: Vec::new(),
    };

    assert!(matches!(app.run(os_args(&["prog"])), Err(1)));
}

#[cfg(unix)]
#[test]
fn application_sigpipe_handler_exits_zero_when_sigpipe_flag_is_set_like_cpp() {
    let status = spawn_sigpipe_child(true);

    assert_eq!(status.code(), Some(0));
}

#[cfg(unix)]
#[test]
fn application_sigpipe_handler_exits_one_without_sigpipe_flag_like_cpp() {
    let status = spawn_sigpipe_child(false);

    assert_eq!(status.code(), Some(1));
}

#[test]
fn application_rejects_profile_with_multiple_threads_like_cpp_run() {
    let mut app = app(
        0,
        vec![Key {
            name: "threads",
            flag: false,
            default_value: "1",
            acronym: None,
        }],
    );

    assert!(matches!(
        app.run(os_args(&["prog", "-profile", "-threads", "2"])),
        Err(1)
    ));
}

#[test]
fn application_accepts_one_required_group_key_like_cpp_run() {
    let mut app = app(
        0,
        vec![
            Key {
                name: "prot",
                flag: false,
                default_value: "",
                acronym: None,
            },
            Key {
                name: "nuc",
                flag: false,
                default_value: "",
                acronym: None,
            },
        ],
    );
    app.set_required_group("prot", "input").unwrap();
    app.set_required_group("nuc", "input").unwrap();

    let run = app.run(os_args(&["prog", "-prot", "in.fa"])).unwrap();

    assert_eq!(run.key_args["prot"], "in.fa");
    assert_eq!(run.key_args["nuc"], "");
}

#[test]
fn application_rejects_missing_required_group_like_cpp_run() {
    let mut app = app(
        0,
        vec![
            Key {
                name: "prot",
                flag: false,
                default_value: "",
                acronym: None,
            },
            Key {
                name: "nuc",
                flag: false,
                default_value: "",
                acronym: None,
            },
        ],
    );
    app.set_required_group("prot", "input").unwrap();
    app.set_required_group("nuc", "input").unwrap();

    assert!(matches!(app.run(os_args(&["prog"])), Err(1)));
}

#[test]
fn application_rejects_conflicting_required_group_keys_like_cpp_run() {
    let mut app = app(
        0,
        vec![
            Key {
                name: "prot",
                flag: false,
                default_value: "",
                acronym: None,
            },
            Key {
                name: "nuc",
                flag: false,
                default_value: "",
                acronym: None,
            },
        ],
    );
    app.set_required_group("prot", "input").unwrap();
    app.set_required_group("nuc", "input").unwrap();

    assert!(matches!(
        app.run(os_args(&["prog", "-prot", "in.fa", "-nuc", "in.fna"])),
        Err(1)
    ));
}

#[test]
fn application_sets_quiet_runtime_flag_when_quiet_key_is_used() {
    let mut quiet_app = app(
        0,
        vec![Key {
            name: "quiet",
            flag: true,
            default_value: "false",
            acronym: None,
        }],
    );

    quiet_app.run(os_args(&["prog", "--quiet"])).unwrap();

    assert!(quiet_enabled());

    let mut plain_app = app(0, Vec::new());
    plain_app.run(os_args(&["prog"])).unwrap();
    assert!(!quiet_enabled());
}

#[test]
fn shell_application_init_var_honors_quiet_flag_like_cpp_stderr() {
    let mut quiet_app = app(
        0,
        vec![Key {
            name: "quiet",
            flag: true,
            default_value: "false",
            acronym: None,
        }],
    );
    quiet_app.run(os_args(&["prog", "--quiet"])).unwrap();
    let mut shell = ShellApplication {
        application: app(0, Vec::new()),
        use_tmp: true,
        tmp: None,
        stderr_quiet: false,
        start_time: None,
    };

    shell.init_var().unwrap();

    assert!(quiet_enabled());
    assert!(shell.tmp.as_ref().unwrap().exists());
}

#[test]
fn shell_application_finds_and_quotes_program_like_cpp_find_prog_full_prog() {
    let dir = tempfile::tempdir().unwrap();
    let bin = dir.path().join("bin");
    std::fs::create_dir(&bin).unwrap();
    let prog = bin.join("prog");
    let tool = bin.join("tool name");
    std::fs::write(&prog, "#!/bin/sh\nexit 0\n").unwrap();
    std::fs::write(&tool, "#!/bin/sh\nexit 0\n").unwrap();
    #[cfg(unix)]
    {
        make_executable(&prog);
        make_executable(&tool);
    }

    let mut common_app = app(0, Vec::new());
    common_app
        .run(vec![prog.as_os_str().to_os_string()])
        .unwrap();
    let shell = ShellApplication {
        application: app(0, Vec::new()),
        use_tmp: false,
        tmp: None,
        stderr_quiet: false,
        start_time: None,
    };

    shell.find_prog("tool name").unwrap();
    assert_eq!(
        shell.full_prog("tool name").unwrap(),
        format!(
            "{} ",
            amrfinder::common::shell_quote(&tool.to_string_lossy())
        )
    );
}

#[test]
fn shell_application_exec2str_requires_one_output_line_like_cpp_exec2str() {
    let tmp = tempfile::tempdir().unwrap();
    let shell = ShellApplication {
        application: app(0, Vec::new()),
        use_tmp: true,
        tmp: Some(tmp.path().to_path_buf()),
        stderr_quiet: false,
        start_time: None,
    };

    assert_eq!(shell.exec2str("printf one", "one_line", "").unwrap(), "one");
    let err = shell
        .exec2str("printf 'one\\ntwo\\n'", "two_lines", "")
        .unwrap_err();
    assert!(err.contains("One line is expected"));
}

#[test]
fn shell_application_uncompress_gz_invokes_gunzip_and_returns_tmp_suffix_like_cpp() {
    let dir = tempfile::tempdir().unwrap();
    let plain = dir.path().join("input.fa");
    let gz = dir.path().join("input.fa.gz");
    std::fs::write(&plain, ">id\nACGT\n").unwrap();
    let output = std::process::Command::new("gzip")
        .arg("-c")
        .arg(&plain)
        .output()
        .unwrap();
    assert!(output.status.success());
    std::fs::write(&gz, output.stdout).unwrap();
    let tmp = tempfile::tempdir().unwrap();
    let shell = ShellApplication {
        application: app(0, Vec::new()),
        use_tmp: true,
        tmp: Some(tmp.path().to_path_buf()),
        stderr_quiet: false,
        start_time: None,
    };

    let out = shell
        .uncompress(
            &amrfinder::common::shell_quote(&gz.to_string_lossy()),
            "flat.fa",
        )
        .unwrap();

    let out_path = PathBuf::from(out.trim_matches('\''));
    assert_eq!(out_path, tmp.path().join("flat.fa"));
    assert_eq!(std::fs::read_to_string(out_path).unwrap(), ">id\nACGT\n");
}

#[test]
fn shell_application_uncompress_leaves_non_gz_quoted_name_unchanged_like_cpp() {
    let tmp = tempfile::tempdir().unwrap();
    let shell = ShellApplication {
        application: app(0, Vec::new()),
        use_tmp: true,
        tmp: Some(tmp.path().to_path_buf()),
        stderr_quiet: false,
        start_time: None,
    };

    assert_eq!(
        shell.uncompress("'plain.fa'", "flat.fa").unwrap(),
        "'plain.fa'"
    );
}

#[test]
fn shell_application_blast_thread_probe_uses_fake_help_without_real_blast() {
    let dir = tempfile::tempdir().unwrap();
    let bin = dir.path().join("bin");
    std::fs::create_dir(&bin).unwrap();
    let prog = bin.join("prog");
    let blast = bin.join("blastp");
    std::fs::write(&prog, "#!/bin/sh\nexit 0\n").unwrap();
    std::fs::write(
        &blast,
        "#!/bin/sh\nif [ \"$1\" = \"-help\" ]; then printf '  -num_threads NUM\\n  -mt_mode MODE\\n'; fi\n",
    )
    .unwrap();
    #[cfg(unix)]
    {
        make_executable(&prog);
        make_executable(&blast);
    }

    let mut common_app = app(
        0,
        vec![Key {
            name: "threads",
            flag: false,
            default_value: "1",
            acronym: None,
        }],
    );
    common_app
        .run(vec![
            prog.as_os_str().to_os_string(),
            OsString::from("-threads"),
            OsString::from("4"),
        ])
        .unwrap();

    let tmp = tempfile::tempdir().unwrap();
    let shell = ShellApplication {
        application: app(0, Vec::new()),
        use_tmp: true,
        tmp: Some(tmp.path().to_path_buf()),
        stderr_quiet: false,
        start_time: None,
    };
    shell.find_prog("blastp").unwrap();

    assert_eq!(
        shell.get_blast_threads_param("blastp", 2).unwrap(),
        "  -num_threads 2  -mt_mode 1"
    );
}

#[test]
fn shell_application_blast_thread_probe_omits_mt_mode_when_help_lacks_it() {
    let dir = tempfile::tempdir().unwrap();
    let bin = dir.path().join("bin");
    std::fs::create_dir(&bin).unwrap();
    let prog = bin.join("prog-no-mt");
    let blast = bin.join("blastn-no-mt");
    std::fs::write(&prog, "#!/bin/sh\nexit 0\n").unwrap();
    std::fs::write(
        &blast,
        "#!/bin/sh\nif [ \"$1\" = \"-help\" ]; then printf '  -num_threads NUM\\n'; fi\n",
    )
    .unwrap();
    #[cfg(unix)]
    {
        make_executable(&prog);
        make_executable(&blast);
    }

    let mut common_app = app(
        0,
        vec![Key {
            name: "threads",
            flag: false,
            default_value: "1",
            acronym: None,
        }],
    );
    common_app
        .run(vec![
            prog.as_os_str().to_os_string(),
            OsString::from("-threads"),
            OsString::from("3"),
        ])
        .unwrap();

    let tmp = tempfile::tempdir().unwrap();
    let shell = ShellApplication {
        application: app(0, Vec::new()),
        use_tmp: true,
        tmp: Some(tmp.path().to_path_buf()),
        stderr_quiet: false,
        start_time: None,
    };
    shell.find_prog("blastn-no-mt").unwrap();

    assert_eq!(
        shell.get_blast_threads_param("blastn-no-mt", 8).unwrap(),
        "  -num_threads 3"
    );
}

#[test]
fn shell_application_blast_thread_probe_returns_empty_without_num_threads() {
    let dir = tempfile::tempdir().unwrap();
    let bin = dir.path().join("bin");
    std::fs::create_dir(&bin).unwrap();
    let prog = bin.join("prog-no-threads");
    let blast = bin.join("blastx-no-threads");
    std::fs::write(&prog, "#!/bin/sh\nexit 0\n").unwrap();
    std::fs::write(
        &blast,
        "#!/bin/sh\nif [ \"$1\" = \"-help\" ]; then printf '  -mt_mode MODE\\n'; fi\n",
    )
    .unwrap();
    #[cfg(unix)]
    {
        make_executable(&prog);
        make_executable(&blast);
    }

    let mut common_app = app(
        0,
        vec![Key {
            name: "threads",
            flag: false,
            default_value: "1",
            acronym: None,
        }],
    );
    common_app
        .run(vec![
            prog.as_os_str().to_os_string(),
            OsString::from("-threads"),
            OsString::from("3"),
        ])
        .unwrap();

    let tmp = tempfile::tempdir().unwrap();
    let shell = ShellApplication {
        application: app(0, Vec::new()),
        use_tmp: true,
        tmp: Some(tmp.path().to_path_buf()),
        stderr_quiet: false,
        start_time: None,
    };
    shell.find_prog("blastx-no-threads").unwrap();

    assert_eq!(
        shell
            .get_blast_threads_param("blastx-no-threads", 8)
            .unwrap(),
        ""
    );
}
