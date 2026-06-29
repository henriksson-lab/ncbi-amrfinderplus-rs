use std::ffi::OsString;
use std::fs;
use std::io::Write;
use std::path::Path;

use amrfinder::amrfinder_index::{main, ThisApplication};

#[test]
fn shell_body_indexes_core_and_mutation_taxgroup_databases() {
    let temp = tempfile::tempdir().unwrap();
    let db = temp.path().join("db");
    let blast = temp.path().join("blast");
    let hmmer = temp.path().join("hmmer");
    let log = temp.path().join("commands.log");
    fs::create_dir(&db).unwrap();
    fs::create_dir(&blast).unwrap();
    fs::create_dir(&hmmer).unwrap();
    fs::write(db.join("AMR.LIB"), b"hmm\n").unwrap();
    fs::write(db.join("AMRProt.fa"), b">p\n").unwrap();
    fs::write(db.join("AMR_CDS.fa"), b">n\n").unwrap();
    fs::write(db.join("AMR_DNA-Escherichia.fa"), b">d\n").unwrap();
    fs::write(
        db.join("taxgroup.tsv"),
        "# comment\nEscherichia\tgp\t2\nSalmonella\tgp\t0\n",
    )
    .unwrap();
    write_fake_tool(
        &hmmer.join("hmmpress"),
        &format!("printf 'hmmpress %s\\n' \"$*\" >> {}", log.display()),
    );
    write_fake_tool(
        &blast.join("makeblastdb"),
        &format!("printf 'makeblastdb %s\\n' \"$*\" >> {}", log.display()),
    );

    ThisApplication::new()
        .shell_body(&db, &blast, &hmmer)
        .unwrap();

    let commands = fs::read_to_string(log).unwrap();
    assert!(commands.contains("hmmpress -f "));
    assert!(commands.contains("AMR.LIB"));
    assert!(commands.contains("makeblastdb -in "));
    assert!(commands.contains("AMRProt.fa -dbtype prot"));
    assert!(commands.contains("AMR_CDS.fa -dbtype nucl"));
    assert!(commands.contains("AMR_DNA-Escherichia.fa -dbtype nucl"));
    assert!(!commands.contains("AMR_DNA-Salmonella.fa"));
}

#[test]
fn shell_body_rejects_missing_database_directory() {
    let temp = tempfile::tempdir().unwrap();

    let err = ThisApplication::new()
        .shell_body(&temp.path().join("missing"), Path::new(""), Path::new(""))
        .unwrap_err();

    assert!(err.to_string().contains("does not exist"));
}

#[test]
fn shell_body_rejects_taxgroup_with_space() {
    let temp = tempfile::tempdir().unwrap();
    let db = temp.path().join("db");
    fs::create_dir(&db).unwrap();
    fs::write(db.join("taxgroup.tsv"), "Bad taxon\tgp\t1\n").unwrap();

    let err = ThisApplication::new()
        .shell_body(&db, Path::new(""), Path::new(""))
        .unwrap_err();

    assert!(err.to_string().contains("space in taxgroup"));
}

#[test]
fn main_parses_original_keys_and_runs_shell_body() {
    let temp = tempfile::tempdir().unwrap();
    let db = temp.path().join("db");
    let blast = temp.path().join("blast");
    let hmmer = temp.path().join("hmmer");
    let log = temp.path().join("commands.log");
    fs::create_dir(&db).unwrap();
    fs::create_dir(&blast).unwrap();
    fs::create_dir(&hmmer).unwrap();
    fs::write(db.join("AMR.LIB"), b"hmm\n").unwrap();
    fs::write(db.join("AMRProt.fa"), b">p\n").unwrap();
    fs::write(db.join("AMR_CDS.fa"), b">n\n").unwrap();
    fs::write(db.join("taxgroup.tsv"), "Escherichia\tgp\t0\n").unwrap();
    write_fake_tool(
        &hmmer.join("hmmpress"),
        &format!("printf 'hmmpress %s\\n' \"$*\" >> {}", log.display()),
    );
    write_fake_tool(
        &blast.join("makeblastdb"),
        &format!("printf 'makeblastdb %s\\n' \"$*\" >> {}", log.display()),
    );

    assert_eq!(
        main(vec![
            OsString::from("amrfinder_index"),
            db.into_os_string(),
            OsString::from("-blast_bin"),
            blast.into_os_string(),
            OsString::from("-hmmer_bin"),
            hmmer.into_os_string(),
        ])
        .unwrap(),
        0
    );
    assert_eq!(
        fs::read_to_string(log)
            .unwrap()
            .lines()
            .filter(|line| line.starts_with("makeblastdb "))
            .count(),
        2
    );
}

fn write_fake_tool(path: &Path, body: &str) {
    let mut file = fs::File::create(path).unwrap();
    writeln!(file, "#!/bin/sh").unwrap();
    writeln!(file, "{body}").unwrap();
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        fs::set_permissions(path, fs::Permissions::from_mode(0o755)).unwrap();
    }
}
