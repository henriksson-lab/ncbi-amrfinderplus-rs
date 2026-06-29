use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{bail, ensure, Context, Result};
use tempfile::TempDir;

use crate::common::{exec, Application, Key};

pub struct ThisApplication {
    application: Application,
}

impl ThisApplication {
    pub fn new() -> Self {
        Self {
            application: Application {
                description: "Index the database for AMRFinder",
                usage: "amrfinder_index DATABASE [-blast_bin BLAST_DIR] [-hmmer_bin HMMER_DIR]",
                positionals: 1,
                needs_arg: false,
                threads_used: false,
                key_args: vec![
                    Key {
                        name: "blast_bin",
                        flag: false,
                        default_value: "",
                        acronym: None,
                    },
                    Key {
                        name: "hmmer_bin",
                        flag: false,
                        default_value: "",
                        acronym: None,
                    },
                ],
            },
        }
    }

    pub fn shell_body(&self, db_dir: &Path, blast_bin: &Path, hmmer_bin: &Path) -> Result<()> {
        let mut db_dir = db_dir.to_string_lossy().into_owned();
        let mut blast_bin = blast_bin.to_string_lossy().into_owned();
        let mut hmmer_bin = hmmer_bin.to_string_lossy().into_owned();
        if !db_dir.is_empty() && !db_dir.ends_with('/') {
            db_dir.push('/');
        }
        if !blast_bin.is_empty() && !blast_bin.ends_with('/') {
            blast_bin.push('/');
        }
        if !hmmer_bin.is_empty() && !hmmer_bin.ends_with('/') {
            hmmer_bin.push('/');
        }

        if !Path::new(&db_dir).is_dir() {
            bail!("Database directory {db_dir} does not exist");
        }

        let exec_dir = std::env::current_exe()
            .ok()
            .and_then(|path| path.parent().map(Path::to_path_buf))
            .and_then(|path| fs::canonicalize(path).ok());

        let mut makeblastdb_dir = blast_bin;
        if makeblastdb_dir.is_empty() {
            if let Some(dir) = exec_dir
                .as_ref()
                .filter(|dir| dir.join("makeblastdb").is_file())
            {
                makeblastdb_dir = format!("{}/", dir.display());
            } else if let Some(dir) = std::env::var_os("PATH").and_then(|path| {
                std::env::split_paths(&path)
                    .find(|dir| !dir.as_os_str().is_empty() && dir.join("makeblastdb").is_file())
            }) {
                makeblastdb_dir = format!("{}/", dir.display());
            } else {
                bail!(
                    "Binary 'makeblastdb' is not found.\nPlease make sure that 'makeblastdb' is in the same directory as 'amrfinder_index' or is in your $PATH."
                );
            }
        }
        let mut hmmpress_dir = hmmer_bin;
        if hmmpress_dir.is_empty() {
            if let Some(dir) = exec_dir
                .as_ref()
                .filter(|dir| dir.join("hmmpress").is_file())
            {
                hmmpress_dir = format!("{}/", dir.display());
            } else if let Some(dir) = std::env::var_os("PATH").and_then(|path| {
                std::env::split_paths(&path)
                    .find(|dir| !dir.as_os_str().is_empty() && dir.join("hmmpress").is_file())
            }) {
                hmmpress_dir = format!("{}/", dir.display());
            } else {
                bail!(
                    "Binary 'hmmpress' is not found.\nPlease make sure that 'hmmpress' is in the same directory as 'amrfinder_index' or is in your $PATH."
                );
            }
        }

        let shell_quote = |s: &str| format!("'{}'", s.replace('\'', "'\"'\"'"));
        let makeblastdb = format!("{} ", shell_quote(&(makeblastdb_dir + "makeblastdb")));
        let hmmpress = format!("{} ", shell_quote(&(hmmpress_dir + "hmmpress")));

        let tmp = TempDir::new().context("failed to create temporary directory")?;
        let tmp_s = tmp.path().to_string_lossy();

        let taxgroup_path = format!("{db_dir}taxgroup.tsv");
        let taxgroup_text = fs::read_to_string(&taxgroup_path)
            .with_context(|| format!("failed to read {taxgroup_path}"))?;
        let mut dna_point_muts = Vec::new();
        for line in taxgroup_text.lines() {
            if line.starts_with('#') {
                continue;
            }
            let mut taxgroup = line.to_string();
            let Some(pos) = taxgroup.rfind('\t') else {
                bail!("bad taxgroup line: {line}");
            };
            let n: i32 = taxgroup[pos + 1..]
                .parse()
                .with_context(|| format!("bad point mutation count in {line}"))?;
            taxgroup.truncate(pos);
            let Some(pos) = taxgroup.rfind('\t') else {
                bail!("bad taxgroup line: {line}");
            };
            taxgroup.truncate(pos);
            ensure!(n >= 0, "negative point mutation count in {line}");
            ensure!(!taxgroup.contains(' '), "space in taxgroup {taxgroup}");
            if n != 0 {
                dna_point_muts.push(taxgroup);
            }
        }

        eprintln!("Indexing");
        let hmmpress_err = format!("{tmp_s}/hmmpress.err");
        exec(
            &format!(
                "{hmmpress}-f {} > /dev/null 2> {hmmpress_err}",
                shell_quote(&(db_dir.clone() + "AMR.LIB"))
            ),
            &hmmpress_err,
        )
        .map_err(anyhow::Error::msg)?;

        let db_link = format!("{tmp_s}/db");
        #[cfg(unix)]
        std::os::unix::fs::symlink(
            fs::canonicalize(&db_dir)
                .with_context(|| format!("failed to canonicalize {db_dir}"))?,
            &db_link,
        )
        .with_context(|| format!("failed to symlink {db_link}"))?;
        #[cfg(not(unix))]
        fs::create_dir_all(&db_link).with_context(|| format!("failed to create {db_link}"))?;

        let amr_prot_log = format!("{tmp_s}/makeblastdb.AMRProt");
        exec(
            &format!(
                "{makeblastdb}-in {db_link}/AMRProt.fa  -dbtype prot  -logfile {amr_prot_log}"
            ),
            &amr_prot_log,
        )
        .map_err(anyhow::Error::msg)?;

        let amr_cds_log = format!("{tmp_s}/makeblastdb.AMR_CDS");
        exec(
            &format!("{makeblastdb}-in {db_link}/AMR_CDS.fa  -dbtype nucl  -logfile {amr_cds_log}"),
            &amr_cds_log,
        )
        .map_err(anyhow::Error::msg)?;

        for dna_point_mut in dna_point_muts {
            let log = format!("{tmp_s}/makeblastdb.AMR_DNA-{dna_point_mut}");
            exec(
                &format!(
                    "{makeblastdb}-in {db_link}/AMR_DNA-{dna_point_mut}.fa  -dbtype nucl  -logfile {log}"
                ),
                &log,
            )
            .map_err(anyhow::Error::msg)?;
        }

        Ok(())
    }
}

pub fn main(argv: Vec<OsString>) -> Result<i32> {
    let mut app = ThisApplication::new();
    let run = match app.application.run(argv) {
        Ok(run) => run,
        Err(code) => return Ok(code),
    };
    let db_dir = PathBuf::from(&run.positional_args[0]);
    let blast_bin = PathBuf::from(run.key_args.get("blast_bin").expect("blast_bin defaulted"));
    let hmmer_bin = PathBuf::from(run.key_args.get("hmmer_bin").expect("hmmer_bin defaulted"));
    app.shell_body(&db_dir, &blast_bin, &hmmer_bin)?;
    Ok(0)
}
