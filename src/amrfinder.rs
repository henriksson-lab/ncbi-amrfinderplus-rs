// Main AMRFinder pipeline orchestrator — port of amrfinder.cpp
// Coordinates BLAST, HMM, and report generation

use std::ffi::OsString;
use std::fs::{self, File, OpenOptions};
use std::io::ErrorKind;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{bail, ensure, Result};

use crate::common::{Application, DataVersion, Key, Run, SoftwareVersion};
use crate::gff::GffType;
use crate::seq::Hsp;
use crate::tsv::{Header, TextTable};

const HMMSEARCH_EFFECTIVE_Z: &str = "10000";
const DATA_VER_MIN: &str = "2025-09-22.2";
const STXTYPER_VERSION: &str = "1.0.45";
const AMRFINDER_VERSION: &str = include_str!("../amr/version.txt");

#[derive(Debug, Clone, Copy)]
enum TblastnSearch {
    Fast,
    Slow,
}

impl TblastnSearch {
    fn base_args(self) -> Vec<&'static str> {
        let args = match self {
            Self::Fast => Hsp::BLASTP_FAST,
            Self::Slow => Hsp::BLASTP_SLOW,
        };
        args.split_whitespace().collect()
    }
}

#[derive(Debug, Default, PartialEq, Eq)]
struct ReportBlastxArgs {
    blastx_path: Option<PathBuf>,
    dna_len_path: Option<PathBuf>,
}

/// Pipeline configuration
#[derive(Debug, Clone)]
#[allow(clippy::struct_excessive_bools)]
pub struct PipelineConfig {
    pub protein: Option<PathBuf>,
    pub nucleotide: Option<PathBuf>,
    pub gff: Option<PathBuf>,
    pub database: PathBuf,
    pub organism: String,
    pub ident_min: f64,
    pub coverage_min: f64,
    pub threads: usize,
    pub plus: bool,
    pub report_common: bool,
    pub report_all_equal: bool,
    pub print_node: bool,
    pub pgap: bool,
    pub gpipe_org: bool,
    pub mutation_all: Option<PathBuf>,
    pub annotation_format: String,
    pub translation_table: u32,
    pub name: String,
    pub blast_bin: String,
    pub hmmer_bin: String,
    pub output: Option<PathBuf>,
    pub protein_output: Option<PathBuf>,
    pub nucleotide_output: Option<PathBuf>,
    pub nucleotide_flank5_output: Option<PathBuf>,
    pub nucleotide_flank5_size: usize,
    pub parm: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DatabaseVersion {
    pub directory: PathBuf,
    pub data_version: DataVersion,
}

pub struct ThisApplication {
    application: Application,
}

impl ThisApplication {
    pub fn new() -> Self {
        Self {
            application: Application {
                description: "Identify AMR genes and point mutations",
                usage: "amrfinder [-p PROT_FASTA] [-n NUCL_FASTA] [-g GFF] [-d DATABASE] [options]",
                positionals: 0,
                needs_arg: false,
                threads_used: true,
                key_args: vec![
                    Key {
                        name: "protein",
                        flag: false,
                        default_value: "",
                        acronym: Some('p'),
                    },
                    Key {
                        name: "nucleotide",
                        flag: false,
                        default_value: "",
                        acronym: Some('n'),
                    },
                    Key {
                        name: "gff",
                        flag: false,
                        default_value: "",
                        acronym: Some('g'),
                    },
                    Key {
                        name: "database",
                        flag: false,
                        default_value: "$BASE/data/latest",
                        acronym: Some('d'),
                    },
                    Key {
                        name: "dir",
                        flag: false,
                        default_value: "",
                        acronym: None,
                    },
                    Key {
                        name: "organism",
                        flag: false,
                        default_value: "",
                        acronym: Some('O'),
                    },
                    Key {
                        name: "ident_min",
                        flag: false,
                        default_value: "-1",
                        acronym: Some('i'),
                    },
                    Key {
                        name: "coverage_min",
                        flag: false,
                        default_value: "0.5",
                        acronym: Some('c'),
                    },
                    Key {
                        name: "threads",
                        flag: false,
                        default_value: "4",
                        acronym: None,
                    },
                    Key {
                        name: "plus",
                        flag: true,
                        default_value: "false",
                        acronym: None,
                    },
                    Key {
                        name: "report_common",
                        flag: true,
                        default_value: "false",
                        acronym: None,
                    },
                    Key {
                        name: "report_all_equal",
                        flag: true,
                        default_value: "false",
                        acronym: None,
                    },
                    Key {
                        name: "print_node",
                        flag: true,
                        default_value: "false",
                        acronym: None,
                    },
                    Key {
                        name: "pgap",
                        flag: true,
                        default_value: "false",
                        acronym: None,
                    },
                    Key {
                        name: "gpipe_org",
                        flag: true,
                        default_value: "false",
                        acronym: None,
                    },
                    Key {
                        name: "mutation_all",
                        flag: false,
                        default_value: "",
                        acronym: None,
                    },
                    Key {
                        name: "annotation_format",
                        flag: false,
                        default_value: "genbank",
                        acronym: Some('a'),
                    },
                    Key {
                        name: "translation_table",
                        flag: false,
                        default_value: "11",
                        acronym: Some('t'),
                    },
                    Key {
                        name: "name",
                        flag: false,
                        default_value: "",
                        acronym: None,
                    },
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
                    Key {
                        name: "output",
                        flag: false,
                        default_value: "",
                        acronym: Some('o'),
                    },
                    Key {
                        name: "protein_output",
                        flag: false,
                        default_value: "",
                        acronym: None,
                    },
                    Key {
                        name: "nucleotide_output",
                        flag: false,
                        default_value: "",
                        acronym: None,
                    },
                    Key {
                        name: "nucleotide_flank5_output",
                        flag: false,
                        default_value: "",
                        acronym: None,
                    },
                    Key {
                        name: "nucleotide_flank5_size",
                        flag: false,
                        default_value: "0",
                        acronym: None,
                    },
                    Key {
                        name: "parm",
                        flag: false,
                        default_value: "",
                        acronym: None,
                    },
                    Key {
                        name: "debug",
                        flag: true,
                        default_value: "false",
                        acronym: None,
                    },
                    Key {
                        name: "database_version",
                        flag: true,
                        default_value: "false",
                        acronym: Some('V'),
                    },
                    Key {
                        name: "list_organisms",
                        flag: true,
                        default_value: "false",
                        acronym: Some('l'),
                    },
                    Key {
                        name: "update",
                        flag: true,
                        default_value: "false",
                        acronym: Some('u'),
                    },
                    Key {
                        name: "force_update",
                        flag: true,
                        default_value: "false",
                        acronym: Some('U'),
                    },
                ],
            },
        }
    }

    pub fn shell_body(&mut self, run: &Run, out: &mut dyn Write) -> Result<()> {
        let dir = run.key_args["dir"].clone();
        let prepend_dir = |value: &str| -> PathBuf {
            if dir.is_empty() {
                PathBuf::from(value)
            } else {
                PathBuf::from(&dir).join(value)
            }
        };
        let protein = run.key_args["protein"]
            .is_empty()
            .then_some(None)
            .unwrap_or_else(|| Some(prepend_dir(&run.key_args["protein"])));
        let nucleotide = run.key_args["nucleotide"]
            .is_empty()
            .then_some(None)
            .unwrap_or_else(|| Some(prepend_dir(&run.key_args["nucleotide"])));
        let gff = run.key_args["gff"]
            .is_empty()
            .then_some(None)
            .unwrap_or_else(|| Some(prepend_dir(&run.key_args["gff"])));
        let _gff_type = GffType::name2type(&run.key_args["annotation_format"])?;
        let database_arg_given = crate::common::get_command_line()
            .split_whitespace()
            .any(|arg| matches!(arg, "-d" | "--database" | "-database"));
        let mut database = if !database_arg_given {
            std::env::var_os("AMRFINDER_DB")
                .map(PathBuf::from)
                .unwrap_or_else(|| PathBuf::from(&run.key_args["database"]))
        } else {
            PathBuf::from(&run.key_args["database"])
        };
        let update = run.key_args["update"] == "true" || run.key_args["force_update"] == "true";
        let software_dir = std::env::current_exe()
            .ok()
            .and_then(|path| path.parent().map(Path::to_path_buf))
            .unwrap_or_else(|| PathBuf::from("."));
        if run.key_args["database_version"] == "true" {
            writeln!(out, "Software directory: '{}'", software_dir.display())?;
            writeln!(out, "Software version: {}", AMRFINDER_VERSION.trim_end())?;
        } else {
            eprintln!("Software directory: {}", software_dir.display());
            eprintln!("Software version: {}", AMRFINDER_VERSION.trim_end());
        }

        if run.key_args["name"].contains('\t') {
            bail!("NAME cannot contain a tab character");
        }
        let ident_min: f64 = run.key_args["ident_min"].parse()?;
        let coverage_min: f64 = run.key_args["coverage_min"].parse()?;
        if ident_min != -1.0 && !(0.0..=1.0).contains(&ident_min) {
            bail!("ident_min must be between 0 and 1");
        }
        if !(0.0..=1.0).contains(&coverage_min) {
            bail!("coverage_min must be between 0 and 1");
        }
        if run.key_args["report_common"] == "true" && run.key_args["organism"].is_empty() {
            bail!("--report_common requires --organism");
        }
        if run.key_args["report_common"] == "true" && run.key_args["plus"] != "true" {
            bail!("--report_common requires --plus");
        }
        if !run.key_args["mutation_all"].is_empty() && run.key_args["organism"].is_empty() {
            eprintln!(
                "WARNING: --mutation_all option used without -O/--organism option. No point mutations will be screened"
            );
        }

        let output = run.key_args["output"]
            .is_empty()
            .then_some(None)
            .unwrap_or_else(|| Some(PathBuf::from(&run.key_args["output"])));
        if let Some(output_path) = &output {
            if File::create(output_path).is_err() {
                bail!("Cannot create output file '{}'", output_path.display());
            }
            fs::remove_file(output_path)?;
        }

        let mut blast_bin = run.key_args["blast_bin"].clone();
        if blast_bin.is_empty() {
            if let Some(env_blast_bin) = std::env::var_os("BLAST_BIN") {
                blast_bin = env_blast_bin.to_string_lossy().into_owned();
            }
        }

        if update {
            if protein.is_some() || nucleotide.is_some() {
                bail!("AMRFinder -u/--update option cannot be run with -n/--nucleotide or -p/--protein options");
            }
            if database_arg_given {
                bail!("AMRFinder update option (-u/--update) only operates on the default database directory. The -d/--database option is not permitted");
            }
            if std::env::var_os("AMRFINDER_DB").is_some() {
                eprintln!("Warning: AMRFINDER_DB is set, but AMRFinder auto-update only downloads to the default database directory");
                database = PathBuf::from(&run.key_args["database"]);
            }
            let mut argv = vec![OsString::from("amrfinder_update")];
            argv.push(OsString::from("-d"));
            let update_database = PathBuf::from(&run.key_args["database"]);
            if update_database
                .file_name()
                .is_some_and(|name| name == "latest")
            {
                argv.push(
                    update_database
                        .parent()
                        .map(PathBuf::from)
                        .unwrap_or(update_database)
                        .into(),
                );
                if run.key_args["force_update"] == "true" {
                    argv.push(OsString::from("--force_update"));
                }
                if !blast_bin.is_empty() {
                    argv.push(OsString::from("-blast_bin"));
                    argv.push(OsString::from(&blast_bin));
                }
                if !run.key_args["hmmer_bin"].is_empty() {
                    argv.push(OsString::from("-hmmer_bin"));
                    argv.push(OsString::from(&run.key_args["hmmer_bin"]));
                }
                if run.key_args["debug"] == "true" {
                    argv.push(OsString::from("--debug"));
                } else if run.key_args["qc"] == "true" {
                    argv.push(OsString::from("--qc"));
                }
                crate::amrfinder_update::main(argv)?;
            } else {
                let parent = update_database.parent().unwrap_or(Path::new(""));
                eprintln!(
                    "Warning: Updating database directory works only for databases with the default data directory format.\n         Please see https://github.com/ncbi/amr/wiki for details.\n         Current database directory is: {}\n         New database directories will be created as subdirectories of: {}",
                    update_database.display(),
                    parent.display()
                );
            }
        }

        if run.key_args["database_version"] == "true" {
            let download_latest_instr =
                "\nTo download the latest version to the default directory run: amrfinder -u";
            if !database.is_dir() {
                bail!(
                    "No valid AMRFinder database is found.\nThis directory (or symbolic link to directory) is not found: {}{}",
                    database.display(),
                    download_latest_instr
                );
            }
            let directory = fs::canonicalize(&database).unwrap_or_else(|_| database.clone());
            writeln!(out, "Database directory: '{}'", directory.display())?;

            let db_test = database.join("AMRProt.fa.phr");
            if !db_test.exists() {
                bail!(
                    "The BLAST database for AMRProt.fa ({}) was not found.\nUse amrfinder -u or amrfinder --force_update to download and prepare database for AMRFinderPlus",
                    db_test.display()
                );
            }

            let software_version = SoftwareVersion::from_text(AMRFINDER_VERSION.trim_end(), false)
                .map_err(anyhow::Error::msg)?;
            let software_version_min =
                SoftwareVersion::from_file(&database.join("database_format_version.txt"))
                    .map_err(anyhow::Error::msg)?;
            let data_version = DataVersion::from_file(&database.join("version.txt"))
                .map_err(anyhow::Error::msg)?;
            let data_version_min =
                DataVersion::from_text(DATA_VER_MIN).map_err(anyhow::Error::msg)?;

            writeln!(out, "Database version: {}", data_version)?;
            if software_version < software_version_min {
                bail!(
                    "Database requires software version at least {}",
                    software_version_min
                );
            }
            if data_version < data_version_min {
                bail!(
                    "Software requires database version at least {}{}",
                    data_version_min,
                    download_latest_instr
                );
            }
            if let Some(protein) = protein {
                bail!("No processing of {} is done", protein.display());
            }
            if let Some(nucleotide) = nucleotide {
                bail!("No processing of {} is done", nucleotide.display());
            }
            return Ok(());
        }

        if run.key_args["list_organisms"] == "true" {
            let database_info = database_version(&database)?;
            eprintln!("Database directory: {}", database_info.directory.display());
            eprintln!("Database version: {}", database_info.data_version);
            let organisms = db2organisms(&database)?;
            writeln!(
                out,
                "\nAvailable --organism options: {}",
                organisms.join(", ")
            )?;
            return Ok(());
        }

        if update {
            let _ = database_version(&database)?;
            if protein.is_none() && nucleotide.is_none() {
                return Ok(());
            }
        }

        let config = PipelineConfig {
            protein,
            nucleotide,
            gff,
            database,
            organism: run.key_args["organism"].clone(),
            ident_min,
            coverage_min,
            threads: run.key_args["threads"].parse()?,
            plus: run.key_args["plus"] == "true",
            report_common: run.key_args["report_common"] == "true",
            report_all_equal: run.key_args["report_all_equal"] == "true",
            print_node: run.key_args["print_node"] == "true",
            pgap: run.key_args["pgap"] == "true",
            gpipe_org: run.key_args["gpipe_org"] == "true",
            mutation_all: run.key_args["mutation_all"]
                .is_empty()
                .then_some(None)
                .unwrap_or_else(|| Some(PathBuf::from(&run.key_args["mutation_all"]))),
            annotation_format: run.key_args["annotation_format"].clone(),
            translation_table: run.key_args["translation_table"].parse()?,
            name: run.key_args["name"].clone(),
            blast_bin,
            hmmer_bin: run.key_args["hmmer_bin"].clone(),
            output: output.clone(),
            protein_output: run.key_args["protein_output"]
                .is_empty()
                .then_some(None)
                .unwrap_or_else(|| Some(PathBuf::from(&run.key_args["protein_output"]))),
            nucleotide_output: run.key_args["nucleotide_output"]
                .is_empty()
                .then_some(None)
                .unwrap_or_else(|| Some(PathBuf::from(&run.key_args["nucleotide_output"]))),
            nucleotide_flank5_output: run.key_args["nucleotide_flank5_output"]
                .is_empty()
                .then_some(None)
                .unwrap_or_else(|| Some(PathBuf::from(&run.key_args["nucleotide_flank5_output"]))),
            nucleotide_flank5_size: run.key_args["nucleotide_flank5_size"].parse()?,
            parm: run.key_args["parm"].clone(),
        };
        let output_text = run_pipeline(&config)?;
        if output.is_none() {
            write!(out, "{output_text}")?;
        }
        Ok(())
    }
}

pub fn main(argv: Vec<OsString>, out: &mut dyn Write) -> Result<i32> {
    let mut app = ThisApplication::new();
    let run = match app.application.run(argv) {
        Ok(run) => run,
        Err(code) => return Ok(code),
    };
    app.shell_body(&run, out)?;
    Ok(0)
}

impl Default for PipelineConfig {
    fn default() -> Self {
        PipelineConfig {
            protein: None,
            nucleotide: None,
            gff: None,
            database: PathBuf::new(),
            organism: String::new(),
            ident_min: -1.0,
            coverage_min: 0.5,
            threads: 4,
            plus: false,
            report_common: false,
            report_all_equal: false,
            print_node: false,
            pgap: false,
            gpipe_org: false,
            mutation_all: None,
            annotation_format: "genbank".to_string(),
            translation_table: 11,
            name: String::new(),
            blast_bin: String::new(),
            hmmer_bin: String::new(),
            output: None,
            protein_output: None,
            nucleotide_output: None,
            nucleotide_flank5_output: None,
            nucleotide_flank5_size: 0,
            parm: String::new(),
        }
    }
}

impl PipelineConfig {
    fn find_prog(&self, name: &str) -> PathBuf {
        let dir = match name {
            "blastp" | "blastn" | "blastx" | "tblastn" | "makeblastdb" => &self.blast_bin,
            "hmmsearch" | "hmmpress" => &self.hmmer_bin,
            _ => "",
        };
        if dir.is_empty() {
            PathBuf::from(name)
        } else {
            Path::new(dir).join(name)
        }
    }
}

fn fasta_check(
    f_name: &Path,
    prot: bool,
    ambig_max: usize,
    len_path: Option<&Path>,
    out_path: Option<&Path>,
) -> Result<(usize, usize, usize)> {
    if prot {
        crate::fasta_check::body(f_name, true, false, false, ambig_max, true, None, out_path)
    } else {
        crate::fasta_check::body(f_name, false, true, true, 0, false, len_path, None)
    }
}

/// Split a FASTA file using the translated fasta2parts body.
fn split_fasta(query: &Path, n: usize, out_dir: &Path) -> Result<Vec<PathBuf>> {
    crate::fasta2parts::body(query, n.max(2), out_dir)?;
    let mut paths = Vec::new();
    for part in 1..=n.max(2) {
        let path = out_dir.join(part.to_string());
        if path.exists() {
            paths.push(path);
        }
    }
    Ok(paths)
}

/// Run hmmsearch the way the C++ amrfinder does (amrfinder.cpp:1226-1243).
///
/// HMMER's internal `--cpu N` threading scales poorly on the AMR workload
/// (small query, ~780 models): measured ~4x slower than `--cpu 0`. So when
/// there is enough work to parallelize, split the query into `threads` chunks
/// and run up to `threads - 1` single-threaded (`--cpu 0`) hmmsearch processes
/// in parallel, then concatenate the per-chunk tblout/domtblout outputs (the
/// report parser skips the `#` comment blocks). Otherwise fall back to a single
/// process with `--cpu (threads-1)`.
#[allow(clippy::too_many_arguments)]
fn run_hmmsearch(
    hmm_prog: &Path,
    hmm_lib: &str,
    query: &Path,
    n_prot: usize,
    threads: usize,
    tblout: &Path,
    domtblout: &Path,
    work_dir: &Path,
) -> Result<()> {
    if threads > 1 && n_prot > threads / 2 {
        let chunk_dir = work_dir.join("hmm_chunk");
        fs::create_dir_all(&chunk_dir)?;
        let chunks = split_fasta(query, threads, &chunk_dir)?;

        let worker_limit = threads.saturating_sub(1).max(1);
        let mut finished = Vec::with_capacity(chunks.len());
        let mut handles = Vec::with_capacity(worker_limit);
        for (i, chunk) in chunks.into_iter().enumerate() {
            let prog = hmm_prog.to_path_buf();
            let lib = hmm_lib.to_string();
            let chunk_tbl = work_dir.join(format!("hmm_tbl_{i}"));
            let chunk_dom = work_dir.join(format!("hmm_dom_{i}"));
            let tbl = chunk_tbl.clone();
            let dom = chunk_dom.clone();
            let handle = std::thread::spawn(move || -> Result<()> {
                let status = Command::new(&prog)
                    .args([
                        "--noali",
                        "--cut_tc",
                        "-Z",
                        HMMSEARCH_EFFECTIVE_Z,
                        "--cpu",
                        "0",
                        "--tblout",
                        chunk_tbl.to_str().unwrap(),
                        "--domtblout",
                        chunk_dom.to_str().unwrap(),
                        &lib,
                        chunk.to_str().unwrap(),
                    ])
                    .stderr(std::process::Stdio::null())
                    .stdout(std::process::Stdio::null())
                    .status()?;
                if status.success() {
                    Ok(())
                } else {
                    bail!("hmmsearch failed")
                }
            });
            handles.push((tbl, dom, handle));
            if handles.len() == worker_limit {
                for (tbl, dom, handle) in handles.drain(..) {
                    handle
                        .join()
                        .map_err(|_| anyhow::anyhow!("hmmsearch thread panicked"))??;
                    finished.push((tbl, dom));
                }
            }
        }

        // Concatenate chunk outputs in order once each process finishes.
        let mut tbl_out = File::create(tblout)?;
        let mut dom_out = File::create(domtblout)?;
        for (tbl, dom) in finished {
            std::io::copy(&mut File::open(&tbl)?, &mut tbl_out)?;
            std::io::copy(&mut File::open(&dom)?, &mut dom_out)?;
        }
        for (tbl, dom, handle) in handles {
            handle
                .join()
                .map_err(|_| anyhow::anyhow!("hmmsearch thread panicked"))??;
            std::io::copy(&mut File::open(&tbl)?, &mut tbl_out)?;
            std::io::copy(&mut File::open(&dom)?, &mut dom_out)?;
        }
        Ok(())
    } else {
        let cpu = threads.saturating_sub(1).to_string();
        let status = Command::new(hmm_prog)
            .args([
                "--noali",
                "--cut_tc",
                "-Z",
                HMMSEARCH_EFFECTIVE_Z,
                "--cpu",
                &cpu,
                "--tblout",
                tblout.to_str().unwrap(),
                "--domtblout",
                domtblout.to_str().unwrap(),
                hmm_lib,
                query.to_str().unwrap(),
            ])
            .stderr(std::process::Stdio::null())
            .stdout(std::process::Stdio::null())
            .status()?;
        if status.success() {
            Ok(())
        } else {
            bail!("hmmsearch failed")
        }
    }
}

/// Run tblastn the way the C++ amrfinder does (amrfinder.cpp:1283-1300).
///
/// `tblastn -subject` does not benefit from `-num_threads`, so when more than
/// one thread is available, split the protein query into `threads` chunks and
/// run up to `threads - 1` tblastn processes in parallel, then concatenate the
/// outputs (outfmt 6 has no comment lines). Each query sequence is scored
/// independently against the fixed `-dbsize`, so splitting does not change
/// results.
fn run_tblastn(
    blast_prog: &Path,
    query: &str,
    subject: &Path,
    gencode: u32,
    threads: usize,
    out: &Path,
    work_dir: &Path,
) -> Result<()> {
    run_tblastn_search(
        TblastnSearch::Fast,
        blast_prog,
        query,
        subject,
        gencode,
        threads,
        out,
        work_dir,
    )
}

#[allow(clippy::too_many_arguments)]
fn run_tblastn_search(
    search: TblastnSearch,
    blast_prog: &Path,
    query: &str,
    subject: &Path,
    gencode: u32,
    threads: usize,
    out: &Path,
    work_dir: &Path,
) -> Result<()> {
    let gencode = gencode.to_string();
    let run_one = |query: &str, out: &Path| -> Result<()> {
        let mut args: Vec<String> = ["-query", query, "-subject", subject.to_str().unwrap()]
            .into_iter()
            .map(str::to_string)
            .chain(search.base_args().into_iter().map(str::to_string))
            .collect();
        if matches!(search, TblastnSearch::Fast) {
            args.extend(["-task".to_string(), "tblastn-fast".to_string()]);
        }
        args.extend([
            "-window_size".to_string(),
            "15".to_string(),
            "-threshold".to_string(),
            "100".to_string(),
            "-db_gencode".to_string(),
            gencode.clone(),
            "-outfmt".to_string(),
            "6 qseqid sseqid qstart qend qlen sstart send slen qseq sseq".to_string(),
            "-out".to_string(),
            out.to_str().unwrap().to_string(),
        ]);
        let status = Command::new(blast_prog)
            .args(&args)
            .stderr(std::process::Stdio::null())
            .stdout(std::process::Stdio::null())
            .status()?;
        if status.success() {
            Ok(())
        } else {
            bail!("tblastn failed")
        }
    };

    if threads <= 1 {
        return run_one(query, out);
    }

    let chunk_dir = work_dir.join("AMRProt_chunk");
    fs::create_dir_all(&chunk_dir)?;
    let chunks = split_fasta(Path::new(query), threads, &chunk_dir)?;

    let worker_limit = threads.saturating_sub(1).max(1);
    let mut finished = Vec::with_capacity(chunks.len());
    let mut handles = Vec::with_capacity(worker_limit);
    for (i, chunk) in chunks.into_iter().enumerate() {
        let prog = blast_prog.to_path_buf();
        let subject = subject.to_path_buf();
        let gencode = gencode.clone();
        let chunk_out = work_dir.join(format!("tblastn_{i}"));
        let out_ret = chunk_out.clone();
        let handle = std::thread::spawn(move || -> Result<()> {
            let mut args: Vec<String> = [
                "-query",
                chunk.to_str().unwrap(),
                "-subject",
                subject.to_str().unwrap(),
            ]
            .into_iter()
            .map(str::to_string)
            .chain(search.base_args().into_iter().map(str::to_string))
            .collect();
            if matches!(search, TblastnSearch::Fast) {
                args.extend(["-task".to_string(), "tblastn-fast".to_string()]);
            }
            args.extend([
                "-window_size".to_string(),
                "15".to_string(),
                "-threshold".to_string(),
                "100".to_string(),
                "-db_gencode".to_string(),
                gencode,
                "-outfmt".to_string(),
                "6 qseqid sseqid qstart qend qlen sstart send slen qseq sseq".to_string(),
                "-out".to_string(),
                chunk_out.to_str().unwrap().to_string(),
            ]);
            let status = Command::new(&prog)
                .args(&args)
                .stderr(std::process::Stdio::null())
                .stdout(std::process::Stdio::null())
                .status()?;
            if status.success() {
                Ok(())
            } else {
                bail!("tblastn failed")
            }
        });
        handles.push((out_ret, handle));
        if handles.len() == worker_limit {
            for (chunk_out, handle) in handles.drain(..) {
                handle
                    .join()
                    .map_err(|_| anyhow::anyhow!("tblastn thread panicked"))??;
                finished.push(chunk_out);
            }
        }
    }

    let mut out_file = File::create(out)?;
    for chunk_out in finished {
        std::io::copy(&mut File::open(&chunk_out)?, &mut out_file)?;
    }
    for (chunk_out, handle) in handles {
        handle
            .join()
            .map_err(|_| anyhow::anyhow!("tblastn thread panicked"))??;
        std::io::copy(&mut File::open(&chunk_out)?, &mut out_file)?;
    }
    Ok(())
}

fn parse_report_blastx_args(args: &str) -> ReportBlastxArgs {
    let mut parsed = ReportBlastxArgs::default();
    let parts: Vec<&str> = args.split_whitespace().collect();
    let mut i = 0;
    while i < parts.len() {
        match parts[i] {
            "-blastx" => {
                i += 1;
                if i < parts.len() {
                    parsed.blastx_path = Some(PathBuf::from(parts[i]));
                }
            }
            "-dna_len" => {
                i += 1;
                if i < parts.len() {
                    parsed.dna_len_path = Some(PathBuf::from(parts[i]));
                }
            }
            _ => {}
        }
        i += 1;
    }
    parsed
}

pub fn db2organisms(db: &Path) -> Result<Vec<String>> {
    let taxgroup = TextTable::from_file(db.join("taxgroup.tsv").to_str().unwrap())?;
    let amrprot_mutation = TextTable::from_file(db.join("AMRProt-mutation.tsv").to_str().unwrap())?;
    let amrprot_susceptible =
        TextTable::from_file(db.join("AMRProt-susceptible.tsv").to_str().unwrap())?;
    taxgroup.qc()?;
    amrprot_mutation.qc()?;
    amrprot_susceptible.qc()?;
    let mut vec = taxgroup.col2values(0);
    vec.extend(amrprot_mutation.col2values(0));
    vec.extend(amrprot_susceptible.col2values(0));
    vec.sort();
    vec.dedup();
    Ok(vec)
}

pub fn database_version(db: &Path) -> Result<DatabaseVersion> {
    let download_latest_instr =
        "\nTo download the latest version to the default directory run: amrfinder -u";
    if !db.is_dir() {
        bail!(
            "No valid AMRFinder database is found.\nThis directory (or symbolic link to directory) is not found: {}{}",
            db.display(),
            download_latest_instr
        );
    }

    let db_test = db.join("AMRProt.fa.phr");
    if !db_test.exists() {
        bail!(
            "The BLAST database for AMRProt.fa ({}) was not found.\nUse amrfinder -u or amrfinder --force_update to download and prepare database for AMRFinderPlus",
            db_test.display()
        );
    }

    let software_version = SoftwareVersion::from_text(AMRFINDER_VERSION.trim_end(), false)
        .map_err(anyhow::Error::msg)?;
    let software_version_min = SoftwareVersion::from_file(&db.join("database_format_version.txt"))
        .map_err(anyhow::Error::msg)?;
    let data_version =
        DataVersion::from_file(&db.join("version.txt")).map_err(anyhow::Error::msg)?;
    let data_version_min = DataVersion::from_text(DATA_VER_MIN).map_err(anyhow::Error::msg)?;

    if software_version < software_version_min {
        bail!(
            "Database requires software version at least {}",
            software_version_min
        );
    }
    if data_version < data_version_min {
        bail!(
            "Software requires database version at least {}{}",
            data_version_min,
            download_latest_instr
        );
    }

    Ok(DatabaseVersion {
        directory: fs::canonicalize(db).unwrap_or_else(|_| db.to_path_buf()),
        data_version,
    })
}

/// Run the AMRFinder pipeline
/// This initially shells out to external tools (BLAST, HMM, amr_report)
/// matching the C++ behavior exactly
pub fn run_pipeline(config: &PipelineConfig) -> Result<String> {
    // Validate inputs
    if config.name.contains('\t') {
        bail!("NAME cannot contain a tab character");
    }
    if config.ident_min != -1.0 && (config.ident_min < 0.0 || config.ident_min > 1.0) {
        bail!("ident_min must be between 0 and 1");
    }
    if config.coverage_min < 0.0 || config.coverage_min > 1.0 {
        bail!("coverage_min must be between 0 and 1");
    }
    if config.report_common && config.organism.is_empty() {
        bail!("--report_common requires --organism");
    }
    if config.report_common && !config.plus {
        bail!("--report_common requires --plus");
    }
    let mut gff_type = GffType::name2type(&config.annotation_format)?;
    if config.database.as_os_str().is_empty() {
        bail!("Database directory must be specified");
    }

    let database = database_version(&config.database)?;
    eprintln!("Database directory: {}", database.directory.display());
    eprintln!("Database version: {}", database.data_version);
    if config.protein.is_none() && config.nucleotide.is_none() {
        bail!("Parameter --protein or --nucleotide must be present");
    }
    if config.protein.is_none() {
        if config.gff.is_some() {
            bail!("Parameter --gff is redundant");
        }
    } else if config.nucleotide.is_some() && config.gff.is_none() {
        bail!("If parameters --protein and --nucleotide are present then parameter --gff must be present");
    }

    if config.protein.is_none() && config.protein_output.is_some() {
        bail!("Parameter --protein must be present for --protein_output");
    }
    if config.nucleotide.is_none() && config.nucleotide_output.is_some() {
        bail!("Parameter --nucleotide must be present for --nucleotide_output");
    }
    if config.nucleotide.is_none() && config.nucleotide_flank5_output.is_some() {
        bail!("Parameter --nucleotide must be present for --nucleotide_flank5_output");
    }
    if config.nucleotide_flank5_output.is_none() && config.nucleotide_flank5_size > 0 {
        bail!("Parameter --nucleotide_flank5_output must be present for --nucleotide_flank5_size");
    }
    if config.nucleotide_flank5_output.is_some() && config.nucleotide_flank5_size == 0 {
        bail!("Parameter --nucleotide_flank5_size must be present with a positive value for --nucleotide_flank5_output");
    }
    let mut search_mode = String::new();
    let mut includes = Vec::new();
    if config.protein.is_none() {
        search_mode.push_str("translated nucleotide");
    } else if config.nucleotide.is_none() {
        search_mode.push_str("protein-only");
        includes.push("-n/--nucleotide and -g/--gff options to add translated searches");
    } else {
        search_mode.push_str("combined translated and protein");
    }
    if config.organism.is_empty() {
        includes.push("-O/--organism option to add mutation searches and suppress common proteins");
    } else {
        search_mode.push_str(" and mutation");
    }
    eprintln!("AMRFinder {search_mode} search");
    for include in includes {
        eprintln!("  - include {include}");
    }
    // Create temp directory
    let tmp_dir = tempfile::tempdir()?;
    let tmp = tmp_dir.path();
    let uncompress = |path: &Path, suffix: &str| -> Result<PathBuf> {
        if path.extension().and_then(|ext| ext.to_str()) != Some("gz") {
            return Ok(path.to_path_buf());
        }
        let out = tmp.join(suffix);
        if out == path {
            bail!("Input and output file names are the same");
        }
        let output = Command::new("gunzip").arg("-c").arg(path).output()?;
        if !output.status.success() {
            bail!(
                "gunzip failed: {}",
                String::from_utf8_lossy(&output.stderr).trim()
            );
        }
        fs::write(&out, output.stdout)?;
        Ok(out)
    };
    let protein = config
        .protein
        .as_deref()
        .map(|path| uncompress(path, "prot_flat"))
        .transpose()?;
    let nucleotide = config
        .nucleotide
        .as_deref()
        .map(|path| uncompress(path, "dna_flat"))
        .transpose()?;
    let gff = config
        .gff
        .as_deref()
        .map(|path| uncompress(path, "gff_flat"))
        .transpose()?;

    let mut empty_files = Vec::new();
    if let (Some(original), Some(flat)) = (config.protein.as_deref(), protein.as_deref()) {
        if fs::metadata(flat)
            .map(|metadata| metadata.len() == 0)
            .unwrap_or(false)
        {
            empty_files.push(original);
        }
    }
    if let (Some(original), Some(flat)) = (config.nucleotide.as_deref(), nucleotide.as_deref()) {
        if fs::metadata(flat)
            .map(|metadata| metadata.len() == 0)
            .unwrap_or(false)
        {
            empty_files.push(original);
        }
    }
    if let (Some(original), Some(flat)) = (config.gff.as_deref(), gff.as_deref()) {
        if fs::metadata(flat)
            .map(|metadata| metadata.len() == 0)
            .unwrap_or(false)
        {
            empty_files.push(original);
        }
    }
    for empty_file in empty_files {
        eprintln!("WARNING: Empty file: {}", empty_file.display());
    }

    let db = database.directory.to_str().unwrap().to_string();
    let organism1 = if config.organism.is_empty() {
        String::new()
    } else {
        let mut organism1 = config.organism.replace(' ', "_");
        if config.gpipe_org {
            let tab = TextTable::from_file(&format!("{db}/taxgroup.tsv"))?;
            tab.qc()?;
            let taxgroup_col = tab.col2num("taxgroup")?;
            let gpipe_taxgroup_col = tab.col2num("gpipe_taxgroup")?;
            let mut found = false;
            for row in &tab.rows {
                let gpipe_org_vec = row[gpipe_taxgroup_col]
                    .split(',')
                    .map(str::trim)
                    .filter(|s| !s.is_empty());
                if gpipe_org_vec.into_iter().any(|s| s == organism1) {
                    organism1 = row[taxgroup_col].clone();
                    found = true;
                    break;
                }
            }
            if !found {
                eprintln!("WARNING: Non-existant GPipe taxgroup: {}", config.organism);
                organism1.clear();
            }
        }
        if !organism1.is_empty() {
            let organisms = db2organisms(Path::new(&db))?;
            if !organisms.contains(&organism1) {
                bail!("Possible organisms: {}", organisms.join(", "));
            }
        }
        organism1
    };

    if config.pgap {
        match gff_type {
            GffType::Genbank => gff_type = GffType::Pgap,
            GffType::Pgap => {}
            _ => bail!(
                "--pgap conflicts with GFF type \"{}\"",
                config.annotation_format
            ),
        }
    }
    let mut lcl = false;
    if matches!(gff_type, GffType::Pgap) {
        if let Some(dna_path) = nucleotide.as_deref() {
            let dna_text = fs::read_to_string(dna_path)?;
            for line in dna_text.lines() {
                if line.starts_with('>') {
                    lcl = line.starts_with(">lcl|");
                    break;
                }
            }
        }
    }

    let blastn = config.nucleotide.is_some()
        && !organism1.is_empty()
        && Path::new(&format!("{}/AMR_DNA-{}.fa", db, organism1)).exists();
    let stxtyper = blastn && organism1 == "Escherichia" && config.plus;
    if blastn {
        let dna_db = format!("{}/AMR_DNA-{}.fa", db, organism1);
        let db_test1 = format!("{dna_db}.ndb");
        let db_test2 = format!("{dna_db}.nin");
        if !Path::new(&db_test1).exists() && !Path::new(&db_test2).exists() {
            bail!(
                "The BLAST database for AMR_DNA-{}.fa was not found.\nUse amrfinder -u or amrfinder --force_update to download and prepare database for AMRFinderPlus",
                organism1
            );
        }
    }
    if stxtyper {
        let Some(stxtyper_path) = find_stxtyper(Path::new(&db)) else {
            bail!("stxtyper is required for Escherichia nucleotide plus reports");
        };
        let output = Command::new(&stxtyper_path).arg("-v").output()?;
        if !output.status.success() {
            bail!(
                "stxtyper failed: {}",
                String::from_utf8_lossy(&output.stderr).trim()
            );
        }
        let version = String::from_utf8_lossy(&output.stdout)
            .lines()
            .next()
            .map(str::to_string)
            .unwrap_or_default();
        if version.is_empty() {
            bail!("Cannot get the version of StxTyper");
        }
        if version != STXTYPER_VERSION {
            bail!(
                "AMRFinder invokes StxTyper version {}. Expected StxTyper version is {}",
                version,
                STXTYPER_VERSION
            );
        }
    }

    let mut amr_report_blastp = String::new();
    let mut amr_report_blastx = String::new();

    // Process protein input
    if let Some(ref prot_path) = protein {
        if !prot_path.exists() {
            bail!(
                "Protein FASTA file is empty or missing: {}",
                prot_path.display()
            );
        }
        let protein_empty = fs::metadata(prot_path)?.len() == 0;

        let blastp_out = tmp.join("blastp");
        let hmmsearch_out = tmp.join("hmmsearch");
        let dom_out = tmp.join("dom");
        let mut gff_prot_match = None;
        let mut gff_dna_match = None;
        if protein_empty {
            File::create(&blastp_out)?;
            File::create(&hmmsearch_out)?;
            File::create(&dom_out)?;
        } else {
            // fasta_check
            let mut prot1 = prot_path.to_path_buf();
            let (n_prot, _prot_len_max, prot_len_total) = match fasta_check(
                prot_path, true, 20, None, None,
            ) {
                Ok(counts) => counts,
                Err(err) => {
                    let err_s = err.to_string();
                    let fixable = if err_s.contains("Hyphen in the sequence") {
                        eprintln!(
                            "WARNING: Ignoring dash '-' characters in the sequences of the protein file {}",
                            prot_path.display()
                        );
                        true
                    } else if err_s.contains("Too many ambiguities") {
                        eprintln!(
                            "WARNING: Removing sequences with >= 20 Xs from the protein file {}",
                            prot_path.display()
                        );
                        true
                    } else {
                        false
                    };
                    if !fixable {
                        return Err(err);
                    }
                    prot1 = tmp.join("prot");
                    fasta_check(prot_path, true, 20, None, Some(&prot1))?
                }
            };

            if let Some(ref gff_path) = gff {
                if !config.parm.split_whitespace().any(|arg| arg == "-bed") {
                    if gff_type == GffType::Pgap
                        && crate::gff_check::body(
                            &gff_path.to_string_lossy(),
                            GffType::Standard,
                            &prot1,
                            nucleotide.as_deref().unwrap_or(Path::new("")),
                            Path::new(""),
                            Path::new(""),
                            false,
                            false,
                            &mut Vec::new(),
                        )
                        .is_ok()
                    {
                        gff_type = GffType::Standard;
                    }

                    let mut gff_prot_match_path = PathBuf::new();
                    match gff_type {
                        GffType::Genbank => {
                            let mut has_locus_tag = false;
                            for line in fs::read_to_string(&prot1)?.lines() {
                                if line.starts_with('>') {
                                    has_locus_tag = line.contains(crate::gff_check::LOCUS_TAG_S);
                                    break;
                                }
                            }
                            if has_locus_tag {
                                gff_prot_match_path = tmp.join("prot_match");
                            }
                        }
                        GffType::Microscope | GffType::Prodigal => {
                            gff_prot_match_path = tmp.join("prot_match");
                        }
                        _ => {}
                    }
                    if nucleotide.is_some() && gff_type == GffType::Pseudomonasdb {
                        gff_dna_match = Some(tmp.join("dna_match"));
                    }
                    if !gff_prot_match_path.as_os_str().is_empty() {
                        gff_prot_match = Some(gff_prot_match_path);
                    }

                    let gff_check_result = crate::gff_check::body(
                        &gff_path.to_string_lossy(),
                        gff_type,
                        &prot1,
                        nucleotide.as_deref().unwrap_or(Path::new("")),
                        gff_prot_match.as_deref().unwrap_or(Path::new("")),
                        gff_dna_match.as_deref().unwrap_or(Path::new("")),
                        lcl,
                        false,
                        &mut Vec::new(),
                    );
                    if let Err(err) = gff_check_result {
                        bail!("GFF file mismatch.\n{}", err);
                    }
                }
            }

            // Run blastp
            let blastp_prog = config.find_prog("blastp");
            let blastp_threads = std::cmp::min(
                config.threads,
                std::cmp::min(n_prot, prot_len_total / 10000),
            );

            let mut blastp_args = vec![
                "-query".to_string(),
                prot1.to_str().unwrap().to_string(),
                "-db".to_string(),
                format!("{}/AMRProt.fa", db),
                "-comp_based_stats".to_string(),
                "0".to_string(),
                "-seg".to_string(),
                "no".to_string(),
                "-max_target_seqs".to_string(),
                "10000".to_string(),
                "-dbsize".to_string(),
                "10000".to_string(),
                "-evalue".to_string(),
                "1e-10".to_string(),
                "-word_size".to_string(),
                "5".to_string(),
                "-task".to_string(),
                "blastp-fast".to_string(),
            ];
            if blastp_threads > 1 {
                blastp_args.extend(["-num_threads".to_string(), blastp_threads.to_string()]);
            }
            blastp_args.extend([
                "-outfmt".to_string(),
                "6 sseqid qseqid sstart send slen qstart qend qlen sseq qseq".to_string(),
                "-out".to_string(),
                blastp_out.to_str().unwrap().to_string(),
            ]);

            let blastp_output = Command::new(&blastp_prog).args(&blastp_args).output()?;

            if !blastp_output.status.success() {
                let stderr = String::from_utf8_lossy(&blastp_output.stderr);
                bail!("blastp failed: {}", stderr);
            }

            let hmm_prog = config.find_prog("hmmsearch");
            let hmm_lib = format!("{}/AMR.LIB", db);
            run_hmmsearch(
                &hmm_prog,
                &hmm_lib,
                &prot1,
                n_prot,
                config.threads,
                &hmmsearch_out,
                &dom_out,
                tmp,
            )?;
        }

        amr_report_blastp = format!(
            "-blastp {} -hmmsearch {} -hmmdom {}",
            blastp_out.to_str().unwrap(),
            hmmsearch_out.to_str().unwrap(),
            dom_out.to_str().unwrap(),
        );

        if let Some(ref gff_path) = gff {
            let gff_type_name = match gff_type {
                GffType::Bakta => "bakta",
                GffType::Genbank => "genbank",
                GffType::Microscope => "microscope",
                GffType::Patric => "patric",
                GffType::Pgap => "pgap",
                GffType::Prodigal => "prodigal",
                GffType::Prokka => "prokka",
                GffType::Pseudomonasdb => "pseudomonasdb",
                GffType::Rast => "rast",
                GffType::Standard => "standard",
            };
            amr_report_blastp += &format!(
                " -gff {} -gfftype {}",
                gff.as_deref()
                    .unwrap_or(gff_path.as_path())
                    .to_str()
                    .unwrap(),
                gff_type_name,
            );
            if let Some(path) = gff_prot_match {
                amr_report_blastp += &format!(" -gff_prot_match {}", path.to_str().unwrap());
            }
            if let Some(path) = gff_dna_match {
                amr_report_blastp += &format!(" -gff_dna_match {}", path.to_str().unwrap());
            }
            if lcl {
                amr_report_blastp.push_str(" -lcl");
            }
        }
    } else {
        // Create empty files
        File::create(tmp.join("blastp"))?;
        File::create(tmp.join("hmmsearch"))?;
        File::create(tmp.join("dom"))?;
    }

    // Process nucleotide input
    if let Some(ref dna_path) = nucleotide {
        if !dna_path.exists() {
            bail!(
                "Nucleotide FASTA file is empty or missing: {}",
                dna_path.display()
            );
        }

        // fasta_check for DNA
        let dna_len_path = tmp.join("len");
        let blastx_out = tmp.join("blastx");
        let nucleotide_empty = fs::metadata(dna_path)?.len() == 0;
        if nucleotide_empty {
            File::create(&blastx_out)?;
            File::create(&dna_len_path)?;
            if blastn {
                File::create(tmp.join("blastn"))?;
            }
        } else {
            let (n_dna, dna_len_max, dna_len_total) =
                fasta_check(dna_path, false, 0, Some(&dna_len_path), None)?;

            // Choose blastx or tblastn based on sequence length
            let use_tblastn = dna_len_max > 100_000;
            let blast_prog_name = if use_tblastn { "tblastn" } else { "blastx" };
            let blast_prog = config.find_prog(blast_prog_name);

            if use_tblastn {
                // tblastn: protein query (AMRProt.fa) vs translated nucleotide subject.
                // `-subject` mode does not parallelize via -num_threads, so the C++
                // amrfinder splits the query across processes instead (amrfinder.cpp
                // :1283-1300). Replicate that here.
                run_tblastn(
                    &blast_prog,
                    &format!("{}/AMRProt.fa", db),
                    dna_path,
                    config.translation_table,
                    config.threads,
                    &blastx_out,
                    tmp,
                )?;
            } else {
                // blastx: nucleotide query vs protein DB
                let blast_threads =
                    std::cmp::min(config.threads, std::cmp::min(n_dna, dna_len_total / 10002));
                let mut blast_args = vec![
                    "-query".to_string(),
                    dna_path.to_str().unwrap().to_string(),
                    "-db".to_string(),
                    format!("{}/AMRProt.fa", db),
                    "-comp_based_stats".to_string(),
                    "0".to_string(),
                    "-seg".to_string(),
                    "no".to_string(),
                    "-max_target_seqs".to_string(),
                    "10000".to_string(),
                    "-dbsize".to_string(),
                    "10000".to_string(),
                    "-evalue".to_string(),
                    "1e-10".to_string(),
                    "-word_size".to_string(),
                    "5".to_string(),
                    "-query_gencode".to_string(),
                    config.translation_table.to_string(),
                ];
                if blast_threads > 1 {
                    blast_args.extend(["-num_threads".to_string(), blast_threads.to_string()]);
                }
                blast_args.extend([
                    "-outfmt".to_string(),
                    "6 sseqid qseqid sstart send slen qstart qend qlen sseq qseq".to_string(),
                    "-out".to_string(),
                    blastx_out.to_str().unwrap().to_string(),
                ]);
                let blast_status = Command::new(&blast_prog)
                    .args(&blast_args)
                    .stderr(std::process::Stdio::null())
                    .stdout(std::process::Stdio::null())
                    .status()?;

                if !blast_status.success() {
                    bail!("{} failed", blast_prog_name);
                }
            }

            if !organism1.is_empty() {
                let tab = TextTable::from_file(&format!("{}/AMRProt-susceptible.tsv", db))?;
                tab.qc()?;
                let taxgroup_col = tab.col2num("taxgroup")?;
                let mut found = false;
                for row in &tab.rows {
                    ensure!(!row[taxgroup_col].is_empty(), "QC assertion failed");
                    if row[taxgroup_col] == organism1 {
                        found = true;
                        break;
                    }
                }
                if found {
                    let susceptible_fasta = format!("{}/AMRProt-susceptible.fa", db);
                    let blastx_slow_out = tmp.join("blastx-slow");
                    run_tblastn_search(
                        TblastnSearch::Slow,
                        &config.find_prog("tblastn"),
                        &susceptible_fasta,
                        dna_path,
                        config.translation_table,
                        1,
                        &blastx_slow_out,
                        tmp,
                    )?;
                    let mut blastx_file = OpenOptions::new().append(true).open(&blastx_out)?;
                    std::io::copy(&mut File::open(&blastx_slow_out)?, &mut blastx_file)?;
                }
            }
        }

        amr_report_blastx = format!(
            "-blastx {} -dna_len {}",
            blastx_out.to_str().unwrap(),
            dna_len_path.to_str().unwrap()
        );

        // DNA mutation search (blastn)
        let blastn_handle = if !nucleotide_empty && blastn {
            let dna_db = format!("{}/AMR_DNA-{}.fa", db, organism1);
            let blastn_out = tmp.join("blastn");
            let blastn_prog = config.find_prog("blastn");
            let dna_str = dna_path.to_str().unwrap().to_string();
            let blastn_out_str = blastn_out.to_str().unwrap().to_string();
            // Use reverse format (sseqid qseqid) so reference is in qseqid
            Some(std::thread::spawn(move || {
                Command::new(&blastn_prog)
                    .args([
                        "-query",
                        &dna_str,
                        "-db",
                        &dna_db,
                        "-evalue",
                        "1e-20",
                        "-dust",
                        "no",
                        "-max_target_seqs",
                        "10000",
                        "-outfmt",
                        "6 sseqid qseqid sstart send slen qstart qend qlen sseq qseq",
                        "-out",
                        &blastn_out_str,
                    ])
                    .stderr(std::process::Stdio::null())
                    .stdout(std::process::Stdio::null())
                    .status()
            }))
        } else {
            None
        };

        // Wait for blastn if it was started
        if let Some(handle) = blastn_handle {
            let result = handle
                .join()
                .map_err(|_| anyhow::anyhow!("blastn thread panicked"))??;
            if !result.success() {
                bail!("blastn failed");
            }
        }
    }

    // Run Rust amr_report + dna_mutation
    let (raw_result, fasta_source) = run_rust_amr_report(
        config,
        tmp,
        &db,
        &organism1,
        &amr_report_blastp,
        &amr_report_blastx,
        nucleotide.as_deref(),
        blastn,
    )?;

    // Post-processing: sort by sort columns
    let result = sort_tsv_output(&raw_result, config, true)?;

    // Handle output
    if let Some(ref output_path) = config.output {
        fs::write(output_path, &result)?;
    }
    write_report_fasta_outputs(
        &fasta_source,
        config,
        tmp,
        protein.as_deref(),
        nucleotide.as_deref(),
    )?;

    Ok(result)
}

/// Run the Rust amr_report implementation
fn run_rust_amr_report(
    config: &PipelineConfig,
    tmp: &Path,
    db: &str,
    organism1: &str,
    amr_report_blastp: &str,
    amr_report_blastx: &str,
    nucleotide_flat: Option<&Path>,
    blastn: bool,
) -> Result<(String, String)> {
    use crate::amr_report::{run_amr_report, AmrReportConfig};

    // Parse blastp/hmmsearch/dom file paths from the accumulated args
    let mut blastp_path = None;
    let mut hmmsearch_path = None;
    let mut hmmdom_path = None;
    let mut gff_path = None;
    let mut gff_type = "genbank";
    let mut gff_lcl = false;

    let parts: Vec<&str> = amr_report_blastp.split_whitespace().collect();
    let mut i = 0;
    while i < parts.len() {
        match parts[i] {
            "-blastp" => {
                i += 1;
                if i < parts.len() {
                    blastp_path = Some(PathBuf::from(parts[i]));
                }
            }
            "-hmmsearch" => {
                i += 1;
                if i < parts.len() {
                    hmmsearch_path = Some(PathBuf::from(parts[i]));
                }
            }
            "-hmmdom" => {
                i += 1;
                if i < parts.len() {
                    hmmdom_path = Some(PathBuf::from(parts[i]));
                }
            }
            "-gff" => {
                i += 1;
                if i < parts.len() {
                    gff_path = Some(PathBuf::from(parts[i]));
                }
            }
            "-gfftype" => {
                i += 1;
                if i < parts.len() {
                    gff_type = parts[i];
                }
            }
            "-lcl" => {
                gff_lcl = true;
            }
            _ => {}
        }
        i += 1;
    }

    if gff_lcl {
        if gff_type != "pgap" {
            bail!("-lcl requires PGAP GFF type");
        }
        if let Some(path) = gff_path.as_deref() {
            let mut normalized = String::new();
            for line in fs::read_to_string(path)?.lines() {
                if line.is_empty() || line.starts_with('#') {
                    normalized.push_str(line);
                    normalized.push('\n');
                    continue;
                }
                if let Some((contig, rest)) = line.split_once('\t') {
                    if contig.contains(':') {
                        bail!(
                            "Accession \"{}\" cannot have \"gnl|\" and \"lcl|\" at the same time",
                            contig
                        );
                    }
                    if contig.starts_with("lcl|") {
                        normalized.push_str(line);
                    } else {
                        normalized.push_str("lcl|");
                        normalized.push_str(contig);
                        normalized.push('\t');
                        normalized.push_str(rest);
                    }
                    normalized.push('\n');
                } else {
                    normalized.push_str(line);
                    normalized.push('\n');
                }
            }
            let normalized_path = tmp.join("gff_lcl");
            fs::write(&normalized_path, normalized)?;
            gff_path = Some(normalized_path);
        }
    }

    let blastx_args = parse_report_blastx_args(amr_report_blastx);
    let blastx_path = blastx_args.blastx_path;
    let dna_len_path = blastx_args.dna_len_path;
    let mut nosame = false;
    let mut noblast = false;
    let mut nohmm = false;
    let mut retain_blasts = false;
    let mut skip_hmm_check = false;
    let mut print_node = config.print_node;
    let mut print_node_raw = false;
    for arg in config.parm.split_whitespace() {
        match arg {
            "-nosame" => nosame = true,
            "-noblast" => noblast = true,
            "-nohmm" => nohmm = true,
            "-retain_blasts" => retain_blasts = true,
            "-skip_hmm_check" => skip_hmm_check = true,
            "-print_node" => print_node = true,
            "-print_node_raw" => print_node_raw = true,
            _ => {}
        }
    }
    if noblast {
        blastp_path = None;
    }

    let cds_exist = gff_path.is_some()
        || (!noblast && blastx_path.is_some())
        || (config.nucleotide.is_some() && !organism1.is_empty());

    let fam_file = PathBuf::from(format!("{}/fam.tsv", db));
    let mutation_file = PathBuf::from(format!("{}/AMRProt-mutation.tsv", db));
    let susceptible_file = PathBuf::from(format!("{}/AMRProt-susceptible.tsv", db));
    let suppress_file = if !organism1.is_empty() && !config.report_common {
        let suppress_source = PathBuf::from(format!("{}/AMRProt-suppress.tsv", db));
        let suppress_filtered = tmp.join("suppress_prot");
        let mut suppress_out = File::create(&suppress_filtered)?;
        for line in fs::read_to_string(&suppress_source)?.lines() {
            if line.starts_with('#') {
                continue;
            }
            let mut fields = line.split_whitespace();
            let Some(org) = fields.next() else {
                continue;
            };
            let Some(accver) = fields.next() else {
                bail!(
                    "Missing suppressed protein accession in {}",
                    suppress_source.display()
                );
            };
            if org == organism1 {
                writeln!(suppress_out, "{accver}")?;
            }
        }
        Some(suppress_filtered)
    } else {
        None
    };

    let report_config = AmrReportConfig {
        fam_file: &fam_file,
        blastp_file: blastp_path.as_deref(),
        blastx_file: if noblast {
            None
        } else {
            blastx_path.as_deref()
        },
        dna_len_file: dna_len_path.as_deref(),
        hmmsearch_file: hmmsearch_path.as_deref(),
        hmmdom_file: hmmdom_path.as_deref(),
        gff_file: gff_path.as_deref(),
        gff_type,
        organism: organism1,
        mutation_file: Some(&mutation_file),
        susceptible_file: Some(&susceptible_file),
        suppress_file: suppress_file.as_deref(),
        coverage_min: config.coverage_min,
        ident_min: config.ident_min,
        print_node,
        print_node_raw,
        mutation_all: config.mutation_all.as_deref(),
        name: &config.name,
        non_reportable: false,
        report_core_only: !config.plus,
        report_all_equal: config.report_all_equal,
        cds_exist,
        nosame,
        noblast,
        nohmm,
        retain_blasts,
        skip_hmm_check,
    };

    let mut output = Vec::new();
    run_amr_report(&report_config, &mut output)?;

    // Add dna_mutation results if applicable
    let mut raw_result = String::from_utf8(output)?;

    // DNA mutation detection using Rust implementation
    if blastn {
        let blastn_file = tmp.join("blastn");
        let dna_tsv_path = PathBuf::from(format!("{}/AMR_DNA-{}.tsv", db, organism1));
        if blastn_file.exists() && dna_tsv_path.exists() {
            let mut dna_output = Vec::new();
            let mutation_all_file = if config.mutation_all.is_some() {
                Some(tempfile::NamedTempFile::new()?)
            } else {
                None
            };
            crate::dna_mutation::body(
                &blastn_file,
                &dna_tsv_path,
                organism1,
                mutation_all_file.as_ref().map(|file| file.path()),
                &config.name,
                config.print_node,
                &mut dna_output,
            )?;
            let snp_output = String::from_utf8_lossy(&dna_output);
            for (i, line) in snp_output.lines().enumerate() {
                if i > 0 && !line.is_empty() {
                    raw_result.push_str(line);
                    raw_result.push('\n');
                }
            }
            if let (Some(path), Some(mutation_all_file)) = (&config.mutation_all, mutation_all_file)
            {
                let has_existing_rows = path.exists() && fs::metadata(path)?.len() > 0;
                let mut file = OpenOptions::new().create(true).append(true).open(path)?;
                let mutation_all_output = fs::read_to_string(mutation_all_file.path())?;
                for (i, line) in mutation_all_output.lines().enumerate() {
                    if line.is_empty() || (has_existing_rows && i == 0) {
                        continue;
                    }
                    writeln!(file, "{line}")?;
                }
            }
        }
    }

    if config.nucleotide.is_some()
        && organism1 == "Escherichia"
        && config.plus
        && Path::new(&format!("{}/AMR_DNA-{}.fa", db, organism1)).exists()
    {
        append_stxtyper_rows_if_available(config, tmp, db, &mut raw_result, nucleotide_flat)?;
    }

    let fasta_source = raw_result.clone();
    {
        let mut amr_tab = text_table_from_tsv(&raw_result)?;
        amr_tab_disruptions(&mut amr_tab, db, nucleotide_flat, config.translation_table)?;
        raw_result = text_table_to_tsv(&amr_tab)?;
    }

    if let Some(path) = &config.mutation_all {
        if path.exists() {
            let mutation_all = fs::read_to_string(path)?;
            let mutation_all = if config.nucleotide.is_some() {
                let mut mutation_all_tab = text_table_from_tsv(&mutation_all)?;
                amr_tab_disruptions(
                    &mut mutation_all_tab,
                    db,
                    nucleotide_flat,
                    config.translation_table,
                )?;
                text_table_to_tsv(&mutation_all_tab)?
            } else {
                mutation_all
            };
            let sorted = sort_tsv_output(&mutation_all, config, false)?;
            fs::write(path, sorted)?;
        }
    }

    Ok((raw_result, fasta_source))
}

fn amr_tab_disruptions(
    amr_tab: &mut TextTable,
    db: &str,
    dna_flat: Option<&Path>,
    gencode: u32,
) -> Result<()> {
    let Some(dna_flat) = dna_flat else {
        return Ok(());
    };
    let contig_col = amr_tab.col2num(crate::columns::CONTIG_COL_NAME)?;
    let prot_col = amr_tab.col2num(crate::columns::CLOSEST_REF_ACCESSION_COL_NAME)?;
    let genesymbol_col = amr_tab.col2num(crate::columns::GENESYMBOL_COL_NAME)?;

    let mut disr_raw = tempfile::NamedTempFile::new()?;
    for row in &amr_tab.rows {
        let contig = row_field(row, contig_col);
        let prot = row_field(row, prot_col);
        let genesymbol = row_field(row, genesymbol_col);
        let Some((_, disruption)) = genesymbol.split_once(crate::columns::DISRUPTION_DELIM) else {
            continue;
        };
        if contig == crate::columns::NA || prot == crate::columns::NA {
            continue;
        }
        writeln!(disr_raw, "{contig}\t{prot}\t{disruption}")?;
    }

    if disr_raw.as_file().metadata()?.len() == 0 {
        return Ok(());
    }

    let mut disr_out = Vec::new();
    crate::disruption2genesymbol::body(
        dna_flat,
        &Path::new(db).join("AMRProt-susceptible.fa"),
        disr_raw.path(),
        gencode.try_into()?,
        1,
        &mut disr_out,
    )?;
    let disr_out = String::from_utf8(disr_out)?;
    let mut disr_rows = Vec::new();
    for line in disr_out.lines().filter(|line| !line.is_empty()) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            bail!("Cannot read disruption2genesymbol output: {line}");
        }
        disr_rows.push((
            fields[0].to_string(),
            fields[1].to_string(),
            fields[2].to_string(),
            fields[3].to_string(),
        ));
    }

    for row in &mut amr_tab.rows {
        let contig = row_field(row, contig_col);
        let prot = row_field(row, prot_col);
        let genesymbol = row_field(row, genesymbol_col);
        let Some((prefix, disruption)) = genesymbol.split_once(crate::columns::DISRUPTION_DELIM)
        else {
            continue;
        };

        let disr = disr_rows
            .iter()
            .find(|(disr_contig, disr_prot, _, disr_raw)| {
                disr_contig == contig && disr_prot == prot && disr_raw == disruption
            })
            .map(|(_, _, disr, _)| disr.clone());
        let Some(disr) = disr else {
            bail!("Disruption is not replaced by gene symbol:\n{contig} {prot} {disruption}");
        };
        row[genesymbol_col] = format!("{prefix}_{disr}");
    }

    Ok(())
}

fn text_table_from_tsv(tsv: &str) -> Result<TextTable> {
    let mut lines = tsv.lines().filter(|line| !line.is_empty());
    let Some(header_line) = lines.next() else {
        bail!("Cannot read the table header");
    };
    let pound = header_line.starts_with('#');
    let header_line = header_line.trim_start_matches('#');
    let header = header_line.split('\t').map(Header::new).collect();
    let rows = lines
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect();
    let table = TextTable {
        name: String::new(),
        pound,
        save_header: true,
        header,
        rows,
    };
    table.qc()?;
    Ok(table)
}

fn text_table_to_tsv(table: &TextTable) -> Result<String> {
    let mut output = Vec::new();
    table.save_text(&mut output)?;
    Ok(String::from_utf8(output)?)
}

#[allow(dead_code)]
fn prepare_fasta_extract(
    amr_tab: &TextTable,
    columns: &[&str],
    save_header: bool,
) -> Result<TextTable> {
    amr_tab.qc()?;
    let col_nums = amr_tab.columns2nums(columns)?;
    let header = col_nums
        .iter()
        .map(|&col| amr_tab.header[col].clone())
        .collect();
    let mut rows: Vec<Vec<String>> = amr_tab
        .rows
        .iter()
        .map(|row| col_nums.iter().map(|&col| row[col].clone()).collect())
        .filter(|row: &Vec<String>| row.first().map(String::as_str) != Some(crate::columns::NA))
        .collect();
    rows.sort();
    rows.dedup();
    let table = TextTable {
        name: String::new(),
        pound: false,
        save_header,
        header,
        rows,
    };
    table.qc()?;
    Ok(table)
}

fn write_report_fasta_outputs(
    tsv: &str,
    config: &PipelineConfig,
    tmp: &Path,
    protein_flat: Option<&Path>,
    nucleotide_flat: Option<&Path>,
) -> Result<()> {
    if config.protein_output.is_none()
        && config.nucleotide_output.is_none()
        && config.nucleotide_flank5_output.is_none()
    {
        return Ok(());
    }

    let amr_tab = text_table_from_tsv(tsv)?;
    if let (Some(prot), Some(prot_out)) = (protein_flat, &config.protein_output) {
        let target = prepare_fasta_extract(
            &amr_tab,
            &[
                crate::columns::PROT_COL_NAME,
                crate::columns::GENESYMBOL_COL_NAME,
                crate::columns::ELEM_NAME_COL_NAME,
            ],
            false,
        )?;
        let target_path = tmp.join("prot_out");
        fs::write(&target_path, text_table_to_tsv(&target)?)?;
        let mut out = File::create(prot_out)?;
        crate::fasta_extract::body(prot, &target_path, true, false, &mut out)?;
    }
    if let (Some(dna), Some(dna_out)) = (nucleotide_flat, &config.nucleotide_output) {
        let target = prepare_fasta_extract(
            &amr_tab,
            &[
                crate::columns::CONTIG_COL_NAME,
                crate::columns::START_COL_NAME,
                crate::columns::STOP_COL_NAME,
                crate::columns::STRAND_COL_NAME,
                crate::columns::GENESYMBOL_COL_NAME,
                crate::columns::ELEM_NAME_COL_NAME,
            ],
            false,
        )?;
        let target_path = tmp.join("dna_out");
        fs::write(&target_path, text_table_to_tsv(&target)?)?;
        let mut out = File::create(dna_out)?;
        crate::fasta_extract::body(dna, &target_path, false, false, &mut out)?;
    }
    if let (Some(dna), Some(dna_flank5_out)) = (nucleotide_flat, &config.nucleotide_flank5_output) {
        let mut target = prepare_fasta_extract(
            &amr_tab,
            &[
                crate::columns::CONTIG_COL_NAME,
                crate::columns::START_COL_NAME,
                crate::columns::STOP_COL_NAME,
                crate::columns::STRAND_COL_NAME,
                crate::columns::GENESYMBOL_COL_NAME,
                crate::columns::ELEM_NAME_COL_NAME,
            ],
            true,
        )?;
        for row in &mut target.rows {
            if row[3] == "+" {
                let start = row[1].parse::<usize>()?;
                row[1] = start
                    .saturating_sub(config.nucleotide_flank5_size)
                    .max(1)
                    .to_string();
            } else {
                let stop = row[2].parse::<usize>()?;
                row[2] = (stop + config.nucleotide_flank5_size).to_string();
            }
        }
        target.save_header = false;
        target.qc()?;
        let target_path = tmp.join("dnaFlank5_out");
        fs::write(&target_path, text_table_to_tsv(&target)?)?;
        let mut out = File::create(dna_flank5_out)?;
        crate::fasta_extract::body(dna, &target_path, false, false, &mut out)?;
    }
    Ok(())
}

fn append_stxtyper_rows_if_available(
    config: &PipelineConfig,
    tmp: &Path,
    db: &str,
    raw_result: &mut String,
    nucleotide_flat: Option<&Path>,
) -> Result<bool> {
    let Some(nucleotide) = nucleotide_flat.or(config.nucleotide.as_deref()) else {
        return Ok(false);
    };
    let stxtyper = find_stxtyper(Path::new(db));
    let Some(stxtyper) = stxtyper else {
        bail!("stxtyper is required for Escherichia nucleotide plus reports");
    };

    let stxtyper_out = tmp.join("stxtyper.tsv");
    let mut command = Command::new(&stxtyper);
    command
        .arg("-n")
        .arg(nucleotide)
        .arg("--amrfinder")
        .arg("-o")
        .arg(&stxtyper_out)
        .arg("--name")
        .arg(&config.name)
        .arg("--threads")
        .arg(config.threads.to_string())
        .arg("-q");
    if config.print_node {
        command.arg("--print_node");
    }
    if !config.blast_bin.is_empty() {
        command.arg("--blast_bin").arg(&config.blast_bin);
    }

    let output = match command.output() {
        Ok(output) => output,
        Err(err) if err.kind() == ErrorKind::NotFound => {
            bail!("stxtyper is required for Escherichia nucleotide plus reports")
        }
        Err(err) => return Err(err.into()),
    };
    if !output.status.success() {
        bail!(
            "stxtyper failed: {}",
            String::from_utf8_lossy(&output.stderr).trim()
        );
    }
    if !stxtyper_out.exists() {
        bail!("stxtyper did not write {}", stxtyper_out.display());
    }
    let stxtyper_tsv = fs::read_to_string(stxtyper_out)?;
    append_report_rows(raw_result, &stxtyper_tsv);
    Ok(true)
}

fn find_stxtyper(db: &Path) -> Option<PathBuf> {
    let _ = db;
    if let Ok(exe) = std::env::current_exe() {
        if let Some(exec_dir) = exe.parent() {
            let candidate = exec_dir.join("stx").join("stxtyper");
            if candidate.is_file() && is_executable(&candidate) {
                return Some(candidate);
            }
        }
    }
    None
}

#[cfg(unix)]
fn is_executable(path: &Path) -> bool {
    use std::os::unix::fs::PermissionsExt;
    fs::metadata(path)
        .map(|metadata| metadata.permissions().mode() & 0o111 != 0)
        .unwrap_or(false)
}

#[cfg(not(unix))]
fn is_executable(path: &Path) -> bool {
    path.is_file()
}

fn append_report_rows(raw_result: &mut String, extra_tsv: &str) {
    let Some(raw_header) = raw_result.lines().next() else {
        return;
    };
    let raw_header = raw_header.trim_start_matches('#').to_string();
    let mut existing: std::collections::HashSet<String> =
        raw_result.lines().skip(1).map(str::to_string).collect();
    for line in extra_tsv.lines() {
        if line.is_empty() {
            continue;
        }
        if line.trim_start_matches('#') == raw_header {
            continue;
        }
        if existing.contains(line) {
            continue;
        }
        raw_result.push_str(line);
        raw_result.push('\n');
        existing.insert(line.to_string());
    }
}

/// Sort the TSV output to match C++ amrfinder's post-processing
/// Sort columns: Contig id, Start, Stop, Strand, Protein id, Element symbol
fn sort_tsv_output(tsv: &str, config: &PipelineConfig, deredundify: bool) -> Result<String> {
    let has_cds = config.nucleotide.is_some() || config.gff.is_some();
    let mut table = text_table_from_tsv(tsv)?;
    table.set_header()?;

    if deredundify && config.nucleotide.is_some() {
        let subtype_col = table.col2num(crate::columns::SUBTYPE_COL_NAME)?;
        let start_col = table.col2num(crate::columns::START_COL_NAME)?;
        let stop_col = table.col2num(crate::columns::STOP_COL_NAME)?;
        let strand_col = table.col2num(crate::columns::STRAND_COL_NAME)?;
        table.deredundify(
            &[
                crate::columns::CONTIG_COL_NAME,
                crate::columns::STRAND_COL_NAME,
                crate::columns::GENESYMBOL_COL_NAME,
            ],
            |better, worse| {
                amr_tab_equiv_better(better, worse, subtype_col, start_col, stop_col, strand_col)
            },
        )?;
    }

    if has_cds {
        table.sort(&[
            crate::columns::CONTIG_COL_NAME,
            crate::columns::START_COL_NAME,
            crate::columns::STOP_COL_NAME,
            crate::columns::STRAND_COL_NAME,
            crate::columns::PROT_COL_NAME,
            crate::columns::GENESYMBOL_COL_NAME,
        ])?;
    } else {
        table.sort(&[
            crate::columns::PROT_COL_NAME,
            crate::columns::GENESYMBOL_COL_NAME,
        ])?;
    }
    table.rows.dedup();
    table.qc()?;

    text_table_to_tsv(&table)
}

fn amr_tab_equiv_better(
    row_better: &[String],
    row_worse: &[String],
    subtype_col: usize,
    start_col: usize,
    stop_col: usize,
    strand_col: usize,
) -> i32 {
    if row_field(row_better, subtype_col) == "POINT"
        && row_field(row_worse, subtype_col) == "POINT_DISRUPT"
        && row_field(row_better, strand_col) == row_field(row_worse, strand_col)
        && ((row_field(row_better, strand_col) == "+"
            && row_field(row_better, start_col) == row_field(row_worse, start_col))
            || (row_field(row_better, strand_col) == "-"
                && row_field(row_better, stop_col) == row_field(row_worse, stop_col)))
    {
        1
    } else {
        0
    }
}

fn row_field(row: &[String], col: usize) -> &str {
    row.get(col).map(String::as_str).unwrap_or("")
}
