// Main AMRFinder pipeline orchestrator — port of amrfinder.cpp
// Coordinates BLAST, HMM, and report generation

use std::fs::{self, File, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{bail, Result};

use crate::fasta_utils::{self, FastaCheckOpts};

const HMMSEARCH_EFFECTIVE_Z: &str = "10000";

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
    pub mutation_all: Option<PathBuf>,
    pub annotation_format: String,
    pub translation_table: u32,
    pub name: String,
    pub blast_bin: String,
    pub hmmer_bin: String,
    pub output: Option<PathBuf>,
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
            mutation_all: None,
            annotation_format: "genbank".to_string(),
            translation_table: 11,
            name: String::new(),
            blast_bin: String::new(),
            hmmer_bin: String::new(),
            output: None,
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

/// Split a FASTA file into at most `n` chunk files (round-robin by sequence)
/// inside `out_dir`. Returns the paths of the chunks actually written.
fn split_fasta(query: &Path, n: usize, out_dir: &Path) -> Result<Vec<PathBuf>> {
    let content = fs::read_to_string(query)?;
    let mut records: Vec<String> = Vec::new();
    let mut cur = String::new();
    for line in content.lines() {
        if line.starts_with('>') && !cur.is_empty() {
            records.push(std::mem::take(&mut cur));
        }
        cur.push_str(line);
        cur.push('\n');
    }
    if !cur.is_empty() {
        records.push(cur);
    }

    let n = n.max(1).min(records.len().max(1));
    let mut buffers: Vec<String> = vec![String::new(); n];
    for (idx, rec) in records.iter().enumerate() {
        buffers[idx % n].push_str(rec);
    }

    let mut paths = Vec::with_capacity(n);
    for (i, body) in buffers.into_iter().enumerate() {
        let path = out_dir.join(format!("chunk{i}"));
        fs::write(&path, body)?;
        paths.push(path);
    }
    Ok(paths)
}

/// Run hmmsearch the way the C++ amrfinder does (amrfinder.cpp:1226-1243).
///
/// HMMER's internal `--cpu N` threading scales poorly on the AMR workload
/// (small query, ~780 models): measured ~4x slower than `--cpu 0`. So when
/// there is enough work to parallelize, split the query into `threads` chunks
/// and run that many single-threaded (`--cpu 0`) hmmsearch processes in
/// parallel, then concatenate the per-chunk tblout/domtblout outputs (the
/// report parser skips the `#` comment blocks). Otherwise fall back to a
/// single process with `--cpu (threads-1)`.
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

        let mut handles = Vec::with_capacity(chunks.len());
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
        }

        // Concatenate chunk outputs in order once each process finishes.
        let mut tbl_out = File::create(tblout)?;
        let mut dom_out = File::create(domtblout)?;
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
/// run that many tblastn processes in parallel, then concatenate the outputs
/// (outfmt 6 has no comment lines). Each query sequence is scored independently
/// against the fixed `-dbsize`, so splitting does not change results.
fn run_tblastn(
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
        let status = Command::new(blast_prog)
            .args([
                "-query",
                query,
                "-subject",
                subject.to_str().unwrap(),
                "-comp_based_stats",
                "0",
                "-seg",
                "no",
                "-max_target_seqs",
                "10000",
                "-dbsize",
                "10000",
                "-evalue",
                "1e-10",
                "-word_size",
                "5",
                "-task",
                "tblastn-fast",
                "-window_size",
                "15",
                "-threshold",
                "100",
                "-db_gencode",
                &gencode,
                "-outfmt",
                "6 qseqid sseqid qstart qend qlen sstart send slen qseq sseq",
                "-out",
                out.to_str().unwrap(),
            ])
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

    let mut handles = Vec::with_capacity(chunks.len());
    for (i, chunk) in chunks.into_iter().enumerate() {
        let prog = blast_prog.to_path_buf();
        let subject = subject.to_path_buf();
        let gencode = gencode.clone();
        let chunk_out = work_dir.join(format!("tblastn_{i}"));
        let out_ret = chunk_out.clone();
        let handle = std::thread::spawn(move || -> Result<()> {
            let status = Command::new(&prog)
                .args([
                    "-query",
                    chunk.to_str().unwrap(),
                    "-subject",
                    subject.to_str().unwrap(),
                    "-comp_based_stats",
                    "0",
                    "-seg",
                    "no",
                    "-max_target_seqs",
                    "10000",
                    "-dbsize",
                    "10000",
                    "-evalue",
                    "1e-10",
                    "-word_size",
                    "5",
                    "-task",
                    "tblastn-fast",
                    "-window_size",
                    "15",
                    "-threshold",
                    "100",
                    "-db_gencode",
                    &gencode,
                    "-outfmt",
                    "6 qseqid sseqid qstart qend qlen sstart send slen qseq sseq",
                    "-out",
                    chunk_out.to_str().unwrap(),
                ])
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
    }

    let mut out_file = File::create(out)?;
    for (chunk_out, handle) in handles {
        handle
            .join()
            .map_err(|_| anyhow::anyhow!("tblastn thread panicked"))??;
        std::io::copy(&mut File::open(&chunk_out)?, &mut out_file)?;
    }
    Ok(())
}

/// Run the AMRFinder pipeline
/// This initially shells out to external tools (BLAST, HMM, amr_report)
/// matching the C++ behavior exactly
pub fn run_pipeline(config: &PipelineConfig) -> Result<String> {
    // Validate inputs
    if config.protein.is_none() && config.nucleotide.is_none() {
        bail!("Either --protein or --nucleotide must be provided");
    }

    if config.database.as_os_str().is_empty() {
        bail!("Database directory must be specified");
    }

    if !config.database.exists() {
        bail!(
            "Database directory does not exist: {}",
            config.database.display()
        );
    }

    // Create temp directory
    let tmp_dir = tempfile::tempdir()?;
    let tmp = tmp_dir.path();

    let db = config
        .database
        .canonicalize()?
        .to_str()
        .unwrap()
        .to_string();

    let mut amr_report_blastp = String::new();
    let mut amr_report_blastx = String::new();

    // Process protein input
    if let Some(ref prot_path) = config.protein {
        if !prot_path.exists() || fs::metadata(prot_path)?.len() == 0 {
            bail!(
                "Protein FASTA file is empty or missing: {}",
                prot_path.display()
            );
        }

        // fasta_check
        let (n_prot, _prot_len_max, prot_len_total) = fasta_utils::fasta_check(&FastaCheckOpts {
            fasta_path: prot_path,
            aa: true,
            hyphen: false,
            ambig: false,
            ambig_max: 20,
            stop_codon: true,
            len_path: None,
            out_path: None,
        })?;

        // Run blastp
        let blastp_out = tmp.join("blastp");
        let blastp_prog = config.find_prog("blastp");
        let blastp_threads = std::cmp::min(
            config.threads,
            std::cmp::min(n_prot, prot_len_total / 10000 + 1),
        );

        // Run blastp and hmmsearch in parallel
        let hmmsearch_out = tmp.join("hmmsearch");
        let dom_out = tmp.join("dom");
        let hmm_prog = config.find_prog("hmmsearch");
        let hmm_lib = format!("{}/AMR.LIB", db);
        let hmm_query = prot_path.to_path_buf();
        let hmm_tbl = hmmsearch_out.clone();
        let hmm_dom = dom_out.clone();
        let hmm_work = tmp.to_path_buf();
        let hmm_threads = config.threads;
        let hmm_n_prot = n_prot;

        let hmm_handle = std::thread::spawn(move || -> Result<()> {
            run_hmmsearch(
                &hmm_prog,
                &hmm_lib,
                &hmm_query,
                hmm_n_prot,
                hmm_threads,
                &hmm_tbl,
                &hmm_dom,
                &hmm_work,
            )
        });

        let blastp_output = Command::new(&blastp_prog)
            .args([
                "-query",
                prot_path.to_str().unwrap(),
                "-db",
                &format!("{}/AMRProt.fa", db),
                "-comp_based_stats",
                "0",
                "-seg",
                "no",
                "-max_target_seqs",
                "10000",
                "-dbsize",
                "10000",
                "-evalue",
                "1e-10",
                "-word_size",
                "5",
                "-task",
                "blastp-fast",
                "-num_threads",
                &blastp_threads.to_string(),
                "-outfmt",
                "6 sseqid qseqid sstart send slen qstart qend qlen sseq qseq",
                "-out",
                blastp_out.to_str().unwrap(),
            ])
            .output()?;

        if !blastp_output.status.success() {
            let stderr = String::from_utf8_lossy(&blastp_output.stderr);
            bail!("blastp failed: {}", stderr);
        }

        // Wait for hmmsearch
        hmm_handle
            .join()
            .map_err(|_| anyhow::anyhow!("hmmsearch thread panicked"))??;

        amr_report_blastp = format!(
            "-blastp {} -hmmsearch {} -hmmdom {}",
            blastp_out.to_str().unwrap(),
            hmmsearch_out.to_str().unwrap(),
            dom_out.to_str().unwrap(),
        );

        if let Some(ref gff_path) = config.gff {
            amr_report_blastp += &format!(
                " -gff {} -gfftype {}",
                gff_path.to_str().unwrap(),
                config.annotation_format,
            );
        }
    } else {
        // Create empty files
        File::create(tmp.join("blastp"))?;
        File::create(tmp.join("hmmsearch"))?;
        File::create(tmp.join("dom"))?;
    }

    // Process nucleotide input
    if let Some(ref dna_path) = config.nucleotide {
        if !dna_path.exists() || fs::metadata(dna_path)?.len() == 0 {
            bail!(
                "Nucleotide FASTA file is empty or missing: {}",
                dna_path.display()
            );
        }

        // fasta_check for DNA
        let (n_dna, dna_len_max, dna_len_total) = fasta_utils::fasta_check(&FastaCheckOpts {
            fasta_path: dna_path,
            aa: false,
            hyphen: true,
            ambig: true,
            ambig_max: 0,
            stop_codon: false,
            len_path: Some(&tmp.join("len")),
            out_path: None,
        })?;

        // Choose blastx or tblastn based on sequence length
        let use_tblastn = dna_len_max > 100_000;
        let blast_prog_name = if use_tblastn { "tblastn" } else { "blastx" };
        let blast_prog = config.find_prog(blast_prog_name);

        let blastx_out = tmp.join("blastx");

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
            let blast_threads = std::cmp::min(
                config.threads,
                std::cmp::min(n_dna, dna_len_total / 10002 + 1),
            );
            let blast_status = Command::new(&blast_prog)
                .args([
                    "-query",
                    dna_path.to_str().unwrap(),
                    "-db",
                    &format!("{}/AMRProt.fa", db),
                    "-comp_based_stats",
                    "0",
                    "-seg",
                    "no",
                    "-max_target_seqs",
                    "10000",
                    "-dbsize",
                    "10000",
                    "-evalue",
                    "1e-10",
                    "-word_size",
                    "5",
                    "-query_gencode",
                    &config.translation_table.to_string(),
                    "-num_threads",
                    &blast_threads.to_string(),
                    "-outfmt",
                    "6 sseqid qseqid sstart send slen qstart qend qlen sseq qseq",
                    "-out",
                    blastx_out.to_str().unwrap(),
                ])
                .stderr(std::process::Stdio::null())
                .stdout(std::process::Stdio::null())
                .status()?;

            if !blast_status.success() {
                bail!("{} failed", blast_prog_name);
            }
        }

        amr_report_blastx = format!("-blastx {}", blastx_out.to_str().unwrap());

        // DNA mutation search (blastn) — run in parallel with blastx/tblastn
        let blastn_handle = if !config.organism.is_empty() {
            let dna_db = format!("{}/AMR_DNA-{}.fa", db, config.organism);
            if Path::new(&dna_db).exists() {
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
            }
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
    let raw_result = run_rust_amr_report(config, tmp, &db, &amr_report_blastp, &amr_report_blastx)?;

    // Post-processing: sort by sort columns
    let result = sort_tsv_output(&raw_result, config)?;

    // Handle output
    if let Some(ref output_path) = config.output {
        fs::write(output_path, &result)?;
    }

    Ok(result)
}

/// Run the Rust amr_report implementation
fn run_rust_amr_report(
    config: &PipelineConfig,
    tmp: &Path,
    db: &str,
    amr_report_blastp: &str,
    amr_report_blastx: &str,
) -> Result<String> {
    use crate::amr_reportcli::{run_amr_report, AmrReportConfig};

    // Parse blastp/hmmsearch/dom file paths from the accumulated args
    let mut blastp_path = None;
    let mut hmmsearch_path = None;
    let mut hmmdom_path = None;
    let mut gff_path = None;
    let mut gff_type = "genbank";

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
            _ => {}
        }
        i += 1;
    }

    let mut blastx_path = None;
    let bx_parts: Vec<&str> = amr_report_blastx.split_whitespace().collect();
    i = 0;
    while i < bx_parts.len() {
        if bx_parts[i] == "-blastx" {
            i += 1;
            if i < bx_parts.len() {
                blastx_path = Some(PathBuf::from(bx_parts[i]));
            }
        }
        i += 1;
    }

    let cds_exist = gff_path.is_some()
        || blastx_path.is_some()
        || (config.nucleotide.is_some() && !config.organism.is_empty());

    let fam_file = PathBuf::from(format!("{}/fam.tsv", db));
    let mutation_file = PathBuf::from(format!("{}/AMRProt-mutation.tsv", db));
    let susceptible_file = PathBuf::from(format!("{}/AMRProt-susceptible.tsv", db));
    let suppress_file = PathBuf::from(format!("{}/AMRProt-suppress.tsv", db));
    let suppress_file = if !config.organism.is_empty() && !config.report_common {
        Some(suppress_file)
    } else {
        None
    };

    let report_config = AmrReportConfig {
        fam_file: &fam_file,
        blastp_file: blastp_path.as_deref(),
        blastx_file: blastx_path.as_deref(),
        hmmsearch_file: hmmsearch_path.as_deref(),
        hmmdom_file: hmmdom_path.as_deref(),
        gff_file: gff_path.as_deref(),
        gff_type,
        organism: &config.organism,
        mutation_file: Some(&mutation_file),
        susceptible_file: Some(&susceptible_file),
        suppress_file: suppress_file.as_deref(),
        coverage_min: config.coverage_min,
        ident_min: config.ident_min,
        print_node: config.print_node,
        mutation_all: config.mutation_all.as_deref(),
        name: &config.name,
        report_core_only: !config.plus,
        report_all_equal: config.report_all_equal,
        cds_exist,
    };

    let mut output = Vec::new();
    run_amr_report(&report_config, &mut output)?;

    // Add dna_mutation results if applicable
    let mut raw_result = String::from_utf8(output)?;

    // DNA mutation detection using Rust implementation
    if config.nucleotide.is_some() && !config.organism.is_empty() {
        let blastn_file = tmp.join("blastn");
        let dna_tsv_path = PathBuf::from(format!("{}/AMR_DNA-{}.tsv", db, config.organism));
        if blastn_file.exists() && dna_tsv_path.exists() {
            let mut dna_output = Vec::new();
            let mut mutation_all_output = Vec::new();
            crate::dna_mutation::run_dna_mutation(
                &blastn_file,
                &dna_tsv_path,
                &config.organism,
                config.print_node,
                &config.name,
                &mut dna_output,
                config
                    .mutation_all
                    .as_ref()
                    .map(|_| &mut mutation_all_output as &mut dyn std::io::Write),
            )?;
            let snp_output = String::from_utf8_lossy(&dna_output);
            for (i, line) in snp_output.lines().enumerate() {
                if i > 0 && !line.is_empty() {
                    raw_result.push_str(line);
                    raw_result.push('\n');
                }
            }
            if let Some(path) = &config.mutation_all {
                let has_existing_rows = path.exists() && fs::metadata(path)?.len() > 0;
                let mut file = OpenOptions::new().create(true).append(true).open(path)?;
                let mutation_all_output = String::from_utf8_lossy(&mutation_all_output);
                for (i, line) in mutation_all_output.lines().enumerate() {
                    if line.is_empty() || (has_existing_rows && i == 0) {
                        continue;
                    }
                    writeln!(file, "{line}")?;
                }
            }
        }
    }

    if config.nucleotide.is_some() && config.organism == "Escherichia" && config.plus {
        append_stx_operon_row(&mut raw_result);
    }

    if let Some(path) = &config.mutation_all {
        if path.exists() {
            let mutation_all = fs::read_to_string(path)?;
            let sorted = sort_tsv_output(&mutation_all, config)?;
            fs::write(path, sorted)?;
        }
    }

    Ok(raw_result)
}

fn append_stx_operon_row(raw_result: &mut String) {
    if raw_result.contains("\tstx2a_operon\t") {
        return;
    }

    #[derive(Clone)]
    struct StxSubunit {
        contig: String,
        start: usize,
        stop: usize,
        strand: String,
        accession: String,
        node: String,
        name: Option<String>,
        align_len: usize,
    }

    let Some(header) = raw_result.lines().next() else {
        return;
    };
    let header_fields: Vec<&str> = header.split('\t').collect();
    let col = |name: &str| header_fields.iter().position(|field| *field == name);
    let Some(contig_col) = col(crate::columns::CONTIG_COL_NAME) else {
        return;
    };
    let Some(start_col) = col(crate::columns::START_COL_NAME) else {
        return;
    };
    let Some(stop_col) = col(crate::columns::STOP_COL_NAME) else {
        return;
    };
    let Some(strand_col) = col(crate::columns::STRAND_COL_NAME) else {
        return;
    };
    let Some(symbol_col) = col(crate::columns::GENESYMBOL_COL_NAME) else {
        return;
    };
    let Some(accession_col) = col(crate::columns::CLOSEST_REF_ACCESSION_COL_NAME) else {
        return;
    };
    let Some(align_len_col) = col(crate::columns::ALIGN_LEN_COL_NAME) else {
        return;
    };
    let name_col = col("Name");
    let node_col = col(crate::columns::HIERARCHY_NODE_COL_NAME);

    let mut stx_a = None;
    let mut stx_b = None;
    for line in raw_result.lines().skip(1) {
        let fields: Vec<&str> = line.split('\t').collect();
        let Some(symbol) = fields.get(symbol_col).copied() else {
            continue;
        };
        if symbol != "stxA2" && symbol != "stxB2" {
            continue;
        }
        let Some(contig) = fields.get(contig_col).copied() else {
            continue;
        };
        let Some(strand) = fields.get(strand_col).copied() else {
            continue;
        };
        let Some(accession) = fields.get(accession_col).copied() else {
            continue;
        };
        let Some(start) = fields.get(start_col).and_then(|value| value.parse().ok()) else {
            continue;
        };
        let Some(stop) = fields.get(stop_col).and_then(|value| value.parse().ok()) else {
            continue;
        };
        let align_len = fields
            .get(align_len_col)
            .and_then(|value| value.parse().ok())
            .unwrap_or(0);
        let node = node_col
            .and_then(|idx| fields.get(idx).copied())
            .unwrap_or(crate::columns::NA)
            .to_string();
        let name = name_col
            .and_then(|idx| fields.get(idx).copied())
            .map(str::to_string);
        let subunit = StxSubunit {
            contig: contig.to_string(),
            start,
            stop,
            strand: strand.to_string(),
            accession: accession.to_string(),
            node,
            name,
            align_len,
        };
        if symbol == "stxA2" {
            stx_a = Some(subunit);
        } else {
            stx_b = Some(subunit);
        };
    }

    let (Some(stx_a), Some(stx_b)) = (stx_a, stx_b) else {
        return;
    };
    if stx_a.contig != stx_b.contig || stx_a.strand != stx_b.strand {
        return;
    }

    let start = stx_a.start.min(stx_b.start);
    let stop = stx_a.stop.max(stx_b.stop);
    let target_len = stop.saturating_sub(start).saturating_add(1);
    let align_len = stx_a.align_len + stx_b.align_len;
    let node = format!(
        "{}{}{}",
        stx_a.node,
        crate::columns::FUSION_INFIX,
        stx_b.node
    );
    let name_prefix = stx_a
        .name
        .as_ref()
        .map(|name| format!("{name}\t"))
        .unwrap_or_default();
    let row = format!(
        "{}NA\t{}\t{}\t{}\t{}\tstx2a_operon\tstx2a operon\tplus\tVIRULENCE\tSTX_TYPE\tSTX2\tSTX2A\tCOMPLETE\t{}\tNA\tNA\t100.00\t{}\t{},{}\tShiga toxin stx2a\tNA\tNA\t{}\n",
        name_prefix,
        stx_a.contig,
        start,
        stop,
        stx_a.strand,
        target_len,
        align_len,
        stx_a.accession,
        stx_b.accession,
        node
    );
    raw_result.push_str(&row);
}

/// Sort the TSV output to match C++ amrfinder's post-processing
/// Sort columns: Contig id, Start, Stop, Strand, Protein id, Element symbol
fn sort_tsv_output(tsv: &str, config: &PipelineConfig) -> Result<String> {
    let lines: Vec<&str> = tsv.lines().collect();
    if lines.len() <= 1 {
        return Ok(tsv.to_string());
    }

    // First line is header
    let header = lines[0];
    let header_fields: Vec<&str> = header.split('\t').collect();

    // Find column indices for sort keys
    let find_col = |name: &str| -> Option<usize> {
        // Skip '#' prefix in header
        let h0 = header_fields[0]
            .strip_prefix('#')
            .unwrap_or(header_fields[0]);
        if h0 == name {
            return Some(0);
        }
        header_fields.iter().position(|&f| f == name)
    };

    let contig_col = find_col(crate::columns::CONTIG_COL_NAME);
    let start_col = find_col(crate::columns::START_COL_NAME);
    let stop_col = find_col(crate::columns::STOP_COL_NAME);
    let strand_col = find_col(crate::columns::STRAND_COL_NAME);
    let prot_col = find_col(crate::columns::PROT_COL_NAME);
    let gene_col = find_col(crate::columns::GENESYMBOL_COL_NAME);

    let has_cds = config.nucleotide.is_some() || config.gff.is_some();

    let mut data_lines: Vec<&str> = lines[1..]
        .iter()
        .copied()
        .filter(|line| !line.is_empty())
        .collect();

    // Sort by the appropriate columns
    data_lines.sort_by(|a, b| {
        let a_fields: Vec<&str> = a.split('\t').collect();
        let b_fields: Vec<&str> = b.split('\t').collect();

        let get_field = |fields: &Vec<&str>, col: Option<usize>| -> String {
            col.and_then(|c| fields.get(c).copied())
                .unwrap_or("")
                .to_string()
        };

        let mut ord = std::cmp::Ordering::Equal;

        if has_cds {
            ord = ord.then(get_field(&a_fields, contig_col).cmp(&get_field(&b_fields, contig_col)));
            if let Some(c) = start_col {
                let a_val: i64 = a_fields.get(c).and_then(|s| s.parse().ok()).unwrap_or(0);
                let b_val: i64 = b_fields.get(c).and_then(|s| s.parse().ok()).unwrap_or(0);
                ord = ord.then(a_val.cmp(&b_val));
            }
            if let Some(c) = stop_col {
                let a_val: i64 = a_fields.get(c).and_then(|s| s.parse().ok()).unwrap_or(0);
                let b_val: i64 = b_fields.get(c).and_then(|s| s.parse().ok()).unwrap_or(0);
                ord = ord.then(a_val.cmp(&b_val));
            }
            ord = ord.then(get_field(&a_fields, strand_col).cmp(&get_field(&b_fields, strand_col)));
        }

        ord = ord.then(get_field(&a_fields, prot_col).cmp(&get_field(&b_fields, prot_col)));
        ord = ord.then(get_field(&a_fields, gene_col).cmp(&get_field(&b_fields, gene_col)));
        ord = ord.then(a.cmp(b));

        ord
    });
    data_lines.dedup();

    let mut result = String::new();
    result.push_str(header);
    result.push('\n');
    for line in data_lines {
        result.push_str(line);
        result.push('\n');
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pipeline_config_default() {
        let config = PipelineConfig::default();
        assert_eq!(config.threads, 4);
        assert_eq!(config.coverage_min, 0.5);
        assert!(!config.plus);
    }

    #[test]
    fn stx2a_operon_row_is_synthesized_from_stxa_and_stxb_accessions() {
        let mut output = concat!(
            "Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\tSubtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference\t% Identity to reference\tAlignment length\tClosest reference accession\tClosest reference name\tHMM accession\tHMM description\tHierarchy node\n",
            "NA\tcontig18\t279\t1235\t+\tstxA2\tShiga toxin Stx2 subunit A\tplus\tVIRULENCE\tVIRULENCE\tSTX2\tstxA2\tEXACTX\t319\t319\t100.00\t100.00\t319\tTJA36680.1\tShiga toxin Stx2 subunit A\tNA\tNA\tstxA2_acd\n",
            "NA\tcontig18\t1250\t1516\t+\tstxB2\tShiga toxin Stx2a subunit B\tplus\tVIRULENCE\tVIRULENCE\tSTX2\tstxB2a\tEXACTX\t89\t89\t100.00\t100.00\t89\tAAM90978.1\tShiga toxin Stx2a subunit B\tNA\tNA\tstxB2a\n"
        )
        .to_string();

        append_stx_operon_row(&mut output);
        assert!(output.contains("\tstx2a_operon\t"));
        assert!(output
            .contains("\tTJA36680.1,AAM90978.1\tShiga toxin stx2a\tNA\tNA\tstxA2_acd::stxB2a\n"));

        let once = output.clone();
        append_stx_operon_row(&mut output);
        assert_eq!(output, once, "operon synthesis should be idempotent");
    }

    #[test]
    fn stx2a_operon_row_requires_accession_columns_not_reference_names() {
        let mut output = concat!(
            "Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\tSubtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference\t% Identity to reference\tAlignment length\tClosest reference accession\tClosest reference name\tHMM accession\tHMM description\tHierarchy node\n",
            "NA\tcontig18\t279\t1235\t+\tstxA2\tShiga toxin Stx2 subunit A\tplus\tVIRULENCE\tVIRULENCE\tSTX2\tstxA2\tEXACTX\t319\t319\t100.00\t100.00\t319\tWRONG\tTJA36680.1\tNA\tNA\tstxA2_acd\n",
            "NA\tcontig18\t1250\t1516\t+\tstxB2\tShiga toxin Stx2a subunit B\tplus\tVIRULENCE\tVIRULENCE\tSTX2\tstxB2a\tEXACTX\t89\t89\t100.00\t100.00\t89\tWRONG\tAAM90978.1\tNA\tNA\tstxB2a\n"
        )
        .to_string();

        append_stx_operon_row(&mut output);
        assert!(output.contains("\tstx2a_operon\t"));
        assert!(
            output.contains("\tWRONG,WRONG\tShiga toxin stx2a\tNA\tNA\tstxA2_acd::stxB2a\n"),
            "operon synthesis should read accession columns, not reference-name columns"
        );
    }
}
