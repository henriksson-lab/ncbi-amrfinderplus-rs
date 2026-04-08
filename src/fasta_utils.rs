use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use anyhow::{bail, Result};

/// Options for fasta_check
pub struct FastaCheckOpts<'a> {
    pub fasta_path: &'a Path,
    pub aa: bool,
    pub hyphen: bool,
    pub ambig: bool,
    pub ambig_max: usize,
    pub stop_codon: bool,
    pub len_path: Option<&'a Path>,
    pub out_path: Option<&'a Path>,
}

/// Check the correctness of a FASTA file.
/// Returns (num_sequences, max_sequence_length, total_sequence_length).
/// Matches the C++ fasta_check binary exactly.
pub fn fasta_check(opts: &FastaCheckOpts) -> Result<(usize, usize, usize)> {
    let FastaCheckOpts {
        fasta_path,
        aa,
        hyphen,
        ambig,
        ambig_max,
        stop_codon,
        len_path,
        out_path,
    } = opts;
    let aa = *aa;
    let hyphen = *hyphen;
    let ambig = *ambig;
    let ambig_max = *ambig_max;
    let stop_codon = *stop_codon;
    if stop_codon && !aa {
        bail!("stop_codon requires aa");
    }

    let mut len_file = len_path
        .map(|p| File::create(p).map(BufWriter::new))
        .transpose()?;
    let mut out_file = out_path
        .map(|p| File::create(p).map(BufWriter::new))
        .transpose()?;

    let mut ids: Vec<String> = Vec::with_capacity(100_000);
    let mut seq_size_max: usize = 0;
    let mut seq_size_sum: usize = 0;

    // Per-sequence state
    let mut xs: usize = 0;
    let mut header = String::new();
    let mut seq = String::new();
    let mut lines: usize = 0;
    let mut nuc: usize = 0;

    // Inline processSeq logic as a macro to match the C++ lambda behavior
    // (captures all local variables by mutable reference)
    macro_rules! process_seq {
        () => {
            if lines > 0 {
                assert!(!header.is_empty());
                assert!(!ids.is_empty());
                let id = ids.last().unwrap().clone();

                if aa && !stop_codon {
                    while !seq.is_empty() && seq.ends_with('*') {
                        if out_file.is_some() {
                            seq.pop();
                        } else {
                            bail!("{}: '*' at the sequence end", id);
                        }
                    }
                }
                if seq.is_empty() {
                    bail!("{}: Empty sequence", id);
                }
                let mut skip = false;
                if !ambig && xs > ambig_max {
                    if out_file.is_some() {
                        skip = true;
                    } else {
                        bail!("{}: Too many ambiguities", id);
                    }
                }
                if skip {
                    eprintln!("Skipping {}", id);
                } else {
                    if let Some(ref mut lf) = len_file {
                        writeln!(lf, "{}\t{}", id, seq.len())?;
                    }
                    if let Some(ref mut of_) = out_file {
                        writeln!(of_, "{}", header)?;
                        writeln!(of_, "{}", seq)?;
                    }
                    if seq.len() > seq_size_max {
                        seq_size_max = seq.len();
                    }
                    seq_size_sum += seq.len();
                }
                #[allow(unused_assignments)]
                { xs = 0; }
                header.clear();
                seq.clear();
            }
        };
    }

    {
        let file = File::open(fasta_path)?;
        let reader = BufReader::new(file);
        let mut line_num: usize = 0;

        for line_result in reader.lines() {
            let mut line = line_result?;
            // trimTrailing: remove trailing whitespace
            let trimmed_len = line.trim_end().len();
            line.truncate(trimmed_len);

            if line.is_empty() {
                line_num += 1;
                continue;
            }

            let error_prefix =
                format!("File {}, line {}: ", fasta_path.display(), line_num);

            if let Some(after_gt) = line.strip_prefix('>') {
                let pos = after_gt
                    .find(|c: char| c.is_ascii_whitespace())
                    .unwrap_or(after_gt.len());
                let id = &after_gt[..pos];
                if id.is_empty() {
                    bail!("{}Empty sequence identifier", error_prefix);
                }
                for c in id.chars() {
                    if !printable(c) {
                        bail!(
                            "{}Non-printable character in the sequence identifier: {}",
                            error_prefix,
                            c as u32
                        );
                    }
                }
                // BLAST: PD-4548
                if !aa {
                    if id.starts_with('?') {
                        bail!("{}Sequence identifier starts with '?'", error_prefix);
                    }
                    for c in [',', ';', '.', '~'] {
                        if id.ends_with(c) {
                            bail!(
                                "{}Sequence identifier ends with \"{}\"",
                                error_prefix,
                                c
                            );
                        }
                    }
                    if id.contains("\\t") {
                        bail!("{}Sequence identifier contains '\\t'", error_prefix);
                    }
                    if id.contains(",,") {
                        bail!("{}Sequence identifier contains ',,'", error_prefix);
                    }
                }

                process_seq!();
                header = line.clone();
                ids.push(id.to_string());
            } else {
                if lines == 0 {
                    bail!("{}FASTA should start with '>'", error_prefix);
                }
                for c in line.chars() {
                    let mut skip = false;
                    if c == '-' {
                        if hyphen {
                            // allowed
                        } else if out_file.is_some() {
                            skip = true;
                        } else {
                            bail!("{}Hyphen in the sequence", error_prefix);
                        }
                    } else {
                        let c1 = c.to_ascii_lowercase();
                        if aa {
                            if !"acdefghiklmnpqrstvwyxbzjuo*".contains(c1) {
                                bail!(
                                    "{}Wrong amino acid character: (code = {}) '{}'",
                                    error_prefix,
                                    c as u32,
                                    c
                                );
                            }
                            if "acgt".contains(c1) {
                                nuc += 1;
                            }
                            if "xbzjuo".contains(c1) {
                                xs += 1;
                            }
                        } else {
                            if !"acgtbdhkmnrsvwy".contains(c1) {
                                bail!(
                                    "{}Wrong nucleotide character: (code = {}) '{}'",
                                    error_prefix,
                                    c as u32,
                                    c
                                );
                            }
                            if "bdhkmnrsvwy".contains(c1) {
                                xs += 1;
                            }
                        }
                    }
                    if !skip {
                        seq.push(c);
                    }
                }
            }
            lines += 1;
        }
    }

    // Process last sequence
    process_seq!();

    if lines == 0 {
        bail!("Empty file");
    }
    if aa && (nuc as f64) / (seq_size_sum as f64) > 0.9 {
        bail!("Protein sequences looks like a nucleotide sequences");
    }

    // Check for duplicate identifiers
    ids.sort();
    for i in 1..ids.len() {
        if ids[i] == ids[i - 1] {
            bail!("Duplicate identifier: {}", ids[i]);
        }
    }

    Ok((ids.len(), seq_size_max, seq_size_sum))
}

// --- fasta_extract ---

struct Segment {
    start: usize,
    stop: usize,
    strand: bool,
    genesymbol: String,
    name: String,
}

impl Segment {
    fn is_dna(&self) -> bool {
        self.stop > 0
    }
}

fn complementary_nucleotide(c: char) -> Result<char> {
    let r = match c.to_ascii_lowercase() {
        'a' => 't', 'c' => 'g', 'g' => 'c', 't' => 'a',
        'm' => 'k', 'r' => 'y', 'w' => 'w', 's' => 's',
        'y' => 'r', 'k' => 'm', 'v' => 'b', 'h' => 'd',
        'd' => 'h', 'b' => 'v', 'n' => 'n', '-' => '-',
        _ => bail!("Bad nucleotide {}", c),
    };
    if c.is_ascii_uppercase() {
        Ok(r.to_ascii_uppercase())
    } else {
        Ok(r)
    }
}

/// Extract sequences from a FASTA file.
/// Matches the C++ fasta_extract binary exactly.
pub fn fasta_extract(
    fasta_path: &Path,
    target_path: &Path,
    aa: bool,
    out: &mut dyn Write,
) -> Result<()> {
    // Parse target file
    let mut id2segments: HashMap<String, Vec<Segment>> = HashMap::new();
    {
        let file = File::open(target_path)?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.is_empty() {
                continue;
            }
            let id = fields[0].to_string();
            let seg = if aa {
                if fields.len() < 3 {
                    bail!("AA target line needs: <id> <gene> <name>");
                }
                Segment {
                    start: 0,
                    stop: 0,
                    strand: true,
                    genesymbol: fields[1].to_string(),
                    name: fields[2..].join(" "),
                }
            } else {
                if fields.len() < 6 {
                    bail!("DNA target line needs: <id> <start> <stop> <strand> <gene> <name>");
                }
                let start: usize = fields[1].parse::<usize>()?;
                let stop: usize = fields[2].parse()?;
                let strand_char = fields[3];
                if start == 0 {
                    bail!("start must be >= 1");
                }
                if start > stop {
                    bail!("start must be <= stop");
                }
                if strand_char != "+" && strand_char != "-" {
                    bail!("strand must be '+' or '-'");
                }
                Segment {
                    start: start - 1, // 1-based to 0-based
                    stop,
                    strand: strand_char == "+",
                    genesymbol: fields[4].to_string(),
                    name: fields[5..].join(" "),
                }
            };
            id2segments.entry(id).or_default().push(seg);
        }
    }

    if id2segments.is_empty() {
        return Ok(());
    }

    let mut processed: usize = 0;
    {
        let file = File::open(fasta_path)?;
        let reader = BufReader::new(file);
        let mut current_id = String::new();
        let mut seq = String::new();

        let process = |id: &str, seq: &mut String, out: &mut dyn Write| -> Result<bool> {
            if id.is_empty() {
                return Ok(false);
            }
            let segments = match id2segments.get(id) {
                Some(s) => s,
                None => return Ok(false),
            };

            // Remove hyphens
            seq.retain(|c| c != '-');
            assert!(!seq.is_empty());

            for seg in segments {
                write!(out, ">{}", id)?;
                if seg.is_dna() {
                    let start = seg.start.min(seq.len());
                    let stop = seg.stop.min(seq.len());
                    assert!(start < stop);
                    write!(
                        out,
                        ":{}−{} strand:{}",
                        start + 1,
                        stop,
                        if seg.strand { '+' } else { '-' }
                    )?;
                }
                writeln!(out, " {} {}", seg.genesymbol, seg.name)?;

                let mut seq1 = seq.clone();
                if seg.is_dna() {
                    assert!(seg.stop <= seq1.len());
                    seq1 = seq1[seg.start..seg.stop].to_string();
                    if !seg.strand {
                        seq1 = seq1
                            .chars()
                            .rev()
                            .map(|c| complementary_nucleotide(c).unwrap())
                            .collect();
                    }
                }
                let line_len = 60;
                for chunk in seq1.as_bytes().chunks(line_len) {
                    writeln!(out, "{}", std::str::from_utf8(chunk).unwrap())?;
                }
            }

            Ok(true)
        };

        for line_result in reader.lines() {
            let mut line = line_result?;
            let trimmed_len = line.trim_end().len();
            line.truncate(trimmed_len);
            if line.is_empty() {
                continue;
            }
            if let Some(after_gt) = line.strip_prefix('>') {
                if process(&current_id, &mut seq, out)? {
                    processed += 1;
                }
                let pos = after_gt
                    .find(|c: char| c.is_ascii_whitespace())
                    .unwrap_or(after_gt.len());
                current_id = after_gt[..pos].to_string();
                seq.clear();
            } else {
                seq.push_str(&line);
            }
        }
        if process(&current_id, &mut seq, out)? {
            processed += 1;
        }
    }

    if processed != id2segments.len() {
        bail!(
            "Requested identifiers: {}, but processed: {}",
            id2segments.len(),
            processed
        );
    }

    Ok(())
}

// --- fasta2parts ---

/// Split FASTA file into parts without breaking sequences.
/// Matches the C++ fasta2parts binary exactly.
pub fn fasta2parts(
    fasta_path: &Path,
    parts_max: usize,
    out_dir: &Path,
) -> Result<()> {
    if parts_max <= 1 {
        bail!("Number of parts must be >= 2");
    }

    let file_size = fs::metadata(fasta_path)?.len() as usize;
    let chunk_min = file_size / parts_max + 1;

    let mut part: usize = 0;
    let mut out: Option<BufWriter<File>> = None;
    let mut seq_size: usize = 0;

    let file = File::open(fasta_path)?;
    let reader = BufReader::new(file);

    for line_result in reader.lines() {
        let mut line = line_result?;
        let trimmed_len = line.trim_end().len();
        line.truncate(trimmed_len);
        if line.is_empty() {
            continue;
        }

        if line.starts_with('>') && seq_size >= chunk_min && part < parts_max {
            out = None;
            seq_size = 0;
        }

        if out.is_none() {
            part += 1;
            assert!(part <= parts_max);
            let out_path = out_dir.join(part.to_string());
            out = Some(BufWriter::new(File::create(out_path)?));
        }

        if let Some(ref mut f) = out {
            writeln!(f, "{}", line)?;
        }
        seq_size += line.len();
    }

    Ok(())
}

fn printable(c: char) -> bool {
    let code = c as u32;
    (32..127).contains(&code)
}

/// Write fasta_check output to stdout in the same format as the C++ binary.
pub fn fasta_check_print(num_seqs: usize, max_len: usize, total_len: usize) {
    println!("{}", num_seqs);
    println!("{}", max_len);
    println!("{}", total_len);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use std::process::Command;

    fn test_data_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amr")
    }

    fn cpp_binary(name: &str) -> PathBuf {
        test_data_dir().join(name)
    }

    #[test]
    fn test_fasta_check_protein() {
        let path = test_data_dir().join("test_prot.fa");
        if !path.exists() {
            return;
        }
        let result = fasta_check(&FastaCheckOpts {
            fasta_path: &path, aa: true, hyphen: false, ambig: false,
            ambig_max: 0, stop_codon: false, len_path: None, out_path: None,
        });
        assert!(result.is_ok(), "fasta_check failed: {:?}", result.err());
    }

    #[test]
    fn test_fasta_check_dna() {
        let path = test_data_dir().join("test_dna.fa");
        if !path.exists() {
            return;
        }
        let result = fasta_check(&FastaCheckOpts {
            fasta_path: &path, aa: false, hyphen: false, ambig: true,
            ambig_max: 0, stop_codon: false, len_path: None, out_path: None,
        });
        assert!(result.is_ok(), "fasta_check failed: {:?}", result.err());
    }

    #[test]
    fn test_fasta_check_protein_matches_cpp() {
        let path = test_data_dir().join("test_prot.fa");
        let cpp_bin = cpp_binary("fasta_check");
        if !path.exists() || !cpp_bin.exists() {
            return;
        }

        let cpp_output = Command::new(&cpp_bin)
            .arg(&path)
            .arg("-aa")
            .output()
            .expect("failed to run C++ fasta_check");
        let cpp_stdout = String::from_utf8_lossy(&cpp_output.stdout);

        let (num_seqs, max_len, total_len) = fasta_check(&FastaCheckOpts {
            fasta_path: &path, aa: true, hyphen: false, ambig: false,
            ambig_max: 0, stop_codon: false, len_path: None, out_path: None,
        }).unwrap();
        let rust_stdout = format!("{}\n{}\n{}\n", num_seqs, max_len, total_len);

        assert_eq!(cpp_stdout, rust_stdout,
            "C++ and Rust fasta_check output differ for protein");
    }

    #[test]
    fn test_fasta_check_dna_matches_cpp() {
        let path = test_data_dir().join("test_dna.fa");
        let cpp_bin = cpp_binary("fasta_check");
        if !path.exists() || !cpp_bin.exists() {
            return;
        }

        let cpp_output = Command::new(&cpp_bin)
            .arg(&path)
            .arg("-ambig")
            .output()
            .expect("failed to run C++ fasta_check");
        let cpp_stdout = String::from_utf8_lossy(&cpp_output.stdout);

        let (num_seqs, max_len, total_len) = fasta_check(&FastaCheckOpts {
            fasta_path: &path, aa: false, hyphen: false, ambig: true,
            ambig_max: 0, stop_codon: false, len_path: None, out_path: None,
        }).unwrap();
        let rust_stdout = format!("{}\n{}\n{}\n", num_seqs, max_len, total_len);

        assert_eq!(cpp_stdout, rust_stdout,
            "C++ and Rust fasta_check output differ for DNA");
    }
}
