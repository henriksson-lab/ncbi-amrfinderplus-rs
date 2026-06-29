use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use anyhow::{bail, Result};

fn process_seq(
    lines: usize,
    aa: bool,
    ambig: bool,
    ambig_max: usize,
    stop_codon: bool,
    ids: &[String],
    xs: &mut usize,
    header: &mut Vec<u8>,
    seq: &mut Vec<u8>,
    len_file: &mut Option<BufWriter<File>>,
    out_file: &mut Option<BufWriter<File>>,
    seq_size_max: &mut usize,
    seq_size_sum: &mut usize,
) -> Result<()> {
    if lines == 0 {
        return Ok(());
    }
    assert!(!header.is_empty());
    assert!(!ids.is_empty());
    let id = ids.last().unwrap().clone();

    if aa && !stop_codon {
        while !seq.is_empty() && seq.last() == Some(&b'*') {
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
    if !ambig && *xs > ambig_max {
        if out_file.is_some() {
            skip = true;
        } else {
            bail!("{}: Too many ambiguities", id);
        }
    }
    if !skip {
        if let Some(lf) = len_file {
            writeln!(lf, "{}\t{}", id, seq.len())?;
        }
        if let Some(of_) = out_file {
            of_.write_all(header)?;
            of_.write_all(b"\n")?;
            of_.write_all(seq)?;
            of_.write_all(b"\n")?;
        }
        if seq.len() > *seq_size_max {
            *seq_size_max = seq.len();
        }
        *seq_size_sum += seq.len();
    }
    *xs = 0;
    header.clear();
    seq.clear();
    Ok(())
}

/// Check the correctness of a FASTA file.
/// Returns (num_sequences, max_sequence_length, total_sequence_length).
/// Matches the C++ fasta_check binary exactly.
pub fn body(
    fasta_path: &Path,
    aa: bool,
    hyphen: bool,
    ambig: bool,
    ambig_max: usize,
    stop_codon: bool,
    len_path: Option<&Path>,
    out_path: Option<&Path>,
) -> Result<(usize, usize, usize)> {
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
    let mut header = Vec::new();
    let mut seq = Vec::new();
    let mut lines: usize = 0;
    let mut nuc: usize = 0;

    {
        let file = File::open(fasta_path)?;
        let mut reader = BufReader::new(file);
        let mut line_num: usize = 0;
        let mut line = Vec::new();

        while reader.read_until(b'\n', &mut line)? != 0 {
            line_num += 1;
            while line
                .last()
                .is_some_and(|c| *c > b'\0' && *c <= b' ' && c.is_ascii_whitespace())
            {
                line.pop();
            }

            if line.is_empty() {
                line.clear();
                continue;
            }

            let error_prefix = format!("File {}, line {}: ", fasta_path.display(), line_num);

            if line[0] == b'>' {
                let after_gt = &line[1..];
                let pos = after_gt
                    .iter()
                    .position(|c| c.is_ascii_whitespace())
                    .unwrap_or(after_gt.len());
                let id_bytes = &after_gt[..pos];
                if id_bytes.is_empty() {
                    bail!("{}Empty sequence identifier", error_prefix);
                }
                for &c in id_bytes {
                    if !(32..127).contains(&c) {
                        bail!(
                            "{}Non-printable character in the sequence identifier: {}",
                            error_prefix,
                            c as i8
                        );
                    }
                }
                let id = String::from_utf8(id_bytes.to_vec())?;
                // BLAST: PD-4548
                if !aa {
                    if id.starts_with('?') {
                        bail!("{}Sequence identifier starts with '?'", error_prefix);
                    }
                    for c in [',', ';', '.', '~'] {
                        if id.ends_with(c) {
                            bail!("{}Sequence identifier ends with \"{}\"", error_prefix, c);
                        }
                    }
                    if id.contains("\\t") {
                        bail!("{}Sequence identifier contains '\\t'", error_prefix);
                    }
                    if id.contains(",,") {
                        bail!("{}Sequence identifier contains ',,'", error_prefix);
                    }
                }

                process_seq(
                    lines,
                    aa,
                    ambig,
                    ambig_max,
                    stop_codon,
                    &ids,
                    &mut xs,
                    &mut header,
                    &mut seq,
                    &mut len_file,
                    &mut out_file,
                    &mut seq_size_max,
                    &mut seq_size_sum,
                )?;
                header = line.clone();
                ids.push(id);
            } else {
                if lines == 0 {
                    bail!("{}FASTA should start with '>'", error_prefix);
                }
                for &c in &line {
                    let mut skip = false;
                    if c == b'-' {
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
                            if !b"acdefghiklmnpqrstvwyxbzjuo*".contains(&c1) {
                                bail!(
                                    "{}Wrong amino acid character: (code = {}) '{}'",
                                    error_prefix,
                                    c as i8,
                                    c as char
                                );
                            }
                            if b"acgt".contains(&c1) {
                                nuc += 1;
                            }
                            if b"xbzjuo".contains(&c1) {
                                xs += 1;
                            }
                        } else {
                            if !b"acgtbdhkmnrsvwy".contains(&c1) {
                                bail!(
                                    "{}Wrong nucleotide character: (code = {}) '{}'",
                                    error_prefix,
                                    c as i8,
                                    c as char
                                );
                            }
                            if b"bdhkmnrsvwy".contains(&c1) {
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
            line.clear();
        }
    }

    // Process last sequence
    process_seq(
        lines,
        aa,
        ambig,
        ambig_max,
        stop_codon,
        &ids,
        &mut xs,
        &mut header,
        &mut seq,
        &mut len_file,
        &mut out_file,
        &mut seq_size_max,
        &mut seq_size_sum,
    )?;

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
