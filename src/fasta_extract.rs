use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use anyhow::{bail, Result};

struct Segment {
    start: usize,
    stop: usize,
    strand: bool,
    gene_symbol: String,
    name: String,
}

impl Segment {
    fn is_dna(&self) -> bool {
        self.stop > 0
    }

    fn size(&self) -> usize {
        self.stop - self.start
    }

    fn save_text(&self, out: &mut dyn Write) -> Result<()> {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            self.start, self.stop, self.strand, self.gene_symbol, self.name
        )?;
        Ok(())
    }
}

fn complementary_nucleotide(c: char) -> Result<char> {
    let r = match c.to_ascii_lowercase() {
        'a' => 't',
        'c' => 'g',
        'g' => 'c',
        't' => 'a',
        'm' => 'k',
        'r' => 'y',
        'w' => 'w',
        's' => 's',
        'y' => 'r',
        'k' => 'm',
        'v' => 'b',
        'h' => 'd',
        'd' => 'h',
        'b' => 'v',
        'n' => 'n',
        '-' => '-',
        _ => bail!("Bad nucleotide {}", c),
    };
    if c.is_ascii_uppercase() {
        Ok(r.to_ascii_uppercase())
    } else {
        Ok(r)
    }
}

pub fn body(
    fasta_path: &Path,
    target_path: &Path,
    aa: bool,
    verbose: bool,
    out: &mut dyn Write,
) -> Result<()> {
    // Parse target file
    let mut id2segments: BTreeMap<String, Vec<Segment>> = BTreeMap::new();
    {
        let file = File::open(target_path)?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.is_empty() {
                bail!("Empty target line");
            }
            let id = fields[0].to_string();
            let tokens_before_name = if aa { 2 } else { 5 };
            let mut name_start = 0;
            for _ in 0..tokens_before_name {
                let token_start = line[name_start..]
                    .find(|c: char| !c.is_whitespace())
                    .map(|offset| name_start + offset)
                    .unwrap_or(line.len());
                let token_stop = line[token_start..]
                    .find(|c: char| c.is_whitespace())
                    .map(|offset| token_start + offset)
                    .unwrap_or(line.len());
                name_start = token_stop;
            }
            let name = line[name_start..].trim().to_string();
            let seg = if aa {
                if fields.len() < 3 {
                    bail!("AA target line needs: <id> <gene> <name>");
                }
                Segment {
                    start: 0,
                    stop: 0,
                    strand: true,
                    gene_symbol: fields[1].to_string(),
                    name,
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
                    gene_symbol: fields[4].to_string(),
                    name,
                }
            };
            id2segments.entry(id).or_default().push(seg);
        }
    }

    if id2segments.is_empty() {
        return Ok(());
    }

    if verbose {
        for (id, segments) in &id2segments {
            writeln!(out, "{}: ", id)?;
            for seg in segments {
                write!(out, "  ")?;
                seg.save_text(out)?;
            }
        }
    }

    let mut processed: usize = 0;
    {
        let file = File::open(fasta_path)?;
        let reader = BufReader::new(file);
        let mut current_id = String::new();
        let mut seq = String::new();

        for line_result in reader.lines() {
            let mut line = line_result?;
            let trimmed_len = line.trim_end().len();
            line.truncate(trimmed_len);
            if line.is_empty() {
                continue;
            }
            if let Some(after_gt) = line.strip_prefix('>') {
                if process(&current_id, &mut seq, &mut id2segments, out)? {
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
        if process(&current_id, &mut seq, &mut id2segments, out)? {
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

fn process(
    id: &str,
    seq: &mut String,
    id2segments: &mut BTreeMap<String, Vec<Segment>>,
    out: &mut dyn Write,
) -> Result<bool> {
    if id.is_empty() {
        return Ok(false);
    }
    let segments = match id2segments.get_mut(id) {
        Some(s) => s,
        None => return Ok(false),
    };

    seq.retain(|c| c != '-');
    assert!(!seq.is_empty());

    for seg in segments {
        write!(out, ">{}", id)?;
        if seg.is_dna() {
            assert!(seg.start <= seq.len());
            seg.stop = seg.stop.min(seq.len());
            assert!(seg.start < seg.stop);
            write!(
                out,
                ":{}-{} strand:{}",
                seg.start + 1,
                seg.stop,
                if seg.strand { '+' } else { '-' }
            )?;
        }
        writeln!(out, " {} {}", seg.gene_symbol, seg.name)?;

        let mut seq1 = seq.clone();
        if seg.is_dna() {
            assert!(seg.stop <= seq1.len());
            seq1 = seq1[seg.start..seg.start + seg.size()].to_string();
            if !seg.strand {
                let mut reverse_complement = String::with_capacity(seq1.len());
                for c in seq1.chars().rev() {
                    reverse_complement.push(complementary_nucleotide(c)?);
                }
                seq1 = reverse_complement;
            }
        }
        let line_len = 60;
        for chunk in seq1.as_bytes().chunks(line_len) {
            writeln!(out, "{}", std::str::from_utf8(chunk).unwrap())?;
        }
    }

    Ok(true)
}
