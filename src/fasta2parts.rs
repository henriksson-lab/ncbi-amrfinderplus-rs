use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use anyhow::{bail, Result};

/// Split FASTA file into parts without breaking sequences.
/// Matches the C++ fasta2parts binary exactly.
pub fn body(fasta_path: &Path, parts_max: usize, out_dir: &Path) -> Result<()> {
    if parts_max <= 1 {
        bail!("Number of parts must be >= 2");
    }

    let file_size = fs::metadata(fasta_path)?.len() as usize;
    let chunk_min = file_size / parts_max + 1;

    let mut part: usize = 0;
    let mut out: Option<BufWriter<File>> = None;
    let mut seq_size: usize = 0;

    let file = File::open(fasta_path)?;
    let mut reader = BufReader::new(file);

    let mut line = Vec::new();
    while reader.read_until(b'\n', &mut line)? != 0 {
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

        if line[0] == b'>' && seq_size >= chunk_min && part < parts_max {
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
            f.write_all(&line)?;
            f.write_all(b"\n")?;
        }
        seq_size += line.len();
        line.clear();
    }

    Ok(())
}
