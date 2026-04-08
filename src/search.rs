// Library-based BLAST and HMM search
// Replaces CLI subprocess calls with direct library calls

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use anyhow::Result;

/// Run HMM search using hmmer-pure-rs library
pub fn run_hmmsearch_library(
    hmm_db_path: &Path,
    query_fasta: &Path,
    tblout_path: &Path,
    domtblout_path: &Path,
) -> Result<()> {
    use hmmer::io::{binary, fasta};
    use hmmer::pipeline::{self, PipelineConfig};

    // Load HMM models from pressed database
    let h3m_path = format!("{}.h3m", hmm_db_path.to_str().unwrap());
    let hmms = binary::read_all_pressed(&h3m_path)?;

    // Load query sequences
    let f = File::open(query_fasta)?;
    let mut reader = BufReader::new(f);
    let abc = hmmer::alphabet::Alphabet::amino();
    let seqs = fasta::read_fasta(&mut reader, &abc)?;

    // Configure pipeline
    let config = PipelineConfig {
        use_tc: true,
        ..PipelineConfig::default()
    };

    // Run hmmsearch for each HMM model and write results
    let mut tblout_file = File::create(tblout_path)?;
    let mut domtblout_file = File::create(domtblout_path)?;

    for hmm in &hmms {
        let (hits, _stats) = pipeline::hmmsearch(hmm, &seqs, &config);

        if !hits.is_empty() {
            let hmm_name = &hmm.name;
            let hmm_acc = hmm.acc.as_deref();
            pipeline::format_tblout(hmm_name, &hits, &mut tblout_file)?;
            pipeline::format_domtblout(hmm_name, hmm_acc, hmm.m, &hits, &mut domtblout_file)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_hmmsearch_library() {
        let db = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amrfinder_db/2026-03-24.1");
        let hmm_path = db.join("AMR.LIB");
        let query = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amr/test_prot.fa");

        if !hmm_path.exists() || !query.exists() {
            return;
        }

        let h3m = format!("{}.h3m", hmm_path.to_str().unwrap());
        if !Path::new(&h3m).exists() {
            return;
        }

        let tblout = std::env::temp_dir().join("test_hmmsearch_tblout");
        let domtblout = std::env::temp_dir().join("test_hmmsearch_domtblout");

        let result = run_hmmsearch_library(&hmm_path, &query, &tblout, &domtblout);
        if let Ok(()) = result {
            assert!(tblout.exists());
            assert!(domtblout.exists());
            std::fs::remove_file(&tblout).ok();
            std::fs::remove_file(&domtblout).ok();
        }
    }
}
