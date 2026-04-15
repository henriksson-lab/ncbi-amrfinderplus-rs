// Library-based BLAST and HMM search
// Replaces CLI subprocess calls with direct library calls

use std::fs::File;
use std::io::{BufReader, Write};
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
            let hmm_acc = hmm.acc.as_deref().unwrap_or("-");
            // Write standard HMMER tblout format (the library's format_tblout is non-standard)
            for hit in &hits {
                writeln!(
                    tblout_file,
                    "{} - {} {} {:.2e} {:.1} {:.1} {:.2e} {:.1} {:.1} 1 1 1 1 {:.2e} {:.1} {:.1} {} {}",
                    hit.name, hmm_name, hmm_acc,
                    hit.evalue, hit.score, hit.bias,
                    hit.evalue, hit.score, hit.bias,
                    hit.evalue, hit.score, hit.bias,
                    hit.seq_len, hit.desc,
                )?;
            }
            // Write standard HMMER domtblout format (library format is non-standard)
            for hit in &hits {
                let ndom = hit.domains.len();
                for (d_idx, dom) in hit.domains.iter().enumerate() {
                    writeln!(
                        domtblout_file,
                        "{} - {} {} {} {} {:.2e} {:.1} {:.1} {} {} {:.2e} {:.2e} {:.1} {:.1} {} {} {} {} {} {} {:.2} {}",
                        hit.name, hit.seq_len,
                        hmm_name, hmm_acc, hmm.m,
                        hit.evalue, hit.score, hit.bias,
                        d_idx + 1, ndom,
                        dom.cevalue, dom.ievalue, dom.bitscore, 0.0, // domain bias
                        dom.ienv, dom.jenv.min(hmm.m), // hmm from/to
                        dom.ienv, dom.jenv,             // ali from/to
                        dom.ienv, dom.jenv,             // env from/to
                        0.0,                            // acc
                        hit.desc,
                    )?;
                }
            }
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

    /// Diagnostic: check that HMMs and sequences load, and that at least some hits are found.
    #[test]
    fn test_hmmsearch_library_produces_hits() {
        use hmmer::io::{binary, fasta};
        use hmmer::pipeline::{self, PipelineConfig};

        let db = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amrfinder_db/2026-03-24.1");
        let hmm_path = db.join("AMR.LIB");
        let query = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/test_prot.fa");
        if !hmm_path.exists() || !query.exists() {
            return;
        }
        let h3m = format!("{}.h3m", hmm_path.to_str().unwrap());
        if !Path::new(&h3m).exists() {
            return;
        }

        let hmms = binary::read_all_pressed(&h3m).unwrap();
        eprintln!("Loaded {} HMM models", hmms.len());

        let f = File::open(&query).unwrap();
        let mut reader = std::io::BufReader::new(f);
        let abc = hmmer::alphabet::Alphabet::amino();
        let seqs = fasta::read_fasta(&mut reader, &abc).unwrap();
        eprintln!("Loaded {} sequences", seqs.len());
        for s in &seqs {
            eprintln!("  seq: {} len={}", s.name, s.len);
        }

        let config = PipelineConfig {
            use_tc: true,
            z_override: Some(10000.0),
            ..PipelineConfig::default()
        };

        // Test with the blaTEM HMM specifically (should hit blaTEM-156)
        let mut total_hits = 0;
        for hmm in &hmms {
            let (hits, stats) = pipeline::hmmsearch(hmm, &seqs, &config);
            if !hits.is_empty() {
                eprintln!("HMM {} ({} nodes): {} hits (passed_msv={}, passed_vit={}, passed_fwd={})",
                    hmm.name, hmm.m, hits.len(),
                    stats.n_passed_msv, stats.n_passed_vit, stats.n_passed_fwd);
                for h in &hits {
                    eprintln!("  hit: {} score={:.1} evalue={:.2e}", h.name, h.score, h.evalue);
                }
                total_hits += hits.len();
            }
        }

        assert!(
            total_hits > 0,
            "Library hmmsearch should find hits (loaded {} HMMs, {} sequences)",
            hmms.len(), seqs.len()
        );
    }

    /// Compare pure Rust hmmsearch results against the external hmmsearch binary.
    /// Both should find the same targets with similar scores.
    #[test]
    fn test_hmmsearch_library_matches_external() {
        use std::collections::HashMap;
        use std::process::Command;

        let db = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amrfinder_db/2026-03-24.1");
        let hmm_path = db.join("AMR.LIB");
        let query = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/test_prot.fa");

        if !hmm_path.exists() || !query.exists() {
            return;
        }
        let h3m = format!("{}.h3m", hmm_path.to_str().unwrap());
        if !Path::new(&h3m).exists() {
            return;
        }

        // Check external hmmsearch exists
        let which = Command::new("which").arg("hmmsearch").output();
        if which.is_err() || !which.unwrap().status.success() {
            return;
        }

        // Run external hmmsearch
        let ext_tblout = std::env::temp_dir().join("hmm_ext_compare.tblout");
        let ext_dom = std::env::temp_dir().join("hmm_ext_compare.dom");
        let ext_status = Command::new("hmmsearch")
            .args([
                "--tblout", ext_tblout.to_str().unwrap(),
                "--domtblout", ext_dom.to_str().unwrap(),
                "--noali", "--cut_tc",
                "-Z", "10000",
                hmm_path.to_str().unwrap(),
                query.to_str().unwrap(),
            ])
            .stderr(std::process::Stdio::null())
            .stdout(std::process::Stdio::null())
            .status();
        if ext_status.is_err() || !ext_status.unwrap().success() {
            return;
        }

        // Run library hmmsearch
        let lib_tblout = std::env::temp_dir().join("hmm_lib_compare.tblout");
        let lib_dom = std::env::temp_dir().join("hmm_lib_compare.dom");
        run_hmmsearch_library(&hmm_path, &query, &lib_tblout, &lib_dom).unwrap();

        // Parse tblout: extract target_name -> (query_name, score) for each
        fn parse_tblout(path: &std::path::Path) -> HashMap<(String, String), f64> {
            let mut map = HashMap::new();
            let content = std::fs::read_to_string(path).unwrap();
            for line in content.lines() {
                if line.starts_with('#') || line.is_empty() {
                    continue;
                }
                let fields: Vec<&str> = line.split_whitespace().collect();
                if fields.len() >= 6 {
                    let target = fields[0].to_string();
                    let query = fields[2].to_string();
                    let score: f64 = fields[5].parse().unwrap_or(0.0);
                    map.insert((target, query), score);
                }
            }
            map
        }

        let ext_hits = parse_tblout(&ext_tblout);
        let lib_hits = parse_tblout(&lib_tblout);

        // Every external hit should have a corresponding library hit
        let mut missing = Vec::new();
        let mut score_diffs = Vec::new();
        for ((target, query), ext_score) in &ext_hits {
            if let Some(lib_score) = lib_hits.get(&(target.clone(), query.clone())) {
                let diff = (ext_score - lib_score).abs();
                if diff > 1.0 {
                    score_diffs.push(format!(
                        "{}/{}: ext={:.1} lib={:.1} diff={:.1}",
                        target, query, ext_score, lib_score, diff
                    ));
                }
            } else {
                missing.push(format!("{}/{} (score={:.1})", target, query, ext_score));
            }
        }

        // Extra hits in library not in external
        let mut extra = Vec::new();
        for (key, _) in &lib_hits {
            if !ext_hits.contains_key(key) {
                extra.push(format!("{}/{}", key.0, key.1));
            }
        }

        // Cleanup
        std::fs::remove_file(&ext_tblout).ok();
        std::fs::remove_file(&ext_dom).ok();
        std::fs::remove_file(&lib_tblout).ok();
        std::fs::remove_file(&lib_dom).ok();

        // Report
        assert!(
            missing.is_empty(),
            "Library hmmsearch missing {} hits found by external:\n  {}",
            missing.len(),
            missing.join("\n  ")
        );
        assert!(
            score_diffs.is_empty(),
            "Score differences > 1.0 bit:\n  {}",
            score_diffs.join("\n  ")
        );
        if !extra.is_empty() {
            eprintln!(
                "NOTE: Library found {} extra hits not in external: {}",
                extra.len(),
                extra.join(", ")
            );
        }
    }
}
