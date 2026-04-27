// Library-based BLAST and HMM search
// Replaces CLI subprocess calls with direct library calls

use std::fs::File;
use std::io::Write;
use std::path::Path;

use anyhow::{anyhow, Result};

/// Run HMM search using hmmer-pure-rs library
pub fn run_hmmsearch_library(
    hmm_db_path: &Path,
    query_fasta: &Path,
    tblout_path: &Path,
    domtblout_path: &Path,
) -> Result<()> {
    use hmmer::pipeline::{BitCutoff, Pipeline, ZSetBy};
    use hmmer::profile::{self, P7_LOCAL};
    use hmmer::{Alphabet, Bg, OProfile, Profile, Sequence, TopHits};

    // Load HMM models from pressed database
    let h3m_path = format!("{}.h3m", hmm_db_path.to_str().unwrap());
    let hmms = hmmer::hmmfile_binary::read_binary_hmm_file(Path::new(&h3m_path))
        .or_else(|_| hmmer::hmmfile::read_hmm_file(hmm_db_path))
        .map_err(|e| anyhow!("failed to read HMM database: {e}"))?;

    // Run hmmsearch for each HMM model and write results
    let mut tblout_file = File::create(tblout_path)?;
    let mut domtblout_file = File::create(domtblout_path)?;

    for hmm in &hmms {
        let abc = Alphabet::new(hmm.abc_type);
        let mut bg = Bg::new(&abc);
        let mut gm = Profile::new(hmm.m, &abc);
        profile::profile_config(hmm, &bg, &mut gm, 100, P7_LOCAL);
        bg.set_filter(hmm.m, &hmm.compo);

        let mut om = OProfile::convert(&gm);
        let mut pli = Pipeline::new();
        pli.new_model(&gm);
        pli.use_bit_cutoffs = BitCutoff::TC;
        if let Err(e) = pli.new_model_thresholds(&hmm.cutoff) {
            return Err(anyhow!("failed to apply TC cutoffs for {}: {e}", hmm.name));
        }

        let mut th = TopHits::new();
        let mut sqf = hmmer::sequence::open_seq_file(query_fasta, &abc)
            .map_err(|e| anyhow!("failed to open query FASTA: {e}"))?;
        let mut sq = Sequence::new();
        let mut n_targets = 0_u64;
        while sqf
            .read(&mut sq)
            .map_err(|e| anyhow!("failed to read query FASTA: {e}"))?
        {
            n_targets += 1;
            bg.set_length(sq.n);
            let mut local_th = TopHits::new();
            if pli.run(&mut gm, &mut om, &bg, hmm, &sq, &mut local_th) {
                th.hits.extend(local_th.hits.into_iter());
            }
            sq.reuse();
        }

        pli.n_targets = n_targets;
        let z = match pli.z_setby {
            ZSetBy::Option => pli.z,
            ZSetBy::Ntargets => pli.n_targets as f64,
        };
        th.sort_by_sortkey();
        th.threshold(&pli, z, z);
        let domz = match pli.domz_setby {
            ZSetBy::Option => pli.domz,
            ZSetBy::Ntargets => th.nreported.max(1) as f64,
        };
        if domz != z {
            th.threshold(&pli, z, domz);
        }

        hmmer::output::write_tblout(&mut tblout_file, &hmm.name, hmm.acc.as_deref(), &th, z);
        hmmer::output::write_domtblout(
            &mut domtblout_file,
            &hmm.name,
            hmm.acc.as_deref(),
            &th,
            z,
            domz,
        );
    }

    tblout_file.flush()?;
    domtblout_file.flush()?;

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

        let tblout = std::env::temp_dir().join("hmm_lib_hits.tblout");
        let domtblout = std::env::temp_dir().join("hmm_lib_hits.domtblout");
        run_hmmsearch_library(&hmm_path, &query, &tblout, &domtblout).unwrap();
        let total_hits = std::fs::read_to_string(&tblout)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
            .count();
        std::fs::remove_file(&tblout).ok();
        std::fs::remove_file(&domtblout).ok();

        assert!(total_hits > 0, "Library hmmsearch should find hits");
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
                "--tblout",
                ext_tblout.to_str().unwrap(),
                "--domtblout",
                ext_dom.to_str().unwrap(),
                "--noali",
                "--cut_tc",
                "-Z",
                "10000",
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
        for key in lib_hits.keys() {
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
