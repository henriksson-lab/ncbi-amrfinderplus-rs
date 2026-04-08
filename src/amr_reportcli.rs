// amr_report CLI — standalone report generation from BLAST/HMM results
// This module provides a Rust implementation of the C++ amr_report binary

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use anyhow::Result;

use crate::gff::{Annot, GffType};
use crate::report::{Batch, BlastAlignment, BlastRule, HmmAlignment, HmmDomain};
use crate::seq::Hsp;

/// Parse BLASTP/BLASTX tabular output and add to batch
pub fn load_blast_results(
    blast_file: &Path,
    batch: &mut Batch,
    is_protein: bool, // true for blastp, false for blastx
    nosame: bool,
) -> Result<()> {
    let file = File::open(blast_file)?;
    let reader = BufReader::new(file);

    let default_complete_br = BlastRule::default();
    let default_partial_br = BlastRule::default();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        let al = BlastAlignment::from_blast_line(
            trimmed,
            true,           // q_prot (reference is always protein)
            is_protein,     // s_prot (target is protein for blastp, DNA for blastx)
            &default_complete_br,
            &default_partial_br,
        )?;

        if nosame && al.ref_accession == al.hsp.sseqid {
            continue;
        }

        batch.add_blast_al(al);
    }

    Ok(())
}

/// Parse HMM domain table output
pub fn load_hmm_domains(
    dom_file: &Path,
    batch: &mut Batch,
) -> Result<()> {
    let file = File::open(dom_file)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() < 23 {
            continue;
        }

        // --domtblout format:
        // 0:target_name 1:accession 2:tlen 3:query_name 4:query_accession 5:qlen
        // 6:E-value 7:score 8:bias 9:# 10:of 11:c-evalue 12:i-evalue 13:score 14:bias
        // 15:hmm_from 16:hmm_to 17:ali_from 18:ali_to 19:env_from 20:env_to 21:acc
        let target_name = fields[0].to_string();
        let query_accession = fields[4].to_string();
        let seq_len: usize = fields[2].parse().unwrap_or(0);
        let hmm_len: usize = fields[5].parse().unwrap_or(0);
        let score: f64 = fields[13].parse().unwrap_or(0.0);
        let hmm_start: usize = fields[15].parse::<usize>().unwrap_or(1).saturating_sub(1);
        let hmm_stop: usize = fields[16].parse().unwrap_or(0);
        let seq_start: usize = fields[17].parse::<usize>().unwrap_or(1).saturating_sub(1);
        let seq_stop: usize = fields[18].parse().unwrap_or(0);

        if score <= 0.0 || hmm_start >= hmm_stop || seq_start >= seq_stop {
            continue;
        }

        // Look up fam by HMM accession
        let fam_id = match batch.hmm2fam.get(&query_accession) {
            Some(id) => id.clone(),
            None => continue,
        };

        let domain = HmmDomain {
            score,
            hmm_len,
            hmm_start,
            hmm_stop,
            seq_len,
            seq_start,
            seq_stop,
        };

        let key = (target_name, fam_id);
        let existing = batch.domains.get(&key);
        if existing.is_none() || existing.unwrap().score < score {
            batch.domains.insert(key, domain);
        }
    }

    Ok(())
}

/// Parse HMM search results
pub fn load_hmm_results(
    hmmsearch_file: &Path,
    batch: &mut Batch,
) -> Result<()> {
    let file = File::open(hmmsearch_file)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() < 9 {
            continue;
        }

        // hmmsearch --tblout format:
        // target_name  accession  query_name  query_accession  E-value  score  bias  E-value  score  ...
        let sseqid = fields[0].to_string();
        let hmm = fields[2].to_string();
        let hmm_accession = fields[3].to_string();
        let score1: f64 = fields[5].parse().unwrap_or(0.0);
        let score2: f64 = fields[8].parse().unwrap_or(0.0);

        // Look up fam
        let fam_id = match batch.hmm2fam.get(&hmm_accession) {
            Some(id) => id.clone(),
            None => match batch.hmm2fam.get(&hmm) {
                Some(id) => id.clone(),
                None => continue,
            },
        };

        let hmm_al = HmmAlignment {
            sseqid: sseqid.clone(),
            score1,
            score2,
            fam_id: fam_id.clone(),
            domain: None,
            blast_al_idx: None,
        };

        // Check if good
        if !hmm_al.good(&batch.fam_map) {
            continue;
        }

        // Get domain info
        let domain = batch.domains.get(&(sseqid.clone(), fam_id.clone())).cloned();

        if domain.is_none() || domain.as_ref().unwrap().hmm_len == 0 {
            continue;
        }

        // Find best BLAST alignment for this target
        let best_blast_idx = batch.target2blast_als
            .get(&sseqid)
            .and_then(|indices| {
                indices.iter()
                    .max_by_key(|&&idx| batch.blast_als[idx].hsp.nident)
                    .copied()
            });

        // Create synthetic BlastAlignment from HMM
        let domain = domain.unwrap();
        let hmm_blast_al = if let Some(blast_idx) = best_blast_idx {
            // Base on best BLAST alignment
            let mut al = batch.blast_als[blast_idx].clone();
            al.from_hmm = true;
            al.hmm_al_idx = Some(batch.hmm_als.len());
            al
        } else {
            // Stand-alone HMM hit — create minimal BlastAlignment
            let fam = batch.fam_map.get(&fam_id);
            let genesymbol = fam.map(|f| f.genesymbol.clone()).unwrap_or_default();
            let family_name = fam.map(|f| f.family_name.clone()).unwrap_or_default();

            let align_len = domain.hmm_stop - domain.hmm_start;
            let hsp = Hsp {
                sseqid: sseqid.clone(),
                qseqid: fam_id.clone(),
                slen: domain.seq_len,
                s_int: crate::seq::Interval::new(domain.seq_start, domain.seq_stop, 0),
                qlen: domain.hmm_len,
                q_int: crate::seq::Interval::new(domain.hmm_start, domain.hmm_stop, 0),
                q_prot: true,
                s_prot: true,
                a_prot: true,
                length: align_len,
                nident: (align_len as f64 * 0.7) as usize,
                ..Hsp::default()
            };

            BlastAlignment {
                hsp,
                from_hmm: true,
                ref_accession: String::new(),
                part: 1,
                parts: 1,
                fam_id: fam_id.clone(),
                gene: fam_id.clone(),
                resistance: String::new(),
                class: fam.map(|f| f.class.clone()).unwrap_or_default(),
                subclass: fam.map(|f| f.subclass.clone()).unwrap_or_default(),
                product: family_name,
                reportable: fam.map(|f| f.reportable).unwrap_or(0),
                genesymbol,
                method: "HMM".to_string(),
                complete_br: BlastRule::default(),
                partial_br: BlastRule::default(),
                cdss: Vec::new(),
                hmm_al_idx: Some(batch.hmm_als.len()),
                susceptible_idx: None,
                seq_changes: Vec::new(),
                ref_mutation: crate::alignment::AmrMutation::default(),
                fusion_ids: Vec::new(),
            }
        };

        // Update domain info on HMM alignment
        let mut hmm_al_with_domain = hmm_al;
        hmm_al_with_domain.domain = Some(domain);

        batch.add_hmm_al(hmm_al_with_domain);
        batch.add_blast_al(hmm_blast_al);
    }

    Ok(())
}

/// Configuration for amr_report
pub struct AmrReportConfig<'a> {
    pub fam_file: &'a Path,
    pub blastp_file: Option<&'a Path>,
    pub blastx_file: Option<&'a Path>,
    pub hmmsearch_file: Option<&'a Path>,
    pub hmmdom_file: Option<&'a Path>,
    pub gff_file: Option<&'a Path>,
    pub gff_type: &'a str,
    pub organism: &'a str,
    pub mutation_file: Option<&'a Path>,
    pub susceptible_file: Option<&'a Path>,
    pub coverage_min: f64,
    pub ident_min: f64,
    pub print_node: bool,
    pub report_core_only: bool,
    pub cds_exist: bool,
}

/// Run the amr_report processing pipeline
pub fn run_amr_report(config: &AmrReportConfig, out: &mut dyn Write) -> Result<()> {
    let reportable_min = if config.report_core_only { 1 } else { 0 };
    let mut batch = Batch::from_fam_file(config.fam_file, reportable_min)?;
    batch.cds_exist = config.cds_exist;

    // Load mutations and susceptible data
    if let Some(mut_file) = config.mutation_file {
        let organism_display = config.organism.replace('_', " ");
        batch.load_mutations(mut_file, &organism_display)?;
    }
    if let Some(sus_file) = config.susceptible_file {
        let organism_display = config.organism.replace('_', " ");
        batch.load_susceptible(sus_file, &organism_display)?;
    }

    // Load BLAST results
    if let Some(bp_file) = config.blastp_file {
        if bp_file.exists() {
            load_blast_results(bp_file, &mut batch, true, false)?;
        }
    }

    // Load HMM results
    if let (Some(dom_file), Some(search_file)) = (config.hmmdom_file, config.hmmsearch_file) {
        if dom_file.exists() && search_file.exists() {
            load_hmm_domains(dom_file, &mut batch)?;
            load_hmm_results(search_file, &mut batch)?;
        }
    }

    // Load BLASTX results
    if let Some(bx_file) = config.blastx_file {
        if bx_file.exists() {
            load_blast_results(bx_file, &mut batch, false, false)?;
        }
    }

    // Load GFF annotations and assign CDSs to BLAST alignments
    if let Some(gff_file) = config.gff_file {
        let gff_type = GffType::from_name(config.gff_type).unwrap_or(GffType::Genbank);
        if let Ok(annot) = Annot::from_gff(gff_file.to_str().unwrap_or(""), gff_type, false, false) {
            // Assign CDS loci to protein BLAST alignments
            for al in &mut batch.blast_als {
                if !al.hsp.s_prot || al.from_hmm {
                    continue;
                }
                if let Ok(loci) = annot.find_loci(&al.hsp.sseqid) {
                    al.cdss = loci.iter().map(|l| crate::gff::Locus {
                        line_num: l.line_num,
                        contig: l.contig.clone(),
                        start: l.start,
                        stop: l.stop,
                        strand: l.strand,
                        partial: l.partial,
                        contig_len: l.contig_len,
                        cross_origin: l.cross_origin,
                        gene: l.gene.clone(),
                        product: l.product.clone(),
                    }).collect();
                }
            }
        }
    }

    // Process
    batch.process();

    // Report
    batch.report(out, config.print_node)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn test_fixtures() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
    }

    fn db_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amrfinder_db/2026-03-24.1")
    }

    #[test]
    fn test_load_blast_results() {
        let blastp_file = test_fixtures().join("blastp");
        let db = db_dir();
        if !blastp_file.exists() || !db.exists() {
            return;
        }

        let mut batch = Batch::from_fam_file(&db.join("fam.tsv"), 0).unwrap();
        load_blast_results(&blastp_file, &mut batch, true, false).unwrap();
        assert!(!batch.blast_als.is_empty(), "Should have loaded BLAST alignments");
        assert!(!batch.target2blast_als.is_empty(), "Should have indexed by target");
    }

    #[test]
    fn test_load_hmm_results() {
        let fixtures = test_fixtures();
        let db = db_dir();
        let dom_file = fixtures.join("dom");
        let search_file = fixtures.join("hmmsearch");
        if !dom_file.exists() || !search_file.exists() || !db.exists() {
            return;
        }

        let mut batch = Batch::from_fam_file(&db.join("fam.tsv"), 0).unwrap();
        // Need to load BLAST first for HMM-BLAST linking
        let blastp_file = fixtures.join("blastp");
        if blastp_file.exists() {
            load_blast_results(&blastp_file, &mut batch, true, false).unwrap();
        }
        load_hmm_domains(&dom_file, &mut batch).unwrap();
        load_hmm_results(&search_file, &mut batch).unwrap();
        assert!(!batch.hmm_als.is_empty(), "Should have loaded HMM alignments");
    }

    #[test]
    fn test_run_amr_report_output() {
        let fixtures = test_fixtures();
        let db = db_dir();
        let blastp_file = fixtures.join("blastp");
        let dom_file = fixtures.join("dom");
        let search_file = fixtures.join("hmmsearch");
        let expected_file = fixtures.join("expected_output");
        let gff_file = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amr/test_prot.gff");

        if !blastp_file.exists() || !db.exists() || !expected_file.exists() {
            return;
        }

        let mut output = Vec::new();
        let config = AmrReportConfig {
            fam_file: &db.join("fam.tsv"),
            blastp_file: Some(blastp_file.as_path()),
            blastx_file: None,
            hmmsearch_file: Some(search_file.as_path()),
            hmmdom_file: Some(dom_file.as_path()),
            gff_file: Some(&gff_file),
            gff_type: "genbank",
            organism: "Escherichia",
            mutation_file: Some(&db.join("AMRProt-mutation.tsv")),
            susceptible_file: Some(&db.join("AMRProt-susceptible.tsv")),
            coverage_min: 0.5,
            ident_min: -1.0,
            print_node: true,
            report_core_only: false,
            cds_exist: true,
        };
        run_amr_report(&config, &mut output).unwrap();

        let rust_output = String::from_utf8(output).unwrap();
        let cpp_output = std::fs::read_to_string(&expected_file).unwrap();

        // Compare line counts first
        let rust_lines: Vec<&str> = rust_output.lines().collect();
        let cpp_lines: Vec<&str> = cpp_output.lines().collect();

        eprintln!("Rust output: {} lines", rust_lines.len());
        eprintln!("C++ output: {} lines", cpp_lines.len());

        assert!(!rust_output.is_empty(), "Rust amr_report should produce output");

        // Save for debugging
        std::fs::write("/tmp/rust_report_debug.tsv", &rust_output).ok();

        // Check header matches
        if !rust_lines.is_empty() && !cpp_lines.is_empty() {
            assert_eq!(rust_lines[0], cpp_lines[0], "Headers should match");
        }

        // Count targets in Rust output
        let mut target_counts: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
        for line in &rust_lines[1..] {
            let target = line.split('\t').next().unwrap_or("");
            *target_counts.entry(target).or_insert(0) += 1;
        }
        eprintln!("Targets with >1 hit:");
        for (t, c) in &target_counts {
            if *c > 1 {
                eprintln!("  {}: {} hits", t, c);
            }
        }
    }
}
