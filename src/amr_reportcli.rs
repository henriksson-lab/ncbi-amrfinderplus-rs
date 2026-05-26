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

    let default_complete_br = BlastRule::new(0.9, 0.0, 0.9);
    let default_partial_br = BlastRule::new(0.9, 0.0, 0.5);

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        let mut al = BlastAlignment::from_blast_line(
            trimmed,
            true,       // q_prot (reference is always protein)
            is_protein, // s_prot (target is protein for blastp, DNA for blastx)
            &default_complete_br,
            &default_partial_br,
        )?;

        if let Some(fam) = batch.fam_map.get(&al.fam_id) {
            if !fam.complete_br.empty() {
                al.complete_br = fam.complete_br.clone();
                al.partial_br = fam.partial_br.clone();
            }
        }

        if nosame && al.ref_accession == al.hsp.sseqid {
            continue;
        }

        batch.add_blast_al(al);
    }

    Ok(())
}

/// Parse HMM domain table output
pub fn load_hmm_domains(dom_file: &Path, batch: &mut Batch) -> Result<()> {
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
pub fn load_hmm_results(hmmsearch_file: &Path, batch: &mut Batch) -> Result<()> {
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
        let domain = batch
            .domains
            .get(&(sseqid.clone(), fam_id.clone()))
            .cloned();

        if domain.is_none() || domain.as_ref().unwrap().hmm_len == 0 {
            continue;
        }

        // Find best BLAST alignment for this target
        let best_blast_idx = batch.target2blast_als.get(&sseqid).and_then(|indices| {
            indices
                .iter()
                .max_by_key(|&&idx| batch.blast_als[idx].hsp.nident)
                .copied()
        });

        // Create synthetic BlastAlignment from HMM
        let domain = domain.unwrap();
        // Update domain info on HMM alignment
        let mut hmm_al_with_domain = hmm_al;
        hmm_al_with_domain.domain = Some(domain.clone());
        let hmm_idx = batch.hmm_als.len();
        batch.add_hmm_al(hmm_al_with_domain);

        if let Some(blast_idx) = best_blast_idx {
            if batch.blast_als[blast_idx].good() {
                let replace = batch.blast_als[blast_idx]
                    .hmm_al_idx
                    .map(|old_idx| batch.hmm_als[hmm_idx].score1 > batch.hmm_als[old_idx].score1)
                    .unwrap_or(true);
                if replace {
                    batch.blast_als[blast_idx].hmm_al_idx = Some(hmm_idx);
                }
                continue;
            }

            let fam = batch.fam_map.get(&fam_id);
            let mut al = batch.blast_als[blast_idx].clone();
            al.from_hmm = true;
            al.hmm_al_idx = Some(hmm_idx);
            al.fam_id = fam_id.clone();
            al.gene = fam_id.clone();
            al.resistance.clear();
            al.class = fam.map(|f| f.class.clone()).unwrap_or_default();
            al.subclass = fam.map(|f| f.subclass.clone()).unwrap_or_default();
            al.product = fam.map(|f| f.family_name.clone()).unwrap_or_default();
            al.reportable = fam.map(|f| f.reportable).unwrap_or(0);
            al.genesymbol = fam.map(|f| f.genesymbol.clone()).unwrap_or_default();
            batch.add_blast_al(al);
            continue;
        }

        let hmm_blast_al = {
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
    pub suppress_file: Option<&'a Path>,
    pub coverage_min: f64,
    pub ident_min: f64,
    pub print_node: bool,
    pub mutation_all: Option<&'a Path>,
    pub name: &'a str,
    pub report_core_only: bool,
    pub report_all_equal: bool,
    pub cds_exist: bool,
}

/// Run the amr_report processing pipeline
pub fn run_amr_report(config: &AmrReportConfig, out: &mut dyn Write) -> Result<()> {
    let reportable_min = if config.report_core_only { 2 } else { 0 };
    let mut batch = Batch::from_fam_file(config.fam_file, reportable_min)?;
    batch.cds_exist = config.cds_exist;
    batch.name = config.name.to_string();
    batch.ident_min = config.ident_min;
    batch.coverage_min = config.coverage_min;
    batch.report_all_equal = config.report_all_equal;

    // Load mutations and susceptible data
    if let Some(mut_file) = config.mutation_file {
        let organism_display = config.organism.replace('_', " ");
        batch.load_mutations(mut_file, &organism_display)?;
    }
    if let Some(sus_file) = config.susceptible_file {
        let organism_display = config.organism.replace('_', " ");
        batch.load_susceptible(sus_file, &organism_display)?;
    }
    if let Some(suppress_file) = config.suppress_file {
        batch.load_suppress(suppress_file, config.organism)?;
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
        let gff_type = GffType::from_name(config.gff_type)?;
        let annot = Annot::from_gff(gff_file.to_str().unwrap_or(""), gff_type, false, false)?;
        // Assign CDS loci to protein BLAST alignments
        for al in &mut batch.blast_als {
            if !al.hsp.s_prot {
                continue;
            }
            if let Ok(loci) = annot.find_loci(&al.hsp.sseqid) {
                al.cdss = loci
                    .iter()
                    .map(|l| crate::gff::Locus {
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
                    })
                    .collect();
            }
        }
    }

    // Process
    batch.process();

    if let Some(path) = config.mutation_all {
        let mut mutation_all = File::create(path)?;
        batch.report_mutation_all(&mut mutation_all, config.print_node)?;
    }

    // Report
    batch.report(out, config.print_node)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::{Path, PathBuf};

    fn test_fixtures() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
    }

    fn db_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amrfinder_db/2026-03-24.1")
    }

    fn require_file(path: &Path) {
        assert!(
            path.exists(),
            "required parity fixture is missing: {}",
            path.display()
        );
    }

    fn run_cpp_golden_report() -> (String, String) {
        let fixtures = test_fixtures();
        let db = db_dir();
        let blastp_file = fixtures.join("blastp");
        let dom_file = fixtures.join("dom");
        let search_file = fixtures.join("hmmsearch");
        let expected_file = fixtures.join("expected_output");
        let gff_file = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amr/test_prot.gff");

        require_file(&db.join("fam.tsv"));
        require_file(&db.join("AMRProt-mutation.tsv"));
        require_file(&db.join("AMRProt-susceptible.tsv"));
        require_file(&blastp_file);
        require_file(&dom_file);
        require_file(&search_file);
        require_file(&expected_file);
        require_file(&gff_file);

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
            suppress_file: None,
            coverage_min: 0.5,
            ident_min: -1.0,
            print_node: true,
            mutation_all: None,
            name: "",
            report_core_only: false,
            report_all_equal: false,
            cds_exist: true,
        };
        run_amr_report(&config, &mut output).unwrap();

        let rust_output = String::from_utf8(output).unwrap();
        let cpp_output = std::fs::read_to_string(&expected_file).unwrap();
        (rust_output, cpp_output)
    }

    fn run_cpp_blastx_golden_report() -> (String, String) {
        let fixtures = test_fixtures();
        let db = db_dir();
        let blastx_file = fixtures.join("blastx");
        let expected_file = fixtures.join("expected_blastx_report");

        require_file(&db.join("fam.tsv"));
        require_file(&db.join("AMRProt-mutation.tsv"));
        require_file(&db.join("AMRProt-susceptible.tsv"));
        require_file(&blastx_file);
        require_file(&expected_file);

        let mut output = Vec::new();
        let config = AmrReportConfig {
            fam_file: &db.join("fam.tsv"),
            blastp_file: None,
            blastx_file: Some(blastx_file.as_path()),
            hmmsearch_file: None,
            hmmdom_file: None,
            gff_file: None,
            gff_type: "genbank",
            organism: "Escherichia",
            mutation_file: Some(&db.join("AMRProt-mutation.tsv")),
            susceptible_file: Some(&db.join("AMRProt-susceptible.tsv")),
            suppress_file: None,
            coverage_min: 0.5,
            ident_min: -1.0,
            print_node: true,
            mutation_all: None,
            name: "",
            report_core_only: false,
            report_all_equal: false,
            cds_exist: true,
        };
        run_amr_report(&config, &mut output).unwrap();

        let rust_output = String::from_utf8(output).unwrap();
        let cpp_output = std::fs::read_to_string(&expected_file).unwrap();
        (rust_output, cpp_output)
    }

    fn run_combined_golden_report() -> String {
        let fixtures = test_fixtures();
        let db = db_dir();
        let blastp_file = fixtures.join("blastp");
        let blastx_file = fixtures.join("blastx");
        let dom_file = fixtures.join("dom");
        let search_file = fixtures.join("hmmsearch");
        let gff_file = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amr/test_prot.gff");

        require_file(&db.join("fam.tsv"));
        require_file(&db.join("AMRProt-mutation.tsv"));
        require_file(&db.join("AMRProt-susceptible.tsv"));
        require_file(&blastp_file);
        require_file(&blastx_file);
        require_file(&dom_file);
        require_file(&search_file);
        require_file(&gff_file);

        let mut output = Vec::new();
        let config = AmrReportConfig {
            fam_file: &db.join("fam.tsv"),
            blastp_file: Some(blastp_file.as_path()),
            blastx_file: Some(blastx_file.as_path()),
            hmmsearch_file: Some(search_file.as_path()),
            hmmdom_file: Some(dom_file.as_path()),
            gff_file: Some(&gff_file),
            gff_type: "genbank",
            organism: "Escherichia",
            mutation_file: Some(&db.join("AMRProt-mutation.tsv")),
            susceptible_file: Some(&db.join("AMRProt-susceptible.tsv")),
            suppress_file: None,
            coverage_min: 0.5,
            ident_min: -1.0,
            print_node: true,
            mutation_all: None,
            name: "",
            report_core_only: false,
            report_all_equal: false,
            cds_exist: true,
        };
        run_amr_report(&config, &mut output).unwrap();

        String::from_utf8(output).unwrap()
    }

    fn run_protein_report_for_organism(organism: &str, suppress: bool) -> String {
        let fixtures = test_fixtures();
        let db = db_dir();
        let blastp_file = fixtures.join("blastp");
        let dom_file = fixtures.join("dom");
        let search_file = fixtures.join("hmmsearch");
        let gff_file = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amr/test_prot.gff");
        let suppress_file = db.join("AMRProt-suppress.tsv");

        require_file(&db.join("fam.tsv"));
        require_file(&db.join("AMRProt-mutation.tsv"));
        require_file(&db.join("AMRProt-susceptible.tsv"));
        require_file(&blastp_file);
        require_file(&dom_file);
        require_file(&search_file);
        require_file(&gff_file);
        if suppress {
            require_file(&suppress_file);
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
            organism,
            mutation_file: Some(&db.join("AMRProt-mutation.tsv")),
            susceptible_file: Some(&db.join("AMRProt-susceptible.tsv")),
            suppress_file: suppress.then_some(suppress_file.as_path()),
            coverage_min: 0.5,
            ident_min: -1.0,
            print_node: true,
            mutation_all: None,
            name: "",
            report_core_only: false,
            report_all_equal: false,
            cds_exist: true,
        };
        run_amr_report(&config, &mut output).unwrap();

        String::from_utf8(output).unwrap()
    }

    fn tsv_row<'a>(output: &'a str, protein_id: &str) -> Vec<&'a str> {
        output
            .lines()
            .skip(1)
            .find(|line| line.split('\t').next() == Some(protein_id))
            .unwrap_or_else(|| panic!("missing report row for {protein_id}"))
            .split('\t')
            .collect()
    }

    fn col<'a>(header: &[&str], row: &'a [&str], name: &str) -> &'a str {
        let idx = header
            .iter()
            .position(|field| *field == name)
            .unwrap_or_else(|| panic!("missing report column {name}"));
        row[idx]
    }

    #[test]
    fn amr_report_cpp_golden_fixture_is_complete() {
        let (_rust_output, cpp_output) = run_cpp_golden_report();
        let cpp_lines: Vec<&str> = cpp_output.lines().collect();

        assert_eq!(
            cpp_lines.len(),
            16,
            "C++ golden should contain header plus 15 rows"
        );
        assert!(
            cpp_output.contains("nimIJ_hmm\tcontigX"),
            "C++ golden should include the standalone HMM-only hit"
        );
        assert!(
            cpp_output.contains("pmrB_C84R\tcontig14"),
            "C++ golden should include protein mutation hits"
        );
    }

    #[test]
    fn amr_report_matches_cpp_golden_byte_for_byte() {
        let (rust_output, cpp_output) = run_cpp_golden_report();

        assert_eq!(
            rust_output, cpp_output,
            "Rust amr_report output must be byte-identical to the C++ golden output"
        );
    }

    #[test]
    fn amr_report_matches_cpp_blastx_golden_byte_for_byte() {
        let (rust_output, cpp_output) = run_cpp_blastx_golden_report();

        assert_eq!(
            rust_output, cpp_output,
            "Rust BLASTX amr_report output must be byte-identical to the C++ golden output"
        );
    }

    #[test]
    fn hmm_loader_links_or_synthesizes_expected_alignment_kinds() {
        let fixtures = test_fixtures();
        let db = db_dir();
        let mut batch = Batch::from_fam_file(&db.join("fam.tsv"), 0).unwrap();

        load_blast_results(&fixtures.join("blastp"), &mut batch, true, false).unwrap();
        let blast_count_before_hmm = batch.blast_als.len();
        load_hmm_domains(&fixtures.join("dom"), &mut batch).unwrap();
        load_hmm_results(&fixtures.join("hmmsearch"), &mut batch).unwrap();

        let bla_tem_rows: Vec<_> = batch
            .blast_als
            .iter()
            .filter(|al| al.hsp.sseqid == "blaTEM-156")
            .collect();
        assert!(
            bla_tem_rows
                .iter()
                .any(|al| al.ref_accession == "WP_061158039.1" && al.hmm_al_idx.is_some()),
            "good BLAST hit should link to its HMM alignment"
        );
        assert!(
            bla_tem_rows.iter().all(|al| !al.from_hmm),
            "good BLAST target should not get duplicate synthetic HMM rows"
        );

        let nimij_hmm = batch
            .blast_als
            .iter()
            .find(|al| al.hsp.sseqid == "nimIJ_hmm" && al.from_hmm)
            .expect("below-threshold BLAST target should get a synthetic HMM row");
        assert_eq!(nimij_hmm.fam_id, "nimIJ");
        assert_eq!(nimij_hmm.ref_accession, "WP_005812825.1");
        assert!(nimij_hmm.hmm_al_idx.is_some());
        assert!(
            batch.blast_als.len() > blast_count_before_hmm,
            "HMM processing should add synthetic rows only where BLAST alone is insufficient"
        );
    }

    #[test]
    fn hmm_loader_creates_standalone_row_without_blast_target() {
        let db = db_dir();
        let mut batch = Batch::from_fam_file(&db.join("fam.tsv"), 0).unwrap();
        let mut dom = tempfile::NamedTempFile::new().unwrap();
        let mut search = tempfile::NamedTempFile::new().unwrap();

        writeln!(
            dom,
            "standalone_hmm - 166 NimIJ-NCBIFAM NF000262.1 151 1.1e-82 274.7 0.3 1 1 1.2e-86 1.2e-82 274.5 0.3 1 151 4 154 4 154 1.00 description"
        )
        .unwrap();
        writeln!(
            search,
            "standalone_hmm - NimIJ-NCBIFAM NF000262.1 1.1e-82 274.7 0.3 1.2e-82 274.5 0.3 1.0 1 0 0 1 1 1 1 description"
        )
        .unwrap();

        load_hmm_domains(dom.path(), &mut batch).unwrap();
        load_hmm_results(search.path(), &mut batch).unwrap();

        assert_eq!(batch.hmm_als.len(), 1);
        assert_eq!(batch.blast_als.len(), 1);
        let al = &batch.blast_als[0];
        assert!(al.from_hmm);
        assert_eq!(al.hsp.sseqid, "standalone_hmm");
        assert_eq!(al.ref_accession, "");
        assert_eq!(al.fam_id, "nimIJ");
        assert_eq!(al.method, "HMM");
    }

    #[test]
    fn report_rows_cover_known_parity_edge_cases() {
        let (rust_output, _cpp_output) = run_cpp_golden_report();
        let header: Vec<&str> = rust_output.lines().next().unwrap().split('\t').collect();

        let aph = tsv_row(&rust_output, "aph3pp-Ib_partial_5p_neg");
        assert_eq!(col(&header, &aph, "% Identity to reference"), "100.00");
        assert_eq!(col(&header, &aph, "Method"), "PARTIAL_CONTIG_ENDP");

        let mutation = tsv_row(&rust_output, "pmrB_C84R");
        assert_eq!(col(&header, &mutation, "Element symbol"), "pmrB_C84R");
        assert_eq!(col(&header, &mutation, "Type"), "AMR");
        assert_eq!(col(&header, &mutation, "Subtype"), "POINT");
        assert_eq!(col(&header, &mutation, "Hierarchy node"), "pmrB");

        let oxa = tsv_row(&rust_output, "blaOXA-436_partial");
        assert_eq!(col(&header, &oxa, "Hierarchy node"), "blaOXA-48_fam");

        let sul2 = tsv_row(&rust_output, "sul2_partial_3p_neg");
        assert_eq!(col(&header, &sul2, "HMM accession"), "NA");
        assert_eq!(col(&header, &sul2, "HMM description"), "NA");
    }

    #[test]
    fn blastx_report_uses_dna_coordinates_and_x_methods() {
        let db = db_dir();
        let mut blastx = tempfile::NamedTempFile::new().unwrap();
        let qseq = format!("M{}", "A".repeat(285));
        let sseq = qseq.clone();
        writeln!(
            blastx,
            "WP_061158039.1|1|1|blaTEM-156|blaTEM|beta-lactamase|2|BETA-LACTAM|BETA-LACTAM|class_A_beta-lactamase_TEM-156\tcontig_blastx\t1\t286\t287\t101\t958\t1200\t{}\t{}",
            qseq, sseq
        )
        .unwrap();

        let mut output = Vec::new();
        let config = AmrReportConfig {
            fam_file: &db.join("fam.tsv"),
            blastp_file: None,
            blastx_file: Some(blastx.path()),
            hmmsearch_file: None,
            hmmdom_file: None,
            gff_file: None,
            gff_type: "genbank",
            organism: "Escherichia",
            mutation_file: Some(&db.join("AMRProt-mutation.tsv")),
            susceptible_file: Some(&db.join("AMRProt-susceptible.tsv")),
            suppress_file: None,
            coverage_min: 0.5,
            ident_min: -1.0,
            print_node: true,
            mutation_all: None,
            name: "",
            report_core_only: false,
            report_all_equal: false,
            cds_exist: true,
        };
        run_amr_report(&config, &mut output).unwrap();

        let output = String::from_utf8(output).unwrap();
        let header: Vec<&str> = output.lines().next().unwrap().split('\t').collect();
        let row = tsv_row(&output, "NA");
        assert_eq!(col(&header, &row, "Contig id"), "contig_blastx");
        assert_eq!(col(&header, &row, "Start"), "101");
        assert_eq!(col(&header, &row, "Stop"), "958");
        assert_eq!(col(&header, &row, "Method"), "ALLELEX");
        assert_eq!(col(&header, &row, "Target length"), "286");
        assert_eq!(col(&header, &row, "Reference sequence length"), "286");
        assert_eq!(col(&header, &row, "% Coverage of reference"), "100.00");
        assert_eq!(col(&header, &row, "% Identity to reference"), "100.00");
    }

    #[test]
    fn blastx_mutation_all_reports_wildtype_and_mutant_rows() {
        let fixtures = test_fixtures();
        let db = db_dir();
        let blastx_file = fixtures.join("blastx");

        if !blastx_file.exists() {
            return;
        }
        require_file(&db.join("fam.tsv"));
        require_file(&db.join("AMRProt-mutation.tsv"));
        require_file(&db.join("AMRProt-susceptible.tsv"));

        let mutation_all = tempfile::NamedTempFile::new().unwrap();
        let mut output = Vec::new();
        let config = AmrReportConfig {
            fam_file: &db.join("fam.tsv"),
            blastp_file: None,
            blastx_file: Some(blastx_file.as_path()),
            hmmsearch_file: None,
            hmmdom_file: None,
            gff_file: None,
            gff_type: "genbank",
            organism: "Escherichia",
            mutation_file: Some(&db.join("AMRProt-mutation.tsv")),
            susceptible_file: Some(&db.join("AMRProt-susceptible.tsv")),
            suppress_file: None,
            coverage_min: 0.5,
            ident_min: -1.0,
            print_node: true,
            mutation_all: Some(mutation_all.path()),
            name: "",
            report_core_only: false,
            report_all_equal: false,
            cds_exist: true,
        };
        run_amr_report(&config, &mut output).unwrap();

        let mutation_all = std::fs::read_to_string(mutation_all.path()).unwrap();
        assert!(
            mutation_all.contains("\nNA\tcontig14\t1\t1089\t+\tpmrB_C84R\t"),
            "BLASTX mutation_all should include mutant pmrB row"
        );
        assert!(
            mutation_all.contains("\nNA\tcontig14\t1\t1089\t+\tpmrB_A159A\t"),
            "BLASTX mutation_all should include wildtype pmrB rows"
        );
        assert!(
            mutation_all.contains("\nNA\tcontig16\t1\t720\t+\tnfsA_K141Ter\t")
                && mutation_all.contains("\nNA\tcontig16\t1\t720\t+\tnfsA_R15C\t"),
            "BLASTX mutation_all should include all observed nfsA mutation rows"
        );
        assert!(
            mutation_all
                .lines()
                .next()
                .is_some_and(|header| header.ends_with("\tHierarchy node")),
            "direct amr_report mutation_all should preserve --print_node"
        );
    }

    #[test]
    fn combined_report_suppresses_blastx_when_protein_cds_covers_region() {
        let output = run_combined_golden_report();

        assert!(
            output.contains("\nblaTEM-156\tcontig01\t101\t961\t+\tblaTEM-156\t"),
            "protein/CDS row should be retained"
        );
        assert!(
            !output.contains("\nNA\tcontig01\t101\t958\t+\tblaTEM-156\t"),
            "overlapping BLASTX row should be suppressed when a protein CDS row covers it"
        );
        assert!(
            output.contains("\nNA\tcontig04\t1261\t2391\t+\tblaEC\t"),
            "non-overlapping BLASTX hit on the same contig should remain reportable"
        );
    }

    #[test]
    fn combined_report_transfers_internal_stop_from_blastx_to_protein_row() {
        let output = run_combined_golden_report();
        let header: Vec<&str> = output.lines().next().unwrap().split('\t').collect();
        let row = tsv_row(&output, "blaTEM-internal_stop");

        assert_eq!(col(&header, &row, "Method"), "INTERNAL_STOP");
        assert!(
            !output.contains("\nNA\tcontig11\t101\t958\t+\tblaTEM\t"),
            "BLASTX row carrying the same internal-stop evidence should be suppressed"
        );
    }

    #[test]
    fn combined_report_prefers_full_blastx_mutations_over_partial_protein_mutation() {
        let output = run_combined_golden_report();

        assert!(
            !output.contains("\nnfsA_R15C_K141STOP\t"),
            "partial protein mutation row should be suppressed by fuller BLASTX evidence"
        );
        assert!(
            output.contains("\nNA\tcontig16\t1\t720\t+\tnfsA_K141Ter\t")
                && output.contains("\nNA\tcontig16\t1\t720\t+\tnfsA_R15C\t"),
            "distinct BLASTX mutation rows should remain reportable"
        );
    }

    #[test]
    fn organism_suppress_file_removes_common_proteins_unless_report_common_is_used() {
        let suppressed = run_protein_report_for_organism("Escherichia", true);
        let report_common = run_protein_report_for_organism("Escherichia", false);

        assert!(
            !suppressed.contains("\narsR-suppressed-in-escherichia\t"),
            "Escherichia suppress file should remove arsR"
        );
        assert!(
            report_common.contains("\narsR-suppressed-in-escherichia\t"),
            "report_common path should omit suppress file and retain arsR"
        );
    }

    #[test]
    fn campylobacter_fixture_reports_campylobacter_specific_mutations() {
        let output = run_protein_report_for_organism("Campylobacter", false);

        assert!(
            output.contains("\ngyrA\tcontig06\t31\t2616\t+\tgyrA_T86A\t")
                || output.contains("\tgyrA_T86A\tCampylobacter quinolone resistant GyrA\t"),
            "Campylobacter organism should report the gyrA_T86A protein mutation"
        );
        assert!(
            output.contains("\t50S_L22:A103V\tCampylobacter macrolide resistant RplV\t")
                || output.contains("\trplV_A103V\tCampylobacter macrolide resistant RplV\t"),
            "Campylobacter organism should report the 50S_L22/rplV mutation fixture"
        );
        assert!(
            !output.contains("\tnfsA_R15C\tEscherichia nitrofurantoin resistant NfsA\t"),
            "Escherichia-specific nfsA mutation should not be reported for Campylobacter"
        );
    }
}
