use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};

use anyhow::{bail, Result};

/// GFF file type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GffType {
    Bakta,
    Genbank,
    Microscope,
    Patric,
    Pgap,
    Prodigal,
    Prokka,
    Pseudomonasdb,
    Rast,
    Standard,
}

impl GffType {
    pub const NAMES: &[&str] = &[
        "bakta", "genbank", "microscope", "patric", "pgap",
        "prodigal", "prokka", "pseudomonasdb", "rast", "standard",
    ];

    pub fn from_name(name: &str) -> Result<GffType> {
        match name {
            "bakta" => Ok(GffType::Bakta),
            "genbank" => Ok(GffType::Genbank),
            "microscope" => Ok(GffType::Microscope),
            "patric" => Ok(GffType::Patric),
            "pgap" => Ok(GffType::Pgap),
            "prodigal" => Ok(GffType::Prodigal),
            "prokka" => Ok(GffType::Prokka),
            "pseudomonasdb" => Ok(GffType::Pseudomonasdb),
            "rast" => Ok(GffType::Rast),
            "standard" => Ok(GffType::Standard),
            _ => bail!("Unknown GFF type: \"{}\"", name),
        }
    }
}

/// Genomic locus
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Locus {
    pub line_num: usize,
    pub contig: String,
    pub start: usize,
    pub stop: usize,
    pub strand: bool,
    pub partial: bool,
    pub contig_len: usize,
    pub cross_origin: bool,
    pub gene: String,
    pub product: String,
}

impl Locus {
    const END_DELTA: usize = 3;

    #[allow(clippy::too_many_arguments)]
    pub fn new(
        line_num: usize,
        contig: &str,
        start: usize,
        stop: usize,
        strand: bool,
        partial: bool,
        cross_origin_seq_len: usize,
        gene: String,
        product: String,
    ) -> Result<Self> {
        let contig = contig.trim().to_string();
        if contig.is_empty() {
            bail!("Empty contig name");
        }

        let (start, stop, cross_origin) = if cross_origin_seq_len > 0 {
            // swap and adjust for cross-origin
            let new_start = stop - 1;
            let new_stop = start + 1;
            assert!(cross_origin_seq_len > 0);
            assert!(new_stop <= cross_origin_seq_len);
            (new_start, new_stop, true)
        } else {
            (start, stop, false)
        };

        assert!(start < stop);

        Ok(Locus {
            line_num,
            contig,
            start,
            stop,
            strand,
            partial,
            contig_len: cross_origin_seq_len,
            cross_origin,
            gene,
            product,
        })
    }

    pub fn empty(&self) -> bool {
        self.contig.is_empty()
    }

    pub fn size(&self) -> usize {
        if self.cross_origin {
            self.contig_len - self.stop + self.start
        } else {
            self.stop - self.start
        }
    }

    pub fn at_contig_start(&self) -> bool {
        self.start <= Self::END_DELTA
    }

    pub fn at_contig_stop(&self) -> bool {
        self.contig_len > 0 && self.contig_len - self.stop <= Self::END_DELTA
    }
}

impl Ord for Locus {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.contig
            .cmp(&other.contig)
            .then(self.start.cmp(&other.start))
            .then(self.stop.cmp(&other.stop))
            .then(self.strand.cmp(&other.strand))
            .then(self.cross_origin.cmp(&other.cross_origin))
    }
}

impl PartialOrd for Locus {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// URL percent-decode and trim
fn unescape(s: &str) -> String {
    let decoded = percent_decode(s);
    decoded.trim().to_string()
}

fn percent_decode(s: &str) -> String {
    let mut result = String::with_capacity(s.len());
    let bytes = s.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%' && i + 2 < bytes.len() {
            if let Ok(val) = u8::from_str_radix(
                std::str::from_utf8(&bytes[i + 1..i + 3]).unwrap_or(""),
                16,
            ) {
                result.push(val as char);
                i += 3;
                continue;
            }
        }
        result.push(bytes[i] as char);
        i += 1;
    }
    result
}

fn pgap_accession(accession: &mut String, lcl: bool) -> Result<()> {
    let gnl_prefix = "gnl|";
    let lcl_prefix = "lcl|";

    if let Some(pos) = accession.rfind(':') {
        if lcl {
            bail!(
                "Accession \"{}\" cannot have \"{}\" and \"{}\" at the same time",
                accession, gnl_prefix, lcl_prefix
            );
        }
        // Replace ':' with '|' at the position
        let mut chars: Vec<char> = accession.chars().collect();
        chars[pos] = '|';
        *accession = format!("{}{}", gnl_prefix, chars.into_iter().collect::<String>());
    } else if lcl {
        *accession = format!("{}{}", lcl_prefix, accession);
    }

    assert!(!accession.is_empty());
    Ok(())
}

fn printable_char(c: char) -> bool {
    let code = c as u32;
    (32..127).contains(&code)
}

/// GFF/BED annotation data
pub struct Annot {
    pub prot2loci: BTreeMap<String, BTreeSet<Locus>>,
    pub fasta2gff_prot: HashMap<String, String>,
}

impl Annot {
    /// Parse a GFF file
    pub fn from_gff(
        fname: &str,
        gff_type: GffType,
        prot_match: bool,
        lcl: bool,
    ) -> Result<Self> {
        assert!(
            !prot_match
                || gff_type == GffType::Genbank
                || gff_type == GffType::Microscope
                || gff_type == GffType::Prodigal
        );
        if gff_type == GffType::Microscope {
            assert!(prot_match);
        }
        if lcl {
            assert!(gff_type == GffType::Pgap);
        }

        if fname.is_empty() {
            bail!("Empty GFF file name");
        }

        let mut prot2loci: BTreeMap<String, BTreeSet<Locus>> = BTreeMap::new();

        let file = File::open(fname)?;
        let reader = BufReader::new(file);

        for (line_num, line_result) in reader.lines().enumerate() {
            let mut line = line_result?;
            let trimmed = line.trim().to_string();
            line = trimmed;

            if (gff_type == GffType::Prokka || gff_type == GffType::Bakta)
                && line == "##FASTA"
            {
                break;
            }

            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let parse_result = (|| -> Result<()> {
                let mut rest = line.as_str();

                let mut contig = unescape(split_tab(&mut rest)?);
                let _source = unescape(split_tab(&mut rest)?);
                let feat_type = unescape(split_tab(&mut rest)?);
                let start_s = unescape(split_tab(&mut rest)?);
                let stop_s = unescape(split_tab(&mut rest)?);
                let _score = unescape(split_tab(&mut rest)?);
                let strand_s = unescape(split_tab(&mut rest)?);
                let _phase = unescape(split_tab(&mut rest)?);
                let attributes = rest.trim().to_string();

                if attributes.is_empty() {
                    bail!("9 fields are expected in each line");
                }

                if contig.is_empty() {
                    bail!("empty sequence indentifier");
                }
                for c in contig.chars() {
                    if !printable_char(c) {
                        bail!(
                            "Non-printable character in the sequence identifier: {}",
                            c as u32
                        );
                    }
                }

                if feat_type != "CDS" && feat_type != "gene" && feat_type != "pseudogene" {
                    return Ok(());
                }

                if gff_type == GffType::Pgap && feat_type != "CDS" {
                    return Ok(());
                }

                let start: i64 = start_s
                    .parse()
                    .map_err(|_| anyhow::anyhow!("Cannot read start"))?;
                if start <= 0 {
                    bail!("start should be >= 1");
                }

                let stop: i64 = stop_s
                    .parse()
                    .map_err(|_| anyhow::anyhow!("Cannot read stop"))?;
                if stop <= 0 {
                    bail!("stop should be >= 1");
                }

                if start > stop {
                    bail!("start cannot be greater than stop");
                }

                let start = (start - 1) as usize;
                let stop = stop as usize;

                if strand_s != "+" && strand_s != "-" {
                    bail!("strand should be '+' or '-'");
                }

                let pseudo = attributes.contains("pseudo=true")
                    || attributes.contains("gene_biotype=pseudogene")
                    || feat_type == "pseudogene";

                let partial = attributes.contains("partial=true")
                    || attributes.contains("partial=01")
                    || attributes.contains("partial=10")
                    || attributes.contains("partial=11");

                let prot_attr = match gff_type {
                    GffType::Bakta => "ID",
                    GffType::Genbank => {
                        if prot_match || pseudo {
                            "locus_tag"
                        } else {
                            "Name"
                        }
                    }
                    GffType::Microscope => "ID",
                    GffType::Patric => "ID",
                    GffType::Prodigal => "ID",
                    GffType::Prokka => "ID",
                    GffType::Pseudomonasdb => "Alias",
                    GffType::Rast => "ID",
                    _ => "Name",
                };
                let prot_attr_eq = format!("{}=", prot_attr);

                let mut prot_ = String::new();
                let mut gene_ = String::new();
                let mut product_ = String::new();
                let mut locus_tag = String::new();

                for attr_part in attributes.split(';') {
                    let s = attr_part.trim();
                    if let Some(val) = s.strip_prefix(&*prot_attr_eq) {
                        prot_ = val.to_string();
                    } else if let Some(val) = s.strip_prefix("gene=") {
                        gene_ = val.to_string();
                    } else if let Some(val) = s.strip_prefix("product=") {
                        product_ = val.to_string();
                    } else if gff_type == GffType::Patric {
                        if let Some(val) = s.strip_prefix("locus_tag=") {
                            locus_tag = val.to_string();
                        }
                    }
                }

                // Trim surrounding quotes
                if prot_.starts_with('"') {
                    prot_ = prot_.trim_start_matches('"').to_string();
                }
                if prot_.ends_with('"') {
                    prot_ = prot_.trim_end_matches('"').to_string();
                }

                if prot_.is_empty() {
                    return Ok(());
                }

                match gff_type {
                    GffType::Genbank => {
                        if let Some(pos) = prot_.find(':') {
                            prot_ = prot_[pos + 1..].to_string();
                        }
                    }
                    GffType::Patric => {
                        if !locus_tag.is_empty() {
                            prot_ = format!("{}|{}", prot_, locus_tag);
                        }
                        if contig.starts_with("accn|") {
                            contig = contig[5..].to_string();
                        }
                    }
                    _ => {}
                }

                assert!(!prot_.is_empty());

                let mut prot = unescape(&prot_);
                let gene = unescape(&gene_);
                let product = unescape(&product_);

                if gff_type == GffType::Pgap {
                    pgap_accession(&mut prot, false)?;
                    pgap_accession(&mut contig, lcl)?;
                }
                assert!(!prot.is_empty());

                let locus = Locus::new(
                    line_num + 1,
                    &contig,
                    start,
                    stop,
                    strand_s == "+",
                    partial,
                    0,
                    gene,
                    product,
                )?;

                prot2loci.entry(prot).or_default().insert(locus);

                Ok(())
            })();

            if let Err(e) = parse_result {
                bail!("File {}, line {}: {}", fname, line_num + 1, e);
            }
        }

        Ok(Annot {
            prot2loci,
            fasta2gff_prot: HashMap::new(),
        })
    }

    /// Parse a BED file
    pub fn from_bed(fname: &str) -> Result<Self> {
        if fname.is_empty() {
            bail!("Empty BED file name");
        }

        let mut prot2loci: BTreeMap<String, BTreeSet<Locus>> = BTreeMap::new();

        let file = File::open(fname)?;
        let reader = BufReader::new(file);

        for (line_num, line_result) in reader.lines().enumerate() {
            let mut line = line_result?;
            line = line.trim().to_string();

            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            // Replace spaces with underscores to use tab as delimiter
            line = line.replace(' ', "_");

            let error_prefix = format!("File {}, line {}: ", fname, line_num + 1);

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 6 {
                bail!("{}at least 5 fields are expected in each line", error_prefix);
            }

            let contig = fields[0];
            let start: usize = fields[1].parse()?;
            let stop: usize = fields[2].parse()?;
            let prot = fields[3].trim_matches('_').to_string();
            let _score: f64 = fields[4].parse()?;
            let strand_c = fields[5].chars().next().unwrap_or(' ');

            for c in contig.chars() {
                if !printable_char(c) {
                    bail!(
                        "{}Non-printable character in the sequence identifier: {}",
                        error_prefix,
                        c as u32
                    );
                }
            }

            if start >= stop {
                bail!("{}start should be less than stop", error_prefix);
            }

            if strand_c != '+' && strand_c != '-' {
                bail!("{}strand should be '+' or '-'", error_prefix);
            }

            assert!(!prot.is_empty());

            let locus = Locus::new(
                line_num + 1,
                contig,
                start,
                stop,
                strand_c == '+',
                false,
                0,
                String::new(),
                String::new(),
            )?;

            prot2loci.entry(prot).or_default().insert(locus);
        }

        Ok(Annot {
            prot2loci,
            fasta2gff_prot: HashMap::new(),
        })
    }

    /// Load FASTA-to-GFF protein ID mapping
    pub fn load_fasta2gff_prot(&mut self, fname: &str) -> Result<()> {
        assert!(self.fasta2gff_prot.is_empty());

        let file = File::open(fname)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line?;
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 2 {
                continue;
            }
            let fasta_prot = fields[0].to_string();
            let gff_prot = fields[1].to_string();
            assert!(!gff_prot.is_empty());
            self.fasta2gff_prot.insert(fasta_prot, gff_prot);
        }

        if self.fasta2gff_prot.is_empty() {
            bail!("File {} is empty", fname);
        }

        Ok(())
    }

    /// Load FASTA-to-GFF DNA contig mapping and update contig names
    pub fn load_fasta2gff_dna(&mut self, fname: &str) -> Result<()> {
        let mut gff2fasta: HashMap<String, String> = HashMap::new();

        let file = File::open(fname)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line?;
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 2 {
                continue;
            }
            let fasta_dna = fields[0].to_string();
            let gff_dna = fields[1].to_string();
            assert!(!gff_dna.is_empty());
            gff2fasta.insert(gff_dna, fasta_dna);
        }

        if gff2fasta.is_empty() {
            bail!("File {} is empty", fname);
        }

        // Update contig names in all loci
        let mut updated: BTreeMap<String, BTreeSet<Locus>> = BTreeMap::new();
        for (prot, loci) in &self.prot2loci {
            let mut new_loci = BTreeSet::new();
            for locus in loci {
                let new_contig = gff2fasta.get(&locus.contig).ok_or_else(|| {
                    anyhow::anyhow!(
                        "FASTA DNA contig \"{}\" is not found in GFF-DNA match file \"{}\"",
                        locus.contig,
                        fname
                    )
                })?;
                let mut new_locus = locus.clone();
                new_locus.contig = new_contig.clone();
                new_loci.insert(new_locus);
            }
            updated.insert(prot.clone(), new_loci);
        }
        self.prot2loci = updated;

        Ok(())
    }

    /// Find loci for a FASTA protein ID
    pub fn find_loci(&self, fasta_prot: &str) -> Result<&BTreeSet<Locus>> {
        assert!(!fasta_prot.is_empty());

        let gff_prot = if !self.fasta2gff_prot.is_empty() {
            self.fasta2gff_prot
                .get(fasta_prot)
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "FASTA protein \"{}\" is not found in GFF-protein match file",
                        fasta_prot
                    )
                })?
                .as_str()
        } else {
            fasta_prot
        };

        let loci = self.prot2loci.get(gff_prot).ok_or_else(|| {
            let suffix = if fasta_prot == gff_prot {
                String::new()
            } else {
                format!(" (converted to GFF protein {})", gff_prot)
            };
            anyhow::anyhow!(
                "FASTA protein {}{} is misssing in .gff-file",
                fasta_prot,
                suffix
            )
        })?;
        assert!(!loci.is_empty());

        Ok(loci)
    }
}

fn split_tab<'a>(rest: &mut &'a str) -> Result<&'a str> {
    if let Some(pos) = rest.find('\t') {
        let field = &rest[..pos];
        *rest = &rest[pos + 1..];
        Ok(field)
    } else {
        bail!("Expected tab-delimited field");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn test_data_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amr")
    }

    #[test]
    fn test_parse_gff() {
        let path = test_data_dir().join("test_prot.gff");
        if !path.exists() {
            return;
        }
        let annot = Annot::from_gff(
            path.to_str().unwrap(),
            GffType::Genbank,
            false,
            false,
        );
        assert!(annot.is_ok(), "GFF parse failed: {:?}", annot.err());
        let annot = annot.unwrap();
        assert!(!annot.prot2loci.is_empty());
    }

    #[test]
    fn test_gff_type_names() {
        for name in GffType::NAMES {
            let result = GffType::from_name(name);
            assert!(result.is_ok(), "Failed to parse GFF type: {}", name);
        }
        assert!(GffType::from_name("invalid").is_err());
    }
}
