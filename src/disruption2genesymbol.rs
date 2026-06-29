use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Result};

use crate::common::{Application, Key};
use crate::seq::{
    codon2aa, reverse_dna, DisruptionType, Gencode, Strand, STOP_SUFFIX, TERMINATOR_WORD,
};

const NO_AA: char = '?';

#[derive(Debug, Clone)]
pub struct SymbolRaw {
    pub contig: String,
    pub prot: String,
    pub dtype: DisruptionType,
    pub qstart: usize,
    pub qend: usize,
    pub sstart: usize,
    pub send: usize,
    pub strand: Strand,
    pub stop: bool,
    pub rest: String,
    pub ref_: String,
    pub allele: String,
}

impl SymbolRaw {
    pub fn new(line: &str) -> Result<Self> {
        let bytes = line.as_bytes();
        let mut spans = Vec::with_capacity(3);
        let mut i = 0;
        while spans.len() < 3 {
            while i < bytes.len() && bytes[i].is_ascii_whitespace() {
                i += 1;
            }
            if i == bytes.len() {
                break;
            }
            let start = i;
            while i < bytes.len() && !bytes[i].is_ascii_whitespace() {
                i += 1;
            }
            spans.push((start, i));
        }
        if spans.len() < 3 {
            bail!("line needs contig, protein, and disruption fields");
        }

        let contig = line[spans[0].0..spans[0].1].to_string();
        let prot = line[spans[1].0..spans[1].1].to_string();
        let mut s = line[spans[2].0..spans[2].1].to_string();
        let rest_start = spans[2].1;
        let rest = line[rest_start..]
            .trim_end_matches(['\r', '\n'])
            .to_string();

        let stop = if s.ends_with(STOP_SUFFIX) {
            s.truncate(s.len() - STOP_SUFFIX.len());
            true
        } else {
            false
        };

        let strand_s = s
            .rsplit_once('_')
            .ok_or_else(|| anyhow::anyhow!("missing strand in disruption: {s}"))?
            .1
            .to_string();
        s.truncate(s.len() - strand_s.len() - 1);
        let strand = match strand_s.as_str() {
            "0" => -1,
            "1" => 1,
            _ => bail!("Unknown strand: {strand_s:?}"),
        };

        let send_s = s
            .rsplit_once('_')
            .ok_or_else(|| anyhow::anyhow!("missing send in disruption"))?
            .1
            .to_string();
        s.truncate(s.len() - send_s.len() - 1);
        let sstart_s = s
            .rsplit_once('_')
            .ok_or_else(|| anyhow::anyhow!("missing sstart in disruption"))?
            .1
            .to_string();
        s.truncate(s.len() - sstart_s.len() - 1);
        let qend_s = s
            .rsplit_once('_')
            .ok_or_else(|| anyhow::anyhow!("missing qend in disruption"))?
            .1
            .to_string();
        s.truncate(s.len() - qend_s.len() - 1);
        let qstart_s = s
            .rsplit_once('_')
            .ok_or_else(|| anyhow::anyhow!("missing qstart in disruption"))?
            .1
            .to_string();
        s.truncate(s.len() - qstart_s.len() - 1);

        let qstart = qstart_s.parse::<usize>()?;
        let qend = qend_s.parse::<usize>()?;
        let sstart = sstart_s.parse::<usize>()?;
        let send = send_s.parse::<usize>()?;
        assert!(qstart <= qend);
        assert!(sstart <= send);

        let dtype = DisruptionType::from_name(&s)
            .ok_or_else(|| anyhow::anyhow!("unknown disruption type: {s}"))?;
        assert_ne!(dtype, DisruptionType::None);
        assert_ne!(dtype, DisruptionType::Smooth);

        Ok(Self {
            contig,
            prot,
            dtype,
            qstart,
            qend,
            sstart,
            send,
            strand,
            stop,
            rest,
            ref_: String::new(),
            allele: String::new(),
        })
    }

    pub fn save_text(&self, os: &mut dyn Write, verbose: bool) -> Result<()> {
        assert!(!self.ref_.is_empty());
        write!(os, "{}\t{}\t", self.contig, self.prot)?;
        if verbose {
            write!(
                os,
                "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
                DisruptionType::NAMES[self.dtype as usize],
                self.qstart,
                self.qend,
                self.sstart,
                self.send,
                self.strand,
                self.stop,
                self.ref_,
                self.allele
            )?;
        }
        assert!(!self.ref_.contains('*'));
        let mut allele = self.allele.clone();
        let allele_stop = if allele.ends_with('*') {
            allele.pop();
            true
        } else {
            false
        };
        let allele_size = allele.len();
        if self.stop {
            assert!(allele_stop);
        }
        const DISPLAY_MAX: usize = 1 + 5;
        if allele_size > DISPLAY_MAX {
            allele = "ins".to_string();
        }
        if allele_stop {
            allele.push_str(TERMINATOR_WORD);
        }
        assert!(!allele.contains('*'));

        if self.ref_.len() > DISPLAY_MAX {
            write!(
                os,
                "{}{}_{}{}",
                self.ref_.chars().next().unwrap(),
                self.qstart + 1,
                self.ref_.chars().last().unwrap(),
                self.qstart + self.ref_.len()
            )?;
        } else {
            write!(os, "{}{}", self.ref_, self.qstart + 1)?;
        }

        match self.dtype {
            DisruptionType::Frameshift => {
                assert_eq!(self.ref_.len(), 1);
                assert!(!self.allele.is_empty());
                if allele_stop && allele_size == 0 {
                    write!(os, "{TERMINATOR_WORD}")?;
                } else {
                    write!(os, "{}", self.allele.chars().next().unwrap())?;
                }
                write!(os, "{}", DisruptionType::NAMES[self.dtype as usize])?;
                if allele_stop {
                    write!(os, "{TERMINATOR_WORD}{allele_size}")?;
                }
            }
            DisruptionType::Deletion => {
                if allele.is_empty() {
                    write!(os, "{}", DisruptionType::NAMES[self.dtype as usize])?;
                } else {
                    write!(os, "{allele}")?;
                    if allele_size > DISPLAY_MAX {
                        write!(os, "{}", allele_size - 1)?;
                    }
                }
            }
            DisruptionType::Insertion => {
                assert_eq!(self.ref_.len(), 1);
                assert!(!allele.is_empty());
                write!(os, "{allele}")?;
                if allele_size > DISPLAY_MAX {
                    write!(os, "{}", allele_size - 1)?;
                }
            }
            _ => {}
        }

        write!(
            os,
            "\t{}_{}_{}_{}_{}_{}",
            DisruptionType::NAMES[self.dtype as usize],
            self.qstart,
            self.qend,
            self.sstart,
            self.send,
            if self.strand == 1 { 1 } else { 0 }
        )?;
        if self.stop {
            write!(os, "{STOP_SUFFIX}")?;
        }
        writeln!(os, "\t{}", self.rest)?;
        Ok(())
    }

    pub fn contig2aa(&self, dna: &str, offset: usize, gencode: Gencode) -> Result<char> {
        assert!(self.send <= dna.len());

        if self.strand == 1 {
            let i = self.sstart + offset * 3;
            if i + 3 > self.send {
                return Ok(NO_AA);
            }
            return codon2aa(&dna[i..i + 3], gencode, false);
        }

        assert_eq!(self.strand, -1);
        if self.send < (offset + 1) * 3 {
            return Ok(NO_AA);
        }
        let i = self.send - (offset + 1) * 3;
        assert!(i + 3 <= dna.len());
        if i < self.sstart {
            return Ok(NO_AA);
        }
        let mut s = dna[i..i + 3].to_string();
        reverse_dna(&mut s);
        codon2aa(&s, gencode, false)
    }
}

pub struct ThisApplication {
    application: Application,
}

impl ThisApplication {
    pub const ID_DELIM: char = '|';

    pub fn new() -> Self {
        Self {
            application: Application {
                description: "Convert Disruption::genesymbol_raw() to standard gene symbols according to https://hgvs-nomenclature.org/stable/recommendations/protein/frameshift/.\nA stop codon is 'Ter'.\nPrint: <tab row> where <genesymbol> is inserted before <Disruption::genesymbol_raw()>",
                usage: "Usage: disruption2genesymbol <nucl> <prot> <tab> [-gencode <n>] [-prot_id_pos <n>]",
                positionals: 3,
                needs_arg: false,
                threads_used: false,
                key_args: vec![
                    Key {
                        name: "gencode",
                        flag: false,
                        default_value: "11",
                        acronym: None,
                    },
                    Key {
                        name: "prot_id_pos",
                        flag: false,
                        default_value: "0",
                        acronym: None,
                    },
                ],
            },
        }
    }
}

pub fn body(
    nucl_fname: &Path,
    prot_fname: &Path,
    tab_fname: &Path,
    gencode: Gencode,
    prot_id_pos: usize,
    out: &mut dyn Write,
) -> Result<()> {
    let tab_file = File::open(tab_fname)?;
    let tab_reader = BufReader::new(tab_file);
    let mut symbol_raws = Vec::new();
    for line in tab_reader.lines() {
        symbol_raws.push(SymbolRaw::new(&line?)?);
    }
    if symbol_raws.is_empty() {
        return Ok(());
    }

    let nucl_file = File::open(nucl_fname)?;
    let nucl_reader = BufReader::new(nucl_file);
    let mut name = String::new();
    let mut id = String::new();
    let mut seq = String::new();
    let mut in_record = false;
    let mut after_blank = false;
    let mut line_num = 0usize;
    for line in nucl_reader.lines() {
        line_num += 1;
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            if in_record {
                if name.trim().is_empty() {
                    bail!("Blank sequence name");
                }
                if let Some((i, c)) = seq
                    .chars()
                    .enumerate()
                    .find(|(_, c)| !"-acgtbdhkmnrsvwyACGTBDHKMNRSVWY".contains(*c))
                {
                    bail!(
                        "{name}\nBad sequence character in {:?}: ASCII={} pos={}",
                        id,
                        c as u32,
                        i + 1
                    );
                }
                for symbol_raw in &mut symbol_raws {
                    if symbol_raw.contig == id {
                        for offset in 0.. {
                            let aa = symbol_raw.contig2aa(&seq, offset, gencode)?;
                            if aa == NO_AA {
                                break;
                            }
                            symbol_raw.allele.push(aa);
                            if aa == '*' {
                                break;
                            }
                        }
                    }
                }
            }
            name = header.replace('\t', " ");
            if name.trim().is_empty() {
                bail!("Blank sequence name");
            }
            id = name.split(' ').next().unwrap_or("").to_string();
            seq.clear();
            in_record = true;
            after_blank = false;
        } else if line.trim().is_empty() {
            if in_record {
                if name.trim().is_empty() {
                    bail!("Blank sequence name");
                }
                if let Some((i, c)) = seq
                    .chars()
                    .enumerate()
                    .find(|(_, c)| !"-acgtbdhkmnrsvwyACGTBDHKMNRSVWY".contains(*c))
                {
                    bail!(
                        "{name}\nBad sequence character in {:?}: ASCII={} pos={}",
                        id,
                        c as u32,
                        i + 1
                    );
                }
                for symbol_raw in &mut symbol_raws {
                    if symbol_raw.contig == id {
                        for offset in 0.. {
                            let aa = symbol_raw.contig2aa(&seq, offset, gencode)?;
                            if aa == NO_AA {
                                break;
                            }
                            symbol_raw.allele.push(aa);
                            if aa == '*' {
                                break;
                            }
                        }
                    }
                }
                seq.clear();
                in_record = false;
                after_blank = true;
            } else if !after_blank {
                bail!("Error in Multifasta, line {}", line_num + 1);
            }
        } else {
            if !in_record {
                bail!("Error in Multifasta, line {}", line_num + 1);
            }
            seq.push_str(
                &line
                    .trim_end_matches(|c: char| c.is_ascii_whitespace())
                    .to_ascii_lowercase(),
            );
        }
    }
    if in_record {
        if name.trim().is_empty() {
            bail!("Blank sequence name");
        }
        if let Some((i, c)) = seq
            .chars()
            .enumerate()
            .find(|(_, c)| !"-acgtbdhkmnrsvwyACGTBDHKMNRSVWY".contains(*c))
        {
            bail!(
                "{name}\nBad sequence character in {:?}: ASCII={} pos={}",
                id,
                c as u32,
                i + 1
            );
        }
        for symbol_raw in &mut symbol_raws {
            if symbol_raw.contig == id {
                for offset in 0.. {
                    let aa = symbol_raw.contig2aa(&seq, offset, gencode)?;
                    if aa == NO_AA {
                        break;
                    }
                    symbol_raw.allele.push(aa);
                    if aa == '*' {
                        break;
                    }
                }
            }
        }
    }

    let prot_file = File::open(prot_fname)?;
    let prot_reader = BufReader::new(prot_file);
    let mut proteins: HashMap<String, String> = HashMap::new();
    let mut prot_name = String::new();
    let mut id_whole = String::new();
    let mut pep = String::new();
    let mut in_record = false;
    let mut after_blank = false;
    let mut line_num = 0usize;
    for line in prot_reader.lines() {
        line_num += 1;
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            if in_record {
                if prot_name.trim().is_empty() {
                    bail!("Blank sequence name");
                }
                if let Some((i, c)) = pep
                    .chars()
                    .enumerate()
                    .find(|(_, c)| !"-ACDEFGHIKLMNPQRSTVWYXBZJUO*".contains(*c))
                {
                    bail!(
                        "{prot_name}\nBad sequence character in {:?}: ASCII={} pos={}",
                        id_whole,
                        c as u32,
                        i + 1
                    );
                }
                if pep.find('*').is_some_and(|pos| pos + 1 != pep.len()) {
                    bail!("{prot_name}\nPeptide has an internal stop codon");
                }
                proteins.insert(id_whole.clone(), pep.clone());
            }
            prot_name = header.replace('\t', " ");
            if prot_name.trim().is_empty() {
                bail!("Blank sequence name");
            }
            id_whole = prot_name.split(' ').next().unwrap_or("").to_string();
            pep.clear();
            in_record = true;
            after_blank = false;
        } else if line.trim().is_empty() {
            if in_record {
                if prot_name.trim().is_empty() {
                    bail!("Blank sequence name");
                }
                if let Some((i, c)) = pep
                    .chars()
                    .enumerate()
                    .find(|(_, c)| !"-ACDEFGHIKLMNPQRSTVWYXBZJUO*".contains(*c))
                {
                    bail!(
                        "{prot_name}\nBad sequence character in {:?}: ASCII={} pos={}",
                        id_whole,
                        c as u32,
                        i + 1
                    );
                }
                if pep.find('*').is_some_and(|pos| pos + 1 != pep.len()) {
                    bail!("{prot_name}\nPeptide has an internal stop codon");
                }
                proteins.insert(id_whole.clone(), pep.clone());
                pep.clear();
                in_record = false;
                after_blank = true;
            } else if !after_blank {
                bail!("Error in Multifasta, line {}", line_num + 1);
            }
        } else {
            if !in_record {
                bail!("Error in Multifasta, line {}", line_num + 1);
            }
            pep.push_str(
                &line
                    .trim_end_matches(|c: char| c.is_ascii_whitespace())
                    .to_ascii_uppercase(),
            );
        }
    }
    if in_record {
        if prot_name.trim().is_empty() {
            bail!("Blank sequence name");
        }
        if let Some((i, c)) = pep
            .chars()
            .enumerate()
            .find(|(_, c)| !"-ACDEFGHIKLMNPQRSTVWYXBZJUO*".contains(*c))
        {
            bail!(
                "{prot_name}\nBad sequence character in {:?}: ASCII={} pos={}",
                id_whole,
                c as u32,
                i + 1
            );
        }
        if pep.find('*').is_some_and(|pos| pos + 1 != pep.len()) {
            bail!("{prot_name}\nPeptide has an internal stop codon");
        }
        proteins.insert(id_whole, pep);
    }

    for symbol_raw in &mut symbol_raws {
        for (id_whole, pep) in &proteins {
            let id = if prot_id_pos != 0 {
                let vec = id_whole
                    .split(ThisApplication::ID_DELIM)
                    .collect::<Vec<_>>();
                if prot_id_pos - 1 >= vec.len() {
                    bail!(
                        "Protein identifier position {} is outside of the list of identifiers: {:?}",
                        prot_id_pos,
                        id_whole
                    );
                }
                vec[prot_id_pos - 1]
            } else {
                id_whole.as_str()
            };
            if symbol_raw.prot == id {
                symbol_raw.ref_ = pep[symbol_raw.qstart..symbol_raw.qend].to_string();
            }
        }
    }

    for symbol_raw in &symbol_raws {
        symbol_raw.save_text(out, false)?;
    }
    Ok(())
}

pub fn main(argv: Vec<OsString>, out: &mut dyn Write) -> Result<i32> {
    let mut app = ThisApplication::new();
    let run = match app.application.run(argv) {
        Ok(run) => run,
        Err(code) => return Ok(code),
    };
    let nucl_fname = PathBuf::from(&run.positional_args[0]);
    let prot_fname = PathBuf::from(&run.positional_args[1]);
    let tab_fname = PathBuf::from(&run.positional_args[2]);
    let gencode = run.key_args["gencode"].parse::<Gencode>()?;
    let prot_id_pos = run.key_args["prot_id_pos"].parse::<usize>()?;

    body(
        &nucl_fname,
        &prot_fname,
        &tab_fname,
        gencode,
        prot_id_pos,
        out,
    )?;
    Ok(0)
}
