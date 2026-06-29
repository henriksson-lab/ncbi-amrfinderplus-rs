use amrfinder::seq::{
    aa2num, aa_match, codon2aa, count_insertions, dna2codons_len, get_intersect_nucleotide,
    get_union_nucleotide, hsp2exons, more_general_aminoacid, nucleotide_match, nucleotides2wild,
    reverse_dna, sparse_seq_len, strand2char, wild2nucleotides, Cds, Disruption, DisruptionType,
    Dna, DnaAnnot, Exon, Hsp, HspMerge, Interval, Intron, Multifasta, Mutation, Peptide,
    PeptideOrf, Seq, SubstMat,
};
use std::io::Write;

#[test]
fn test_nucleotide_match() {
    assert!(nucleotide_match('A', 'A'));
    assert!(nucleotide_match('a', 'A'));
    assert!(!nucleotide_match('a', 'c'));
    assert!(nucleotide_match('A', 'C'));
    assert!(nucleotide_match('N', 'A'));
    assert!(nucleotide_match('r', 'a'));
    assert!(nucleotide_match('r', 'g'));
    assert!(!nucleotide_match('r', 'c'));
    assert!(nucleotide_match('R', 'C'));
    assert!(nucleotide_match('-', 'A'));
}

#[test]
fn nucleotide_set_conversions_match_original_sparse_semantics() {
    let mut acgtb = [false; 5];

    assert_eq!(wild2nucleotides('r', &mut acgtb), 2);
    assert_eq!(acgtb, [true, false, true, false, false]);
    assert_eq!(nucleotides2wild(&acgtb), 'r');

    assert_eq!(wild2nucleotides('R', &mut acgtb), 3);
    assert_eq!(acgtb, [true, false, true, false, true]);
    assert_eq!(nucleotides2wild(&acgtb), 'R');

    assert_eq!(wild2nucleotides('-', &mut acgtb), 1);
    assert_eq!(acgtb, [false, false, false, false, true]);
    assert_eq!(nucleotides2wild(&acgtb), '-');
}

#[test]
fn union_and_intersection_use_original_wildcard_rules() {
    assert_eq!(get_union_nucleotide("ag"), 'r');
    assert_eq!(get_union_nucleotide("ct"), 'y');
    assert_eq!(get_union_nucleotide(""), ' ');

    assert_eq!(get_intersect_nucleotide("rn"), 'r');
    assert_eq!(get_intersect_nucleotide("ry"), ' ');
    assert_eq!(get_intersect_nucleotide("-A"), '-');
}

#[test]
fn reverse_dna_complements_in_place_and_preserves_gaps() {
    let mut seq = String::from("aC-gTRYN");

    assert_eq!(reverse_dna(&mut seq), "NRYAc-Gt");
    assert_eq!(seq, "NRYAc-Gt");
}

#[test]
fn codon2aa_translates_and_lowercases_possible_start_codon() {
    assert_eq!(codon2aa("ATG", 11, false).unwrap(), 'M');
    assert_eq!(codon2aa("ATG", 11, true).unwrap(), 'm');
    assert_eq!(codon2aa("TGA", 11, false).unwrap(), '*');
    assert!(codon2aa("ATG", 2, true).is_err());
}

#[test]
fn test_aa_match() {
    assert!(aa_match('A', 'A'));
    assert!(aa_match('x', 'M'));
    assert!(!aa_match('X', 'M'));
    assert!(aa_match('b', 'd'));
    assert!(aa_match('b', 'n'));
    assert!(!aa_match('B', 'D'));
    assert!(!aa_match('b', 'a'));
    assert!(aa_match('z', 'e'));
    assert!(aa_match('j', 'i'));
    assert!(!aa_match('-', 'A'));
}

#[test]
fn original_sequence_free_functions_cover_sparse_and_peptide_edges() {
    assert_eq!(count_insertions("A-C--D"), 3);
    assert_eq!(sparse_seq_len("A-C--D"), 3);
    assert_eq!(dna2codons_len(0), 0);
    assert_eq!(dna2codons_len(1), 1);
    assert_eq!(dna2codons_len(3), 1);
    assert_eq!(dna2codons_len(4), 2);

    assert_eq!(aa2num('A'), 0);
    assert_eq!(aa2num('y'), 19);
    assert_eq!(aa2num('*'), 20);
    assert_eq!(aa2num('X'), 21);

    assert!(more_general_aminoacid('x', 'M'));
    assert!(!more_general_aminoacid('X', 'M'));
    assert!(more_general_aminoacid('b', 'd'));
    assert!(!more_general_aminoacid('B', 'D'));
    assert!(!more_general_aminoacid('D', 'B'));
    assert!(aa_match('u', '*'));
    assert!(aa_match('*', 'o'));
    assert!(aa_match('-', '-'));
    assert!(!aa_match('-', 'A'));
}

#[test]
fn subst_mat_parses_scores_and_matches_original_special_cases() {
    let mut matrix = tempfile::NamedTempFile::new().unwrap();
    write!(matrix, "  {}\n", SubstMat::CHARS).unwrap();
    for c1 in SubstMat::CHARS.chars() {
        write!(matrix, "{c1}").unwrap();
        for c2 in SubstMat::CHARS.chars() {
            let score = if c1 == c2 { 5 } else { -4 };
            write!(matrix, " {score}").unwrap();
        }
        writeln!(matrix).unwrap();
    }

    let mat = SubstMat::new(matrix.path().to_str().unwrap()).unwrap();
    mat.qc();

    assert!(mat.good_index('A' as usize));
    assert_eq!(mat.char2score('A', 'A').unwrap(), 5.0);
    assert_eq!(mat.char2score('A', 'R').unwrap(), -4.0);
    assert_eq!(
        mat.char2score('B', 'D').unwrap(),
        mat.char2score('X', 'D').unwrap()
    );
    assert_eq!(mat.char2score('*', 'A').unwrap(), -10.0);
    assert_eq!(mat.char2score('-', 'A').unwrap(), -1.0);
    assert_eq!(mat.get_substitution_dist('A' as usize, 'R' as usize), 18.0);
    assert_eq!(mat.get_deletion_dist('A' as usize, -2.0), 9.0);
    assert_eq!(SubstMat::char2score_static(None, 'b', 'd').unwrap(), 1.0);
    assert_eq!(SubstMat::char2score_static(None, 'B', 'D').unwrap(), 0.0);

    let mut out = Vec::new();
    mat.save_text(&mut out).unwrap();
    let text = String::from_utf8(out).unwrap();
    assert!(text.starts_with(" * A B C D E F G H I K L M N P Q R S T V W X Y Z\n* 5 -4"));
    assert!(text
        .ends_with("Z -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 5\n"));
}

#[test]
fn seq_members_match_id_description_sparse_and_counts() {
    let mut seq = Seq::new("gi|123|abc description [Taxon name]", "A-C--NN", false);

    assert_eq!(seq.seq, "ACNN");
    assert_eq!(seq.get_id_size(), 10);
    assert_eq!(seq.get_id(), "gi|123|abc");
    assert_eq!(seq.get_gi().unwrap(), 123);
    assert_eq!(seq.get_taxon_start(), Some(23));
    assert_eq!(seq.get_description(false), "description [Taxon name]");
    assert_eq!(seq.get_description(true), "description");
    assert_eq!(seq.get_xs(|c| c == 'N'), 2);
    assert_eq!(seq.get_contiguous_xs(|c| c == 'N'), 2);

    seq.append_id("sfx");
    assert_eq!(seq.name, "gi|123|abc-sfx description [Taxon name]");

    seq.append_id("");
    assert_eq!(seq.name, "gi|123|abc-sfx description [Taxon name]");

    let counts = seq.get_char_count();
    assert_eq!(counts.get(&'N'), Some(&2));
}

#[test]
fn sequence_name_constructors_validate_like_upstream_seq_constructor() {
    assert!(std::panic::catch_unwind(|| Seq::new("", "acgt", false)).is_err());
    assert!(std::panic::catch_unwind(|| Seq::with_len("bad\tname", 4, false)).is_err());
    assert!(std::panic::catch_unwind(|| Dna::with_len("", 4, false)).is_err());
    assert!(std::panic::catch_unwind(|| Peptide::new("bad\tname", "ACD", false)).is_err());
}

#[test]
fn multifasta_reads_dna_records_like_upstream_seq_constructor() {
    let mut fasta_file = tempfile::NamedTempFile::new().unwrap();
    write!(fasta_file, ">seq1\tdesc\nacg-\nTT \n\n>seq2\nNN--aa\n").unwrap();

    let mut fasta = Multifasta::new(fasta_file.path(), false, 1000).unwrap();
    assert!(fasta.next().unwrap());

    let seq1 = Dna::from_multifasta(&mut fasta, 500, false).unwrap();
    assert_eq!(seq1.name, "seq1 desc");
    assert_eq!(seq1.seq, "acgtt");
    assert!(!seq1.sparse);

    assert!(fasta.next().unwrap());
    let seq2 = Dna::from_multifasta(&mut fasta, 500, true).unwrap();
    assert_eq!(seq2.name, "seq2");
    assert_eq!(seq2.seq, "nn--aa");
    assert!(seq2.sparse);

    assert!(!fasta.next().unwrap());
}

#[test]
fn multifasta_reads_peptide_records_uppercase_and_checks_alphabet() {
    let mut fasta_file = tempfile::NamedTempFile::new().unwrap();
    write!(fasta_file, ">prot\nmk-x*\n").unwrap();

    let mut fasta = Multifasta::new(fasta_file.path(), true, 1000).unwrap();
    let pep = Peptide::from_multifasta(&mut fasta, Peptide::STD_AVE_LEN, true).unwrap();

    assert_eq!(pep.name, "prot");
    assert_eq!(pep.seq, "MK-X*");
    assert!(pep.sparse);
    assert!(!fasta.next().unwrap());
}

#[test]
fn multifasta_rejects_non_header_after_blank_boundary() {
    let mut fasta_file = tempfile::NamedTempFile::new().unwrap();
    write!(fasta_file, ">seq1\nacgt\n\norphan\n").unwrap();

    let mut fasta = Multifasta::new(fasta_file.path(), false, 1000).unwrap();
    let seq1 = Dna::from_multifasta(&mut fasta, 500, false).unwrap();
    assert_eq!(seq1.seq, "acgt");

    let err = fasta.next().unwrap_err();
    assert!(err.to_string().contains("Error in Multifasta, line 4"));
}

#[test]
fn dna_members_cover_complexity_repeats_reverse_and_polya() {
    let dna = Dna::new("contig one [Taxon]", "acgtacgt", false);
    assert_eq!(dna.get_id(), "contig");
    assert_eq!(dna.get_description(true), "one");
    assert_eq!(dna.get_seq_alphabet(), "acgtbdhkmnrsvwyACGTBDHKMNRSVWY");
    assert!(Dna::is_ambiguous('n'));
    assert_eq!(dna.get_xs(), 0);
    assert!(dna.get_complexity() > 0.0);

    let rev = dna.make_complementary();
    assert_eq!(rev.name, "contig.rev");
    assert_eq!(rev.seq, "acgtacgt");

    let mut repeat = Dna::new("repeat", "aaaatttc", false);
    assert_eq!(repeat.mono_nuc2n(3), 7);
    assert_eq!(repeat.seq, "nnnnnnnc");

    let mut tail = Dna::new("tail", "ggccaaaaa", false);
    assert!(tail.poly_a_window(4, 5));
    assert_eq!(tail.find_poly_a(4), 4);
    assert_eq!(tail.remove_poly_a(4), 5);
    assert_eq!(tail.seq, "ggcc");

    let mut all_tail = Dna::new("all_tail", "aaaa", false);
    assert_eq!(all_tail.find_poly_a(4), 0);
    assert_eq!(all_tail.remove_poly_a(4), 4);
    assert!(all_tail.seq.is_empty());
}

#[test]
fn dna_translation_cds_and_orf_peptides_follow_reference_rules() {
    let dna = Dna::new(
        "gene",
        "atgaaataatttatgccc taa".replace(' ', "").as_str(),
        false,
    );

    let (pep, translation_start) = dna.make_peptide(1, 11, true, true).unwrap();
    assert_eq!(translation_start, 0);
    assert_eq!(pep.name, "gene.fr1");
    assert_eq!(pep.seq, "mK*FmP*");
    assert!(pep.pseudo);

    let cds = Dna::new("cds", "atgaaataa", false);
    let prot = cds.cds2prot(11, false, false, true, false).unwrap();
    assert_eq!(prot.name, "cds");
    assert_eq!(prot.seq, "mK");

    let extra_stop = Dna::new("extra", "atgtaaaaataa", false);
    assert!(extra_stop.cds2prot(11, false, false, true, false).is_err());
    assert_eq!(
        extra_stop
            .cds2prot(11, false, false, true, true)
            .unwrap()
            .seq,
        "mK*"
    );

    let peptides = dna.get_peptides(1, 11, 2).unwrap();
    assert_eq!(peptides.len(), 2);
    assert_eq!(peptides[0].name, "gene:1..6");
    assert_eq!(peptides[0].seq, "MK");
    assert_eq!(peptides[1].name, "gene:13..18");
    assert_eq!(peptides[1].seq, "MP");

    let minus = Dna::new("minus", "tttcat", false);
    let (pep_minus, translation_start) = minus.make_peptide(-1, 11, true, true).unwrap();
    assert_eq!(translation_start, 6);
    assert_eq!(pep_minus.seq, "mK");
}

#[test]
fn dna_ambiguity_prefix_suffix_and_trim_n_follow_original_rules() {
    let dna = Dna::new("ambig", "NNRacgtBD", true);
    assert!(dna.contains_ambiguity());
    assert_eq!(dna.get_ambiguous_prefix_end(), 3);
    assert_eq!(dna.get_ambiguous_suffix_start(), 7);

    let plain = Dna::new("plain", "acgt", false);
    assert!(!plain.contains_ambiguity());
    assert_eq!(plain.get_ambiguous_prefix_end(), 0);
    assert_eq!(plain.get_ambiguous_suffix_start(), 4);

    let mut prefix = Dna::new("prefix", "NNRacgt", true);
    assert!(prefix.delete_ambiguous_prefix());
    assert_eq!(prefix.seq, "acgt");
    assert!(!prefix.delete_ambiguous_prefix());

    let mut suffix = Dna::new("suffix", "acgtBD", true);
    assert!(suffix.delete_ambiguous_suffix());
    assert_eq!(suffix.seq, "acgt");
    assert!(!suffix.delete_ambiguous_suffix());

    let mut not_gap_coded = Dna::new("trim", "nnNacgNnn", true);
    not_gap_coded.trim_n(false);
    assert_eq!(not_gap_coded.seq, "acg");

    let mut gap_coded = Dna::new("trim", "nnNacgNnn", true);
    gap_coded.trim_n(true);
    assert_eq!(gap_coded.seq, "nnNacgNnn");

    let mut upper_only = Dna::new("trim", "NNacgNN", true);
    upper_only.trim_n(true);
    assert_eq!(upper_only.seq, "acg");
}

#[test]
fn mutation_parse_save_qc_order_and_replace_follow_original_rules() {
    let prot = Mutation::from_line(true, " gyrA-S83L ").unwrap();
    assert!(prot.prot);
    assert_eq!(prot.gene_name, "gyrA");
    assert_eq!(prot.pos, 82);
    assert_eq!(prot.ref_, "S");
    assert_eq!(prot.allele, "L");
    assert!(!prot.frameshift);
    assert!(!prot.ambig);
    prot.qc().unwrap();

    let mut out = Vec::new();
    prot.save_text(&mut out).unwrap();
    assert_eq!(String::from_utf8(out).unwrap(), "gyrA-S83L");

    let insertion = Mutation::from_line(true, "gyrA-ins83A").unwrap();
    assert_eq!(insertion.ref_, "");
    let deletion = Mutation::from_line(true, "gyrA-S83del").unwrap();
    assert_eq!(deletion.allele, "");
    let frameshift = Mutation::from_line(true, "gyrA-S83fs").unwrap();
    assert!(frameshift.frameshift);
    assert_eq!(frameshift.allele, "");
    let ambiguous = Mutation::new("gyrA".to_string(), 82, "S".to_string(), "X".to_string());
    assert!(ambiguous.ambig);

    let dna = Mutation::from_line(false, "a3t").unwrap();
    assert!(!dna.prot);
    assert_eq!(dna.pos, 2);
    assert_eq!(dna.ref_, "a");
    assert_eq!(dna.allele, "t");
    assert_eq!(dna.stop(), 3);
    dna.qc().unwrap();

    let mut ref_dna = Dna::new("ref", "ccacc", false);
    dna.replace(&mut ref_dna);
    assert_eq!(ref_dna.seq, "cctcc");

    let dna_insertion = Mutation::from_line(false, "INS3t").unwrap();
    assert_eq!(dna_insertion.ref_, "");
    let dna_deletion = Mutation::from_line(false, "a3DEL").unwrap();
    assert_eq!(dna_deletion.allele, "");

    let mut muts = vec![
        Mutation::from_line(true, "gyrB-S83L").unwrap(),
        prot.clone(),
        Mutation::from_line(false, "a2t").unwrap(),
    ];
    muts.sort();
    assert_eq!(muts[0], Mutation::from_line(false, "a2t").unwrap());
    assert_eq!(muts[1], prot);

    let bad_protein_ref = Mutation::new("gyrA".to_string(), 82, "X".to_string(), "L".to_string());
    assert!(bad_protein_ref.qc().is_err());
    let bad_dna_ref = Mutation::new("".to_string(), 2, "R".to_string(), "t".to_string());
    assert!(bad_dna_ref.qc().is_err());
}

#[test]
fn peptide_members_cover_stop_partial_complexity_and_reduced_alphabet() {
    let mut peptide = Peptide::new("prot partial protein", "ABZJUOX*", false);

    assert_eq!(peptide.get_seq_alphabet(), "ACDEFGHIKLMNPQRSTVWYXBZJUO*");
    assert!(Peptide::is_ambiguous('X'));
    assert!(!peptide.has_inside_stop());
    assert!(peptide.is_description_partial());
    assert!(Peptide::is_start_aa('m'));
    assert!(!Peptide::is_start_aa('M'));
    assert!(peptide.get_complexity() >= 0.0);

    assert!(peptide.trim_stop());
    assert_eq!(peptide.seq, "ABZJUOX");
    assert_eq!(peptide.ambig2x(), 5);
    assert_eq!(peptide.seq, "AXXXXXX");

    let mut gbmr = Peptide::new("gbmr", "ADKYFILVMCWGXP*", false);
    gbmr.to_gbmr4();
    assert_eq!(gbmr.seq, "AAAYYYYYYYYGXP*");
}

#[test]
fn peptide_similarity_methods_follow_reference_gap_scoring() {
    let mut matrix = tempfile::NamedTempFile::new().unwrap();
    write!(matrix, "  {}\n", SubstMat::CHARS).unwrap();
    for c1 in SubstMat::CHARS.chars() {
        write!(matrix, "{c1}").unwrap();
        for c2 in SubstMat::CHARS.chars() {
            let score = if c1 == c2 { 5 } else { -4 };
            write!(matrix, " {score}").unwrap();
        }
        writeln!(matrix).unwrap();
    }
    let mat = SubstMat::new(matrix.path().to_str().unwrap()).unwrap();

    let left = Peptide::new("left", "AC--DE", true);
    let right = Peptide::new("right", "A--CDE", true);

    assert_eq!(left.get_self_similarity(&mat, 0, 0), 20.0);
    assert_eq!(left.get_self_similarity(&mat, 1, 5), 10.0);
    assert_eq!(left.get_similarity(&right, &mat, 2.0, 0.5), 12.0);
}

#[test]
fn peptide_orf_members_follow_reference_start_stop_and_coordinates() {
    let peptide = Peptide::new("prot", "AKlm*QQ", false);

    let orf = PeptideOrf::new(12, 1, &peptide, 2);
    assert_eq!(orf.translation_start, 12);
    assert_eq!(orf.strand, 1);
    assert_eq!(orf.start, 2);
    assert!(orf.start_m);
    assert_eq!(orf.stop, 4);
    assert!(orf.stop_terminator);
    assert_eq!(orf.size(), 2);
    assert!(orf.good(2));
    assert!(!orf.good(3));
    assert_eq!(orf.dna_pos(3), 21);
    assert_eq!(orf.cds_start(), 18);
    assert_eq!(orf.cds_stop(), 27);
    assert_eq!(orf.to_peptide(&peptide).seq, "lm*");

    let minus = PeptideOrf::new_fields(30, -1, 2, true, 4, true);
    minus.qc();
    assert_eq!(minus.dna_pos(3), 21);
    assert_eq!(minus.cds_start(), 24);
    assert_eq!(minus.cds_stop(), 15);

    let mut text = Vec::new();
    orf.save_text(&mut text).unwrap();
    assert_eq!(String::from_utf8(text).unwrap(), "12\t1\t2\t4\t1\t1");
    assert_eq!(PeptideOrf::from_text("12 1 2 4 1 1").unwrap(), orf);

    let orfs = peptide.get_peptide_orfs(12, 1, true, true, 1);
    assert_eq!(orfs.len(), 2);
    assert_eq!(orfs[0], orf);
    assert_eq!(orfs[1].start, 0);
    assert!(!orfs[1].start_m);
    assert_eq!(orfs[1].to_peptide(&peptide).seq, "AKlm*");
}

#[test]
fn test_bacterial_start_codon_liv_matches_reference_methionine() {
    let line = concat!(
        "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
        "target\t1\t3\t3\t1\t3\t3\tMKA\tLKA"
    );

    let hsp = Hsp::from_blast_line(line, true, true, true, false, true).unwrap();

    assert_eq!(hsp.sseq, "MKA");
    assert_eq!(hsp.nident, 3);
    assert_eq!(hsp.rel_identity(), 1.0);
}

#[test]
fn hsp_from_blast_line_accepts_whitespace_like_stream_extraction() {
    let line = concat!(
        "ref|1|1|fam|gene|resistance|2|subclass|class|product ",
        "target 1 3 3 1 3 3 MKA MRA"
    );

    let hsp = Hsp::from_blast_line(line, true, true, true, false, false).unwrap();

    assert_eq!(
        hsp.qseqid,
        "ref|1|1|fam|gene|resistance|2|subclass|class|product"
    );
    assert_eq!(hsp.sseqid, "target");
    assert_eq!(hsp.qseq, "MKA");
    assert_eq!(hsp.sseq, "MRA");
    assert_eq!(hsp.nident, 2);
}

#[test]
fn hsp_from_blast_line_uses_stream_extraction_for_tabular_ids_with_spaces() {
    let line = concat!(
        "ref|1|1|fam|gene|resistance|2|subclass|class|product with spaces\t",
        "target\t1\t3\t3\t1\t3\t3\tMKA\tMRA"
    );

    assert!(Hsp::from_blast_line(line, true, true, true, false, false).is_err());
}

#[test]
fn finish_hsp_moves_dashes_right_then_trims_terminal_gaps() {
    let line = concat!(
        "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
        "target\t1\t4\t4\t1\t3\t3\tAAAA\t-AAA"
    );

    let hsp = Hsp::from_blast_line(line, true, true, true, false, false).unwrap();

    assert_eq!(hsp.qseq, "AAA");
    assert_eq!(hsp.sseq, "AAA");
    assert_eq!(hsp.q_int, Interval::new(0, 3, 0));
    assert_eq!(hsp.s_int, Interval::new(0, 3, 0));
    assert_eq!(hsp.length, 3);
    assert_eq!(hsp.nident, 3);
    assert_eq!(hsp.sgap, 0);
}

#[test]
fn hsp_qc_and_save_text_use_original_summary_fields() {
    let hsp = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "target\t1\t3\t3\t1\t3\t3\tMKA\tMRA"
        ),
        true,
        true,
        true,
        false,
        false,
    )
    .unwrap();
    hsp.qc();

    let mut out = Vec::new();
    hsp.save_text(&mut out).unwrap();
    let text = String::from_utf8(out).unwrap();
    assert!(text.starts_with(
        "Hsp: merged=0 ref|1|1|fam|gene|resistance|2|subclass|class|product(3) 1-3 target(3) 1-3"
    ));
    assert!(text.contains(" qLen=3 sLen=3 length=3 nident=2"));
    assert!(text.ends_with(" #disrs=0"));
}

#[test]
fn finish_hsp_does_not_count_gap_columns_as_identity() {
    let mut hsp = Hsp::default();
    hsp.qseqid = "dna_ref".to_string();
    hsp.sseqid = "dna_target".to_string();
    hsp.q_int = Interval::new(0, 2, 0);
    hsp.s_int = Interval::new(0, 3, 1);
    hsp.qlen = 2;
    hsp.slen = 3;
    hsp.qseq = "A-C".to_string();
    hsp.sseq = "AAA".to_string();

    hsp.finish_hsp(false, false);

    assert_eq!(hsp.qgap, 1);
    assert_eq!(hsp.sgap, 0);
    assert_eq!(hsp.nident, 2);
    hsp.qc();
}

#[test]
fn finish_hsp_counts_ambiguity_with_original_case_sensitive_alphabets() {
    let mut dna_hsp = Hsp::default();
    dna_hsp.qseqid = "dna_ref".to_string();
    dna_hsp.sseqid = "dna_target".to_string();
    dna_hsp.q_int = Interval::new(0, 4, 0);
    dna_hsp.s_int = Interval::new(0, 4, 1);
    dna_hsp.qlen = 4;
    dna_hsp.slen = 4;
    dna_hsp.qseq = "ACgt".to_string();
    dna_hsp.sseq = "acGT".to_string();

    dna_hsp.finish_hsp(false, false);

    assert_eq!(dna_hsp.qx, 2);
    assert_eq!(dna_hsp.sx, 2);

    let mut peptide_hsp = Hsp::default();
    peptide_hsp.q_prot = true;
    peptide_hsp.s_prot = true;
    peptide_hsp.a_prot = true;
    peptide_hsp.qseqid = "prot_ref".to_string();
    peptide_hsp.sseqid = "prot_target".to_string();
    peptide_hsp.q_int = Interval::new(0, 4, 0);
    peptide_hsp.s_int = Interval::new(0, 4, 0);
    peptide_hsp.qlen = 4;
    peptide_hsp.slen = 4;
    peptide_hsp.qseq = "XxBZ".to_string();
    peptide_hsp.sseq = "xXbz".to_string();

    peptide_hsp.finish_hsp(false, false);

    assert_eq!(peptide_hsp.qx, 3);
    assert_eq!(peptide_hsp.sx, 1);
}

#[test]
fn hsp_position_maps_follow_query_gaps_like_upstream() {
    let hsp = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "target\t1\t4\t5\t10\t14\t20\tA-BCD\tAXBCD"
        ),
        true,
        true,
        true,
        false,
        false,
    )
    .unwrap();

    assert_eq!(hsp.pos2real_q(1, true), 2);
    assert_eq!(hsp.pos2real_q(1, false), 0);
    assert_eq!(hsp.pos2q(1, true), 1);
    assert_eq!(hsp.q2pos(1, true), 2);
    assert_eq!(hsp.q2pos(1, false), 0);
    assert_eq!(hsp.q_map(5), "ABCD-");
}

#[test]
fn hsp_position_maps_follow_subject_gaps_and_minus_strand() {
    let plus = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t4\t4\t10\t18\t30\tABCD\tA-CD"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    assert_eq!(plus.pos2real_s(1, true), 2);
    assert_eq!(plus.pos2real_s(1, false), 0);
    assert_eq!(plus.pos2s(1, true), 12);
    assert_eq!(plus.s2pos(12, true), 2);
    assert_eq!(plus.s2pos(12, false), 0);
    assert_eq!(plus.q_map(4), "A-CD");

    let minus = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t4\t4\t18\t10\t30\tABCD\tA-CD"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    assert_eq!(minus.s_int, Interval::new(9, 18, -1));
    assert_eq!(minus.pos2s(1, true), 15);
    assert_eq!(minus.s2pos(15, true), 2);
    assert_eq!(minus.s2pos(15, false), 0);
}

#[test]
fn cds_methods_follow_original_interval_rules() {
    let cds = Cds::new_fields(300, 90, "refA".to_string(), 0.75);
    cds.qc();

    assert!(!cds.strand());
    assert_eq!(cds.left(), 90);
    assert_eq!(cds.right(), 300);
    assert_eq!(cds.start_human(), 91);
    assert_eq!(cds.stop_human(), 300);
    assert_eq!(cds.size(), 210);
    assert_eq!(cds.size_effective(), 150);

    let overlapping = Cds::new_fields(260, 410, "refB".to_string(), 0.8);
    assert_eq!(cds.get_overlap(&overlapping), 40);
    assert!(!cds.coexists(&overlapping));

    let downstream_plus = Cds::new_fields(360, 510, "refC".to_string(), 0.8);
    assert!(cds.coexists(&downstream_plus));

    let mut out = Vec::new();
    cds.save_text(&mut out).unwrap();
    assert_eq!(String::from_utf8(out).unwrap(), "90\t300\t0\t69");
}

#[test]
fn dna_annot_run_sorts_and_links_best_cds_chain() {
    let mut annot = DnaAnnot::new();
    annot
        .cdss
        .push(Cds::new_fields(240, 390, "refC".to_string(), 0.92));
    annot
        .cdss
        .push(Cds::new_fields(0, 150, "refA".to_string(), 0.95));
    annot
        .cdss
        .push(Cds::new_fields(30, 180, "refB".to_string(), 0.9));

    let last = annot.run().unwrap();

    assert_eq!(annot.cdss[0].ref_prot, "refA");
    assert_eq!(annot.cdss[1].ref_prot, "refB");
    assert_eq!(annot.cdss[2].ref_prot, "refC");
    assert_eq!(last, 2);
    assert_eq!(annot.cdss[last].best_prev, Some(0));
    assert_eq!(
        annot.cdss[last].sum_len,
        annot.cdss[0].size_effective() + annot.cdss[last].size_effective()
    );
    assert!(annot.cdss[1].prevs.is_empty());
}

#[test]
fn test_interval() {
    let i = Interval::new(10, 20, 1);
    assert_eq!(i.len(), 10);
    assert!(!i.empty());
    assert!(i.valid());

    let j = Interval::new(12, 18, 1);
    assert!(i.contains(&j));
    assert!(i.overlaps(&j));
    assert!(!j.contains(&i));
}

#[test]
fn interval_rest_matches_upstream_strand_logic() {
    let plus = Interval::new(10, 20, 1);
    assert_eq!(plus.rest(100, true), 10);
    assert_eq!(plus.rest(100, false), 80);

    let minus = Interval::new(10, 20, -1);
    assert_eq!(minus.rest(100, true), 80);
    assert_eq!(minus.rest(100, false), 10);
}

#[test]
fn test_strand2char() {
    assert_eq!(strand2char(1), '+');
    assert_eq!(strand2char(-1), '-');
    assert_eq!(strand2char(0), '?');
}

#[test]
fn disruption_reportable_excludes_none_and_smooth_only() {
    let smooth = Disruption {
        prev_hsp_idx: Some(0),
        next_hsp_idx: Some(1),
        prev_start: 0,
        next_stop: 0,
        prev_slen: 100,
        intron: false,
        s_stop_codon: false,
        q_interval: Interval::new(10, 10, 1),
        s_interval: Interval::new(20, 20, 1),
    };
    let frameshift = Disruption {
        s_interval: Interval::new(20, 22, 1),
        ..smooth.clone()
    };
    let deletion = Disruption {
        q_interval: Interval::new(10, 12, 1),
        s_interval: Interval::new(20, 23, 1),
        ..smooth.clone()
    };
    let insertion = Disruption {
        s_interval: Interval::new(20, 23, 1),
        ..smooth.clone()
    };

    assert_eq!(smooth.disruption_type(), DisruptionType::Smooth);
    assert!(!smooth.reportable());
    assert_eq!(frameshift.disruption_type(), DisruptionType::Frameshift);
    assert!(frameshift.reportable());
    assert_eq!(deletion.disruption_type(), DisruptionType::Deletion);
    assert!(deletion.reportable());
    assert_eq!(insertion.disruption_type(), DisruptionType::Insertion);
    assert!(insertion.reportable());
}

#[test]
fn same_hsp_non_codon_disruption_is_not_subject_stop_codon() {
    let disruption = Disruption {
        prev_hsp_idx: Some(0),
        next_hsp_idx: Some(0),
        prev_start: 2,
        next_stop: 4,
        prev_slen: 100,
        intron: false,
        s_stop_codon: false,
        q_interval: Interval::new(10, 12, 0),
        s_interval: Interval::new(30, 36, 1),
    };

    assert_eq!(disruption.disruption_type(), DisruptionType::Deletion);
    assert!(!disruption.s_stop_codon());
    assert_eq!(disruption.genesymbol_raw(), "del_10_12_30_36_1");

    let mut text = Vec::new();
    disruption.save_text(&mut text).unwrap();
    assert_eq!(String::from_utf8(text).unwrap(), "Disruption:11-12:31-36/+");
}

#[test]
fn same_hsp_broad_slice_containing_stop_is_subject_stop_codon() {
    let disruption = Disruption {
        prev_hsp_idx: Some(0),
        next_hsp_idx: Some(0),
        prev_start: 2,
        next_stop: 5,
        prev_slen: 100,
        intron: false,
        s_stop_codon: true,
        q_interval: Interval::new(10, 12, 0),
        s_interval: Interval::new(30, 39, 1),
    };

    assert_eq!(disruption.disruption_type(), DisruptionType::Deletion);
    assert!(disruption.s_stop_codon());
    assert_eq!(disruption.genesymbol_raw(), "del_10_12_30_39_1_STOP");

    let mut text = Vec::new();
    disruption.save_text(&mut text).unwrap();
    assert_eq!(
        String::from_utf8(text).unwrap(),
        "Disruption:11-12:31-39/+/*"
    );
}

#[test]
fn reverse_strand_insertion_genesymbol_extends_logical_start_like_upstream() {
    let disruption = Disruption {
        prev_hsp_idx: Some(0),
        next_hsp_idx: Some(1),
        prev_start: 3,
        next_stop: 4,
        prev_slen: 100,
        intron: true,
        s_stop_codon: false,
        q_interval: Interval::new(5, 5, 0),
        s_interval: Interval::new(30, 33, -1),
    };

    assert_eq!(disruption.disruption_type(), DisruptionType::Insertion);
    assert_eq!(disruption.genesymbol_raw(), "ins_4_5_30_36_0");
}

#[test]
fn plus_strand_frameshift_genesymbol_extends_subject_stop_to_prev_slen_like_upstream() {
    let disruption = Disruption {
        prev_hsp_idx: Some(0),
        next_hsp_idx: Some(1),
        prev_start: 3,
        next_stop: 4,
        prev_slen: 100,
        intron: true,
        s_stop_codon: false,
        q_interval: Interval::new(5, 6, 0),
        s_interval: Interval::new(30, 32, 1),
    };

    assert_eq!(disruption.disruption_type(), DisruptionType::Frameshift);
    assert_eq!(disruption.genesymbol_raw(), "fs_5_6_30_100_1");
}

#[test]
fn merge_blastx_chain_records_inter_hsp_disruption() {
    let left = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t3\t6\t1\t9\t30\tAAA\tAAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let right = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t5\t6\t6\t13\t18\t30\tCC\tCC"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    let merged = Hsp::merge_blastx_chain(&[left, right]).unwrap();

    assert!(merged.merged);
    assert_eq!(merged.q_int, Interval::new(0, 6, 0));
    assert_eq!(merged.s_int, Interval::new(0, 18, 1));
    assert_eq!(merged.qseq, "AAACC");
    assert_eq!(merged.sseq, "AAACC");
    assert_eq!(merged.disrs.len(), 1);
    assert_eq!(merged.disrs[0].q_interval, Interval::new(3, 4, 0));
    assert_eq!(merged.disrs[0].s_interval, Interval::new(9, 12, 1));
    assert!(merged.disrs[0].intron);
    assert_eq!(merged.disrs[0].disruption_type(), DisruptionType::Deletion);
}

#[test]
fn exon_arcable_and_intron_disruption_match_blastx_gap() {
    let left = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t3\t6\t1\t9\t30\tAAA\tAAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let right = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t5\t6\t6\t13\t18\t30\tCC\tCC"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    let prev_exon = Exon::new(0, 0, left.length);
    let next_exon = Exon::new(1, 0, right.length);
    assert!(prev_exon.arcable(&next_exon, &[&left, &right], true));
    assert!(prev_exon.arcable(&next_exon, &[&left, &right], false));

    let hsps = vec![left, right];
    let intron = Intron::new(
        Exon::new(0, 0, hsps[0].length),
        Exon::new(1, 0, hsps[1].length),
        &hsps,
        None,
    );
    let disruption = intron.disruption(&hsps);
    assert_eq!(disruption.prev_hsp_idx, Some(0));
    assert_eq!(disruption.next_hsp_idx, Some(1));
    assert_eq!(disruption.prev_start, 3);
    assert_eq!(disruption.next_stop, 0);
    assert_eq!(disruption.q_interval, Interval::new(3, 4, 0));
    assert_eq!(disruption.s_interval, Interval::new(9, 12, 1));
    assert!(disruption.intron);
    assert_eq!(disruption.disruption_type(), DisruptionType::Deletion);
}

#[test]
fn exon_arcable_uses_exon_local_rules_not_whole_hsp_rules() {
    let hsp = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t6\t6\t1\t18\t30\tAAAAAA\tAAAAAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let other = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t7\t9\t9\t19\t27\t30\tCCC\tCCC"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    let left = Exon::new(0, 0, 3);
    let adjacent = Exon::new(0, 3, 6);
    assert!(left.arcable(&adjacent, &[&hsp, &other], true));

    let mut insertion = Exon::new(1, 0, other.length);
    insertion.is_insertion = true;
    assert!(!left.arcable(&insertion, &[&hsp, &other], true));
}

#[test]
fn intron_disruption_uses_original_minus_strand_orientation() {
    let left = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t3\t6\t30\t22\t40\tAAA\tAAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let right = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t5\t6\t6\t18\t13\t40\tCC\tCC"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    let prev_exon = Exon::new(0, 0, left.length);
    let next_exon = Exon::new(1, 0, right.length);
    assert!(prev_exon.arcable(&next_exon, &[&left, &right], true));

    let hsps = vec![left, right];
    let intron = Intron::new(
        Exon::new(0, 0, hsps[0].length),
        Exon::new(1, 0, hsps[1].length),
        &hsps,
        None,
    );
    let disruption = intron.disruption(&hsps);
    assert_eq!(disruption.q_interval, Interval::new(3, 4, 0));
    assert_eq!(disruption.s_interval, Interval::new(18, 21, -1));
    assert_eq!(disruption.disruption_type(), DisruptionType::Deletion);
}

#[test]
fn hsp2exons_keeps_uniform_high_density_hsp_as_one_exon() {
    let hsp = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t12\t12\t1\t36\t60\tAAAAAAAAAAAA\tAAAAAAAAAAAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    let exons = hsp2exons(0, &hsp, None);

    assert_eq!(exons.len(), 1);
    assert_eq!(exons[0].hsp_idx, 0);
    assert_eq!(exons[0].start, 0);
    assert_eq!(exons[0].stop, hsp.length);
    assert!(!exons[0].is_insertion);
    assert_eq!(exons[0].q_interval(&[hsp.clone()]), Interval::new(0, 12, 0));
    assert_eq!(exons[0].s_interval(&[hsp]), Interval::new(0, 36, 1));
}

#[test]
fn hsp2exons_marks_low_density_stretch_as_insertion_exon() {
    let hsp = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t24\t24\t1\t72\t90\tAAAAAAAAAFFFFFFAAAAAAAAA\tAAAAAAAAAYYYYYYAAAAAAAAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    let exons = hsp2exons(0, &hsp, None);

    assert!(exons.iter().any(|exon| exon.is_insertion));
    assert!(exons.iter().any(|exon| !exon.is_insertion));
    assert_eq!(exons.first().unwrap().stop, hsp.length);
    assert_eq!(exons.last().unwrap().start, 0);
    for exon in &exons {
        assert!(exon.start < exon.stop);
        assert_eq!(exon.hsp_idx, 0);
    }
}

#[test]
fn merge_blastx_chain_trims_overlapping_hsps_at_rightmost_split() {
    let left = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t5\t10\t1\t15\t60\tABCDE\tABCDE"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let right = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t4\t10\t10\t16\t36\t60\tDEFGHIJ\tDEFGHIJ"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    assert!(left.can_merge_blastx_with(&right));
    let merged = Hsp::merge_blastx_chain(&[left, right]).unwrap();

    assert!(merged.merged);
    assert_eq!(merged.q_int, Interval::new(0, 10, 0));
    assert_eq!(merged.s_int, Interval::new(0, 36, 1));
    assert_eq!(merged.qseq, "ABCDEFGHIJ");
    assert_eq!(merged.sseq, "ABCDEFGHIJ");
    assert_eq!(merged.disrs.len(), 1);
    assert_eq!(merged.disrs[0].prev_start, 5);
    assert_eq!(merged.disrs[0].next_stop, 2);
    assert_eq!(merged.disrs[0].q_interval, Interval::new(5, 5, 0));
    assert_eq!(merged.disrs[0].s_interval, Interval::new(15, 21, 1));
    assert_eq!(merged.disrs[0].disruption_type(), DisruptionType::Insertion);
}

#[test]
fn blastx_hsp_records_subject_frame_like_reference() {
    let plus = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t3\t6\t1\t9\t30\tAAA\tAAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    plus.qc();
    assert_eq!(plus.sframe, 1);

    let minus = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t2\t6\t18\t13\t30\tAA\tAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    minus.qc();
    assert_eq!(minus.sframe, -1);
}

#[test]
fn hsp_q_better_eq_matches_original_containment_and_tie_order() {
    let broad = Hsp::from_blast_line(
        concat!(
            "refB|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t6\t10\t1\t18\t100\tAAAAAA\tAAAAAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let contained_worse_identity = Hsp::from_blast_line(
        concat!(
            "refA|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t2\t5\t10\t4\t15\t100\tAAAA\tAATA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let outside = Hsp::from_blast_line(
        concat!(
            "refA|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t7\t10\t1\t21\t100\tAAAAAAA\tAAAAAAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let same_score_smaller_query = Hsp::from_blast_line(
        concat!(
            "refA|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t6\t10\t1\t18\t100\tAAAAAA\tAAAAAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    assert!(broad.q_better_eq(&contained_worse_identity));
    assert!(!contained_worse_identity.q_better_eq(&broad));
    assert!(!broad.q_better_eq(&outside));
    assert!(same_score_smaller_query.q_better_eq(&broad));
    assert!(!broad.q_better_eq(&same_score_smaller_query));
}

#[test]
fn hsp_s_start_global_uses_original_strand_formula() {
    let plus = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t3\t5\t10\t20\t28\t100\tAAA\tAAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let minus = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t3\t5\t10\t28\t20\t100\tAAA\tAAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    assert_eq!(plus.s_start_global(), 13);
    assert_eq!(minus.s_start_global(), 34);
}

#[test]
fn blastx_merge_rejects_subject_gaps_past_intron_max() {
    let left = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t3\t6\t1\t9\t8000\tAAA\tAAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let right = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t5\t6\t6\t5011\t5016\t8000\tCC\tCC"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    assert!(!left.can_merge_blastx_with(&right));
    assert!(Hsp::merge_blastx_chain(&[left, right]).is_none());
}

#[test]
fn hsp_merge_get_returns_best_scoring_chain_then_empty() {
    let weak = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "other_contig\t1\t2\t10\t1\t6\t100\tAA\tAT"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let left = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t3\t5\t10\t10\t18\t100\tCCC\tCCC"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let right = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t6\t8\t10\t25\t33\t100\tDDD\tDDD"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    let mut merge = HspMerge::new(vec![weak, left, right], None, 10.0, true);

    let (best, orig_idx, score) = merge.get();
    assert_eq!(orig_idx, Some(1));
    assert!(best.merged);
    assert_eq!(best.q_int, Interval::new(2, 8, 0));
    assert_eq!(best.s_int, Interval::new(9, 33, 1));
    assert_eq!(best.qseq, "CCCDDD");
    assert_eq!(best.sseq, "CCCDDD");
    assert!(score > 0.0);

    let (next, orig_idx, _) = merge.get();
    assert_eq!(orig_idx, Some(0));
    assert!(next.merged);
    assert_eq!(next.qseq, "AA");

    let (empty, orig_idx, score) = merge.get();
    assert!(orig_idx.is_none());
    assert!(empty.empty());
    assert!(score.is_infinite() && score.is_sign_negative());
}

#[test]
fn hsp_merge_builds_arcs_independent_of_input_order() {
    let left = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t3\t10\t1\t9\t100\tCCC\tCCC"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let right = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t4\t6\t10\t25\t33\t100\tDDD\tDDD"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    let mut merge = HspMerge::new(vec![right, left], None, 1.0, true);
    let (best, orig_idx, score) = merge.get();

    assert_eq!(orig_idx, Some(1));
    assert!(best.merged);
    assert_eq!(best.q_int, Interval::new(0, 6, 0));
    assert_eq!(best.s_int, Interval::new(0, 33, 1));
    assert_eq!(best.qseq, "CCCDDD");
    assert_eq!(best.sseq, "CCCDDD");
    assert!(score > 3.0);
}

#[test]
fn hsp_merge_tail_skips_non_arcable_intervening_hsp() {
    let left = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t3\t10\t1\t9\t100\tCCC\tCCC"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let blocker = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "other_contig\t4\t4\t10\t10\t12\t100\tA\tA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let right = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t4\t6\t10\t25\t33\t100\tDDD\tDDD"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    let mut merge = HspMerge::new(vec![left, blocker, right], None, 1.0, true);

    let (best, orig_idx, score) = merge.get();

    assert_eq!(orig_idx, Some(0));
    assert!(best.merged);
    assert_eq!(best.q_int, Interval::new(0, 6, 0));
    assert_eq!(best.s_int, Interval::new(0, 33, 1));
    assert_eq!(best.qseq, "CCCDDD");
    assert_eq!(best.sseq, "CCCDDD");
    assert!(score > 3.0);
}

#[test]
fn exon_merge_tail_uses_explicit_intron_edges() {
    let left = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t3\t10\t1\t9\t100\tCCC\tCCC"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let right = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t4\t6\t10\t25\t33\t100\tDDD\tDDD"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let hsps = vec![left.clone(), right.clone()];
    let exons = vec![Exon::new(0, 0, left.length), Exon::new(1, 0, right.length)];
    let introns = vec![Intron::new(exons[0].clone(), exons[1].clone(), &hsps, None)];

    let mut no_edge_memo = vec![None; exons.len()];
    let (no_edge_score, no_edge_indices) =
        exons[0].merge_tail(&exons, &[], &hsps, &mut no_edge_memo, true, 1.0, None);
    assert_eq!(
        no_edge_indices
            .iter()
            .map(|exon| exon.hsp_idx)
            .collect::<Vec<_>>(),
        vec![0]
    );

    let mut edge_memo = vec![None; exons.len()];
    let (edge_score, edge_indices) =
        exons[0].merge_tail(&exons, &introns, &hsps, &mut edge_memo, true, 1.0, None);
    assert_eq!(
        edge_indices
            .iter()
            .map(|exon| exon.hsp_idx)
            .collect::<Vec<_>>(),
        vec![0, 1]
    );
    assert!(edge_score > no_edge_score);
}

#[test]
fn hsp_merge_keeps_insertion_exons_in_graph_scoring() {
    let mut matrix = tempfile::NamedTempFile::new().unwrap();
    write!(matrix, "  {}\n", SubstMat::CHARS).unwrap();
    for c1 in SubstMat::CHARS.chars() {
        write!(matrix, "{c1}").unwrap();
        for c2 in SubstMat::CHARS.chars() {
            let score = if c1 == 'F' && c2 == 'Y' {
                20
            } else if c1 == c2 {
                1
            } else {
                -10
            };
            write!(matrix, " {score}").unwrap();
        }
        writeln!(matrix).unwrap();
    }
    let mat = SubstMat::new(matrix.path().to_str().unwrap()).unwrap();

    let hsp = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t24\t24\t1\t72\t90\tAAAAAAAAAFFFFFFAAAAAAAAA\tAAAAAAAAAYYYYYYAAAAAAAAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    let exons = hsp2exons(0, &hsp, None);
    assert!(exons.iter().any(|exon| exon.is_insertion));
    assert!(exons
        .iter()
        .any(|exon| exon.is_insertion && &hsp.qseq[exon.start..exon.stop] == "FFFFFF"));

    let mut merge = HspMerge::new(vec![hsp], Some(&mat), 10.0, true);
    let (merged, orig_idx, score) = merge.get();

    assert_eq!(orig_idx, Some(0));
    assert!(merged.merged);
    assert_eq!(merged.qseq, "AAAAAAAAAFFFFFFAAAAAAAAA");
    assert_eq!(merged.sseq, "AAAAAAAAAYYYYYYAAAAAAAAA");
    assert_eq!(merged.q_int, Interval::new(0, 24, 0));
    assert_eq!(merged.s_int, Interval::new(0, 72, 1));
    assert!(score > 120.0);
}

#[test]
fn hsp_merge_get_keeps_unselected_exons_from_same_hsp() {
    let mut matrix = tempfile::NamedTempFile::new().unwrap();
    write!(matrix, "  {}\n", SubstMat::CHARS).unwrap();
    for c1 in SubstMat::CHARS.chars() {
        write!(matrix, "{c1}").unwrap();
        for c2 in SubstMat::CHARS.chars() {
            let score = if c1 == c2 { 1 } else { -10 };
            write!(matrix, " {score}").unwrap();
        }
        writeln!(matrix).unwrap();
    }
    let mat = SubstMat::new(matrix.path().to_str().unwrap()).unwrap();

    let hsp = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t24\t24\t1\t72\t90\tAAAAAAAAAFFFFFFAAAAAAAAA\tAAAAAAAAAYYYYYYAAAAAAAAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    let mut merge = HspMerge::new(vec![hsp], Some(&mat), 10.0, true);
    let (first, first_idx, first_score) = merge.get();
    let (second, second_idx, second_score) = merge.get();
    let (empty, empty_idx, empty_score) = merge.get();

    assert_eq!(first_idx, Some(0));
    assert_eq!(second_idx, Some(0));
    assert_eq!(first.qseq, "AAAAAAAAA");
    assert_eq!(second.qseq, "AAAAAAAAA");
    assert_eq!(first_score, 9.0);
    assert_eq!(second_score, 9.0);

    let mut intervals = vec![(first.q_int, first.s_int), (second.q_int, second.s_int)];
    intervals.sort();
    assert_eq!(
        intervals,
        vec![
            (Interval::new(0, 9, 0), Interval::new(0, 27, 1)),
            (Interval::new(15, 24, 0), Interval::new(45, 72, 1)),
        ]
    );

    assert!(empty_idx.is_none());
    assert!(empty.empty());
    assert!(empty_score.is_infinite() && empty_score.is_sign_negative());
}

#[test]
fn hsp_merge_preserves_exon_stop_codon_disruptions() {
    let hsp = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t3\t3\t1\t9\t9\tAAA\tA*A"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    let exons = hsp2exons(0, &hsp, None);
    assert_eq!(exons.len(), 1);
    assert_eq!(exons[0].disrs.len(), 1);
    assert_eq!(
        exons[0].disrs[0].disruption_type(),
        DisruptionType::Deletion
    );
    assert_eq!(exons[0].disrs[0].q_interval, Interval::new(1, 2, 0));
    assert_eq!(exons[0].disrs[0].s_interval, Interval::new(3, 6, 1));
    assert!(!exons[0].disrs[0].intron);
    assert!(exons[0].disrs[0].s_stop_codon());
    assert_eq!(exons[0].disrs[0].genesymbol_raw(), "del_1_2_3_6_1_STOP");
    let mut disr_text = Vec::new();
    exons[0].disrs[0].save_text(&mut disr_text).unwrap();
    assert_eq!(
        String::from_utf8(disr_text).unwrap(),
        "Disruption:2-2:4-6/+/*"
    );

    let mut merge = HspMerge::new(vec![hsp], None, 10.0, true);
    let (merged, orig_idx, score) = merge.get();

    assert_eq!(orig_idx, Some(0));
    assert!(merged.merged);
    assert_eq!(merged.disrs.len(), 1);
    assert_eq!(merged.disrs[0].disruption_type(), DisruptionType::Deletion);
    assert_eq!(merged.disrs[0].q_interval, Interval::new(1, 2, 0));
    assert_eq!(merged.disrs[0].s_interval, Interval::new(3, 6, 1));
    assert!(!merged.disrs[0].intron);
    assert!(merged.disrs[0].s_stop_codon());
    assert_eq!(merged.disrs[0].genesymbol_raw(), "del_1_2_3_6_1_STOP");
    assert!(score > 0.0);
}

#[test]
fn hsp_merge_clips_stop_codon_disruption_removed_by_overlap_split() {
    let left = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t5\t12\t1\t15\t60\tAAAAA\tAAAAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let right = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t4\t12\t12\t16\t42\t60\tAAAAAAAAA\tA*AAAAAAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    let right_exons = hsp2exons(1, &right, None);
    assert!(right_exons.iter().any(|exon| !exon.disrs.is_empty()));

    let mut merge = HspMerge::new(vec![left, right], None, 0.0, true);
    let (merged, orig_idx, score) = merge.get();

    assert_eq!(orig_idx, Some(0));
    assert_eq!(merged.qseq, "AAAAAAAAAAAA");
    assert_eq!(merged.sseq, "AAAAAAAAAAAA");
    assert!(!merged.s_internal_stop);
    assert!(merged.disrs.iter().all(|disr| !disr.s_stop_codon()));
    assert!(score > 0.0);
}

#[test]
fn hsp_merge_uses_supplied_substitution_matrix_for_exon_scoring() {
    let mut matrix = tempfile::NamedTempFile::new().unwrap();
    write!(matrix, "  {}\n", SubstMat::CHARS).unwrap();
    for c1 in SubstMat::CHARS.chars() {
        write!(matrix, "{c1}").unwrap();
        for c2 in SubstMat::CHARS.chars() {
            let score = if c1 == c2 {
                if c1 == 'C' {
                    20
                } else {
                    1
                }
            } else {
                -10
            };
            write!(matrix, " {score}").unwrap();
        }
        writeln!(matrix).unwrap();
    }
    let mat = SubstMat::new(matrix.path().to_str().unwrap()).unwrap();

    let weak_default = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t2\t10\t1\t6\t100\tAA\tAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let strong_matrix = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "other_contig\t3\t4\t10\t10\t15\t100\tCC\tCC"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    let mut default_merge = HspMerge::new(
        vec![weak_default.clone(), strong_matrix.clone()],
        None,
        1.0,
        true,
    );
    let (default_best, default_idx, default_score) = default_merge.get();
    assert_eq!(default_idx, Some(0));
    assert_eq!(default_best.qseq, "AA");
    assert_eq!(default_score, 2.0);

    let mut matrix_merge = HspMerge::new(vec![weak_default, strong_matrix], Some(&mat), 1.0, true);
    let (matrix_best, matrix_idx, matrix_score) = matrix_merge.get();
    assert_eq!(matrix_idx, Some(1));
    assert_eq!(matrix_best.qseq, "CC");
    assert_eq!(matrix_score, 40.0);
}

#[test]
fn hsp_merge_skips_single_residue_hsps_like_reference_graph() {
    let one_residue = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t1\t10\t1\t3\t100\tA\tA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let two_residue = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "other_contig\t2\t3\t10\t10\t15\t100\tCC\tCC"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    assert!(hsp2exons(0, &one_residue, None).is_empty());

    let mut merge = HspMerge::new(vec![one_residue, two_residue], None, 1.0, true);
    let (best, orig_idx, score) = merge.get();
    let (empty, empty_idx, empty_score) = merge.get();

    assert_eq!(orig_idx, Some(1));
    assert_eq!(best.qseq, "CC");
    assert_eq!(score, 2.0);
    assert!(empty_idx.is_none());
    assert!(empty.empty());
    assert!(empty_score.is_infinite() && empty_score.is_sign_negative());
}

#[test]
fn hsp_merge_subtracts_overlap_split_score_like_reference() {
    let left = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t1\t10\t11\t1\t30\t100\tAAAAAAAAAA\tAAAAAAAAAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let overlapping_right = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "contig\t2\t11\t11\t31\t60\t100\tAAAAAAAAAA\tAAAAAAAAAA"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();
    let standalone = Hsp::from_blast_line(
        concat!(
            "ref|1|1|fam|gene|resistance|2|subclass|class|product\t",
            "other_contig\t1\t12\t12\t1\t36\t100\tCCCCCCCCCCCC\tCCCCCCCCCCCC"
        ),
        true,
        false,
        true,
        false,
        false,
    )
    .unwrap();

    let mut merge = HspMerge::new(vec![left, overlapping_right, standalone], None, 0.0, true);
    let (best, orig_idx, score) = merge.get();

    assert_eq!(orig_idx, Some(2));
    assert_eq!(best.qseq, "CCCCCCCCCCCC");
    assert_eq!(score, 12.0);
}
