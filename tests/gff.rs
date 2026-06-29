use std::path::PathBuf;

use amrfinder::gff::{Annot, GffType};

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amr")
}

#[test]
fn parse_gff() {
    let path = test_data_dir().join("test_prot.gff");
    if !path.exists() {
        return;
    }
    let annot = Annot::from_gff(path.to_str().unwrap(), GffType::Genbank, false, false);
    assert!(annot.is_ok(), "GFF parse failed: {:?}", annot.err());
    let annot = annot.unwrap();
    assert!(!annot.prot2loci.is_empty());
}

#[test]
fn gff_type_names() {
    for name in GffType::NAMES {
        let result = GffType::name2type(name);
        assert!(result.is_ok(), "Failed to parse GFF type: {}", name);
    }
    assert!(GffType::name2type("invalid").is_err());
}
