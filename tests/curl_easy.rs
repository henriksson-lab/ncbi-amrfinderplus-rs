use std::fs;

use amrfinder::curl_easy::{get_lib_version, Curl};

#[test]
fn get_lib_version_reads_linked_libcurl_version() {
    assert_ne!(get_lib_version().str(), "0.0.0");
}

#[test]
fn read_collects_file_url_contents() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("body.txt");
    fs::write(&src, "alpha\nbeta\n").unwrap();

    let mut curl = Curl::new().unwrap();
    let body = curl.read(&format!("file://{}", src.display())).unwrap();

    assert_eq!(body, "alpha\nbeta\n");
}

#[test]
fn download_writes_file_url_contents() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("body.txt");
    let dst = dir.path().join("download.txt");
    fs::write(&src, "payload\n").unwrap();

    let mut curl = Curl::new().unwrap();
    curl.download(&format!("file://{}", src.display()), dst.to_str().unwrap())
        .unwrap();

    assert_eq!(fs::read_to_string(dst).unwrap(), "payload\n");
}

#[test]
fn download_rejects_xml_error_document() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("body.xml");
    let dst = dir.path().join("download.xml");
    fs::write(&src, "<?xml version=\"1.0\"?><Error/>\n").unwrap();

    let mut curl = Curl::new().unwrap();
    let err = curl
        .download(&format!("file://{}", src.display()), dst.to_str().unwrap())
        .unwrap_err()
        .to_string();

    assert!(err.contains("Cannot download"));
}

#[test]
fn read_preserves_utf8_file_url_contents() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("utf8.txt");
    fs::write(&src, "cafe: \u{e9}\n").unwrap();

    let mut curl = Curl::new().unwrap();
    let body = curl.read(&format!("file://{}", src.display())).unwrap();

    assert_eq!(body, "cafe: \u{e9}\n");
}

#[test]
fn download_allows_non_utf8_payloads() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("body.bin");
    let dst = dir.path().join("download.bin");
    fs::write(&src, [0xff, 0x00, b'x']).unwrap();

    let mut curl = Curl::new().unwrap();
    curl.download(&format!("file://{}", src.display()), dst.to_str().unwrap())
        .unwrap();

    assert_eq!(fs::read(dst).unwrap(), [0xff, 0x00, b'x']);
}
