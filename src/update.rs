// Database download and update — port of amrfinder_update.cpp
// Replaces libcurl with ureq

use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

use anyhow::Result;

const NCBI_FTP_BASE: &str = "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database";

/// Download a file from a URL to a local path
pub fn download_file(url: &str, dest: &Path) -> Result<()> {
    let body = ureq::get(url).call()?.into_body().read_to_vec()?;
    let mut file = File::create(dest)?;
    file.write_all(&body)?;
    Ok(())
}

/// Read a URL to a string
pub fn read_url(url: &str) -> Result<String> {
    let body = ureq::get(url).call()?.into_body().read_to_string()?;
    Ok(body)
}

/// Get the latest database version from NCBI
pub fn get_latest_version() -> Result<String> {
    let url = format!("{}/latest/VERSION.txt", NCBI_FTP_BASE);
    let version = read_url(&url)?;
    Ok(version.trim().to_string())
}

/// Check if update is needed
pub fn check_update(db_dir: &Path) -> Result<Option<String>> {
    let latest = get_latest_version()?;

    // Check current version
    let version_file = db_dir.join("version.txt");
    if version_file.exists() {
        let current = fs::read_to_string(&version_file)?;
        if current.trim() == latest {
            return Ok(None); // Up to date
        }
    }

    Ok(Some(latest))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_download_url_parsing() {
        // Just verify the URL constant is valid
        assert!(NCBI_FTP_BASE.starts_with("https://"));
    }
}
