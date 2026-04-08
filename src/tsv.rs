use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

use anyhow::{bail, Result};

// --- TsvOut ---

/// Tab-separated output writer. Matches the C++ TsvOut class.
pub struct TsvOut<'a> {
    os: Option<&'a mut dyn Write>,
    lines: usize,
    fields_max: usize,
    fields: usize,
    pub use_pound: bool,
}

impl<'a> TsvOut<'a> {
    pub fn new(os: Option<&'a mut dyn Write>) -> Self {
        TsvOut {
            os,
            lines: 0,
            fields_max: 0,
            fields: 0,
            use_pound: true,
        }
    }

    pub fn live(&self) -> bool {
        self.os.is_some()
    }

    pub fn empty(&self) -> bool {
        self.lines == 0 && self.fields_max == 0 && self.fields == 0
    }

    pub fn write_field<T: fmt::Display>(&mut self, field: &T) -> Result<()> {
        if let Some(ref mut os) = self.os {
            if self.lines > 0 && self.fields >= self.fields_max {
                bail!("TsvOut: fields_max = {}", self.fields_max);
            }
            if self.fields > 0 {
                write!(os, "\t")?;
            } else if self.lines == 0 && self.use_pound {
                write!(os, "#")?;
            }
            write!(os, "{}", field)?;
            self.fields += 1;
        }
        Ok(())
    }

    pub fn new_line(&mut self) -> Result<()> {
        if let Some(ref mut os) = self.os {
            writeln!(os)?;
            if self.lines == 0 {
                self.fields_max = self.fields;
            }
            self.lines += 1;
            if self.fields != self.fields_max {
                bail!(
                    "TsvOut: fields_max = {}, but fields = {}",
                    self.fields_max,
                    self.fields
                );
            }
            self.fields = 0;
        }
        Ok(())
    }
}

// --- TextTable ---

/// Column header metadata
#[derive(Debug, Clone)]
pub struct Header {
    pub name: String,
    pub len_max: usize,
    pub numeric: bool,
    pub scientific: bool,
    pub decimals: usize,
    pub null: bool,
}

impl Header {
    pub fn new(name: &str) -> Self {
        Header {
            name: name.to_string(),
            len_max: 0,
            numeric: true,
            scientific: false,
            decimals: 0,
            null: false,
        }
    }
}

/// Tab-separated value table with header.
/// Matches the C++ TextTable class.
pub struct TextTable {
    pub name: String,
    pub pound: bool,
    pub save_header: bool,
    pub header: Vec<Header>,
    pub rows: Vec<Vec<String>>,
}

impl TextTable {
    /// Parse a TSV file
    pub fn from_file(fname: &str) -> Result<Self> {
        let file = File::open(fname)?;
        let reader = BufReader::new(file);
        let mut header: Vec<Header> = Vec::new();
        let mut rows: Vec<Vec<String>> = Vec::new();
        let mut pound = false;
        let mut header_found = false;

        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim().to_string();
            if trimmed.is_empty() {
                continue;
            }

            if !header_found {
                // First non-empty line is the header
                let mut hdr_line = trimmed.as_str();
                if hdr_line.starts_with('#') {
                    pound = true;
                    hdr_line = &hdr_line[1..];
                }
                for name in hdr_line.split('\t') {
                    header.push(Header::new(name.trim()));
                }
                header_found = true;
                continue;
            }

            // Skip comment lines
            if trimmed.starts_with('#') {
                continue;
            }

            let fields: Vec<String> = trimmed.split('\t').map(|s| s.trim().to_string()).collect();
            // Pad with empty strings if needed
            let mut row = fields;
            while row.len() < header.len() {
                row.push(String::new());
            }
            rows.push(row);
        }

        Ok(TextTable {
            name: fname.to_string(),
            pound,
            save_header: true,
            header,
            rows,
        })
    }

    /// Create an empty table with the given header
    pub fn with_header(header: Vec<Header>) -> Self {
        TextTable {
            name: String::new(),
            pound: false,
            save_header: true,
            header,
            rows: Vec::new(),
        }
    }

    /// Find column index by name, returns None if not found
    pub fn col2num(&self, column_name: &str) -> Option<usize> {
        self.header.iter().position(|h| h.name == column_name)
    }

    /// Find column index by name, returns error if not found
    pub fn col2num_required(&self, column_name: &str) -> Result<usize> {
        self.col2num(column_name).ok_or_else(|| {
            anyhow::anyhow!(
                "Table has no column \"{}\"\nIn table file: {}",
                column_name,
                self.name
            )
        })
    }

    /// Check if table has a column
    pub fn has_column(&self, column_name: &str) -> bool {
        self.col2num(column_name).is_some()
    }

    /// Sort rows by given columns
    pub fn sort(&mut self, by: &[&str]) -> Result<()> {
        let col_nums: Vec<usize> = by
            .iter()
            .map(|name| self.col2num_required(name))
            .collect::<Result<Vec<_>>>()?;

        self.rows.sort_by(|a, b| {
            for &col in &col_nums {
                let cmp = a.get(col).map(|s| s.as_str()).unwrap_or("")
                    .cmp(b.get(col).map(|s| s.as_str()).unwrap_or(""));
                if cmp != std::cmp::Ordering::Equal {
                    return cmp;
                }
            }
            std::cmp::Ordering::Equal
        });

        Ok(())
    }

    /// Filter to only the specified columns (in order)
    pub fn filter_columns(&mut self, new_column_names: &[&str]) -> Result<()> {
        let col_nums: Vec<usize> = new_column_names
            .iter()
            .map(|name| self.col2num_required(name))
            .collect::<Result<Vec<_>>>()?;

        let new_header: Vec<Header> = col_nums
            .iter()
            .map(|&i| self.header[i].clone())
            .collect();

        let new_rows: Vec<Vec<String>> = self
            .rows
            .iter()
            .map(|row| {
                col_nums
                    .iter()
                    .map(|&i| row.get(i).cloned().unwrap_or_default())
                    .collect()
            })
            .collect();

        self.header = new_header;
        self.rows = new_rows;

        Ok(())
    }

    /// Write to a stream
    pub fn write_to(&self, out: &mut dyn Write) -> Result<()> {
        if self.save_header {
            if self.pound {
                write!(out, "#")?;
            }
            let names: Vec<&str> = self.header.iter().map(|h| h.name.as_str()).collect();
            writeln!(out, "{}", names.join("\t"))?;
        }
        for row in &self.rows {
            writeln!(out, "{}", row.join("\t"))?;
        }
        Ok(())
    }

    /// Build a key index on given columns
    pub fn build_key(&self, columns: &[&str]) -> Result<TextTableKey> {
        let col_nums: Vec<usize> = columns
            .iter()
            .map(|name| self.col2num_required(name))
            .collect::<Result<Vec<_>>>()?;

        let mut data: HashMap<Vec<String>, usize> = HashMap::new();
        for (row_num, row) in self.rows.iter().enumerate() {
            let key: Vec<String> = col_nums
                .iter()
                .map(|&i| row.get(i).cloned().unwrap_or_default())
                .collect();
            data.insert(key, row_num);
        }

        Ok(TextTableKey { col_nums, data })
    }
}

/// Key index for fast row lookup
pub struct TextTableKey {
    pub col_nums: Vec<usize>,
    pub data: HashMap<Vec<String>, usize>,
}

impl TextTableKey {
    pub fn find(&self, values: &[String]) -> Option<usize> {
        self.data.get(values).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tsv_out() {
        let mut buf = Vec::new();
        {
            let mut tsv = TsvOut::new(Some(&mut buf));
            tsv.write_field(&"col1").unwrap();
            tsv.write_field(&"col2").unwrap();
            tsv.write_field(&"col3").unwrap();
            tsv.new_line().unwrap();
            tsv.write_field(&"a").unwrap();
            tsv.write_field(&"b").unwrap();
            tsv.write_field(&"c").unwrap();
            tsv.new_line().unwrap();
        }
        let output = String::from_utf8(buf).unwrap();
        assert_eq!(output, "#col1\tcol2\tcol3\na\tb\tc\n");
    }

    #[test]
    fn test_text_table_basic() {
        let table = TextTable::with_header(vec![
            Header::new("Name"),
            Header::new("Value"),
        ]);
        assert_eq!(table.header.len(), 2);
        assert_eq!(table.col2num("Name"), Some(0));
        assert_eq!(table.col2num("Value"), Some(1));
        assert_eq!(table.col2num("Missing"), None);
    }
}
