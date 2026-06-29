use std::path::PathBuf;
use std::process;

use amrfinder::common::{Application, Key};

fn main() {
    let mut app = Application {
        description: "Check the correctness of a FASTA file. Exit with an error if it is incorrect. Print the number of sequences, max. sequence length and total sequence length",
        usage: "Usage: fasta_check <in> [-aa] [-hyphen] [-ambig] [-ambig_max <n>] [-stop_codon] [-len <file>] [-out <file>]",
        positionals: 1,
        needs_arg: false,
        threads_used: false,
        key_args: vec![
            Key {
                name: "aa",
                flag: true,
                default_value: "false",
                acronym: None,
            },
            Key {
                name: "hyphen",
                flag: true,
                default_value: "false",
                acronym: None,
            },
            Key {
                name: "ambig",
                flag: true,
                default_value: "false",
                acronym: None,
            },
            Key {
                name: "ambig_max",
                flag: false,
                default_value: "0",
                acronym: None,
            },
            Key {
                name: "stop_codon",
                flag: true,
                default_value: "false",
                acronym: None,
            },
            Key {
                name: "len",
                flag: false,
                default_value: "",
                acronym: None,
            },
            Key {
                name: "out",
                flag: false,
                default_value: "",
                acronym: None,
            },
        ],
    };
    let run = match app.run(std::env::args_os().collect()) {
        Ok(run) => run,
        Err(code) => process::exit(code),
    };

    let fasta_path = PathBuf::from(&run.positional_args[0]);
    let ambig_max = match run.key_args["ambig_max"].parse::<usize>() {
        Ok(value) => value,
        Err(e) => {
            eprintln!("{}", e);
            process::exit(1);
        }
    };
    let len_path = (!run.key_args["len"].is_empty()).then(|| PathBuf::from(&run.key_args["len"]));
    let out_path = (!run.key_args["out"].is_empty()).then(|| PathBuf::from(&run.key_args["out"]));

    match amrfinder::fasta_check::body(
        &fasta_path,
        run.key_args["aa"] == "true",
        run.key_args["hyphen"] == "true",
        run.key_args["ambig"] == "true",
        ambig_max,
        run.key_args["stop_codon"] == "true",
        len_path.as_deref(),
        out_path.as_deref(),
    ) {
        Ok((num_seqs, max_len, total_len)) => {
            println!("{}", num_seqs);
            println!("{}", max_len);
            println!("{}", total_len);
            if !run.key_args["log"].is_empty() {
                let _ = std::fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(&run.key_args["log"]);
            }
            if !run.key_args["json"].is_empty() {
                if let Err(e) = std::fs::write(&run.key_args["json"], "{}\n") {
                    eprintln!("{}", e);
                    process::exit(1);
                }
            }
        }
        Err(e) => {
            eprintln!("{}", e);
            process::exit(1);
        }
    }
}
