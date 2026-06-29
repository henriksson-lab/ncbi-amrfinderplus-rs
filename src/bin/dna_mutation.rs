use std::io;
use std::process;

fn main() {
    let stdout = io::stdout();
    let mut out = stdout.lock();
    match amrfinder::dna_mutation::main(std::env::args_os().collect(), &mut out) {
        Ok(code) => process::exit(code),
        Err(e) => {
            eprintln!("{}", e);
            process::exit(1);
        }
    }
}
