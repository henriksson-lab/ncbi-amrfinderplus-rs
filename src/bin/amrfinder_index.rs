use std::process;

fn main() {
    match amrfinder::amrfinder_index::main(std::env::args_os().collect()) {
        Ok(code) => process::exit(code),
        Err(e) => {
            eprintln!("{}", e);
            process::exit(1);
        }
    }
}
