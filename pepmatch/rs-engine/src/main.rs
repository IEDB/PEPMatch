mod preprocess;
#[path = "match.rs"]
mod matching;

use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        print_usage();
        process::exit(1);
    }

    match args[1].as_str() {
        "preprocess" => {
            if args.len() < 5 {
                eprintln!("Usage: pepmatch-rs preprocess <fasta> <k> <output.pepidx>");
                process::exit(1);
            }
            let k: usize = args[3].parse().unwrap_or_else(|_| {
                eprintln!("Error: k must be a positive integer");
                process::exit(1);
            });
            preprocess::run(&args[2], k, &args[4]);
        }
        "match" => {
            if args.len() < 6 {
                eprintln!("Usage: pepmatch-rs match <pepidx> <peptides> <k> <max_mismatches>");
                process::exit(1);
            }
            let k: usize = args[4].parse().unwrap_or_else(|_| {
                eprintln!("Error: k must be a positive integer");
                process::exit(1);
            });
            let mm: usize = args[5].parse().unwrap_or_else(|_| {
                eprintln!("Error: max_mismatches must be a non-negative integer");
                process::exit(1);
            });
            matching::run(&args[2], &args[3], k, mm);
        }
        _ => {
            eprintln!("Unknown command: {}", args[1]);
            print_usage();
            process::exit(1);
        }
    }
}

fn print_usage() {
    eprintln!("PEPMatch Rust Engine");
    eprintln!();
    eprintln!("Usage:");
    eprintln!("  pepmatch-rs preprocess <fasta> <k> <output.pepidx>");
    eprintln!("  pepmatch-rs match <pepidx> <peptides> <k> <max_mismatches>");
}
