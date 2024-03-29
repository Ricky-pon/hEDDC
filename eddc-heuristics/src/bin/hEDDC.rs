use clap;
use eddc_heuristics::fasta;
use string_decomposer;

use clap::Parser;
use std::path::PathBuf;
#[derive(Parser)]
#[clap(author, version = "0.1", about = "All vs All alirngmnet between tandem repeats", long_about = None)]
struct Args {
    #[clap(short, long, value_parser, value_name = "FASTA")]
    reads: PathBuf,
    #[clap(short, long, value_parser, value_name = "FASTA")]
    units: PathBuf,
}
fn main() -> std::io::Result<()> {
    let args = Args::parse();
    let reads = fasta::parse_into_vec(&args.reads)?;
    let reads: Vec<_> = reads.iter().map(|r| r.seq()).collect();
    let units = fasta::parse_into_vec(&args.units)?;
    let units: Vec<_> = units.iter().map(|u| u.seq()).collect();
    let params = string_decomposer::StringDecompParameters::new(1, -1, -1, -1);
    use rayon::prelude::*;
    let encoded_reads: Vec<Vec<_>> = reads
        .par_iter()
        .map(|read| {
            string_decomposer::DecomposedSeq::new(read, &units, &params)
                .encoding()
                .iter()
                .map(|enc| enc.id())
                .collect()
        })
        .collect();
    let score_matrix = eddc_heuristics::eddc_heuristic_parallel(&encoded_reads, &units);
    for (i, v) in score_matrix.iter().enumerate() {
        let scores: Vec<_> = v
            .iter()
            .enumerate()
            .map(|(j, &score)| {
                let len = ((reads[i].len() * reads[j].len()) as f64).sqrt();
                // Avoid zero div.
                let normalized_score = if 0.000001 < len { score / len } else { score };
                format!("{normalized_score}")
            })
            .collect();
        println!("{}", scores.join(" "));
    }
    Ok(())
}
