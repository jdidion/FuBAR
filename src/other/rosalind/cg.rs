use clap::Parser;
use noodles::fasta;
use std::{
    fs::File,
    io::{self, BufReader},
    path::PathBuf,
};

#[derive(clap::Parser)]
struct Args {
    #[arg(short = 'i')]
    input: PathBuf,
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    let mut reader = File::open(args.input).map(BufReader::new).map(fasta::Reader::new)?;

    if let Some((name, gc_frac)) = reader
        .records()
        .map(|result| -> io::Result<(String, f64)> {
            let record = result?;
            let sequence = record.sequence();
            let seq_len = sequence.len() as f64;
            let gc_frac = if seq_len == 0.0 {
                0.0
            } else {
                sequence.as_ref().iter().filter(|base| base == &&b'C' || base == &&b'G').count()
                    as f64
                    / seq_len
            };
            Ok((record.name().to_owned(), gc_frac))
        })
        .collect::<io::Result<Vec<_>>>()?
        .into_iter()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
    {
        println!("{}\n{}", name, gc_frac);
    } else {
        println!("Empty FASTA file")
    }
    Ok(())
}
