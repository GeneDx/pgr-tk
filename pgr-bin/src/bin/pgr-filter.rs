const VERSION_STRING: &'static str = env!("VERSION_STRING");
use clap::{self, IntoApp, Parser};
use flate2::bufread::MultiGzDecoder;
use pgr_db::fasta_io::{FastaReader, FastqStreamReader, SeqRec, FastaStreamReader};
use pgr_db::kmer_filter::KmerFilter;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufReader, Read};

enum GZFastaReader {
    GZFile(FastaReader<BufReader<MultiGzDecoder<BufReader<File>>>>),
    RegularFile(FastaReader<BufReader<BufReader<File>>>),
}

#[derive(Parser, Debug)]
#[clap(name = "pgr-filter")]
#[clap(author, version)]
#[clap(about = "using Cuckoo Filter for Matching Reads To A Reference Set of Sequences", long_about = None)]
struct CmdOptions {
    ref_fasta_path: String,
    #[clap(long, short)]
    query_fastx_path: Option<String>,
    /// k-mer size
    #[clap(long, short, default_value_t = 32)]
    k: usize,
    /// count threshold
    #[clap(long, short, default_value_t = 0.8)]
    threshold: f32,
    #[clap(long)]
    fasta_stdin: bool,
}

fn get_fastx_reader(filepath: String) -> Result<GZFastaReader, std::io::Error> {
    let file = File::open(&filepath)?;
    let mut reader = BufReader::new(file);
    let mut is_gzfile = false;
    {
        let r = reader.by_ref();
        let mut buf = Vec::<u8>::new();
        let _ = r.take(2).read_to_end(&mut buf);
        if buf == [0x1F_u8, 0x8B_u8] {
            log::info!("input file: {} detected as gz-compressed file", filepath);
            is_gzfile = true;
        }
    }
    drop(reader);

    let file = File::open(&filepath)?;
    let reader = BufReader::new(file);
    let gz_buf = BufReader::new(MultiGzDecoder::new(reader));

    let file = File::open(&filepath)?;
    let reader = BufReader::new(file);
    let std_buf = BufReader::new(reader);

    if is_gzfile {
        drop(std_buf);
        Ok(GZFastaReader::GZFile(
            FastaReader::new(gz_buf, &filepath, 256, false).unwrap(),
        ))
    } else {
        drop(gz_buf);
        Ok(GZFastaReader::RegularFile(
            FastaReader::new(std_buf, &filepath, 256, false).unwrap(),
        ))
    }
}

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();
    let mut filter = KmerFilter::with_capacity(args.k, 1_usize << 24);
    let mut add_seqs = |seq_iter: &mut dyn Iterator<Item = io::Result<SeqRec>>| {
        seq_iter.into_iter().for_each(|r| {
            if let Ok(r) = r {
                filter.add_seq_mmers(&r.seq);
            };
        });
    };

    match get_fastx_reader(args.ref_fasta_path)? {
        GZFastaReader::GZFile(reader) => add_seqs(&mut reader.into_iter()),

        GZFastaReader::RegularFile(reader) => add_seqs(&mut reader.into_iter()),
    };

    let check_seqs = |seq_iter: &mut dyn Iterator<Item = io::Result<SeqRec>>| {
        let mut seq_data = Vec::<SeqRec>::new();
        for r in seq_iter {
            if let Ok(r) = r {
                seq_data.push(r);
            };
            if seq_data.len() == 1024 {
                seq_data
                    .par_iter()
                    .map(|r| {
                        let (total, c) = filter.check_seq_mmers(&r.seq);
                        (r.clone(), total, c)
                    })
                    .collect::<Vec<(SeqRec, usize, usize)>>()
                    .iter()
                    .for_each(|(r, total, c)| {
                        if *total > 0 {
                            if (*c as f32) / (*total as f32) > args.threshold {
                                println!(">{} {} {}", String::from_utf8_lossy(&r.id), total, c);
                                println!("{}", String::from_utf8_lossy(&r.seq[..]));
                            }
                        }
                    });
                seq_data.clear();
            }
        }

        seq_data
            .into_par_iter()
            .map(|r| {
                let (total, c) = filter.check_seq_mmers(&r.seq);
                (r, total, c)
            })
            .collect::<Vec<(SeqRec, usize, usize)>>()
            .iter()
            .for_each(|(r, total, c)| {
                if *total > 0 {
                    if (*c as f32) / (*total as f32) > args.threshold {
                        println!(">{} {} {}", String::from_utf8_lossy(&r.id), total, c);
                        println!("{}", String::from_utf8_lossy(&r.seq[..]));
                    }
                }
            });
    };

    if args.query_fastx_path.is_some() {
        match get_fastx_reader(args.query_fastx_path.unwrap())? {
            GZFastaReader::GZFile(reader) => check_seqs(&mut reader.into_iter()),
            GZFastaReader::RegularFile(reader) => check_seqs(&mut reader.into_iter()),
        }
    } else {
        if args.fasta_stdin {
            let reader = FastaStreamReader::new(256);
            check_seqs(&mut reader.into_iter());
        } else {
            let reader = FastqStreamReader::new(256);
            check_seqs(&mut reader.into_iter());
        }
    }

    Ok(())
}
