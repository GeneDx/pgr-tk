const VERSION_STRING: &'static str = env!("VERSION_STRING");

//use std::path::PathBuf;
use clap::{self, CommandFactory, Parser};

use pgr_db::shmmrutils::ShmmrSpec;
use std::fs::File;
use std::io::{BufWriter, Write};

use pgr_db::seq_db;

#[derive(Parser, Debug)]
#[clap(name = "pgr-seq-smp-count")]
#[clap(author, version)]
#[clap(about = "count shimmer pairs in a shimmer database", long_about = None)]
struct CmdOptions {
    #[clap(long, short)]
    in_fasta: String,
    #[clap(long, short)]
    output_path: String,
    //max_unique_count
    #[clap(long, short, default_value_t = 4)]
    min_count: usize,
    #[clap(long, short, default_value_t = 31)]
    w: u32,
    #[clap(long, short, default_value_t = 31)]
    k: u32,
    #[clap(long, short, default_value_t = 1)]
    r: u32,
    #[clap(long, default_value_t = 0)]
    min_span: u32,
}

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();
    let filepath = args.in_fasta;
    let spec = ShmmrSpec {
        w: args.w,
        k: args.k,
        r: args.r,
        min_span: args.min_span,
        sketch: false,
    };
    let mut sdb = seq_db::CompactSeqDB::new(spec.clone());
    sdb.load_seqs_from_fastx(filepath)?;
    let mut out_buf = BufWriter::new(File::create(args.output_path)?);
    sdb.frag_map
        .into_iter()
        .try_for_each(|(k, v)| -> Result<(), std::io::Error> {
            let c = v.len(); 
            if c >= args.min_count {
                out_buf.write_fmt(format_args!("{:016x} {:016x} {}\n", k.0, k.1, c))?;
            };
            Ok(())
        })?;
    Ok(())
}
