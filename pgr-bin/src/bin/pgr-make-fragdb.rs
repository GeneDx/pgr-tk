const VERSION_STRING: &str = env!("VERSION_STRING");

//use std::path::PathBuf;
use clap::{self, CommandFactory, Parser};

use pgr_bin::SeqIndexDB;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Parser, Debug)]
#[clap(name = "pgr-mdb")]
#[clap(author, version)]
#[clap(about = "create pgr fragment minimizer db", long_about = None)]
struct CmdOptions {
    // file contains the paths to the fastx files to load
    filepath: String,
    prefix: String,
    /// minimizer window size
    #[clap(long, short, default_value_t = 80)]
    w: u32,
    /// minimizer k-mer size
    #[clap(long, short, default_value_t = 56)]
    k: u32,
    /// sparse minimizer (shimmer) reduction factor
    #[clap(long, short, default_value_t = 4)]
    r: u32,
    /// min span for neighboring minimiers
    #[clap(long, short, default_value_t = 64)]
    min_span: u32,
}

fn main() {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();
    // TODO: to log file
    //println!("read data from files in {:?}", args.filepath);
    //println!("output prefix {:?}", args.prefix);
    let _shmmr_spec = pgr_db::shmmrutils::ShmmrSpec {
        w: args.w,
        k: args.k,
        r: args.r,
        min_span: args.min_span,
        sketch: false,
    };
    let mut sdb = SeqIndexDB::new();
    let input_files = BufReader::new(
        File::open(Path::new(&args.filepath))
            .expect("can't open the input file that contains the paths to the fastx files"),
    );
    input_files
        .lines()
        .into_iter()
        .enumerate()
        .for_each(|(fid, filename)| {
            let filepath = filename
                .expect("can't get fastx file name")
                .trim()
                .to_string();
            if fid == 0 {
                sdb.load_from_fastx(filepath.clone(), args.w, args.k, args.r, args.min_span)
                    .expect(format!("fail to read the fastx file: {}", filepath).as_str());
            } else {
                sdb.append_from_fastx(filepath.clone())
                    .expect(format!("fail to read the fastx file: {}", filepath).as_str());
            }
        });

    sdb.write_frag_and_index_files(args.prefix);
}
