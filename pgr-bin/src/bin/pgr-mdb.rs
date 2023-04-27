const VERSION_STRING: &str = env!("VERSION_STRING");

//use std::path::PathBuf;
use clap::{self, CommandFactory, Parser};

#[cfg(feature = "with_agc")]
use pgr_db::agc_io::AGCFile;

#[cfg(feature = "with_agc")]
use pgr_db::shmmrutils::ShmmrSpec;

#[cfg(feature = "with_agc")]
use std::fs::File;

#[cfg(feature = "with_agc")]
use std::io::{BufRead, BufReader};

#[cfg(feature = "with_agc")]
use pgr_db::seq_db;

/// Create pgr minimizer database with AGC backend
#[derive(Parser, Debug)]
#[clap(name = "pgr-mdb")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
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
    /// using sketch k-mer than minimizer
    #[clap(short, long)]
    sketch: bool,
    /// set to use agc prefecting feature (more memory usage but faster, useful for agcfile with many small contigs)
    #[clap(short, long)]
    prefetching: bool,
    /// number of parallel agc reader threads (more memory usage)
    #[clap(long, short, default_value_t = 4)]
    number_of_readers: usize,
}

#[cfg(feature = "with_agc")]
fn load_write_index_from_agcfile(
    path: String,
    prefix: String,
    shmmr_spec: &ShmmrSpec,
    prefetching: bool,
    number_of_readers: usize,
) -> Result<(), std::io::Error> {
    let mut sdb = seq_db::CompactSeqDB::new(shmmr_spec.clone());
    let filelist = File::open(path)?;

    BufReader::new(filelist).lines().into_iter().try_for_each(
        |fp| -> Result<(), std::io::Error> {
            let fp = fp.unwrap();
            //println!("load file {}", fp);
            let mut agcfile: AGCFile = AGCFile::new(fp)?;
            agcfile.set_iter_thread(number_of_readers);
            agcfile.set_prefetching(prefetching);
            //println!("start to load index");
            let _ = sdb.load_index_from_agcfile(agcfile);
            Ok(())
        },
    )?;

    //seq_db::write_shmr_map_file(&sdb.frag_map, "test.db".to_string());
    sdb.write_shmmr_map_index(prefix)?;
    Ok(())
}

fn main() {
    CmdOptions::command().version(VERSION_STRING).get_matches();

    #[cfg(feature = "with_agc")]
    let args = CmdOptions::parse();
    // TODO: to log file
    //println!("read data from files in {:?}", args.filepath);
    //println!("output prefix {:?}", args.prefix);

    #[cfg(feature = "with_agc")]
    let shmmr_spec = pgr_db::shmmrutils::ShmmrSpec {
        w: args.w,
        k: args.k,
        r: args.r,
        min_span: args.min_span,
        sketch: args.sketch,
    };

    #[cfg(feature = "with_agc")]
    load_write_index_from_agcfile(
        args.filepath,
        args.prefix.clone(),
        &shmmr_spec,
        args.prefetching,
        args.number_of_readers,
    )
    .unwrap();

    #[cfg(not(feature = "with_agc"))]
    panic!("the command is not compiled with `with_agc` feature")

}
