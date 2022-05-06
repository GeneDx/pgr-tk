const VERSION_STRING: &'static str = env!("VERSION_STRING");

//use std::path::PathBuf;
use clap::{self, IntoApp, Parser};

use pgr_db::agc_io::AGCFile;
use pgr_db::shmmrutils::ShmmrSpec;
use std::fs::File;
use std::io::{BufRead, BufReader};

use pgr_db::seq_db;

#[derive(Parser, Debug)]
#[clap(name = "pgr-mdb")]
#[clap(author, version)]
#[clap(about = "create pgr minimizer db", long_about = None)]
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
}

fn load_write_index_from_agcfile(
    path: String,
    prefix: String,
    shmmr_spec: &ShmmrSpec,
) -> Result<(), std::io::Error> {
    let mut sdb = seq_db::CompactSeqDB::new(shmmr_spec.clone());
    let filelist = File::open(path)?;

    BufReader::new(filelist).lines().into_iter().try_for_each(|fp| -> Result<(), std::io::Error> {
        let fp = fp.unwrap();
        let agcfile = AGCFile::new(fp)?;
        let _ = sdb.load_index_from_agcfile(agcfile);
        Ok(())
    });

    //seq_db::write_shmr_map_file(&sdb.frag_map, "test.db".to_string());
    sdb.write_shmr_map_index(prefix)?;
    Ok(())
}

fn _dump_index_mdb(prefix: String) -> Result<(), std::io::Error> {
    let (shmmr_spec, new_map) = seq_db::read_mdb_file(prefix + ".mdb")?;
    println!(
        "# {} {} {} {} {}",
        shmmr_spec.w, shmmr_spec.k, shmmr_spec.r, shmmr_spec.min_span, shmmr_spec.sketch
    );
    new_map.iter().for_each(|v| {
        for vv in v.1.iter() {
            println!(
                "{:016x} {:016x} {} {} {} {} {}",
                v.0 .0, v.0 .1, vv.0, vv.1, vv.2, vv.3, vv.4
            );
        }
    });
    Ok(())
}

fn main() {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();
    // TODO: to log file
    // println!("read data from files in {:?}", args.filepath);
    // println!("output prefix {:?}", args.prefix);
    let shmmr_spec = pgr_db::shmmrutils::ShmmrSpec {
        w: args.w,
        k: args.k,
        r: args.r,
        min_span: args.min_span,
        sketch: args.sketch,
    };
    load_write_index_from_agcfile(args.filepath, args.prefix.clone(), &shmmr_spec).unwrap();
    // dump_index_mdb(args.prefix).unwrap();

    //load_write_index_from_agcfile();
}
