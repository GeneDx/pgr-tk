const VERSION_STRING: &'static str = "abc";

use std::path::PathBuf;
use clap::{self,Parser, IntoApp};

use pgr_db::agc_io::AGCFile;
use std::fs::File;
use std::io::{BufRead, BufReader};

use pgr_db::seq_db;

#[derive(Parser, Debug)]
#[clap(name = "pgr-mdb")]
#[clap(author, version, about, long_about = None)]
struct CmdOptions {
    filepath: String,
    prefix: String,
}

fn load_write_index_from_agcfile(path: String, prefix: String) -> Result<(), std::io::Error> {
    let mut sdb = seq_db::CompactSeqDB::new();
    let filelist = File::open(path)?;

    BufReader::new(filelist).lines().into_iter().for_each(|fp| {
        let fp = fp.unwrap();
        let agcfile = AGCFile::new(fp);
        let _ = sdb.load_index_from_agcfile(agcfile);
    });

    //seq_db::write_shmr_map_file(&sdb.frag_map, "test.db".to_string());
    sdb.write_shmr_map_index(prefix)?;
    Ok(())
}

fn dump_index_mdb(prefix: String) -> Result<(), std::io::Error> {
    let (shmmr_spec, new_map) = seq_db::read_shmr_map_file(prefix+".mdb")?;
    println!("# {} {} {} {} {}", shmmr_spec.w, shmmr_spec.k, shmmr_spec.r, shmmr_spec.min_span, shmmr_spec.sketch);
    new_map.iter().for_each(|v| {
        for vv in v.1.iter() {
            println!("{:016x} {:016x} {} {} {} {} {}", v.0.0, v.0.1, vv.0, vv.1, vv.2, vv.3, vv.4);
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
    load_write_index_from_agcfile(args.filepath, args.prefix.clone()).unwrap();
    // dump_index_mdb(args.prefix).unwrap();
    

    
    //load_write_index_from_agcfile();
}