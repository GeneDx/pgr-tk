const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use pgr_bin::SeqIndexDB;
use pgr_db::fasta_io;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// List or fetch sequences from a PGR-TK database
#[derive(Parser, Debug)]
#[clap(name = "pgr-fetch-seqs")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
    /// the prefix to a PGR-TK sequence database
    pgr_db_prefix: String,

    /// using the frg format for the sequence database (default to the AGC backend database if not specified)
    #[clap(long, default_value_t = false)]
    frg_file: bool,

    /// the regions file path
    #[clap(short, long, default_value=None)]
    region_file: Option<String>,

    /// output file name
    #[clap(short, long, default_value=None)]
    output_file: Option<String>,

    /// list all sequence source, contig names in the database
    #[clap(long, default_value_t = false)]
    list: bool,
}

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();

    let mut seq_index_db = SeqIndexDB::new();
    if args.frg_file {
        let _ = seq_index_db.load_from_frg_index(args.pgr_db_prefix);
    } else {
        let _ = seq_index_db.load_from_agc_index(args.pgr_db_prefix);
    }

    if args.list {
        let mut out = if args.output_file.is_some() {
            let f = File::open(args.output_file.unwrap()).expect("can't open the ouptfile");
            Box::new(f) as Box<dyn Write>
        } else {
            Box::new(io::stdout())
        };
        seq_index_db
            .seq_info
            .unwrap()
            .into_iter()
            .for_each(|(sid, (ctg, src, length))| {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}",
                    sid,
                    src.unwrap_or_else(|| "None".to_string()),
                    ctg,
                    length
                )
                .expect("can't write output file")
            });
        return Ok(());
    }

    let region_file = args.region_file.expect("region file not specified");
    let region_file =
        BufReader::new(File::open(Path::new(&region_file)).expect("can't open the region file"));

    let mut out = if args.output_file.is_some() {
        let f = BufWriter::new(
            File::create(args.output_file.clone().unwrap()).expect("can't open the ouptfile"),
        );
        Box::new(f) as Box<dyn Write>
    } else {
        Box::new(io::stdout())
    };
    
    region_file.lines().into_iter().for_each(|line| {
        let line = line.expect("fail to get a line in the region file");
        let fields = line.split('\t').collect::<Vec<&str>>();
        let label = fields[0].to_string();
        let src = fields[1].to_string();
        let ctg = fields[2].to_string();
        let bgn: usize = fields[3].parse().expect("can't parse bgn");
        let end: usize = fields[4].parse().expect("can't parse end");
        let reversed: bool = fields[4].parse::<u32>().expect("can't parse strand") == 1;
        let mut seq = seq_index_db
            .get_sub_seq(src, ctg, bgn, end)
            .expect("fail to fetch sequence");
        if reversed {
            seq = fasta_io::reverse_complement(&seq);
        }

        writeln!(out, ">{}", label).expect("fail to write the sequences");
        writeln!(out, "{}", String::from_utf8_lossy(&seq[..]))
            .expect("fail to write the sequences");
    });

    Ok(())
}
