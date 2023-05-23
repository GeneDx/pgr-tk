const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use pgr_db::ext::SeqIndexDB;
use pgr_db::fasta_io;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::hash::Hash;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use flate2::bufread::MultiGzDecoder;
use pgr_db::fasta_io::{FastaReader, FastaStreamReader, FastqStreamReader, SeqRec};
enum GZFastaReader {
    GZFile(FastaReader<BufReader<MultiGzDecoder<BufReader<File>>>>),
    RegularFile(FastaReader<BufReader<BufReader<File>>>),
}

/// List or fetch sequences from a PGR-TK database
#[derive(Parser, Debug)]
#[clap(name = "pgr-shmmr-count")]
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

    /// ref_fasta
    #[clap(short, long, default_value=None)]
    ref_fastx: Option<String>,

    /// read_fasta
    #[clap(short, long, default_value=None)]
    read_fastx: Option<String>,

    /// output file name
    #[clap(short, long, default_value=None)]
    output_file: Option<String>,

    /// list all sequence source, contig names in the database
    #[clap(long, default_value_t = false)]
    list: bool,

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

    let mut seq_index_db = SeqIndexDB::new();

    #[cfg(feature = "with_agc")]
    if args.frg_file {
        let _ = seq_index_db.load_from_frg_index(args.pgr_db_prefix);
    } else {
        let _ = seq_index_db.load_from_agc_index(args.pgr_db_prefix);
    }
    #[cfg(not(feature = "with_agc"))]
    if args.frg_file {
        let _ = seq_index_db.load_from_frg_index(args.pgr_db_prefix);
    } else {
        panic!("This command is compiled with only frg file support, please specify `--frg-file");
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

    let ref_fastx = args.ref_fastx.expect("ref_fastx file not specified");

    let read_fastx = args.read_fastx.expect("read_fastx file not specified");

    let mut shmmr_count = FxHashMap::<u64, (usize, usize)>::default();

    region_file.lines().into_iter().for_each(|line| {
        let line = line.expect("fail to get a line in the region file");
        let fields = line.split('\t').collect::<Vec<&str>>();
        let _label = fields[0].to_string();
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
        };

        let shmmrs = pgr_db::shmmrutils::sequence_to_shmmrs1(
            0,
            &seq,
            args.w,
            args.k,
            args.r,
            args.min_span,
            false,
        );
        shmmrs.into_iter().for_each(|mmer| {
            shmmr_count.insert(mmer.hash(), (0, 0));
        });
    });

    let mut ref_seqs: Vec<SeqRec> = vec![];
    let mut ref_shmmr_location = Vec::<(usize, usize, u64)>::new();
    let mut sid_to_ctg = FxHashMap::<usize, Vec<u8>>::default();
    let count_ref_seq_shmmrs = |seq_iter: &mut dyn Iterator<Item = io::Result<SeqRec>>| {
        seq_iter.into_iter().for_each(|r| {
            if let Ok(r) = r {
                ref_seqs.push(r);
            };
        });

        let ctg_to_sid = ref_seqs
            .iter()
            .enumerate()
            .map(|(sid, seq_rec)| (seq_rec.id.clone(), sid))
            .collect::<FxHashMap<Vec<u8>, usize>>();

        ctg_to_sid.iter().for_each(|(id, sid)| {
            sid_to_ctg.insert(*sid, id.clone());
        });

        ref_seqs
            .into_par_iter()
            .map(|seq_rec| {
                let mut partial_shmmr_count = FxHashMap::<u64, usize>::default();
                let mut partial_ref_shmmr_location = Vec::<(usize, usize, u64)>::new();
                let sid = ctg_to_sid.get(&seq_rec.id).unwrap();
                let shmmrs = pgr_db::shmmrutils::sequence_to_shmmrs1(
                    *sid as u32,
                    &seq_rec.seq,
                    args.w,
                    args.k,
                    args.r,
                    args.min_span,
                    false,
                );
                shmmrs.into_iter().for_each(|mmer| {
                    let hash = mmer.hash();
                    if shmmr_count.contains_key(&hash) {
                        let e = partial_shmmr_count.entry(mmer.hash()).or_insert_with(|| 0);
                        *e += 1;
                        partial_ref_shmmr_location.push((*sid, mmer.pos() as usize, hash));
                    }
                });
                (
                    partial_shmmr_count.into_iter(),
                    partial_ref_shmmr_location.into_iter(),
                )
            })
            .collect::<Vec<_>>()
            .into_iter()
            .for_each(|(counts, locations)| {
                counts.for_each(|(k, v)| {
                    let mut e = *shmmr_count.get(&k).unwrap();
                    e.0 += v;
                });
                ref_shmmr_location.extend(locations);
            });
    };

    match get_fastx_reader(ref_fastx)? {
        #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
        GZFastaReader::GZFile(reader) => count_ref_seq_shmmrs(&mut reader.into_iter()),

        #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
        GZFastaReader::RegularFile(reader) => count_ref_seq_shmmrs(&mut reader.into_iter()),
    };

    let count_seq_group_shmmrs = |seqs: &Vec<SeqRec>| {
        seqs.into_par_iter()
            .map(|seq_rec| {
                let mut partial_shmmr_count = FxHashMap::<u64, usize>::default();
                let shmmrs = pgr_db::shmmrutils::sequence_to_shmmrs1(
                    0,
                    &seq_rec.seq,
                    args.w,
                    args.k,
                    args.r,
                    args.min_span,
                    false,
                );
                shmmrs.into_iter().for_each(|mmer| {
                    let hash = mmer.hash();
                    if shmmr_count.contains_key(&hash) {
                        let e = partial_shmmr_count.entry(mmer.hash()).or_insert_with(|| 0);
                        *e += 1;
                    }
                });
                partial_shmmr_count.into_iter()
            })
            .collect::<Vec<_>>()
            .into_iter()
            .for_each(|counts| {
                counts.for_each(|(k, v)| {
                    let mut e = *shmmr_count.get(&k).unwrap();
                    e.1 += v;
                });
            });
    };

    let count_read_seq_shmmrs = |seq_iter: &mut dyn Iterator<Item = io::Result<SeqRec>>| {
        let mut read_seqs: Vec<SeqRec> = vec![];
        seq_iter.into_iter().for_each(|r| {
            if let Ok(r) = r {
                read_seqs.push(r);
                if read_seqs.len() >= 128 {
                    count_seq_group_shmmrs(&read_seqs);
                    read_seqs.clear();
                }
            };
        });
    };

    match get_fastx_reader(read_fastx)? {
        #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
        GZFastaReader::GZFile(reader) => count_read_seq_shmmrs(&mut reader.into_iter()),

        #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
        GZFastaReader::RegularFile(reader) => count_read_seq_shmmrs(&mut reader.into_iter()),
    };

    ref_shmmr_location.par_sort();

    let mut out = if args.output_file.is_some() {
        let f = BufWriter::new(
            File::create(args.output_file.clone().unwrap()).expect("can't open the output file"),
        );
        Box::new(f) as Box<dyn Write>
    } else {
        Box::new(io::stdout())
    };

    ref_shmmr_location.into_iter().for_each(|(sid, pos, hash)| {
        let ctg = String::from_utf8_lossy(sid_to_ctg.get(&sid).unwrap());
        let (c0, c1) = shmmr_count.get(&hash).unwrap();
        writeln!(out, "{} {} {} {}", ctg, pos, c0, c1).expect("writing output error");
    });

    Ok(())
}
