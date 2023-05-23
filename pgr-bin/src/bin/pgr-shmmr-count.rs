const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};

use flate2::bufread::MultiGzDecoder;
use pgr_db::fasta_io::{FastaReader, SeqRec};
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
    /// the target fasta file path
    shmmr_target_fastx: String,

    /// ref_fasta
    ref_fastx: String,

    /// read_fasta
    read_fastx: String,

    /// output file name
    #[clap(short, long, default_value=None)]
    output_file: Option<String>,

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

    let mut shmmr_count = FxHashMap::<u64, (usize, usize)>::default();

    let mut set_shmmr_count = |seq_iter: &mut dyn Iterator<Item = io::Result<SeqRec>>| {
        let mut target_seqs: Vec<SeqRec> = vec![];
        seq_iter.into_iter().for_each(|r| {
            if let Ok(r) = r {
                target_seqs.push(r);
            };
        });
        target_seqs
            .into_par_iter()
            .map(|seq_rec| {
                pgr_db::shmmrutils::sequence_to_shmmrs1(
                    0,
                    &seq_rec.seq,
                    args.w,
                    args.k,
                    args.r,
                    args.min_span,
                    false,
                )
                .iter()
                .map(|&mmer| mmer.hash())
                .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
            .into_iter()
            .for_each(|hashes| {
                hashes.into_iter().for_each(|hash| {
                    shmmr_count.insert(hash, (0, 0));
                })
            });
    };

    match get_fastx_reader(args.shmmr_target_fastx)? {
        #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
        GZFastaReader::GZFile(reader) => set_shmmr_count(&mut reader.into_iter()),

        #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
        GZFastaReader::RegularFile(reader) => set_shmmr_count(&mut reader.into_iter()),
    };

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
                    let mut e = shmmr_count.entry(k).or_default();
                    (*e).0 += v;
                });
                ref_shmmr_location.extend(locations);
            });
    };

    match get_fastx_reader(args.ref_fastx)? {
        #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
        GZFastaReader::GZFile(reader) => count_ref_seq_shmmrs(&mut reader.into_iter()),

        #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
        GZFastaReader::RegularFile(reader) => count_ref_seq_shmmrs(&mut reader.into_iter()),
    };

    let mut count_seq_group_shmmrs = |seqs: &Vec<SeqRec>| {
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
                    let mut e = shmmr_count.entry(k).or_default();
                    (*e).1 += v;
                });
            });
    };

    let mut count_read_seq_shmmrs = |seq_iter: &mut dyn Iterator<Item = io::Result<SeqRec>>| {
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
        count_seq_group_shmmrs(&read_seqs);
        read_seqs.clear();
    };

    match get_fastx_reader(args.read_fastx)? {
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
        let (c0, c1) = *shmmr_count.get(&hash).unwrap();
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}",
            ctg,
            pos,
            pos + args.k as usize,
            c1 as f32 / c0 as f32,
            c1,
            c0
        )
        .expect("writing output error");
    });

    Ok(())
}
