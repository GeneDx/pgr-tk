const VERSION_STRING: &'static str = env!("VERSION_STRING");
use clap::{self, IntoApp, Parser};
use flate2::bufread::MultiGzDecoder;
use pgr_db::fasta_io::{reverse_complement, FastaReader, FastqStreamReader, SeqRec};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};

#[derive(Parser, Debug)]
#[clap(name = "pgr-probe-match")]
#[clap(author, version)]
#[clap(about = "using Cuckoo Filter for Matching Reads To A Reference Set of Sequences", long_about = None)]
struct CmdOptions {
    probe_file_path: String,
    #[clap(long, short)]
    query_fastx_path: Option<String>,
}
enum GZFastaReader {
    GZFile(FastaReader<BufReader<MultiGzDecoder<BufReader<File>>>>),
    RegularFile(FastaReader<BufReader<BufReader<File>>>),
}

#[derive(Clone)]
struct ProbeInfo {
    vname: String,
    vprobe: Vec<u8>,
    vprobe_r: Vec<u8>,
    t1name: String,
    t1probe: Vec<u8>,
    t1probe_r: Vec<u8>,
    t2name: String,
    t2probe: Vec<u8>,
    t2probe_r: Vec<u8>,
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
    let probe_reader = BufReader::new(File::open(args.probe_file_path)?);
    let mut all_probes = FxHashMap::<String, ProbeInfo>::default();
    probe_reader
        .lines()
        .into_iter()
        .for_each(|line| match line {
            Ok(line) => {
                let line = line.trim_end();
                let mut fields = line.split("\t");
                let vname = fields.next().expect("error parsing").to_string();
                let tmp = fields.next().expect("error parsing");
                let vprobe = tmp.as_bytes().to_vec();
                let vprobe_r = reverse_complement(&vprobe);
                let t1name = fields.next().expect("error parsing").to_string();
                let tmp = fields.next().expect("error parsing");
                let t1probe = tmp.as_bytes().to_vec();
                let t1probe_r = reverse_complement(&t1probe);
                let t2name = fields.next().expect("error parsing").to_string();
                let tmp = fields.next().expect("error parsing");
                let t2probe = tmp.as_bytes().to_vec();
                let t2probe_r = reverse_complement(&t2probe);
                let probeset = ProbeInfo {
                    vname: vname.clone(),
                    vprobe,
                    vprobe_r,
                    t1name,
                    t1probe,
                    t1probe_r,
                    t2name,
                    t2probe,
                    t2probe_r,
                };
                all_probes.insert(vname, probeset);
            }
            _ => {}
        });

    let match_probe = |seq: &Vec<u8>, probe: &Vec<u8>| -> bool {
        let plen = probe.len();
        let mut flag = false;
        for i in 0..seq.len() - plen {
            if seq[i..i + plen] == probe[..] || seq[i..i + plen] == probe[..] {
                flag = true;
                break;
            }
        }
        flag
    };

    let check_seqs = |seq_iter: &mut dyn Iterator<Item = io::Result<SeqRec>>| {
        let mut seq_data = Vec::<SeqRec>::new();
        for r in seq_iter {
            if let Ok(r) = r {
                seq_data.push(r);
            }
        }

        all_probes.into_iter().for_each(|(_vname, probe_info)| {
            let mut count = (0_usize, 0_usize, 0_usize);
            (&seq_data)
                .into_par_iter()
                .filter(|&r| {
                    match_probe(&r.seq, &probe_info.vprobe)
                        || match_probe(&r.seq, &probe_info.vprobe_r)
                })
                .collect::<Vec<&SeqRec>>()
                .into_iter()
                .for_each(|r| {
                    count.0 += 1;
                    if match_probe(&r.seq, &probe_info.t1probe)
                        || match_probe(&r.seq, &probe_info.t1probe_r)
                    {
                        count.1 += 1;
                    }
                    if match_probe(&r.seq, &probe_info.t2probe)
                        || match_probe(&r.seq, &probe_info.t2probe_r)
                    {
                        count.2 += 1;
                    }
                });
            println!(
                "{} {} {} {} {} {}",
                probe_info.vname, count.0, probe_info.t1name, count.1, probe_info.t2name, count.2
            );
        });
    };

    if args.query_fastx_path.is_some() {
        match get_fastx_reader(args.query_fastx_path.unwrap())? {
            GZFastaReader::GZFile(reader) => check_seqs(&mut reader.into_iter()),
            GZFastaReader::RegularFile(reader) => check_seqs(&mut reader.into_iter()),
        }
    } else {
        let reader = FastqStreamReader::new(128);
        check_seqs(&mut reader.into_iter());
    }

    Ok(())
}
