const VERSION_STRING: &'static str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use pgr_bin::{get_fastx_reader, GZFastaReader, SeqIndexDB};
use pgr_db::fasta_io::SeqRec;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

#[derive(Parser, Debug)]
#[clap(name = "pgr-query")]
#[clap(author, version)]
#[clap(about = "query a pgr pangenome database and ouput the hits", long_about = None)]
struct CmdOptions {
    pgr_db_prefix: String,
    query_fastx_path: String,
    output_prfix: String,

    #[clap(long, default_value_t=false)]
    frg_file: bool,

    #[clap(long, short, default_value_t = 0.025)]
    gap_penality_factor: f32,

    #[clap(long, short, default_value_t = 100000)]
    merge_range_tol: usize,

    #[clap(long, default_value_t = 128)]
    max_count: u32,

    #[clap(long, default_value_t = 128)]
    max_query_count: u32,

    #[clap(long, default_value_t = 128)]
    max_target_count: u32,

    #[clap(long, default_value_t = 8)]
    max_aln_chain_span: u32,
}

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();

    let mut query_seqs: Vec<SeqRec> = vec![];
    let mut add_seqs = |seq_iter: &mut dyn Iterator<Item = io::Result<SeqRec>>| {
        seq_iter.into_iter().for_each(|r| {
            if let Ok(r) = r {
                query_seqs.push(r);
            };
        });
    };

    match get_fastx_reader(args.query_fastx_path)? {
        GZFastaReader::GZFile(reader) => add_seqs(&mut reader.into_iter()),

        GZFastaReader::RegularFile(reader) => add_seqs(&mut reader.into_iter()),
    };

    let mut seq_index_db = SeqIndexDB::new();
    if args.frg_file {
        let _ = seq_index_db.load_from_frg_index(args.pgr_db_prefix);
    } else {
        let _ = seq_index_db.load_from_agc_index(args.pgr_db_prefix);
    }
    let prefix = Path::new(&args.output_prfix);
    let mut hit_file = BufWriter::new(File::create(prefix.with_extension("hit")).unwrap());
    writeln!(
        hit_file,
        "#{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        "q_idx",
        "q_name",
        "query_bgn",
        "query_end",
        "aln_anchor_count",
        "src",
        "ctg",
        "ctg_bgn",
        "ctg_end",
        "orientation",
        "out_seq_name"
    )?;

    query_seqs
        .into_iter()
        .enumerate()
        .for_each(|(idx, seq_rec)| {
            let q_name = String::from_utf8_lossy(&seq_rec.id);
            let query_seq = seq_rec.seq;

            let query_results = seq_index_db.query_fragment_to_hps(
                query_seq,
                args.gap_penality_factor,
                Some(args.max_count),
                Some(args.max_query_count),
                Some(args.max_target_count),
                Some(args.max_aln_chain_span),
            );

            if let Some(qr) = query_results {
                let mut sid_to_alns = FxHashMap::default();
                qr.into_iter().for_each(|(sid, alns)| {
                    let mut aln_lens = vec![];
                    let mut f_count = 0_usize;
                    let mut r_count = 0_usize;
                    alns.into_iter().for_each(|(_score, aln)| {
                        if aln.len() > 2 {
                            aln_lens.push(aln.len());
                            for hp in &aln {
                                if hp.0 .2 == hp.1 .2 {
                                    f_count += 1;
                                } else {
                                    r_count += 1;
                                }
                            }
                            let orientation = if f_count > r_count { 0_u32 } else { 1_u32 };
                            let e = sid_to_alns.entry(sid).or_insert(vec![]);
                            e.push((aln, orientation))
                        }
                    })
                });

                let mut aln_range = FxHashMap::default();
                sid_to_alns.into_iter().for_each(|(sid, alns)| {
                    alns.into_iter().for_each(|(aln, orientation)| {
                        let mut target_coordiantes = aln
                            .iter()
                            .map(|v| (v.1 .0, v.1 .1))
                            .collect::<Vec<(u32, u32)>>();
                        target_coordiantes.sort();
                        let bgn = target_coordiantes[0].0;
                        let end = target_coordiantes[target_coordiantes.len() - 1].1;
                        let e = aln_range.entry(sid).or_insert(vec![]);
                        e.push((bgn, end, end - bgn, orientation, aln));
                    })
                });

                // TODO: merge aln_range
                let aln_range = aln_range
                    .into_iter()
                    .map(|(sid, rgns)| {
                        let mut f_rgns = rgns
                            .iter()
                            .filter(|&v| v.3 == 0)
                            .map(|v| v.clone())
                            .collect::<Vec<_>>();

                        let mut r_rgns = rgns
                            .iter()
                            .filter(|&v| v.3 == 1)
                            .map(|v| v.clone())
                            .collect::<Vec<_>>();

                        f_rgns.sort();
                        r_rgns.sort();

                        let mut out_rgns = vec![];
                        let mut last_rgn: (u32, u32, u32, u32, Vec<_>) = (0, 0, 0, 0, vec![]);
                        f_rgns.into_iter().for_each(|r| {
                            if last_rgn.4.len() == 0 {
                                last_rgn = r.clone();
                                return;
                            } else {
                                let l_bgn = last_rgn.0;
                                let l_end = last_rgn.1;
                                assert!(l_end > l_bgn);
                                let r_bgn = r.0;
                                let r_end = r.1;
                                if (r_bgn as i64) - (l_end as i64) < args.merge_range_tol as i64 {
                                    let bgn = l_bgn;
                                    let end = if r_end > l_end { r_end } else { l_end };
                                    let len = end - bgn;
                                    let orientation = last_rgn.3;
                                    let mut aln = last_rgn.4.clone();
                                    aln.extend(r.4);
                                    last_rgn = (bgn, end, len, orientation, aln);
                                } else {
                                    out_rgns.push(last_rgn.clone());
                                    last_rgn = r.clone();
                                }
                            }
                        });
                        if last_rgn.2 > 0 {
                            //not empty
                            out_rgns.push(last_rgn.clone());
                        };

                        let mut last_rgn: (u32, u32, u32, u32, Vec<_>) = (0, 0, 0, 0, vec![]);
                        r_rgns.into_iter().for_each(|r| {
                            if last_rgn.4.len() == 0 {
                                last_rgn = r.clone();
                                return;
                            } else {
                                let l_bgn = last_rgn.0;
                                let l_end = last_rgn.1;
                                assert!(l_end > l_bgn);
                                let r_bgn = r.0;
                                let r_end = r.1;
                                if (r_bgn as i64) - (l_end as i64) < args.merge_range_tol as i64 {
                                    let bgn = l_bgn;
                                    let end = if r_end > l_end { r_end } else { l_end };
                                    let len = end - bgn;
                                    let orientation = last_rgn.3;
                                    let mut aln = last_rgn.4.clone();
                                    aln.extend(r.4);
                                    last_rgn = (bgn, end, len, orientation, aln);
                                } else {
                                    out_rgns.push(last_rgn.clone());
                                    last_rgn = r.clone();
                                }
                            }
                        });
                        if last_rgn.2 > 0 {
                            //not empty
                            out_rgns.push(last_rgn.clone());
                        };

                        (sid, out_rgns)
                    })
                    .collect::<FxHashMap<_, _>>();

                let ext = format!("{:03}.fa", idx);
                let mut fasta_out =
                    BufWriter::new(File::create(prefix.with_extension(ext)).unwrap());

                aln_range.into_iter().for_each(|(sid, rgns)| {
                    let (ctg, src, _ctg_len) =
                        seq_index_db.seq_info.as_ref().unwrap().get(&sid).unwrap();
                    //let src = *src.unwrap_or("N/A".to_string()).to_string();
                    let src = (*src).as_ref().unwrap_or(&"N/A".to_string()).clone();
                    rgns.into_iter()
                        .for_each(|(b, e, _, orientation, mut aln)| {
                            aln.sort();
                            let q_bgn = aln[0].0 .0;
                            let q_end = aln[aln.len() - 1].0 .1;
                            let base = Path::new(&src).file_stem().unwrap().to_string_lossy();
                            
                            let target_seq_name = format!("{}::{}_{}_{}_{}", base, ctg, b, e, orientation);
                            let _ = writeln!(
                                hit_file,
                                "{:03}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                idx,
                                q_name,
                                q_bgn,
                                q_end,
                                aln.len(),
                                src,
                                ctg,
                                b,
                                e,
                                orientation,
                                target_seq_name
                            );
                            //println!("DBG: {}", seq_id);
                            let target_seq = seq_index_db
                                .get_sub_seq_by_id(sid, b as usize, e as usize)
                                .unwrap();
                            let target_seq = if orientation == 1 {
                                pgr_db::fasta_io::reverse_complement(&target_seq)                               
                            } else {
                                target_seq
                            };
                            let _ = writeln!(fasta_out, ">{}", target_seq_name);
                            let _ = writeln!(fasta_out, "{}", String::from_utf8_lossy(&target_seq));
                        });
                });
            };
        });
    Ok(())
}
