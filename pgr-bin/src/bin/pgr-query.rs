const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use pgr_bin::{get_fastx_reader, GZFastaReader, SeqIndexDB};
use pgr_db::fasta_io::SeqRec;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;


/// Query a PGR-TK pangenome sequence database,
/// output the hit summary and generate fasta files from the target sequences
#[derive(Parser, Debug)]
#[clap(name = "pgr-query")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
    /// the prefix to a PGR-TK sequence database
    pgr_db_prefix: String,
    /// the path to the query fasta file
    query_fastx_path: String,
    /// the prefix of the output file
    output_prefix: String,

    /// using the frg format for the sequence database (default to the AGC backend database if not specified)
    #[clap(long, default_value_t = false)]
    frg_file: bool,

    #[clap(long, default_value_t = false)]
    fastx_file: bool,

    #[clap(long, short, default_value_t = 80)]
    w: u32,
    /// minimizer k-mer size
    #[clap(long, short, default_value_t = 56)]
    k: u32,
    /// sparse minimizer (shimmer) reduction factor
    #[clap(long, short, default_value_t = 4)]
    r: u32,
    /// min span for neighboring minimizers
    #[clap(long, short, default_value_t = 64)]
    min_span: u32,

    /// the gap penalty factor for sparse alignments in the SHIMMER space
    #[clap(long, short, default_value_t = 0.025)]
    gap_penalty_factor: f32,

    /// merge hits with the specified distance
    #[clap(long, short, default_value_t = 100000)]
    merge_range_tol: usize,

    /// the max count of SHIMMER used for the sparse alignment
    #[clap(long, default_value_t = 128)]
    max_count: u32,

    /// the max count of SHIMMER in the query sequences used for the sparse alignment
    #[clap(long, default_value_t = 128)]
    max_query_count: u32,

    /// the max count of SHIMMER in the targets sequences used for the sparse alignment
    #[clap(long, default_value_t = 128)]
    max_target_count: u32,

    /// the span of the chain for building the sparse alignment directed acyclic graph
    #[clap(long, default_value_t = 8)]
    max_aln_chain_span: u32,

    /// option only to output summaries
    #[clap(long, default_value_t = false)]
    only_summary: bool,

    /// output summaries in the bed format
    #[clap(long, default_value_t = false)]
    bed_summary: bool,
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
        #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
        GZFastaReader::GZFile(reader) => add_seqs(&mut reader.into_iter()),

        #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
        GZFastaReader::RegularFile(reader) => add_seqs(&mut reader.into_iter()),
    };

    let mut seq_index_db = SeqIndexDB::new();
    if args.frg_file {
        let stderr = io::stderr();
        let mut handle = stderr.lock();
        let _ = handle.write_all(b"the option `--frg_file` is specified, read the input file as a FRG backed index database files.");
        let _ = seq_index_db.load_from_frg_index(args.pgr_db_prefix);
    } else if args.fastx_file {
        let stderr = io::stderr();
        let mut handle = stderr.lock();
        let _ = handle.write_all(b"the option `--fastx_file` is specified, read the input file as a file file.");
        let _ =
            seq_index_db.load_from_fastx(args.pgr_db_prefix, args.w, args.k, args.r, args.min_span);
    } else {
        let stderr = io::stderr();
        let mut handle = stderr.lock();
        let _ = handle.write_all(b"Read the input as a AGC backed index database files.");
        let _ = seq_index_db.load_from_agc_index(args.pgr_db_prefix);
    }
    let prefix = Path::new(&args.output_prefix);

    let mut hit_file = if args.bed_summary {
        BufWriter::new(File::create(prefix.with_extension("hit.bed")).unwrap())
    } else {
        BufWriter::new(File::create(prefix.with_extension("hit")).unwrap())
    };
    if args.bed_summary {
        writeln!(
            hit_file,
            "#{}",
            [
                "target",
                "bgn",
                "end",
                "query",
                "color",
                "orientation",
                "q_len",
                "aln_anchor_count",
                "q_idx",
                "src",
                "ctg_bgn",
                "ctg_end",
            ]
            .join("\t")
        )?;
    } else {
        writeln!(
            hit_file,
            "#{}",
            [
                "out_seq_name",
                "ctg_bgn",
                "ctg_end",
                "color",
                "q_name",
                "orientation",
                "idx",
                "q_idx",
                "query_bgn",
                "query_end",
                "q_len",
                "aln_anchor_count",
            ]
            .join("\t")
        )?;
    };
    query_seqs
        .into_iter()
        .enumerate()
        .for_each(|(idx, seq_rec)| {
            let q_name = String::from_utf8_lossy(&seq_rec.id);
            let query_seq = seq_rec.seq;
            let q_len = query_seq.len();

            let query_results = seq_index_db.query_fragment_to_hps(
                query_seq,
                args.gap_penalty_factor,
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
                            let e = sid_to_alns.entry(sid).or_insert_with(Vec::new);
                            e.push((aln, orientation))
                        }
                    })
                });

                let mut aln_range = FxHashMap::default();
                sid_to_alns.into_iter().for_each(|(sid, alns)| {
                    alns.into_iter().for_each(|(aln, orientation)| {
                        let mut target_coordinates = aln
                            .iter()
                            .map(|v| (v.1 .0, v.1 .1))
                            .collect::<Vec<(u32, u32)>>();
                        target_coordinates.sort();
                        let bgn = target_coordinates[0].0;
                        let end = target_coordinates[target_coordinates.len() - 1].1;
                        let e = aln_range.entry(sid).or_insert_with(Vec::new);
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
                            .cloned()
                            .collect::<Vec<_>>();

                        let mut r_rgns = rgns
                            .iter()
                            .filter(|&v| v.3 == 1)
                            .cloned()
                            .collect::<Vec<_>>();

                        f_rgns.sort();
                        r_rgns.sort();

                        let mut out_rgns = vec![];
                        let mut last_rgn: (u32, u32, u32, u32, Vec<_>) = (0, 0, 0, 0, vec![]);
                        f_rgns.into_iter().for_each(|r| {
                            if last_rgn.4.is_empty() {
                                last_rgn = r;
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
                                    last_rgn = r;
                                }
                            }
                        });
                        if last_rgn.2 > 0 {
                            //not empty
                            out_rgns.push(last_rgn);
                        };

                        let mut last_rgn: (u32, u32, u32, u32, Vec<_>) = (0, 0, 0, 0, vec![]);
                        r_rgns.into_iter().for_each(|r| {
                            if last_rgn.4.is_empty() {
                                last_rgn = r;
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
                                    last_rgn = r;
                                }
                            }
                        });
                        if last_rgn.2 > 0 {
                            //not empty
                            out_rgns.push(last_rgn);
                        };

                        (sid, out_rgns)
                    })
                    .collect::<FxHashMap<_, _>>();

                let mut fasta_out = None;
                let mut fasta_buf: BufWriter<File>;
                if !args.only_summary {
                    let ext = format!("{:03}.fa", idx);
                    fasta_buf = BufWriter::new(File::create(prefix.with_extension(ext)).unwrap());
                    fasta_out = Some(&mut fasta_buf);
                };
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
                            let target_seq_name =
                            if args.fastx_file && args.bed_summary {
                                format!("{}", ctg)
                            } else {
                                format!("{}::{}_{}_{}_{}", base, ctg, b, e, orientation)
                            };

                            if args.bed_summary {
                                let _ = writeln!(
                                    hit_file,
                                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                    target_seq_name,
                                    b,
                                    e,
                                    q_name,
                                    "#AAAAAA",
                                    orientation,
                                    q_len,
                                    aln.len(),
                                    idx,
                                    src,
                                    q_bgn,
                                    q_end,
                                );
                            } else {
                                let _ = writeln!(
                                    hit_file,
                                    "{:03}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                    idx,
                                    q_name,
                                    q_bgn,
                                    q_end,
                                    q_len,
                                    aln.len(),
                                    src,
                                    ctg,
                                    b,
                                    e,
                                    orientation,
                                    target_seq_name
                                );
                            }
                            //println!("DBG: {}", seq_id);
                            if let Some(fasta_out) = fasta_out.as_mut() {
                                let target_seq = seq_index_db
                                    .get_sub_seq_by_id(sid, b as usize, e as usize)
                                    .unwrap();
                                let target_seq = if orientation == 1 {
                                    pgr_db::fasta_io::reverse_complement(&target_seq)
                                } else {
                                    target_seq
                                };
                                let _ = writeln!(fasta_out, ">{}", target_seq_name);
                                let _ =
                                    writeln!(fasta_out, "{}", String::from_utf8_lossy(&target_seq));
                            };
                        });
                });
            };
        });
    Ok(())
}
