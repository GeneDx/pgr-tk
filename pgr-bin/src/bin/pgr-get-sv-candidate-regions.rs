const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use pgr_db::aln;
use pgr_db::ext::{get_fastx_reader, GZFastaReader, SeqIndexDB};
use pgr_db::fasta_io::{reverse_complement, SeqRec};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

/// Align long contigs and identify potential SV regions with respect to the reference fasta file
#[derive(Parser, Debug)]
#[clap(name = "pgr-get-sv-candidate-regions")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
    /// path to the reference fasta file
    reference_fasta_path: String,
    /// the path to the query assembly contig file
    assembly_contig_path: String,
    /// the prefix of the output file
    output_prefix: String,
    /// number of threads used in parallel (more memory usage), default to "0" using all CPUs available or the number set by RAYON_NUM_THREADS
    #[clap(long, default_value_t = 0)]
    number_of_thread: usize,
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

    /// the span of the chain for building the sparse alignment directed acyclic graph
    #[clap(long, default_value_t = 8)]
    max_aln_chain_span: u32,
}

// ((q_smp_start, q_smp_end, q_smp_orientation), (t_smp_start, t_smp_end, t_smp_orientation))
type AlignSegment = ((u32, u32, u8), (u32, u32, u8));

type AlignSegments = Vec<AlignSegment>;

fn filter_aln(aln_segs: &AlignSegments) -> Vec<((u32, u32), (u32, u32))> {
    let mut aln_segs = aln_segs.clone();

    let mut last_ts = aln_segs[0].1 .0;
    let mut last_te = aln_segs[0].1 .1;

    let mut last_qs = aln_segs[0].0 .0;
    let mut last_qe = aln_segs[0].0 .1;

    let mut rtn = Vec::<((u32, u32), (u32, u32))>::new();
    //rtn.push(((last_ts, last_te), (last_qs, last_qe)));
    for ((qs, qe, qo), (ts, te, to)) in aln_segs {
        if te < ts {
            continue;
        };
        if qo != to {
            continue;
        };
        if ts >= last_te {
            last_ts = last_te;
            last_te = te;

            last_qs = last_qe;
            last_qe = qe;
            if last_ts == last_te {
                continue;
            }
        }
        rtn.push(((last_ts, last_te), (last_qs, last_qe)));
    }
    rtn
}

fn filter_aln_rev(aln_segs: &AlignSegments) -> Vec<((u32, u32), (u32, u32))> {
    let mut aln_segs = aln_segs.clone();
    aln_segs.reverse();

    let mut last_ts = aln_segs[0].1 .0;
    let mut last_te = aln_segs[0].1 .1;

    let mut last_qs = aln_segs[0].0 .0;
    let mut last_qe = aln_segs[0].0 .1;

    let mut rtn = Vec::<((u32, u32), (u32, u32))>::new();
    //rtn.push(((last_ts, last_te), (last_qs, last_qe)));
    for ((qs, qe, qo), (ts, te, to)) in aln_segs {
        if te < ts {
            continue;
        };
        if qo == to {
            continue;
        };
        if ts >= last_te {
            last_ts = last_te;
            last_te = te;

            last_qe = last_qs;
            last_qs = qs;
            if last_ts == last_te {
                continue;
            }
        }
        rtn.push(((last_ts, last_te), (last_qs, last_qe)));
    }
    rtn
}

type AlignmentResult = (Vec<(u32, u32, char, String, String)>, Vec<(u32, u32, char)>);
pub fn get_variant_segments(
    target_str: &str,
    query_str: &str,
    max_wf_length: Option<u32>,
    mismatch_penalty: i32,
    open_penalty: i32,
    extension_penalty: i32,
    max_diff_percent: f32,
) -> Option<AlignmentResult> {
    let set_len_diff = (query_str.len() as i64 - target_str.len() as i64).unsigned_abs() as u32;
    let max_wf_length = if let Some(max_wf_length) = max_wf_length {
        max_wf_length
    } else {
        std::cmp::max(2 * set_len_diff, 128_u32)
    };

    if max_wf_length > 128
        && (max_wf_length as f32 / std::cmp::min(target_str.len(), query_str.len()) as f32)
            > max_diff_percent
    {
        return None;
    };

    if let Some((aln_target_str, aln_query_str)) = aln::wfa_align_bases(
        target_str,
        query_str,
        max_wf_length,
        mismatch_penalty,
        open_penalty,
        extension_penalty,
    ) {
        let aln_pairs = aln::wfa_aln_pair_map(&aln_target_str, &aln_query_str);
        Some((
            aln::get_variants_from_aln_pair_map(&aln_pairs, target_str, query_str),
            aln_pairs,
        ))
    } else {
        None
    }
}

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.number_of_thread)
        .build_global()
        .unwrap();

    let mut ref_seq_index_db = SeqIndexDB::new();
    ref_seq_index_db.load_from_fastx(
        args.reference_fasta_path,
        args.w,
        args.k,
        args.r,
        args.min_span,
    )?;

    let mut query_seqs: Vec<SeqRec> = vec![];
    let mut add_seqs = |seq_iter: &mut dyn Iterator<Item = io::Result<SeqRec>>| {
        seq_iter.into_iter().for_each(|r| {
            if let Ok(r) = r {
                query_seqs.push(r);
            };
        });
    };

    match get_fastx_reader(args.assembly_contig_path)? {
        #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
        GZFastaReader::GZFile(reader) => add_seqs(&mut reader.into_iter()),

        #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
        GZFastaReader::RegularFile(reader) => add_seqs(&mut reader.into_iter()),
    };

    let kmer_size = args.k;

    query_seqs
        .into_par_iter()
        .enumerate()
        .map(|(idx, seq_rec)| {
            // let q_name = String::from_utf8_lossy(&seq_rec.id);
            let query_seq = seq_rec.seq.clone();
            //let q_len = query_seq.len();

            let query_results = ref_seq_index_db.query_fragment_to_hps(
                &query_seq,
                args.gap_penalty_factor,
                Some(1),
                Some(1),
                Some(1),
                Some(args.max_aln_chain_span),
            );
            (idx, seq_rec, query_results)
        })
        .for_each(|(idx, seq_rec, query_results)| {
            if let Some(qr) = query_results {
                let q_name = String::from_utf8_lossy(&seq_rec.id);
                let query_seq = seq_rec.seq;
                let q_len = query_seq.len();
                let mut sid_to_mapped_regions = FxHashMap::default();
                qr.into_iter().for_each(|(sid, mapped_segments)| {
                    let mut aln_lens = vec![];
                    let mut f_count = 0_usize;
                    let mut r_count = 0_usize;
                    mapped_segments.into_iter().for_each(|(_score, aln)| {
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
                            let e = sid_to_mapped_regions.entry(sid).or_insert_with(Vec::new);
                            e.push((aln, orientation))
                        }
                    })
                });

                sid_to_mapped_regions
                    .into_iter()
                    .for_each(|(sid, mapped_regions)| {
                        let ref_seq = ref_seq_index_db.get_seq_by_id(sid).unwrap();
                        let mapped_region_aln = mapped_regions
                            .into_iter()
                            .map(|(aln_segs, orientation)| {
                                let aln_segs = if orientation == 0 {
                                    filter_aln(&aln_segs)
                                } else {
                                    filter_aln_rev(&aln_segs)
                                };

                                aln_segs
                                    .into_iter()
                                    .map(|((ts, te), (qs, qe))| {
                                        let s0str = ref_seq[ts as usize..te as usize].to_vec();
                                        let s1str = if orientation == 0 {
                                            query_seq[qs as usize..qe as usize].to_vec()
                                        } else {
                                            reverse_complement(
                                                &query_seq[(qs - kmer_size) as usize
                                                    ..(qe - kmer_size) as usize],
                                            )
                                        };
                                        // println!("XX0: {}", String::from_utf8_lossy(&s0str[..16]));
                                        // println!("XX1: {}", String::from_utf8_lossy(&s1str[..16]));
                                        // println!("YY0: {}", String::from_utf8_lossy(&s0str[s0str.len() - 16..]));
                                        // println!("YY1: {}", String::from_utf8_lossy(&s1str[s1str.len() - 16..]));

                                        let diff = if s0str.len() > 16 && s1str.len() > 16 {
                                            if s0str[..16] == s1str[..16]
                                                && s0str[s0str.len() - 16..]
                                                    == s1str[s1str.len() - 16..]
                                                && (s0str.len() as isize - s1str.len() as isize)
                                                    .abs()
                                                    < 256
                                            {
                                                let s0str = String::from_utf8_lossy(&s0str[..]);
                                                let s1str = String::from_utf8_lossy(&s1str[..]);
                                                get_variant_segments(
                                                    &s0str,
                                                    &s1str,
                                                    Some(64),
                                                    4,
                                                    3,
                                                    1,
                                                    1.0,
                                                )
                                            } else {
                                                None
                                            }
                                        } else {
                                            None
                                        };

                                        // println!("{:?} {:?}",  ((ts, te), (qs, qe), orientation), diff);
                                        // println!();
                                        ((ts, te), (qs, qe), orientation, diff)
                                    })
                                    .collect::<Vec<_>>()
                            })
                            .collect::<Vec<_>>();

                        let ref_ctg_name = ref_seq_index_db
                            .seq_info
                            .as_ref()
                            .unwrap()
                            .get(&sid)
                            .unwrap()
                            .0
                            .clone();

                        mapped_region_aln.into_iter().for_each(|v| {
                            println!("B\t{}\t{}\t{}", ref_ctg_name, q_name, q_len); 
                            v.into_iter()
                                .for_each(|((ts, te), (qs, qe), orientation, diff)| {
                                    let qs = if orientation == 0 { qs } else { qs - kmer_size };
                                    let qe = if orientation == 0 { qe } else { qe - kmer_size };
                                    if let Some(diff) = diff {
                                        if diff.0.is_empty() {
                                            println!(
                                                "M\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                                ref_ctg_name, ts, te, q_name, qs, qe, orientation
                                            );
                                        } else {
                                            diff.0.into_iter().for_each(
                                                |(td, qd, vt, t_str, q_str)| {
                                                    println!(
                                                        "V\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                                        ref_ctg_name,
                                                        ts,
                                                        te,
                                                        q_name,
                                                        qs,
                                                        qe,
                                                        orientation,
                                                        ts + td,
                                                        qs + qd,
                                                        vt,
                                                        t_str,
                                                        q_str
                                                    );
                                                },
                                            )
                                        }
                                    } else {
                                        println!(
                                            "S\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                            ref_ctg_name, ts, te, q_name, qs, qe, orientation
                                        );
                                    }
                                });
                            println!("E\t{}\t{}\t{}", ref_ctg_name, q_name, q_len);
                        });
                    })
            };
        });
    Ok(())
}
