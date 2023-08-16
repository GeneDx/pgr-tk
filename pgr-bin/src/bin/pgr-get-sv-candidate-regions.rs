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
    /// the prefix of the output files
    output_prefix: String,
    /// number of threads used in parallel (more memory usage), default to "0" using all CPUs available or the number set by RAYON_NUM_THREADS
    #[clap(long, default_value_t = 0)]
    number_of_thread: usize,
    #[clap(long, short, default_value_t = 80)]
    w: u32,
    /// minimizer k-mer size
    #[clap(long, short, default_value_t = 55)]
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

#[derive(Clone)]
enum Record {
    Bgn(u32, u32, u32, u32, u32, u32, u32, u32),
    End(u32, u32, u32, u32, u32, u32, u32, u32),
    Match(u32, u32, u32, u32, u32, u32, u32),
    SvCnd(u32, u32, u32, u32, u32, u32, u32),
    Variant(
        u32,
        u32,
        u32,
        u32,
        u32,
        u32,
        u32,
        u32,
        u32,
        u32,
        char,
        String,
        String,
    ),
}

// ((q_smp_start, q_smp_end, q_smp_orientation), (t_smp_start, t_smp_end, t_smp_orientation))
type AlignSegment = ((u32, u32, u8), (u32, u32, u8));

type AlignSegments = Vec<AlignSegment>;

fn filter_aln(aln_segs: &AlignSegments) -> Vec<((u32, u32), (u32, u32))> {
    // the aln_segs should be sorted already
    let aln_segs = aln_segs.clone();

    let mut last_ts = aln_segs[0].1 .0;
    let mut last_te = aln_segs[0].1 .1;

    let mut last_qs = aln_segs[0].0 .0;
    let mut last_qe = aln_segs[0].0 .1;

    let mut rtn = Vec::<((u32, u32), (u32, u32))>::new();
    rtn.push(((last_ts, last_te), (last_qs, last_qe)));
    for ((_qs, qe, qo), (ts, te, to)) in aln_segs {
        if te < ts {
            continue;
        };
        if qo != to {
            continue;
        };
        if ts > last_te {
            last_ts = last_te;
            last_te = te;

            last_qs = last_qe;
            last_qe = qe;
            if last_ts == last_te {
                continue;
            }
            rtn.push(((last_ts, last_te), (last_qs, last_qe)));
        }
    }
    rtn
}

fn filter_aln_rev(aln_segs: &AlignSegments) -> Vec<((u32, u32), (u32, u32))> {
    // the aln_segs should be sorted already
    let mut aln_segs = aln_segs.clone();
    aln_segs.reverse();

    let mut last_ts = aln_segs[0].1 .0;
    let mut last_te = aln_segs[0].1 .1;

    let mut last_qs = aln_segs[0].0 .0;
    let mut last_qe = aln_segs[0].0 .1;

    let mut rtn = Vec::<((u32, u32), (u32, u32))>::new();
    rtn.push(((last_ts, last_te), (last_qs, last_qe)));
    for ((qs, _qe, qo), (ts, te, to)) in aln_segs {
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
            rtn.push(((last_ts, last_te), (last_qs, last_qe)));
        }
    }
    rtn
}

type AlignmentResult = Vec<(u32, u32, char, String, String)>;
pub fn get_variant_segments(
    target_str: &[u8],
    query_str: &[u8],
    left_padding: usize,
    max_wf_length: Option<u32>,
    mismatch_penalty: i32,
    open_penalty: i32,
    extension_penalty: i32,
) -> Option<AlignmentResult> {
    let set_len_diff = (query_str.len() as i64 - target_str.len() as i64).unsigned_abs() as u32;
    let max_wf_length = if let Some(max_wf_length) = max_wf_length {
        max_wf_length
    } else {
        std::cmp::max(2 * set_len_diff, 128_u32)
    };

    // we need to reverse the string for alignment such the the gaps are on the left
    // maybe we can do this in the stack for performance in the future
    // we assume the left_padding base on the left side are identical
    let mut r_t_str = target_str[left_padding..].to_vec();
    let mut r_q_str = query_str[left_padding..].to_vec();
    r_t_str.reverse();
    r_q_str.reverse();
    let r_t_str = String::from_utf8_lossy(&r_t_str[..]);
    let r_q_str = String::from_utf8_lossy(&r_q_str[..]);
    let t_len_minus_one = left_padding as u32 + r_t_str.len() as u32 - 1;
    let q_len_minus_one = left_padding as u32 + r_q_str.len() as u32 - 1;

    if let Some((aln_target_str, aln_query_str)) = aln::wfa_align_bases(
        &r_t_str,
        &r_q_str,
        max_wf_length,
        mismatch_penalty,
        open_penalty,
        extension_penalty,
    ) {
        let mut aln_pairs = aln::wfa_aln_pair_map(&aln_target_str, &aln_query_str);
        // assume the base on the left are identical  ( # of base = left_padding)
        (0..left_padding).for_each(|delta| {
            aln_pairs.push((
                (r_t_str.len() + delta) as u32,
                (r_q_str.len() + delta) as u32,
                'M',
            ));
        });
        // convert the coordinate from the reverse to the forward sequence
        aln_pairs.iter_mut().for_each(|(t_pos, q_pos, _c)| {
            *t_pos = t_len_minus_one - *t_pos;
            *q_pos = q_len_minus_one - *q_pos;
        });
        aln_pairs.reverse();

        // compute the VCF like variant representation
        let target_str = String::from_utf8_lossy(target_str);
        let query_str = String::from_utf8_lossy(query_str);
        Some(
            aln::get_variants_from_aln_pair_map(&aln_pairs, &target_str, &query_str)
        )
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

    let mut out_alnmap = BufWriter::new(
        File::create(Path::new(&args.output_prefix).with_extension("alnmap")).unwrap(),
    );

    let mut out_vcf =
        BufWriter::new(File::create(Path::new(&args.output_prefix).with_extension("vcf")).unwrap());

    let mut out_svcnd = BufWriter::new(
        File::create(Path::new(&args.output_prefix).with_extension("svcnd.bed")).unwrap(),
    );

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

    let query_name =  query_seqs.iter().enumerate().map(|(idx, seq_rec)| 
        (idx as u32, String::from_utf8_lossy(&seq_rec.id[..]).to_string()) ).collect::<FxHashMap<_,_>>();
    let target_name = ref_seq_index_db
        .seq_info
        .as_ref()
        .unwrap().iter().map(|(k, v)| (*k, v.0.clone())).collect::<FxHashMap<_,_>>();

    let all_records = query_seqs
        .into_par_iter()
        .enumerate()
        .map(|(q_idx, seq_rec)| {
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
            (q_idx, seq_rec, query_results)
        })
        .flat_map(|(q_idx, seq_rec, query_results)| {
            if let Some(qr) = query_results {
                let query_seq = seq_rec.seq;
                let q_len = query_seq.len();
                let mut sid_to_mapped_regions = FxHashMap::default();
                qr.into_iter().for_each(|(t_idx, mapped_segments)| {
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
                            let e = sid_to_mapped_regions.entry(t_idx).or_insert_with(Vec::new);
                            e.push((aln, orientation))
                        }
                    })
                });

                let rtn = sid_to_mapped_regions
                    .into_iter()
                    .flat_map(|(t_idx, mapped_regions)| {
                        let ref_seq = ref_seq_index_db.get_seq_by_id(t_idx).unwrap();
                        let mapped_region_aln = mapped_regions
                            .into_par_iter()
                            .map(|(aln_segs, orientation)| {
                                let aln_segs = if orientation == 0 {
                                    filter_aln(&aln_segs)
                                } else {
                                    filter_aln_rev(&aln_segs)
                                };

                                aln_segs
                                    .into_iter()
                                    .map(|((ts, te), (qs, qe))| {
                                        let ts = ts - kmer_size; // add one to ensure a match base if the first call is deletion
                                        let te = te;
                                        let qs = if orientation == 0 { qs - kmer_size } else { qs };
                                        let qe = if orientation == 0 { qe } else { qe + kmer_size };
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
                                                //let s0str = String::from_utf8_lossy(&s0str[..]);
                                                //let s1str = String::from_utf8_lossy(&s1str[..]);
                                                get_variant_segments(
                                                    &s0str,
                                                    &s1str,
                                                    1,
                                                    Some(128),
                                                    3,
                                                    2,
                                                    1,
                                                )
                                            } else {
                                                None
                                            }
                                        } else {
                                            None
                                        };

                                        // if diff.is_some() {
                                        //     let diff = diff.clone().unwrap().0;
                                        //     if !diff.is_empty() && diff[0].0 == 0 && diff[0].2 == 'I' {
                                        //         let s0str = String::from_utf8_lossy(&s0str[..]);
                                        //         let s1str = String::from_utf8_lossy(&s1str[..]);
                                        //         println!("XX0: {}",s0str);
                                        //         println!("XX1: {}",s1str);
                                        //         println!("XX2: {} {}", diff[0].3, diff[0].4);
                                        //     }
                                        // }
                                        // println!("{:?} {:?}",  ((ts, te), (qs, qe), orientation), diff);
                                        // println!();

                                        ((ts, te), (qs, qe), orientation, diff)
                                    })
                                    .collect::<Vec<_>>()
                            })
                            .filter(|v| !v.is_empty())
                            .collect::<Vec<_>>();

                        mapped_region_aln
                            .into_iter()
                            .map(|v| {
                                let mut output_records = Vec::<Record>::new();
                                let bgn_rec = v[0].clone();
                                output_records.push(Record::Bgn(
                                    t_idx,
                                    bgn_rec.0 .0,
                                    bgn_rec.0 .1,
                                    q_idx as u32,
                                    bgn_rec.1 .0,
                                    bgn_rec.1 .1,
                                    bgn_rec.2,
                                    q_len as u32,
                                ));
                                let v_last = v.last().unwrap().clone();
                                v.into_iter().for_each(
                                    |((ts, te), (qs, qe), orientation, diff)| {
                                        let qs = if orientation == 0 { qs } else { qs - kmer_size };
                                        let qe = if orientation == 0 { qe } else { qe - kmer_size };
                                        if let Some(diff) = diff {
                                            if diff.is_empty() {
                                                output_records.push(Record::Match(
                                                    t_idx,
                                                    ts,
                                                    te,
                                                    q_idx as u32,
                                                    qs,
                                                    qe,
                                                    orientation,
                                                ))
                                            } else {
                                                diff.into_iter().for_each(
                                                    |(td, qd, vt, t_str, q_str)| {
                                                        output_records.push(Record::Variant(
                                                            t_idx,
                                                            ts,
                                                            te,
                                                            q_idx as u32,
                                                            qs,
                                                            qe,
                                                            orientation,
                                                            td,
                                                            qd,
                                                            ts + td,
                                                            vt,
                                                            t_str,
                                                            q_str,
                                                        ));
                                                    },
                                                )
                                            }
                                        } else {
                                            output_records.push(Record::SvCnd(
                                                t_idx,
                                                ts,
                                                te,
                                                q_idx as u32,
                                                qs,
                                                qe,
                                                orientation,
                                            ));
                                        }
                                    },
                                );
                                output_records.push(Record::End(
                                    t_idx,
                                    v_last.0 .0,
                                    v_last.0 .1,
                                    q_idx as u32,
                                    v_last.1 .0,
                                    v_last.1 .1,
                                    v_last.2,
                                    q_len as u32,
                                ));
                                output_records
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>();
                Some(rtn)
            } else {
                None
            }
        })
        .collect::<Vec<_>>();

    writeln!(out_vcf, "##fileformat=VCFv4.2").expect("fail to write the vcf file");
    writeln!(
        out_vcf,
        r#"##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">"#
    )
    .expect("fail to write the vcf file");
    writeln!(
        out_vcf,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS000"
    )
    .expect("fail to write the vcf file");

    let mut vcf_records = Vec::<(u32, u32, String, String)>::new();
    let mut bed_records = Vec::<(u32, u32, u32)>::new();

    all_records
        .into_iter()
        .flatten()
        .enumerate()
        .for_each(|(idx, vr)| {
            vr.into_iter().for_each(|r| {
                let rec_out = match r.clone() {
                    Record::Bgn(t_idx, ts, te, q_idx, qs, qe, orientation, q_len) => {
                        let tn = target_name.get(&t_idx).unwrap();
                        let qn = query_name.get(&q_idx).unwrap(); 
                        format!(
                        "{:06}\tB\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        idx, tn, ts, te, qn, qs, qe, orientation, q_len
                    )},
                    Record::End(t_idx, ts, te, q_idx, qs, qe, orientation, q_len) => {
                        let tn = target_name.get(&t_idx).unwrap();
                        let qn = query_name.get(&q_idx).unwrap(); 
                        format!(
                        "{:06}\tE\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        idx, tn, ts, te, qn, qs, qe, orientation, q_len
                    )},
                    Record::Match(t_idx, ts, te, q_idx, qs, qe, orientation) => {
                        let tn = target_name.get(&t_idx).unwrap();
                        let qn = query_name.get(&q_idx).unwrap();
                        format!(
                        "{:06}\tM\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        idx, tn, ts, te, qn, qs, qe, orientation
                    )},
                    Record::SvCnd(t_idx, ts, te, q_idx, qs, qe, orientation) => {
                        bed_records.push((t_idx, ts + 1, te + 1));
                        let tn = target_name.get(&t_idx).unwrap();
                        let qn = query_name.get(&q_idx).unwrap();

                        format!(
                            "{:06}\tS\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            idx, tn, ts, te, qn, qs, qe, orientation
                        )
                    }
                    Record::Variant(
                        t_idx,
                        ts,
                        te,
                        q_idx,
                        qs,
                        qe,
                        orientation,
                        td,
                        qd,
                        tc,
                        vt,
                        tvs,
                        qvs,
                    ) => {
                        vcf_records.push((t_idx, tc + 1, tvs.clone(), qvs.clone()));
                        let tn = target_name.get(&t_idx).unwrap();
                        let qn = query_name.get(&q_idx).unwrap();
                        format!(
                            "{:06}\tV\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            idx, tn, ts, te, qn, qs, qe, orientation, td, qd, tc, vt, tvs, qvs
                        )
                    }
                };
                writeln!(out_alnmap, "{}", rec_out).expect("fail to write the output file");
            })
        });

    vcf_records.sort();
    vcf_records.into_iter().for_each(|(t_idx, tc, tvs, qvs)| {
        let tn =  target_name.get(&t_idx).unwrap();
        writeln!(
            out_vcf,
            "{}\t{}\t.\t{}\t{}\t60\tPASS\tGT\t1/1",
            tn, tc, tvs, qvs
        )
        .expect("fail to write the vcf file");
    });

    bed_records.sort();
    bed_records.into_iter().for_each(|(t_idx, ts, te)| {
        let tn =  target_name.get(&t_idx).unwrap();
        writeln!(out_svcnd, "{}\t{}\t{}", tn, ts + 1, te + 1)
            .expect("fail to write the 'in-alignment' sv candidate bed file");
    });

    Ok(())
}
