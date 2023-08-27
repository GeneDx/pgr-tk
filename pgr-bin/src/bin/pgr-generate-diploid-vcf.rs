const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use iset::set::IntervalSet;
// use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

/// Generate diploid VCF field from paired alnmap file from two haplotype assembly
#[derive(Parser, Debug)]
#[clap(name = "pgr-generate-diploid-vcf")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
    /// path to the first haplotype alnmap file
    hap0_path: String,
    /// path to the second haplotype alnmap file
    hap1_path: String,
    /// path to a ctgmap.json file
    target_len_json_path: String,
    /// the prefix of the output files
    output_path: String,
    /// number of threads used in parallel (more memory usage), default to "0" using all CPUs available or the number set by RAYON_NUM_THREADS
    #[clap(long, default_value_t = 0)]
    number_of_thread: usize,
}

type TargetSeqLength = Vec<(u32, String, u32)>;

type ShimmerMatchBlock = (String, u32, u32, String, u32, u32, u32);
type VariantRecord = (String, u32, u32, u8, String, String);

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.number_of_thread)
        .build_global()
        .unwrap();

    let mut target_length_json_file = BufReader::new(
        File::open(Path::new(&args.target_len_json_path)).expect("can't open the input file"),
    );
    let mut buffer = Vec::new();
    target_length_json_file.read_to_end(&mut buffer)?;
    let mut target_length: TargetSeqLength =
        serde_json::from_str(&String::from_utf8_lossy(&buffer[..]))
            .expect("can't parse the target_len.json file");

    target_length.sort();

    let hap0_alnmap_file = BufReader::new(File::open(Path::new(&args.hap0_path)).unwrap());

    let hap1_alnmap_file = BufReader::new(File::open(Path::new(&args.hap1_path)).unwrap());

    #[allow(clippy::type_complexity)]
    let get_variant_recs = |f: BufReader<File>,
                            hap_type: u8|
     -> (
        Vec<(String, u32, u32, u8, String, String)>,
        FxHashMap<u32, (Option<ShimmerMatchBlock>, Option<ShimmerMatchBlock>)>,
    ) {
        let mut out = Vec::<(String, u32, u32, u8, String, String)>::new();
        let mut aln_block =
            FxHashMap::<u32, (Option<ShimmerMatchBlock>, Option<ShimmerMatchBlock>)>::default();

        f.lines().for_each(|line| {
            if let Ok(line) = line {
                let fields = line.split('\t').collect::<Vec<&str>>();
                assert!(fields.len() > 3);
                if fields[1] == "V" {
                    assert!(fields.len() == 15);
                    let err_msg = format!("fail to parse on {}", line);
                    let t_name = fields[2];
                    // let ts = fields[3].parse::<u32>().expect(&err_msg);
                    // let te = fields[4].parse::<u32>().expect(&err_msg);
                    // let q_name = fields[5];
                    // let qs = fields[6].parse::<u32>().expect(&err_msg);
                    // let qe = fields[7].parse::<u32>().expect(&err_msg);
                    // let orientation = fields[8].parse::<u32>().expect(&err_msg);
                    // let td = fields[9].parse::<u32>().expect(&err_msg);
                    // let qd = fields[10].parse::<u32>().expect(&err_msg);
                    let tc = fields[11].parse::<u32>().expect(&err_msg);
                    // let tt = fields[12].chars().next().expect(&err_msg);
                    let tvs = fields[13];
                    let qvs = fields[14];
                    out.push((
                        t_name.to_string(),
                        tc,
                        tvs.len() as u32,
                        hap_type,
                        tvs.to_string(),
                        qvs.to_string(),
                    ));
                } else if fields[1] == "B" || fields[1] == "E" {
                    let err_msg = format!("fail to parse on {}", line);
                    let aln_block_id = fields[0].parse::<u32>().expect(&err_msg);
                    let t_name = fields[2];
                    let ts = fields[3].parse::<u32>().expect(&err_msg);
                    let te = fields[4].parse::<u32>().expect(&err_msg);
                    let q_name = fields[5];
                    let qs = fields[6].parse::<u32>().expect(&err_msg);
                    let qe = fields[7].parse::<u32>().expect(&err_msg);
                    let orientation = fields[8].parse::<u32>().expect(&err_msg);
                    let e = aln_block.entry(aln_block_id).or_default();
                    match fields[1] {
                        "B" => {
                            e.0 = Some((
                                t_name.to_string(),
                                ts,
                                te,
                                q_name.to_string(),
                                qs,
                                qe,
                                orientation,
                            ));
                        }
                        "E" => {
                            e.1 = Some((
                                t_name.to_string(),
                                ts,
                                te,
                                q_name.to_string(),
                                qs,
                                qe,
                                orientation,
                            ));
                        }
                        _ => (),
                    }
                }
            }
        });
        (out, aln_block)
    };
    let (hap0_recs, hap0_aln_block) = get_variant_recs(hap0_alnmap_file, 0);
    let (hap1_recs, hap1_aln_block) = get_variant_recs(hap1_alnmap_file, 1);

    let mut hap0_aln_interval = FxHashMap::<String, IntervalSet<u32>>::default();
    hap0_aln_block.into_iter().for_each(|(_id, rec_pair)| {
        if let (Some(b_rec), Some(e_rec)) = rec_pair {
            let t_name = b_rec.0;
            let bgn = b_rec.1;
            let end = e_rec.2;
            let interval_set = hap0_aln_interval.entry(t_name).or_default();
            interval_set.insert(bgn..end);
        }
    });

    let mut hap1_aln_interval = FxHashMap::<String, IntervalSet<u32>>::default();
    hap1_aln_block.into_iter().for_each(|(_id, rec_pair)| {
        if let (Some(b_rec), Some(e_rec)) = rec_pair {
            let t_name = b_rec.0;
            let bgn = b_rec.1;
            let end = e_rec.2;
            let interval_set = hap1_aln_interval.entry(t_name).or_default();
            interval_set.insert(bgn..end);
        }
    });

    let mut out_vcf = BufWriter::new(File::create(Path::new(&args.output_path)).unwrap());
    writeln!(out_vcf, "##fileformat=VCFv4.2").expect("fail to write the vcf file");
    target_length.into_iter().for_each(|(_, t_name, t_len)| {
        writeln!(out_vcf, r#"##contig=<ID={},length={}>"#, t_name, t_len)
            .expect("fail to write the vcf file");
    });
    writeln!(
        out_vcf,
        r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#
    )
    .expect("fail to write the vcf file");
    writeln!(
        out_vcf,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    )
    .expect("fail to write the vcf file");

    let convert_to_vcf_record = |records: &Vec<VariantRecord>| {
        let mut ref_bases = FxHashSet::<(u32, char)>::default();
        let mut h0alleles = Vec::<(u32, VariantRecord)>::new();
        let mut h1alleles = Vec::<(u32, VariantRecord)>::new();
        let mut al_idx_map = FxHashMap::<(u32, String, String), u32>::default();
        let mut al_idx = 0_u32;

        let ref_name = records.first().unwrap().0.clone();
        records.iter().for_each(|rec| {
            let (_t_name, ts, tl, ht, vts, vqs) = rec;
            (0..*tl).for_each(|t_pos| {
                let vts = vts.chars().collect::<Vec<_>>();
                ref_bases.insert((*ts + t_pos, vts[t_pos as usize]));
            });

            let key = (*ts, vts.clone(), vqs.clone());

            al_idx_map.entry(key).or_insert_with(|| {
                al_idx += 1;
                al_idx
            });

            if *ht == 0 {
                h0alleles.push((al_idx, rec.clone()));
            };
            if *ht == 1 {
                h1alleles.push((al_idx, rec.clone()));
            };
        });
        let mut ref_bases = ref_bases.into_iter().collect::<Vec<_>>();
        ref_bases.sort();
        let ref_str = String::from_iter(ref_bases.iter().map(|(_, c)| *c).collect::<Vec<_>>());
        assert!(ref_str.len() == ref_bases.len()); // make sure all bases at the same t_pos are the same, if not the vectors will have different lengths
        let ts0 = ref_bases.first().unwrap().0;
        let tl0 = ref_str.len() as u32;

        let mut query_alleles = al_idx_map
            .iter()
            .map(|((ts, tvs, qvs), &al_idx)| {
                let t_prefix = ref_str[0..(ts - ts0) as usize].to_string();
                let t_suffix = ref_str[(ts + tvs.len() as u32 - ts0) as usize..].to_string();
                (al_idx, [t_prefix, qvs.clone(), t_suffix].join(""))
            })
            .collect::<Vec<_>>();

        query_alleles.sort();
        let query_alleles = query_alleles
            .iter()
            .map(|(_, qs)| qs.clone())
            .collect::<Vec<_>>()
            .join(",");

        let h0_al_idx = if let Some(i_set) = hap0_aln_interval.get(&ref_name) {
            if i_set.has_overlap(ts0..ts0 + tl0) {
                if h0alleles.is_empty() {
                    "0".to_string()
                } else {
                    format!("{}", h0alleles.last().unwrap().0)
                }
            } else {
                ".".to_string()
            }
        } else {
            ".".to_string()
        };
        let h1_al_idx = if let Some(i_set) = hap1_aln_interval.get(&ref_name) {
            if i_set.has_overlap(ts0..ts0 + tl0) {
                if h1alleles.is_empty() {
                    "0".to_string()
                } else {
                    format!("{}", h1alleles.last().unwrap().0)
                }
            } else {
                ".".to_string()
            }
        } else {
            ".".to_string()
        };
        let gt = [h0_al_idx, h1_al_idx].join("|");
        (ref_name, ts0, ref_str, query_alleles, gt)
    };

    let mut variant_records = Vec::<VariantRecord>::new();
    variant_records.extend(hap0_recs);
    variant_records.extend(hap1_recs);

    let mut variant_records = variant_records.into_iter().collect::<Vec<_>>();
    // variant_group: represent a group of overlapped variants, ref_id, ref_start, len, REF, ALT
    let mut variant_group = Vec::<VariantRecord>::new();
    // currrent_vg_end: represent the end coordinate of the current variant group
    let mut currrent_vg_end = Option::<(String, u32)>::None;
    variant_records.sort();
    variant_records
        .into_iter()
        .for_each(|(ref_name, ts, tl, ht, vts, vqs)| {
            if let Some(currrent_vg_end) = currrent_vg_end.clone() {
                //println!("{} {} {} {} {:?} {}", ref_name, ts, tl, ts + tl ,  currrent_vg_end, variant_group.len()  );
                if ref_name == currrent_vg_end.0 && ts < currrent_vg_end.1 {
                    variant_group.push((ref_name.clone(), ts, tl, ht, vts, vqs));
                } else if !variant_group.is_empty() {
                    // println!("X {} {} {} {}", ref_name, ts, tl, variant_group.len());
                    let (vcfrec_ref_name, ts0, ref_str, query_alleles, gt) =
                        convert_to_vcf_record(&variant_group);
                    writeln!(
                        out_vcf,
                        "{}\t{}\t.\t{}\t{}\t60\tPASS\t.\tGT\t{}",
                        vcfrec_ref_name,
                        ts0 + 1,
                        ref_str,
                        query_alleles,
                        gt,
                    )
                    .expect("fail to write the vcf file");
                    variant_group.clear();
                    variant_group.push((ref_name.clone(), ts, tl, ht, vts, vqs));
                }
            } else {
                variant_group.push((ref_name.clone(), ts, tl, ht, vts, vqs));
                currrent_vg_end = Some((ref_name, ts + tl));
                return;
            }
            currrent_vg_end = Some((ref_name.clone(), ts + tl));
        });
    if !variant_group.is_empty() {
        // println!("X {} {} {} {}", ref_name, ts, tl, variant_group.len());
        let (vcfrec_ref_name, ts0, ref_str, query_alleles, gt) =
            convert_to_vcf_record(&variant_group);
        writeln!(
            out_vcf,
            "{}\t{}\t.\t{}\t{}\t60\tPASS\t.\tGT\t{}",
            vcfrec_ref_name,
            ts0 + 1,
            ref_str,
            query_alleles,
            gt,
        )
        .expect("fail to write the vcf file");
    };
    Ok(())
}
