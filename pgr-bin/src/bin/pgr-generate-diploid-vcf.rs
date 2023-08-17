const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Align long contigs and identify potential SV regions with respect to the reference fasta file
#[derive(Parser, Debug)]
#[clap(name = "pgr-generate-diploid-vcf")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
    /// path to the the first haplotype alnmap file
    hap0_path: String,
    /// path to the the second haplotype alnmap file
    hap1_path: String,
    /// the prefix of the output files
    output_prefix: String,
    /// number of threads used in parallel (more memory usage), default to "0" using all CPUs available or the number set by RAYON_NUM_THREADS
    #[clap(long, default_value_t = 0)]
    number_of_thread: usize,
}

type ShimmerMatchBlock = (String, u32, u32, String, u32, u32, u32);
type VariantRecord = (ShimmerMatchBlock, u32, u32, u32, char, String, String);

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.number_of_thread)
        .build_global()
        .unwrap();

    let mut hap0_alnmap_file = BufReader::new(File::open(Path::new(&args.hap0_path)).unwrap());

    let mut hap1_alnmap_file = BufReader::new(File::open(Path::new(&args.hap1_path)).unwrap());

    let get_variant_recs = |f: BufReader<File>| -> FxHashMap<(String, u32), Vec<VariantRecord>> {
        let mut out = FxHashMap::<(String, u32), Vec<VariantRecord>>::default();
        f.lines().for_each(|line| {
            if let Ok(line) = line {
                let fields = line.split('\t').collect::<Vec<&str>>();
                assert!(fields.len() > 3);
                if fields[1] == "V" {
                    assert!(fields.len() == 15);
                    let err_msg = format!("faile to parse on {}", line);
                    let t_name = fields[2];
                    let ts = fields[3].parse::<u32>().expect(&err_msg);
                    let te = fields[4].parse::<u32>().expect(&err_msg);
                    let q_name = fields[5];
                    let qs = fields[6].parse::<u32>().expect(&err_msg);
                    let qe = fields[7].parse::<u32>().expect(&err_msg);
                    let orientation = fields[8].parse::<u32>().expect(&err_msg);
                    let td = fields[9].parse::<u32>().expect(&err_msg);
                    let qd = fields[10].parse::<u32>().expect(&err_msg);
                    let tc = fields[11].parse::<u32>().expect(&err_msg);
                    let tt = fields[12].chars().next().expect(&err_msg);
                    let tvs = fields[13];
                    let qvs = fields[14];
                    let key = (t_name.to_string(), tc);
                    let e = out.entry(key).or_insert_with(Vec::new);
                    e.push((
                        (
                            t_name.to_string(),
                            ts,
                            te,
                            q_name.to_string(),
                            qs,
                            qe,
                            orientation,
                        ),
                        td,
                        qd,
                        tc,
                        tt,
                        tvs.to_string(),
                        qvs.to_string(),
                    ))
                }
            }
        });
        out
    };
    let hap0_recs = get_variant_recs(hap0_alnmap_file);
    let hap1_recs = get_variant_recs(hap1_alnmap_file);

    let mut out_vcf =
        BufWriter::new(File::create(Path::new(&args.output_prefix).with_extension("vcf")).unwrap());
    writeln!(out_vcf, "##fileformat=VCFv4.2").expect("fail to write the vcf file");
    writeln!(
        out_vcf,
        r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#
    )
    .expect("fail to write the vcf file");
    writeln!(
        out_vcf,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tGT\tSAMPLE"
    )
    .expect("fail to write the vcf file");

    let mut rec_keys = FxHashSet::<(String, u32)>::default();
    rec_keys.extend(hap0_recs.keys().cloned());
    rec_keys.extend(hap1_recs.keys().cloned());

    let mut rec_keys = rec_keys.into_iter().collect::<Vec<_>>();
    rec_keys.sort();
    rec_keys.into_iter().for_each(|k| {
        let mut vt_count = FxHashMap::<String, FxHashMap<String, (u32, u32)>>::default();
        if let Some(v0) = hap0_recs.get(&k) {
            v0.iter().for_each(|(_, _, _, _, _, vts, vqs)| {
                let vts = vts.trim_end_matches('-').to_string();
                let e = vt_count
                    .entry(vts)
                    .or_insert(FxHashMap::<String, (u32, u32)>::default());
                let e2 = e.entry(vqs.clone()).or_insert((0, 0));
                e2.0 += 1;
            });
        };

        if let Some(v1) = hap1_recs.get(&k) {
            v1.iter().for_each(|(_, _, _, _, _, vts, vqs)| {
                let vts = vts.trim_end_matches('-').to_string();
                let e = vt_count
                    .entry(vts)
                    .or_insert(FxHashMap::<String, (u32, u32)>::default());
                let e2 = e.entry(vqs.clone()).or_insert((0, 0));
                e2.1 += 1;
            });
        };
        // TODO: we still have to check missing call
        vt_count.into_iter().for_each(|(vts, gt_count)| {
            
            gt_count.into_iter().for_each(| (qvs, (c0, c1)) | {
                let qvs = qvs.trim_end_matches('-').clone();
                let gt0 = if c0 > 0 { 1 } else { 0 };
                let gt1 = if c1 > 0 { 1 } else { 0 };
                let gt = format!("{}|{}", gt0, gt1);
                writeln!(
                    out_vcf,
                    "{}\t{}\t.\t{}\t{}\t60\tPASS\t.\tGT\t{}",
                    k.0,
                    k.1+1,
                    vts,
                    qvs,
                    gt,
                )
                .expect("fail to write the vcf file");
            });

        });
    });
    Ok(())
}
