#![allow(dead_code)]

use super::seq_db::SeqDB;
use super::shmmrutils::MM128;
use core::ops::Range;
use intervaltree::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
pub type MapIntervalRecord = [u32; 8];
pub type MapIntervals = FxHashMap<u32, IntervalTree<u32, MapIntervalRecord>>;
pub type Shmmrs = Vec<Vec<MM128>>;
use petgraph::graphmap::DiGraphMap;
use petgraph::unionfind::UnionFind;
use rayon::{iter::IntoParallelRefIterator, prelude::*};


fn load_shmmr_map(mi: &mut MapIntervals, filename: &String) -> Result<(), std::io::Error> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut itvl_vecs = FxHashMap::<u32, Vec<(Range<u32>, MapIntervalRecord)>>::default();

    for r in reader.lines() {
        let r = r?;
        let r: Vec<&str> = r.split_whitespace().collect();
        if r[0] != "M" {
            continue;
        }
        let mut v: MapIntervalRecord = [0; 8];
        for i in 1..v.len() + 1 {
            v[i - 1] = r[i].parse::<u32>().unwrap();
        }
        itvl_vecs
            .entry(v[0])
            .or_insert_with(|| vec![])
            .push((v[1]..v[2], v.clone()));
    }
    for (sid, itvl_vec) in itvl_vecs {
        let itvl: IntervalTree<u32, MapIntervalRecord> = itvl_vec.iter().cloned().collect();
        mi.insert(sid, itvl);
    }
    Ok(())
}

fn write_shmmr_map(mi: &MapIntervals, filename: &String) -> Result<(), std::io::Error> {
    let mut file = BufWriter::new(File::create(filename)?);

    for (_rid, itvls) in mi.iter() {
        for itvl in itvls.iter_sorted() {
            let v = itvl.clone().value;
            writeln!(
                file,
                "M {} {} {} {} {} {} {} {}",
                v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]
            )?;
        }
    }
    Ok(())
}

fn generate_shmmr_map(shmmrs0: &Shmmrs, shmmrs1: &Shmmrs, max_hits: usize) -> MapIntervals {
    let out_vec = shmmrs1
        .par_chunks(250)
        .into_par_iter()
        .flat_map_iter(|mmers_v| {
            mmers_v
                .iter()
                .filter(|mmers| mmers.len() >= 2)
                .flat_map(|mmers| {
                    let mlen = mmers.len();
                    mmers[0..mlen - 1]
                        .iter()
                        .zip(mmers[1..mlen].iter())
                        .filter(|m| (m.0.x >> 8) != (m.1.x >> 8))
                        .map(|m| {
                            let (m0, m1) = m;
                            let hash0 = (m0.x >> 8) as u128;
                            let hash1 = (m1.x >> 8) as u128;
                            let h128 = (hash0 << 64) | hash1;
                            (h128, m0.y, m1.y)
                        })
                })
        })
        .collect::<Vec<(u128, u64, u64)>>();

    const INDEX_CHUNKS: usize = 16;

    let mut smp_index = Vec::<(u32, FxHashMap<u128, Vec<(u64, u64)>>)>::new();
    for i in 0..INDEX_CHUNKS {
        smp_index.push((i as u32, FxHashMap::<u128, Vec<(u64, u64)>>::default()));
    }
    smp_index.par_iter_mut().for_each(|x| {
        let (idx, idxhash) = x;
        out_vec
            .iter()
            .filter(|&h| ((h.0 & 0xFF) % INDEX_CHUNKS as u128) as u8 == *idx as u8)
            .for_each(|v| {
                idxhash
                    .entry(v.0)
                    .or_insert_with(|| vec![])
                    .push((v.1, v.2))
            })
    });

    let mut all_itvl = MapIntervals::default();

    for (rid0, mmers) in shmmrs0.iter().enumerate() {
        if mmers.len() < 2 {
            continue;
        }

        let mut itvl_vec = Vec::<(Range<u32>, MapIntervalRecord)>::new();

        let mlen = mmers.len();
        if mlen < 2 {
            continue;
        }
        let out_v = mmers[0..mlen - 1]
            .par_iter()
            .zip(mmers[1..mlen].par_iter())
            .filter(|m| (m.0.x >> 8) != (m.1.x >> 8))
            .map(|m| -> [(u128, u64, u64, u32); 2] {
                let (m0, m1) = m;
                let hash0 = (m0.x >> 8) as u128;
                let hash1 = (m1.x >> 8) as u128;
                let h128 = (hash0 << 64) | hash1;
                let rh128 = (hash1 << 64) | hash0;
                [(h128, m0.y, m1.y, 0), (rh128, m0.y, m1.y, 1)]
            })
            .collect::<Vec<[(u128, u64, u64, u32); 2]>>();

        //.collect_into_vec(&mut out_v);

        let mut range_matches = Vec::<Vec<(Range<u32>, MapIntervalRecord)>>::new();
        out_v
            .par_iter()
            .map(|pv| {
                let mut rmv = Vec::<(Range<u32>, MapIntervalRecord)>::new();
                for v in pv {
                    let (h128, y0, y1, orientation) = *v;
                    let shard = ((h128 & 0x0F) % INDEX_CHUNKS as u128) as usize;
                    if !smp_index[shard].1.contains_key(&h128) {
                        continue;
                    }
                    let matches = smp_index[shard].1.get(&h128).unwrap();
                    if (*matches).len() > max_hits {
                        continue;
                    }
                    for (yy0, yy1) in matches.iter() {
                        let rid0 = (y0 >> 32) as u32;
                        let pos0 = ((y0 & 0xFFFFFFFF) >> 1) as u32;
                        let pos1 = ((y1 & 0xFFFFFFFF) >> 1) as u32;

                        let rid1 = (*yy0 >> 32) as u32;
                        let ppos0 = ((*yy0 & 0xFFFFFFFF) >> 1) as u32;
                        let ppos1 = ((*yy1 & 0xFFFFFFFF) >> 1) as u32;
                        let group = 0;
                        rmv.push((
                            pos0..pos1,
                            [rid0, pos0, pos1, rid1, ppos0, ppos1, orientation, group],
                        ));
                        //DEBUG
                        //println!("{} {} {} {} {} {} {} {}", rid0, pos0, pos1, rid1, ppos0, ppos1, orientation, group);
                    }
                }
                rmv
            })
            .collect_into_vec(&mut range_matches);

        for rmv in range_matches {
            itvl_vec.extend(rmv);
        }

        let itvl: IntervalTree<u32, MapIntervalRecord> = itvl_vec.iter().cloned().collect();
        all_itvl.insert(rid0 as u32, itvl);
    }
    all_itvl
}

fn build_shmmer_map_from_query_results(mqr: &Vec<MapIntervalRecord>) -> MapIntervals {
    let mut id2itvls = FxHashMap::<u32, Vec<(Range<u32>, MapIntervalRecord)>>::default();

    mqr.iter().for_each(|&r| {
        id2itvls
            .entry(r[0])
            .or_insert_with(|| vec![])
            .push((r[1]..r[2], r));
    });

    let mut all_itvl = MapIntervals::default();
    id2itvls.iter().for_each(|(sid, recs)| {
        let itvl: IntervalTree<u32, MapIntervalRecord> = recs.iter().cloned().collect();
        all_itvl.insert(*sid, itvl);
    });

    all_itvl
}

fn map_interval_query(
    ivtl: &MapIntervals,
    sid: u32,
    bgn: u32,
    end: u32,
) -> Vec<MapIntervalRecord> {
    let mut q_res: Vec<_> = ivtl
        .get(&sid)
        .unwrap()
        .query(bgn..end)
        .map(|x| x.value)
        .collect();
    q_res.sort();
    let mut mq = Vec::<MapIntervalRecord>::new();
    for e in q_res {
        mq.push(e);
        // DBEUG
        // println!("q: {} {} {} r:{} {} {} {} {} {} {}", sid, bgn, end, e[0], e[1], e[2], e[3], e[4], e[5], e[6]);
    }
    mq
}

pub fn find_match_chain(matches: &Vec<MapIntervalRecord>) -> Vec<MapIntervalRecord> {
    let mut seqpair_count = FxHashMap::<(u32, u32), u32>::default();
    let mut out = Vec::<MapIntervalRecord>::new();
    matches.iter().for_each(|v| {
        let sid0 = v[0];
        let sid1 = v[3];
        *seqpair_count.entry((sid0, sid1)).or_insert(0) += 1;
    });

    for ((sid0, sid1), count) in seqpair_count.iter() {
        if *count < 4 {
            continue;
        }
        let mut aln_g = DiGraphMap::<u32, u32>::new();
        let filtered_matches = matches
            .iter()
            .filter(|x| x[0] == *sid0 && x[3] == *sid1)
            .collect::<Vec<&MapIntervalRecord>>();
        for i in 0..filtered_matches.len() {
            let mut j = i + 1;
            loop {
                if j >= filtered_matches.len() || j > i + 25 {
                    break;
                }
                let r1 = filtered_matches[i];
                let r2 = filtered_matches[j];
                if r1[1] == r2[1] {
                    j += 1;
                    continue;
                }
                let d1: i64;
                let d2: i64;
                if r1[6] == 1 {
                    d1 = r1[1] as i64 + r1[4] as i64;
                } else {
                    d1 = r1[1] as i64 - r1[4] as i64;
                }
                if r2[6] == 1 {
                    d2 = r2[1] as i64 + r2[4] as i64;
                } else {
                    d2 = r2[1] as i64 - r2[4] as i64;
                }
                let d: f32 =
                    ((d1 - d2).abs() as f32) / ((r2[1] as i64 - r1[1] as i64).abs() as f32);
                if d < 0.2 && r1[6] == r2[6] && r2[1] - r1[1] < 100000 {
                    aln_g.add_edge(i as u32, j as u32, 1);
                }
                j += 1;
            }
        }
        let mut vertex_sets = UnionFind::new(filtered_matches.len());
        for e in aln_g.all_edges() {
            vertex_sets.union(e.0, e.1);
        }

        let labels = vertex_sets.into_labeling();
        for i in 0..filtered_matches.len() {
            if !aln_g.contains_node(i as u32) {
                continue;
            }
            let label = labels[i];
            //println!("{} {} {:?} {} {}", sid0, sid1, filtered_matches[i], i, label);
            let mut out_r: [u32; 8] = filtered_matches[i].clone();
            out_r[7] = label;
            out.push(out_r);
        }
    }
    out
}

pub struct ExtDeltas {
    pub sid0: u32,
    pub pos0: u32,
    pub sid1: u32,
    pub pos1: u32,
    pub dk: i32,
    pub strand1: u8,
    pub base: char,
}

fn rc_dna_seq(seq: Vec<u8>) -> Vec<u8> {
    // support ACGTacgtNn -> TGCAtcgaNn (revesred)
    let rev_map: [u8; 256] = [
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
        25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
        48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 84, 66, 71, 68, 69, 70,
        67, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 65, 85, 86, 87, 88, 89, 90, 91, 92, 93,
        94, 95, 96, 116, 98, 103, 100, 101, 102, 99, 104, 105, 106, 107, 108, 109, 110, 111, 112,
        113, 114, 115, 97, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130,
        131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148,
        149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166,
        167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184,
        185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202,
        203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220,
        221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238,
        239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255,
    ];

    let mut rtn = Vec::<u8>::new();
    seq.iter()
        .rev()
        .for_each(|c| rtn.push(rev_map[*c as usize]));
    rtn
}

fn get_deltas(
    subseq0: &Vec<u8>,
    subseq1: &Vec<u8>,
    sid0: u32,
    sid1: u32,
    offset0: u32,
    offset1: u32,
    strand: u8,
) -> Vec<ExtDeltas> {
    let mut out_vec = Vec::<ExtDeltas>::new();
    if subseq0.len() < 32 || subseq1.len() < 32 {
        return out_vec;
    }
    //let hs0 = super::shmmrutils::get_hpc_seq(&subseq0).s;
    //let hs1 = super::shmmrutils::get_hpc_seq(&subseq1).s;
    let hs0 = subseq0;
    let hs1 = subseq1;

    if let Some(ovlp) = super::shmmrutils::match_reads(hs0, hs1, true, 0.02, 32, 24) {
        if let Some(mut deltas) = ovlp.deltas {
            if deltas.len() > 0 {
                deltas.reverse();
                for dpt in deltas {
                    let base = if dpt.dk > 0 {
                        '-'
                    } else {
                        subseq1[dpt.y as usize - 1] as char
                    };
                    out_vec.push(ExtDeltas {
                        sid0: sid0,
                        pos0: offset0 + dpt.x,
                        sid1: sid1,
                        pos1: offset1 + dpt.y,
                        strand1: strand,
                        dk: dpt.dk,
                        base: base,
                    });
                }
            }
        }
    }
    out_vec
}

fn generate_deltas(seqdb0: &SeqDB, seqdb1: &SeqDB, m: MapIntervalRecord, k: u32) -> usize {
    let mut rev = false;
    if m[6] == 1 {
        //strand
        rev = true;
    }
    let sid0 = m[0] as usize;
    let b0 = m[1] as usize;
    let e0 = m[2] as usize;
    let sid1 = m[3] as usize;
    let b1 = m[4] as usize;
    let e1 = m[5] as usize;
    let subseq0;
    let subseq1;
    let strand: u8;

    // we need to pad the sequence with the minimize on both end. the match_reads() need at least 8 base match at the beginning
    match rev {
        true => {
            subseq0 = seqdb0.seqs[sid0 as usize][(b0 - (k - 1) as usize)..e0].to_vec();
            
            let bb1 = b1 - (k - 2) as usize;
            let ee1 = e1 - (k - 2) as usize;

            let subseq_tmp = seqdb1.seqs[sid1 as usize][bb1..(ee1 + (k - 1) as usize)].to_vec();
            subseq1 = rc_dna_seq(subseq_tmp);
            strand = 1;
            /*
            println!("DEBUB s0 {}", String::from_utf8(subseq0.clone()).unwrap());
            println!("DEBUB s1 {}", String::from_utf8(subseq1.clone()).unwrap());
            */
        }
        false => {
            subseq0 = seqdb0.seqs[sid0 as usize][(b0 - (k - 1) as usize)..e0].to_vec();
            subseq1 = seqdb1.seqs[sid1 as usize][(b1 - (k - 1) as usize)..e1].to_vec();
            strand = 0;
        }
    }

 
    let v = get_deltas(
        &subseq0,
        &subseq1,
        sid0 as u32,
        sid1 as u32,
        m[1],
        m[4],
        strand,
    );

    //let sname0 = &seqdb0.id2seqname[&(sid0 as u32)];
    //let sname1 = &seqdb1.id2seqname[&(sid1 as u32)];
    /*
    v.iter().for_each(|e| {
        println!(
            "D {} {} {} {} {} {} {} {} {}",
            sname0, e.sid0, e.pos0, sname1, e.sid1, e.pos1, e.strand1, e.dk, e.base
        );
    });
    */

    v.len()
}

fn linearize_hits(mut v: Vec<[u32; 8]>) -> Vec<[u32; 8]> {
    let mut vv = Vec::<[u32; 8]>::new();

    v.sort_by(|a, b| (a[1], a[3], a[4]).partial_cmp(&(b[1], b[3], b[4])).unwrap());
    let mut strand_count = FxHashMap::<u32, [u32; 2]>::default();
    v.iter().for_each(|x| {
        let sid = x[3];
        let strand = x[6];
        let sc = strand_count.entry(sid).or_insert([0_u32; 2]);
        if strand == 0 {
            sc[0] += 1;
        } else {
            sc[1] += 1;
        }
    });

    // identify the major alignment strand between two contigs
    let mut strands = FxHashMap::<u32, u32>::default();
    strand_count.iter().for_each(|(k, v)| {
        if v[0] < v[1] {
            strands.entry(*k).or_insert(1);
        } else {
            strands.entry(*k).or_insert(0);
        }
    });

    let mut last_postions = FxHashMap::<u32, [u32; 2]>::default();
    v.iter().for_each(|x| {
        let sid = x[3];
        if !last_postions.contains_key(&sid) {
            vv.push(x.clone());
            let strand = *strands.get(&sid).unwrap();
            if x[6] != strand {
                return;
            }
            if strand == 0 {
                last_postions.insert(sid, [x[2], x[5]]);
            } else {
                last_postions.insert(sid, [x[2], x[4]]);
            }
        } else {
            let last_p = last_postions.get(&sid).unwrap();
            let strand = *strands.get(&sid).unwrap();
            if x[6] != strand {
                return;
            }
            if x[2] <= last_p[0] {
                return;
            }
            if strand == 0 {
                if x[5] <= last_p[1] {
                    return;
                } else {
                    vv.push(x.clone());
                    last_postions.insert(sid, [x[2], x[5]]);
                }
            } else {
                if x[4] >= last_p[1] {
                    return;
                } else {
                    vv.push(x.clone());
                    last_postions.insert(sid, [x[2], x[4]]);
                }
            }
        }
    });
    vv
}

fn filter_chain_group(chain_groups: FxHashMap<(u32, u32, u32), Vec<[u32; 8]>>) -> Vec<[u32; 8]> {
    let mut chain_list = chain_groups
        .iter()
        .collect::<Vec<(&(u32, u32, u32), &Vec<[u32; 8]>)>>();

    chain_list.sort_by_key(|v| -(v.1.len() as i64));

    let mut r_intervals = FxHashSet::<(u32, u32, u32)>::default();
    let mut t_intervals = FxHashSet::<(u32, u32, u32)>::default();
    let mut aln_itvls = Vec::<[u32; 8]>::new();
    for (_k, v) in chain_list.iter() {
        //DEBUG
        //println!("# {:?} {:?}", _k, v.len());
        if v.len() < 4 {
            continue;
        }
        let mut aln_itvls0 = Vec::<[u32; 8]>::new();
        for w in v.iter() {
            let r_ivtl = (w[3], w[1], w[2]); // tagged with target id for both
            let t_ivtl = (w[3], w[4], w[5]);
            if !r_intervals.contains(&r_ivtl) && !t_intervals.contains(&t_ivtl) {
                aln_itvls0.push(*w);
                r_intervals.insert(r_ivtl);
                t_intervals.insert(t_ivtl);
            }
        }
        if aln_itvls0.len() > 0 && aln_itvls0.len() as f32 > 0.5 * v.len() as f32 {
            //ad hoc rule to avoid repeat
            aln_itvls.extend(aln_itvls0);
        }
    }
    aln_itvls
}

pub fn multi_map_seqs(
    ref_fasta_file: &String,
    target_fasta_file: &String,
    output_prefix: &String,
    w: u32,
    k: u32,
    r: u32,
    max_hits: usize,
) -> Result<(), std::io::Error> {
    let mut sdb0 = SeqDB::new(ref_fasta_file.clone());
    let mut sdb1 = SeqDB::new(target_fasta_file.clone());

    sdb0.load_sequences()?;
    sdb1.load_sequences()?;
    sdb0.build_shmmrs(w, k, r);
    sdb1.build_shmmrs(w, k, r);

    let shmmrmap = generate_shmmr_map(&sdb0.shmmrs, &sdb1.shmmrs, max_hits);
    let mut out_file = BufWriter::new(File::create(format!("{}.aln", output_prefix))?);

    (0..sdb0.seqs.len())
        .into_par_iter()
        .map(|sid| {
            let seqname = &sdb0.id2name[&(sid as u32)];
            let seqlen = sdb0.lengths[seqname];
            let mq = map_interval_query(&shmmrmap, sid as u32, 0, seqlen as u32);
            let mq_filtered = find_match_chain(&mq);
            let mut chain_group = FxHashMap::<(u32, u32, u32), Vec<[u32; 8]>>::default();
            mq_filtered.iter().for_each(|r| {
                let (sid0, sid1, chain_label) = (r[0], r[3], r[7]);
                let key = (sid0, sid1, chain_label);
                //DEBUG
                /*
                println!(
                    "{} {} {} {} {} {} {} {} {}",
                    r[0],
                    r[1],
                    r[2],
                    r[3],
                    r[4],
                    r[5],
                    r[6],
                    r[7],
                    sdb0.id2seqname.get(&r[3]).unwrap()
                );
                */
                //println!("# {} {} {}", sid0, sid1, chain_label );
                chain_group.entry(key).or_insert_with(|| vec![]).push(*r);
            });

            let new_chains = filter_chain_group(chain_group);
            /*
            let new_chains = chain_group
                .values()
                .into_iter()
                .filter(|chain| chain.len() >= 0)
                .map(|chain| chain.clone())
                .flatten()
                .collect::<Vec<[u32; 8]>>();
                */
            (sid, seqlen, new_chains)
        })
        .filter(|(_sid, _seqlen, new_chains)| new_chains.len() > 0)
        .map(|(sid, seqlen, new_chains)| {
            let shmmrmap = build_shmmer_map_from_query_results(&new_chains);
            let mq = map_interval_query(&shmmrmap, sid as u32, 0, seqlen as u32);

            //TODO: we need to build a new chain group and for overlapped chains, we only use the longest one

            let v = mq
                .into_iter()
                .map(|m| {
                    let dcount = generate_deltas(&sdb0, &sdb1, m, k);
                    let mut out = [0_u32; 8];
                    out[0..7].copy_from_slice(&m[0..7]);
                    out[7] = dcount as u32;
                    out
                })
                .collect::<Vec<[u32; 8]>>();

            linearize_hits(v)
        })
        .flatten()
        .collect::<Vec<[u32; 8]>>()
        .into_iter()
        .try_for_each(|v| -> Result<(), std::io::Error> {
            writeln!(
                out_file,
                "{} {} {} {} {} {} {} {}",
                v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]
            )?;
            Ok(())
        })?;

    let mut idx0_file = BufWriter::new(File::create(format!("{}.0.path", output_prefix))?);
    let mut idx1_file = BufWriter::new(File::create(format!("{}.1.path", output_prefix))?);
    sdb0.shmmrs.into_iter().flatten().try_for_each(|shmmr| {
        let rid = shmmr.y >> 32;
        let pos = (shmmr.y & 0xFFFFFFFF) >> 1;
        let strand = shmmr.y & 0x1;
        writeln!(idx0_file, "{} {} {}", rid, pos, strand)
    })?;

    sdb1.shmmrs.into_iter().flatten().try_for_each(|shmmr| {
        let rid = shmmr.y >> 32;
        let pos = (shmmr.y & 0xFFFFFFFF) >> 1;
        let strand = shmmr.y & 0x1;
        writeln!(idx1_file, "{} {} {}", rid, pos, strand)
    })?;
    Ok(())
}
