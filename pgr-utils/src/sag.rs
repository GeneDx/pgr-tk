#![allow(dead_code)]

//se core::ops::Range;
use petgraph::graphmap::DiGraphMap;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, BufReader}; //{BufWriter, Write};

pub type MapIntervalRecord = [u32; 8];
pub type IntervalPaths = FxHashMap<u32, Vec<(u32, u32)>>;
pub type IntervalAlnGraph = DiGraphMap<(u32, u32, u32, u32), u32>;

//use petgraph::unionfind::UnionFind;
//use rayon::prelude::*;

pub fn read_aln_records(
    filename: &String,
    remove_self: bool,
) -> Result<Vec<MapIntervalRecord>, std::io::Error> {
    let reader = BufReader::new(File::open(filename)?);
    let out = reader
        .lines()
        .map(|line| {
            line.unwrap()
                .split_whitespace()
                .map(|c| c.parse::<u32>().unwrap())
                .collect::<Vec<u32>>()
        })
        .into_iter()
        .filter(|x| !remove_self || x[0] != x[3])
        .map(|x| {
            let mut out = [0_u32; 8];
            out[..].copy_from_slice(&x[..]);
            out
        })
        .collect::<Vec<MapIntervalRecord>>();
    Ok(out)
}

pub fn read_path(filename: &String) -> Result<IntervalPaths, std::io::Error> {
    let reader = BufReader::new(File::open(filename)?);
    let mut path = FxHashMap::<u32, Vec<u32>>::default();
    reader
        .lines()
        .map(|line| {
            line.unwrap()
                .split_whitespace()
                .map(|c| c.parse::<u32>().unwrap())
                .collect::<Vec<u32>>()
        })
        .for_each(|r| {
            path.entry(r[0]).or_insert_with(|| vec![]).push(r[1]);
        });
    //let mut out = FxHashMap::<u32, Vec<(u32, u32)>>::default();

    let out = path
        .into_iter()
        .map(|(sid, path0)| {
            //let mut last = 0_u32;

            //filter small intervals
            let path = path0.into_iter().collect::<Vec<u32>>();

            // convert path to intervals
            let itvls = path[0..path.len() - 1]
                .iter()
                .zip(path[1..path.len()].iter())
                .map(|x| (*x.0, *x.1))
                .collect::<Vec<(u32, u32)>>();
            (sid, itvls)
        })
        .collect::<IntervalPaths>();

    Ok(out)
}

pub fn get_pairs<T: Copy>(p: &Vec<T>) -> Vec<(T, T)> {
    p[0..p.len() - 1]
        .iter()
        .zip(p[1..p.len()].iter())
        .map(|(t0, t1)| (*t0, *t1))
        .collect::<Vec<(T, T)>>()
}

pub fn get_aln_best(aln_recs: &Vec<[u32; 8]>) -> FxHashMap<u32, (u32, u32, u32, u32)> {
    let mut aln_count = FxHashMap::<u32, FxHashMap<(u32, u32), u32>>::default(); // id0 -> (id1, strand), count

    aln_recs.iter().for_each(|&v| {
        let mut trgn: [u32; 4] = Default::default();
        trgn[0..3].copy_from_slice(&v[0..3]);
        trgn[3] = 0;
        let mut qrgn: [u32; 4] = Default::default();
        qrgn.copy_from_slice(&v[3..7]);
        let k = qrgn[0];
        let v = (trgn[0], qrgn[3]);

        let m = aln_count
            .entry(k)
            .or_insert_with(|| FxHashMap::<(u32, u32), u32>::default());
        *m.entry(v).or_insert(0) += 1;
    });

    let mut aln_best = FxHashMap::<u32, (u32, u32, u32, u32)>::default();

    aln_count.iter().for_each(|(s, m)| {
        let mut c = m
            .iter()
            .map(|(k, v)| (*v, *k))
            .collect::<Vec<(u32, (u32, u32))>>();

        c.sort_by(|a, b| b.partial_cmp(a).unwrap());
        let total_count = c.iter().map(|x| x.0).sum();
        let (count, (s1, strand)) = c[0];
        aln_best.insert(*s, (s1, strand, count, total_count));
    });
    aln_best
}

pub fn get_sag_from_aln_records(
    path0: &IntervalPaths,
    aln_records: &Vec<MapIntervalRecord>,
    contig: Option<u32>,
) -> (
    Vec<u32>,
    IntervalAlnGraph,
    FxHashMap<u32, (u32, u32, u32, u32)>,
) {
    let mut sag = DiGraphMap::<(u32, u32, u32, u32), u32>::new();

    let aln_best = get_aln_best(&aln_records);

    let mut a_ctgs = FxHashSet::<u32>::default();

    aln_best.iter().for_each(|(s0, bm)| {
        // ad hoc rule to determined if a contig is associated
        if (bm.2 as f32) / (bm.3 as f32) > 0.5
            && path0.get(&bm.0).unwrap().len() > path0.get(&s0).unwrap().len()
        {
            a_ctgs.insert(*s0);
            //println!("BM {} {:?} {}",s0, bm, (bm.2 as f32) / (bm.3 as f32));
        }
    });

    let mut pctg_ids = Vec::<u32>::new();
    path0.into_iter().for_each(|(&s0, _)| {
        if !a_ctgs.contains(&s0) {
            pctg_ids.push(s0);
        }
    });

    path0
        .into_iter()
        .filter(|(&s0, _path0)| {
            if let Some(the_contig) = contig {
                s0 == the_contig
            } else {
                !a_ctgs.contains(&s0) //not associate contigs
            }
        })
        .for_each(|(s0, path0)| {
            get_pairs(path0).iter().for_each(|e| {
                sag.add_edge((*s0, e.0 .0, e.0 .1, 0), (*s0, e.1 .0, e.1 .1, 0), 0);
                assert!(e.0 .1 == e.1 .0); // make sure there is no gap between nodes
            });
        });

    let aln_map = aln_records
        .into_iter()
        .filter(|&&r| {
            if let Some(the_contig) = contig {
                r[0] == the_contig
            } else {
                true
            }
        })
        .filter(|&&r| aln_best.contains_key(&r[3]))
        .filter(|&&r| {
            let a = *aln_best.get(&r[3]).unwrap();
            a.0 == r[0] && a.1 == r[6] && !a_ctgs.contains(&r[0])
        })
        .filter(|&&r| r[7] == 0)
        .map(|&r| ((r[3], r[4], r[5], r[6]), (r[0], r[1], r[2], 0)))
        .collect::<FxHashMap<(u32, u32, u32, u32), (u32, u32, u32, u32)>>();

    path0
        .iter()
        .filter(|(&s1, _path1)| a_ctgs.contains(&s1))
        .map(|(&s1, path1)| {
            let tag: u32;
            if aln_best.contains_key(&s1) {
                let (s0, _strand, count, _total_count) = aln_best.get(&s1).unwrap();
                if path1.len() > path0.get(s0).unwrap().len() {
                    tag = 2_u32;
                } else if *count > 1 && {
                    if let Some(the_contig) = contig {
                        *s0 == the_contig
                    } else {
                        true
                    }
                } {
                    tag = 1_u32;
                } else {
                    tag = 0_u32;
                }
            } else {
                tag = 0_u32;
            }
            (s1, path1, tag)
        })
        .for_each(|(s1, path1, tag)| {
            get_pairs(path1).iter().for_each(|e| {
                let strand: u32;
                let n0;
                let n1;
                match tag {
                    2 => return,
                    0 => {
                        if contig.is_some() {
                            return;
                        } else {
                            strand = 0;
                            n0 = (s1, e.0 .0, e.0 .1, strand);
                            n1 = (s1, e.1 .0, e.1 .1, strand);
                            sag.add_edge(n0, n1, 0);
                        }
                    }
                    1 => {
                        let (_s0, strand_, _count, _total_count) = aln_best.get(&s1).unwrap();
                        strand = *strand_;
                        let k0 = (s1, e.0 .0, e.0 .1, strand);
                        let k1 = (s1, e.1 .0, e.1 .1, strand);
                        assert!(k0.2 == k1.1); // make sure there is no gap between nodes
                        n0 = *aln_map.get(&k0).unwrap_or(&k0);
                        n1 = *aln_map.get(&k1).unwrap_or(&k1);
                        if strand == 0 {
                            sag.add_edge(k0, k1, 0);
                            sag.add_edge(k0, n1, 0);
                            sag.add_edge(n0, k1, 0);
                            sag.add_edge(n0, n1, 0);
                        } else {
                            sag.add_edge(k1, k0, 0);
                            sag.add_edge(k1, n0, 0);
                            sag.add_edge(n1, k0, 0);
                            sag.add_edge(n1, n0, 0);
                        }
                        //println!("DEBUG {} {:?} {:?} {:?} {:?}", s1, k0, k1, n0, n1);
                    }
                    _ => {
                        strand = 0;
                        n0 = (s1, e.0 .0, e.0 .1, strand);
                        n1 = (s1, e.1 .0, e.1 .1, strand);
                        sag.add_edge(n0, n1, 0);
                    }
                }
            });
        });
    (pctg_ids, sag, aln_best)
}
