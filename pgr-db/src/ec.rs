#![warn(missing_docs)]
//! function for error correction

use crate::aln::query_fragment_to_hps;
use crate::fasta_io::reverse_complement;
use crate::graph_utils::{ShmmrGraphNode, WeightedNode};
use crate::seq_db;
use crate::shmmrutils::{sequence_to_shmmrs, ShmmrSpec};
use petgraph::algo::toposort;
use petgraph::EdgeDirection::Outgoing;
use petgraph::{graphmap::DiGraphMap, EdgeDirection::Incoming};
use rustc_hash::{FxHashMap, FxHashSet};

/// perform error correction using de Bruijn graph
/// just a naive approach for now
/// each input sequence is expected to be starting and ending at the "same" position
///
/// this methods can ignore haplotype specific signals
///
pub fn naive_dbg_consensus(
    seqs: Vec<Vec<u8>>,
    kmer_size: usize,
    min_cov: usize,
) -> Result<Vec<u8>, &'static str> {
    let mut db_g = DiGraphMap::<usize, u32>::new();
    let mut kmer_idx = FxHashMap::<Vec<u8>, usize>::default();
    let mut idx_kmer = Vec::<Vec<u8>>::new();
    let mut kmer_count = FxHashMap::<usize, usize>::default();
    let mut kmer_max_idx = 0;

    let tgt_seq = seqs[0].clone();
    for seq in seqs.into_iter() {
        if seq.len() < kmer_size {
            panic!("sequence needs to be longer than the k-mer size");
        }
        let kmer0 = seq[0..kmer_size].to_vec();
        let mut kidx0 = *kmer_idx.entry(kmer0.clone()).or_insert_with(|| {
            let m = kmer_max_idx;
            idx_kmer.push(kmer0);
            kmer_max_idx += 1;
            m
        });
        *kmer_count.entry(kidx0).or_insert(0) += 1;
        let mut kidx1 = 0;
        (1..seq.len() - kmer_size + 1).into_iter().for_each(|p| {
            let kmer1 = seq[p..p + kmer_size].to_vec();
            kidx1 = *kmer_idx.entry(kmer1.clone()).or_insert_with(|| {
                let m = kmer_max_idx;
                idx_kmer.push(kmer1);
                kmer_max_idx += 1;
                m
            });
            *kmer_count.entry(kidx1).or_insert(0) += 1;
            db_g.add_edge(kidx0, kidx1, 1);
            kidx0 = kidx1;
        });
    }

    let get_best_path = |kmers: Vec<usize>| -> Vec<u8> {
        let mut best_score = 0;
        let mut best_node = 0;

        let mut node_score = FxHashMap::<usize, u64>::default();
        let mut track_back = FxHashMap::<usize, Option<usize>>::default();

        kmers.into_iter().for_each(|m| {
            let in_edges = db_g.edges_directed(m, Incoming);
            let mut bs = 0;
            let mut bn: Option<usize> = None;
            let ms = *kmer_count.get(&m).unwrap();
            in_edges.into_iter().for_each(|(v, _w, _)| {
                if bn.is_none() {
                    bs = *node_score.get(&v).unwrap();
                    bn = Some(v);
                } else {
                    let s = *node_score.get(&v).unwrap();
                    if s > bs {
                        bs = s;
                        bn = Some(v);
                    }
                }
            });
            let ns = bs + ms as u64;
            node_score.insert(m, ns);
            track_back.insert(m, bn);

            if ns > best_score {
                best_score = ns;
                best_node = m;
            }
        });

        let mut tgt_rev_path = FxHashMap::<usize, Option<usize>>::default();
        (0..tgt_seq.len() - kmer_size + 1)
            .into_iter()
            .for_each(|p| {
                if p != 0 {
                    let kmer0 = tgt_seq[p..p + kmer_size].to_vec();
                    let idx0 = *kmer_idx.get(&kmer0).unwrap();
                    let kmer1 = tgt_seq[p - 1..p + kmer_size - 1].to_vec();
                    let idx1 = *kmer_idx.get(&kmer1).unwrap();
                    // println!("{:?} {:?} {} {}", kmer0, kmer1, idx0, idx1);
                    tgt_rev_path.insert(idx0, Some(idx1));
                } else {
                    let kmer0 = tgt_seq[p..p + kmer_size].to_vec();
                    let idx0 = *kmer_idx.get(&kmer0).unwrap();
                    tgt_rev_path.insert(idx0, None);
                }
            });

        let last_kmer = tgt_seq[tgt_seq.len() - kmer_size..tgt_seq.len()].to_vec();
        // println!("{:?}", last_kmer);
        let last_tgt_idx = *kmer_idx.get(&last_kmer.to_vec()).unwrap();
        let mut rev_path = Vec::<usize>::new();
        let mut cur_idx = last_tgt_idx;
        rev_path.push(cur_idx);
        loop {
            if let Some(p_node) = tgt_rev_path.get(&cur_idx) {
                if let Some(p_idx) = p_node {
                    if *kmer_count.get(&p_idx).unwrap() >= min_cov {
                        cur_idx = *p_idx;
                        rev_path.push(cur_idx);
                        continue;
                    }
                }
            }

            if let Some(p_node) = track_back.get(&cur_idx) {
                if let Some(p_idx) = p_node {
                    cur_idx = *p_idx;
                    rev_path.push(cur_idx);
                } else {
                    break;
                }
            } else {
                break;
            }
        }

        rev_path.reverse();
        let path = rev_path;

        let mut bases = Vec::<u8>::new();
        bases.extend(idx_kmer[path[0]].iter());
        path[1..].iter().for_each(|&p| {
            bases.push(idx_kmer[p][kmer_size - 1]);
        });
        bases
    };

    match toposort(&db_g, None) {
        Ok(kmers) => Ok(get_best_path(kmers)),
        Err(_) => Err("circle found"),
    }
}

/// perform error correction using shimmer de Bruijn graph
///
/// this methods can ignore haplotype specific signals
///
pub fn shmmr_dbg_consensus(
    seqs: Vec<Vec<u8>>,
    shmmr_spec: &Option<ShmmrSpec>,
) -> Result<Vec<(Vec<u8>, Vec<u32>)>, &'static str> {
    let shmmr_spec = shmmr_spec.as_ref().unwrap_or(&ShmmrSpec {
        w: 31,
        k: 31,
        r: 1,
        min_span: 0,
        sketch: false,
    });
    assert!(shmmr_spec.k % 2 == 1); // the k needs to odd to break symmetry
    assert!(shmmr_spec.min_span == 0); // if min_span != 0, we don't get consistent path
    let mut sdb = seq_db::CompactSeqDB::new(shmmr_spec.clone());
    let seqs = (0..seqs.len())
        .into_iter()
        .map(|sid| {
            (
                sid as u32,
                Some("Memory".to_string()),
                format!("{}", sid),
                seqs[sid].clone(),
            )
        })
        .collect::<Vec<(u32, Option<String>, String, Vec<u8>)>>();
    sdb.load_index_from_seq_vec(&seqs);

    let mut frg_seqs = FxHashMap::<ShmmrGraphNode, Vec<u8>>::default();

    let mut score = FxHashMap::<ShmmrGraphNode, u32>::default();
    sdb.frag_map.iter().for_each(|(k, v)| {
        let (_, sid, b, e, strand) = v[0];
        let b = (b - shmmr_spec.k) as usize;
        let e = e as usize;
        let seq = seqs[sid as usize].3[b..e].to_vec();
        let node = ShmmrGraphNode(k.0, k.1, strand);
        score.insert(node, v.len() as u32);
        if !frg_seqs.contains_key(&node) {
            frg_seqs.insert(node, seq.clone());
        };
        let seq = reverse_complement(&seq);
        let node = ShmmrGraphNode(k.0, k.1, 1 - strand);
        score.insert(node, v.len() as u32);
        if !frg_seqs.contains_key(&node) {
            frg_seqs.insert(node, seq.clone());
        };
    });

    let adj_list = sdb.generate_smp_adj_list(0);
    let s0 = adj_list[0];

    //println!("s0:{:?}", s0);

    let start = ShmmrGraphNode(s0.1 .0, s0.1 .1, s0.1 .2);

    use crate::graph_utils::BiDiGraphWeightedDfs;

    let mut g = DiGraphMap::<ShmmrGraphNode, ()>::new();
    adj_list.into_iter().for_each(|(_sid, v, w)| {
        let v = ShmmrGraphNode(v.0, v.1, v.2);
        let w = ShmmrGraphNode(w.0, w.1, w.2);
        g.add_edge(v, w, ());

        //println!("DBG: add_edge {:?} {:?}", v, w);
        //*score.entry(v).or_insert(0) += 1;
        //*score.entry(w).or_insert(0) += 1;
    });

    //println!("DBG: node_count {:?} {:?}", g.node_count(), g.edge_count());
    //println!("DBG: {} {}", g.node_count(), g.edge_count());

    let mut wdfs_walker = BiDiGraphWeightedDfs::new(&g, start, &score);
    let mut out = vec![];
    loop {
        if let Some((node, p_node, is_leaf, rank, branch_id, branch_rank)) = wdfs_walker.next(&g) {
            let node_count = *score.get(&node).unwrap();
            let p_node = match p_node {
                Some(pnode) => Some((pnode.0, pnode.1, pnode.2)),
                None => None,
            };
            out.push((
                node,
                p_node,
                node_count,
                is_leaf,
                rank,
                branch_id,
                branch_rank,
            ));
        } else {
            break;
        }
    }

    let mut out_seqs = Vec::<(Vec<u8>, Vec<u32>)>::new();

    let mut out_seq = vec![];
    let mut out_cov = vec![];
    //let mut head_orientation = 0_u8;
    for (node, _p_node, node_count, is_leaf, _rank, _branch_id, _branch_rank) in out {
        if out_seq.len() == 0 {
            let seq = frg_seqs.get(&node).unwrap().clone();
            for _ in 0..seq.len() {
                out_cov.push(node_count);
            }
            out_seq.extend(seq);
        } else {
            let k = shmmr_spec.k as usize;
            let seq = frg_seqs.get(&node).unwrap().clone();
            assert!(out_seq[out_seq.len() - k..] == seq[..k]);
            let seq = seq[k..].to_vec();
            for _ in 0..seq.len() {
                out_cov.push(node_count);
            }
            out_seq.extend(seq);
        }
        if is_leaf {
            out_seqs.push((out_seq.clone(), out_cov.clone()));
            out_seq.clear();
            out_cov.clear()
        }
    }
    Ok(out_seqs)
}

/// perform error correction using shimmer de Bruijn graph
///
/// this methods try to perseve SNP specific to a guide read (the first one in the list)
/// if there is more or equal to the "min_cov"
///
pub fn guided_shmmr_dbg_consensus(
    seqs: Vec<Vec<u8>>,
    shmmr_spec: &Option<ShmmrSpec>,
    min_cov: u32,
) -> Result<(Vec<u8>, Vec<u32>), &'static str> {
    let shmmr_spec = shmmr_spec.as_ref().unwrap_or(&ShmmrSpec {
        w: 31,
        k: 31,
        r: 1,
        min_span: 0,
        sketch: false,
    });
    assert!(shmmr_spec.k % 2 == 1); // the k needs to odd to break symmetry
    assert!(shmmr_spec.min_span == 0); // if min_span != 0, we don't get consistent path
    let mut sdb = seq_db::CompactSeqDB::new(shmmr_spec.clone());
    let seqs = (0..seqs.len())
        .into_iter()
        .map(|sid| {
            (
                sid as u32,
                Some("Memory".to_string()),
                format!("{}", sid),
                seqs[sid].clone(),
            )
        })
        .collect::<Vec<(u32, Option<String>, String, Vec<u8>)>>();
    sdb.load_index_from_seq_vec(&seqs);

    let mut frg_seqs = FxHashMap::<ShmmrGraphNode, Vec<u8>>::default();
    let mut score = FxHashMap::<ShmmrGraphNode, u32>::default();
    sdb.frag_map.iter().for_each(|(k, v)| {
        let (_, sid, b, e, strand) = v[0];
        let b = (b - shmmr_spec.k) as usize;
        let e = e as usize;
        let seq = seqs[sid as usize].3[b..e].to_vec();
        let node = ShmmrGraphNode(k.0, k.1, strand);
        score.insert(node, v.len() as u32);
        if !frg_seqs.contains_key(&node) {
            frg_seqs.insert(node, seq.clone());
        };
        let seq = reverse_complement(&seq);
        let node = ShmmrGraphNode(k.0, k.1, 1 - strand);
        score.insert(node, v.len() as u32);
        if !frg_seqs.contains_key(&node) {
            frg_seqs.insert(node, seq.clone());
        };
    });

    let adj_list = sdb.generate_smp_adj_list(0);
    let s0 = adj_list[0];

    //println!("s0:{:?}", s0);

    let mut g = DiGraphMap::<ShmmrGraphNode, ()>::new();
    adj_list.into_iter().for_each(|(_sid, v, w)| {
        let v = ShmmrGraphNode(v.0, v.1, v.2);
        let w = ShmmrGraphNode(w.0, w.1, w.2);
        g.add_edge(v, w, ());

        //println!("DBG: add_edge {:?} {:?}", v, w);
    });

    let get_shmmr_nodes_from_seq = |seq: &Vec<u8>| -> Vec<((u64, u64, u8), u32)> {
        let shmmrs = sequence_to_shmmrs(0, &seq, &shmmr_spec, false);
        seq_db::pair_shmmrs(&shmmrs)
            .iter()
            .map(|(s0, s1)| {
                let p0 = s0.pos() + 1;
                //let p1 = s1.pos() + 1;
                let s0 = s0.hash();
                let s1 = s1.hash();
                if s0 < s1 {
                    ((s0, s1, 0_u8), p0)
                } else {
                    ((s1, s0, 1_u8), p0)
                }
            })
            .collect::<Vec<((u64, u64, u8), u32)>>()
    };

    let mut guide_nodes = FxHashMap::<ShmmrGraphNode, u32>::default();

    get_shmmr_nodes_from_seq(&seqs[0].3)
        .into_iter()
        .for_each(|(n, p)| {
            let node = ShmmrGraphNode(n.0, n.1, n.2);
            let s = *score.get(&node).unwrap();
            if s >= min_cov {
                guide_nodes.insert(node, p);
            }
        });
    //println!("DBG: node_count {:?} {:?}", g.node_count(), g.edge_count());
    //println!("DBG: {} {}", g.node_count(), g.edge_count());

    //let mut wdfs_walker = BiDiGraphWeightedDfs::new(&g, start, &score);

    let start = ShmmrGraphNode(s0.1 .0, s0.1 .1, s0.1 .2);
    let w = *score.get(&start).unwrap();

    let mut out = vec![];
    let mut next_node: WeightedNode<ShmmrGraphNode> = WeightedNode(w, start);
    let mut visited = FxHashSet::<ShmmrGraphNode>::default();
    let mut last_in_guide_nodes: Option<ShmmrGraphNode> = None;
    loop {
        let node = next_node;
        if visited.insert(node.1) {
            // the node is not visited before if insert() return true
            let mut out_count = 0_usize;
            let mut succ_list_f = Vec::<WeightedNode<ShmmrGraphNode>>::new();
            let mut next_guide_node: Option<WeightedNode<ShmmrGraphNode>> = None;
            let mut min_dist: Option<u32> = None;
            let current_node_position = guide_nodes.get(&node.1);
            for succ in g.neighbors_directed(node.1, Outgoing) {
                //println!("DBG: succ: {:?} {:?}", node.1, succ);
                if !visited.contains(&succ) {
                    //println!("DBG: pushing0: {:?}", succ);
                    out_count += 1;
                    let s = *score.get(&succ).unwrap();
                    if guide_nodes.contains_key(&succ) {
                        // choose the closest one for repetitive cases
                        if let Some(pos) = current_node_position {
                            let pos2 = *guide_nodes.get(&succ).unwrap();
                            if pos2 > *pos {
                                if min_dist.is_none() {
                                    let dist = pos2 - *pos;
                                    min_dist = Some(dist);
                                    next_guide_node = Some(WeightedNode(s, succ));
                                } else {
                                    let dist = pos2 - *pos;
                                    if dist < min_dist.unwrap() {
                                        next_guide_node = Some(WeightedNode(s, succ));
                                    }
                                }
                            }
                        } else {
                            next_guide_node = Some(WeightedNode(s, succ));
                        }
                    } else {
                        succ_list_f.push(WeightedNode(s, succ));
                    }
                }
            }
            if out_count == 0 {
                break;
            } else {
                if next_guide_node.is_some() {
                    next_node = next_guide_node.unwrap();
                    last_in_guide_nodes = Some(next_node.1.clone());
                } else {
                    if succ_list_f.len() > 0 {
                        succ_list_f.sort();
                        next_node = succ_list_f.pop().unwrap();
                    } else {
                        break;
                    }
                }
            }
            out.push((node.1, *score.get(&node.1).unwrap()));
        }
    }

    let mut out_seq = vec![];
    let mut out_cov = vec![];
    //let mut head_orientation = 0_u8;
    for (node, node_count) in out {
        if out_seq.len() == 0 {
            let seq = frg_seqs.get(&node).unwrap().clone();
            for _ in 0..seq.len() {
                out_cov.push(node_count);
            }
            out_seq.extend(seq);
        } else {
            let k = shmmr_spec.k as usize;
            let seq = frg_seqs.get(&node).unwrap().clone();
            /*
            println!(
                "{} {} {} {} {}",
                node.0,
                node.1,
                node.2,
                String::from_utf8_lossy(&out_seq[out_seq.len() - k..]),
                String::from_utf8_lossy(&seq[..k])
            );
            */
            assert!(out_seq[out_seq.len() - k..] == seq[..k]);
            let seq = seq[k..].to_vec();
            for _ in 0..seq.len() {
                out_cov.push(node_count);
            }
            out_seq.extend(seq);
        }
        /*
        if guide_nodes.contains(&node) {
            println!("XX");
        } else {
            println!("YY");
        }
        */
        if last_in_guide_nodes.is_some() {
            if node == last_in_guide_nodes.unwrap() {
                break;
            }
        }
    }
    Ok((out_seq, out_cov))
}

/// perform error correction using shimmer alignment
///
/// this methods try to perseve SNP specific to a guide read (the first one in the list)
/// if there is more or equal to the "min_cov"
///
/// 
pub fn shmmr_sparse_aln_consensus(
    seqs: Vec<Vec<u8>>,
    shmmr_spec: &Option<ShmmrSpec>,
    min_cov: u32,
) -> Result<Vec<(Vec<u8>, Vec<u32>)>, &'static str> {
    let shmmr_spec = shmmr_spec.as_ref().unwrap_or(&ShmmrSpec {
        w: 33,
        k: 33,
        r: 1,
        min_span: 0,
        sketch: false,
    });
    assert!(shmmr_spec.k % 2 == 1); // the k needs to odd to break symmetry
    assert!(shmmr_spec.min_span == 0); // if min_span != 0, we don't get consistent path
    let mut sdb = seq_db::CompactSeqDB::new(shmmr_spec.clone());
    let seqs = (0..seqs.len())
        .into_iter()
        .map(|sid| {
            (
                sid as u32,
                Some("Memory".to_string()),
                format!("{}", sid),
                seqs[sid].clone(),
            )
        })
        .collect::<Vec<(u32, Option<String>, String, Vec<u8>)>>();
    sdb.load_index_from_seq_vec(&seqs);
    let hit_pairs = query_fragment_to_hps(
        &sdb.frag_map,
        &seqs[0].3,
        &shmmr_spec,
        0.1,
        Some(32),
        Some(32),
        Some(32),
        Some(33),
    );

    let mut hit_map = FxHashMap::<(u32, u32, u8), Vec<(u32, (u32, u32, u8))>>::default();
    hit_pairs.into_iter().for_each(|(sid, hits)| {
        if hits.len() > 0 {
            // only use the main chian
            let hps = &hits[0].1;
            hps.into_iter().for_each(|&(v, w)| {
                hit_map.entry(v).or_insert(vec![]).push((sid, w));
            })
        }
    });

    let mut keys = hit_map.keys().map(|v| *v).collect::<Vec<(u32, u32, u8)>>();
    let mut reliable_regions = Vec::<((u32, u32, u8), u32)>::new();
    keys.sort();
    keys.into_iter().for_each(|k| {
        let m = hit_map.get(&k).unwrap();
        //println!("{:?} {:?}",k ,m.len());
        if m.len() >= min_cov as usize {
            reliable_regions.push((k, m.len() as u32));
        };
        // m.into_iter().for_each(|v| {
        //     println!("M : {:?} {:?}", k, v);
        // })
    });

    let mut out_seqs = vec![];
    let mut seq = vec![];
    let mut cov = vec![];
    let mut p_region: Option<((u32, u32, u8), u32)> = None;
    reliable_regions.into_iter().for_each(|(r, c)| {
        if p_region.is_none() {
            p_region = Some((r, c));
            seq.extend(seqs[0].3[0..r.1 as usize].to_vec());
            (0..r.1).into_iter().for_each(|_| {
                cov.push(c);
            });
        } else {
            // println!("R : {:?} {:?}", r, p_region);
            if r.0 == p_region.unwrap().0 .1 {
                seq.extend(seqs[0].3[r.0 as usize..r.1 as usize].to_vec());
                (r.0..r.1).into_iter().for_each(|_| {
                    cov.push(c);
                });
            } else {
                // println!("X : {:?} {}", r, seq.len());
                let p_hit = hit_map.get(&p_region.unwrap().0).unwrap();
                let c_hit = hit_map.get(&r).unwrap();
                let p_hit = p_hit
                    .into_iter()
                    .map(|v| *v)
                    .collect::<FxHashMap<u32, (u32, u32, u8)>>();
                let c_hit = c_hit
                    .into_iter()
                    .map(|v| *v)
                    .collect::<FxHashMap<u32, (u32, u32, u8)>>();

                let mut s = vec![];
                //let mut s0 = vec![];
                let mut c2 = 0_u32;
                let k = shmmr_spec.k as usize;
                for (sid, v) in p_hit {
                    if sid == 0 {
                        //let w = *c_hit.get(&sid).unwrap();
                        //s0 = seqs[sid as usize].3[v.1 as usize..w.0 as usize].to_vec();
                        continue;
                    }
                    if c_hit.contains_key(&sid) {
                        let w = *c_hit.get(&sid).unwrap();
                        // println!("S: {} {:?} {:?}", sid, v, w);
                        if s.len() == 0 {
                            // patch in, TODO: we will need to do some sub-consensus
                            if v.0 < w.0 && v.1 < w.1 && v.1 < w.0 {
                                s = seqs[sid as usize].3[v.1 as usize..w.0 as usize].to_vec();
                            } else if w.0 < v.0 && w.1 < v.1 && w.1 < v.0 {
                                s = seqs[sid as usize].3[w.1 as usize - k..v.0 as usize - k]
                                    .to_vec();
                                s = reverse_complement(&s);
                            } else {
                                continue;
                            }
                        }
                        //println!("S2: {}", String::from_utf8_lossy(&s[..]));
                        c2 += 1;
                    }
                }

                if c2 >= min_cov {
                    (0..s.len()).into_iter().for_each(|_| cov.push(c2));
                    seq.extend(s);
                    seq.extend(seqs[0].3[r.0 as usize..r.1 as usize].to_vec());
                    (r.0..r.1).into_iter().for_each(|_| {
                        cov.push(c);
                    });
                } else {
                    out_seqs.push((seq.clone(), cov.clone()));
                    seq.clear();
                    cov.clear();
                    // seq.extend(s0);
                    seq.extend(seqs[0].3[r.0 as usize..r.1 as usize].to_vec());
                    (r.0..r.1).into_iter().for_each(|_| {
                        cov.push(c);
                    });
                }
            }
            p_region = Some((r, c));
        }
    });

    out_seqs.push((seq.clone(), cov.clone()));

    Ok(out_seqs)
}

#[cfg(test)]
mod test {
    use crate::ec::guided_shmmr_dbg_consensus;
    use crate::ec::naive_dbg_consensus;
    use crate::ec::shmmr_dbg_consensus;
    use crate::ec::shmmr_sparse_aln_consensus;
    use crate::seq_db::CompactSeqDB;
    use crate::shmmrutils::ShmmrSpec;
    #[test]
    fn test_naive_dbg_consensus() {
        let spec = ShmmrSpec {
            w: 24,
            k: 24,
            r: 12,
            min_span: 12,
            sketch: false,
        };
        let mut sdb = CompactSeqDB::new(spec);
        let _ = sdb.load_seqs_from_fastx("test/test_data/consensus_test.fa".to_string());
        let seqs = (0..sdb.seqs.len())
            .into_iter()
            .map(|sid| sdb.get_seq_by_id(sid as u32))
            .collect::<Vec<Vec<u8>>>();

        let r = naive_dbg_consensus(seqs, 48, 2).unwrap();
        println!("{}", String::from_utf8_lossy(&r[..]));
    }

    #[test]
    fn test_shmmr_dbg_consensus() {
        let spec = ShmmrSpec {
            w: 24,
            k: 24,
            r: 12,
            min_span: 12,
            sketch: false,
        };
        let mut sdb = CompactSeqDB::new(spec);
        let _ = sdb.load_seqs_from_fastx("test/test_data/consensus_test3.fa".to_string());
        let seqs = (0..sdb.seqs.len())
            .into_iter()
            .map(|sid| sdb.get_seq_by_id(sid as u32))
            .collect::<Vec<Vec<u8>>>();

        let r = shmmr_dbg_consensus(seqs, &None).unwrap();
        for (s, c) in r {
            println!("{}", String::from_utf8_lossy(&s[..]));
            println!("{:?}", c);
        }
    }

    #[test]
    fn test_guided_shmmr_dbg_consensus() {
        let spec = ShmmrSpec {
            w: 24,
            k: 24,
            r: 12,
            min_span: 12,
            sketch: false,
        };
        let mut sdb = CompactSeqDB::new(spec);
        let _ = sdb.load_seqs_from_fastx("test/test_data/consensus_test.fa".to_string());
        let seqs = (0..sdb.seqs.len())
            .into_iter()
            .map(|sid| sdb.get_seq_by_id(sid as u32))
            .collect::<Vec<Vec<u8>>>();

        let (s, c) = guided_shmmr_dbg_consensus(seqs, &None, 2).unwrap();
        println!("{}", String::from_utf8_lossy(&s[..]));
        println!("{:?}", c);
    }

    #[test]
    fn test_shmmr_sparse_aln_consensus() {
        let spec = ShmmrSpec {
            w: 24,
            k: 24,
            r: 12,
            min_span: 12,
            sketch: false,
        };
        let mut sdb = CompactSeqDB::new(spec);
        let _ = sdb.load_seqs_from_fastx("test/test_data/consensus_test5.fa".to_string());
        let seqs = (0..sdb.seqs.len())
            .into_iter()
            .map(|sid| sdb.get_seq_by_id(sid as u32))
            .collect::<Vec<Vec<u8>>>();

        let r = shmmr_sparse_aln_consensus(seqs, &None, 2).unwrap();
        for (s, c) in r {
            println!("{}", String::from_utf8_lossy(&s[..]));
            println!("{:?}", c);
        }
    }
}
