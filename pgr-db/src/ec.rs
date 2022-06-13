#![warn(missing_docs)]
//! function for error correction

use crate::fasta_io::reverse_complement;
use crate::graph_utils::ShmmrGraphNode;
use crate::seq_db::pair_shmmrs;
use crate::shmmrutils::{sequence_to_shmmrs, ShmmrSpec};
use petgraph::algo::toposort;
use petgraph::{graphmap::DiGraphMap, EdgeDirection::Incoming};
use rayon::prelude::*;
use rustc_hash::FxHashMap;

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
) -> Result<Vec<u8>, &'static str> {
    let shmmr_spec = shmmr_spec.as_ref().unwrap_or(&ShmmrSpec {
        w: 12,
        k: 32,
        r: 1,
        min_span: 0,
        sketch: false,
    });

    let pairs = seqs
        .par_iter()
        .enumerate()
        .flat_map(|(sid, seq)| {
            //println!("len {}", seq.len());
            let shmmrs = sequence_to_shmmrs(sid as u32, seq, shmmr_spec, true);

            let pairs = pair_shmmrs(&shmmrs)
                .iter()
                .map(|(&shmmr0, &shmmr1)| {
                    let s0 = shmmr0.x >> 8;
                    let s1 = shmmr1.x >> 8;
                    let p0 = ((shmmr0.y & 0xFFFF_FFFF) >> 1) as u32;
                    let p1 = ((shmmr1.y & 0xFFFF_FFFF) >> 1) as u32;
                    let shmmr_pair = if s0 <= s1 {
                        (s0, s1, 0_u8, sid as u32, p0, p1)
                    } else {
                        (s1, s0, 1_u8, sid as u32, p0, p1)
                    };
                    shmmr_pair
                })
                .collect::<Vec<(u64, u64, u8, u32, u32, u32)>>();
            pairs
        })
        .collect::<Vec<(u64, u64, u8, u32, u32, u32)>>();

    let mut frg_seqs = FxHashMap::<(u64, u64, u8), Vec<u8>>::default();
    pairs.iter().for_each(|&v| {
        let (s0, s1, strand, sid, p0, p1) = v;
        let b = (p0 - shmmr_spec.k) as usize;
        let e = (p1 + 1) as usize;
        let seq = seqs[sid as usize][b..e].to_vec();
        if !frg_seqs.contains_key(&(s0, s1, strand)) {
            frg_seqs.insert((s0, s1, strand), seq.clone());
        };
        let seq = reverse_complement(&seq);
        if !frg_seqs.contains_key(&(s0, s1, 1 - strand)) {
            frg_seqs.insert((s0, s1, 1 - strand), seq.clone());
        };
    });

    let adj_list = (0..pairs.len() - 1)
        .into_par_iter()
        .flat_map(|i| {
            let v = pairs[i];
            let w = pairs[i + 1];
            vec![
                (v.3 as u32, (v.0, v.1, v.2), (w.0, w.1, w.2)),
                (v.3 as u32, (w.0, w.1, 1 - w.2), (v.0, v.1, 1 - v.2)),
            ]
        })
        .collect::<Vec<(u32, (u64, u64, u8), (u64, u64, u8))>>();

    //adj_list.iter().for_each(|f|println!("{:?}", f));

    let s0 = adj_list[0];

    //println!("s0:{:?}", s0);

    let start = ShmmrGraphNode(s0.1 .0, s0.1 .1, s0.1 .2);

    use crate::graph_utils::BiDiGraphWeightedDfs;

    let mut g = DiGraphMap::<ShmmrGraphNode, ()>::new();
    let mut score = FxHashMap::<ShmmrGraphNode, u32>::default();
    adj_list.into_iter().for_each(|(_sid, v, w)| {
        let v = ShmmrGraphNode(v.0, v.1, v.2);
        let w = ShmmrGraphNode(w.0, w.1, w.2);
        g.add_edge(v, w, ());

        //println!("DBG: add_edge {:?} {:?}", v, w);
        *score.entry(v).or_insert(0) += 1;
        *score.entry(w).or_insert(0) += 1;
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
                (node.0, node.1, node.2),
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

    let mut out_seq = vec![];
    let n = out[0];
    let seq = frg_seqs.get(&n.0).unwrap().clone();

    out_seq.extend(seq[..1 + shmmr_spec.k as usize].to_vec());
    for (node, _p_node, _node_count, is_leaf, _rank, _branch_id, _branch_rank) in out {
        out_seq.extend(
            frg_seqs.get(&(node.0, node.1, node.2)).unwrap()[1 + shmmr_spec.k as usize..].to_vec(),
        );
        if is_leaf {
            break;
        }
    }

    Ok(out_seq)
}

#[cfg(test)]
mod test {
    use crate::ec::naive_dbg_consensus;
    use crate::ec::shmmr_dbg_consensus;
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
        let _ = sdb.load_seqs_from_fastx("test/test_data/consensus_test.fa".to_string());
        let seqs = (0..sdb.seqs.len())
            .into_iter()
            .map(|sid| sdb.get_seq_by_id(sid as u32))
            .collect::<Vec<Vec<u8>>>();

        let r = shmmr_dbg_consensus(seqs, &None).unwrap();
        println!("{}", String::from_utf8_lossy(&r[..]));
    }
}
