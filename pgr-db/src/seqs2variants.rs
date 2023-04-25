use rs_libwfa2::wfa;
use std::result::Result;

#[derive(Debug, Clone)]
pub struct SeqLocus {
    pub id: u32,
    pub bgn: u32,
    pub len: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AlnSegType {
    Match,
    Mismatch,
    Insertion,
    Deletion,
    Unspecified,
}

#[derive(Debug, Clone)]
pub struct AlnSegment {
    pub ref_loc: SeqLocus,
    pub tgt_loc: SeqLocus,
    pub t: AlnSegType,
}

#[derive(Debug)]
pub struct AlnMap {
    pub pmap: Vec<(u32, u32)>,
    pub ref_a_seq: Vec<u8>,
    pub tgt_a_seq: Vec<u8>,
    pub aln_seq: Vec<u8>,
}

pub type AlnSegments = Vec<AlnSegment>;

pub fn get_cigar(seq0: &String, seq1: &String) -> Result<(i32, Vec<u8>), &'static str> {
    let a = seq0;
    let b = seq1;
    unsafe {
        let mut attributes = wfa::wavefront_aligner_attr_default;
        // Do not use a heuristic (enabled by default).
        //attributes.heuristic.strategy = wfa::wf_heuristic_strategy_wf_heuristic_none;
        // Only compute the score (not a path).
        //attributes.alignment_scope = wfa::alignment_scope_t_compute_score;

        // Set the cost model and parameters.
        attributes.distance_metric = wfa::distance_metric_t_gap_affine;
        attributes.affine_penalties.mismatch = 4_i32;
        attributes.affine_penalties.gap_opening = 4_i32;
        attributes.affine_penalties.gap_extension = 1_i32;

        // Initialize the aligner object.
        // This should be reused for multiple queries.
        let wf_aligner = wfa::wavefront_aligner_new(&mut attributes);

        // Do the alignment.
        let status = wfa::wavefront_align(
            wf_aligner,
            a.as_ptr() as *const i8,
            a.len() as i32,
            b.as_ptr() as *const i8,
            b.len() as i32,
        );
        assert_eq!(status, 0);
        if status != 0 {
            return Err("wfa align failed");
        };

        let cigar = *(*wf_aligner).cigar;
        //println!("{:?}", cigar);
        let cigar_vec = (cigar.begin_offset..cigar.end_offset)
            .into_iter()
            .map(|offset| *cigar.operations.add(offset as usize) as u8)
            .collect::<Vec<u8>>();
        let score = cigar.score;
        // Clean up memory.
        wfa::wavefront_aligner_delete(wf_aligner);
        Ok((score, cigar_vec))
    }
}

pub fn get_aln_segments(
    ref_id: u32,
    ref_seq: &String,
    tgt_id: u32,
    tgt_seq: &String,
) -> Result<AlnSegments, &'static str> {
    let cigar = get_cigar(ref_seq, tgt_seq);
    if cigar.is_err() {
        return Err("Fail to find variants");
    };
    let cigar = cigar.unwrap();
    let mut v = AlnSegments::new();
    let mut p0: u32 = 0;
    let mut p1: u32 = 0;

    let mut cigar_seg: Vec<Vec<u8>> = vec![];

    cigar.1.iter().for_each(|&c| {
        if cigar_seg.is_empty() {
            cigar_seg.push(vec![c])
        } else {
            let l = cigar_seg.len();
            let last = &mut cigar_seg[l - 1];
            if last[last.len() - 1] == c {
                last.push(c)
            } else {
                cigar_seg.push(vec![c]);
            }
        }
    });

    cigar_seg
        .iter()
        .map(|v| {
            let tag = v[0];
            let adv = v.len() as u32;
            match tag {
                b'M' => Some((AlnSegType::Match, adv, adv)),
                b'X' => Some((AlnSegType::Mismatch, adv, adv)),
                b'I' => Some((AlnSegType::Insertion, 0, adv)),
                b'D' => Some((AlnSegType::Deletion, adv, 0)),
                _ => None,
            }
            .unwrap()
        })
        .for_each(|advs| {
            let sl_r = SeqLocus {
                id: ref_id,
                bgn: p0,
                len: advs.1,
            };
            let sl_t = SeqLocus {
                id: tgt_id,
                bgn: p1,
                len: advs.2,
            };
            let vr = AlnSegment {
                ref_loc: sl_r,
                tgt_loc: sl_t,
                t: advs.0,
            };
            v.push(vr);
            p0 += advs.1;
            p1 += advs.2;
        });
    Ok(v)
}

pub fn get_aln_map(aln_segs: &[AlnSegment], s0: &str, s1: &str) -> Result<AlnMap, &'static str> {
    let mut pmap = Vec::<(u32, u32)>::with_capacity(1 << 14);
    let mut ref_a_seq = Vec::<u8>::with_capacity(1 << 14);
    let mut tgt_a_seq = Vec::<u8>::with_capacity(1 << 14);
    let mut aln_seq = Vec::<u8>::with_capacity(1 << 14);
    let mut aln_p = 0;
    aln_segs.iter().for_each(|f| {
        let ref_bgn = f.ref_loc.bgn;
        let ref_len = f.ref_loc.len;
        let tgt_bgn = f.tgt_loc.bgn;
        let tgt_len = f.tgt_loc.len;
        match f.t {
            AlnSegType::Match => {
                for idx in 0..ref_len {
                    let p0 = ref_bgn + idx;
                    let p1 = tgt_bgn + idx;
                    assert!(
                        (p0 as usize) < s0.len(),
                        "s0: index out of buond, WFA bug:\n{}\n{}\n",
                        s0,
                        s1
                    );
                    assert!(
                        (p1 as usize) < s1.len(),
                        "s1: index out of buond, WFA bug:\n{}\n{}\n",
                        s0,
                        s1
                    );
                    pmap.push((p1, aln_p + idx));
                    let b0 = s0.as_bytes()[p0 as usize];
                    let b1 = s1.as_bytes()[p1 as usize];
                    ref_a_seq.push(b0);
                    tgt_a_seq.push(b1);
                    aln_seq.push(b'|');
                }
                aln_p += ref_len;
            }
            AlnSegType::Mismatch => {
                for idx in 0..ref_len {
                    let p0 = ref_bgn + idx;
                    let p1 = tgt_bgn + idx;
                    assert!(
                        (p0 as usize) < s0.len(),
                        "s0: index out of buond, WFA bug:\n{}\n{}\n",
                        s0,
                        s1
                    );
                    assert!(
                        (p1 as usize) < s1.len(),
                        "s1: index out of buond, WFA bug:\n{}\n{}\n",
                        s0,
                        s1
                    );
                    pmap.push((p1, aln_p + idx));
                    let b0 = s0.as_bytes()[p0 as usize];
                    let b1 = s1.as_bytes()[p1 as usize];
                    ref_a_seq.push(b0);
                    tgt_a_seq.push(b1);
                    aln_seq.push(b'.');
                }
                aln_p += ref_len;
            }
            AlnSegType::Insertion => {
                for idx in 0..tgt_len {
                    let p1 = tgt_bgn + idx;
                    let b1 = s1.as_bytes()[p1 as usize];
                    ref_a_seq.push(b'-');
                    tgt_a_seq.push(b1);
                    aln_seq.push(b' ');
                }
                aln_p += tgt_len;
            }
            AlnSegType::Deletion => {
                for idx in 0..ref_len {
                    pmap.push((tgt_bgn, aln_p + idx));
                    let p0 = ref_bgn + idx;
                    let b0 = s0.as_bytes()[p0 as usize];
                    ref_a_seq.push(b0);
                    tgt_a_seq.push(b'-');
                    aln_seq.push(b' ');
                }
                aln_p += ref_len;
            }
            AlnSegType::Unspecified => {}
        };
    });
    Ok(AlnMap {
        pmap,
        ref_a_seq,
        tgt_a_seq,
        aln_seq,
    })
}

pub fn get_aln_fragment(
    ref_loc: &SeqLocus,
    aln_map: &AlnMap,
    ref_len: u32,
) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    let ref_a_seq = &aln_map.ref_a_seq;
    let aln_seq = &aln_map.aln_seq;
    let tgt_a_seq = &aln_map.tgt_a_seq;
    let ref_bgn = if ref_loc.bgn > 5 { ref_loc.bgn - 5 } else { 0 };
    let ref_end = if ref_loc.bgn + ref_loc.len + 5 < ref_len as u32 {
        ref_loc.bgn + ref_loc.len + 5
    } else {
        ref_len - 1
    };
    let bgn = aln_map.pmap.get(ref_bgn as usize).unwrap().1 as usize;
    let end = aln_map.pmap.get(ref_end as usize).unwrap().1 as usize;
    (
        ref_a_seq[bgn..end].to_vec(),
        aln_seq[bgn..end].to_vec(),
        tgt_a_seq[bgn..end].to_vec(),
    )
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_get_cigar2() {
        let seq0 = include_str!("../test/test_data/seq0")
            .trim_end()
            .to_string();
        let seq1 = include_str!("../test/test_data/seq1")
            .trim_end()
            .to_string();
        let cigar = get_cigar(&seq0, &seq1);
        println!("{:?}", cigar);
    }

    #[test]
    fn test_get_aln_map() {
        let seq0 = include_str!("../test/test_data/seq0")
            .trim_end()
            .to_string();
        let seq1 = include_str!("../test/test_data/seq1")
            .trim_end()
            .to_string();
        let v = super::get_aln_segments(0, &seq0, 1, &seq1).unwrap();
        let aln_map = get_aln_map(&v, &seq0, &seq1).unwrap();
        let ref_a_seq = aln_map.ref_a_seq;
        let aln_seq = aln_map.aln_seq;
        let tgt_a_seq = aln_map.tgt_a_seq;
        println!("{}", String::from_utf8_lossy(&ref_a_seq));
        println!("{}", String::from_utf8_lossy(&aln_seq));
        println!("{}", String::from_utf8_lossy(&tgt_a_seq));
        assert!(ref_a_seq == b"TCCATTCCCACCAGCAGTGTGTGAAAGTCTGGTACTGGTTCAGCCTGCCGTACTTTAATGATTATTGGTGTCACTCTTTCAAGTAACTTGTTGGTAATA--------AGAAGTCAATTA");
        assert!(  aln_seq == b"|||||||||||||||||||||||||           ||||||||||||||||||||||||||||||||||.||||||||||||||||||||||||||||        ||||||||||||");
        assert!(tgt_a_seq == b"TCCATTCCCACCAGCAGTGTGTGAA-----------GGTTCAGCCTGCCGTACTTTAATGATTATTGGTGACACTCTTTCAAGTAACTTGTTGGTAATATTTATCTAAGAAGTCAATTA");
    }
    #[test]
    fn test_get_aln_fragment() {
        let seq0 = include_str!("../test/test_data/seq0")
            .trim_end()
            .to_string();
        let seq1 = include_str!("../test/test_data/seq1")
            .trim_end()
            .to_string();
        let aln_segs = super::get_aln_segments(0, &seq0, 1, &seq1).unwrap();
        let aln_map = get_aln_map(&aln_segs, &seq0, &seq1).unwrap();
        let aln_segs = aln_segs
            .iter()
            .filter(|&s| s.t != AlnSegType::Match)
            .map(|s| s.clone())
            .collect::<Vec<AlnSegment>>();

        let out = get_aln_fragment(&aln_segs[0].ref_loc, &aln_map, seq0.len() as u32);
        assert!(String::from_utf8_lossy(&out.0) == "GTGAAAGTCTGGTACTGGTTC".to_string());
        assert!(String::from_utf8_lossy(&out.1) == "|||||           |||||".to_string());
        assert!(String::from_utf8_lossy(&out.2) == "GTGAA-----------GGTTC".to_string());

        let out = get_aln_fragment(&aln_segs[1].ref_loc, &aln_map, seq0.len() as u32);
        assert!(String::from_utf8_lossy(&out.0) == "TGGTGTCACTC".to_string());
        assert!(String::from_utf8_lossy(&out.1) == "|||||.|||||".to_string());
        assert!(String::from_utf8_lossy(&out.2) == "TGGTGACACTC".to_string());

        let out = get_aln_fragment(&aln_segs[2].ref_loc, &aln_map, seq0.len() as u32);
        assert!(String::from_utf8_lossy(&out.0) == "TAATA--------AGAAG".to_string());
        assert!(String::from_utf8_lossy(&out.1) == "|||||        |||||".to_string());
        assert!(String::from_utf8_lossy(&out.2) == "TAATATTTATCTAAGAAG".to_string());

        let ref_loc = SeqLocus {
            id: 0,
            bgn: 80,
            len: 1,
        };

        let out = get_aln_fragment(&ref_loc, &aln_map, seq0.len() as u32);

        assert!(String::from_utf8_lossy(&out.0) == "CTTTCAAGTAA".to_string());
        assert!(String::from_utf8_lossy(&out.1) == "|||||||||||".to_string());
        assert!(String::from_utf8_lossy(&out.2) == "CTTTCAAGTAA".to_string());
    }
}
