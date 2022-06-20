use libwfa::{affine_wavefront::*, bindings::*, mm_allocator::*, penalties::*};
use regex::Regex;
use std::result::Result;

#[derive(Debug, Clone)]
pub struct SeqLocus {
    pub id: u32,
    pub bgn: u32,
    pub len: u32,
}

#[derive(Debug, Clone, PartialEq)]
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

pub fn get_cigar(seq0: &String, seq1: &String) -> Result<(isize, String, Vec<u8>), &'static str> {
    let alloc = MMAllocator::new(BUFFER_SIZE_8M as u64);
    let pattern = seq0.to_string();
    let text = seq1.to_string();
    let mut penalties = AffinePenalties {
        match_: 0,
        mismatch: 4,
        gap_opening: 6,
        gap_extension: 1,
    };
    let pat_len = pattern.as_bytes().len();
    let text_len = text.as_bytes().len();

    let mut wavefronts =
        AffineWavefronts::new_reduced(pat_len, text_len, &mut penalties, 100, 100, &alloc);
    let m = wavefronts.align(pattern.as_bytes(), text.as_bytes());
    match m {
        Ok(()) => {
            let score = wavefronts.edit_cigar_score(&mut penalties);
            //let cigar = wavefronts.cigar_bytes_raw();
            let cigar = wavefronts.cigar_bytes();
            let cg_str = std::str::from_utf8(&cigar).unwrap();

            Ok((score, cg_str.to_string(), wavefronts.cigar_bytes_raw()))
        }
        Err(_) => Err("align failed"),
    }
}

pub fn get_aln_segements(
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
    let unit = Regex::new(r"([0-9]+[MIDX])").unwrap();
    unit.captures_iter(cigar.1.as_str())
        .map(|m| {
            let u = &m[0];
            let tag = u.as_bytes().get(u.len() - 1);
            let adv = u[0..u.len() - 1].parse::<u32>().unwrap();
            let advs = match tag {
                Some(b'M') => Some((AlnSegType::Match, adv, adv)),
                Some(b'X') => Some((AlnSegType::Mismatch, adv, adv)),
                Some(b'I') => Some((AlnSegType::Insertion, 0, adv)),
                Some(b'D') => Some((AlnSegType::Deletion, adv, 0)),
                _ => None,
            }
            .unwrap();
            advs
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

pub fn get_aln_map(aln_segs: &Vec<AlnSegment>, s0: &str, s1: &str) -> Result<AlnMap, &'static str> {
    let mut pmap = Vec::<(u32, u32)>::with_capacity(1 << 14);
    let mut ref_a_seq = Vec::<u8>::with_capacity(1 << 14);
    let mut tgt_a_seq = Vec::<u8>::with_capacity(1 << 14);
    let mut aln_seq = Vec::<u8>::with_capacity(1 << 14 );
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
    fn test_get_aln_map() {
        let seq0 = include_str!("../test/test_data/seq0")
            .trim_end()
            .to_string();
        let seq1 = include_str!("../test/test_data/seq1")
            .trim_end()
            .to_string();
        let v = super::get_aln_segements(0, &seq0, 1, &seq1).unwrap();
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
        let aln_segs = super::get_aln_segements(0, &seq0, 1, &seq1).unwrap();
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
