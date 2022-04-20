use crate::fasta_io::{reverse_complement, FastaReader};
use crate::shmmrutils::{match_reads, sequence_to_shmmrs, DeltaPoint, MM128};
use flate2::bufread::MultiGzDecoder;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

pub const KMERSIZE: u32 = 56;

pub type Bases = Vec<u8>;
pub type AlnSegments = (u32, bool, Vec<AlnSegment>);

#[derive(Debug, Clone)]
pub enum AlnSegment {
    // this still use a lot of space, we will find way to reduce the memory footprint later
    FullMatch,
    Match(u32, u32),
    Insertion(u8),
}

#[derive(Clone)]
pub enum Fragment {
    AlnSegments(AlnSegments),
    Prefix(Bases),
    Internal(Bases),
    Suffix(Bases),
}

impl<'a> fmt::Display for Fragment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Fragment::AlnSegments(d) => write!(f, "Seq:{} AlignSegs:{:?}", d.0, d.1),
            Fragment::Prefix(b) => write!(f, "Seq:{} AlignSegs:None", String::from_utf8_lossy(b)),
            Fragment::Internal(b) => write!(f, "Seq:{} AlignSegs:None", String::from_utf8_lossy(b)),
            Fragment::Suffix(b) => write!(f, "Seq:{} AlignSegs:None", String::from_utf8_lossy(b)),
        }
    }
}

pub type ShmmrPair = u128;
pub type Fragments = Vec<Fragment>;
pub type ShmmrToFrags = FxHashMap<ShmmrPair, Vec<(u32, u32)>>;

pub struct CompressedSeq {
    pub name: String,
    pub id: u32,
    pub shmmrs: Vec<MM128>,
    pub seq_frags: Vec<u32>,
}

pub struct CompressedSeqDB {
    pub filepath: String,
    pub seqs: Vec<CompressedSeq>,
    pub frag_map: ShmmrToFrags,
    pub frags: Fragments,
}

pub fn deltas_to_aln_segs(
    deltas: &Vec<DeltaPoint>,
    base_frg: &Vec<u8>,
    frg: &Vec<u8>,
) -> Vec<AlnSegment> {
    let mut aln_segs = Vec::<AlnSegment>::new();
    if deltas.len() == 0 {
        aln_segs.push(AlnSegment::FullMatch);
        //println!("aln_segs: {:?}", aln_segs);
        return aln_segs;
    }

    let mut x = base_frg.len();
    let mut y;
    // note: x - y = k
    for d in deltas.iter() {
        let x1 = d.x as usize;
        let y1 = d.y as usize;

        if x1 < x {
            aln_segs.push(AlnSegment::Match(x1 as u32, x as u32));
        }
        x = x1;
        y = y1;
        if d.dk > 0 {
            x -= d.dk as usize; // deletion from the base_frg
        } else {
            for yy in 0..(-d.dk) as usize {
                aln_segs.push(AlnSegment::Insertion(frg[y - yy - 1]));
            }
        }
    }
    if x != 0 {
        aln_segs.push(AlnSegment::Match(0, x as u32));
    };
    aln_segs.reverse();
    //println!("aln_segs: {:?}", aln_segs);
    aln_segs
}

pub fn reconstruct_seq_from_aln_segs(base_seq: &Vec<u8>, aln_segs: &Vec<AlnSegment>) -> Vec<u8> {
    let mut seq = Vec::<u8>::new();
    for s in aln_segs.iter() {
        match s {
            AlnSegment::FullMatch => {
                seq.extend_from_slice(&base_seq[..]);
            }
            AlnSegment::Match(x1, x2) => {
                seq.extend_from_slice(&base_seq[*x1 as usize..*x2 as usize]);
            }
            AlnSegment::Insertion(c) => {
                seq.push(*c);
            }
        }
    }
    seq
}

impl CompressedSeqDB {
    pub fn new(filepath: String) -> Self {
        let seqs = Vec::<CompressedSeq>::new();
        let frag_map = ShmmrToFrags::default();
        let frags = Vec::<Fragment>::new();
        CompressedSeqDB {
            filepath,
            seqs,
            frag_map,
            frags,
        }
    }

    pub fn seq_to_compressed(
        &mut self,
        name: String,
        id: u32,
        seq: &Vec<u8>,
        shmmrs: Vec<MM128>,
        try_compress: bool,
    ) -> CompressedSeq {
        //let shmmrs = sequence_to_shmmrs(id, &seq, 80, KMERSIZE, 4);
        let mut pos = 0;
        let mut seq_frags = Vec::<u32>::new();
        let mut frg_id = self.frags.len() as u32;
        let mut px: u128 = 0;

        for shmmr in shmmrs.iter() { // TODO: parallelize this
            let next_pos = shmmr.pos() + 1;
            if pos == 0 {
                let frg = seq[pos as usize..next_pos as usize].to_vec();
                self.frags.push(Fragment::Prefix(frg));
                seq_frags.push(frg_id);
                frg_id += 1;
                px = (shmmr.x >> 8) as u128;
                pos = next_pos;
                continue;
            }
            let mut aligned = false;
            let shmmr_pair = px << 64 | (shmmr.x >> 8) as u128;
            //println!("shmmr_pair: {} {} {:?}",px,  shmmr.x >> 8, shmmr_pair);

            if try_compress && self.frag_map.contains_key(&shmmr_pair) {
                let e = self.frag_map.get_mut(&shmmr_pair).unwrap();
                for t_frg_id in e.iter() {
                    let base_frg = self.frags.get(t_frg_id.0 as usize).unwrap();
                    if let Fragment::Internal(b) = base_frg {
                        let base_frg = b;
                        let frg = seq[(pos - KMERSIZE) as usize..next_pos as usize].to_vec();
                        let m = match_reads(base_frg, &frg, true, 0.1, 0, 0, 32);
                        if let Some(m) = m {
                            let deltas: Vec<DeltaPoint> = m.deltas.unwrap();
                            let aln_segs = deltas_to_aln_segs(&deltas, base_frg, &frg);
                            self.frags
                                .push(Fragment::AlnSegments((t_frg_id.0, false, aln_segs))); // false for the original strand
                            seq_frags.push(frg_id);
                            e.push((frg_id, id));
                            frg_id += 1;
                            aligned = true;
                            break; // we aligned to the first one of the fragments
                        } else {
                            continue;
                        }
                    }
                }
            };

            if try_compress && !aligned {
                // try reverse complement
                let shmmr_pair = ((shmmr.x >> 8) as u128) << 64 | px;
                if self.frag_map.contains_key(&shmmr_pair) {
                    let e = self.frag_map.get_mut(&shmmr_pair).unwrap();
                    for t_frg_id in e.iter() {
                        let base_frg = self.frags.get(t_frg_id.0 as usize).unwrap();
                        if let Fragment::Internal(b) = base_frg {
                            let base_frg = &reverse_complement(&b);
                            let frg = seq[(pos - KMERSIZE) as usize..next_pos as usize].to_vec();
                            let m = match_reads(base_frg, &frg, true, 0.1, 0, 0, 32);
                            if let Some(m) = m {
                                let deltas: Vec<DeltaPoint> = m.deltas.unwrap();
                                let aln_segs = deltas_to_aln_segs(&deltas, base_frg, &frg);
                                self.frags
                                    .push(Fragment::AlnSegments((t_frg_id.0, true, aln_segs))); // true for reverse complement
                                seq_frags.push(frg_id);
                                e.push((frg_id, id));
                                frg_id += 1;
                                aligned = true;
                                break; // we aligned to the first one of the fragments
                            } else {
                                continue;
                            }
                        }
                    }
                }
            };

            if !aligned || !try_compress {
                let frg = seq[(pos - KMERSIZE) as usize..next_pos as usize].to_vec();
                self.frags.push(Fragment::Internal(frg));
                seq_frags.push(frg_id);
                let shmmr_pair = px << 64 | (shmmr.x >> 8) as u128;
                let e = self
                    .frag_map
                    .entry(shmmr_pair)
                    .or_insert(Vec::<(u32, u32)>::new());
                e.push((frg_id, id));
                frg_id += 1;
            };

            px = (shmmr.x >> 8) as u128;
            pos = next_pos;
        }

        let frg = seq[pos as usize..].to_vec();
        self.frags.push(Fragment::Suffix(frg));
        seq_frags.push(frg_id);

        CompressedSeq {
            name,
            id,
            shmmrs,
            seq_frags,
        }
    }

    pub fn load_seqs(&mut self) -> Result<(), std::io::Error> {
        let file = File::open(&self.filepath)?;
        let mut reader = BufReader::new(file);
        let mut is_gzfile = false;
        {
            let r = reader.by_ref();
            let mut buf = Vec::<u8>::new();
            let _ = r.take(2).read_to_end(&mut buf);
            if buf == [0x1F_u8, 0x8B_u8] {
                log::info!(
                    "input file: {} detected as gz-compressed file",
                    self.filepath
                );
                is_gzfile = true;
            }
        }
        drop(reader);

        let file = File::open(&self.filepath)?;
        let mut reader = BufReader::new(file);
        let gz_buf = &mut BufReader::new(MultiGzDecoder::new(&mut reader));

        let file = File::open(&self.filepath)?;
        let reader = BufReader::new(file);
        let std_buf = &mut BufReader::new(reader);

        let fastx_buf: &mut dyn BufRead = if is_gzfile {
            drop(std_buf);
            gz_buf
        } else {
            drop(gz_buf);
            std_buf
        };

        let mut fastx_reader = FastaReader::new(fastx_buf, &self.filepath)?;
        let mut sid = 0;

        let mut seqs: Vec<(u32, String, Vec<u8>)> = Vec::new();
        while let Some(rec) = fastx_reader.next_rec() {
            let rec = rec.unwrap();
            let seqname = String::from_utf8_lossy(&rec.id).into_owned();
            seqs.push((sid, seqname, rec.seq));
            sid += 1;
        }

        let all_shmmers = seqs
            .par_iter()
            .map(|(sid, seqname, seq)| {
                let shmmrs = sequence_to_shmmrs(*sid, &seq, 80, KMERSIZE, 4);
                (*sid, shmmrs)
            })
            .collect::<Vec<(u32, Vec<MM128>)>>();

        seqs.iter()
            .zip(all_shmmers)
            .for_each(|((sid, seqname, seq), (_sid, shmmrs))| {
                let compress_seq = self.seq_to_compressed(seqname.clone(), *sid, seq, shmmrs, true);
                self.seqs.push(compress_seq);
            });

        Ok(())
    }
}
