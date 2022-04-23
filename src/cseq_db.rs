use crate::agc_io::AGCFile;
use crate::fasta_io::{reverse_complement, FastaReader, SeqRec};
use crate::shmmrutils::{match_reads, sequence_to_shmmrs, DeltaPoint, ShmmrSpec, MM128};
use byteorder::{LittleEndian, WriteBytesExt};
use flate2::bufread::MultiGzDecoder;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fmt;
use std::fs::File;
use std::io::{self, BufReader, Read, Write};


pub const KMERSIZE: u32 = 56;
pub const SHMMRSPEC: ShmmrSpec = ShmmrSpec {
    w: 80,
    k: KMERSIZE,
    r: 4,
    min_span: 64,
    sketch: true,
};

pub type Bases = Vec<u8>;
pub type AlnSegments = (u32, bool, Vec<AlnSegment>);

#[derive(Debug, Clone)]
pub enum AlnSegment {
    // this still use a lot of space, we will find way to reduce the memory footprint later
    FullMatch,
    Match(u32, u32),
    Insertion(u8),
}

enum GZFastaReader {
    GZFile(FastaReader<BufReader<MultiGzDecoder<BufReader<File>>>>),
    RegularFile(FastaReader<BufReader<BufReader<File>>>),
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

pub type ShmmrPair = (u64, u64);
pub type Fragments = Vec<Fragment>;
pub type FragmentSignature = (u32, u32, u32, u32, u8); //frg_id, seq_id, bgn, end, orientation(to shimmer pair)
pub type ShmmrToFrags = FxHashMap<ShmmrPair, Vec<FragmentSignature>>;
pub struct CompressedSeq {
    pub name: String,
    pub id: u32,
    pub shmmrs: Vec<MM128>,
    pub seq_frags: Vec<u32>,
    pub len: usize,
}

pub struct CompressedSeqDB {
    pub seqs: Vec<CompressedSeq>,
    pub frag_map: ShmmrToFrags,
    pub frags: Fragments,
}

pub fn deltas_to_aln_segs(
    deltas: &Vec<DeltaPoint>,
    endx: usize,
    endy: usize,
    frg: &Vec<u8>,
) -> Vec<AlnSegment> {
    let mut aln_segs = Vec::<AlnSegment>::new();
    if deltas.len() == 0 {
        aln_segs.push(AlnSegment::FullMatch);
        //println!("aln_segs: {:?}", aln_segs);
        return aln_segs;
    }

    let mut x = endx;
    let mut y = endy;

    for yy in (y..frg.len()).rev() {
        aln_segs.push(AlnSegment::Insertion(frg[yy]));
    }

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
    pub fn new() -> Self {
        let seqs = Vec::<CompressedSeq>::new();
        let frag_map = ShmmrToFrags::default();
        let frags = Vec::<Fragment>::new();
        CompressedSeqDB {
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
        let mut seq_frags = Vec::<u32>::new();
        let mut frg_id = self.frags.len() as u32;

        //assert!(shmmrs.len() > 0);
        if shmmrs.len() == 0 {
            let frg = seq[..].to_vec();
            self.frags.push(Fragment::Prefix(frg));
            seq_frags.push(frg_id);
            frg_id += 1;

            let frg = Vec::<u8>::new();
            self.frags.push(Fragment::Suffix(frg));
            seq_frags.push(frg_id);

            return CompressedSeq {
                name,
                id,
                shmmrs,
                seq_frags,
                len: seq.len(),
            };
        }
        // prefix
        let end = (shmmrs[0].pos() + 1) as usize;
        let frg = seq[..end].to_vec();
        self.frags.push(Fragment::Prefix(frg));
        seq_frags.push(frg_id);
        frg_id += 1;

        let shmmr_pairs = shmmrs[0..shmmrs.len() - 1]
            .iter()
            .zip(shmmrs[1..shmmrs.len()].iter())
            .collect::<Vec<_>>();

        let internal_frags = shmmr_pairs
            .par_iter()
            .map(|(shmmr0, shmmr1)| {
                let s0 = shmmr0.x >> 8;
                let s1 = shmmr1.x >> 8;
                let (shmmr_pair, orientation) = if s0 <= s1 {
                    ((s0, s1), 0_u8)
                } else {
                    ((s1, s0), 1_u8)
                };
                //let shmmr_pair = (shmmr0.x >> 8, shmmr1.x >> 8);
                let bgn = shmmr0.pos() + 1;
                let end = shmmr1.pos() + 1;
                let frg_len = end - bgn;
                let mut aligned = false;
                let mut out_frag = None;

                if frg_len > 64 && try_compress && self.frag_map.contains_key(&shmmr_pair) {
                    let e = self.frag_map.get(&shmmr_pair).unwrap();
                    for t_frg_id in e.iter() {
                        let base_frg = self.frags.get(t_frg_id.0 as usize).unwrap();
                        if let Fragment::Internal(b) = base_frg {
                            let base_frg = b;
                            //assert!(base_frg.len() > KMERSIZE as usize);
                            let frg;
                            let rc;
                            if orientation != t_frg_id.4 {
                                frg = reverse_complement(
                                    &seq[(bgn - KMERSIZE) as usize..end as usize].to_vec(),
                                );
                                rc = true;
                            } else {
                                frg = seq[(bgn - KMERSIZE) as usize..end as usize].to_vec();
                                rc = false;
                            }
                            //assert!(frg.len() > KMERSIZE as usize);
                            let m = match_reads(base_frg, &frg, true, 0.1, 0, 0, 32);
                            if let Some(m) = m {
                                let deltas: Vec<DeltaPoint> = m.deltas.unwrap();
                                let aln_segs = deltas_to_aln_segs(
                                    &deltas,
                                    m.end0 as usize,
                                    m.end1 as usize,
                                    &frg,
                                );

                                out_frag = Some((
                                    shmmr_pair,
                                    Fragment::AlnSegments((t_frg_id.0, rc, aln_segs)),
                                    bgn,
                                    end,
                                    orientation,
                                ));
                                aligned = true;
                                break; // we aligned to the first one of the fragments
                            } else {
                                continue;
                            }
                        }
                    }
                };

                if !aligned || !try_compress {
                    let frg = seq[(bgn - KMERSIZE) as usize..end as usize].to_vec();
                    out_frag = Some((shmmr_pair, Fragment::Internal(frg), bgn, end, orientation));
                };
                out_frag
            })
            .collect::<Vec<_>>();

        // TODO: parallize by sharding the key
        internal_frags.iter().for_each(|v| match v {
            Some((shmmr, frg, bgn, end, orientation)) => {
                if !self.frag_map.contains_key(shmmr) {
                    self.frag_map
                        .insert(*shmmr, Vec::<(u32, u32, u32, u32, u8)>::new());
                }
                let e = self.frag_map.get_mut(shmmr).unwrap();
                e.push((frg_id, id, *bgn, *end, *orientation));
                self.frags.push(frg.clone());
                seq_frags.push(frg_id);
                frg_id += 1;
            }
            None => {}
        });

        // suffix
        let bgn = (shmmrs[shmmrs.len() - 1].pos() + 1) as usize;
        let frg = seq[bgn..].to_vec();
        self.frags.push(Fragment::Suffix(frg));
        seq_frags.push(frg_id);

        CompressedSeq {
            name,
            id,
            shmmrs,
            seq_frags,
            len: seq.len(),
        }
    }

    pub fn seq_to_index(
        &mut self,
        name: String,
        id: u32,
        seqlen: usize,
        shmmrs: Vec<MM128>,
    ) -> CompressedSeq {
        let mut seq_frags = Vec::<u32>::new();
        let mut frg_id = self.frags.len() as u32;

        //assert!(shmmrs.len() > 0);
        if shmmrs.len() == 0 {
            return CompressedSeq {
                name,
                id,
                shmmrs,
                seq_frags,
                len: seqlen,
            };
        }

        seq_frags.push(frg_id);
        self.frags.push(Fragment::Prefix(vec![]));
        frg_id += 1;

        let shmmr_pairs = shmmrs[0..shmmrs.len() - 1]
            .iter()
            .zip(shmmrs[1..shmmrs.len()].iter())
            .collect::<Vec<_>>();

        let internal_frags = shmmr_pairs
            .par_iter()
            .map(|(shmmr0, shmmr1)| {
                let s0 = shmmr0.x >> 8;
                let s1 = shmmr1.x >> 8;
                let (shmmr_pair, orientation) = if s0 <= s1 {
                    ((s0, s1), 0_u8)
                } else {
                    ((s1, s0), 1_u8)
                };
                let bgn = shmmr0.pos() + 1;
                let end = shmmr1.pos() + 1;
                (shmmr_pair, bgn, end, orientation)
            })
            .collect::<Vec<_>>();

        internal_frags
            .iter()
            .for_each(|(shmmr, bgn, end, orientation)| {
                let e = self.frag_map.entry(*shmmr).or_insert(vec![]);
                e.push((frg_id, id, *bgn, *end, *orientation));
                seq_frags.push(frg_id);
                self.frags.push(Fragment::Internal(vec![]));
                frg_id += 1;
            });

        seq_frags.push(frg_id);
        self.frags.push(Fragment::Suffix(vec![]));

        CompressedSeq {
            name,
            id,
            shmmrs,
            seq_frags,
            len: seqlen,
        }
    }

    fn get_fastx_reader(&mut self, filepath: String) -> Result<GZFastaReader, std::io::Error> {
        let file = File::open(&filepath)?;
        let mut reader = BufReader::new(file);
        let mut is_gzfile = false;
        {
            let r = reader.by_ref();
            let mut buf = Vec::<u8>::new();
            let _ = r.take(2).read_to_end(&mut buf);
            if buf == [0x1F_u8, 0x8B_u8] {
                log::info!("input file: {} detected as gz-compressed file", filepath);
                is_gzfile = true;
            }
        }
        drop(reader);

        let file = File::open(&filepath)?;
        let reader = BufReader::new(file);
        let gz_buf = BufReader::new(MultiGzDecoder::new(reader));

        let file = File::open(&filepath)?;
        let reader = BufReader::new(file);
        let std_buf = BufReader::new(reader);

        if is_gzfile {
            drop(std_buf);
            Ok(GZFastaReader::GZFile(
                FastaReader::new(gz_buf, &filepath).unwrap(),
            ))
        } else {
            drop(gz_buf);
            Ok(GZFastaReader::RegularFile(
                FastaReader::new(std_buf, &filepath).unwrap(),
            ))
        }
    }

    fn get_shmmrs_from_seqs(
        &mut self,
        seqs: &Vec<(u32, String, Vec<u8>)>,
    ) -> Vec<(u32, Vec<MM128>)> {
        let all_shmmers = seqs
            .par_iter()
            .map(|(sid, _seqname, seq)| {
                let shmmrs = sequence_to_shmmrs(*sid, &seq, SHMMRSPEC);
                //let shmmrs = sequence_to_shmmrs2(*sid, &seq, 80, KMERSIZE, 4);
                (*sid, shmmrs)
            })
            .collect::<Vec<(u32, Vec<MM128>)>>();
        all_shmmers
    }

    fn load_seq_from_reader(&mut self, reader: &mut dyn Iterator<Item = io::Result<SeqRec>>) {
        let mut seqs = <Vec<(u32, String, Vec<u8>)>>::new();
        let mut sid = 0;
        loop {
            let mut count = 0;
            let mut end_ext_loop = false;
            seqs.clear();

            loop {
                if let Some(rec) = reader.next() {
                    let rec = rec.unwrap();
                    let seqname = String::from_utf8_lossy(&rec.id).into_owned();
                    seqs.push((sid, seqname, rec.seq));
                    sid += 1;
                } else {
                    end_ext_loop = true;
                    break;
                }
                count += 1;
                if count > 128 {
                    break;
                }
            }

            self.load_seqs_from_seq_vec(&seqs);
            if end_ext_loop {
                break;
            }
        }
        ();
    }

    fn load_seqs_from_seq_vec(&mut self, seqs: &Vec<(u32, String, Vec<u8>)>) {
        let all_shmmers = self.get_shmmrs_from_seqs(seqs);
        seqs.iter()
            .zip(all_shmmers)
            .for_each(|((sid, seqname, seq), (_sid, shmmrs))| {
                let compress_seq = self.seq_to_compressed(seqname.clone(), *sid, seq, shmmrs, true);
                self.seqs.push(compress_seq);
            });
    }

    pub fn load_seqs(&mut self, filepath: String) -> Result<(), std::io::Error> {
        match self.get_fastx_reader(filepath)? {
            GZFastaReader::GZFile(reader) => self.load_seq_from_reader(&mut reader.into_iter()),

            GZFastaReader::RegularFile(reader) => {
                self.load_seq_from_reader(&mut reader.into_iter())
            }
        };

        Ok(())
    }

    fn load_index_from_reader(&mut self, reader: &mut dyn Iterator<Item = io::Result<SeqRec>>) {
        let mut seqs = <Vec<(u32, String, Vec<u8>)>>::new();
        let mut sid = 0;
        loop {
            let mut count = 0;
            let mut end_ext_loop = false;
            seqs.clear();

            loop {
                if let Some(rec) = reader.next() {
                    let rec = rec.unwrap();
                    let seqname = String::from_utf8_lossy(&rec.id).into_owned();
                    seqs.push((sid, seqname, rec.seq));
                    sid += 1;
                } else {
                    end_ext_loop = true;
                    break;
                }
                count += 1;
                if count > 128 {
                    break;
                }
            }

            self.load_index_from_seq_vec(&seqs);
            if end_ext_loop {
                break;
            }
        }
        ();
    }

    fn load_index_from_seq_vec(&mut self, seqs: &Vec<(u32, String, Vec<u8>)>) {
        let all_shmmers = self.get_shmmrs_from_seqs(seqs);
        let seq_names = seqs
            .iter()
            .map(|(_sid, n, s)| (n.clone(), s.len()))
            .collect::<Vec<(String, usize)>>();
        seq_names
            .iter()
            .zip(all_shmmers)
            .for_each(|((seq_name, seqlen), (sid, shmmrs))| {
                let compress_seq = self.seq_to_index(seq_name.clone(), sid, *seqlen, shmmrs);
                self.seqs.push(compress_seq);
            });
    }

    fn _write_shmmr_vec_from_reader(
        &mut self,
        reader: &mut dyn Iterator<Item = io::Result<SeqRec>>,
        writer: &mut Vec<u8>,
    ) {
        let mut seqs = <Vec<(u32, String, Vec<u8>)>>::new();
        let mut sid = 0;
        loop {
            let mut count = 0;
            let mut end_ext_loop = false;
            seqs.clear();

            loop {
                if let Some(rec) = reader.next() {
                    let rec = rec.unwrap();
                    let seqname = String::from_utf8_lossy(&rec.id).into_owned();
                    seqs.push((sid, seqname, rec.seq));
                    sid += 1;
                } else {
                    end_ext_loop = true;
                    break;
                }
                count += 1;
                if count > 128 {
                    break;
                }
            }

            self.get_shmmrs_from_seqs(&seqs)
                .iter()
                .map(|(_, v)| v)
                .for_each(|v| {
                    v.iter().for_each(|m| {
                        let _ = writer.write_u64::<LittleEndian>(m.x);
                        let _ = writer.write_u64::<LittleEndian>(m.y);
                    });
                });

            if end_ext_loop {
                break;
            }
        }
        ();
    }

    pub fn load_index_from_fastx(&mut self, filepath: String) -> Result<(), std::io::Error> {
        match self.get_fastx_reader(filepath)? {
            GZFastaReader::GZFile(reader) => self.load_index_from_reader(&mut reader.into_iter()),

            GZFastaReader::RegularFile(reader) => {
                self.load_index_from_reader(&mut reader.into_iter())
            }
        };

        Ok(())
    }

    pub fn load_index_from_agcfile(&mut self, agcfile: AGCFile) -> Result<(), std::io::Error> {
        //let agcfile = AGCFile::new(filepath);

        self.load_index_from_reader(&mut agcfile.into_iter());
        Ok(())
    }
}

pub fn query_fragment(
    shmmr_map: &ShmmrToFrags,
    frag: &Vec<u8>,
) -> Vec<((u64, u64), (u32, u32, u8), Vec<FragmentSignature>)> {
    let shmmrs = sequence_to_shmmrs(0, &frag, SHMMRSPEC);
    let query_results = shmmrs[0..shmmrs.len() - 1]
        .iter()
        .zip(shmmrs[1..shmmrs.len()].iter())
        .map(|(s0, s1)| {
            let p0 = s0.pos();
            let p1 = s1.pos();
            let s0 = s0.x >> 8;
            let s1 = s1.x >> 8;
            if s0 < s1 {
                (s0, s1, p0, p1, 0_u8)
            } else {
                (s1, s0, p0, p1, 0_u8)
            }
        })
        .map(|(s0, s1, p0, p1, orientation)| {
            if let Some(m) = shmmr_map.get(&(s0, s1)) {
                ((s0, s1), (p0, p1, orientation), m.clone())
            } else {
                ((s0, s1), (p0, p1, orientation), vec![])
            }
        })
        .collect::<Vec<_>>();
    query_results
}

pub fn write_shmr_map_file(shmmr_map: &ShmmrToFrags, filepath: String) {
    let mut out_file = File::create(filepath).expect("open fail");
    let mut buf = Vec::<u8>::new();
    buf.write_u64::<LittleEndian>(shmmr_map.len() as u64);
    shmmr_map.into_iter().for_each(|(k,v) | {
        buf.write_u64::<LittleEndian>(k.0);
        buf.write_u64::<LittleEndian>(k.1);
        buf.write_u64::<LittleEndian>(v.len() as u64);
        v.iter().for_each(|r| {
            buf.write_u32::<LittleEndian>(r.0);
            buf.write_u32::<LittleEndian>(r.1);
            buf.write_u32::<LittleEndian>(r.2);
            buf.write_u32::<LittleEndian>(r.3);
            buf.write_u8(r.4);
        })
    } );
    let _ = out_file.write_all(&buf);
}


pub fn read_shmr_map_file(filepath: String) -> ShmmrToFrags {

    let mut in_file = File::open(filepath).expect("open fail");
    let mut buf = Vec::<u8>::new();
    let mut u64bytes = [0_u8;8];
    let mut u32bytes = [0_u8;4];
    in_file.read_to_end(&mut buf);
    let mut cursor = 0_usize;
    u64bytes.clone_from_slice(&buf[cursor..cursor+8]);
    let shmmr_key_len = usize::from_le_bytes(u64bytes);
    cursor += 8;
    let mut shmmr_map = ShmmrToFrags::default();
    (0..shmmr_key_len).into_iter().for_each( |_| {

        u64bytes.clone_from_slice(&buf[cursor..cursor+8]);
        let k1 = u64::from_le_bytes(u64bytes);
        cursor += 8;
        
        u64bytes.clone_from_slice(&buf[cursor..cursor+8]);
        let k2 = u64::from_le_bytes(u64bytes);
        cursor += 8;

        u64bytes.clone_from_slice(&buf[cursor..cursor+8]);
        let vec_len = usize::from_le_bytes(u64bytes);
        cursor += 8;

        let value = (0..vec_len).into_iter().map(|_| {
            let mut v = (0_u32, 0_u32, 0_u32, 0_u32, 0_u8);

            u32bytes.clone_from_slice(&buf[cursor..cursor+4]);
            v.0 = u32::from_le_bytes(u32bytes);
            cursor += 4;
            
            u32bytes.clone_from_slice(&buf[cursor..cursor+4]);
            v.1 = u32::from_le_bytes(u32bytes);
            cursor += 4;
            
            u32bytes.clone_from_slice(&buf[cursor..cursor+4]);
            v.2 = u32::from_le_bytes(u32bytes);
            cursor += 4;
            
            u32bytes.clone_from_slice(&buf[cursor..cursor+4]);
            v.3 = u32::from_le_bytes(u32bytes);
            cursor += 4;

            v.4 = buf[cursor..cursor+1][0];
            cursor += 1;

           v 
        }).collect::<Vec<(u32, u32, u32, u32, u8)>>();

        shmmr_map.insert((k1,  k2), value);
        
    });

    shmmr_map
}
