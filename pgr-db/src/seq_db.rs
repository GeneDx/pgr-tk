use crate::agc_io::AGCFile;
use crate::fasta_io::{reverse_complement, FastaReader, SeqRec};
use crate::graph_utils::{AdjList, AdjPair, ShmmrGraphNode};
use crate::shmmrutils::{match_reads, sequence_to_shmmrs, DeltaPoint, ShmmrSpec, MM128};
use bincode::{config, Decode, Encode};
use byteorder::{ByteOrder, LittleEndian, WriteBytesExt};
use flate2::bufread::MultiGzDecoder;
use flate2::write::DeflateEncoder;
use flate2::Compression;
use memmap2::Mmap;
use petgraph::graphmap::DiGraphMap;
use petgraph::visit::Dfs;
use petgraph::EdgeDirection::{Incoming, Outgoing};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use std::fmt;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};

pub const KMERSIZE: u32 = 56;
pub const SHMMRSPEC: ShmmrSpec = ShmmrSpec {
    w: 80,
    k: KMERSIZE,
    r: 4,
    min_span: 64,
    sketch: true,
};

pub type Bases = Vec<u8>;
pub type AlnSegments = (u32, bool, u32, Vec<AlnSegment>); //(refFragID, orientation, SeqLength, AlnSegments)

#[derive(Debug, Clone, Decode, Encode)]
pub enum AlnSegment {
    // this still use a lot of space, we will find way to reduce the memory footprint later
    FullMatch,
    // u16 should be enough, the max span should be less than 128 * 144 = 18423 * 2 < 2**16
    Match(u32, u32),
    Insertion(u8),
}
#[allow(clippy::large_enum_variant)]
enum GZFastaReader {
    GZFile(FastaReader<BufReader<MultiGzDecoder<BufReader<File>>>>),
    RegularFile(FastaReader<BufReader<BufReader<File>>>),
}

#[derive(Debug, Clone, Decode, Encode)]
pub enum Fragment {
    // size = 40, align = 8
    AlnSegments(AlnSegments),
    Prefix(Bases),
    Internal(Bases),
    Suffix(Bases),
}

impl fmt::Display for Fragment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Fragment::AlnSegments(d) => write!(
                f,
                "RefFragId:{} Orientation:{:?} len:{:?} AlignSegs:{:?}",
                d.0, d.1, d.2, d.3
            ),
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
pub type ShmmrIndexFileLocation = Vec<(ShmmrPair, (usize, usize))>;
pub type ShmmrToIndexFileLocation = FxHashMap<ShmmrPair, (usize, usize)>;

pub trait GetSeq {
    fn get_seq_by_id(&self, sid: u32) -> Vec<u8>;
    fn get_sub_seq_by_id(&self, sid: u32, bgn: u32, end: u32) -> Vec<u8>;
}

#[derive(Debug, Clone, Decode, Encode)]
pub struct CompactSeq {
    pub source: Option<String>,
    pub name: String,
    pub id: u32,
    pub seq_frag_range: (u32, u32), // (start, len)
    pub len: usize,
}

#[derive(Debug, Clone)]
pub struct CompactSeqDB {
    pub shmmr_spec: ShmmrSpec,
    pub seqs: Vec<CompactSeq>,
    pub frag_map: ShmmrToFrags,
    pub frags: Option<Fragments>,
}

pub fn pair_shmmrs(shmmrs: &Vec<MM128>) -> Vec<(&MM128, &MM128)> {
    if shmmrs.len() < 2 {
        return vec![];
    }
    let shmmr_pairs = shmmrs[0..shmmrs.len() - 1]
        .iter()
        .zip(shmmrs[1..shmmrs.len()].iter())
        .collect::<Vec<_>>();
    shmmr_pairs
}

pub fn deltas_to_aln_segs(
    deltas: &Vec<DeltaPoint>,
    endx: usize,
    endy: usize,
    base_frg: &Vec<u8>,
    frg: &Vec<u8>,
) -> Vec<AlnSegment> {
    let mut aln_segs = Vec::<AlnSegment>::new();
    if deltas.is_empty() && base_frg.len() == frg.len() {
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

pub fn reconstruct_seq_from_aln_segs(base_seq: &[u8], aln_segs: &[AlnSegment]) -> Vec<u8> {
    let mut seq = Vec::<u8>::new();
    for s in aln_segs.iter() {
        match s {
            AlnSegment::FullMatch => {
                seq.extend_from_slice(base_seq);
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

impl CompactSeqDB {
    pub fn new(shmmr_spec: ShmmrSpec) -> Self {
        let seqs = Vec::<CompactSeq>::new();
        let frag_map = ShmmrToFrags::default();
        let frags = None;
        CompactSeqDB {
            shmmr_spec,
            seqs,
            frag_map,
            frags,
        }
    }

    pub fn seq_to_compressed(
        &mut self,
        source: Option<String>,
        name: String,
        id: u32,
        seq: &Vec<u8>,
        shmmrs: Vec<MM128>,
        try_compress: bool,
    ) -> CompactSeq {
        let mut seq_frags = Vec::<u32>::new();

        assert!(self.frags.is_some());
        let frags: &mut Vec<Fragment> = self.frags.as_mut().unwrap();

        let mut frg_id = frags.len() as u32;
        let mut seq_len = 0_usize;

        //assert!(shmmrs.len() > 0);
        if shmmrs.is_empty() {
            let frg = seq[..].to_vec();
            frags.push(Fragment::Prefix(frg));
            seq_frags.push(frg_id);
            // frg_id += 1;

            let frg = Vec::<u8>::new();
            frags.push(Fragment::Suffix(frg));
            seq_frags.push(frg_id);

            return CompactSeq {
                source,
                name,
                id,
                seq_frag_range: (seq_frags[0] as u32, seq_frags.len() as u32),
                len: seq.len(),
            };
        }
        // prefix
        let end = (shmmrs[0].pos() + 1) as usize;
        let frg = seq[..end].to_vec();
        seq_len += frg.len();
        frags.push(Fragment::Prefix(frg));
        seq_frags.push(frg_id);
        frg_id += 1;

        let internal_frags = pair_shmmrs(&shmmrs)
            .par_iter()
            .map(|(shmmr0, shmmr1)| {
                let s0 = shmmr0.hash();
                let s1 = shmmr1.hash();
                let (shmmr_pair, orientation) = if s0 <= s1 {
                    ((s0, s1), 0_u8)
                } else {
                    ((s1, s0), 1_u8)
                };
                let bgn = shmmr0.pos() + 1;
                let end = shmmr1.pos() + 1;
                let frg_len = end - bgn;
                let mut aligned = false;
                let mut out_frag = None;

                if frg_len > 128 && try_compress && self.frag_map.contains_key(&shmmr_pair) {
                    let e = self.frag_map.get(&shmmr_pair).unwrap();
                    for t_frg_id in e.iter() {
                        let base_frg = frags.get(t_frg_id.0 as usize).unwrap();
                        if let Fragment::Internal(b) = base_frg {
                            let base_frg = b;
                            //assert!(base_frg.len() > KMERSIZE as usize);
                            let frg;
                            let rc;
                            if orientation != t_frg_id.4 {
                                frg = reverse_complement(
                                    &seq[(bgn - self.shmmr_spec.k) as usize..end as usize],
                                );
                                rc = true;
                            } else {
                                frg =
                                    seq[(bgn - self.shmmr_spec.k) as usize..end as usize].to_vec();
                                rc = false;
                            }
                            //assert!(frg.len() > KMERSIZE as usize);
                            //the max span should be less than 128 * 144 = 18423 * 2 < 2**16
                            assert!(base_frg.len() < (1 << 32) - 1);
                            let m = match_reads(base_frg, &frg, true, 0.1, 0, 0, 32);
                            if let Some(m) = m {
                                let deltas: Vec<DeltaPoint> = m.deltas.unwrap();
                                let aln_segs = deltas_to_aln_segs(
                                    &deltas,
                                    m.end0 as usize,
                                    m.end1 as usize,
                                    &base_frg,
                                    &frg,
                                );

                                /*  // For debugging
                                let out = reconstruct_seq_from_aln_segs(base_frg, &aln_segs);
                                if out != frg {
                                    println!("DBG: {:?} {} {}", deltas, m.end0, m.end1);
                                    println!("DBG: {:?} {:?} {:?} {} {}", String::from_utf8_lossy(base_frg), aln_segs,
                                    String::from_utf8_lossy(&frg), base_frg.len(), frg.len());
                                };
                                assert_eq!(out, frg);
                                */

                                if std::mem::align_of_val(&aln_segs) > (frg.len() >> 2) {
                                    continue;
                                }

                                out_frag = Some((
                                    shmmr_pair,
                                    Fragment::AlnSegments((
                                        t_frg_id.0,
                                        rc,
                                        frg.len() as u32,
                                        aln_segs,
                                    )),
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
                    let frg = seq[(bgn - self.shmmr_spec.k) as usize..end as usize].to_vec();
                    out_frag = Some((shmmr_pair, Fragment::Internal(frg), bgn, end, orientation));
                };
                out_frag
            })
            .collect::<Vec<_>>();

        // TODO: parallelize by sharding the key
        internal_frags.iter().for_each(|v| match v {
            Some((shmmr, frg, bgn, end, orientation)) => {
                if !self.frag_map.contains_key(shmmr) {
                    self.frag_map
                        .insert(*shmmr, Vec::<(u32, u32, u32, u32, u8)>::new());
                }
                let e = self.frag_map.get_mut(shmmr).unwrap();
                e.push((frg_id, id, *bgn, *end, *orientation));
                seq_len += (*end - *bgn) as usize;
                frags.push(frg.clone());
                seq_frags.push(frg_id);
                frg_id += 1;
            }
            None => {}
        });

        // suffix
        let bgn = (shmmrs[shmmrs.len() - 1].pos() + 1) as usize;
        let frg = seq[bgn..].to_vec();
        seq_len += frg.len();
        frags.push(Fragment::Suffix(frg));
        seq_frags.push(frg_id);

        assert_eq!(seq_len, seq.len());
        CompactSeq {
            source,
            name,
            id,
            seq_frag_range: (seq_frags[0] as u32, seq_frags.len() as u32),
            len: seq.len(),
        }
    }

    #[allow(clippy::type_complexity)]
    pub fn seq_to_index(
        source: Option<String>,
        name: String,
        id: u32,
        seqlen: usize,
        shmmrs: Vec<MM128>,
    ) -> (CompactSeq, Vec<((u64, u64), u32, u32, u8)>) {
        //assert!(shmmrs.len() > 0);
        if shmmrs.is_empty() {
            return (
                CompactSeq {
                    source,
                    name,
                    id,
                    seq_frag_range: (0, 0),
                    len: seqlen,
                },
                vec![],
            );
        }

        let shmmr_pairs = shmmrs[0..shmmrs.len() - 1]
            .iter()
            .zip(shmmrs[1..shmmrs.len()].iter())
            .collect::<Vec<_>>();

        let internal_frags: Vec<((u64, u64), u32, u32, u8)> = shmmr_pairs
            .par_iter()
            .map(|(shmmr0, shmmr1)| {
                let s0 = shmmr0.hash();
                let s1 = shmmr1.hash();
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

        let seq_frags: Vec<u32> = (0..(shmmr_pairs.len() as u32)).collect();
        let seq_frag_range = if seq_frags.is_empty() {
            (0, 0)
        } else {
            (seq_frags[0] as u32, seq_frags.len() as u32)
        };
        (
            CompactSeq {
                source,
                name,
                id,
                seq_frag_range,
                len: seqlen,
            },
            internal_frags,
        )
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
                FastaReader::new(gz_buf, &filepath, 1 << 14, true).unwrap(),
            ))
        } else {
            drop(gz_buf);
            Ok(GZFastaReader::RegularFile(
                FastaReader::new(std_buf, &filepath, 1 << 14, true).unwrap(),
            ))
        }
    }

    fn get_shmmrs_from_seqs(
        &mut self,
        seqs: &Vec<(u32, Option<String>, String, Vec<u8>)>,
    ) -> Vec<(u32, Vec<MM128>)> {
        let all_shmmrs = seqs
            .par_iter()
            .map(|(sid, _, _, seq)| {
                let shmmrs = sequence_to_shmmrs(*sid, seq, &self.shmmr_spec, false);
                //let shmmrs = sequence_to_shmmrs2(*sid, &seq, 80, KMERSIZE, 4);
                (*sid, shmmrs)
            })
            .collect::<Vec<(u32, Vec<MM128>)>>();
        all_shmmrs
    }

    fn load_seq_from_reader(&mut self, reader: &mut dyn Iterator<Item = io::Result<SeqRec>>) {
        let mut seqs = <Vec<(u32, Option<String>, String, Vec<u8>)>>::new();
        let mut sid = self.seqs.len() as u32;
        if self.frags.is_none() {
            self.frags = Some(Fragments::new());
        };

        loop {
            let mut count = 0;
            let mut end_ext_loop = false;
            seqs.clear();

            loop {
                if let Some(rec) = reader.next() {
                    let rec = rec.unwrap();
                    let source = rec.source.clone();
                    let seqname = String::from_utf8_lossy(&rec.id).into_owned();
                    seqs.push((sid, source, seqname, rec.seq));
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
    }

    pub fn load_seqs_from_seq_vec(&mut self, seqs: &Vec<(u32, Option<String>, String, Vec<u8>)>) {
        if self.frags.is_none() {
            self.frags = Some(Fragments::new());
        }
        let all_shmmrs = self.get_shmmrs_from_seqs(seqs);
        seqs.iter()
            .zip(all_shmmrs)
            .for_each(|((sid, source, seqname, seq), (_sid, shmmrs))| {
                let compress_seq = self.seq_to_compressed(
                    source.clone(),
                    seqname.clone(),
                    *sid,
                    seq,
                    shmmrs,
                    true,
                );
                self.seqs.push(compress_seq);
            });
    }

    pub fn load_seqs_from_fastx(&mut self, filepath: String) -> Result<(), std::io::Error> {
        match self.get_fastx_reader(filepath)? {
            #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
            GZFastaReader::GZFile(reader) => self.load_seq_from_reader(&mut reader.into_iter()),

            #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
            GZFastaReader::RegularFile(reader) => {
                self.load_seq_from_reader(&mut reader.into_iter())
            }
        };

        Ok(())
    }

    fn load_index_from_reader(&mut self, reader: &mut dyn Iterator<Item = io::Result<SeqRec>>) {
        let mut seqs = <Vec<(u32, Option<String>, String, Vec<u8>)>>::new();
        let mut sid = 0;
        loop {
            let mut count = 0;
            let mut end_ext_loop = false;
            seqs.clear();

            loop {
                if let Some(rec) = reader.next() {
                    let rec = rec.unwrap();
                    let source = rec.source;
                    let seqname = String::from_utf8_lossy(&rec.id).into_owned();
                    seqs.push((sid, source, seqname, rec.seq));
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
    }

    pub fn load_index_from_seq_vec(&mut self, seqs: &Vec<(u32, Option<String>, String, Vec<u8>)>) {
        let all_shmmrs = self.get_shmmrs_from_seqs(seqs);
        let seq_names = seqs
            .iter()
            .map(|(_sid, src, n, s)| (src.clone(), n.clone(), s.len()))
            .collect::<Vec<(Option<String>, String, usize)>>();

        /*
        seq_names.iter().zip(all_shmmrs).for_each(
            |((source, seq_name, seqlen), (sid, shmmrs))| {
                let compress_seq =
                    self._seq_to_index(source.clone(), seq_name.clone(), sid, *seqlen, shmmrs);
                self.seqs.push(compress_seq);
            },
        );
        */

        seq_names
            .par_iter()
            .zip(all_shmmrs)
            .map(|((source, seq_name, seqlen), (sid, shmmrs))| {
                let tmp = self::CompactSeqDB::seq_to_index(
                    source.clone(),
                    seq_name.clone(),
                    sid,
                    *seqlen,
                    shmmrs,
                );
                (sid, tmp.0, tmp.1)
            })
            .collect::<Vec<(u32, CompactSeq, Vec<_>)>>()
            .into_iter()
            .for_each(|(sid, cs, internal_frags)| {
                internal_frags
                    .iter()
                    .zip(cs.seq_frag_range.0..cs.seq_frag_range.0 + cs.seq_frag_range.1)
                    .for_each(|((shmmr, bgn, end, orientation), frg_id)| {
                        let e = self.frag_map.entry(*shmmr).or_default();
                        e.push((frg_id, sid, *bgn, *end, *orientation));
                    });
                self.seqs.push(cs);
            });
    }

    fn _write_shmmr_vec_from_reader(
        &mut self,
        reader: &mut dyn Iterator<Item = io::Result<SeqRec>>,
        writer: &mut Vec<u8>,
    ) {
        let mut seqs = <Vec<(u32, Option<String>, String, Vec<u8>)>>::new();
        let mut sid = 0;
        loop {
            let mut count = 0;
            let mut end_ext_loop = false;
            seqs.clear();

            loop {
                if let Some(rec) = reader.next() {
                    let rec = rec.unwrap();
                    let source = rec.source;
                    let seqname = String::from_utf8_lossy(&rec.id).into_owned();
                    seqs.push((sid, source, seqname, rec.seq));
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
    }

    pub fn load_index_from_fastx(&mut self, filepath: String) -> Result<(), std::io::Error> {
        match self.get_fastx_reader(filepath)? {
            #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
            GZFastaReader::GZFile(reader) => self.load_index_from_reader(&mut reader.into_iter()),

            #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
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

impl CompactSeqDB {
    fn reconstruct_seq_from_frags<I: Iterator<Item = u32>>(&self, frag_ids: I) -> Vec<u8> {
        let mut reconstructed_seq = <Vec<u8>>::new();
        let frags: &Vec<Fragment> = self.frags.as_ref().unwrap();
        // let mut _p = 0;
        frag_ids.for_each(|frag_id| {
            //println!("{}:{}", frg_id, sdb.frags[*frg_id as usize]);
            match frags.get(frag_id as usize).unwrap() {
                Fragment::Prefix(b) => {
                    reconstructed_seq.extend_from_slice(&b[..]);
                    //println!("P p: {} {} {}", frag_id, _p, _p + b.len());
                    //_p += b.len();
                }
                Fragment::Suffix(b) => {
                    reconstructed_seq.extend_from_slice(&b[..]);
                    //println!("S p: {} {} {}", frag_id, _p, _p + b.len());
                    //_p += b.len();
                }
                Fragment::Internal(b) => {
                    reconstructed_seq.extend_from_slice(&b[self.shmmr_spec.k as usize..]);
                    //println!("I p: {} {} {}", frag_id, _p, _p + b.len()-self.shmmr_spec.k as usize);
                    //_p += b.len()-self.shmmr_spec.k as usize;
                }
                Fragment::AlnSegments((frg_id, reversed, _length, a)) => {
                    if let Fragment::Internal(base_seq) = frags.get(*frg_id as usize).unwrap() {
                        let mut seq = reconstruct_seq_from_aln_segs(base_seq, a);
                        /*  // for debugging
                        if *_length as usize != seq.len() {
                            println!("DBG X: {:?} {:?}", String::from_utf8_lossy(base_seq), a);
                        }
                        */

                        assert_eq!(*_length as usize, seq.len());
                        if *reversed {
                            seq = reverse_complement(&seq);
                        }
                        reconstructed_seq.extend_from_slice(&seq[self.shmmr_spec.k as usize..]);
                        // println!("A p: {} {} {}", frag_id, _p, _p + seq.len()-self.shmmr_spec.k as usize);
                        // _p += seq.len()-self.shmmr_spec.k as usize;
                    }
                }
            }
        });

        reconstructed_seq
    }

    pub fn get_seq(&self, seq: &CompactSeq) -> Vec<u8> {
        self.reconstruct_seq_from_frags(
            (seq.seq_frag_range.0..seq.seq_frag_range.0 + seq.seq_frag_range.1).into_iter(),
        )
    }

    /* TODO */
    /*
    pub fn get_sub_seq(&self, seq: &CompactSeq, b: usize, e:usize) -> Vec<u8> {
        vec![]
    }
    */
}

impl GetSeq for CompactSeqDB {
    fn get_seq_by_id(&self, sid: u32) -> Vec<u8> {
        let seq = self.seqs.get(sid as usize).unwrap();
        self.reconstruct_seq_from_frags(
            (seq.seq_frag_range.0..seq.seq_frag_range.0 + seq.seq_frag_range.1).into_iter(),
        )
    }

    fn get_sub_seq_by_id(&self, sid: u32, bgn: u32, end: u32) -> Vec<u8> {
        assert!((sid as usize) < self.seqs.len());
        let frag_range = &self.seqs[sid as usize].seq_frag_range;

        let mut _p = 0;
        let mut base_offset = 0_u32;
        let mut sub_seq_frag = vec![];
        let frags: &Vec<Fragment> = self.frags.as_ref().unwrap();
        for frag_id in frag_range.0..frag_range.0 + frag_range.1 {
            let f = &frags[frag_id as usize];
            let frag_len = match f {
                Fragment::AlnSegments(d) => d.2 - self.shmmr_spec.k,
                Fragment::Prefix(b) => b.len() as u32,
                Fragment::Internal(b) => b.len() as u32 - self.shmmr_spec.k,
                Fragment::Suffix(b) => b.len() as u32,
            };
            if (base_offset <= bgn && bgn < base_offset + frag_len)
                || (base_offset <= end && end < base_offset + frag_len)
                || (bgn <= base_offset && base_offset + frag_len <= end)
            {
                sub_seq_frag.push((frag_id, base_offset));
            }

            base_offset += frag_len;
        }

        // println!("DBG0: {} {}" , sid,  sub_seq_frag.len() );
        let reconstructed_seq = self.reconstruct_seq_from_frags(sub_seq_frag.iter().map(|v| v.0));

        // println!("DBG: {} {} {} {} {}", sid, frag_range.1, sub_seq_frag.len(), end, reconstructed_seq.len() );

        let offset = bgn - sub_seq_frag[0].1;
        reconstructed_seq[(offset as usize)..((offset + end - bgn) as usize)].to_vec()
    }
}

impl CompactSeqDB {
    pub fn write_shmmr_map_index(&self, fp_prefix: String) -> Result<(), std::io::Error> {
        let seq_idx_fp = fp_prefix.clone() + ".midx";
        let data_fp = fp_prefix + ".mdb";
        write_shmmr_map_file(&self.shmmr_spec, &self.frag_map, data_fp)?;
        let mut idx_file = BufWriter::new(File::create(seq_idx_fp).expect("file create error"));
        self.seqs
            .iter()
            .try_for_each(|s| -> Result<(), std::io::Error> {
                writeln!(
                    idx_file,
                    "{}\t{}\t{}\t{}",
                    s.id,
                    s.len,
                    s.name,
                    s.source.clone().unwrap_or_else(|| "-".to_string())
                )?;
                Ok(())
            })?;

        Ok(())
    }
}

impl CompactSeqDB {
    pub fn write_to_frag_files(&self, file_prefix: String) {
        let mut sdx_file = BufWriter::new(
            File::create(file_prefix.clone() + ".sdx").expect("sdx file creating fail\n"),
        );
        let mut frg_file =
            BufWriter::new(File::create(file_prefix + ".frg").expect("frg file creating fail\n"));

        let config = config::standard();

        //self.seqs.iter().for_each(|s| {
        //    println!("{:?} {:?} {} {}", s.id, s.seq_frags.len(), s.seq_frags[0], s.seq_frags[s.seq_frags.len()-1]);
        //});

        let compressed_frags = self
            .frags
            .as_ref()
            .unwrap()
            .par_iter()
            .map(|f| {
                let frag_len = match f {
                    Fragment::AlnSegments(d) => d.2,
                    Fragment::Prefix(b) => b.len() as u32,
                    Fragment::Internal(b) => b.len() as u32,
                    Fragment::Suffix(b) => b.len() as u32,
                };
                let w = bincode::encode_to_vec(f, config).unwrap();
                let mut compressor = DeflateEncoder::new(Vec::new(), Compression::default());
                compressor.write_all(&w).unwrap();
                let compress_frag = compressor.finish().unwrap();
                (frag_len, compress_frag)
            })
            .collect::<Vec<(u32, Vec<u8>)>>();

        let mut frag_addr_offset = vec![];
        let mut offset = 0_usize;
        compressed_frags.iter().for_each(|(frag_len, v)| {
            let l = v.len();
            frag_addr_offset.push((offset, v.len(), *frag_len));
            offset += l;
            frg_file.write_all(v).expect("frag file writing error\n");
        });

        bincode::encode_into_std_write((frag_addr_offset, &self.seqs), &mut sdx_file, config)
            .expect("sdx file writing error\n");
        //bincode::encode_into_std_write(compressed_frags, &mut frg_file, config)
        //    .expect(" frag file writing error");
    }
}

pub fn frag_map_to_adj_list(
    frag_map: &ShmmrToFrags,
    min_count: usize,
    keeps: Option<Vec<u32>>, // a list of sequence id that we like to keep the sequence in the adj list regardless the coverage
) -> AdjList {
    let mut out = frag_map
        .par_iter()
        .flat_map(|v| {
            v.1.iter()
                .map(|vv| (vv.1, vv.2, vv.3, ShmmrGraphNode(v.0 .0, v.0 .1, vv.4)))
                .collect::<Vec<(u32, u32, u32, ShmmrGraphNode)>>() //(seq_id, bgn, end, (hash0, hash1, orientation))
        })
        .collect::<Vec<(u32, u32, u32, ShmmrGraphNode)>>();
    if out.len() < 2 {
        return vec![];
    }
    out.par_sort();

    let out = if let Some(keeps) = keeps {
        let keeps = FxHashSet::<u32>::from_iter(keeps.into_iter());

        // more or less duplicate code, but this takes the hashset check out of the loop if keeps is None.
        out.into_par_iter()
            .map(|v| {
                if frag_map.get(&(v.3 .0, v.3 .1)).unwrap().len() >= min_count
                    || keeps.contains(&v.0)
                {
                    Some(v)
                } else {
                    None
                }
            })
            .collect::<Vec<Option<(u32, u32, u32, ShmmrGraphNode)>>>()
    } else {
        out.into_par_iter()
            .map(|v| {
                if frag_map.get(&(v.3 .0, v.3 .1)).unwrap().len() >= min_count {
                    Some(v)
                } else {
                    None
                }
            })
            .collect::<Vec<Option<(u32, u32, u32, ShmmrGraphNode)>>>()
    };

    (0..out.len() - 1)
        //.into_par_iter()
        .into_iter()
        .flat_map(|i| {
            if let (Some(v), Some(w)) = (out[i], out[i + 1]) {
                // println!("DBG v: {} {} {} {:?} w: {} {} {} {:?}", v.0, v.1, v.2, v.3, w.0, w.1, w.2, w.3); // XXX
                if v.0 != w.0 || v.2 != w.1 {
                    vec![None]
                } else {
                    vec![
                        Some((v.0, v.3, w.3)),
                        Some((
                            v.0,
                            ShmmrGraphNode(w.3 .0, w.3 .1, 1 - w.3 .2),
                            ShmmrGraphNode(v.3 .0, v.3 .1, 1 - v.3 .2),
                        )),
                    ]
                }
            } else {
                vec![None]
            }
        })
        .filter(|v| v.is_some())
        .map(|v| v.unwrap())
        .collect::<AdjList>() // seq_id, node0, node1
}

pub fn generate_smp_adj_list_for_seq(
    seq: &Vec<u8>,
    sid: u32,
    frag_map: &ShmmrToFrags,
    shmmr_spec: &ShmmrSpec,
    min_count: usize,
) -> AdjList {
    let shmmrs = sequence_to_shmmrs(0, seq, shmmr_spec, false);
    let res = pair_shmmrs(&shmmrs)
        .iter()
        .map(|(s0, s1)| {
            let p0 = s0.pos() + 1;
            let p1 = s1.pos() + 1;
            let s0 = s0.x >> 8;
            let s1 = s1.x >> 8;
            if s0 < s1 {
                (s0, s1, p0, p1, 0_u8)
            } else {
                (s1, s0, p0, p1, 1_u8)
            }
        })
        .collect::<Vec<(u64, u64, u32, u32, u8)>>();

    if res.len() < 2 {
        vec![]
    } else {
        (0..res.len() - 1)
            .into_iter()
            .flat_map(|i| {
                let v = res[i];
                let w = res[i + 1];
                if (frag_map.get(&(v.0, v.1)).is_none() || frag_map.get(&(w.0, w.1)).is_none())
                    || (frag_map.get(&(v.0, v.1)).unwrap().len() < min_count
                        || frag_map.get(&(w.0, w.1)).unwrap().len() < min_count)
                    || v.3 != w.2
                {
                    vec![None]
                } else {
                    vec![
                        Some((
                            sid,
                            ShmmrGraphNode(v.0, v.1, v.4),
                            ShmmrGraphNode(w.0, w.1, w.4),
                        )),
                        Some((
                            sid,
                            ShmmrGraphNode(w.0, w.1, 1 - w.4),
                            ShmmrGraphNode(v.0, v.1, 1 - v.4),
                        )),
                    ]
                }
            })
            .flatten()
            .collect::<AdjList>()
    }
}

type PBundleNode = (
    // node, Option<previous_node>, node_weight, is_leaf, global_rank, branch, branch_rank
    ShmmrGraphNode,
    Option<ShmmrGraphNode>,
    u32,
    bool,
    u32,
    u32,
    u32,
);

pub fn sort_adj_list_by_weighted_dfs(
    frag_map: &ShmmrToFrags,
    adj_list: &[AdjPair],
    start: ShmmrGraphNode,
) -> Vec<PBundleNode> {
    use crate::graph_utils::BiDiGraphWeightedDfs;

    let mut g = DiGraphMap::<ShmmrGraphNode, ()>::new();
    let mut score = FxHashMap::<ShmmrGraphNode, u32>::default();
    adj_list.iter().for_each(|&(_sid, v, w)| {
        let vv = (v.0, v.1);
        let ww = (w.0, w.1);
        let v = ShmmrGraphNode(v.0, v.1, v.2);
        let w = ShmmrGraphNode(w.0, w.1, w.2);
        g.add_edge(v, w, ());

        // println!("DBG: add_edge {:?} {:?}", v, w);
        score
            .entry(v)
            .or_insert_with(|| frag_map.get(&vv).unwrap().len() as u32);
        score
            .entry(w)
            .or_insert_with(|| frag_map.get(&ww).unwrap().len() as u32);
    });

    // println!("DBG: # node: {}, # edge: {}", g.node_count(), g.edge_count());

    let start = ShmmrGraphNode(start.0, start.1, start.2);

    let mut weighted_dfs_walker = BiDiGraphWeightedDfs::new(&g, start, &score);
    let mut out = vec![];
    while let Some((node, p_node, is_leaf, rank, branch_id, branch_rank)) =
        weighted_dfs_walker.next(&g)
    {
        let node_count = *score.get(&node).unwrap();
        let p_node = p_node.map(|pnode| ShmmrGraphNode(pnode.0, pnode.1, pnode.2));
        out.push((
            ShmmrGraphNode(node.0, node.1, node.2),
            p_node,
            node_count,
            is_leaf,
            rank,
            branch_id,
            branch_rank,
        ));
        //println!("DBG, next node: {:?}", node);
    }
    out
}

pub fn get_principal_bundles_from_adj_list(
    frag_map: &ShmmrToFrags,
    adj_list: &[AdjPair],
    path_len_cutoff: usize,
) -> (Vec<Vec<ShmmrGraphNode>>, AdjList) {
    assert!(!adj_list.is_empty());
    // println!("DBG: adj_list[0]: {:?}", adj_list[0]);
    let s = adj_list[0].1;
    let sorted_adj_list = sort_adj_list_by_weighted_dfs(frag_map, adj_list, s);

    // println!("DGB: sorted_adj_list len: {}", sorted_adj_list.len());

    let mut paths: Vec<Vec<ShmmrGraphNode>> = vec![];
    let mut path: Vec<ShmmrGraphNode> = vec![];
    for v in sorted_adj_list.into_iter() {
        path.push(v.0);
        if v.3 {
            // it is a leaf node
            paths.push(path.clone());
            path.clear()
        }
    }

    let long_paths = paths
        .into_iter()
        .filter(|p| p.len() > path_len_cutoff as usize);

    let mut main_bundle_path_vertices = FxHashSet::<(u64, u64)>::default();

    long_paths.for_each(|p| {
        p.into_iter().for_each(|v| {
            main_bundle_path_vertices.insert((v.0, v.1));
        })
    });

    let mut g0 = DiGraphMap::<ShmmrGraphNode, ()>::new();
    let mut filtered_adj_list = AdjList::new();
    adj_list.iter().for_each(|&(sid, v, w)| {
        if main_bundle_path_vertices.contains(&(v.0, v.1))
            && main_bundle_path_vertices.contains(&(w.0, w.1))
        {
            g0.add_edge(
                ShmmrGraphNode(v.0, v.1, v.2),
                ShmmrGraphNode(w.0, w.1, w.2),
                (),
            );
            filtered_adj_list.push((sid, v, w));
        }
    });

    let mut g1 = g0.clone();
    let mut terminal_vertices = FxHashSet::<ShmmrGraphNode>::default();

    for (v, w, _) in g0.all_edges() {
        if g0.neighbors_directed(v, Outgoing).count() > 1 {
            terminal_vertices.insert(v);
        };
        if g0.neighbors_directed(w, Incoming).count() > 1 {
            terminal_vertices.insert(v);
        };
    }

    let mut starts = Vec::<ShmmrGraphNode>::default();
    for v in g1.nodes() {
        if g1.neighbors_directed(v, Incoming).count() == 0 {
            starts.push(v);
        }
    }
    // if the whole graph is a loop
    if starts.is_empty() {
        if let Some(v) = g1.nodes().next() {
            starts.push(v);
        }
    };

    let mut principal_bundles = Vec::<Vec<ShmmrGraphNode>>::new();

    while !starts.is_empty() {
        let s = starts.pop().unwrap();
        let mut dfs = Dfs::new(&g1, s);
        let mut path = Vec::<ShmmrGraphNode>::new();
        while let Some(v) = dfs.next(&g1) {
            if terminal_vertices.contains(&v) {
                path.push(v);
                break;
            } else {
                path.push(v);
            }
        }
        if !path.is_empty() {
            path.iter().for_each(|&v| {
                g1.remove_node(v);
                g1.remove_node(ShmmrGraphNode(v.0, v.1, 1 - v.2));
            });

            /*
            let v = path[path.len()-1];

            for w in g1.neighbors_directed(v, Outgoing) {
                if g1.neighbors_directed(w, Incoming).count() == 0 {
                    starts.push(w);
                }
            }
            */
            starts.clear();
            for v in g1.nodes() {
                if g1.neighbors_directed(v, Incoming).count() == 0 {
                    starts.push(v);
                }
            }

            principal_bundles.push(path);
        }

        // if the whole graph is a loop
        if starts.is_empty() {
            if let Some(v) = g1.nodes().next() {
                starts.push(v);
            }
        };
    }
    principal_bundles.sort_by(|a, b| b.len().partial_cmp(&(a.len())).unwrap());
    (principal_bundles, filtered_adj_list)
}

impl CompactSeqDB {
    pub fn generate_smp_adj_list_from_frag_map(
        &self,
        min_count: usize,
        keeps: Option<Vec<u32>>,
    ) -> AdjList {
        frag_map_to_adj_list(&self.frag_map, min_count, keeps)
    }
}

pub type FragmentHit = ((u64, u64), (u32, u32, u8), Vec<FragmentSignature>); // ((hash0, hash1), (pos0, pos1, orientation), fragments)

pub fn raw_query_fragment(
    frag_map: &ShmmrToFrags,
    query_frag: &Vec<u8>,
    shmmr_spec: &ShmmrSpec,
) -> Vec<FragmentHit> {
    let shmmrs = sequence_to_shmmrs(0, query_frag, shmmr_spec, false);
    let query_results = pair_shmmrs(&shmmrs)
        .par_iter()
        .map(|(s0, s1)| {
            let p0 = s0.pos() + 1;
            let p1 = s1.pos() + 1;
            let s0 = s0.hash();
            let s1 = s1.hash();
            if s0 < s1 {
                (s0, s1, p0, p1, 0_u8)
            } else {
                (s1, s0, p0, p1, 1_u8)
            }
        })
        .map(|(s0, s1, p0, p1, orientation)| {
            if let Some(m) = frag_map.get(&(s0, s1)) {
                ((s0, s1), (p0, p1, orientation), m.clone())
            } else {
                ((s0, s1), (p0, p1, orientation), vec![])
            }
        })
        .collect::<Vec<_>>();
    query_results
}

pub fn raw_query_fragment_from_mmap_midx(
    frag_map_location: &ShmmrToIndexFileLocation,
    frag_map_mmap_file: &Mmap,
    query_frag: &Vec<u8>,
    shmmr_spec: &ShmmrSpec,
) -> Vec<FragmentHit> {
    let shmmrs = sequence_to_shmmrs(0, query_frag, shmmr_spec, false);
    let query_results = pair_shmmrs(&shmmrs)
        .par_iter()
        .map(|(s0, s1)| {
            let p0 = s0.pos() + 1;
            let p1 = s1.pos() + 1;
            let s0 = s0.hash();
            let s1 = s1.hash();
            if s0 < s1 {
                (s0, s1, p0, p1, 0_u8)
            } else {
                (s1, s0, p0, p1, 1_u8)
            }
        })
        .map(|(s0, s1, p0, p1, orientation)| {
            if let Some(&(start, vec_len)) = frag_map_location.get(&(s0, s1)) {
                let m = get_fragment_signatures_from_mmap_file(&frag_map_mmap_file, start, vec_len);
                ((s0, s1), (p0, p1, orientation), m)
            } else {
                ((s0, s1), (p0, p1, orientation), vec![])
            }
        })
        .collect::<Vec<_>>();
    query_results
}

pub fn get_match_positions_with_fragment(
    shmmr_map: &ShmmrToFrags,
    frag: &Vec<u8>,
    shmmr_spec: &ShmmrSpec,
) -> FxHashMap<u32, Vec<(u32, u32, u8)>> {
    let mut res = FxHashMap::<u32, Vec<(u32, u32, u8)>>::default();
    raw_query_fragment(shmmr_map, frag, shmmr_spec)
        .into_iter()
        .for_each(|v| {
            let q_direction = v.1 .2;
            v.2.into_iter().for_each(|w| {
                let (_, sid, p0, p1, direction) = w;
                let direction = if direction == q_direction { 0 } else { 1 };
                res.entry(sid).or_default().push((p0, p1, direction));
            });
        });
    res.iter_mut().for_each(|(_k, v)| v.sort());
    res
}

pub fn write_shmmr_map_file(
    shmmr_spec: &ShmmrSpec,
    shmmr_map: &ShmmrToFrags,
    filepath: String,
) -> Result<(), std::io::Error> {
    let mut out_file =
        File::create(filepath).expect("open fail while writing the SHIMMER map (.mdb) file\n");
    let mut buf = Vec::<u8>::new();

    buf.extend("mdb".to_string().into_bytes());

    buf.write_u32::<LittleEndian>(shmmr_spec.w as u32)?;
    buf.write_u32::<LittleEndian>(shmmr_spec.k as u32)?;
    buf.write_u32::<LittleEndian>(shmmr_spec.r as u32)?;
    buf.write_u32::<LittleEndian>(shmmr_spec.min_span as u32)?;
    buf.write_u32::<LittleEndian>(shmmr_spec.sketch as u32)?;

    buf.write_u64::<LittleEndian>(shmmr_map.len() as u64)?;
    shmmr_map
        .iter()
        .try_for_each(|(k, v)| -> Result<(), std::io::Error> {
            buf.write_u64::<LittleEndian>(k.0)?;
            buf.write_u64::<LittleEndian>(k.1)?;
            buf.write_u64::<LittleEndian>(v.len() as u64)?;
            v.iter().try_for_each(|r| -> Result<(), std::io::Error> {
                buf.write_u32::<LittleEndian>(r.0)?;
                buf.write_u32::<LittleEndian>(r.1)?;
                buf.write_u32::<LittleEndian>(r.2)?;
                buf.write_u32::<LittleEndian>(r.3)?;
                buf.write_u8(r.4)?;
                Ok(())
            })
        })?;
    let _ = out_file.write_all(&buf);
    Ok(())
}

pub fn read_mdb_file(filepath: String) -> Result<(ShmmrSpec, ShmmrToFrags), io::Error> {
    let mut in_file =
        File::open(filepath).expect("Error while opening the SHIMMER map file (.mdb) file");
    let mut buf = Vec::<u8>::new();

    let mut u64bytes = [0_u8; 8];
    let mut u32bytes = [0_u8; 4];
    in_file.read_to_end(&mut buf)?;
    let mut cursor = 0_usize;
    assert!(buf[0..3] == "mdb".to_string().into_bytes());
    cursor += 3; // skip "mdb"

    let w = LittleEndian::read_u32(&buf[cursor..cursor + 4]);
    cursor += 4;
    let k = LittleEndian::read_u32(&buf[cursor..cursor + 4]);
    cursor += 4;
    let r = LittleEndian::read_u32(&buf[cursor..cursor + 4]);
    cursor += 4;
    let min_span = LittleEndian::read_u32(&buf[cursor..cursor + 4]);
    cursor += 4;
    let flag = LittleEndian::read_u32(&buf[cursor..cursor + 4]);
    cursor += 4;
    let sketch = (flag & 0b01) == 0b01;

    let shmmr_spec = ShmmrSpec {
        w,
        k,
        r,
        min_span,
        sketch,
    };
    u64bytes.clone_from_slice(&buf[cursor..cursor + 8]);
    let shmmr_key_len = usize::from_le_bytes(u64bytes);
    cursor += 8;
    let mut shmmr_map = ShmmrToFrags::default();
    (0..shmmr_key_len).into_iter().for_each(|_| {
        u64bytes.clone_from_slice(&buf[cursor..cursor + 8]);
        let k1 = u64::from_le_bytes(u64bytes);
        cursor += 8;

        u64bytes.clone_from_slice(&buf[cursor..cursor + 8]);
        let k2 = u64::from_le_bytes(u64bytes);
        cursor += 8;

        u64bytes.clone_from_slice(&buf[cursor..cursor + 8]);
        let vec_len = usize::from_le_bytes(u64bytes);
        cursor += 8;

        let value = (0..vec_len)
            .into_iter()
            .map(|_| {
                let mut v = (0_u32, 0_u32, 0_u32, 0_u32, 0_u8);

                u32bytes.clone_from_slice(&buf[cursor..cursor + 4]);
                v.0 = u32::from_le_bytes(u32bytes);
                cursor += 4;

                u32bytes.clone_from_slice(&buf[cursor..cursor + 4]);
                v.1 = u32::from_le_bytes(u32bytes);
                cursor += 4;

                u32bytes.clone_from_slice(&buf[cursor..cursor + 4]);
                v.2 = u32::from_le_bytes(u32bytes);
                cursor += 4;

                u32bytes.clone_from_slice(&buf[cursor..cursor + 4]);
                v.3 = u32::from_le_bytes(u32bytes);
                cursor += 4;

                v.4 = buf[cursor..cursor + 1][0];
                cursor += 1;

                v
            })
            .collect::<Vec<FragmentSignature>>();

        shmmr_map.insert((k1, k2), value);
    });

    Ok((shmmr_spec, shmmr_map))
}

pub fn read_mdb_file_to_frag_locations(
    filepath: String,
) -> Result<(ShmmrSpec, ShmmrIndexFileLocation), io::Error> {
    let mut in_file =
        File::open(filepath).expect("open fail while reading the SHIMMER map (.mdb) file");
    let mut tag_buf = [0_u8; 3];

    let mut u32bytes = [0_u8; 4];
    let mut u64bytes = [0_u8; 8];

    in_file.read_exact(&mut tag_buf)?;
    let mut cursor = 0_usize;
    assert!(tag_buf[0..3] == "mdb".to_string().into_bytes());
    cursor += 3; // skip "mdb"

    in_file.read_exact(&mut u32bytes)?;
    let w = LittleEndian::read_u32(&u32bytes);

    in_file.read_exact(&mut u32bytes)?;
    let k = LittleEndian::read_u32(&u32bytes);

    in_file.read_exact(&mut u32bytes)?;
    let r = LittleEndian::read_u32(&u32bytes);

    in_file.read_exact(&mut u32bytes)?;
    let min_span = LittleEndian::read_u32(&u32bytes);

    in_file.read_exact(&mut u32bytes)?;
    let flag = LittleEndian::read_u32(&u32bytes);
    let sketch = (flag & 0b01) == 0b01;

    cursor += 4 * 5;

    let shmmr_spec = ShmmrSpec {
        w,
        k,
        r,
        min_span,
        sketch,
    };

    in_file.read_exact(&mut u64bytes)?;
    let shmmr_key_len = usize::from_le_bytes(u64bytes);
    cursor += 8;
    let mut rec_loc = Vec::<((u64, u64), (usize, usize))>::new();
    for _ in 0..shmmr_key_len {
        in_file.read_exact(&mut u64bytes)?;
        let k1 = u64::from_le_bytes(u64bytes);

        in_file.read_exact(&mut u64bytes)?;
        let k2 = u64::from_le_bytes(u64bytes);

        in_file.read_exact(&mut u64bytes)?;
        let vec_len = usize::from_le_bytes(u64bytes);
        cursor += 8 * 3;
        let start = cursor;
        let advance = 17 * vec_len;
        cursor += advance;
        in_file.seek(SeekFrom::Current(advance as i64))?;
        rec_loc.push(((k1, k2), (start, vec_len)));
    }
    Ok((shmmr_spec, rec_loc))
}

pub fn get_fragment_signatures_from_mmap_file(
    frag_map_file: &Mmap,
    start: usize,
    vec_len: usize,
) -> Vec<FragmentSignature> {
    let mut cursor = start;
    (0..vec_len)
        .into_iter()
        .map(|_| {
            let mut u32bytes = [0_u8; 4];
            let mut v = (0_u32, 0_u32, 0_u32, 0_u32, 0_u8);
            u32bytes.clone_from_slice(&frag_map_file[cursor..cursor + 4]);
            v.0 = u32::from_le_bytes(u32bytes);
            cursor += 4;

            u32bytes.clone_from_slice(&frag_map_file[cursor..cursor + 4]);
            v.1 = u32::from_le_bytes(u32bytes);
            cursor += 4;

            u32bytes.clone_from_slice(&frag_map_file[cursor..cursor + 4]);
            v.2 = u32::from_le_bytes(u32bytes);
            cursor += 4;

            u32bytes.clone_from_slice(&frag_map_file[cursor..cursor + 4]);
            v.3 = u32::from_le_bytes(u32bytes);
            cursor += 4;

            v.4 = frag_map_file[cursor..cursor + 1][0];
            cursor += 1;
            v
        })
        .collect::<Vec<FragmentSignature>>()
}

pub fn read_mdb_file_parallel(filepath: String) -> Result<(ShmmrSpec, ShmmrToFrags), io::Error> {
    let in_file =
        File::open(filepath.clone()).expect("open fail while reading the SHIMMER map (.mdb) file");
    let frag_map_file = unsafe {
        Mmap::map(&in_file).expect("open fail while reading the SHIMMER map (.mdb) file")
    };

    let (shmmr_spec, rec_loc) = read_mdb_file_to_frag_locations(filepath)?;

    let shmmr_map = rec_loc
        .par_iter()
        .map(|&((k1, k2), (start, vec_len))| {
            let value = get_fragment_signatures_from_mmap_file(&frag_map_file, start, vec_len);
            ((k1, k2), value)
        })
        .collect::<FxHashMap<ShmmrPair, Vec<FragmentSignature>>>();
    Ok((shmmr_spec, shmmr_map))
}
