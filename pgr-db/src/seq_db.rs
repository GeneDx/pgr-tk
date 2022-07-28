use crate::agc_io::AGCFile;
use crate::fasta_io::{reverse_complement, FastaReader, SeqRec};
use crate::graph_utils::ShmmrGraphNode;
use crate::shmmrutils::{match_reads, sequence_to_shmmrs, DeltaPoint, ShmmrSpec, MM128};
use byteorder::{ByteOrder, LittleEndian, WriteBytesExt};
use flate2::bufread::MultiGzDecoder;
use petgraph::graphmap::DiGraphMap;
use petgraph::visit::Dfs;
use petgraph::EdgeDirection::{Incoming, Outgoing};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fmt;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};

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

#[derive(Debug, Clone)]
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
#[derive(Debug, Clone)]
pub struct CompactSeq {
    pub source: Option<String>,
    pub name: String,
    pub id: u32,
    pub seq_frags: Vec<u32>,
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

        //assert!(shmmrs.len() > 0);
        if shmmrs.len() == 0 {
            let frg = seq[..].to_vec();
            frags.push(Fragment::Prefix(frg));
            seq_frags.push(frg_id);
            frg_id += 1;

            let frg = Vec::<u8>::new();
            frags.push(Fragment::Suffix(frg));
            seq_frags.push(frg_id);

            return CompactSeq {
                source,
                name,
                id,
                seq_frags,
                len: seq.len(),
            };
        }
        // prefix
        let end = (shmmrs[0].pos() + 1) as usize;
        let frg = seq[..end].to_vec();
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

                if frg_len > 64 && try_compress && self.frag_map.contains_key(&shmmr_pair) {
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
                                    &seq[(bgn - self.shmmr_spec.k) as usize..end as usize].to_vec(),
                                );
                                rc = true;
                            } else {
                                frg =
                                    seq[(bgn - self.shmmr_spec.k) as usize..end as usize].to_vec();
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
                    let frg = seq[(bgn - self.shmmr_spec.k) as usize..end as usize].to_vec();
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
                frags.push(frg.clone());
                seq_frags.push(frg_id);
                frg_id += 1;
            }
            None => {}
        });

        // suffix
        let bgn = (shmmrs[shmmrs.len() - 1].pos() + 1) as usize;
        let frg = seq[bgn..].to_vec();
        frags.push(Fragment::Suffix(frg));
        seq_frags.push(frg_id);

        CompactSeq {
            source,
            name,
            id,
            seq_frags,
            len: seq.len(),
        }
    }

    pub fn _seq_to_index(
        &mut self,
        source: Option<String>,
        name: String,
        id: u32,
        seqlen: usize,
        shmmrs: Vec<MM128>,
    ) -> CompactSeq {
        assert!(self.frags.is_none());
        let mut seq_frags = Vec::<u32>::new();
        let mut frg_id = 0;

        //assert!(shmmrs.len() > 0);
        if shmmrs.len() == 0 {
            return CompactSeq {
                source,
                name,
                id,
                seq_frags,
                len: seqlen,
            };
        }

        seq_frags.push(frg_id);
        frg_id += 1;

        let shmmr_pairs = shmmrs[0..shmmrs.len() - 1]
            .iter()
            .zip(shmmrs[1..shmmrs.len()].iter())
            .collect::<Vec<_>>();

        let internal_frags = shmmr_pairs
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

        internal_frags
            .iter()
            .for_each(|(shmmr, bgn, end, orientation)| {
                let e = self.frag_map.entry(*shmmr).or_insert(vec![]);
                e.push((frg_id, id, *bgn, *end, *orientation));
                seq_frags.push(frg_id);
                frg_id += 1;
            });

        seq_frags.push(frg_id);

        CompactSeq {
            source,
            name,
            id,
            seq_frags,
            len: seqlen,
        }
    }

    pub fn seq_to_index(
        source: Option<String>,
        name: String,
        id: u32,
        seqlen: usize,
        shmmrs: Vec<MM128>,
    ) -> (CompactSeq, Vec<((u64, u64), u32, u32, u8)>) {
        //assert!(shmmrs.len() > 0);
        if shmmrs.len() == 0 {
            return (
                CompactSeq {
                    source,
                    name,
                    id,
                    seq_frags: vec![],
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

        (
            CompactSeq {
                source,
                name,
                id,
                seq_frags,
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
        let all_shmmers = seqs
            .par_iter()
            .map(|(sid, _, _, seq)| {
                let shmmrs = sequence_to_shmmrs(*sid, &seq, &self.shmmr_spec, false);
                //let shmmrs = sequence_to_shmmrs2(*sid, &seq, 80, KMERSIZE, 4);
                (*sid, shmmrs)
            })
            .collect::<Vec<(u32, Vec<MM128>)>>();
        all_shmmers
    }

    fn load_seq_from_reader(&mut self, reader: &mut dyn Iterator<Item = io::Result<SeqRec>>) {
        let mut seqs = <Vec<(u32, Option<String>, String, Vec<u8>)>>::new();
        let mut sid = 0;
        self.frags = Some(Fragments::new());
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
        ();
    }

    pub fn load_seqs_from_seq_vec(&mut self, seqs: &Vec<(u32, Option<String>, String, Vec<u8>)>) {
        if self.frags.is_none() {
            self.frags = Some(Fragments::new());
        }
        let all_shmmers = self.get_shmmrs_from_seqs(seqs);
        seqs.iter()
            .zip(all_shmmers)
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
            GZFastaReader::GZFile(reader) => self.load_seq_from_reader(&mut reader.into_iter()),

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
        ();
    }

    pub fn load_index_from_seq_vec(&mut self, seqs: &Vec<(u32, Option<String>, String, Vec<u8>)>) {
        let all_shmmers = self.get_shmmrs_from_seqs(seqs);
        let seq_names = seqs
            .iter()
            .map(|(_sid, src, n, s)| (src.clone(), n.clone(), s.len()))
            .collect::<Vec<(Option<String>, String, usize)>>();

        /*
        seq_names.iter().zip(all_shmmers).for_each(
            |((source, seq_name, seqlen), (sid, shmmrs))| {
                let compress_seq =
                    self._seq_to_index(source.clone(), seq_name.clone(), sid, *seqlen, shmmrs);
                self.seqs.push(compress_seq);
            },
        );
        */

        seq_names
            .par_iter()
            .zip(all_shmmers)
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
                internal_frags.iter().zip(cs.seq_frags.clone()).for_each(
                    |((shmmr, bgn, end, orientation), frg_id)| {
                        let e = self.frag_map.entry(*shmmr).or_insert(vec![]);
                        e.push((frg_id, sid, *bgn, *end, *orientation));
                    },
                );
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

impl CompactSeqDB {
    pub fn get_seq(&self, seq: &CompactSeq) -> Vec<u8> {
        let mut reconstructed_seq = <Vec<u8>>::new();
        let frags: &Vec<Fragment> = self.frags.as_ref().unwrap();
        let mut _p = 0;
        for frg_id in seq.seq_frags.iter() {
            //println!("{}:{}", frg_id, sdb.frags[*frg_id as usize]);
            match frags.get(*frg_id as usize).unwrap() {
                Fragment::Prefix(b) => {
                    reconstructed_seq.extend_from_slice(&b[..]);
                    //println!("p: {} {}", p, p + b.len());
                    _p += b.len();
                }
                Fragment::Suffix(b) => {
                    reconstructed_seq.extend_from_slice(&b[..]);
                    //println!("p: {} {}", p, p + b.len());
                    _p += b.len();
                }
                Fragment::Internal(b) => {
                    reconstructed_seq.extend_from_slice(&b[self.shmmr_spec.k as usize..]);
                    //println!("p: {} {}", p, p + b.len());
                    _p += b.len();
                }
                Fragment::AlnSegments((frg_id, reverse, a)) => {
                    if let Fragment::Internal(base_seq) = frags.get(*frg_id as usize).unwrap() {
                        let mut seq = reconstruct_seq_from_aln_segs(&base_seq, a);
                        if *reverse == true {
                            seq = reverse_complement(&seq);
                        }
                        reconstructed_seq.extend_from_slice(&seq[self.shmmr_spec.k as usize..]);
                        //println!("p: {} {}", p, p + seq.len());
                        _p += seq.len();
                    }
                }
            }
        }

        reconstructed_seq
    }

    /* TODO */
    /*
    pub fn get_sub_seq(&self, seq: &CompactSeq, b: usize, e:usize) -> Vec<u8> {
        vec![]
    }
    */

    pub fn get_seq_by_id(&self, sid: u32) -> Vec<u8> {
        let seq = self.seqs.get(sid as usize).unwrap();
        self.get_seq(seq)
    }
}

impl CompactSeqDB {
    pub fn write_shmr_map_index(&self, fp_prefix: String) -> Result<(), std::io::Error> {
        let seq_idx_fp = fp_prefix.clone() + ".midx";
        let data_fp = fp_prefix + ".mdb";
        write_shmr_map_file(&self.shmmr_spec, &self.frag_map, data_fp)?;
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
                    s.source.clone().unwrap_or("-".to_string())
                )?;
                Ok(())
            })?;

        Ok(())
    }
}

pub fn frag_map_to_adj_list(
    frag_map: &ShmmrToFrags,
    min_count: usize,
) -> Vec<(u32, (u64, u64, u8), (u64, u64, u8))> {
    let mut out = frag_map
        .par_iter()
        .flat_map(|v| {
            v.1.iter()
                .map(|vv| (vv.1, vv.2, vv.3, (v.0 .0, v.0 .1, vv.4)))
                .collect::<Vec<(u32, u32, u32, (u64, u64, u8))>>()
        })
        .collect::<Vec<(u32, u32, u32, (u64, u64, u8))>>();
    if out.len() < 2 {
        return vec![];
    }
    out.par_sort();
    let out = out
        .into_par_iter()
        .map(|v| {
            if frag_map.get(&(v.3 .0, v.3 .1)).unwrap().len() >= min_count {
                Some(v)
            } else {
                None
            }
        })
        .collect::<Vec<Option<(u32, u32, u32, (u64, u64, u8))>>>();

    (0..out.len() - 1)
        .into_par_iter()
        .flat_map(|i| {
            let v = out[i];
            let w = out[i + 1];
            if v.is_none() || w.is_none() {
                vec![None]
            } else {
                let v = v.unwrap();
                let w = w.unwrap();
                if v.0 != w.0 {
                    vec![None]
                } else {
                    vec![
                        Some((v.0, v.3, w.3)),
                        Some((
                            v.0,
                            (w.3 .0, w.3 .1, 1 - w.3 .2),
                            (v.3 .0, v.3 .1, 1 - v.3 .2),
                        )),
                    ]
                }
            }
        })
        .filter(|v| v.is_some())
        .map(|v| v.unwrap())
        .collect::<Vec<(u32, (u64, u64, u8), (u64, u64, u8))>>() // seq_id, node0, node1
}

pub fn sort_adj_list_by_weighted_dfs(
    frag_map: &ShmmrToFrags,
    adj_list: &Vec<(u32, (u64, u64, u8), (u64, u64, u8))>,
    start: (u64, u64, u8),
) -> Vec<(
    (u64, u64, u8),
    Option<(u64, u64, u8)>,
    u32,
    bool,
    u32,
    u32,
    u32,
)> {
    // node, node_weight, is_leaf, global_rank, branch, branch_rank
    use crate::graph_utils::BiDiGraphWeightedDfs;

    let mut g = DiGraphMap::<ShmmrGraphNode, ()>::new();
    let mut score = FxHashMap::<ShmmrGraphNode, u32>::default();
    adj_list.into_iter().for_each(|&(_sid, v, w)| {
        let vv = (v.0, v.1);
        let ww = (w.0, w.1);
        let v = ShmmrGraphNode(v.0, v.1, v.2);
        let w = ShmmrGraphNode(w.0, w.1, w.2);
        g.add_edge(v, w, ());

        //println!("DBG: add_edge {:?} {:?}", v, w);
        score
            .entry(v)
            .or_insert_with(|| frag_map.get(&vv).unwrap().len() as u32);
        score
            .entry(w)
            .or_insert_with(|| frag_map.get(&ww).unwrap().len() as u32);
    });

    //println!("DBG: {} {}", g.node_count(), g.edge_count());

    let start = ShmmrGraphNode(start.0, start.1, start.2);

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
            //println!("{:?}", node);
        } else {
            break;
        }
    }
    out
}

pub fn get_principal_bundles_from_adj_list(
    frag_map: &ShmmrToFrags,
    adj_list: &Vec<(u32, (u64, u64, u8), (u64, u64, u8))>,
    path_len_cutoff: usize,
) -> (Vec<Vec<ShmmrGraphNode>>, Vec<(u32, (u64, u64, u8), (u64, u64, u8))>) {
    assert!(adj_list.len() > 0);
    let s = adj_list[0].1;
    let sorted_adj_list = sort_adj_list_by_weighted_dfs(frag_map, adj_list, s);

    let mut paths: Vec<Vec<(u64, u64, u8)>> = vec![];
    let mut path: Vec<(u64, u64, u8)> = vec![];
    for v in sorted_adj_list.into_iter() {
        path.push(v.0);
        if v.3 == true {
            // it is a leaf node
            paths.push(path.clone());
            path.clear()
        }
    }

    let long_paths: Vec<Vec<(u64, u64, u8)>> = paths
        .into_iter()
        .filter(|p| p.len() > path_len_cutoff as usize)
        .collect();

    let mut main_bundle_path_vertices = FxHashSet::<(u64, u64)>::default();

    long_paths.into_iter().for_each(|p| {
        p.into_iter().for_each(|v| {
            main_bundle_path_vertices.insert((v.0, v.1));
        })
    });

    let mut g0 = DiGraphMap::<ShmmrGraphNode, ()>::new();
    let mut filtered_adj_list = Vec::<(u32, (u64, u64, u8), (u64, u64, u8))>::new();
    adj_list.into_iter().for_each(|&(sid, v, w)| {
        if main_bundle_path_vertices.contains(&(v.0, v.1))
            && main_bundle_path_vertices.contains(&(w.0, w.1))
        {
            g0.add_edge(
                ShmmrGraphNode(v.0, v.1, v.2),
                ShmmrGraphNode(w.0, w.1, w.2),
                (),
            );
            filtered_adj_list.push( (sid, v, w) ); 
        }
    });

    let mut g1 = DiGraphMap::<ShmmrGraphNode, ()>::new();
    let mut terminal_vertices = FxHashSet::<ShmmrGraphNode>::default();
    for (v, w, _) in g0.all_edges() {
        let mut remove_edge = false;
        if g0.neighbors_directed(v, Outgoing).count() > 1
            || g0.neighbors_directed(v, Incoming).count() > 1
        {
            terminal_vertices.insert(v);
            remove_edge = true;
        };
        if g0.neighbors_directed(w, Outgoing).count() > 1
            || g0.neighbors_directed(w, Incoming).count() > 1
        {
            terminal_vertices.insert(w);
            remove_edge = true;
        };
        if !remove_edge {
            g1.add_edge(
                ShmmrGraphNode(v.0, v.1, v.2),
                ShmmrGraphNode(w.0, w.1, w.2),
                (),
            );
        }
    }

    let mut starts = Vec::<ShmmrGraphNode>::default();
    for v in g1.nodes() {
        if g1.neighbors_directed(v, Incoming).count() == 0 {
            starts.push(v);
        }
    }

    let mut principal_bundles = Vec::<Vec<ShmmrGraphNode>>::new();

    while starts.len() != 0 {
        let s = starts.pop().unwrap();
        let mut dfs = Dfs::new(&g1, s);
        let mut path = Vec::<ShmmrGraphNode>::new();
        while let Some(v) = dfs.next(&g1) {
            path.push(v);
        }
        if path.len() > 0 {
            path.iter().for_each(|&v| {
                g1.remove_node(v);
                g1.remove_node(ShmmrGraphNode(v.0, v.1, 1 - v.2));
            });
            principal_bundles.push(path);
        }
    }
    principal_bundles.sort_by(|a, b| b.len().partial_cmp(&(a.len())).unwrap());
    (principal_bundles, filtered_adj_list)
}

impl CompactSeqDB {
    pub fn generate_smp_adj_list(
        &self,
        min_count: usize,
    ) -> Vec<(u32, (u64, u64, u8), (u64, u64, u8))> {
        frag_map_to_adj_list(&self.frag_map, min_count)
    }
}

pub fn query_fragment(
    shmmr_map: &ShmmrToFrags,
    frag: &Vec<u8>,
    shmmr_spec: &ShmmrSpec,
) -> Vec<((u64, u64), (u32, u32, u8), Vec<FragmentSignature>)> {
    let shmmrs = sequence_to_shmmrs(0, &frag, &shmmr_spec, false);
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
            if let Some(m) = shmmr_map.get(&(s0, s1)) {
                ((s0, s1), (p0, p1, orientation), m.clone())
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
    query_fragment(shmmr_map, &frag, shmmr_spec)
        .into_iter()
        .for_each(|v| {
            let q_direction = v.1 .2;
            v.2.into_iter().for_each(|w| {
                let (_, sid, p0, p1, direction) = w;
                let direction = if direction == q_direction { 0 } else { 1 };
                res.entry(sid).or_insert(vec![]).push((p0, p1, direction));
            });
        });
    res.iter_mut().for_each(|(_k, v)| v.sort());
    res
}

pub fn write_shmr_map_file(
    shmmr_spec: &ShmmrSpec,
    shmmr_map: &ShmmrToFrags,
    filepath: String,
) -> Result<(), std::io::Error> {
    let mut out_file = File::create(filepath).expect("open fail");
    let mut buf = Vec::<u8>::new();

    buf.extend("mdb".to_string().into_bytes());

    buf.write_u32::<LittleEndian>(shmmr_spec.w as u32)?;
    buf.write_u32::<LittleEndian>(shmmr_spec.k as u32)?;
    buf.write_u32::<LittleEndian>(shmmr_spec.r as u32)?;
    buf.write_u32::<LittleEndian>(shmmr_spec.min_span as u32)?;
    buf.write_u32::<LittleEndian>(if shmmr_spec.sketch { 1 } else { 0 } as u32)?;

    buf.write_u64::<LittleEndian>(shmmr_map.len() as u64)?;
    shmmr_map
        .into_iter()
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
    let mut in_file = File::open(filepath).expect("open fail");
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
            .collect::<Vec<(u32, u32, u32, u32, u8)>>();

        shmmr_map.insert((k1, k2), value);
    });

    Ok((shmmr_spec, shmmr_map))
}

pub fn read_mdb_file_parallel(filepath: String) -> Result<(ShmmrSpec, ShmmrToFrags), io::Error> {
    let mut in_file = File::open(filepath).expect("open fail");
    let mut buf = Vec::<u8>::new();

    let mut u64bytes = [0_u8; 8];

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
    ShmmrToFrags::default();
    let mut rec_loc = Vec::<(u64, u64, usize, usize)>::new();
    for _ in 0..shmmr_key_len {
        u64bytes.clone_from_slice(&buf[cursor..cursor + 8]);
        let k1 = u64::from_le_bytes(u64bytes);
        cursor += 8;

        u64bytes.clone_from_slice(&buf[cursor..cursor + 8]);
        let k2 = u64::from_le_bytes(u64bytes);
        cursor += 8;

        u64bytes.clone_from_slice(&buf[cursor..cursor + 8]);
        let vec_len = usize::from_le_bytes(u64bytes);
        cursor += 8;

        let start = cursor;
        cursor += vec_len * 17;
        rec_loc.push((k1, k2, start, vec_len))
    }

    let shmmr_map = rec_loc
        .par_iter()
        .map(|&(k1, k2, start, vec_len)| {
            let mut cursor = start;
            let value = (0..vec_len)
                .into_iter()
                .map(|_| {
                    let mut u32bytes = [0_u8; 4];
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
                .collect::<Vec<(u32, u32, u32, u32, u8)>>();
            ((k1, k2), value)
        })
        .collect::<FxHashMap<(u64, u64), Vec<(u32, u32, u32, u32, u8)>>>();
    Ok((shmmr_spec, shmmr_map))
}
