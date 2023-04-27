use crate::seq_db::{
    self, read_mdb_file_to_frag_locations, CompactSeq, Fragment, Fragments, GetSeq,
};
use crate::shmmrutils::ShmmrSpec;
use bincode::config;
use flate2::read::DeflateDecoder;
use memmap2::Mmap;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
pub type ShmmrToFragMapLocation = FxHashMap<(u64, u64), (usize, usize)>;

pub struct CompactSeqFragFileStorage {
    pub shmmr_spec: ShmmrSpec,
    pub seqs: Vec<CompactSeq>,
    pub frag_location_map: ShmmrToFragMapLocation,
    pub frag_map_file: Mmap,
    pub frag_file_prefix: String,
    pub frag_file: Mmap,
    pub frag_addr_offsets: Vec<(usize, usize, u32)>, //offset, compress_chunk_size, frag_len_in_bases
    pub frag_compress_chunk_size: usize,
    pub seq_index: FxHashMap<(String, Option<String>), (u32, u32)>,
    /// a dictionary maps id -> (ctg_name, source, len)
    pub seq_info: FxHashMap<u32, (String, Option<String>, u32)>,
}

impl CompactSeqFragFileStorage {
    pub fn new(prefix: String) -> Self {
        let frag_file_prefix = prefix;

        let fmap_file =
            File::open(frag_file_prefix.clone() + ".mdb").expect("frag map file open fail");

        let frag_map_file =
            unsafe { Mmap::map(&fmap_file).expect("frag map file memory map creation fail") };

        let (shmmr_spec, frag_location_map) =
            read_mdb_file_to_frag_locations(frag_file_prefix.clone() + ".mdb").unwrap();

        let frag_location_map =
            FxHashMap::<(u64, u64), (usize, usize)>::from_iter(frag_location_map);

        let mut sdx_file = BufReader::new(
            File::open(frag_file_prefix.clone() + ".sdx").expect("sdx file open error"),
        );
        let mut sdx_version_string = [0_u8;7];
        sdx_file.read_exact(&mut sdx_version_string).expect("sdx file reading error");
        let config = config::standard();
        let (frag_compress_chunk_size, frag_addr_offsets, seqs): (
            usize,
            Vec<(usize, usize, u32)>,
            Vec<CompactSeq>,
        ) = bincode::decode_from_std_read(&mut sdx_file, config).expect("read sdx file error");
        let f_file = File::open(frag_file_prefix.clone() + ".frg").expect("frag file open fail");
        let frag_file = unsafe { Mmap::map(&f_file).expect("frag file memory map creation fail") };

        let mut seq_index = FxHashMap::<(String, Option<String>), (u32, u32)>::default();
        let mut seq_info = FxHashMap::<u32, (String, Option<String>, u32)>::default();

        let midx_file = BufReader::new(
            File::open(frag_file_prefix.clone() + ".midx").expect("open midx file fail"),
        );
        midx_file
            .lines()
            .into_iter()
            .try_for_each(|line| -> Result<(), std::io::Error> {
                let line = line.unwrap();
                let mut line = line.as_str().split('\t');
                let sid = line.next().unwrap().parse::<u32>().unwrap();
                let len = line.next().unwrap().parse::<u32>().unwrap();
                let ctg_name = line.next().unwrap().to_string();
                let source = line.next().unwrap().to_string();
                seq_index.insert((ctg_name.clone(), Some(source.clone())), (sid, len));
                seq_info.insert(sid, (ctg_name, Some(source), len));
                Ok(())
            })
            .expect("read midx file fail");

        Self {
            shmmr_spec,
            seqs,
            frag_location_map,
            frag_map_file,
            frag_file_prefix,
            frag_file,
            frag_addr_offsets,
            frag_compress_chunk_size,
            seq_index,
            seq_info,
        }
    }

    fn reconstruct_sequence_from_frags(&self, frags: Fragments) -> Vec<u8> {
        let mut reconstructed_seq = <Vec<u8>>::new();
        let sub_seqs = frags
            .chunks(32)
            .collect::<Vec<&[Fragment]>>()
            .par_iter()
            .flat_map(|&frags| {
                let mut frag_group_cache = FxHashMap::<u32, Fragments>::default();
                frags.into_iter().map(|frag| {
                    let mut reconstructed_seq = <Vec<u8>>::new();
                    let mut _p = 0;
                    match frag {
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
                        Fragment::AlnSegments((frag_id, reversed, _length, a)) => {
                            let frag_group_id = *frag_id / self.frag_compress_chunk_size as u32;
                            let frag_group =
                                frag_group_cache.entry(frag_group_id).or_insert_with(|| {
                                    fetch_frag_group(
                                        frag_group_id,
                                        &self.frag_addr_offsets,
                                        &self.frag_file,
                                    )
                                });

                            if let Fragment::Internal(base_seq) = frag_group
                                [*frag_id as usize % self.frag_compress_chunk_size]
                                .clone()
                            {
                                let mut seq = seq_db::reconstruct_seq_from_aln_segs(&base_seq, &a);
                                if *reversed {
                                    seq = crate::fasta_io::reverse_complement(&seq);
                                }
                                reconstructed_seq
                                    .extend_from_slice(&seq[self.shmmr_spec.k as usize..]);
                                //println!("p: {} {}", p, p + seq.len());
                                _p += seq.len();
                            }
                        }
                        
                    }
                    reconstructed_seq
                }).collect::<Vec<Vec<u8>>>()
            })
            .collect::<Vec<Vec<u8>>>();
        sub_seqs
            .into_iter()
            .for_each(|s| reconstructed_seq.extend(s));
        reconstructed_seq
    }

    fn get_seq_from_frag_ids<I: Iterator<Item = u32>>(&self, frag_ids: I) -> Vec<u8> {
        let mut _p = 0;
        let mut frag_group_cache = FxHashMap::<u32, Fragments>::default();
        let frags = frag_ids
            .map(|frag_id| {
                let frag_group_id = frag_id / self.frag_compress_chunk_size as u32;
                let frag_group = frag_group_cache.entry(frag_group_id).or_insert_with(|| {
                    fetch_frag_group(frag_group_id, &self.frag_addr_offsets, &self.frag_file)
                });

                let frag = frag_group[frag_id as usize % self.frag_compress_chunk_size].clone();
                frag
            })
            .collect::<Fragments>();

        self.reconstruct_sequence_from_frags(frags)
    }
}

impl GetSeq for CompactSeqFragFileStorage {
    fn get_seq_by_id(&self, sid: u32) -> Vec<u8> {
        assert!((sid as usize) < self.seqs.len());
        let frag_range = &self.seqs[sid as usize].seq_frag_range;
        self.get_seq_from_frag_ids(frag_range.0..frag_range.0 + frag_range.1)
    }

    fn get_sub_seq_by_id(&self, sid: u32, bgn: u32, end: u32) -> Vec<u8> {
        assert!((sid as usize) < self.seqs.len());
        // get these fragment of the first group
        let frag_range = self.seqs[sid as usize].seq_frag_range;
        let frag_range = (frag_range.0, frag_range.0 + frag_range.1); // original frag_range.0:start, frag_range.1: length
        let frag_group_ids: Vec<(u32, u32)> = (frag_range.0..frag_range.1)
            .map(|v| (v / self.frag_compress_chunk_size as u32, v))
            .collect();
        let first_group_ids = frag_group_ids
            .iter()
            .filter(|(gid, _id)| *gid == frag_group_ids[0].0)
            .map(|v| v.1);

        let first_group_seq = self.get_seq_from_frag_ids(first_group_ids.clone());

        let mut current_chunk_bgn;
        let mut current_chunk_end = first_group_seq.len() as u32;
        let mut sub_seqs = Vec::<(u32, Vec<u8>)>::new();
        if bgn < current_chunk_end as u32 {
            sub_seqs.push((0, first_group_seq));
        };

        let group_ids = (frag_group_ids[0].0..=frag_group_ids[frag_group_ids.len() - 1].0)
            .map(|v| v)
            .collect::<Vec<u32>>();

        if group_ids.len() > 1 {
            for &group_id in group_ids[1..].into_iter()  {
                let (_, _, frag_seq_len) = self.frag_addr_offsets[group_id as usize];
                current_chunk_bgn = current_chunk_end;
                current_chunk_end = current_chunk_bgn + frag_seq_len as u32;
                if (current_chunk_bgn <= bgn && bgn < current_chunk_end)
                    || (current_chunk_bgn <= end && end < current_chunk_end)
                    || (bgn <= current_chunk_bgn && current_chunk_end <= end)
                {
                    let frags =
                        fetch_frag_group(group_id, &self.frag_addr_offsets, &self.frag_file);
                    let sub_seq = self.reconstruct_sequence_from_frags(frags);
                    sub_seqs.push((current_chunk_bgn, sub_seq));
                }
            }
        }
        let mut seq = Vec::<u8>::new();
        println!("{:?} {} {} {} {} {}", frag_range, sid, group_ids.len(), bgn, end, end-bgn);
        let offset = (bgn - sub_seqs[0].0) as usize;
        sub_seqs.into_iter().for_each(|ss| seq.extend(ss.1));
        return seq[offset..offset + (end - bgn) as usize].to_vec();
    }
}

fn fetch_frag_group(
    frag_group_id: u32,
    frag_addr_offsets: &[(usize, usize, u32)],
    frag_file: &Mmap,
) -> Fragments {
    let config = config::standard();
    let (offset, size, _) = frag_addr_offsets[frag_group_id as usize];
    let version_string_offset = 7;
    let offset = offset + version_string_offset; 
    let compress_chunk = frag_file[offset..(offset + size)].to_vec();
    let mut deflater = DeflateDecoder::new(&compress_chunk[..]);
    let mut s: Vec<u8> = vec![];
    deflater.read_to_end(&mut s).expect("decompression error");
    let (frags, _size): (Fragments, usize) =
        bincode::decode_from_slice::<Fragments, bincode::config::Configuration>(&s[..], config)
            .unwrap();
    frags
}
