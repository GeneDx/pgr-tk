use crate::seq_db::{self, read_mdb_file_parallel, CompactSeq, Fragment, ShmmrToFrags, GetSeq, read_mdb_file_to_frag_locations};
use crate::shmmrutils::ShmmrSpec;
use bincode::config;
use flate2::read::DeflateDecoder;
use memmap2::Mmap;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

type ShmmrToFragMapLocation = FxHashMap<(u64, u64), (usize, usize)>;

pub struct CompactSeqDBStorage {
    pub shmmr_spec: ShmmrSpec,
    pub seqs: Vec<CompactSeq>,
    pub frag_location_map: ShmmrToFragMapLocation,
    pub frag_map_file: Mmap,
    pub frag_map: ShmmrToFrags,
    pub frag_file_prefix: String,
    pub frag_file: Mmap,
    pub frag_addr_offsets: Vec<(usize, usize, u32)>, //offset, compress_chunk_size, frag_len_in_bases
    pub seq_index: FxHashMap<(String, Option<String>), (u32, u32)>,
    /// a dictionary maps id -> (ctg_name, source, len)
    pub seq_info: FxHashMap<u32, (String, Option<String>, u32)>,
}

impl CompactSeqDBStorage {
    pub fn new(prefix: String) -> Self {
        let frag_file_prefix = prefix;

        let fmap_file = File::open(frag_file_prefix.clone() + ".mdb").expect("frag map file open fail");
        let frag_map_file = unsafe { Mmap::map(&fmap_file).expect("frag map file memory map creation fail") };

        let (shmmr_spec, frag_map) =
            read_mdb_file_parallel(frag_file_prefix.clone() + ".mdb").unwrap();

        let (_shmmr_spec, frag_location_map) =
            read_mdb_file_to_frag_locations(frag_file_prefix.clone() + ".mdb").unwrap();

        let frag_location_map = FxHashMap::<(u64, u64),(usize, usize)>::from_iter(frag_location_map);

        let mut sdx_file = BufReader::new(
            File::open(frag_file_prefix.clone() + ".sdx").expect("sdx file open error"),
        );
        let config = config::standard();
        let (frag_addr_offsets, seqs): (Vec<(usize, usize, u32)>, Vec<CompactSeq>) =
            bincode::decode_from_std_read(&mut sdx_file, config).expect("read sdx file error");
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
            frag_map,
            frag_file_prefix,
            frag_file,
            frag_addr_offsets,
            seq_index,
            seq_info,
        }
    }



    fn get_seq_from_frag_ids<I: Iterator<Item = u32>>(&self, frag_ids: I) -> Vec<u8> {
        let mut reconstructed_seq = <Vec<u8>>::new();

        let mut _p = 0;
        frag_ids.for_each(|frag_id| {
            let frag = fetch_frag(frag_id, &self.frag_addr_offsets, &self.frag_file);
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
                    if let Fragment::Internal(base_seq) =
                        fetch_frag(frag_id, &self.frag_addr_offsets, &self.frag_file)
                    {
                        let mut seq = seq_db::reconstruct_seq_from_aln_segs(&base_seq, &a);
                        if reversed {
                            seq = crate::fasta_io::reverse_complement(&seq);
                        }
                        reconstructed_seq.extend_from_slice(&seq[self.shmmr_spec.k as usize..]);
                        //println!("p: {} {}", p, p + seq.len());
                        _p += seq.len();
                    }
                }
            }
        });

        reconstructed_seq
    }
}

impl GetSeq for CompactSeqDBStorage {
    fn get_seq_by_id(&self, sid: u32) -> Vec<u8> {
        assert!((sid as usize) < self.seqs.len());
        let frag_range = &self.seqs[sid as usize].seq_frag_range;
        self.get_seq_from_frag_ids(frag_range.0..frag_range.0 + frag_range.1)
    }

    fn get_sub_seq_by_id(&self, sid: u32, bgn: u32, end: u32) -> Vec<u8> {
        assert!((sid as usize) < self.seqs.len());
        let frag_range = &self.seqs[sid as usize].seq_frag_range;

        let mut _p = 0;
        let mut base_offset = 0_u32;
        let mut sub_seq_frag = vec![];
        for frag_id in frag_range.0..frag_range.0 + frag_range.1 {
            let (_, _, mut frag_len) = self.frag_addr_offsets[frag_id as usize];
            if frag_id != frag_range.0 && frag_id != frag_range.0 + frag_range.1 - 1 {
                // for internal segments
                frag_len -= self.shmmr_spec.k;
            }
            if (base_offset <= bgn && bgn < base_offset + frag_len)
                || (base_offset <= end && end < base_offset + frag_len)
                || (bgn <= base_offset && base_offset + frag_len <= end)
            {
                sub_seq_frag.push((frag_id, base_offset));
            }

            base_offset += frag_len;
        }

        let reconstructed_seq = self.get_seq_from_frag_ids(sub_seq_frag.iter().map(|v| v.0));

        let offset = bgn - sub_seq_frag[0].1;

        reconstructed_seq[(offset as usize)..((offset + end - bgn) as usize)].to_vec()
    }
    
}

fn fetch_frag(
    frag_id: u32,
    frag_addr_offsets: &[(usize, usize, u32)],
    frag_file: &Mmap,
) -> Fragment {
    let config = config::standard();
    let (offset, size, _) = frag_addr_offsets[frag_id as usize];
    let compress_chunk = frag_file[offset..(offset + size as usize)].to_vec();
    let mut deflater = DeflateDecoder::new(&compress_chunk[..]);
    let mut s: Vec<u8> = vec![];
    deflater.read_to_end(&mut s).expect("decompression error");
    let (frag, _size): (Fragment, usize) =
        bincode::decode_from_slice::<Fragment, bincode::config::Configuration>(&s[..], config)
            .unwrap();
    frag
}
