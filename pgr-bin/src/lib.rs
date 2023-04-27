use flate2::bufread::MultiGzDecoder;
use memmap2::Mmap;
use pgr_db::aln;
use pgr_db::fasta_io::FastaReader;
use pgr_db::graph_utils::{AdjList, ShmmrGraphNode};
pub use pgr_db::seq_db::pair_shmmrs;
use pgr_db::seq_db::{self, GetSeq, raw_query_fragment, raw_query_fragment_from_mmap_midx};
pub use pgr_db::shmmrutils::{sequence_to_shmmrs, ShmmrSpec};
use pgr_db::{agc_io, frag_file_io};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

#[allow(clippy::large_enum_variant)]
pub enum GZFastaReader {
    GZFile(FastaReader<BufReader<MultiGzDecoder<BufReader<File>>>>),
    RegularFile(FastaReader<BufReader<BufReader<File>>>),
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Backend {
    AGC,
    FRG,
    FASTX,
    MEMORY,
    UNKNOWN,
}

pub struct SeqIndexDB {
    /// Rust internal: store the specification of the shmmr_spec
    pub shmmr_spec: Option<ShmmrSpec>,
    /// Rust internal: store the sequences
    pub seq_db: Option<seq_db::CompactSeqDB>,
    /// Rust internal: store the agc file and the index
    pub agc_db: Option<agc_io::AGCSeqDB>,
    pub frg_db: Option<frag_file_io::CompactSeqFragFileStorage>,
    /// a dictionary maps (ctg_name, source) -> (id, len)
    #[allow(clippy::type_complexity)]
    pub seq_index: Option<FxHashMap<(String, Option<String>), (u32, u32)>>,
    /// a dictionary maps id -> (ctg_name, source, len)
    #[allow(clippy::type_complexity)]
    pub seq_info: Option<FxHashMap<u32, (String, Option<String>, u32)>>,
    pub backend: Backend,
}

impl Default for SeqIndexDB {
    fn default() -> Self {
        Self::new()
    }
}

impl SeqIndexDB {
    pub fn new() -> Self {
        SeqIndexDB {
            seq_db: None,
            frg_db: None,
            agc_db: None,
            shmmr_spec: None,
            seq_index: None,
            seq_info: None,
            backend: Backend::UNKNOWN,
        }
    }

    pub fn load_from_agc_index(&mut self, prefix: String) -> Result<(), std::io::Error> {
        let (shmmr_spec, frag_location_map) =
            seq_db::read_mdb_file_to_frag_locations(prefix.to_string() + ".mdb").unwrap();

        let frag_location_map =
            FxHashMap::<(u64, u64), (usize, usize)>::from_iter(frag_location_map);

        let agc_file = agc_io::AGCFile::new(prefix.to_string() + ".agc")?;

        let fmap_file =
            File::open(prefix.clone() + ".mdb").expect("frag map file open fail");
        let frag_map_file =
            unsafe { Mmap::map(&fmap_file).expect("frag map file memory map creation fail") };

        self.agc_db = Some(agc_io::AGCSeqDB {agc_file, frag_location_map, frag_map_file });
        self.backend = Backend::AGC;
        self.shmmr_spec = Some(shmmr_spec);

        let mut seq_index = FxHashMap::<(String, Option<String>), (u32, u32)>::default();
        let mut seq_info = FxHashMap::<u32, (String, Option<String>, u32)>::default();

        let midx_file = BufReader::new(File::open(prefix + ".midx")?);
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
            })?;

        self.seq_index = Some(seq_index);
        self.seq_info = Some(seq_info);
        Ok(())
    }

    pub fn load_from_frg_index(&mut self, prefix: String) -> Result<(), std::io::Error> {
        let mut frag_db = pgr_db::frag_file_io::CompactSeqFragFileStorage::new(prefix);

        let seq_index = frag_db.seq_index.into_iter().map(|(k, v)| (k, v)).collect();

        let seq_info = frag_db.seq_info.into_iter().map(|(k, v)| (k, v)).collect();

        frag_db.seq_index = FxHashMap::<(String, Option<String>), (u32, u32)>::default();
        frag_db.seq_info = FxHashMap::<u32, (String, Option<String>, u32)>::default();

        let shmmr_spec = frag_db.shmmr_spec.clone();

        self.frg_db = Some(frag_db);
        self.backend = Backend::FRG;
        self.shmmr_spec = Some(shmmr_spec);

        self.seq_index = Some(seq_index);
        self.seq_info = Some(seq_info);
        Ok(())
    }

    pub fn load_from_fastx(
        &mut self,
        filepath: String,
        w: u32,
        k: u32,
        r: u32,
        min_span: u32,
    ) -> Result<(), std::io::Error> {
        let spec = ShmmrSpec {
            w,
            k,
            r,
            min_span,
            sketch: false,
        };
        let mut sdb = seq_db::CompactSeqDB::new(spec.clone());
        sdb.load_seqs_from_fastx(filepath)?;
        self.shmmr_spec = Some(spec);
        let mut seq_index = FxHashMap::<(String, Option<String>), (u32, u32)>::default();
        let mut seq_info = FxHashMap::<u32, (String, Option<String>, u32)>::default();
        sdb.seqs.iter().for_each(|v| {
            seq_index.insert((v.name.clone(), v.source.clone()), (v.id, v.len as u32));
            seq_info.insert(v.id, (v.name.clone(), v.source.clone(), v.len as u32));
        });
        self.seq_index = Some(seq_index);
        self.seq_info = Some(seq_info);
        self.seq_db = Some(sdb);
        self.backend = Backend::FASTX;
        Ok(())
    }

    pub fn append_from_fastx(&mut self, filepath: String) -> Result<(), std::io::Error> {
        assert!(
            self.backend == Backend::FASTX,
            "Only DB created with load_from_fastx() can add data from another fastx file"
        );
        let sdb = self.seq_db.as_mut().unwrap();
        sdb.load_seqs_from_fastx(filepath)?;
        let mut seq_index = FxHashMap::<(String, Option<String>), (u32, u32)>::default();
        let mut seq_info = FxHashMap::<u32, (String, Option<String>, u32)>::default();
        sdb.seqs.iter().for_each(|v| {
            seq_index.insert((v.name.clone(), v.source.clone()), (v.id, v.len as u32));
            seq_info.insert(v.id, (v.name.clone(), v.source.clone(), v.len as u32));
        });
        self.seq_index = Some(seq_index);
        self.seq_info = Some(seq_info);
        Ok(())
    }

    pub fn write_frag_and_index_files(&self, file_prefix: String) {
        if self.seq_db.is_some() {
            let internal = self.seq_db.as_ref().unwrap();

            internal.write_to_frag_files(file_prefix.clone(), None);
            internal
                .write_shmmr_map_index(file_prefix)
                .expect("write mdb file fail");
        };
    }

    pub fn load_from_seq_list(
        &mut self,
        seq_list: Vec<(String, Vec<u8>)>,
        source: Option<&str>,
        w: u32,
        k: u32,
        r: u32,
        min_span: u32,
    ) -> Result<(), std::io::Error> {
        let spec = ShmmrSpec {
            w,
            k,
            r,
            min_span,
            sketch: false,
        };
        self.backend = Backend::MEMORY;
        let source = Some(source.unwrap().to_string());
        let mut sdb = seq_db::CompactSeqDB::new(spec.clone());
        let seq_vec = seq_list
            .into_iter()
            .enumerate()
            .map(|(sid, v)| (sid as u32, source.clone(), v.0, v.1))
            .collect::<Vec<(u32, Option<String>, String, Vec<u8>)>>();
        sdb.load_seqs_from_seq_vec(&seq_vec);

        self.shmmr_spec = Some(spec);
        let mut seq_index = FxHashMap::<(String, Option<String>), (u32, u32)>::default();
        let mut seq_info = FxHashMap::<u32, (String, Option<String>, u32)>::default();
        sdb.seqs.iter().for_each(|v| {
            seq_index.insert((v.name.clone(), v.source.clone()), (v.id, v.len as u32));
            seq_info.insert(v.id, (v.name.clone(), v.source.clone(), v.len as u32));
        });
        self.seq_index = Some(seq_index);
        self.seq_info = Some(seq_info);
        self.seq_db = Some(sdb);
        Ok(())
    }

    #[allow(clippy::type_complexity)]
    pub fn query_fragment_to_hps(
        &self,
        seq: Vec<u8>,
        penalty: f32,
        max_count: Option<u32>,
        max_count_query: Option<u32>,
        max_count_target: Option<u32>,
        max_aln_span: Option<u32>,
    ) -> Option<Vec<(u32, Vec<(f32, Vec<aln::HitPair>)>)>> {
        let shmmr_spec = &self.shmmr_spec.as_ref().unwrap();
        if let Some(frag_map) = self.get_shmmr_map_internal() {
            let raw_query_hits = raw_query_fragment(frag_map, &seq, shmmr_spec);
            let res = aln::query_fragment_to_hps(
                raw_query_hits,
                &seq,
                shmmr_spec,
                penalty,
                max_count,
                max_count_query,
                max_count_target,
                max_aln_span,
            );
            Some(res)
        } else {
            None
        }
    }

    #[allow(clippy::type_complexity)]
    pub fn query_fragment_to_hps_from_mmap_file(
        &self,
        seq: Vec<u8>,
        penalty: f32,
        max_count: Option<u32>,
        max_count_query: Option<u32>,
        max_count_target: Option<u32>,
        max_aln_span: Option<u32>,
    ) -> Option<Vec<(u32, Vec<(f32, Vec<aln::HitPair>)>)>> {
        let shmmr_spec = self.shmmr_spec.as_ref().unwrap();
            let (frag_location_map, frag_map_file) = if self.backend == Backend::AGC {
                (&self.agc_db.as_ref().unwrap().frag_location_map, &self.agc_db.as_ref().unwrap().frag_map_file) 
            } else if self.backend == Backend::FRG {
                (&self.frg_db.as_ref().unwrap().frag_location_map, &self.frg_db.as_ref().unwrap().frag_map_file)  
            } else {
                panic!("the call query_fragment_to_hps_from_mmap_file() needs AGC or FRAG backend file");
            };
          
            let raw_query_hits = raw_query_fragment_from_mmap_midx(frag_location_map,
                frag_map_file, &seq, shmmr_spec);
            let res = aln::query_fragment_to_hps(
                raw_query_hits,
                &seq,
                shmmr_spec,
                penalty,
                max_count,
                max_count_query,
                max_count_target,
                max_aln_span,
            );
            Some(res)
    }

    pub fn get_sub_seq(
        &self,
        sample_name: String,
        ctg_name: String,
        bgn: usize,
        end: usize,
    ) -> Result<Vec<u8>, std::io::Error> {
        match self.backend {
            Backend::AGC => {
                Ok(self
                    .agc_db
                    .as_ref()
                    .unwrap()
                    .agc_file
                    .get_sub_seq(sample_name, ctg_name, bgn, end))
            }
            Backend::MEMORY | Backend::FASTX => {
                let &(sid, _) = self
                    .seq_index
                    .as_ref()
                    .unwrap()
                    .get(&(ctg_name, Some(sample_name)))
                    .unwrap();
                Ok(self
                    .seq_db
                    .as_ref()
                    .unwrap()
                    .get_sub_seq_by_id(sid, bgn as u32, end as u32))
            }
            Backend::FRG => {
                let &(sid, _) = self
                    .seq_index
                    .as_ref()
                    .unwrap()
                    .get(&(ctg_name, Some(sample_name)))
                    .unwrap();
                Ok(self
                    .frg_db
                    .as_ref()
                    .unwrap()
                    .get_sub_seq_by_id(sid, bgn as u32, end as u32))
            }
            Backend::UNKNOWN => Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "fetching sequence fail, database type in not determined",
            )),
        }
    }

    pub fn get_seq(
        &self,
        sample_name: String,
        ctg_name: String,
    ) -> Result<Vec<u8>, std::io::Error> {
        match self.backend {
            Backend::AGC => Ok(self
                .agc_db
                .as_ref()
                .unwrap()
                .agc_file
                .get_seq(sample_name, ctg_name)),
            Backend::MEMORY | Backend::FASTX => {
                let &(sid, _) = self
                    .seq_index
                    .as_ref()
                    .unwrap()
                    .get(&(ctg_name, Some(sample_name)))
                    .unwrap();
                Ok(self.seq_db.as_ref().unwrap().get_seq_by_id(sid))
            }
            Backend::FRG => {
                let &(sid, _) = self
                    .seq_index
                    .as_ref()
                    .unwrap()
                    .get(&(ctg_name, Some(sample_name)))
                    .unwrap();
                Ok(self.frg_db.as_ref().unwrap().get_seq_by_id(sid))
            }
            Backend::UNKNOWN => Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "fetching sequence fail, database type in not determined",
            )),
        }
    }

    pub fn get_seq_by_id(&self, sid: u32) -> Result<Vec<u8>, std::io::Error> {
        match self.backend {
            Backend::AGC => {
                let (ctg_name, sample_name, _) = self.seq_info.as_ref().unwrap().get(&sid).unwrap(); //TODO: handle Option unwrap properly
                let ctg_name = ctg_name.clone();
                let sample_name = sample_name.as_ref().unwrap().clone();
                Ok(self
                    .agc_db
                    .as_ref()
                    .unwrap()
                    .agc_file
                    .get_seq(sample_name, ctg_name))
            }
            Backend::MEMORY | Backend::FASTX => {
                Ok(self.seq_db.as_ref().unwrap().get_seq_by_id(sid))
            }
            Backend::FRG => Ok(self.frg_db.as_ref().unwrap().get_seq_by_id(sid)),
            Backend::UNKNOWN => Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "fetching sequence fail, database type in not determined",
            )),
        }
    }

    pub fn get_sub_seq_by_id(
        &self,
        sid: u32,
        bgn: usize,
        end: usize,
    ) -> Result<Vec<u8>, std::io::Error> {
        match self.backend {
            Backend::AGC => {
                let (ctg_name, sample_name, _) = self.seq_info.as_ref().unwrap().get(&sid).unwrap(); //TODO: handle Option unwrap properly
                let ctg_name = ctg_name.clone();
                let sample_name = sample_name.as_ref().unwrap().clone();
                Ok(self
                    .agc_db
                    .as_ref()
                    .unwrap()
                    .agc_file
                    .get_sub_seq(sample_name, ctg_name, bgn, end))
            }
            Backend::MEMORY | Backend::FASTX => Ok(self
                .seq_db
                .as_ref()
                .unwrap()
                .get_sub_seq_by_id(sid, bgn as u32, end as u32)),
            Backend::FRG => Ok(self
                .frg_db
                .as_ref()
                .unwrap()
                .get_sub_seq_by_id(sid, bgn as u32, end as u32)),
            Backend::UNKNOWN => Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "fetching sequence fail, database type in not determined",
            )),
        }
    }

    pub fn get_principal_bundles(
        &self,
        min_count: usize,
        path_len_cutoff: usize,
        keeps: Option<Vec<u32>>,
    ) -> Vec<Vec<(u64, u64, u8)>> {
        if let Some(frag_map) = self.get_shmmr_map_internal() {
            let adj_list = seq_db::frag_map_to_adj_list(frag_map, min_count as usize, keeps);

            seq_db::get_principal_bundles_from_adj_list(frag_map, &adj_list, path_len_cutoff)
                .0
                .into_iter()
                .map(|p| p.into_iter().map(|v| (v.0, v.1, v.2)).collect())
                .collect::<Vec<Vec<(u64, u64, u8)>>>()
        } else {
            vec![]
        }
    }

    fn get_vertex_map_from_principal_bundles(
        &self,
        pb: Vec<Vec<(u64, u64, u8)>>,
    ) -> FxHashMap<(u64, u64), (usize, u8, usize)> {
        // count segment for filtering, some unidirectional seg may have both forward and reverse in the principle bundles
        // let mut seg_count = FxHashMap::<(u64, u64), usize>::default();
        // pb.iter().for_each(|bundle| {
        //    bundle.iter().for_each(|v| {
        //        *seg_count.entry((v.0, v.1)).or_insert(0) += 1;
        //    })
        // });

        pb.iter()
            .enumerate()
            .flat_map(|(bundle_id, path)| {
                path.iter()
                    .enumerate()
                    //.filter(|(_, &v)| *seg_count.get(&(v.0, v.1)).unwrap_or(&0) == 1)
                    .map(|(p, v)| ((v.0, v.1), (bundle_id, v.2, p)))
                    .collect::<Vec<((u64, u64), (usize, u8, usize))>>()
            })
            .collect()
    }

    #[allow(clippy::type_complexity)] // TODO: Define the type for readability
    pub fn get_principal_bundle_decomposition(
        &self,
        min_count: usize,
        path_len_cutoff: usize,
        decomp_fasta_db: Option<&SeqIndexDB>,
        keeps: Option<Vec<u32>>,
    ) -> (
        Vec<(usize, usize, Vec<(u64, u64, u8)>)>,
        Vec<(
            u32,
            Vec<((u64, u64, u32, u32, u8), Option<(usize, u8, usize)>)>,
        )>,
    ) {
        fn get_smps(seq: Vec<u8>, shmmr_spec: &ShmmrSpec) -> Vec<(u64, u64, u32, u32, u8)> {
            let shmmrs = sequence_to_shmmrs(0, &seq, shmmr_spec, false);
            seq_db::pair_shmmrs(&shmmrs)
                .par_iter()
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
                .collect::<Vec<(u64, u64, u32, u32, u8)>>()
        }

        let pb = self.get_principal_bundles(min_count, path_len_cutoff, keeps);
        //println!("DBG: # bundles {}", pb.len());

        let mut vertex_to_bundle_id_direction_pos =
            self.get_vertex_map_from_principal_bundles(pb.clone()); //not efficient but it is PyO3 limit now

        let mut seqid_smps: Vec<(u32, Vec<(u64, u64, u32, u32, u8)>)> = self
            .seq_info
            .clone()
            .unwrap_or_default()
            .iter()
            .map(|(sid, data)| {
                let (ctg_name, source, _) = data;
                let source = source.clone().unwrap();
                let seq = self.get_seq(source, ctg_name.clone()).unwrap();
                (*sid, get_smps(seq, &self.shmmr_spec.clone().unwrap()))
            })
            .collect();

        // data for reordering the bundles and for re-ordering them along the sequences
        let mut bundle_id_to_directions = FxHashMap::<usize, Vec<u32>>::default();
        let mut bundle_id_to_orders = FxHashMap::<usize, Vec<f32>>::default();
        seqid_smps.iter().for_each(|(_sid, smps)| {
            let mut bundle_visited = FxHashSet::<usize>::default();
            smps.iter().enumerate().for_each(|(order, v)| {
                if let Some(bid) = vertex_to_bundle_id_direction_pos.get(&(v.0, v.1)) {
                    if !bundle_visited.contains(&bid.0) {
                        bundle_id_to_orders
                            .entry(bid.0)
                            .or_default()
                            .push(order as f32);
                        bundle_visited.insert(bid.0);
                    }
                    let direction = match bid.1 == v.4 {
                        true => 0,
                        false => 1,
                    };
                    bundle_id_to_directions
                        .entry(bid.0)
                        .or_default()
                        .push(direction);
                }
            })
        });

        // determine the bundles' overall orders and directions by consensus voting
        let mut bundle_mean_order_direction = (0..pb.len())
            .into_iter()
            .map(|bid| {
                if let Some(orders) = bundle_id_to_orders.get(&bid) {
                    let sum: f32 = orders.iter().sum();
                    let mean_ord = sum / (orders.len() as f32);
                    let mean_ord = mean_ord as usize;
                    let directions = bundle_id_to_directions.get(&bid).unwrap();
                    let dir_sum = directions.iter().sum::<u32>() as usize;
                    let direction = if dir_sum < (directions.len() >> 1) {
                        0_u8
                    } else {
                        1_u8
                    };
                    (mean_ord, bid, direction)
                } else {
                    let mean_ord = usize::MAX;
                    (mean_ord, bid, 0)
                }
            })
            .collect::<Vec<(usize, usize, u8)>>();

        //println!("DBG: length of bundle_mean_order_direction: {}", bundle_mean_order_direction.len());

        bundle_mean_order_direction.sort();
        // re-order the principal bundles
        let principal_bundles = bundle_mean_order_direction
            .iter()
            .map(|(ord, bid, direction)| {
                let bundle = if *direction == 1 {
                    let rpb = pb[*bid]
                        .iter()
                        .rev()
                        .map(|v| (v.0, v.1, 1 - v.2))
                        .collect::<Vec<(u64, u64, u8)>>();
                    rpb.iter().enumerate().for_each(|(p, v)| {
                        vertex_to_bundle_id_direction_pos.insert((v.0, v.1), (*bid, v.2, p));
                        // override what in the hashmap
                    });
                    rpb
                } else {
                    pb[*bid].clone()
                };

                (*bid, *ord, bundle)
            })
            .collect::<Vec<(usize, usize, Vec<(u64, u64, u8)>)>>();

        if let Some(seq_db) = decomp_fasta_db {
            seqid_smps = seq_db
                .seq_info
                .clone()
                .unwrap_or_default()
                .iter()
                .map(|(sid, data)| {
                    let (ctg_name, source, _) = data;
                    let source = source.clone().unwrap();
                    let seq = seq_db.get_seq(source, ctg_name.clone()).unwrap();
                    (*sid, get_smps(seq, &self.shmmr_spec.clone().unwrap()))
                })
                .collect();
        }

        // loop through each sequence and generate the decomposition for the sequence
        let seqid_smps_with_bundle_id_seg_direction = seqid_smps
            .iter()
            .map(|(sid, smps)| {
                let smps = smps
                    .iter()
                    .map(|v| {
                        let seg_match = vertex_to_bundle_id_direction_pos.get(&(v.0, v.1)).copied();
                        (*v, seg_match)
                    })
                    .collect::<Vec<((u64, u64, u32, u32, u8), Option<(usize, u8, usize)>)>>();
                (*sid, smps)
            })
            .collect::<Vec<(
                u32,
                Vec<((u64, u64, u32, u32, u8), Option<(usize, u8, usize)>)>,
            )>>();

        (principal_bundles, seqid_smps_with_bundle_id_seg_direction)
    }

    pub fn generate_mapg_gfa(
        &self,
        min_count: usize,
        filepath: &str,
        method: &str,
        keeps: Option<Vec<u32>>,
    ) -> Result<(), std::io::Error> {
        let get_seq_by_id = |sid| -> Vec<u8> {
            match self.backend {
                Backend::AGC => {
                    let (ctg_name, sample_name, _) =
                        self.seq_info.as_ref().unwrap().get(&sid).unwrap(); //TODO: handle Option unwrap properly
                    let ctg_name = ctg_name.clone();
                    let sample_name = sample_name.as_ref().unwrap().clone();
                    self.agc_db
                        .as_ref()
                        .unwrap()
                        .agc_file
                        .get_seq(sample_name, ctg_name)
                }
                Backend::MEMORY => self.seq_db.as_ref().unwrap().get_seq_by_id(sid),
                Backend::FASTX => self.seq_db.as_ref().unwrap().get_seq_by_id(sid),
                Backend::FRG => self.frg_db.as_ref().unwrap().get_seq_by_id(sid),
                Backend::UNKNOWN => vec![],
            }
        };

        let frag_map = self.get_shmmr_map_internal();
        if frag_map.is_none() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "can get frag_map",
            ));
        }
        let mut overlaps =
            FxHashMap::<(ShmmrGraphNode, ShmmrGraphNode), Vec<(u32, u8, u8)>>::default();
        let mut frag_id = FxHashMap::<(u64, u64), usize>::default();
        let mut id = 0_usize;

        let frag_map = frag_map.unwrap();

        let adj_list = if method == "from_fragmap" {
            seq_db::frag_map_to_adj_list(frag_map, min_count, keeps)
        } else {
            let keeps = keeps.map(FxHashSet::<u32>::from_iter);

            self.seq_info
                .as_ref()
                .unwrap()
                .keys()
                .into_iter()
                .copied()
                .collect::<Vec<u32>>()
                .into_par_iter()
                .flat_map(|sid| {
                    let seq = get_seq_by_id(sid);
                    let mc = if let Some(keeps) = &keeps {
                        if keeps.contains(&sid) {
                            0
                        } else {
                            min_count
                        }
                    } else {
                        min_count
                    };
                    seq_db::generate_smp_adj_list_for_seq(
                        &seq,
                        sid,
                        frag_map,
                        self.shmmr_spec.as_ref().unwrap(),
                        mc,
                    )
                })
                .collect::<AdjList>()
        };

        adj_list.iter().for_each(|(k, v, w)| {
            if v.0 <= w.0 {
                let key = (*v, *w);
                let val = (*k, v.2, w.2);
                overlaps.entry(key).or_insert_with(Vec::new).push(val);
                frag_id.entry((v.0, v.1)).or_insert_with(|| {
                    let c_id = id;
                    id += 1;
                    c_id
                });
                frag_id.entry((w.0, w.1)).or_insert_with(|| {
                    let c_id = id;
                    id += 1;
                    c_id
                });
            }
        });

        let mut out_file = BufWriter::new(File::create(filepath).unwrap());

        let kmer_size = self.shmmr_spec.as_ref().unwrap().k;
        out_file
            .write_all("H\tVN:Z:1.0\tCM:Z:Sparse Genome Graph Generated By pgr-tk\n".as_bytes())?;
        frag_id
            .iter()
            .try_for_each(|(smp, id)| -> Result<(), std::io::Error> {
                let hits = frag_map.get(smp).unwrap();
                let ave_len =
                    hits.iter().fold(0_u32, |len_sum, &s| len_sum + s.3 - s.2) / hits.len() as u32;
                let seg_line = format!(
                    "S\t{}\t*\tLN:i:{}\tSN:Z:{:016x}_{:016x}\n",
                    id,
                    ave_len + kmer_size,
                    smp.0,
                    smp.1
                );
                out_file.write_all(seg_line.as_bytes())?;
                Ok(())
            })?;

        overlaps
            .into_iter()
            .try_for_each(|(op, vs)| -> Result<(), std::io::Error> {
                let o1 = if op.0 .2 == 0 { "+" } else { "-" };
                let o2 = if op.1 .2 == 0 { "+" } else { "-" };
                let id0 = frag_id.get(&(op.0 .0, op.0 .1)).unwrap();
                let id1 = frag_id.get(&(op.1 .0, op.1 .1)).unwrap();
                let overlap_line = format!(
                    "L\t{}\t{}\t{}\t{}\t{}M\tSC:i:{}\n",
                    id0,
                    o1,
                    id1,
                    o2,
                    kmer_size,
                    vs.len()
                );
                out_file.write_all(overlap_line.as_bytes())?;
                Ok(())
            })?;

        Ok(())
    }

    pub fn write_mapg_idx(&self, filepath: &str) -> Result<(), std::io::Error> {
        let mut writer = BufWriter::new(File::create(filepath)?);

        if let Some(shmmr_spec) = self.shmmr_spec.clone() {
            writer.write_all(
                format!(
                    "K\t{}\t{}\t{}\t{}\t{}\n",
                    shmmr_spec.w,
                    shmmr_spec.k,
                    shmmr_spec.r,
                    shmmr_spec.min_span,
                    shmmr_spec.sketch
                )
                .as_bytes(),
            )?;
        }

        self.seq_info.as_ref().unwrap().iter().try_for_each(
            |(k, v)| -> Result<(), std::io::Error> {
                let line = format!(
                    "C\t{}\t{}\t{}\t{}\n",
                    k,
                    v.0,
                    v.1.clone().unwrap_or_else(|| "NA".to_string()),
                    v.2
                );
                writer.write_all(line.as_bytes())?;
                Ok(())
            },
        )?;

        let frag_map = self.get_shmmr_map_internal();
        if frag_map.is_none() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "fail to load index",
            ));
        };
        let frag_map = frag_map.unwrap();
        frag_map
            .iter()
            .try_for_each(|v| -> Result<(), std::io::Error> {
                v.1.iter()
                    .try_for_each(|vv| -> Result<(), std::io::Error> {
                        writer.write_all(
                            format!(
                                "F\t{:016x}_{:016x}\t{}\t{}\t{}\t{}\t{}\n",
                                v.0 .0, v.0 .1, vv.0, vv.1, vv.2, vv.3, vv.4
                            )
                            .as_bytes(),
                        )?;
                        Ok(())
                    })?;
                Ok(())
            })?;
        Ok(())
    }

    pub fn generate_principal_mapg_gfa(
        &self,
        min_count: usize,
        path_len_cutoff: usize,
        filepath: &str,
        keeps: Option<Vec<u32>>,
    ) -> Result<(), std::io::Error> {
        let frag_map = self.get_shmmr_map_internal();
        if frag_map.is_none() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "can't load index",
            ));
        };
        let frag_map = frag_map.unwrap();
        let adj_list = seq_db::frag_map_to_adj_list(frag_map, min_count, keeps);

        // println!("DBG: adj_list len {:?}", adj_list.len());

        let mut overlaps =
            FxHashMap::<(ShmmrGraphNode, ShmmrGraphNode), Vec<(u32, u8, u8)>>::default();
        let mut frag_id = FxHashMap::<(u64, u64), usize>::default();
        let mut id = 0_usize;
        let (pb, filtered_adj_list) =
            seq_db::get_principal_bundles_from_adj_list(frag_map, &adj_list, path_len_cutoff);

        // println!("DBG: pb len {:?}, filtered_adj_list len: {:?} ", pb.len(), filtered_adj_list.len());

        // TODO: we will remove this redundant conversion in the future
        let pb = pb
            .into_iter()
            .map(|p| p.into_iter().map(|v| (v.0, v.1, v.2)).collect())
            .collect::<Vec<Vec<(u64, u64, u8)>>>();

        let vertex_to_bundle_id_direction_pos = self.get_vertex_map_from_principal_bundles(pb);

        filtered_adj_list.iter().for_each(|(k, v, w)| {
            if v.0 <= w.0 {
                let key = (*v, *w);
                let val = (*k, v.2, w.2);
                overlaps.entry(key).or_insert_with(Vec::new).push(val);
                frag_id.entry((v.0, v.1)).or_insert_with(|| {
                    let c_id = id;
                    id += 1;
                    c_id
                });
                frag_id.entry((w.0, w.1)).or_insert_with(|| {
                    let c_id = id;
                    id += 1;
                    c_id
                });
            }
        });

        let mut out_file = BufWriter::new(File::create(filepath).unwrap());

        let kmer_size = self.shmmr_spec.as_ref().unwrap().k;
        out_file
            .write_all("H\tVN:Z:1.0\tCM:Z:Sparse Genome Graph Generated By pgr-tk\n".as_bytes())?;
        frag_id
            .iter()
            .try_for_each(|(smp, id)| -> Result<(), std::io::Error> {
                let hits = frag_map.get(smp).unwrap();
                let ave_len =
                    hits.iter().fold(0_u32, |len_sum, &s| len_sum + s.3 - s.2) / hits.len() as u32;
                let seg_line;
                if let Some(bundle_id) = vertex_to_bundle_id_direction_pos.get(smp) {
                    seg_line = format!(
                        "S\t{}\t*\tLN:i:{}\tSN:Z:{:016x}_{:016x}\tBN:i:{}\tBP:i:{}\n",
                        id,
                        ave_len + kmer_size,
                        smp.0,
                        smp.1,
                        bundle_id.0,
                        bundle_id.2
                    );
                } else {
                    seg_line = format!(
                        "S\t{}\t*\tLN:i:{}\tSN:Z:{:016x}_{:016x}\n",
                        id,
                        ave_len + kmer_size,
                        smp.0,
                        smp.1
                    );
                }
                out_file.write_all(seg_line.as_bytes())?;
                Ok(())
            })?;

        overlaps
            .into_iter()
            .try_for_each(|(op, vs)| -> Result<(), std::io::Error> {
                let o1 = if op.0 .2 == 0 { "+" } else { "-" };
                let o2 = if op.1 .2 == 0 { "+" } else { "-" };
                let id0 = frag_id.get(&(op.0 .0, op.0 .1)).unwrap();
                let id1 = frag_id.get(&(op.1 .0, op.1 .1)).unwrap();
                let overlap_line = format!(
                    "L\t{}\t{}\t{}\t{}\t{}M\tSC:i:{}\n",
                    id0,
                    o1,
                    id1,
                    o2,
                    kmer_size,
                    vs.len()
                );
                out_file.write_all(overlap_line.as_bytes())?;
                Ok(())
            })?;

        Ok(())
    }
}

impl SeqIndexDB {
    // depending on the storage type, return the corresponded index
    pub fn get_shmmr_map_internal(&self) -> Option<&seq_db::ShmmrToFrags> {
        match self.backend {
            Backend::AGC => None,
            Backend::FASTX => Some(&self.seq_db.as_ref().unwrap().frag_map),
            Backend::MEMORY => Some(&self.seq_db.as_ref().unwrap().frag_map),
            Backend::FRG => None,
            Backend::UNKNOWN => None,
        }
    }
}

pub fn get_fastx_reader(filepath: String) -> Result<GZFastaReader, std::io::Error> {
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
            FastaReader::new(gz_buf, &filepath, 256, false).unwrap(),
        ))
    } else {
        drop(gz_buf);
        Ok(GZFastaReader::RegularFile(
            FastaReader::new(std_buf, &filepath, 256, false).unwrap(),
        ))
    }
}
