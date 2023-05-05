// src/lib.rs
pub const VERSION_STRING: &'static str = env!("VERSION_STRING");
use pgr_db::aln::{self, HitPair};
use pgr_db::graph_utils::{AdjList, ShmmrGraphNode};
use pgr_db::seq_db;
//use pgr_db::seqs2variants;
use pgr_db::shmmrutils::{sequence_to_shmmrs, DeltaPoint, ShmmrSpec};

#[cfg(feature = "with_agc")]
use pgr_db::agc_io;

use pgr_db::fasta_io;
use pyo3::exceptions;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use pyo3::Python;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use pgr_db::ext::Backend;

/// Get the revision (git-hashtag) of the build
#[pyfunction]
pub fn pgr_lib_version() -> PyResult<String> {
    Ok(VERSION_STRING.to_string())
}

/// A class that stores pangenome indices and sequences with multiple backend storage options (AGC, fasta file, memory)
/// Large set of genomic sequences, a user should use AGC backend. A binary file provides the command ``pgr-mdb``
/// which can read an AGC to create the index file. For example, we can create the index files from an AGC file::
///
///     # create a file that contains a list of file that contains a set of files from which we want to build the indices
///  
///     $ echo HPRC-y1-rebuild-04252022.agc > filelist
///  
///     # using pgr-mdb to create the index files, for 97 haplotyped genome assembly from HPRC year one release,
///     # it takes about 30 to 40 min to create the index files
///
///     $ pgr-mdb filelist HPRC-y1-rebuild-04252022
///
///     # two index files will be created by the pgr-mdb command
///     # one with a suffix .mdb and another one with a suffix .midx
///     # when we use the load_from_agc_index() method, all three files, e.g., genomes.agc, genomes.mdb and
///     # genomes.midx should have the same prefix as the parameter used to call  load_from_agc_index() method
///
/// One can also create index and load the sequences from a fasta file using ```load_from_fastx()``` methods.
/// Currently, this might be a good option for mid-size dataset (up to a couple of hundred megabases).
///
/// Or, a user can load the sequence from memory using a Python list. This is convenient when one needs to
/// rebuild the SHIMMER index with different parameters for a different resolution.
///  
/// Once the index is built, the database can be queried quickly by using the ``query_fragment()`` or
/// the ``query_fragment_to_hps()`` method.
///  
#[pyclass]
struct SeqIndexDB {
    /// Rust internal:
    db_internal: pgr_db::ext::SeqIndexDB,
    pub principal_bundles: Option<(usize, usize, Vec<Vec<(u64, u64, u8)>>)>,
}

#[pymethods]
impl SeqIndexDB {
    /// constructor, take no argument
    #[new]
    pub fn new() -> Self {
        SeqIndexDB {
            db_internal: pgr_db::ext::SeqIndexDB {
                seq_db: None,
                frg_db: None,
                #[cfg(feature = "with_agc")]
                agc_db: None,
                shmmr_spec: None,
                seq_index: None,
                seq_info: None,
                backend: Backend::UNKNOWN,
            },
            principal_bundles: None,
        }
    }

    /// use AGC file for sequences and the index created from an AGC file
    ///
    /// Parameters
    /// ----------
    ///
    /// prefix: string
    ///     the prefix to the `.agc`, `.mdb` and `.midx` files
    ///
    /// Returns
    /// -------
    ///
    /// None or I/O Error
    ///
    #[cfg(feature = "with_agc")]
    #[pyo3(text_signature = "($self, prefix)")]
    pub fn load_from_agc_index(&mut self, prefix: String) -> PyResult<()> {
        self.db_internal.load_from_agc_index(prefix)?;
        Ok(())
    }

    #[pyo3(text_signature = "($self, prefix)")]
    pub fn load_from_frg_index(&mut self, prefix: String) -> PyResult<()> {
        self.db_internal.load_from_frg_index(prefix)?;
        Ok(())
    }

    /// load and create the index created from a fasta / fastq file
    ///
    /// Parameters
    /// ----------
    ///
    ///filepath : string
    ///     the path the fasta or fastq file
    ///
    /// w : int
    ///     the window size of the shimmer index, default to 80
    ///
    /// k : int
    ///     the k-mer size of the shimmer index, default to 56
    ///
    /// r : int
    ///     the reduction factor of the shimmer index, default to 4
    ///
    /// min_span : int
    ///     the min_span ofr the shimmer index, default to 8
    ///
    /// Returns
    /// -------
    ///
    /// None or I/O Error
    ///     None
    ///
    #[pyo3(signature = (filepath, w=80, k=56, r=4, min_span=64))]
    pub fn load_from_fastx(
        &mut self,
        filepath: String,
        w: u32,
        k: u32,
        r: u32,
        min_span: u32,
    ) -> PyResult<()> {
        self.db_internal
            .load_from_fastx(filepath, w, k, r, min_span)?;
        Ok(())
    }

    #[pyo3(text_signature = "($self)")]
    pub fn append_from_fastx(&mut self, filepath: String) -> PyResult<()> {
        assert!(
            self.db_internal.backend == Backend::FASTX,
            "Only DB created with load_from_fastx() can add data from another fastx file"
        );
        let sdb = self.db_internal.seq_db.as_mut().unwrap();
        sdb.load_seqs_from_fastx(filepath)?;
        Ok(())
    }

    /// load and create the index created from a python list
    ///
    /// Parameters
    /// ----------
    ///
    /// seq_list : list
    ///     a list of tuple of the form (sequence_id : int, sequence_name : string, sequence: list of bytes)
    ///
    /// source : string
    ///     a string indicating the source of the sequence, default to "Memory"
    ///
    /// w : int
    ///     the window size of the shimmer index, default to 80
    ///
    /// k : int
    ///     the k-mer size of the shimmer index, default to 56
    ///
    /// r : int
    ///     the reduction factor of the shimmer index, default to 4
    ///
    /// min_span : int
    ///     the min_span ofr the shimmer index, default to 8
    ///
    /// Returns
    /// -------
    ///
    /// None or I/O Error
    ///     None
    ///
    #[pyo3(signature = (seq_list, source="Memory", w=80, k=56, r=4, min_span=8))]
    pub fn load_from_seq_list(
        &mut self,
        seq_list: Vec<(String, Vec<u8>)>,
        source: Option<&str>,
        w: u32,
        k: u32,
        r: u32,
        min_span: u32,
    ) -> PyResult<()> {
        self.db_internal
            .load_from_seq_list(seq_list, source, w, k, r, min_span)?;

        Ok(())
    }

    /// get a dictionary that maps (ctg_name, source) -> (id, len)
    #[getter]
    pub fn get_seq_index(
        &self,
    ) -> PyResult<Option<FxHashMap<(String, Option<String>), (u32, u32)>>> {
        Ok(self.db_internal.seq_index.clone())
    }

    /// a dictionary that maps id -> (ctg_name, source, len)
    #[getter]
    pub fn get_seq_info(&self) -> PyResult<Option<FxHashMap<u32, (String, Option<String>, u32)>>> {
        Ok(self.db_internal.seq_info.clone())
    }

    /// use a fragment of sequence to query the database to get all hits
    ///
    /// Parameters
    /// ----------
    ///
    /// seq : list
    ///     the sequence in bytes used for query
    ///
    /// Returns
    /// -------
    ///
    /// list
    ///   a list of hits in the format (shimmer_pair, query_fragment, target_fragments), where
    ///     - shimmer_pair: (int, int), tuple of the shimmer_pair
    ///     - query_fragment: (int, int, int) = (start_coordinate, end_coordinate, orientation)
    ///     - target_fragments: a list of ``FragmentSignature``: (frg_id, seq_id, bgn, end,
    ///       orientation(to the shimmer pair)) defined as::
    ///
    ///           pub type FragmentSignature = (u32, u32, u32, u32, u8);
    ///       
    ///
    #[pyo3(text_signature = "($self, seq)")]
    pub fn query_fragment(
        &self,
        seq: Vec<u8>,
    ) -> PyResult<Vec<((u64, u64), (u32, u32, u8), Vec<seq_db::FragmentSignature>)>> {
        match self.db_internal.backend {
            #[cfg(feature = "with_agc")]
            Backend::AGC => {
                let (frag_location_map, frag_map_file) = (
                    &self.db_internal.agc_db.as_ref().unwrap().frag_location_map,
                    &self.db_internal.agc_db.as_ref().unwrap().frag_map_file,
                );
                let shmmr_spec = self.db_internal.shmmr_spec.as_ref().unwrap().clone();
                Ok(pgr_db::seq_db::raw_query_fragment_from_mmap_midx(
                    frag_location_map,
                    frag_map_file,
                    &seq,
                    &shmmr_spec,
                ))
            }
            Backend::FRG => {
                let (frag_location_map, frag_map_file) = (
                    &self.db_internal.frg_db.as_ref().unwrap().frag_location_map,
                    &self.db_internal.frg_db.as_ref().unwrap().frag_map_file,
                );
                let shmmr_spec = self.db_internal.shmmr_spec.as_ref().unwrap().clone();
                Ok(pgr_db::seq_db::raw_query_fragment_from_mmap_midx(
                    frag_location_map,
                    frag_map_file,
                    &seq,
                    &shmmr_spec,
                ))
            }
            Backend::MEMORY | Backend::FASTX => {
                let shmmr_spec = &self.db_internal.shmmr_spec.as_ref().unwrap();
                let shmmr_to_frags = self.get_shmmr_map_internal().unwrap();
                let res: Vec<((u64, u64), (u32, u32, u8), Vec<seq_db::FragmentSignature>)> =
                    seq_db::raw_query_fragment(shmmr_to_frags, &seq, shmmr_spec);
                Ok(res)
            }
            Backend::UNKNOWN => Ok(vec![]),
        }
    }

    /// use a fragment of sequence to query the database to get all hits and sort it by the data base sequence id
    ///
    /// Parameters
    /// ----------
    ///
    /// seq : list
    ///     the sequence in bytes used for query
    ///
    /// Returns
    /// -------
    ///
    /// dict
    ///   a dictionary maps sequence id in the database to a list of tuple (position0, position1, direction)
    ///       
    ///
    #[pyo3(text_signature = "($self, seq)")]
    pub fn get_match_positions_with_fragment(
        &self,
        seq: Vec<u8>,
    ) -> PyResult<FxHashMap<u32, Vec<(u32, u32, u8)>>> {
        let shmmr_spec = &self.db_internal.shmmr_spec.as_ref().unwrap();
        match self.db_internal.backend {
            Backend::FASTX | Backend::MEMORY => {
                let shmmr_to_frags = self.get_shmmr_map_internal().unwrap();
                let res =
                    seq_db::get_match_positions_with_fragment(shmmr_to_frags, &seq, shmmr_spec);
                Ok(res)
            }
            _ => Err(exceptions::PyException::new_err(
                "This method only support FASTX or MEMORY backend.",
            )),
        }
    }

    /// use a fragment of sequence to query the database to get all hits
    ///
    /// sparse dynamic programming is performed to long chain of alignment
    ///  
    /// Parameters
    /// ----------
    /// seq : list of bytes
    ///    a list of bytes representing the DNA sequence
    ///
    /// penalty : float
    ///    the gap penalty factor used in sparse dynamic programming for finding the hits
    ///
    /// merge_range_tol : int
    ///    a parameter used to merge the alignment ranges
    ///
    /// max_count : int
    ///    only use the shimmer pairs that less than the ``max_count`` for sparse dynamic programming
    ///
    /// max_query_count : int
    ///    only use the shimmer pairs that less than the ``max_count`` in the query sequence for sparse dynamic programming
    ///
    /// max_query_count : int
    ///    only use the shimmer pairs that less than the ``max_count`` in the target sequence for sparse dynamic programming
    ///
    /// max_aln_span : int
    ///    the size of span used in the sparse dynamic alignment for finding the hits
    ///
    /// Returns
    /// -------
    ///
    /// list
    ///     a list of tuples of
    ///     (``target_sequence_id``, (``score``, ``list_of_the_hit_pairs``)), where
    ///     the ``list_of_the_hit_pairs`` is a list of tuples of
    ///     ((``query_start``, ``query_end``, ``query_orientation``),
    ///     (``target_start``, ``target_end``, ``target_orientation``))
    #[pyo3(
        text_signature = "($self, seq, penalty, max_count, max_query_count, max_target_count, max_aln_span)"
    )]
    pub fn query_fragment_to_hps(
        &self,
        seq: Vec<u8>,
        penalty: f32,
        max_count: Option<u32>,
        max_count_query: Option<u32>,
        max_count_target: Option<u32>,
        max_aln_span: Option<u32>,
    ) -> PyResult<Vec<(u32, Vec<(f32, Vec<aln::HitPair>)>)>> {
        match self.db_internal.backend {
            #[cfg(feature = "with_agc")]
            Backend::AGC => Ok(self
                .db_internal
                .query_fragment_to_hps_from_mmap_file(
                    seq,
                    penalty,
                    max_count,
                    max_count_query,
                    max_count_target,
                    max_aln_span,
                )
                .unwrap()),
            Backend::FRG => Ok(self
                .db_internal
                .query_fragment_to_hps_from_mmap_file(
                    seq,
                    penalty,
                    max_count,
                    max_count_query,
                    max_count_target,
                    max_aln_span,
                )
                .unwrap()),
            Backend::MEMORY | Backend::FASTX => Ok(self
                .db_internal
                .query_fragment_to_hps(
                    seq,
                    penalty,
                    max_count,
                    max_count_query,
                    max_count_target,
                    max_aln_span,
                )
                .unwrap()),
            Backend::UNKNOWN => Ok(vec![]),
        }
    }

    /// Given a sequence context, this function maps the specific positions in the context
    /// to the sequences in the database. The context sequence is aligned to the sequences
    /// in the database with sparse dynamic programming, then the regions include the
    /// positions of interest are identified. A wavefront alignment is performed to pin
    /// down the exact mapped positions in the sequences in the database.
    ///
    /// Parameters
    /// ----------
    /// positions : list of integer
    ///    a list of integers of the position to map
    ///  
    /// seq : list of bytes
    ///    a list of bytes representing the DNA sequence providing to context
    ///
    /// penalty : float
    ///    the gap penalty factor used in sparse dynamic programming for finding the hits
    ///
    /// merge_range_tol : int
    ///    a parameter used to merge the alignment ranges
    ///
    /// max_count : int
    ///    only use the shimmer pairs that less than the ``max_count`` for sparse dynamic programming
    ///
    /// max_query_count : int
    ///    only use the shimmer pairs that less than the ``max_count`` in the query sequence for sparse dynamic programming
    ///
    /// max_query_count : int
    ///    only use the shimmer pairs that less than the ``max_count`` in the target sequence for sparse dynamic programming
    ///
    /// max_aln_span : int
    ///    the size of span used in the sparse dynamic alignment for finding the hits
    ///
    /// Returns
    /// -------
    ///
    /// list
    ///     a list of tuples of
    ///     (``position_in_the_context``,
    ///     (``target_seq_id``, ``target_position``, ``orientation``),
    ///     (``context_end``, ``context_end``),
    ///     (``target_end``, ``target_end``))
    ///
    ///     the sequences from (``context_end``, ``context_end``) in the context sequence and
    ///     the sequences from (``target_end``, ``target_end``) in the target sequence are
    ///     used for the detailed alignment to pin down the exact mapped positions.
    ///
    #[pyo3(
        text_signature = "($self, positions, seq, penalty, max_count, max_query_count, max_target_count, max_aln_span)"
    )]
    pub fn map_positions_in_seq(
        &self,
        positions: Vec<u32>,
        seq: Vec<u8>,
        penalty: f32,
        max_count: Option<u32>,
        max_count_query: Option<u32>,
        max_count_target: Option<u32>,
        max_aln_span: Option<u32>,
    ) -> PyResult<Vec<(u32, (u32, u32, u8), (u32, u32), (u32, u32))>> {
        let shmmr_spec = self.db_internal.shmmr_spec.as_ref().unwrap();
        let mut all_alns = {
            let raw_query_hits = self.query_fragment(seq.clone()).unwrap();
            aln::query_fragment_to_hps(
                raw_query_hits,
                &seq,
                shmmr_spec,
                penalty,
                max_count,
                max_count_query,
                max_count_target,
                max_aln_span,
            )
        };

        // for reach position, we find the left_match and right_match shimmer pair that sandwiched the
        // positions. the algorithm is based on brute force search and should be optimized in the future
        let mut pos2hits = FxHashMap::<u32, Vec<(u32, f32, HitPair, HitPair)>>::default();
        all_alns.iter_mut().for_each(|(t_id, alns)| {
            alns.iter_mut().for_each(|(score, hits)| {
                hits.sort();
                positions.iter().for_each(|&pos| {
                    let mut out: Vec<(u32, f32, HitPair, HitPair)> = vec![];
                    let mut left_match = None;
                    let mut right_match = None;
                    hits.iter().for_each(|&(v, w)| {
                        // println!("{:?} {:?} {:?}", pos, v, w);
                        if v.0 < pos {
                            left_match = Some((v, w));
                            // println!("left set: {:?} {:?}", v, w);
                        }
                        if right_match.is_none() && pos < v.1 {
                            right_match = Some((v, w));
                            // println!("right set: {:?} {:?}", v, w);
                        }
                    });
                    if left_match.is_some() && right_match.is_some() {
                        out.push((*t_id, *score, left_match.unwrap(), right_match.unwrap()));
                    };
                    pos2hits.entry(pos).or_insert(vec![]).extend(out);
                });
            });
        });

        // fetch the sequence for each match if possible
        let mut out = vec![];
        if self.db_internal.seq_info.is_none() {
            return Ok(out);
        };
        pos2hits.iter().for_each(|(pos, hits)| {
            hits.iter()
                .for_each(|(seq_id, _score, left_match, right_match)| {
                    let (ctg, src, t_len) = self
                        .db_internal
                        .seq_info
                        .as_ref()
                        .unwrap()
                        .get(&seq_id)
                        .unwrap(); //TODO, check if seq_info is None
                    let same_orientation = if left_match.0 .2 == left_match.1 .2 {
                        true
                    } else {
                        false
                    };

                    let qb = left_match.0 .0;
                    let qe = right_match.0 .1;
                    let tb;
                    let te;

                    match same_orientation {
                        true => {
                            tb = left_match.1 .0;
                            te = right_match.1 .1;
                        }
                        false => {
                            tb = right_match.1 .0 - shmmr_spec.k;
                            te = left_match.1 .1 - shmmr_spec.k;
                        }
                    };
                    if tb >= te {
                        // println!("{:?} {:?} {} {} {} {}", left_match, right_match, qb, qe, tb, te);
                        // TBD: raise an warning? or error? The coordinates are not consistent wit the shimmer alignment orientation
                        return;
                    }
                    let mut t_seq = self
                        .get_sub_seq(
                            src.clone().unwrap().to_string(),
                            ctg.clone(),
                            tb as usize,
                            te as usize,
                        )
                        .unwrap();

                    if !same_orientation {
                        t_seq = fasta_io::reverse_complement(&t_seq);
                    }
                    let q_seq = seq[qb as usize..qe as usize].to_vec();
                    let ovlp =
                        pgr_db::shmmrutils::match_reads(&q_seq, &t_seq, true, 0.10, 1, 1, 1000);
                    // if ovlp.is_none() {
                    //    println!("aln fail for pos: {:?} {:?} {:?}", pos, left_match, right_match);
                    //    println!("qseq: {}", String::from_utf8_lossy(&q_seq[..]));
                    //    println!("tseq: {}", String::from_utf8_lossy(&t_seq[..]));
                    // }
                    if ovlp.is_some() {
                        let dpos = pos - qb;

                        let mut delta = ovlp.unwrap().deltas.unwrap();

                        delta.push(DeltaPoint { x: 0, y: 0, dk: 0 });

                        let mut dref = None;

                        for dp in delta.iter() {
                            if dp.x <= dpos {
                                dref = Some((dp.x, dp.y));
                                break;
                            };
                        }

                        let dref = dref.unwrap();

                        let orientation = if same_orientation { 0_u8 } else { 1_u8 };
                        let dpos = dpos + dref.1 - dref.0;
                        let (tb, te, tpos) = if same_orientation {
                            (tb, te, tb + dpos)
                        } else {
                            (*t_len - te, *t_len - tb, *t_len - (te - dpos))
                        };

                        out.push((*pos, (*seq_id, tpos, orientation), (qb, qe), (tb, te)));
                    }
                });
        });

        Ok(out)
    }

    /// count the number of shimmer hits in the database
    ///
    /// Parameters
    /// ----------
    ///
    /// shmmr_pair : tuple
    ///     a shimmer pair used for query
    ///
    /// Returns
    /// -------
    ///
    /// int
    ///     number of hits
    #[pyo3(text_signature = "($self, shmmr_pair)")]
    pub fn get_shmmr_pair_count(&self, shmmr_pair: (u64, u64)) -> PyResult<usize> {
        if let Some(shmmr_to_frags) = self.get_shmmr_map_internal() {
            if shmmr_to_frags.contains_key(&shmmr_pair) {
                Ok(shmmr_to_frags.get(&shmmr_pair).unwrap().len())
            } else {
                Ok(0)
            }
        } else {
            Err(exceptions::PyException::new_err(
                "This method only support FASTX or MEMORY backend.",
            ))
        }
    }

    /// count the number of shimmer hits partitioned by the source file in the database
    ///
    /// Parameters
    /// ----------
    ///
    /// shmmr_pair : tuple
    ///     a shimmer pair used for query
    ///
    /// max_unique_count : int
    ///     a integer to filter out shimmer pairs with count that are greater
    ///     than the `max_unique_count`  
    ///
    /// Returns
    /// -------
    ///
    /// list
    ///     a list of the tuple (source_name : string, count : int)
    ///
    #[pyo3(signature = (shmmr_pair, max_unique_count))]
    pub fn get_shmmr_pair_source_count(
        &self,
        shmmr_pair: (u64, u64),
        max_unique_count: Option<usize>,
    ) -> PyResult<Vec<(String, usize)>> {
        let mut count = FxHashMap::<String, usize>::default();
        let shmmr_to_frags = self.get_shmmr_map_internal();
        if shmmr_to_frags.is_none() {
            return Err(exceptions::PyException::new_err(
                "This method only support FASTX or MEMORY backend.",
            ));
        };
        let shmmr_to_frags = shmmr_to_frags.unwrap();

        if shmmr_to_frags.contains_key(&shmmr_pair) {
            shmmr_to_frags
                .get(&shmmr_pair)
                .unwrap()
                .iter()
                .for_each(|v| {
                    let sid = v.1;
                    let source = self
                        .db_internal
                        .seq_info
                        .as_ref()
                        .unwrap()
                        .get(&sid)
                        .unwrap()
                        .1
                        .as_ref()
                        .unwrap_or(&"".to_string())
                        .clone();
                    *count.entry(source).or_insert(0) += 1;
                });

            let out = count
                .into_par_iter()
                .filter(|(_k, v)| {
                    if let Some(muc) = max_unique_count {
                        if *v > muc {
                            false
                        } else {
                            true
                        }
                    } else {
                        true
                    }
                })
                .map(|(k, v)| (k, v))
                .collect::<Vec<(String, usize)>>();
            Ok(out)
        } else {
            Ok(vec![])
        }
    }

    /// Output the specific of the shimmer used to build the index
    ///
    /// Returns
    /// -------
    ///
    /// tuple
    ///     (window_size, k_mer_size, reduction_factor, min_space, use_sketch)
    ///
    pub fn get_shmmr_spec(&self) -> PyResult<Option<(u32, u32, u32, u32, bool)>> {
        if let Some(spec) = self.db_internal.shmmr_spec.as_ref() {
            Ok(Some((spec.w, spec.k, spec.r, spec.min_span, spec.sketch)))
        } else {
            Ok(None)
        }
    }

    /// get the ``shmmr_pair`` to ``fragment_id`` map in Python
    ///
    /// this can be very expensive to generate the Python objects of a large hashmap in Rust
    ///
    /// Returns
    /// -------
    ///
    /// dict
    ///     the ``shmmr_pair`` to ``fragments`` map
    ///
    ///     fragments: a list of ``FragmentSignature``: (frg_id, seq_id, bgn, end,
    ///     orientation(to the shimmer pair)) defined as::
    ///
    ///         pub type FragmentSignature = (u32, u32, u32, u32, u8);
    ///
    pub fn get_shmmr_map(&self) -> PyResult<PyObject> {
        // very expansive as the Rust FxHashMap will be converted to Python's dictionary
        // maybe limit the size that can be converted to avoid OOM
        if let Some(shmmr_to_frags) = self.get_shmmr_map_internal() {
            pyo3::Python::with_gil(|py| Ok(shmmr_to_frags.to_object(py)))
        } else {
            Err(exceptions::PyException::new_err(
                "This method only support FASTX or MEMORY backend.",
            ))
        }
    }

    /// get the ``shmmr_pair`` to ``fragment_id`` map in Python as a list
    ///
    /// this can be very expensive to generate the Python objects of a large hashmap in Rust
    ///
    /// Returns
    /// -------
    ///
    /// list
    ///     list of the tuple (shmmr0, shmmr1, seq_id, position0, position1, orientation)
    ///   
    pub fn get_shmmr_pair_list(&mut self) -> PyResult<Vec<(u64, u64, u32, u32, u32, u8)>> {
        if let Some(shmmr_to_frags) = self.get_shmmr_map_internal() {
            let py_out = shmmr_to_frags
                .par_iter()
                .flat_map(|v| {
                    v.1.iter()
                        .map(|vv| (v.0 .0, v.0 .1, vv.1, vv.2, vv.3, vv.4))
                        .collect::<Vec<(u64, u64, u32, u32, u32, u8)>>()
                })
                .collect::<Vec<(u64, u64, u32, u32, u32, u8)>>();
            Ok(py_out)
        } else {
            Err(exceptions::PyException::new_err(
                "This method only support FASTX or MEMORY backend.",
            ))
        }
    }

    /// fetch a contiguous sub-sequence
    ///
    /// Parameters
    /// ----------
    /// sample_name : string
    ///     the sample name stored in the AGC file
    /// ctg_name : string
    ///     the contig name stored in the AGC file
    /// bgn : int
    ///     the starting coordinate (0-based)
    /// end : int
    ///     the ending coordinate (exclusive)  
    ///
    /// Returns
    /// -------
    /// list
    ///     a list of bytes representing the sequence
    #[pyo3(text_signature = "($self, sample_name, ctg_name, bgn, end)")]
    pub fn get_sub_seq(
        &self,
        sample_name: String,
        ctg_name: String,
        bgn: usize,
        end: usize,
    ) -> PyResult<Vec<u8>> {
        Ok(self
            .db_internal
            .get_sub_seq(sample_name, ctg_name, bgn, end)?)
    }

    /// fetch a contiguous sub-sequence by a sequence id
    ///
    /// Parameters
    /// ----------
    /// sid : int
    ///     sequence id in the database
    ///
    /// bgn : int
    ///     the starting coordinate (0-based)
    /// end : int
    ///     the ending coordinate (exclusive)  
    ///
    /// Returns
    /// -------
    /// list
    ///     a list of bytes representing the sequence
    #[pyo3(text_signature = "($self, sample_name, ctg_name, bgn, end)")]
    pub fn get_sub_seq_by_id(&self, sid: u32, bgn: usize, end: usize) -> PyResult<Vec<u8>> {
        Ok(self.db_internal.get_sub_seq_by_id(sid, bgn, end)?)
    }

    /// fetch a sequence
    ///
    /// Parameters
    /// ----------
    /// sample_name : string
    ///     the sample name stored in the AGC file
    ///
    /// ctg_name : string
    ///     the contig name stored in the AGC file
    ///
    /// Returns
    /// -------
    /// list
    ///     a list of bytes representing the sequence
    #[pyo3(text_signature = "($self, sample_name, ctg_name)")]
    pub fn get_seq(&self, sample_name: String, ctg_name: String) -> PyResult<Vec<u8>> {
        Ok(self.db_internal.get_seq(sample_name, ctg_name)?)
    }

    /// fetch a sequence by the sequence id in the database
    ///
    /// Parameters
    /// ----------
    /// sid : int
    ///     sequence id in the database
    ///
    /// ctg_name : string
    ///     the contig name stored in the AGC file
    ///
    /// Returns
    /// -------
    /// list
    ///     a list of bytes representing the sequence
    #[pyo3(text_signature = "($self, sample_name, ctg_name)")]
    pub fn get_seq_by_id(&self, sid: u32) -> PyResult<Vec<u8>> {
        Ok(self.db_internal.get_seq_by_id(sid).unwrap())
    }

    /// Get adjacent list of the shimmer graph shimmer_pair -> shimmer_pair
    ///
    /// Parameters
    /// ----------
    /// min_count : int
    ///     the minimum number of times a pair of shimmers must be observed to be included in the graph
    ///
    /// Returns
    /// -------
    /// list
    ///     list of pairs of shimmer pairs ((h00, h01, orientation0),(h10, h11, orientation1))  
    ///
    #[pyo3(signature = (min_count, keeps=None))]
    pub fn get_smp_adj_list(
        &self,
        min_count: usize,
        keeps: Option<Vec<u32>>,
    ) -> PyResult<Vec<(u32, (u64, u64, u8), (u64, u64, u8))>> {
        let frag_map = self.get_shmmr_map_internal();
        if frag_map.is_none() {
            return Err(exceptions::PyException::new_err(
                "This method only support FASTX or MEMORY backend.",
            ));
        } else {
            let frag_map = frag_map.unwrap();
            let out = seq_db::frag_map_to_adj_list(frag_map, min_count, keeps)
                .iter()
                .map(|adj_pair| {
                    (
                        adj_pair.0,
                        (adj_pair.1 .0, adj_pair.1 .1, adj_pair.1 .2),
                        (adj_pair.2 .0, adj_pair.2 .1, adj_pair.2 .2),
                    )
                })
                .collect();
            Ok(out)
        }
    }

    /// Sort the adjacent list of the shimmer graph
    ///
    /// Parameters
    /// ----------
    /// adj_list : list
    ///     the list from the output from ``get_smp_adj_list``
    ///
    /// sort_by : (u64, u64)
    ///     the starting node signature
    ///
    /// Returns
    /// -------
    /// list
    ///     list of node in the tuple (node, parent_node, node_weight, is_leaf, global_rank, branch, branch_rank)
    ///
    ///

    pub fn sort_adj_list_by_weighted_dfs(
        &self,
        adj_list: Vec<(u32, (u64, u64, u8), (u64, u64, u8))>,
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
        let adj_list = adj_list
            .iter()
            .map(|&(v0, v1, v2)| {
                (
                    v0,
                    ShmmrGraphNode(v1.0, v1.1, v1.2),
                    ShmmrGraphNode(v2.0, v2.1, v2.2),
                )
            })
            .collect::<AdjList>();

        let start = ShmmrGraphNode(start.0, start.1, start.2);

        if let Some(frag_map) = self.get_shmmr_map_internal() {
            seq_db::sort_adj_list_by_weighted_dfs(&frag_map, &adj_list, start)
                .iter()
                .map(|v| {
                    (
                        (v.0 .0, v.0 .1, v.0 .2),
                        v.1.map(|w| (w.0, w.1, w.2)),
                        v.2,
                        v.3,
                        v.4,
                        v.5,
                        v.6,
                    )
                })
                .collect::<Vec<_>>()
        } else {
            vec![]
        }
    }

    /// Get the principal bundles in MAPG
    ///
    /// Parameters
    /// ----------
    /// min_count : int
    ///     minimum coverage count to be included in the graph
    ///
    /// path_len_cut_off : int
    ///     remove short path less than path_len_cut_off when generating the principal path
    ///     
    ///     if the number is small, the generated principal paths will be more fragmented.
    ///  
    /// Returns
    /// -------
    /// list
    ///     list of paths, each path is a list of nodes
    ///
    #[pyo3(signature = (min_count, path_len_cutoff, keeps=None))]
    pub fn get_principal_bundles(
        &mut self,
        min_count: usize,
        path_len_cutoff: usize,
        keeps: Option<Vec<u32>>,
    ) -> Vec<Vec<(u64, u64, u8)>> {
        let pb = self
            .db_internal
            .get_principal_bundles(min_count, path_len_cutoff, keeps);
        self.principal_bundles = Some((min_count, path_len_cutoff, pb.clone()));
        pb
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

    /// Get the principal bundles and bundle decomposition of all sequences
    ///
    /// Parameters
    /// ----------
    /// min_count : int
    ///     minimum coverage count to be included in the graph
    ///
    /// path_len_cut_off : int
    ///     remove short path less than path_len_cut_off when generating the principal path
    ///     
    ///     if the number is small, the generated principal paths will be more fragmented.
    ///  
    /// Returns
    /// -------
    /// tuple
    ///     a tuple consist of two lists: (principal_bundles, seqid_smps_with_bundle_id_seg_direction)
    ///  
    ///     principal_bundles = list of (principal_bundle_id, ave_bundle_position, list_bundle_vertex)
    ///    
    ///     list_of_bundle_vertex = list of (hash0:u64, hash0:u64, direction:u8)
    ///
    ///     seqid_smps_with_bundle_id_seg_direction = list of shimmer pairs in the database annotated with principal bundle id and direction
    ///     
    ///     the elements of the list are ((hash0:u64, hash1:u64, pos0:u32, pos0:u32, direction:0),
    ///                                   (principal_bundle_id, direction, order_in_the_bundle))
    ///
    #[pyo3(signature = (min_count, path_len_cutoff, keeps=None))]
    pub fn get_principal_bundle_decomposition(
        &mut self,
        min_count: usize,
        path_len_cutoff: usize,
        keeps: Option<Vec<u32>>,
    ) -> (
        Vec<(usize, usize, Vec<(u64, u64, u8)>)>,
        Vec<(
            u32,
            Vec<((u64, u64, u32, u32, u8), Option<(usize, u8, usize)>)>,
        )>,
    ) {
        let pb = self.get_principal_bundles(min_count, path_len_cutoff, keeps);
        //println!("DBG: # bundles {}", pb.len());

        let seqid_seq_list: Vec<(u32, Vec<u8>)> = self
            .db_internal
            .seq_info
            .clone()
            .unwrap_or_default()
            .iter()
            .map(|(sid, data)| {
                let (ctg_name, source, _) = data;
                let source = source.clone().unwrap();
                let seq = self.get_seq(source.clone(), ctg_name.clone()).unwrap();
                (*sid, seq)
            })
            .collect();

        self._get_principal_bundle_projection_internal(pb, seqid_seq_list)
    }

    /// Project sequences outside the sequence database on to a principal bundle decomposition  
    ///
    /// Parameters
    /// ----------
    /// min_count : int
    ///     minimum coverage count to be included in the graph
    ///
    /// path_len_cut_off : int
    ///     remove short path less than path_len_cut_off when generating the principal path
    ///     
    ///     if the number is small, the generated principal paths will be more fragmented.
    ///  
    /// sequences : (contig_id: int, list of sequences)
    ///
    /// Returns
    /// -------
    /// tuple
    ///     a tuple consist of two lists: (principal_bundles, seqid_smps_with_bundle_id_seg_direction)
    ///  
    ///     principal_bundles = list of (principal_bundle_id, ave_bundle_position, list_bundle_vertex)
    ///    
    ///     list_of_bundle_vertex = list of (hash0:u64, hash0:u64, direction:u8)
    ///
    ///     seqid_smps_with_bundle_id_seg_direction = list of shimmer pairs in the database annotated with principal bundle id and direction
    ///     
    ///     the elements of the list are ((hash0:u64, hash1:u64, pos0:u32, pos0:u32, direction:0),
    ///                                   (principal_bundle_id, direction, order_in_the_bundle))
    ///
    ///
    #[pyo3(signature = (min_count, path_len_cutoff, sequence, keeps=None))]
    pub fn get_principal_bundle_projection(
        &mut self,
        min_count: usize,
        path_len_cutoff: usize,
        sequence: Vec<(u32, Vec<u8>)>,
        keeps: Option<Vec<u32>>,
    ) -> (
        Vec<(usize, usize, Vec<(u64, u64, u8)>)>,
        Vec<(
            u32,
            Vec<((u64, u64, u32, u32, u8), Option<(usize, u8, usize)>)>,
        )>,
    ) {
        let pb = self.get_principal_bundles(min_count, path_len_cutoff, keeps);
        //println!("DBG: # bundles {}", pb.len());
        self._get_principal_bundle_projection_internal(pb, sequence)
    }

    fn _get_principal_bundle_projection_internal(
        &self,
        pb: Vec<Vec<(u64, u64, u8)>>,
        sequences: Vec<(u32, Vec<u8>)>,
    ) -> (
        Vec<(usize, usize, Vec<(u64, u64, u8)>)>,
        Vec<(
            u32,
            Vec<((u64, u64, u32, u32, u8), Option<(usize, u8, usize)>)>,
        )>,
    ) {
        fn get_smps(seq: Vec<u8>, shmmr_spec: &ShmmrSpec) -> Vec<(u64, u64, u32, u32, u8)> {
            let shmmrs = sequence_to_shmmrs(0, &seq, &shmmr_spec, false);
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

        //let pb = self.get_principal_bundles(min_count, path_len_cutoff, keeps);
        //println!("DBG: # bundles {}", pb.len());

        let mut vertex_to_bundle_id_direction_pos =
            self.get_vertex_map_from_principal_bundles(pb.clone()); //not efficient but it is PyO3 limit now

        let seqid_smps: Vec<(u32, Vec<(u64, u64, u32, u32, u8)>)> = sequences
            .into_iter()
            .map(|(sid, seq)| {
                (
                    sid as u32,
                    get_smps(seq, &self.db_internal.shmmr_spec.clone().unwrap()),
                )
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
                            .or_insert(vec![])
                            .push(order as f32);
                        bundle_visited.insert(bid.0);
                    }
                    let direction = match bid.1 == v.4 {
                        true => 0,
                        false => 1,
                    };
                    bundle_id_to_directions
                        .entry(bid.0)
                        .or_insert(vec![])
                        .push(direction);
                }
            })
        });

        // determine the bundles' overall orders and directions by consensus voting
        let mut bundle_mean_order_direction = (0..pb.len())
            .into_iter()
            .map(|bid| {
                if let Some(orders) = bundle_id_to_orders.get(&bid) {
                    let sum: f32 = orders.into_iter().sum();
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

        // loop through each sequence and generate the decomposition for the sequence
        let seqid_smps_with_bundle_id_seg_direction = seqid_smps
            .iter()
            .map(|(sid, smps)| {
                let smps = smps
                    .into_iter()
                    .map(|v| {
                        let seg_match =
                            if let Some(m) = vertex_to_bundle_id_direction_pos.get(&(v.0, v.1)) {
                                Some(*m)
                            } else {
                                None
                            };
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

    /// Convert the adjacent list of the shimmer graph shimmer_pair -> GFA
    ///
    /// Parameters
    /// ----------
    /// min_count : int
    ///     the minimum number of times a pair of shimmers must be observed to be included in the graph
    ///
    /// filepath : string
    ///     the path to the output file
    ///
    /// Returns
    /// -------
    ///
    /// None
    ///     The data is written into the file at filepath
    ///
    #[pyo3(signature = (min_count, filepath, method="from_fragmap", keeps=None))]
    pub fn generate_mapg_gfa(
        &self,
        min_count: usize,
        filepath: &str,
        method: &str,
        keeps: Option<Vec<u32>>,
    ) -> PyResult<()> {
        self.db_internal
            .generate_mapg_gfa(min_count, filepath, method, keeps)?;
        Ok(())
    }

    /// Write additional meta data for GFA into a file
    ///
    /// Parameters
    /// ----------
    /// filepath : string
    ///     the path to the output file
    ///
    /// Returns
    /// -------
    ///
    /// None
    ///     The data is written into the file at filepath
    ///

    fn write_mapg_idx(&self, filepath: &str) -> Result<(), std::io::Error> {
        self.db_internal.write_mapg_idx(filepath)?;
        Ok(())
    }

    // for backward compatibility
    pub fn write_midx_to_text_file(&self, filepath: &str) -> Result<(), std::io::Error> {
        self.write_mapg_idx(filepath)
    }

    /// Convert the adjacent list of the shimmer graph shimmer_pair -> GFA
    ///
    /// Parameters
    /// ----------
    /// min_count : int
    ///     the minimum number of times a pair of shimmers must be observed to be included in the graph
    ///
    /// filepath : string
    ///     the path to the output file
    ///
    /// Returns
    /// -------
    ///
    /// None
    ///     The data is written into the file at filepath
    ///     
    #[pyo3(signature = (min_count, path_len_cutoff, filepath, keeps=None))]
    pub fn generate_principal_mapg_gfa(
        &self,
        min_count: usize,
        path_len_cutoff: usize,
        filepath: &str,
        keeps: Option<Vec<u32>>,
    ) -> PyResult<()> {
        self.db_internal.generate_principal_mapg_gfa(
            min_count,
            path_len_cutoff,
            filepath,
            keeps,
        )?;
        Ok(())
    }

    fn write_frag_and_index_files(&self, file_prefix: String) -> () {
        if self.db_internal.seq_db.is_some() {
            let internal = self.db_internal.seq_db.as_ref().unwrap();

            internal.write_to_frag_files(file_prefix.clone(), None);
            internal
                .write_shmmr_map_index(file_prefix.clone())
                .expect("write mdb file fail");
        };
    }

    /// generate consensus sequence for one sequence in the database
    #[pyo3(signature = (sids, min_cov))]
    pub fn shmmr_sparse_aln_consensus(
        &self,
        sids: Vec<u32>,
        min_cov: u32,
    ) -> PyResult<Vec<(u32, Vec<(Vec<u8>, Vec<u32>)>)>> {
        assert!(
            self.db_internal.backend == Backend::FASTX
                || self.db_internal.backend == Backend::MEMORY,
            "Only DB created with load_from_fastx() can add data from another fastx file"
        );
        let sdb = &self.db_internal.seq_db.as_ref().unwrap();
        let consensus = pgr_db::ec::shmmr_sparse_aln_consensus_with_sdb(sids, sdb, min_cov);
        match consensus {
            Ok(seq) => Ok(seq),
            Err(_) => Err(exceptions::PyException::new_err("consensus failed")),
        }
    }
}

impl SeqIndexDB {
    // depending on the storage type, return the corresponded index
    fn get_shmmr_map_internal(&self) -> Option<&seq_db::ShmmrToFrags> {
        match self.db_internal.backend {
            #[cfg(feature = "with_agc")]
            Backend::AGC => None,
            Backend::FASTX => Some(&self.db_internal.seq_db.as_ref().unwrap().frag_map),
            Backend::MEMORY => Some(&self.db_internal.seq_db.as_ref().unwrap().frag_map),
            Backend::FRG => None,
            Backend::UNKNOWN => None,
        }
    }
}

/// A PyO3 class wrapping an existing AGC file for reading
///
/// Example::
///
///      >>> agc_file = AGCFile("/path/to/genomes.agc")
///
#[cfg(feature = "with_agc")]
#[pyclass(unsendable)] // lock in one thread (see https://github.com/PyO3/pyo3/blob/main/guide/src/class.md)
struct AGCFile {
    /// internal agc_file handle
    agc_file: agc_io::AGCFile,

    /// A hashmap mapping (source, ctg_name) to sequence length.
    /// It is more efficient to make a copy in Python to access it to avoid
    /// expensive Rust struct ot Python object conversion than using it directly. For example::
    ///
    ///     >>> agc_file = AGCFile("genomes.agc")
    ///     >>> agc_ctg_lens = agc_file.ctg_lens.copy()
    ///  
    #[pyo3(get)]
    pub ctg_lens: FxHashMap<(String, String), usize>,
}

#[cfg(feature = "with_agc")]
#[pymethods]
impl AGCFile {
    /// constructor
    ///
    /// Parameters
    /// ----------
    /// filepath: string
    ///     the path to a AGC file
    #[pyo3(signature=(filepath))]
    #[new]
    pub fn new(filepath: String) -> PyResult<Self> {
        let agc_file = agc_io::AGCFile::new(filepath)?;
        let mut ctg_lens = FxHashMap::<(String, String), usize>::default();
        agc_file.ctg_lens.iter().for_each(|(k, v)| {
            ctg_lens.insert((k.0.clone(), k.1.clone()), *v);
        });
        Ok(AGCFile { agc_file, ctg_lens })
    }

    /// fetch a contiguous sub-sequence from an AGC file
    ///
    /// Parameters
    /// ----------
    /// sample_name : string
    ///     the sample name stored in the AGC file
    /// ctg_name : string
    ///     the contig name stored in the AGC file
    /// bgn : int
    ///     the starting coordinate (0-based)
    /// end : int
    ///     the ending coordinate (exclusive)  
    ///
    /// Returns
    /// -------
    /// list
    ///     a list of bytes representing the sequence
    #[pyo3(signature = (sample_name, ctg_name, bgn, end))]
    pub fn get_sub_seq(
        &self,
        sample_name: String,
        ctg_name: String,
        bgn: usize,
        end: usize,
    ) -> PyResult<Vec<u8>> {
        Ok(self.agc_file.get_sub_seq(sample_name, ctg_name, bgn, end))
    }

    /// fetch a full contig sequence from an AGC file
    ///
    /// Parameters
    /// ----------
    /// sample_name : string
    ///     the sample name stored in the AGC file
    /// ctg_name : string
    ///     the contig name stored in the AGC file
    ///
    /// Returns
    /// -------
    /// list
    ///     a list of bytes representing the sequence
    #[pyo3(signature = (sample_name, ctg_name))]
    pub fn get_seq(&self, sample_name: String, ctg_name: String) -> PyResult<Vec<u8>> {
        Ok(self.agc_file.get_seq(sample_name, ctg_name))
    }
}

/// Perform sparse dynamic programming to identify alignment between sequence
/// using matched shimmer pairs
///
/// Parameters
/// ----------
/// sp_hits : list
///     a list of tuple of ``HitPair`` defined as
///     ``pub type HitPair = ((u32, u32, u8), (u32, u32, u8))``
///     This represents the hits as the position of matched Shimmer Pairs from
///     the two sequence. For example, if there two shimmers at position 2342
///     and 4322 of the query sequence that matches the shimmers at positions
///     6125465 and 6127445, them the HitPair will be ``(2342, 4322, 0)`` and
///     ``(6125465, 6127445, 0)``. The third number would be 1 if shimmer paired
///     are reversed matched to the sequence orientation and 0 otherwise.
///       
/// max_span : int
///     For a give hit, the max_span defines how many other following hits are
///     considered for the next aligned position. This will limit the search
///     space for the best alignment. If the two sequences are very repetitive,
///     then one needs to use larger ``max_span`` to ensure capturing the right
///     alignment path
//
/// penalty : float
///     this parameter will determine when to break the alignment if there are big
///     gaps between the alignment segment. One can set it to zero to catch large
///     chunk alignment ignoring the gaps. Typically, a number between 0.1 to 0.5 should
///     be used.
///
#[pyfunction(signature = (sp_hits, max_span, penalty))]
pub fn sparse_aln(
    sp_hits: Vec<HitPair>,
    max_span: u32,
    penalty: f32,
) -> PyResult<Vec<(f32, Vec<HitPair>)>> {
    let mut hp = sp_hits.clone();
    Ok(aln::sparse_aln(&mut hp, max_span, penalty))
}

/// Generate a list of shimmer pair from a sequence
///
/// Parameters
/// ----------
/// w : int
///     window size, default to 80, max allowed is 128
///
/// k : int
///     k-mer size, default to 56, max allowed is 56
///
/// r : int
///     reduction factor for generate sparse hierarchical minimizers (shimmer),
///     default to 4, max allowed is 12
///
/// min_span : int
///     - a parameter to remove close-by shimmer pairs
///     - if not zero, shimmer pairs which the distance between them are less
///       than ``the min_span`` will be removed
///     - default to 16
///
/// padding : bool
///     - for short fragment that segmented by using shimmer, set ``padding`` to true
///       to preserve the first and last shimmers
///     - default to false
///
/// Returns
/// -------
/// list of tuple
///     a list fo tuple of ``(shmmr0, shmmr1, position0, position1, orientation)``  
///
#[pyfunction(signature = (seq, w=80, k=56, r=4, min_span=16, padding=false))]
fn get_shmmr_pairs_from_seq(
    seq: Vec<u8>,
    w: u32,
    k: u32,
    r: u32,
    min_span: u32,
    padding: bool,
) -> PyResult<Vec<(u64, u64, u32, u32, u8)>> {
    let shmmr_spec = ShmmrSpec {
        w,
        k,
        r,
        min_span,
        sketch: false,
    };
    let shmmrs = sequence_to_shmmrs(0, &seq, &shmmr_spec, padding);
    let res = seq_db::pair_shmmrs(&shmmrs)
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
        .collect::<Vec<(u64, u64, u32, u32, u8)>>();
    Ok(res)
}

/// Generate a list of shimmer matches for creating a dot plot between two sequences
///
/// Parameters
/// ----------
///
/// seq0 : list
///     a list of bytes representing the first sequences
///
/// seq1 : list
///     a list of bytes representing the second sequences
///
/// w : int
///     window size, default to 80, max allowed is 128
///
/// k : int
///     k-mer size, default to 56, max allowed is 56
///
/// r : int
///     reduction factor for generate sparse hierarchical minimizers (shimmer),
///     default to 4, max allowed is 12
///
/// min_span : int
///     -  a parameter to remove close-by shimmer pairs
///     -  if not zero, shimmer pairs which the distance between them are less
///        than ``the min_span`` will be removed
///
/// Returns
/// -------
/// tuple of two lists
///     ``(x, y)``:
///
///     -  ``x``: the matched shimmer positions in sequence 0
///     -  ``y``: the matched shimmer positions in sequence 1
///
#[pyfunction(signature = (seq0, seq1, w = 80, k = 56, r = 4, min_span = 16))]
fn get_shmmr_dots(
    seq0: Vec<u8>,
    seq1: Vec<u8>,
    w: u32,
    k: u32,
    r: u32,
    min_span: u32,
) -> (Vec<u32>, Vec<u32>) {
    let mut x = Vec::<u32>::new();
    let mut y = Vec::<u32>::new();
    //let seq0v = seq0.to_string_lossy().as_bytes().to_vec();
    //let seq1v = seq1.to_string_lossy().as_bytes().to_vec();
    let shmmr_spec = ShmmrSpec {
        w,
        k,
        r,
        min_span,
        sketch: false,
    };

    let shmmr0 = sequence_to_shmmrs(0, &seq0, &shmmr_spec, false);
    let shmmr1 = sequence_to_shmmrs(1, &seq1, &shmmr_spec, false);
    let mut base_mmer_x = FxHashMap::<u64, Vec<u32>>::default();

    for m in shmmr0 {
        let hash = m.x >> 8;
        let pos = ((m.y & 0xFFFFFFFF) >> 1) as u32;
        base_mmer_x.entry(hash).or_insert_with(|| vec![]).push(pos);
    }

    for m in shmmr1 {
        let hash = m.x >> 8;
        let py = ((m.y & 0xFFFFFFFF) >> 1) as u32;
        if base_mmer_x.contains_key(&hash) {
            for px in base_mmer_x.get(&hash).unwrap() {
                x.push(*px);
                y.push(py);
            }
        }
    }

    (x, y)
}

/// A wrapper class to represent alignment segment for python
///
/// This wraps the Rust struct ``seq2variants::AlnSegment`` mapping
/// the enum ``seq2variants::AlnSegType`` to integers
///
#[pyclass]
#[derive(Clone)]
struct AlnSegment {
    /// alignment type: value =  ``b'M'``, ``b'I'``, ``b'D'``, ``b'X'``, or ``b'?'``
    #[pyo3(get, set)]
    t: u8,
    /// segment coordinate in the reference: (begin, end, length)
    #[pyo3(get, set)]
    ref_loc: (u32, u32, u32),
    /// segment coordinate in the target: (begin, end, length)
    #[pyo3(get, set)]
    tgt_loc: (u32, u32, u32),
}

/// A  class to represent alignment mapping between the two sequences
#[pyclass]
#[derive(Clone)]
pub struct AlnMap {
    /// mapped location between two sequence as a list of paired coordinate
    #[pyo3(get, set)]
    pmap: Vec<(u32, u32)>,
    /// alignment string of the reference
    #[pyo3(get, set)]
    ref_a_seq: Vec<u8>,
    /// alignment string of the target
    #[pyo3(get, set)]
    tgt_a_seq: Vec<u8>,
    /// alignment string of alignment symbols
    #[pyo3(get, set)]
    aln_seq: Vec<u8>,
}

/// Generate the CIGAR string from two sequences with WFA
///
/// Parameters
/// ----------
///
/// seq0 : string
///     a string representing the first sequence
///
/// seq1 : string
///     a string representing the second sequence
///
/// Returns
/// -------
/// tuple
///     tuple of (alignment_score, CIGAR_list)
///
// #[pyfunction(seq0, seq1)]
// #[pyo3(text_signature = "($self, seq0, seq1)")]
// fn get_cigar(seq0: &PyString, seq1: &PyString) -> PyResult<(i32, Vec<u8>)> {
//     if let Ok((score, cigar)) = seqs2variants::get_cigar(&seq0.to_string(), &seq1.to_string()) {
//         Ok((score, cigar))
//     } else {
//         Err(PyValueError::new_err("fail to align"))
//     }
// }

/// Get alignment segments from two sequences
///
/// Parameters
/// ----------
/// ref_id : int
///     a integer id for the reference sequence
///
/// ref_seq : string
///     a python string of the reference sequence
///
/// tgt_id : int
///     a integer id for the target sequence
///
/// tgt_seq : string
///     a python string of the target sequence
///
/// Returns
/// -------
/// list
///     a list of ``AlnSegment``
///
///     the ``AlnSegment`` is a Rust struct defined as::
///
///         pub struct SeqLocus {
///             pub id: u32,
///             pub bgn: u32,
///             pub len: u32,
///         }
///
///         pub enum AlnSegType {
///             Match,
///             Mismatch,
///             Insertion,
///             Deletion,
///             Unspecified,
///         }
///
///         pub struct AlnSegment {
///             pub ref_loc: SeqLocus,
///             pub tgt_loc: SeqLocus,
///             pub t: AlnSegType,
///         }
///
// #[pyfunction(ref_id, ref_seq, tgt_id, tgt_seq)]
// #[pyo3(text_signature = "($self, ref_id, ref_seq, tgt_id, tgt_seq)")]
// fn get_aln_segments(
//     ref_id: u32,
//     ref_seq: &PyString,
//     tgt_id: u32,
//     tgt_seq: &PyString,
// ) -> PyResult<Vec<AlnSegment>> {
//     let ref_seq = ref_seq.to_string();
//     let tgt_seq = tgt_seq.to_string();
//     let aln_segs = seqs2variants::get_aln_segments(ref_id, &ref_seq, tgt_id, &tgt_seq);

//     match aln_segs {
//         Ok(segs) => Ok(segs
//             .par_iter()
//             .map(|seg| {
//                 let t = match seg.t {
//                     seqs2variants::AlnSegType::Match => b'M',
//                     seqs2variants::AlnSegType::Mismatch => b'X',
//                     seqs2variants::AlnSegType::Insertion => b'I',
//                     seqs2variants::AlnSegType::Deletion => b'D',
//                     seqs2variants::AlnSegType::Unspecified => b'?',
//                 };
//                 AlnSegment {
//                     t: t,
//                     ref_loc: (seg.ref_loc.id, seg.ref_loc.bgn, seg.ref_loc.len),
//                     tgt_loc: (seg.tgt_loc.id, seg.tgt_loc.bgn, seg.tgt_loc.len),
//                 }
//             })
//             .collect()),
//         Err(_) => Err(exceptions::PyException::new_err("alignment failed")),
//     }
// }

/// Get alignment map from a list of alignment segments
///
/// Parameters
/// ----------
/// aln_segs : list
///     a list of the ``AlnSegment``
///
/// s0: string
///     a python string of the reference sequence
///
/// s1: int
///     a integer id for the target sequence
///
/// Returns
/// -------
/// list
///     a list of ``AlnSegment``
// #[pyfunction(aln_segs, s0, s1)]
// #[pyo3(text_signature = "($self, aln_segs, s0, s1)")]
// fn get_aln_map(
//     aln_segs: Vec<AlnSegment>,
//     ref_seq: &PyString,
//     tgt_seq: &PyString,
// ) -> PyResult<AlnMap> {
//     let s0 = ref_seq.to_string();
//     let s1 = tgt_seq.to_string();

//     let aln_segs = aln_segs
//         .par_iter()
//         .map(|s| seqs2variants::AlnSegment {
//             ref_loc: seqs2variants::SeqLocus {
//                 id: s.ref_loc.0,
//                 bgn: s.ref_loc.1,
//                 len: s.ref_loc.2,
//             },
//             tgt_loc: seqs2variants::SeqLocus {
//                 id: s.tgt_loc.0,
//                 bgn: s.tgt_loc.1,
//                 len: s.tgt_loc.2,
//             },
//             t: match s.t {
//                 b'M' => seqs2variants::AlnSegType::Match,
//                 b'X' => seqs2variants::AlnSegType::Mismatch,
//                 b'I' => seqs2variants::AlnSegType::Insertion,
//                 b'D' => seqs2variants::AlnSegType::Deletion,
//                 _ => seqs2variants::AlnSegType::Unspecified,
//             },
//         })
//         .collect::<seqs2variants::AlnSegments>();

//     let aln_map = seqs2variants::get_aln_map(&aln_segs, &s0, &s1).unwrap();

//     Ok(AlnMap {
//         pmap: aln_map.pmap,
//         ref_a_seq: aln_map.ref_a_seq,
//         tgt_a_seq: aln_map.tgt_a_seq,
//         aln_seq: aln_map.aln_seq,
//     })
// }

/// Perform a naive de Bruijn graph consensus
///
/// Parameters
/// ----------
/// aln_segs : list
///     a list of the list of bytes representing the bases of each sequence
///
/// kmer_size : int  
///     the size of kmers used for constructing the de Bruijn graph
///
/// min_cov : int
///     to keep haplotype specific consensus, if a kmer has coverage more or equal to min_cov, it will be kept
///
/// Returns
/// -------
/// list
///     a list of bytes representing the consensus sequence
///
#[pyfunction(signature = (seqs, kmer_size=33, min_cov=2))]
pub fn naive_dbg_consensus(
    seqs: Vec<Vec<u8>>,
    kmer_size: usize,
    min_cov: usize,
) -> PyResult<Vec<u8>> {
    let consensus = pgr_db::ec::naive_dbg_consensus(seqs, kmer_size, min_cov);
    match consensus {
        Ok(seq) => Ok(seq),
        Err(_) => Err(exceptions::PyException::new_err(
            "consensus failed, trying bigger kmer size",
        )),
    }
}

/// Perform a shimmer de Bruijn graph consensus
///
/// Parameters
/// ----------
/// aln_segs : list
///     a list of the list of bytes representing the bases of each sequence
///
/// k, w, r, min_span : int  
///     specification of the shimmers for construting graph
///
/// Returns
/// -------
/// list
///     a list of a set of bytes representing the consensus sequences of all branches in the graph
///
#[pyfunction(signature= (seqs, w = 33, k = 33, r = 1, min_span = 0))]
pub fn shmmr_dbg_consensus(
    seqs: Vec<Vec<u8>>,
    w: u32,
    k: u32,
    r: u32,
    min_span: u32,
) -> PyResult<Vec<(Vec<u8>, Vec<u32>)>> {
    let spec = ShmmrSpec {
        w,
        k,
        r,
        min_span,
        sketch: false,
    };
    let consensus = pgr_db::ec::shmmr_dbg_consensus(seqs, &Some(spec));
    match consensus {
        Ok(seq) => Ok(seq),
        Err(_) => Err(exceptions::PyException::new_err(
            "consensus failed, trying bigger kmer size",
        )),
    }
}

/// Perform a guided shimmer de Bruijn graph consensus
///
/// Parameters
/// ----------
/// aln_segs : list
///     a list of the list of bytes representing the bases of each sequence
///
/// k, w, r, min_span : int  
///     specification of the shimmers for construting graph
///
/// min_cov : int
///     to keep hyplotype specific consensus, if a kmer has coverage more or equal to min_cov, it will be kept
///
/// Returns
/// -------
/// list
///     a list of a set of bytes representing the consensus sequences of all branches in the graph
///
#[pyfunction(signature = (seqs, w = 33, k = 33, r = 1, min_span = 0, min_cov = 2))]
pub fn guided_shmmr_dbg_consensus(
    seqs: Vec<Vec<u8>>,
    w: u32,
    k: u32,
    r: u32,
    min_span: u32,
    min_cov: u32,
) -> PyResult<(Vec<u8>, Vec<u32>)> {
    let spec = ShmmrSpec {
        w,
        k,
        r,
        min_span,
        sketch: false,
    };
    let consensus = pgr_db::ec::guided_shmmr_dbg_consensus(seqs, &Some(spec), min_cov);
    match consensus {
        Ok(seq) => Ok(seq),
        Err(_) => Err(exceptions::PyException::new_err(
            "consensus failed, trying bigger kmer size",
        )),
    }
}

/// Perform a shimmer de Bruijn graph consensus
///
/// Parameters
/// ----------
/// aln_segs : list
///     a list of the list of bytes representing the bases of each sequence
///
/// k, w, r, min_span : int  
///     specification of the shimmers for construting graph
///
/// Returns
/// -------
/// list
///     a list of a set of bytes representing the consensus sequences of all branches in the graph
///
#[pyfunction(signature = (seqs, w = 33, k = 33, r = 1, min_span = 0, min_cov = 2))]
pub fn shmmr_sparse_aln_consensus(
    seqs: Vec<Vec<u8>>,
    w: u32,
    k: u32,
    r: u32,
    min_span: u32,
    min_cov: u32,
) -> PyResult<Vec<(Vec<u8>, Vec<u32>)>> {
    let spec = ShmmrSpec {
        w,
        k,
        r,
        min_span,
        sketch: false,
    };
    let consensus = pgr_db::ec::shmmr_sparse_aln_consensus(seqs, &Some(spec), min_cov);
    match consensus {
        Ok(seq) => Ok(seq),
        Err(_) => Err(exceptions::PyException::new_err(
            "consensus failed, trying bigger kmer size",
        )),
    }
}

/// The internal `pgrtk` modules implemented with Rust.
/// These classes and fucntion are re-exported as `pgrtk.*`
/// so `import pgrtk` will bring these classes and function
/// into `pgrtk.*` scope to avoid using the verbose
/// `pgrtk.pgrtk.*`.
#[pymodule]
fn pgrtk(_: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<SeqIndexDB>()?;
    #[cfg(feature = "with_agc")]
    m.add_class::<AGCFile>()?;
    m.add_function(wrap_pyfunction!(sparse_aln, m)?)?;
    m.add_function(wrap_pyfunction!(get_shmmr_dots, m)?)?;
    //m.add_function(wrap_pyfunction!(get_cigar, m)?)?;
    //m.add_function(wrap_pyfunction!(get_aln_segments, m)?)?;
    //m.add_function(wrap_pyfunction!(get_aln_map, m)?)?;
    m.add_function(wrap_pyfunction!(pgr_lib_version, m)?)?;
    m.add_function(wrap_pyfunction!(get_shmmr_pairs_from_seq, m)?)?;
    m.add_function(wrap_pyfunction!(naive_dbg_consensus, m)?)?;
    m.add_function(wrap_pyfunction!(shmmr_dbg_consensus, m)?)?;
    m.add_function(wrap_pyfunction!(guided_shmmr_dbg_consensus, m)?)?;
    m.add_function(wrap_pyfunction!(shmmr_sparse_aln_consensus, m)?)?;
    Ok(())
}
