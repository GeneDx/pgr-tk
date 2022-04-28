
pub mod agc_io;
pub mod seq_db;
pub mod shmmrutils;
pub mod bindings;
pub mod aln;

#[cfg(test)]
mod tests {
    use pgr_utils::fasta_io::{reverse_complement, FastaReader};
    use crate::shmmrutils::{match_reads, DeltaPoint};
    use flate2::bufread::MultiGzDecoder;
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::{BufRead, BufReader, Read};

    use crate::seq_db::{self, deltas_to_aln_segs, reconstruct_seq_from_aln_segs};
    use crate::seq_db::{Fragment, KMERSIZE};

    pub fn load_seqs() -> HashMap<String, Vec<u8>> {
        let mut seqs = HashMap::<String, Vec<u8>>::new();
        let filepath = "test/test_data/test_seqs.fa";
        let file = File::open(filepath.to_string()).unwrap();
        let mut reader = BufReader::new(file);
        let mut is_gzfile = false;
        {
            let r = reader.by_ref();
            let mut buf = Vec::<u8>::new();
            let _ = r.take(2).read_to_end(&mut buf);
            if buf == [0x1F_u8, 0x8B_u8] {
                log::info!("input file detected as gz-compressed file",);
                is_gzfile = true;
            }
        }
        drop(reader);

        let file = File::open(&filepath).unwrap();
        let mut reader = BufReader::new(file);
        let gz_buf = &mut BufReader::new(MultiGzDecoder::new(&mut reader));

        let file = File::open(&filepath).unwrap();
        let reader = BufReader::new(file);
        let std_buf = &mut BufReader::new(reader);

        let fastx_buf: &mut dyn BufRead = if is_gzfile {
            drop(std_buf);
            gz_buf
        } else {
            drop(gz_buf);
            std_buf
        };

        let mut fastx_reader = FastaReader::new(fastx_buf, &filepath.to_string()).unwrap();
        while let Some(rec) = fastx_reader.next_rec() {
            let rec = rec.unwrap();
            let seqname = String::from_utf8_lossy(&rec.id).into_owned();
            seqs.insert(seqname, rec.seq.clone());
        }
        seqs
    }

    #[test]
    pub fn gz_file_read_test() {
        let mut sdb = seq_db::CompactSeqDB::new();
        let _ = sdb.load_seqs("test/test_data/test_seqs2.fa.gz".to_string());
        println!("{:?}", sdb.seqs[0].seq_frags);
    }

    #[test]
    fn load_seq_test() {
        let seqs = load_seqs();
        let mut sdb = seq_db::CompactSeqDB::new();
        let _ = sdb.load_seqs("test/test_data/test_seqs2.fa.gz".to_string());
        //println!("test");
        for seq in sdb.seqs.iter() {
            //println!();
            //println!("{}", seq.name);
            let mut reconstruct_seq = <Vec<u8>>::new();
            let mut _p = 0;
            for frg_id in seq.seq_frags.iter() {
                //println!("{}:{}", frg_id, sdb.frags[*frg_id as usize]);
                match sdb.frags.get(*frg_id as usize).unwrap() {
                    Fragment::Prefix(b) => {
                        reconstruct_seq.extend_from_slice(&b[..]);
                        //println!("p: {} {}", p, p + b.len());
                        _p += b.len();
                    }
                    Fragment::Suffix(b) => {
                        reconstruct_seq.extend_from_slice(&b[..]);
                        //println!("p: {} {}", p, p + b.len());
                        _p += b.len();
                    }
                    Fragment::Internal(b) => {
                        reconstruct_seq.extend_from_slice(&b[KMERSIZE as usize..]);
                        //println!("p: {} {}", p, p + b.len());
                        _p += b.len();
                    }
                    Fragment::AlnSegments((frg_id, reverse, a)) => {
                        if let Fragment::Internal(base_seq) =
                            sdb.frags.get(*frg_id as usize).unwrap()
                        {
                            let bs = base_seq.clone();
                            let mut seq = seq_db::reconstruct_seq_from_aln_segs(&bs, a);
                            if *reverse == true {
                                seq = reverse_complement(&seq);
                            }
                            reconstruct_seq.extend_from_slice(&seq[KMERSIZE as usize..]);
                            //println!("p: {} {}", p, p + seq.len());
                            _p += seq.len();
                        }
                    }
                }
            }
            let orig_seq = seqs.get(&seq.name).unwrap();
            if reconstruct_seq != *orig_seq {
                //println!("{}", seq.name);
                //println!("{:?}", reconstruct_seq);
                //println!("{:?}", orig_seq);
                for i in 0..reconstruct_seq.len() {
                    if orig_seq[i] != reconstruct_seq[i] {
                        println!("{} {} {} X", i, orig_seq[i], reconstruct_seq[i]);
                    } else {
                        println!("{} {} {}  ", i, orig_seq[i], reconstruct_seq[i]);
                    }
                }
            };
            assert_eq!(reconstruct_seq, *orig_seq);
        }
    }

    #[test]
    fn reconstruct_test1() {
        let base_frg = "TATTTATATTTATTTATATATATTTATATATTTATATATATATTTATATATAAATAT"
            .as_bytes()
            .to_vec();
        let frg = "TTTTTATTTTTTTAATTAATTAATTATTTATTTATTTATTTATTTATTTATTTATTT"
            .as_bytes()
            .to_vec();
        //let frg = "TTATATTTATTTATATATATTTATATAGTTTATATATATATTTATATATAAATATATA".as_bytes().to_vec();
        let m = match_reads(&base_frg, &frg, true, 0.1, 0, 0, 32);
        if let Some(m) = m {
            let deltas: Vec<DeltaPoint> = m.deltas.unwrap();
            let aln_segs = deltas_to_aln_segs(&deltas, m.end0 as usize, m.end1 as usize, &frg);
            let re_seq = reconstruct_seq_from_aln_segs(&base_frg, &aln_segs);
            if frg != re_seq || true {
                println!("{} {}", String::from_utf8_lossy(&base_frg), base_frg.len());
                println!("{} {}", String::from_utf8_lossy(&frg), frg.len());
                println!("{} {} {} {}", m.bgn0, m.end0, m.bgn1, m.end1);
                println!("{:?}", deltas);
                println!(
                    "{}",
                    String::from_utf8_lossy(&reconstruct_seq_from_aln_segs(&base_frg, &aln_segs))
                );
                println!("{:?}", aln_segs);
            }
            assert_eq!(frg, reconstruct_seq_from_aln_segs(&base_frg, &aln_segs));
        }
    }

    #[test]
    fn reconstruct_test2() {
        let base_frg = "TATTTATATTTATTTATATATATTTATATATTTATATATATATTTATATATAAATAT"
            .as_bytes()
            .to_vec();
        let frg = "TTTTTTATTTTTTTAATTAATTAATTATTTATTTATTTATTTATTTATTTATTTATT"
            .as_bytes()
            .to_vec();
        //let frg = "TTATATTTATTTATATATATTTATATAGTTTATATATATATTTATATATAAATATATA".as_bytes().to_vec();
        let m = match_reads(&base_frg, &frg, true, 0.1, 0, 0, 32);
        if let Some(m) = m {
            let deltas: Vec<DeltaPoint> = m.deltas.unwrap();
            let aln_segs = deltas_to_aln_segs(&deltas, m.end0 as usize, m.end1 as usize, &frg);
            let re_seq = reconstruct_seq_from_aln_segs(&base_frg, &aln_segs);
            if frg != re_seq || true {
                println!("{} {}", String::from_utf8_lossy(&base_frg), base_frg.len());
                println!("{} {}", String::from_utf8_lossy(&frg), frg.len());
                println!("{} {} {} {}", m.bgn0, m.end0, m.bgn1, m.end1);
                println!("{:?}", deltas);
                println!(
                    "{}",
                    String::from_utf8_lossy(&reconstruct_seq_from_aln_segs(&base_frg, &aln_segs))
                );
                println!("{:?}", aln_segs);
            }
            assert_eq!(frg, reconstruct_seq_from_aln_segs(&base_frg, &aln_segs));
        }
    }
    #[test]
    fn rc_match() {
        let mut sdb = seq_db::CompactSeqDB::new();
        let _ = sdb.load_seqs("test/test_data/test_rev.fa".to_string());
        let cs0 = sdb.seqs.get(0).unwrap();
        let cs1 = sdb.seqs.get(1).unwrap();
        let shmmr0 = cs0.shmmrs.iter().map(|m| m.x >> 8).collect::<Vec<u64>>();
        let shmmr1 = cs1
            .shmmrs
            .iter()
            .rev()
            .map(|m| m.x >> 8)
            .collect::<Vec<u64>>();
        assert!(shmmr0.len() > 0);
        assert_eq!(shmmr0, shmmr1);
    }
    #[test]
    fn raw_agc_test() {
        use crate::bindings::{
            agc_get_ctg_len, agc_get_ctg_seq, agc_list_ctg, agc_list_sample, agc_n_ctg,
            agc_n_sample, agc_open,
        };
        use libc::strlen;

        use std::ffi::CString;
        use std::{mem, slice, str};
        let c: i32 = 0_i32;
        unsafe {
            let agc_file = agc_open(
                CString::new("test/test_data/test.agc").unwrap().into_raw(),
                c,
            );
            let n_samples = agc_n_sample(agc_file);
            println!("agc_test n_sample: {:?}", n_samples);
            let mut n_samples: i32 = n_samples;
            let samples: *mut *mut ::std::os::raw::c_char =
                agc_list_sample(agc_file, &mut n_samples);
            let sample = *(samples.add(1));
            println!(
                "agc_test list_sample: {:?}",
                str::from_utf8_unchecked(slice::from_raw_parts(
                    sample as *const u8,
                    strlen(sample) + 1
                ))
            );
            let mut n_contig = agc_n_ctg(agc_file, sample);
            println!("agc_test n_contig: {:?}", n_contig);
            let contigs: *mut *mut ::std::os::raw::c_char =
                agc_list_ctg(agc_file, sample, &mut n_contig);
            let ctg = *(contigs.add(1));
            println!(
                "agc_test list_ctg: {:?}",
                str::from_utf8_unchecked(slice::from_raw_parts(ctg as *const u8, strlen(ctg) + 1))
            );
            let ctg_len = agc_get_ctg_len(agc_file, sample, ctg);
            println!("agc_test ctg_len: {:?}", ctg_len);
            let seq_buf: *mut i8 = libc::malloc(mem::size_of::<i8>() * ctg_len as usize) as *mut i8;
            agc_get_ctg_seq(agc_file, sample, ctg, 0, ctg_len, seq_buf);
            let seq = str::from_utf8_unchecked(slice::from_raw_parts(
                seq_buf as *const u8,
                strlen(seq_buf),
            ));
            println!("agc_test seq: {}", seq);
        }
    }

    #[test]
    fn act_io_test() {
        use crate::agc_io::AGCFile;

        let agcfile = AGCFile::new(String::from("test/test_data/test.agc"));
        let seq = agcfile.get_sub_seq(
            "test_agc_ref".to_string(),
            "NA21309#1#JAHEPC010000026.1:3279880-3319873".to_string(),
            500,
            1000,
        );
        assert!(seq.len() == 500);
        //println!("seq_read_test: {}", String::from_utf8_lossy(&seq[..]));

        agcfile.into_iter().for_each(|s| {
            if let Ok(r) = s {
                println!("{} {}", String::from_utf8_lossy(&r.id[..]), r.seq.len());
            }
        });

        let agcfile = AGCFile::new(String::from("test/test_data/test.agc"));
        let mut sdb = seq_db::CompactSeqDB::new();
        let _ = sdb.load_index_from_agcfile(agcfile);
        println!("index size: {}", sdb.frag_map.len());
    }


    #[test]
    fn query_frag_test() {
        use crate::agc_io::AGCFile;
        use seq_db::{query_fragment};
        let agcfile = AGCFile::new(String::from("test/test_data/test.agc"));
        
        let mut sdb = seq_db::CompactSeqDB::new();
        let _ = sdb.load_index_from_agcfile(agcfile);

        let mut agcfile = AGCFile::new(String::from("test/test_data/test.agc"));
        let seq = agcfile.next();
        let r_frags = query_fragment(&sdb.frag_map,&seq.unwrap().unwrap().seq);
        let mut out = vec![];
        for res in r_frags {
            for v in res.2 {
                //println!("Q {:?} {:?} {:?}", res.0, res.1, v);
                out.push( (v, res.1, res.0))
            }
        }
        out.sort();
        for v in out {
            println!("Q {:?}", v);
        }
        //write_shmr_map_bincode(&sdb.frag_map, "test_shmmr.db".to_string());


    }

    #[test]
    fn test_shmmrmap_read_write () -> Result<(), std::io::Error> {
        use crate::agc_io::AGCFile;
        use seq_db::{ query_fragment, write_shmr_map_file, read_shmr_map_file};
        let agcfile = AGCFile::new(String::from("test/test_data/test.agc"));
        let mut sdb = seq_db::CompactSeqDB::new();
        let _ = sdb.load_index_from_agcfile(agcfile);
        write_shmr_map_file(&sdb.shmmr_spec, &sdb.frag_map, "test/test_data/test_shmmr.db".to_string())?;
        let (_shmmr_spec, new_map) = read_shmr_map_file("test/test_data/test_shmmr.db".to_string()).unwrap();

        let mut agcfile = AGCFile::new(String::from("test/test_data/test.agc"));
        let seq = agcfile.next();
        let r_frags = query_fragment(&new_map,&seq.unwrap().unwrap().seq);
        let mut out = vec![];
        for res in r_frags {
            for v in res.2 {
                //println!("Q {:?} {:?} {:?}", res.0, res.1, v);
                out.push( (v, res.1, res.0))
            }
        }
        out.sort();
        for v in out {
            println!("Q {:?}", v);
        }
        Ok(())
    }
}
