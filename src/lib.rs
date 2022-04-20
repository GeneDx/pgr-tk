pub mod cseq_db;
pub mod fasta_io;
pub mod shmmrutils;

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::{BufRead, BufReader, Read};
    use std::collections::HashMap;
    use crate::fasta_io::{FastaReader, reverse_complement};
    use flate2::bufread::MultiGzDecoder;

    use crate::cseq_db;
    use crate::cseq_db::{Fragment, KMERSIZE};
    

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
                log::info!(
                    "input file detected as gz-compressed file",
                );
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
        };
        seqs
    }

    #[test]
    fn load_seq_test() {
        let seqs = load_seqs();
        let mut csdb = cseq_db::CompressedSeqDB::new("test/test_data/test_seqs.fa".to_string());
        let _ = csdb.load_seqs();
        print!("test");
        for seq in csdb.seqs.iter() {
            println!();
            println!("{}", seq.name);
            for shmr in seq.shmmrs.iter() {
                println!("{}", shmr);
            }
            let mut reconstruct_seq = <Vec<u8>>::new();
            let mut p = 0;
            for frg_id in seq.seq_frags.iter() {
                println!("{}:{}", frg_id, csdb.frags[*frg_id as usize]);
                match csdb.frags.get(*frg_id as usize).unwrap() {
                    Fragment::Prefix(b) => {
                        reconstruct_seq.extend_from_slice(&b[..]);
                        println!("p: {} {}", p, p+b.len());
                        p += b.len();
                    }
                    Fragment::Suffix(b) => {
                        reconstruct_seq.extend_from_slice(&b[..]);
                        println!("p: {} {}", p, p+b.len());
                        p += b.len();
                    }
                    Fragment::Internal(b) => {
                        reconstruct_seq.extend_from_slice(&b[KMERSIZE as usize..]);
                        println!("p: {} {}", p, p+b.len());
                        p += b.len();
                    }
                    Fragment::AlnSegments((frg_id, reverse, a)) => {
                        if let Fragment::Internal(base_seq) = csdb.frags.get(*frg_id as usize).unwrap()
                        {
                            let mut bs = base_seq.clone();
                            if *reverse == true {
                                bs = reverse_complement(&bs);
                            }
                            let seq = cseq_db::reconstruct_seq_from_aln_segs(&bs, a);
                            reconstruct_seq.extend_from_slice(&seq[KMERSIZE as usize..]);
                        println!("p: {} {}", p, p+seq.len());
                        p += seq.len();
                        }
                    }
                }
            }
            let orig_seq = seqs.get(&seq.name).unwrap();
            if reconstruct_seq != *orig_seq  {
                println!("{}", seq.name);
                println!("{:?}", reconstruct_seq);
                println!("{:?}", orig_seq);
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
}
