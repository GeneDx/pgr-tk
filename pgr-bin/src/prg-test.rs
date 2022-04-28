use flate2::bufread::MultiGzDecoder;
use pgr_db::agc_io::AGCFile;
use pgr_utils::fasta_io::{reverse_complement, FastaReader};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

use pgr_db::seq_db::{self, Fragment, KMERSIZE, read_shmr_map_file, query_fragment};

pub fn load_seqs() -> HashMap<String, Vec<u8>> {
    let mut seqs = HashMap::<String, Vec<u8>>::new();
    //let filepath = "test/test_data/test_seqs.fa";
    let filepath = "/wd/peregrine-r-ext/phasing_test/PanMHCgraph/HPRCy1.MHC.fa";
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

fn load_seq_test() {
    let seqs = load_seqs();
    let mut sdb = seq_db::CompressedSeqDB::new();
    let _ =
        sdb.load_seqs("/wd/peregrine-r-ext/phasing_test/PanMHCgraph/HPRCy1.MHC.fa".to_string());
    //println!("test");
    for seq in sdb.seqs.iter() {
        println!("S {} {} {}", seq.name, seq.id, seq.len);
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
                    if let Fragment::Internal(base_seq) = sdb.frags.get(*frg_id as usize).unwrap()
                    {
                        let bs = base_seq.clone();
                        let mut seq = seq_db::reconstruct_seq_from_aln_segs(&bs, a);
                        if *reverse == true {
                            seq = reverse_complement(&seq);
                        }
                        assert!(base_seq.len() > KMERSIZE as usize);
                        if seq.len() < KMERSIZE as usize {
                            println!("{} {:?} {:?}", seq.len(), seq, a);
                        }
                        assert!(seq.len() > KMERSIZE as usize);
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
        } else {
            println!("{} matched", seq.name);
        };
        assert_eq!(reconstruct_seq, *orig_seq);
    }
    for (shmmr_pair, frg_ids) in sdb.frag_map.into_iter() {
        for ids in frg_ids {
            println!(
                "M {:016X} {:016X} {} {} {} {} {}",
                shmmr_pair.0, shmmr_pair.1, ids.0, ids.1, ids.2, ids.3, ids.4
            );
        }
    }
}

fn load_index_from_fastx() -> Result<(), std::io::Error> {
    let mut sdb = seq_db::CompressedSeqDB::new();
    let filelist = File::open("./filelist").unwrap();

    BufReader::new(filelist).lines().into_iter().for_each(|fp| {
        let fp = fp.unwrap();
        let _ = sdb.load_index_from_fastx(fp);
    });

    seq_db::write_shmr_map_file(&sdb.frag_map, "test.db".to_string())?;

    for seq in sdb.seqs.iter() {
        println!("S {} {} {}", seq.name, seq.id, seq.len);
    }
    for (shmmr_pair, frg_ids) in sdb.frag_map.into_iter() {
        for ids in frg_ids {
            println!(
                "M {:016X} {:016X} {} {} {} {} {}",
                shmmr_pair.0, shmmr_pair.1, ids.0, ids.1, ids.2, ids.3, ids.4
            );
        }
    };
    Ok(())
}


fn load_index_from_agcfile() {
    let mut sdb = seq_db::CompressedSeqDB::new();
    let filelist = File::open("./filelist").unwrap();

    BufReader::new(filelist).lines().into_iter().for_each(|fp| {
        let fp = fp.unwrap();
        let agcfile = AGCFile::new(fp);
        let _ = sdb.load_index_from_agcfile(agcfile);
    });

    seq_db::write_shmr_map_file(&sdb.frag_map, "test.db".to_string());

    for seq in sdb.seqs.iter() {
        println!("S {} {} {}", seq.name, seq.id, seq.len);
    }
    for (shmmr_pair, frg_ids) in sdb.frag_map.into_iter() {
        for ids in frg_ids {
            println!(
                "M {:016X} {:016X} {} {} {} {} {}",
                shmmr_pair.0, shmmr_pair.1, ids.0, ids.1, ids.2, ids.3, ids.4
            );
        }
    }
}

fn load_index_sb() {
    let agcfile = AGCFile::new(String::from("grch38.agc"));
    for sample in agcfile.samples.iter() {
        for contig in sample.contigs.iter() {
            let (n, t) = contig;
            //println!("{}:{}:{}", sample.name, n, t);
        }
    }
    let seq_mhc = agcfile.get_sub_seq("GCA_000001405.15_GRCh38_no_alt_analysis_set".to_string(), 
    "chr6  AC:CM000668.2  gi:568336018  LN:170805979  rl:Chromosome  M5:5691468a67c7e7a7b5f2a3a683792c29  AS:GRCh38".to_string(), 
    28510120, 33480577);
    // println!("MHC seq len: {}", MHCseq.len());
    let new_map = read_shmr_map_file("test.db".to_string()).unwrap();
    
    let r_frags = query_fragment(&new_map,&seq_mhc);
    let mut out = vec![];
    for res in r_frags {
        for v in res.2 {
            //println!("Q {:?} {:?} {:?}", res.0, res.1, v);
            out.push( (v, res.1, res.0))
        }
    }
    out.sort();
    for (v0, v1, _) in out {
        println!("Q {} {} {} {} {} {} {} {}", v0.0, v0.1, v0.2, v0.3, v0.4, v1.0, v1.1, v1.2);
    }
}

fn main() {
    //load_seq_test();
    //load_index_from_fastx();
    //load_index_from_agcfile();
    load_index_sb();
}
