use flate2::bufread::MultiGzDecoder;
use pgr_db::agc_io::AGCFile;
use pgr_db::fasta_io::FastaReader;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

use pgr_db::seq_db::{self, query_fragment, read_mdb_file};

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

    let mut fastx_reader =
        FastaReader::new(fastx_buf, &filepath.to_string(), 1 << 14, true).unwrap();
    while let Some(rec) = fastx_reader.next_rec() {
        let rec = rec.unwrap();
        let seqname = String::from_utf8_lossy(&rec.id).into_owned();
        seqs.insert(seqname, rec.seq.clone());
    }
    seqs
}

fn _load_seq_test() {
    let seqs = load_seqs();
    let mut sdb = seq_db::CompactSeqDB::new(seq_db::SHMMRSPEC);
    let _shmmr_spec = &pgr_db::seq_db::SHMMRSPEC;
    let _ = sdb.load_seqs_from_fastx(
        "/wd/peregrine-r-ext/phasing_test/PanMHCgraph/HPRCy1.MHC.fa".to_string(),
    );
    //println!("test");
    for seq in sdb.seqs.iter() {
        println!("S {} {} {}", seq.name, seq.id, seq.len);
        //println!();
        //println!("{}", seq.name);
        let reconstruct_seq = sdb.get_seq(&seq);
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

fn _load_index_from_fastx() -> Result<(), std::io::Error> {
    let mut sdb = seq_db::CompactSeqDB::new(seq_db::SHMMRSPEC);
    let filelist = File::open("./filelist").unwrap();

    BufReader::new(filelist).lines().into_iter().for_each(|fp| {
        let fp = fp.unwrap();
        let _ = sdb.load_index_from_fastx(fp);
    });

    seq_db::write_shmr_map_file(&sdb.shmmr_spec, &sdb.frag_map, "test.db".to_string())?;

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
    Ok(())
}

fn load_index_from_agcfile() -> Result<(), std::io::Error> {
    let mut sdb = seq_db::CompactSeqDB::new(seq_db::SHMMRSPEC);
    let filelist = File::open("./filelist").unwrap();

    BufReader::new(filelist).lines().into_iter().try_for_each(
        |fp| -> Result<(), std::io::Error> {
            let fp = fp.unwrap();
            let agcfile = AGCFile::new(fp)?;
            let _ = sdb.load_index_from_agcfile(agcfile);
            Ok(())
        },
    )?;

    //seq_db::write_shmr_map_file(&sdb.frag_map, "test.db".to_string());
    sdb.write_shmr_map_index("test".to_string())?;
    Ok(())
}

fn _load_index_mdb() -> Result<(), std::io::Error> {
    let agcfile = AGCFile::new(String::from("grch38.agc"))?;
    for sample in agcfile.samples.iter() {
        for contig in sample.contigs.iter() {
            let (_n, _t) = contig;
            //println!("{}:{}:{}", sample.name, n, t);
        }
    }
    let seq_mhc = agcfile.get_sub_seq("GCA_000001405.15_GRCh38_no_alt_analysis_set".to_string(), 
    "chr6  AC:CM000668.2  gi:568336018  LN:170805979  rl:Chromosome  M5:5691468a67c7e7a7b5f2a3a683792c29  AS:GRCh38".to_string(), 
    28510120, 33480577);
    // println!("MHC seq len: {}", MHCseq.len());
    let (_shmmr_spec, new_map) = read_mdb_file("test.db".to_string()).unwrap();
    let shmmr_spec = &pgr_db::seq_db::SHMMRSPEC;
    let r_frags = query_fragment(&new_map, &seq_mhc, shmmr_spec);
    let mut out = vec![];
    for res in r_frags {
        for v in res.2 {
            //println!("Q {:?} {:?} {:?}", res.0, res.1, v);
            out.push((v, res.1, res.0))
        }
    }
    out.sort();
    for (v0, v1, _) in out {
        println!(
            "Q {} {} {} {} {} {} {} {}",
            v0.0, v0.1, v0.2, v0.3, v0.4, v1.0, v1.1, v1.2
        );
    }
    Ok(())
}

fn main() -> Result<(), std::io::Error> {
    //load_seq_test();
    //load_index_from_fastx();
    load_index_from_agcfile()?;
    //load_index_mdb();
    Ok(())
}
