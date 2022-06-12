#![allow(dead_code)]

use flate2::bufread::MultiGzDecoder;
use std::fs::File;
use std::io::prelude::*;
use std::io::{self, BufReader, SeekFrom};
#[derive(Debug, Clone)]
pub struct SeqRec {
    pub source: Option<String>,
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
}

enum Fastx {
    FastQ,
    FastA,
}
pub struct FastaReader<R> {
    inner: R,
    t: Fastx,
    filename: String,
    seq_capacity: usize,
    keep_source: bool,
}

pub fn reverse_complement(seq: &Vec<u8>) -> Vec<u8> {
    let mut rev_seq = Vec::new();
    for b in seq.iter().rev() {
        match b {
            b'A' => rev_seq.push(b'T'),
            b'C' => rev_seq.push(b'G'),
            b'G' => rev_seq.push(b'C'),
            b'T' => rev_seq.push(b'A'),
            _ => rev_seq.push(*b),
        }
    }
    rev_seq
}

impl<R: BufRead> FastaReader<R> {
    pub fn new(
        mut inner: R,
        filename: &String,
        seq_capacity: usize,
        keep_source: bool,
    ) -> Result<Self, io::Error> {
        let t: Fastx;
        {
            let r = inner.by_ref();
            let mut buf = Vec::<u8>::new();
            r.take(1).read_to_end(&mut buf)?;
            if buf.len() < 1 {
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    format!("empty file: {}", filename),
                ));
            }
            match buf[0] {
                b'>' => t = Fastx::FastA,
                b'@' => t = Fastx::FastQ,
                _ => t = Fastx::FastA,
            }
        }
        Ok(Self {
            inner,
            t,
            filename: filename.to_string(),
            seq_capacity,
            keep_source,
        })
    }

    pub fn next_rec(&mut self) -> Option<io::Result<SeqRec>> {
        match self.t {
            Fastx::FastA => self.fasta_next_rec(),
            Fastx::FastQ => self.fastq_next_rec(),
        }
    }

    pub fn fasta_next_rec(&mut self) -> Option<io::Result<SeqRec>> {
        let mut id_tmp = Vec::<u8>::with_capacity(128);
        let mut seq = Vec::<u8>::with_capacity(self.seq_capacity);

        let res = self.inner.read_until(b'\n', &mut id_tmp);
        if res.is_err() {
            Some(res);
        } else if res.ok() == Some(0) {
            return None;
        }
        let mut r = BufReader::new(&id_tmp[..]);
        let mut id = Vec::<u8>::with_capacity(128);
        let _res = r.read_until(b' ', &mut id);

        let id = id
            .into_iter()
            .filter(|c| *c != b'\n' && *c != b' ' && *c != b'\r')
            .collect::<Vec<u8>>();
        let _x = self.inner.read_until(b'>', &mut seq);
        let mut seq = seq
            .drain(..)
            .filter(|c| *c != b'\n' && *c != b'>' && *c != b'\r')
            .collect::<Vec<u8>>();
        if seq.capacity() as f32 > seq.len() as f32 * 1.2 {
            seq.shrink_to_fit();
        }
        let source = if self.keep_source {Some(self.filename.to_string())} else {None};
        let rec = SeqRec {
            source,
            id: id,
            seq: seq,
        };

        Some(Ok(rec))
    }

    pub fn fastq_next_rec(&mut self) -> Option<io::Result<SeqRec>> {
        let mut id_tmp = Vec::<u8>::with_capacity(128);
        let mut seq = Vec::<u8>::with_capacity(self.seq_capacity);

        let _res = self.inner.read_until(b'\n', &mut id_tmp); //read id
        // fetch the first id up to the first space, strip '\n'
        let mut r = BufReader::new(&id_tmp[..]);
        let mut id = Vec::<u8>::with_capacity(128);
        let _res = r.read_until(b' ', &mut id);
        let id = id
            .into_iter()
            .filter(|c| *c != b'\n' && *c != b' ' && *c != b'\r')
            .collect::<Vec<u8>>();
        // get the seq
        let _res = self.inner.read_until(b'\n', &mut seq);

        let mut seq = seq
            .drain(..)
            .filter(|c| *c != b'\n' && *c != b'\r')
            .collect::<Vec<u8>>();
        if seq.capacity() as f32 > seq.len() as f32 * 1.2 {
            seq.shrink_to_fit();
        }
        let source = if self.keep_source {Some(self.filename.to_string())} else {None};
        let rec = SeqRec {
            source,
            id: id,
            seq: seq,
        };
        // ignore QV
        let mut buf = Vec::<u8>::with_capacity(1024);
        let _res = self.inner.read_until(b'+', &mut buf);
        let _res = self.inner.read_until(b'\n', &mut buf);
        let _res = self.inner.read_until(b'\n', &mut buf);
        let res = self.inner.read_until(b'@', &mut buf); //get to id line
        if res.is_err() {
            Some(res);
        } else if res.ok() == Some(0) {
            return None;
        }
        Some(Ok(rec))
    }
}

impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = io::Result<SeqRec>;
    fn next(&mut self) -> Option<Self::Item> {
        self.next_rec()
    }
}


pub struct FastqStreamReader {
    inner: std::io::Stdin,
    seq_capacity: usize 
}

impl FastqStreamReader {
    pub fn new(seq_capacity: usize) -> Self {
        FastqStreamReader {
            inner: std::io::stdin(),
            seq_capacity
        } 
    }
}

impl Iterator for FastqStreamReader {
    type Item = io::Result<SeqRec>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut tmp = String::with_capacity(128);
        match self.inner.read_line(&mut tmp) {
            Ok(n) => {
                if n == 0 {
                    return None
                }
                let tmp = tmp.trim();
                if &tmp[0..1] == "@" {
                    let id = tmp[1..].as_bytes().to_vec();
                    let mut seq = String::with_capacity(self.seq_capacity);
                    if self.inner.read_line(&mut seq).unwrap_or(0) == 0 {
                        return None
                    }
                    let seq = seq.trim();
                    let seq = seq[..].as_bytes().to_vec();
                    let mut buf = String::with_capacity(self.seq_capacity);
                    if self.inner.read_line(&mut buf).unwrap_or(0) == 0 {
                        return None
                    };
                    if self.inner.read_line(&mut buf).unwrap_or(0) == 0 {
                        return None
                    };
                    let source = None;
                    let rec = SeqRec {
                        source,
                        id: id,
                        seq: seq,
                    }; 
                    Some(Ok(rec))
                } else {
                    None
                }
            },
            Err(_) => {
                None
            }
        }
    }
}






fn encode_biseq(s: Vec<u8>) -> Vec<u8> {
    let fourbit_map_f: [u8; 256] = [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ];
    let fourbit_map_r: [u8; 256] = [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 8, 0, 4, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 8, 0, 4, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ];
    let len = s.len();
    let mut out_s = Vec::<u8>::with_capacity(len);
    for p in 0..len {
        let rp = len - 1 - p;
        let code = (fourbit_map_r[s[rp] as usize] << 4) | fourbit_map_f[s[p] as usize];
        out_s.push(code);
    }
    out_s
}

pub fn build(seq_list_file: &String, out_prefix: &String) -> Result<usize, io::Error> {
    let seqdb_name = format!("{}.seqdb", out_prefix);
    log::info!("create seq db: {}", seqdb_name);
    let mut out_db_file = File::create(seqdb_name)?;
    let seqidx_name = format!("{}.idx", out_prefix);
    log::info!("create seq index: {}", seqidx_name);
    let mut out_idx_file = File::create(seqidx_name)?;
    let mut start = 0_usize;
    let mut seq_id = 0_u32;

    log::info!("get input files from: {}", seq_list_file);
    let f = File::open(seq_list_file)?;
    let seq_list_buf = BufReader::new(f);

    for fastx_file in seq_list_buf.lines() {
        let input_fn = fastx_file.unwrap();
        log::info!("input file: {}", input_fn);
        let metadata = std::fs::metadata(&input_fn)?;
        if !metadata.is_file() || metadata.len() < (1 << 16) {
            log::info!(
                "input file: {} may not be proper input file (filesize = {}), ignore",
                input_fn,
                metadata.len()
            );
            continue;
        }
        let input_file = File::open(&input_fn)?;
        let mut reader = BufReader::new(input_file);
        let mut is_gzfile = false;
        {
            let r = reader.by_ref();
            let mut buf = Vec::<u8>::new();
            let _ = r.take(2).read_to_end(&mut buf);
            if buf == [0x1F_u8, 0x8B_u8] {
                log::info!("input file: {} detected as gz-compressed file", input_fn);
                is_gzfile = true;
            }
        }

        let _ = reader.seek(SeekFrom::Start(0));
        if is_gzfile {
            let fastx_buf = BufReader::new(MultiGzDecoder::new(&mut reader));
            let mut fastx_reader = FastaReader::new(fastx_buf, &input_fn, 1 << 14, true)?;
            while let Some(r) = fastx_reader.next_rec() {
                let r = r.unwrap();
                if r.seq.len() < 500 {
                    //ignore very short reads
                    continue;
                }
                //println!("{}", String::from_utf8_lossy(&r.id));
                //println!("{}", String::from_utf8_lossy(&r.seq));
                let biseq = encode_biseq(r.seq);
                let _ = out_db_file.write(&biseq);
                let _ = writeln!(
                    out_idx_file,
                    "{:09} {} {} {}",
                    seq_id,
                    String::from_utf8_lossy(&r.id),
                    biseq.len(),
                    start
                );
                start += biseq.len();
                seq_id += 1;
            }
        } else {
            let mut fastx_reader = FastaReader::new(reader, &input_fn, 1 << 14, true)?;
            while let Some(r) = fastx_reader.next_rec() {
                let r = r.unwrap();
                if r.seq.len() < 500 {
                    //ignore very short reads
                    continue;
                }
                //println!("{}", String::from_utf8_lossy(&r.id));
                //println!("{}", String::from_utf8_lossy(&r.seq));
                let biseq = encode_biseq(r.seq);
                let _ = out_db_file.write(&biseq);
                let _ = writeln!(
                    out_idx_file,
                    "{:09} {} {} {}",
                    seq_id,
                    String::from_utf8_lossy(&r.id),
                    biseq.len(),
                    start
                );
                start += biseq.len();
                seq_id += 1;
            }
        }
    }
    log::info!("total number of reads indexed: {}", seq_id);
    log::info!("total number of bases: {}", start);
    log::info!("average read length: {}", start as f32 / seq_id as f32);
    Ok(start)
}
