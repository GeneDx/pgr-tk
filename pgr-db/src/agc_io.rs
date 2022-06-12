use crate::bindings::{
    agc_close, agc_get_ctg_len, agc_get_ctg_seq, agc_list_ctg, agc_list_destroy, agc_list_sample,
    agc_n_ctg, agc_n_sample, agc_open, agc_t,
};
use crate::fasta_io::SeqRec;
use libc::strlen;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;
use rayon::ThreadPool;
use rayon::ThreadPoolBuilder;
use std::cell::RefCell;
use std::collections::HashMap;
use std::ffi::CString;
use std::io;
use std::mem;
//use std::sync::Arc;

#[derive(Debug, Clone)]
struct AGCHandle(*mut agc_t);

unsafe impl Send for AGCHandle {}
unsafe impl Sync for AGCHandle {}

#[derive(Debug, Clone)]
pub struct AGCSample {
    pub name: String,
    pub contigs: Vec<(String, usize)>, //name, len count
}

#[derive(Debug, Clone)]
pub struct AGCFile {
    pub filepath: String,
    agc_handle: AGCHandle,
    pub samples: Vec<AGCSample>,
    pub ctg_lens: HashMap<(String, String), usize>,
    sample_ctg: Vec<(String, String)>,
}

pub struct AGCFileIter<'a> {
    agc_file: &'a AGCFile,
    agc_thread_pool: ThreadPool,
    current_ctg: usize,
    seq_buf: RefCell<Option<(usize, usize, Vec<SeqRec>)>>,
}

fn cstr_to_string(cstr_ptr: *mut i8) -> String {
    unsafe { String::from_raw_parts(cstr_ptr as *mut u8, strlen(cstr_ptr), strlen(cstr_ptr)) }
}

impl AGCFile {
    pub fn new(filepath: String) -> Result<Self, std::io::Error> {
        if !std::path::Path::new(&filepath).exists() {
            return Err(std::io::Error::new(std::io::ErrorKind::NotFound, filepath));
        }

        let agc_handle;
        let mut samples = vec![];
        let mut ctg_lens = HashMap::new();
        let mut sample_ctg = vec![];

        unsafe {
            agc_handle = AGCHandle(agc_open(
                CString::new(filepath.clone()).unwrap().into_raw(),
                0_i32,
            ));
            let mut n_samples = agc_n_sample(agc_handle.0);
            let samples_ptr: *mut *mut ::std::os::raw::c_char =
                agc_list_sample(agc_handle.0, &mut n_samples);

            for i in 0..n_samples as usize {
                let s_ptr = *(samples_ptr.add(i));
                let sample_name = cstr_to_string(s_ptr);
                //println!("sample: {}", sample_name);
                let mut n_contig = agc_n_ctg(agc_handle.0, s_ptr);
                let ctg_ptr = agc_list_ctg(agc_handle.0, s_ptr, &mut n_contig);
                let mut ctgs: Vec<(String, usize)> = Vec::new();
                for j in 0..n_contig as usize {
                    let c_ptr = *(ctg_ptr.add(j));
                    let ctg_name = cstr_to_string(c_ptr);
                    //println!("ctg: {}", ctg_name);
                    let ctg_len = agc_get_ctg_len(agc_handle.0, s_ptr, c_ptr);
                    ctg_lens.insert((sample_name.clone(), ctg_name.clone()), ctg_len as usize);
                    sample_ctg.push((sample_name.clone(), ctg_name.clone()));
                    ctgs.push((ctg_name, ctg_len as usize));
                }
                agc_list_destroy(ctg_ptr);
                samples.push(AGCSample {
                    name: sample_name,
                    contigs: ctgs,
                });
            }
            agc_list_destroy(samples_ptr);
        }
        Ok(Self {
            filepath,
            agc_handle,
            samples,
            ctg_lens,
            sample_ctg,
        })
    }

    pub fn get_sub_seq(
        &self,
        sample_name: String,
        ctg_name: String,
        bgn: usize,
        end: usize,
    ) -> Vec<u8> {
        let key = (sample_name.clone(), ctg_name.clone());
        assert!(self.ctg_lens.contains_key(&key));
        assert!(*self.ctg_lens.get(&key).unwrap() >= end);
        assert!(*self.ctg_lens.get(&key).unwrap() >= bgn);
        assert!(bgn < end);

        let c_sample_name: *mut i8 = CString::new(sample_name).unwrap().into_raw();
        let c_ctg_name: *mut i8 = CString::new(ctg_name).unwrap().into_raw();
        let seq;
        let ctg_len = end - bgn + 1;

        unsafe {
            let seq_buf: *mut i8 = libc::malloc(mem::size_of::<i8>() * ctg_len as usize) as *mut i8;
            agc_get_ctg_seq(
                self.agc_handle.0,
                c_sample_name,
                c_ctg_name,
                bgn as i32,
                end as i32 - 1,
                seq_buf,
            );
            seq = <Vec<u8>>::from_raw_parts(seq_buf as *mut u8, ctg_len - 1, ctg_len);
            //check this, it takes over the pointer? we don't need to free the point manually?
        }
        seq
    }

    pub fn get_seq(&self, sample_name: String, ctg_name: String) -> Vec<u8> {
        let key = (sample_name.clone(), ctg_name.clone());
        assert!(self.ctg_lens.contains_key(&key));
        let bgn = 0;
        let end = *self.ctg_lens.get(&key).unwrap();
        let seq = self.get_sub_seq(sample_name, ctg_name, bgn, end);
        assert!(seq.len() == end - bgn);
        seq
    }
}

impl Drop for AGCFile {
    fn drop(&mut self) {
        unsafe {
            agc_close(self.agc_handle.0);
        }
    }
}

impl<'a> IntoIterator for &'a AGCFile {
    // can we parallelized this?
    type Item = io::Result<SeqRec>;
    type IntoIter = AGCFileIter<'a>;
    fn into_iter(self) -> Self::IntoIter {
        AGCFileIter::new(self)
    }
}

impl<'a> AGCFileIter<'a> {
    pub fn new(agc_file: &'a AGCFile) -> Self {
        let agc_thread_pool = ThreadPoolBuilder::new().num_threads(36).build().unwrap();
        AGCFileIter {
            agc_file,
            agc_thread_pool,
            current_ctg: 0,
            seq_buf: RefCell::new(None),
        }
    }
}

thread_local! {
    pub static TL_AGCHANDLE: RefCell<Option<AGCFile>> = RefCell::new(None);
}

impl<'a> Iterator for AGCFileIter<'a> {
    // can we parallelized this?
    type Item = io::Result<SeqRec>;

    fn next(&mut self) -> Option<Self::Item> {
        let number_decoder = 128_usize;

        if self.current_ctg == self.agc_file.sample_ctg.len() {
            return None;
        }

        if self.seq_buf.borrow().is_none() {
            self.seq_buf.replace(Some((0, 0, vec![])));
        }

        let seq_buf: (usize, usize, Vec<SeqRec>) = self.seq_buf.take().unwrap();

        let buf_e = seq_buf.1;
        self.seq_buf.replace(Some(seq_buf));

        if self.current_ctg == buf_e {
            // buffer exhausted
            let mut next_batch = <Vec<(String, String, String, usize, usize)>>::new();
            for i in self.current_ctg..self.current_ctg + number_decoder {
                if i == self.agc_file.sample_ctg.len() {
                    break;
                }
                let (sample_name, ctg_name) = self.agc_file.sample_ctg.get(i).unwrap();
                let bgn = 0;
                let end = *self
                    .agc_file
                    .ctg_lens
                    .get(&(sample_name.clone(), ctg_name.clone()))
                    .unwrap();
                next_batch.push((
                    self.agc_file.filepath.clone(),
                    sample_name.clone(),
                    ctg_name.clone(),
                    bgn,
                    end,
                ));
            }

            let mut seq_buf: (usize, usize, Vec<SeqRec>) = self.seq_buf.take().unwrap();

            let v_seq_rec = self.agc_thread_pool.install(|| {
                let seq_buf = next_batch
                    .par_iter()
                    .map(|(filepath, s, c, _bgn, _end)| {
                        TL_AGCHANDLE.with(|tl_agc_handle| {
                            if (*tl_agc_handle.borrow_mut()).is_none() {
                                *tl_agc_handle.borrow_mut() =
                                    Some(AGCFile::new(filepath.clone()).unwrap());
                            }
                            let t = tl_agc_handle.borrow_mut();
                            let agc_handle = (t.as_ref()).unwrap();
                            let seq = agc_handle.get_seq(s.clone(), c.clone());
                            SeqRec {
                                source: Some(s.clone()),
                                id: c.as_bytes().to_vec(),
                                seq: seq,
                            }
                        })
                        //let seq = self.get_seq(s.clone(), c.clone());
                    })
                    .collect::<Vec<SeqRec>>();
                seq_buf
            });
            seq_buf.2 = v_seq_rec;

            seq_buf.0 = self.current_ctg;
            seq_buf.1 = self.current_ctg + seq_buf.2.len();
            self.seq_buf.replace(Some(seq_buf));
        }

        let seq_buf: (usize, usize, Vec<SeqRec>) = self.seq_buf.take().unwrap();
        let rtn = Some(Ok(seq_buf.2[self.current_ctg - seq_buf.0].clone()));
        self.seq_buf.replace(Some(seq_buf));
        self.current_ctg += 1;
        rtn
    }
}
