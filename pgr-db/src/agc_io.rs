use crate::bindings::{
    agc_close, agc_get_ctg_len, agc_get_ctg_seq, agc_list_ctg, agc_list_destroy, agc_list_sample,
    agc_n_ctg, agc_n_sample, agc_open, agc_t,
};
use libc::strlen;
use pgr_utils::fasta_io::SeqRec;
use std::collections::HashMap;
use std::ffi::CString;
use std::io;
use std::mem;

#[derive(Debug, Clone)]
pub struct AGCSample {
    pub name: String,
    pub contigs: Vec<(String, usize)>, //name, len count
}

#[derive(Debug)]
pub struct AGCFile {
    pub filepath: String,
    agc_handle: *mut agc_t,
    pub samples: Vec<AGCSample>,
    pub ctg_lens: HashMap<(String, String), usize>,
    sample_ctg: Vec<(String, String)>,
    current_ctg: usize,
}

fn cstr_to_string(cstr_ptr: *mut i8) -> String {
    unsafe { String::from_raw_parts(cstr_ptr as *mut u8, strlen(cstr_ptr), strlen(cstr_ptr)) }
}

impl AGCFile {
    pub fn new(filepath: String) -> Self {
        let agc_handle;
        let mut samples = vec![];
        let mut ctg_lens = HashMap::new();
        let mut sample_ctg = vec![];
        unsafe {
            agc_handle = agc_open(CString::new(filepath.clone()).unwrap().into_raw(), 0_i32);
            let mut n_samples = agc_n_sample(agc_handle);
            let samples_ptr: *mut *mut ::std::os::raw::c_char =
                agc_list_sample(agc_handle, &mut n_samples);

            for i in 0..n_samples as usize {
                let s_ptr = *(samples_ptr.add(i));
                let sample_name = cstr_to_string(s_ptr);
                //println!("sample: {}", sample_name);
                let mut n_contig = agc_n_ctg(agc_handle, s_ptr);
                let ctg_ptr = agc_list_ctg(agc_handle, s_ptr, &mut n_contig);
                let mut ctgs: Vec<(String, usize)> = Vec::new();
                for j in 0..n_contig as usize {
                    let c_ptr = *(ctg_ptr.add(j));
                    let ctg_name = cstr_to_string(c_ptr);
                    //println!("ctg: {}", ctg_name);
                    let ctg_len = agc_get_ctg_len(agc_handle, s_ptr, c_ptr);
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
        Self {
            filepath,
            agc_handle,
            samples,
            ctg_lens,
            sample_ctg,
            current_ctg: 0,
        }
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
                self.agc_handle,
                c_sample_name,
                c_ctg_name,
                bgn as i32,
                end as i32 - 1,
                seq_buf,
            );
            seq = String::from_raw_parts(seq_buf as *mut u8, ctg_len - 1, ctg_len);
            //check this, it takes over the pointer? we don't need to free the point manually?
        }
        seq.as_bytes().to_vec()
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
            agc_close(self.agc_handle);
        }
    }
}

impl Iterator for AGCFile {
    type Item = io::Result<SeqRec>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.current_ctg < self.sample_ctg.len() {
            let (sample_name, ctg_name) = self.sample_ctg.get(self.current_ctg).unwrap();
            let id = ctg_name.as_bytes().to_vec();
            let seq = self.get_seq(sample_name.clone(), ctg_name.clone());
            self.current_ctg += 1;
            Some(Ok(SeqRec {
                source: Some(sample_name.clone()),
                id,
                seq,
            }))
        } else {
            None
        }
    }
}
