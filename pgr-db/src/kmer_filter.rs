use cuckoofilter::CuckooFilter;
use std::collections::hash_map::DefaultHasher;

pub struct KmerFilter {
    filter: CuckooFilter<DefaultHasher>,
    kmer_size: usize,
}

impl KmerFilter {
    pub fn new(kmer_size: usize) -> Self {
        let filter = CuckooFilter::new();
        KmerFilter { filter, kmer_size }
    }

    pub fn with_capacity(kmer_size: usize, capacity: usize) -> Self {
        let filter = CuckooFilter::with_capacity(capacity);
        KmerFilter { filter, kmer_size }
    }
}

impl KmerFilter {
    pub fn add_seq(&mut self, seq: &Vec<u8>) {
        (0..seq.len() - self.kmer_size).into_iter().for_each(|pos| {
            self.filter.test_and_add(&seq[pos..pos + self.kmer_size]).unwrap();
        })
    }

    pub fn check_seq(&self, seq: &Vec<u8>) -> usize {
        let mut count = 0_usize;
        (0..seq.len() - self.kmer_size).into_iter().for_each(|pos| {
            if self.filter.contains(&seq[pos..pos + self.kmer_size]) {
                count += 1
            };
        });
        count
    }
}
