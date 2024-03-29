use cuckoofilter::CuckooFilter;
use rustc_hash::FxHashSet;
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
        (0..seq.len() - self.kmer_size).for_each(|pos| {
            self.filter.test_and_add(&seq[pos..pos + self.kmer_size]).unwrap();
        })
    }

    pub fn check_seq(&self, seq: &Vec<u8>) -> usize {
        let mut count = 0_usize;
        (0..seq.len() - self.kmer_size).for_each(|pos| {
            if self.filter.contains(&seq[pos..pos + self.kmer_size]) {
                count += 1
            };
        });
        count
    }
    
    pub fn add_seq_mmers(&mut self, seq: &Vec<u8>) {
        let k = self.kmer_size as u32;
        let w = k >> 1;
        let shmmrs = crate::shmmrutils::sequence_to_shmmrs1(0, seq, w, k, 1, 0, false);
        shmmrs.into_iter().for_each(|mmer| {
            self.filter.test_and_add(&mmer.x).unwrap();
        })
    }

    pub fn check_seq_mmers(&self, seq: &Vec<u8>) -> (usize, usize) {
        let mut count = 0_usize;
        let k = self.kmer_size as u32;
        let w = k >> 1;
        let shmmrs = crate::shmmrutils::sequence_to_shmmrs1(0, seq, w, k, 1, 0, false);
        shmmrs.iter().for_each(|mmer| {
            if self.filter.contains(&mmer.x) {
                count += 1
            };
        });
        (shmmrs.len(), count)
    }
}

pub struct MinimizerFilter {
    filter: FxHashSet<u64>,
    kmer_size: usize,

}

impl MinimizerFilter {
    pub fn new(kmer_size: usize) -> Self {
        let filter = FxHashSet::default();
        MinimizerFilter { filter, kmer_size }
    }
}

impl MinimizerFilter {
    
    pub fn add_seq_mmers(&mut self, seq: &Vec<u8>) {
        let k = self.kmer_size as u32;
        let w = k >> 1;
        let shmmrs = crate::shmmrutils::sequence_to_shmmrs1(0, seq, w, k, 1, 0, false);
        shmmrs.into_iter().for_each(|mmer| {
            self.filter.insert(mmer.x);
        })
    }

    pub fn check_seq_mmers(&self, seq: &Vec<u8>) -> (usize, usize) {
        let mut count = 0_usize;
        let k = self.kmer_size as u32;
        let w = k >> 1;
        let shmmrs = crate::shmmrutils::sequence_to_shmmrs1(0, seq, w, k, 1, 0, false);
        shmmrs.iter().for_each(|mmer| {
            if self.filter.contains(&mmer.x) {
                count += 1
            };
        });
        (shmmrs.len(), count)
    }
}
