pub mod fasta_io;
pub mod gff_db;
pub mod multi_seqmap;
pub mod read_db;
pub mod sag;
pub mod seq_db;
pub mod seqmap;
pub mod seqs2variants;
pub mod shmmrutils;
pub use core::mem::MaybeUninit;
pub use libc::{getrusage, rusage, RUSAGE_SELF, RUSAGE_THREAD};

#[derive(Copy, Clone)]
pub struct Parameters {
    pub nthreads: u32,
    pub nchunks: u32,
    pub k: u32,
    pub w: u32,
    pub r: u32,
    pub tol: f64,
}

#[allow(dead_code)]
pub fn log_resource(msg: &str, data: &mut rusage) -> (u64, u64, u64) {
    let _res = unsafe { getrusage(RUSAGE_SELF, data) };
    log::info!(
        "{} : (maxRSS, utime, stime): {} {} {}",
        msg,
        data.ru_maxrss,
        data.ru_utime.tv_sec,
        data.ru_stime.tv_sec
    );

    (
        data.ru_maxrss as u64,
        data.ru_utime.tv_sec as u64,
        data.ru_stime.tv_sec as u64,
    )
}
