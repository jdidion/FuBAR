use bair::ch1::{ch1_3::*, ch1_4::*};
use criterion::{criterion_group, criterion_main, Criterion};
use lazy_static::lazy_static;
use pulp::Arch;
use reqwest::blocking::get;
use std::io::{copy, BufWriter};

const GENOME_URL: &str =
    "http://bioinformaticsalgorithms.com/data/realdatasets/Replication/Vibrio_cholerae.txt";
const GENOME_SIZE: usize = 1_108_250;

lazy_static! {
    static ref TEST_GENOME: Vec<u8> = {
        let mut resp = get(GENOME_URL).expect("request failed");
        let mut buf = BufWriter::new(Vec::with_capacity(GENOME_SIZE));
        copy(&mut resp, &mut buf).expect("failed to copy content");
        buf.into_inner().unwrap()
    };
}

fn benchmark_find_pattern(c: &mut Criterion) {
    let pattern: &[u8] = b"CTTGATCAT";
    c.bench_function("find_pattern", |b| {
        b.iter(|| find_pattern(TEST_GENOME.as_ref(), pattern));
    });
    c.bench_function("find_pattern_func", |b| {
        b.iter(|| find_pattern_func(TEST_GENOME.as_ref(), pattern));
    });
    let pattern_str = "CTTGATCAT";
    c.bench_function("find_pattern_re", |b| {
        b.iter(|| find_pattern_re(TEST_GENOME.as_ref(), pattern_str));
    });
}

fn benchmark_reverse_complement(c: &mut Criterion) {
    let mut seq = Seq::new(b"ATGATCAAG");
    c.bench_function("rc", |b| b.iter(|| seq.reverse_complement()));
    let mut seq = Seq::new(b"ATGATCAAG");
    let arch = Arch::new();
    c.bench_function("rc_simd", |b| b.iter(|| seq.reverse_complement_simd(&arch)));
}

fn benchmark_find_clumps(c: &mut Criterion) {
    let k: usize = 9;
    let window_size: usize = 500;
    let min_clump_size: usize = 3;
    c.bench_function("fast", |b| {
        b.iter(|| find_clumps_fast(TEST_GENOME.as_ref(), k, window_size, min_clump_size))
    });
    c.bench_function("optimal", |b| {
        b.iter(|| find_clumps_optimal(TEST_GENOME.as_ref(), k, window_size, min_clump_size))
    });
}

criterion_group!(
    benches,
    benchmark_find_pattern,
    benchmark_reverse_complement,
    benchmark_find_clumps
);
criterion_main!(benches);
