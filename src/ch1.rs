pub mod ch1_2 {
    use pcre2::bytes::Regex;
    use std::{collections::HashMap, convert::Into};

    pub fn pattern_count(text: &[u8], pattern: &[u8]) -> usize {
        // assert!(
        //     text.len() >= pattern.len(),
        //     "text length must be >= pattern length"
        // );
        let pattern_len = pattern.len();
        if text.len() < pattern_len {
            return 0;
        }
        let mut matches = 0;
        for i in 0..(text.len() - pattern_len) {
            if &text[i..(i + pattern_len)] == pattern {
                matches += 1;
            }
        }
        matches
    }

    pub fn pattern_count_func<'a, T: Into<&'a [u8]>, P: Into<&'a [u8]>>(
        text: T,
        pattern: P,
    ) -> usize {
        let text: &[u8] = text.into();
        let pattern: &[u8] = pattern.into();
        if text.len() < pattern.len() {
            return 0;
        }
        text.windows(pattern.len())
            .enumerate()
            .filter(|(_, kmer)| *kmer == pattern)
            .map(|(i, _)| i)
            .count()
    }

    pub fn pattern_count_re<'a, T: Into<&'a [u8]>, P: Into<&'a str>>(text: T, pattern: P) -> usize {
        let text: &[u8] = text.into();
        let pattern: &str = pattern.into();
        if text.len() < pattern.len() {
            return 0;
        }
        let re = Regex::new(format!("(?={})", pattern).as_str()).unwrap();
        re.find_iter(text).count()
    }

    // Builds a table of kmer -> Vec(offset)
    pub fn offset_table<'a, T: Into<&'a [u8]>>(text: T, k: usize) -> HashMap<&'a [u8], Vec<usize>> {
        let text: &[u8] = text.into();
        if text.len() < k {
            return HashMap::new();
        }
        if text.len() == k {
            return HashMap::from([(text, vec![0])]);
        }
        let mut offsets: HashMap<&'a [u8], Vec<usize>> = HashMap::new();
        for (i, kmer) in text.windows(k).enumerate() {
            offsets
                .entry(kmer)
                .and_modify(|e| e.push(i))
                .or_insert(Vec::from([i]));
        }
        offsets
    }

    pub fn most_frequent_kmers<'a, T: Into<&'a [u8]>>(text: T, k: usize) -> (Vec<&'a [u8]>, usize) {
        let counter = offset_table(text, k);
        if counter.is_empty() {
            return (Vec::new(), 0);
        }
        let max_frequency = counter.values().map(|v| v.len()).max().unwrap();
        let most_frequent = counter
            .into_iter()
            .filter(|(_, v)| v.len() == max_frequency)
            .map(|(k, _)| k)
            .collect();
        (most_frequent, max_frequency)
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use lazy_static::lazy_static;

        const TEXT: &[u8] = b"ACAACTATGCATACTATCGGGAACTATCCT";
        const PATTERN: &[u8] = b"ACTAT";
        const PATTERN_STR: &str = "ACTAT";
        const EXPECTED_COUNT: usize = 3;

        #[test]
        fn test_pattern_count() {
            assert_eq!(pattern_count(TEXT, PATTERN), EXPECTED_COUNT);
        }

        #[test]
        fn test_pattern_count_func() {
            assert_eq!(pattern_count_func(TEXT, PATTERN), EXPECTED_COUNT);
        }

        #[test]
        fn test_pattern_count_re() {
            assert_eq!(pattern_count_re(TEXT, PATTERN_STR), EXPECTED_COUNT);
        }

        const K1: usize = 5;
        const TEXT2: &[u8] = b"CGATATATCCATAG";
        const K2: usize = 3;
        const TEXT3: &[u8] = b"ACACACA";
        const K3: usize = 2;
        const TEXT_ORI: &[u8] =
            b"atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaacctgagtggatgacatcaagatagg\
            tcgttgtatctccttcctctcgtactctcatgaccacggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaata\
            cttgtgacttgtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggattacgaaagcatgatcat\
            ggctgttgttctgtttatcttgttttgactgagacttgttaggatagacggtttttcatcactgactagccaaagccttactctgc\
            ctgacatcgaccgtaaattgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaagatcttcaa\
            ttgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtttccttaaccctctattttttacggaagaatgatc\
            aagctgctgctcttgatcatcgtttc";

        lazy_static! {
            static ref EXPECTED1: (Vec<&'static [u8]>, usize) = (vec![b"ACTAT"], 3);
            static ref EXPECTED2: (Vec<&'static [u8]>, usize) = (vec![b"ATA"], 3);
            static ref EXPECTED3: (Vec<&'static [u8]>, usize) = (vec![b"AC", b"CA"], 3);
            static ref EXPECTED_ORI: [(usize, (Vec<&'static [u8]>, usize)); 7] = [
                (3, (vec![b"tga"], 25)),
                (4, (vec![b"atga"], 12)),
                (5, (vec![b"gatca", b"tgatc"], 8)),
                (6, (vec![b"tgatca"], 8)),
                (7, (vec![b"atgatca"], 5)),
                (8, (vec![b"atgatcaa"], 4)),
                (
                    9,
                    (
                        vec![b"atgatcaag", b"ctcttgatc", b"cttgatcat", b"tcttgatca"],
                        3
                    )
                )
            ];
        }

        #[test]
        fn test_most_frequent_kmers() {
            assert_eq!(most_frequent_kmers(TEXT, K1), *EXPECTED1);
            assert_eq!(most_frequent_kmers(TEXT2, K2), *EXPECTED2);
            let mut result3 = most_frequent_kmers(TEXT3, K3);
            result3.0.sort();
            assert_eq!(result3, *EXPECTED3);
        }

        #[test]
        fn test_most_frequent_kmers_ori() {
            for (k, expected) in (*EXPECTED_ORI).iter() {
                let mut result = most_frequent_kmers(TEXT_ORI, *k);
                result.0.sort();
                assert_eq!(result, *expected);
            }
        }
    }
}

pub mod ch1_3 {
    use lazy_static::lazy_static;
    use pcre2::bytes::Regex;
    use pulp::Arch;
    use std::collections::HashMap;

    lazy_static! {
        static ref COMPLEMENT: HashMap<u8, u8> =
            HashMap::from([(b'A', b'T'), (b'C', b'G'), (b'G', b'C'), (b'T', b'A')]);
    }

    /// Complement a mutable byte. Assumes that `c` is one of `{A,C,G,T,N}`. Uses XOR to flip the
    /// bits that differ between A<->T (21) or C<->G (4). `N` is its own complement.
    fn complement_mut(c: &mut u8) {
        let val = *c;
        if val != b'N' {
            if val & 2 != 0 {
                *c = val ^ 4
            } else {
                *c = val ^ 21
            }
        }
    }

    #[derive(Debug, PartialEq, Eq)]
    pub struct Seq(Vec<u8>);

    impl Seq {
        pub fn new(bytes: &[u8]) -> Self {
            Seq(Vec::from(bytes))
        }

        pub fn from_str(s: &str) -> Self {
            Seq(Vec::from(s.as_bytes()))
        }

        /// Reverse-complements this sequence
        pub fn reverse_complement(&mut self) {
            self.0.reverse();
            self.0.iter_mut().for_each(complement_mut);
        }

        // Note: currently pulp only supports x86 intrinsics, so this has no effect on ARM
        pub fn reverse_complement_simd(&mut self, arch: &Arch) {
            arch.dispatch(|| {
                self.0.reverse();
                self.0.iter_mut().for_each(complement_mut);
            });
        }

        pub fn as_reverse_complement(&self) -> Self {
            Seq(self
                .0
                .iter()
                .rev()
                .map(|nuc| *COMPLEMENT.get(nuc).unwrap())
                .collect())
        }
    }

    pub fn find_pattern(text: &[u8], pattern: &[u8]) -> Vec<usize> {
        let pattern_len = pattern.len();
        if text.len() < pattern_len {
            return Vec::new();
        }
        let mut matches: Vec<usize> = Vec::new();
        for i in 0..(text.len() - pattern_len) {
            if &text[i..(i + pattern_len)] == pattern {
                matches.push(i);
            }
        }
        matches
    }

    pub fn find_pattern_func<'a, T: Into<&'a [u8]>, P: Into<&'a [u8]>>(
        text: T,
        pattern: P,
    ) -> Vec<usize> {
        let text: &[u8] = text.into();
        let pattern: &[u8] = pattern.into();
        if text.len() < pattern.len() {
            return Vec::new();
        }
        text.windows(pattern.len())
            .enumerate()
            .filter(|(_, kmer)| *kmer == pattern)
            .map(|(i, _)| i)
            .collect()
    }

    pub fn find_pattern_re<'a, T: Into<&'a [u8]>, P: Into<&'a str>>(
        text: T,
        pattern: P,
    ) -> Vec<usize> {
        let text: &[u8] = text.into();
        let pattern: &str = pattern.into();
        if text.len() < pattern.len() {
            return Vec::new();
        }
        let re = Regex::new(format!("(?={})", pattern).as_str()).unwrap();
        re.find_iter(text).map(|m| m.unwrap().start()).collect()
    }

    #[cfg(test)]
    mod tests {
        use super::{super::tests::TEST_GENOME, *};

        const SEQ1: &[u8] = b"ATGATCAAG";
        const SEQ2: &[u8] = b"CTTGATCAT";
        const PATTERN: &[u8] = b"ATGATCAAG";
        const POSITIONS: [usize; 17] = [
            116556, 149355, 151913, 152013, 152394, 186189, 194276, 200076, 224527, 307692, 479770,
            610980, 653338, 679985, 768828, 878903, 985368,
        ];

        #[test]
        fn test_revcomp() {
            let seq1 = Seq::new(SEQ1);
            let seq2 = Seq::new(SEQ2);
            assert_eq!(seq1.as_reverse_complement(), seq2);
        }

        #[test]
        fn test_find_pattern() {
            let positions = find_pattern(TEST_GENOME.as_ref(), PATTERN);
            assert_eq!(positions, POSITIONS);
        }
    }
}

pub mod ch1_4 {
    use super::ch1_2::offset_table;
    use std::collections::{HashMap, VecDeque};

    fn find_offset_clumps(
        offsets: Vec<usize>,
        k: usize,
        window_size: usize,
        min_clump_size: usize,
    ) -> Vec<(usize, usize)> {
        let max_end_offset = window_size - k;
        let mut clumps = Vec::new();
        for i in 0..(offsets.len() - min_clump_size) {
            let start = offsets[i];
            for j in i + 1..offsets.len() {
                if offsets[j] - start > max_end_offset {
                    let clump_size = j - i;
                    if clump_size >= min_clump_size {
                        clumps.push((start, clump_size));
                    }
                    break;
                }
            }
        }
        clumps
    }

    // T = text size
    // K = number of unique kmers
    // W = window size
    // Assumption: T >> K

    /// Computational complexity: O(T)
    /// Memory complexity: O(T)
    pub fn find_clumps_fast<'a, T: Into<&'a [u8]>>(
        text: T,
        k: usize,
        window_size: usize,
        min_clump_size: usize,
    ) -> HashMap<&'a [u8], Vec<(usize, usize)>> {
        let kmers = offset_table(text, k);
        kmers
            .into_iter()
            // ignore infrequent kmers
            .filter(|(_, v)| v.len() >= min_clump_size)
            // look for clumps in the offsets
            .map(|(kmer, v)| (kmer, find_offset_clumps(v, k, window_size, min_clump_size)))
            // filter out kmers with no clumps
            .filter(|(_, v)| !v.is_empty())
            .collect()
    }

    /// Computational complexity: O(T * W)
    /// Memory complexity: O(K)
    pub fn find_clumps_frugal<'a, T: Into<&'a [u8]>>(
        text: T,
        k: usize,
        window_size: usize,
        min_clump_size: usize,
    ) -> HashMap<&'a [u8], Vec<(usize, usize)>> {
        let text: &[u8] = text.into();
        let mut clumps: HashMap<&'a [u8], Vec<(usize, usize)>> = HashMap::new();
        for window in text.windows(window_size) {
            let kmers = offset_table(window, k);
            for (kmer, v) in kmers.into_iter() {
                if v.len() >= min_clump_size {
                    clumps
                        .entry(kmer)
                        .and_modify(|c| c.push((v[0], v.len())))
                        .or_insert(Vec::new());
                }
            }
        }
        clumps
    }

    // Computational complexity: O(T)
    /// Memory complexity: O(K)
    pub fn find_clumps_optimal<'a, T: Into<&'a [u8]>>(
        text: T,
        k: usize,
        window_size: usize,
        min_clump_size: usize,
    ) -> HashMap<&'a [u8], Vec<(usize, usize)>> {
        let text: &[u8] = text.into();
        let max_end_offset = window_size - k;
        // for each kmer, store the current clump and the completed clumps
        let mut clumps: HashMap<&'a [u8], (VecDeque<usize>, Vec<(usize, usize)>)> = HashMap::new();
        for (i, kmer) in text.windows(k).enumerate() {
            clumps
                .entry(kmer)
                .and_modify(|(cur_clump, found_clumps)| {
                    // add a separate clump for each set of offsets where the first and the last
                    // offset are within `window_size` of each other and the number of offsets
                    // is >= `min_clump_size`
                    while !cur_clump.is_empty() && i - cur_clump.front().unwrap() > max_end_offset {
                        let clump_size = cur_clump.len();
                        let start = cur_clump.pop_front().unwrap();
                        if clump_size >= min_clump_size {
                            found_clumps.push((start, clump_size))
                        }
                    }
                    cur_clump.push_back(i);
                })
                .or_insert((VecDeque::from([i]), Vec::new()));
        }
        clumps
            .into_iter()
            .map(|(kmer, (cur_clump, mut found_clumps))| {
                // add the last clump if it is long enough
                if cur_clump.len() >= min_clump_size {
                    found_clumps.push((*cur_clump.front().unwrap(), cur_clump.len()))
                }
                (kmer, found_clumps)
            })
            // filter out kmers with no clumps
            .filter(|(_, v)| !v.is_empty())
            .collect()
    }

    #[cfg(test)]
    mod tests {
        use super::{super::tests::TEST_GENOME, *};

        const PATTERN: &[u8] = b"ATGATCAAG";
        const K: usize = 9;
        const WINDOW_SIZE: usize = 500;
        const MIN_CLUMP_SIZE: usize = 3;
        const TEST_CLUMP_SIZE: usize = 3;

        #[test]
        fn test_find_clumps_fast() {
            let clumps = find_clumps_fast(TEST_GENOME.as_ref(), K, WINDOW_SIZE, MIN_CLUMP_SIZE);
            if let Some(test_clumps) = clumps.get(PATTERN) {
                assert_eq!(test_clumps.len(), 1);
                assert_eq!(test_clumps[0].1, TEST_CLUMP_SIZE);
            } else {
                panic!("no clump found for {:?}", PATTERN);
            }
        }

        #[ignore] // really really slow
        #[test]
        fn test_find_clumps_frugal() {
            let clumps = find_clumps_frugal(TEST_GENOME.as_ref(), K, WINDOW_SIZE, MIN_CLUMP_SIZE);
            if let Some(test_clumps) = clumps.get(PATTERN) {
                assert_eq!(test_clumps.len(), 1);
                assert_eq!(test_clumps[0].1, TEST_CLUMP_SIZE);
            } else {
                panic!("no clump found for {:?}", PATTERN);
            }
        }

        #[test]
        fn test_find_clumps_optimal() {
            let clumps = find_clumps_optimal(TEST_GENOME.as_ref(), K, WINDOW_SIZE, MIN_CLUMP_SIZE);
            if let Some(test_clumps) = clumps.get(PATTERN) {
                assert_eq!(test_clumps.len(), 1);
                assert_eq!(test_clumps[0].1, TEST_CLUMP_SIZE);
            } else {
                panic!("no clump found for {:?}", PATTERN);
            }
        }
    }
}

pub mod ch1_7 {
    #[derive(Default)]
    pub struct SkewLimits {
        pub min_i: usize,
        pub min_val: isize,
        pub max_i: usize,
        pub max_val: isize,
    }

    pub fn skew_limits<'a, T: Into<&'a [u8]>>(text: T) -> SkewLimits {
        let text: &[u8] = text.into();
        let mut val = 0isize;
        let mut skew = SkewLimits::default();
        for (i, x) in text
            .iter()
            .map(|c| match c {
                b'c' | b'C' => -1isize,
                b'g' | b'G' => 1isize,
                _ => 0isize,
            })
            .enumerate()
        {
            if x != 0 {
                val += x;
                if val < skew.min_val {
                    skew.min_i = i;
                    skew.min_val = val;
                } else if val > skew.max_val {
                    skew.max_i = i;
                    skew.max_val = val;
                }
            }
        }
        skew
    }
}

pub mod ch1_8 {
    use bitvec::prelude::*;
    use std::collections::{HashMap, VecDeque};

    pub fn dist(a: &[u8], b: &[u8]) -> usize {
        assert_eq!(a.len(), b.len());
        a.iter().zip(b).filter(|(x, y)| x != y).count()
    }

    const A: u8 = 0;
    const C: u8 = 1;
    const G: u8 = 2;
    const T: u8 = 3;

    fn kmer_to_bitvec(kmer: &[u8]) -> BitVec {
        let n = kmer.len() * 2;
        let mut b = BitVec::with_capacity(n);
        for (c, i) in kmer.iter().zip((0..n).step_by(2)) {
            match c {
                b'A' | b'a' => b[i..i + 1].store::<u8>(A),
                b'C' | b'c' => b[i..i + 1].store::<u8>(C),
                b'G' | b'g' => b[i..i + 1].store::<u8>(G),
                b'T' | b't' => b[i..i + 1].store::<u8>(T),
                _ => panic!("invalid character {}", c),
            }
        }
        b
    }

    pub fn find_pattern_approx<'a, T: Into<&'a [u8]>>(
        text: T,
        pattern: &[u8],
        d: usize,
    ) -> Vec<usize> {
        let text: &[u8] = text.into();
        let mut matches = Vec::new();
        for (i, kmer) in text.windows(pattern.len()).enumerate() {
            if pattern
                .iter()
                .zip(kmer)
                .scan(0usize, |mismatches, (a, b)| {
                    if *mismatches > d {
                        None
                    } else {
                        if a != b {
                            *mismatches += 1;
                        }
                        Some(*mismatches <= d)
                    }
                })
                .last()
                .unwrap_or(true)
            {
                matches.push(i)
            }
        }
        matches
    }

    fn complement(c: u8) -> u8 {
        if c == b'N' {
            c
        } else if c & 2 != 0 {
            c ^ 4
        } else {
            c ^ 21
        }
    }

    pub fn find_clumps<'a, T: Into<&'a [u8]>>(
        text: T,
        k: usize,
        window_size: usize,
        min_clump_size: usize,
    ) -> HashMap<Vec<u8>, Vec<(usize, usize)>> {
        let text: &[u8] = text.into();
        let max_end_offset = window_size - k;
        // for each kmer, store the current clump and the completed clumps
        let mut clumps: HashMap<Vec<u8>, (VecDeque<usize>, Vec<(usize, usize)>)> = HashMap::new();
        for (i, kmer) in text.windows(k).enumerate() {
            // use the lexicographical minimum between the kmer and its reverse-complement
            let mut key = Vec::from(kmer);
            let rc = kmer
                .iter()
                .rev()
                .map(|c| complement(*c))
                .collect::<Vec<_>>();
            if rc < key {
                key = rc;
            }
            clumps
                .entry(key)
                .and_modify(|(cur_clump, found_clumps)| {
                    // add a separate clump for each set of offsets where the first and the last
                    // offset are within `window_size` of each other and the number of offsets
                    // is >= `min_clump_size`
                    while !cur_clump.is_empty() && i - cur_clump.front().unwrap() > max_end_offset {
                        let clump_size = cur_clump.len();
                        let start = cur_clump.pop_front().unwrap();
                        if clump_size >= min_clump_size {
                            found_clumps.push((start, clump_size))
                        }
                    }
                    cur_clump.push_back(i);
                })
                .or_insert((VecDeque::from([i]), Vec::new()));
        }
        clumps
            .into_iter()
            .map(|(kmer, (cur_clump, mut found_clumps))| {
                // add the last clump if it is long enough
                if cur_clump.len() >= min_clump_size {
                    found_clumps.push((*cur_clump.front().unwrap(), cur_clump.len()))
                }
                (kmer, found_clumps)
            })
            // filter out kmers with no clumps
            .filter(|(_, v)| !v.is_empty())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use lazy_static::lazy_static;
    use reqwest::blocking::get;
    use std::io::{copy, BufWriter};

    const GENOME_URL: &str =
        "http://bioinformaticsalgorithms.com/data/realdatasets/Replication/Vibrio_cholerae.txt";
    const GENOME_SIZE: usize = 1_108_250;

    lazy_static! {
        pub static ref TEST_GENOME: Vec<u8> = {
            let mut resp = get(GENOME_URL).expect("request failed");
            let mut buf = BufWriter::new(Vec::with_capacity(GENOME_SIZE));
            copy(&mut resp, &mut buf).expect("failed to copy content");
            buf.into_inner().unwrap()
        };
    }
}
