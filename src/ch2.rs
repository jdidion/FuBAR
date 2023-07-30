pub mod ch2_2 {
    use lazy_static::lazy_static;
    use many_to_many::ManyToMany;
    use std::{
        collections::{hash_map::DefaultHasher, HashMap, HashSet},
        hash::{Hash, Hasher},
        sync::{Arc, Mutex, MutexGuard, RwLock},
    };

    /// Generates all nucleotide combinations of length `d` with replacement.
    pub struct MutationGenerator {
        d: usize,
        indices: Vec<usize>,
        nucs: Vec<u8>,
    }

    impl MutationGenerator {
        const NUCS: [u8; 4] = [b'A', b'C', b'G', b'T'];

        fn reset(&mut self) {
            self.indices.clear();
            self.nucs.clear();
        }

        /// Returns the next combination of nucleotides, or `None` if the generator is exhausted.
        pub fn next(&mut self) -> Option<&Vec<u8>> {
            if self.indices.is_empty() {
                // first iteration
                self.indices.resize(self.d, 0);
                self.nucs.resize(self.d, b'A');
            } else {
                match self
                    .indices
                    .iter()
                    .enumerate()
                    .rev()
                    .skip_while(|(_, idx)| **idx == 3)
                    .next()
                {
                    Some((i, idx)) => {
                        let new_idx = idx + 1;
                        self.indices[i] = new_idx;
                        self.nucs[i] = Self::NUCS[new_idx];
                        for j in i + 1..self.d {
                            self.indices[j] = 0;
                            self.nucs[j] = b'A';
                        }
                    }
                    None => return None,
                }
            }
            Some(&self.nucs)
        }
    }

    /// Wrapper for a MutationGenerator that enables a single instance to be cached and reused.
    pub struct MutationGeneratorSingleton(Mutex<MutationGenerator>);

    lazy_static! {
        /// Cache of MutationGeneratorSingleton instances by value of `d`.
        static ref MUTATION_GENERATORS: RwLock<HashMap<usize, Arc<MutationGeneratorSingleton>>> =
            RwLock::new(HashMap::new());
    }

    impl MutationGeneratorSingleton {
        fn new(d: usize) -> Self {
            MutationGeneratorSingleton(Mutex::new(MutationGenerator {
                d,
                indices: Vec::with_capacity(d),
                nucs: Vec::with_capacity(d),
            }))
        }

        /// Returns the `MutationGeneratorSingleton` for the given number of mutations `d`.
        pub fn get(d: usize) -> Arc<MutationGeneratorSingleton> {
            MUTATION_GENERATORS
                .write()
                .unwrap()
                .entry(d)
                .or_insert_with(|| Arc::new(MutationGeneratorSingleton::new(d)))
                .clone()
        }

        /// Locks the contained `MutationGenerator`, calls `reset()`, and returns the guard.
        pub fn lock_and_reset(&self) -> MutexGuard<MutationGenerator> {
            let mut guard = self.0.lock().unwrap();
            guard.reset();
            guard
        }
    }

    /// Generates all possible combinations of size 1..=d from a range 0..k with no replacement.
    pub struct MotifGenerator {
        k: usize,
        d: usize,
        dmax: usize,
        indices: Vec<usize>,
    }

    impl MotifGenerator {
        pub fn new(k: usize, dmax: usize) -> Self {
            assert!(dmax > 0);
            assert!(k >= dmax);
            MotifGenerator {
                k,
                d: 0,
                dmax,
                indices: Vec::with_capacity(dmax),
            }
        }

        fn reset(&mut self) {
            self.d = 0;
            self.indices.clear();
        }

        /// Returns the next motif or `None` if there are no more combinations.
        fn next(&mut self) -> Option<&[usize]> {
            if self.d == 0 {
                self.d = 1;
                self.indices.push(0);
            } else {
                // scan from the end looking for an index to increment
                match (0..self.d)
                    .zip(self.k - self.d..self.k)
                    .rev()
                    .skip_while(|(i, k)| self.indices[*i] == *k)
                    .next()
                {
                    Some((i, _)) => {
                        // increment the index at `i` and reset the ones to the right
                        self.indices[i] += 1;
                        for j in i + 1..self.d {
                            self.indices[j] = self.indices[j - 1] + 1;
                        }
                    }
                    None if self.d == self.dmax => {
                        // reached the last combination of kmax
                        return None;
                    }
                    None => {
                        // reached the last combination of the current k but there are still more
                        // k <= kmax
                        self.d += 1;
                        self.indices.clear();
                        self.indices.extend(0..self.d);
                    }
                }
            }

            Some(&self.indices[..self.d])
        }

        /// Adds all motifs - i.e. all kmers with edit distance from `kmer` <= `kmax` - to `accu`,
        /// which is a `ManyToMany` data structure that preserves the relationships between the
        /// source kmers and their motifs. The source kmer (on the "left" side of the map) will
        /// always be among the motifs (on the "right" side of the map). Note that we store the
        /// hash value for each motif, rather than the motif itself, to take advantage of the fact
        /// that the hash of a `&[u8]` and a `Vec<u8>` containing the same bytes is identical.
        /// This simplifies the lookup (otherwise we'd have to allocate a new `Vec<u8>` from a kmer
        /// (`&[u8]`) to look it up in the right side). There is a small chance of collision, but
        /// we could look at using an alternative Hasher specialized to sequences drawn from small
        /// alphabets (e.g. DNA).
        pub fn add_all_motifs<'kmer>(
            &mut self,
            kmer: &'kmer [u8],
            accu: &mut ManyToMany<&'kmer [u8], u64>,
        ) {
            self.reset();
            while let Some(indices) = self.next() {
                let mut motif = Vec::from(kmer);
                let mutex = MutationGeneratorSingleton::get(indices.len());
                let mut mutations = mutex.lock_and_reset();
                while let Some(nucs) = mutations.next() {
                    for (i, nuc) in indices.iter().zip(nucs) {
                        motif[*i] = *nuc;
                    }
                    let mut hasher = DefaultHasher::new();
                    motif.hash(&mut hasher);
                    accu.insert(kmer, hasher.finish());
                }
            }
        }
    }

    /// Generates all motifs of size `k` with edit distance <= `d` from the first sequence in
    /// `dna`, then scans all subsequent sequences to find motifs that appear in all of them.
    /// Returns a set of motifs, each represented by the first kmer found in the first sequence.
    pub fn find_all_motifs<'a, Dna: AsRef<[&'a [u8]]>>(
        dna: Dna,
        k: usize,
        d: usize,
    ) -> HashSet<&'a [u8]> {
        let dna = dna.as_ref();
        let mut c = MotifGenerator::new(k, d);
        let mut motifs: ManyToMany<&[u8], u64> = ManyToMany::new();
        // iterate first dna in windows of size k
        for kmer in dna[0].windows(k) {
            // for each kmer, iterate through all possible mutations with edit distance <= d
            // and add them to the set of all kmers
            c.add_all_motifs(kmer, &mut motifs);
        }
        // reduce the set of motifs by intersecting with kmers in each other dna sequence
        let final_motifs = dna[1..].iter().fold(motifs, |accu, seq| {
            let mut new_accu: ManyToMany<&[u8], u64> = ManyToMany::new();
            for motif in seq.windows(k) {
                let mut hasher = DefaultHasher::new();
                motif.hash(&mut hasher);
                if let Some(kmers) = accu.get_right(&hasher.finish()) {
                    for kmer in kmers {
                        for kmer_motif in accu.get_left(&kmer).unwrap() {
                            new_accu.insert(kmer, kmer_motif);
                        }
                    }
                }
            }
            new_accu
        });
        // convert final motifs to owned values so they can be returned
        final_motifs
            .get_left_keys()
            .into_iter()
            .map(|k| *k)
            .collect()
    }

    #[cfg(test)]
    mod tests {
        use super::{find_all_motifs, MotifGenerator, MutationGeneratorSingleton};

        #[test]
        fn test_mutation_generator() {
            let mutex = MutationGeneratorSingleton::get(2);
            let mut m = mutex.lock_and_reset();
            let mut actual: Vec<Vec<u8>> = Vec::with_capacity(16);
            while let Some(i) = m.next() {
                actual.push(i.to_owned())
            }
            let expected: Vec<&[u8]> = vec![
                b"AA", b"AC", b"AG", b"AT", b"CA", b"CC", b"CG", b"CT", b"GA", b"GC", b"GG", b"GT",
                b"TA", b"TC", b"TG", b"TT",
            ];
            assert_eq!(actual.len(), expected.len());
            for (a, e) in actual.iter().zip(expected) {
                assert_eq!(a.as_slice(), e);
            }
        }

        #[test]
        fn test_combinations() {
            let mut c = MotifGenerator::new(5, 3);
            // (5,3) + (5,2) + (5,1) == 25
            let mut actual: Vec<Vec<usize>> = Vec::with_capacity(25);
            while let Some(i) = c.next() {
                actual.push(i.to_owned());
            }
            assert_eq!(actual.len(), 25);
            assert_eq!(
                actual,
                vec![
                    vec![0],
                    vec![1],
                    vec![2],
                    vec![3],
                    vec![4],
                    vec![0, 1],
                    vec![0, 2],
                    vec![0, 3],
                    vec![0, 4],
                    vec![1, 2],
                    vec![1, 3],
                    vec![1, 4],
                    vec![2, 3],
                    vec![2, 4],
                    vec![3, 4],
                    vec![0, 1, 2],
                    vec![0, 1, 3],
                    vec![0, 1, 4],
                    vec![0, 2, 3],
                    vec![0, 2, 4],
                    vec![0, 3, 4],
                    vec![1, 2, 3],
                    vec![1, 2, 4],
                    vec![1, 3, 4],
                    vec![2, 3, 4],
                ]
            )
        }

        #[ignore] // ignoring expensive test; to run do `cargo test -- --ignore`
        #[test]
        fn test_motif_enumeration_big() {
            let dna: Vec<&[u8]> = vec![
                b"ATGACCGGGATACTGATAAAAAAAAGGGGGGGGGCGTACACATTAGATAAACGTATGAAGTACGTTAGACTCGGCGCCGCCG",
                b"ACCCCTATTTTTTGAGCAGATTTAGTGACCTGGAAAAAAAATTTGAGTACAAAACTTTTCCGAATACAATAAAACGGCGGGA",
                b"TGAGTATCCCTGGGATGACTTAAAATAATGGAGTGGTGCTCTCCCGATTTTTGAATATGTAGGATCATTCGCCAGGGTCCGA",
                b"GCTGAGAATTGGATGCAAAAAAAGGGATTGTCCACGCAATCGCGAACCAACGCGGACCCAAAGGCAAGACCGATAAAGGAGA",
                b"TCCCTTTTGCGGTAATGTGCCGGGAGGCTGGTTACGTAGGGAAGCCCTAACGGACTTAATATAATAAAGGAAGGGCTTATAG",
                b"GTCAATCATGTTCTTGTGAATGGATTTAACAATAAGGGCTGGGACCGCTTGGCGCACCCAAATTCAGTGTGGGCGAGCGCAA",
                b"CGGTTTTGGCCCTTGTTAGAGGCCCCCGTATAAACAAGGAGGGCCAATTATGAGAGAGCTAATCTATCGCGTGCGTGTTCAT",
                b"AACTTGAGTTAAAAAATAGGGAGCCCTGGGGCACATACAAGAGGAGTCTTCCTTATCAGTTAATGCTGTATGACACTATGTA",
                b"TTGGCCCATTGGCTAAAAGCCCAACTTGACAAATGGAAGATAGAATCCTTGCATACTAAAAAGGAGCGGACCGAAAGGGAAG",
                b"CTGGTGAGCAACGACAGATTCTTACGTGCATTAGCTCGCTTCCGGGGATCTAATAGCACGAAGCTTACTAAAAAGGAGCGGA"
            ];
            let motifs = find_all_motifs(dna, 15, 4);
            assert_eq!(motifs.len(), 1);
            assert_eq!(motifs.iter().next().unwrap(), &b"AAAAAAAAGGGGGGG");
        }
    }
}

pub mod ch2_4 {
    use super::ch2_2::MutationGeneratorSingleton;

    pub fn find_median_string<'a, Dna: AsRef<[&'a [u8]]>>(dna: Dna, k: usize) -> (Vec<u8>, usize) {
        let dna = dna.as_ref();
        let mutex = MutationGeneratorSingleton::get(k);
        let mut generator = mutex.lock_and_reset();
        let mut median_string: Option<Vec<u8>> = None;
        let mut distance: usize = usize::MAX;
        while let Some(pattern) = generator.next() {
            let pattern_distance = dna
                .iter()
                .map(|seq| {
                    let mut seq_distance: usize = usize::MAX;
                    for kmer in seq.windows(k) {
                        let kmer_distance =
                            pattern.iter().zip(kmer).filter(|(a, b)| a != b).count();
                        if kmer_distance < seq_distance {
                            seq_distance = kmer_distance;
                            if kmer_distance == 0 {
                                break;
                            }
                        }
                    }
                    seq_distance
                })
                .sum();
            if pattern_distance < distance {
                distance = pattern_distance;
                median_string.replace(pattern.clone());
            }
        }
        (median_string.unwrap(), distance)
    }

    #[cfg(test)]
    mod tests {
        use super::find_median_string;

        #[test]
        fn test_find_median_string() {
            let dna: Vec<&[u8]> = vec![
                b"ACCAGTAGGAGCA",
                b"TGGCGACCGATTG",
                b"AGGATCAGGATAA",
                b"CCCGGATTACTGA",
                b"TAAGTAAATCGAT",
            ];
            let (median_string, distance) = find_median_string(dna, 3);
            assert_eq!(median_string, b"GAT");
            assert_eq!(distance, 1);
        }
    }
}

pub mod ch2_5 {
    use rand::distributions::WeightedIndex;
    use rand::prelude::*;

    pub struct Profile<'a> {
        k: usize,
        n: usize,
        pub motifs: Vec<&'a [u8]>,
        counts: Vec<Vec<usize>>,
    }

    impl<'a> Profile<'a> {
        const NUCS: [u8; 4] = [b'A', b'C', b'G', b'T'];

        fn nuc_to_row(nuc: &u8) -> usize {
            match nuc {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => panic!("invalid character {}", nuc),
            }
        }

        pub fn add(&mut self, motif: &'a [u8]) {
            assert_eq!(motif.len(), self.k, "all motifs must have the same length");
            for (i, nuc) in motif.iter().enumerate() {
                // TODO: test if a lookup table is faster
                let row: usize = Self::nuc_to_row(nuc);
                self.counts[row][i] += 1;
            }
            self.motifs.push(motif);
        }

        pub fn new<Dna: AsRef<[&'a [u8]]>>(motifs: Dna) -> Self {
            let motifs = motifs.as_ref();
            assert!(!motifs.is_empty());
            let n = motifs.len();
            let k = motifs[0].len();
            let mut profile = Self {
                k,
                n,
                motifs: Vec::with_capacity(n),
                counts: vec![vec![0; k]; 4],
            };
            for motif in motifs.iter() {
                profile.add(motif);
            }
            profile
        }

        pub fn empty(k: usize) -> Self {
            Self {
                k,
                n: 0,
                motifs: Vec::new(),
                counts: vec![vec![0; k]; 4],
            }
        }

        pub fn with_one(kmer: &'a [u8]) -> Self {
            let k = kmer.len();
            let mut profile = Self {
                k,
                n: 1,
                motifs: Vec::with_capacity(1),
                counts: vec![vec![0; k]; 4],
            };
            profile.add(kmer);
            profile
        }

        pub fn consensus(&self) -> Vec<u8> {
            (0..self.k)
                .map(|i| {
                    let (_, nuc) = self
                        .counts
                        .iter()
                        .enumerate()
                        .map(|(row, counts)| (counts[i], row))
                        .max()
                        .unwrap();
                    Self::NUCS[nuc]
                })
                .collect()
        }

        pub fn score(&self) -> usize {
            let matches: usize = (0..self.k)
                .map(|i| self.counts.iter().map(|v| v[i]).max().unwrap())
                .sum();
            (self.k * self.n) - matches
        }

        pub fn probability_of(&self, kmer: &[u8]) -> f64 {
            kmer.iter()
                .enumerate()
                .map(|(i, nuc)| self.counts[Self::nuc_to_row(nuc)][i] as f64 / self.n as f64)
                .product()
        }

        pub fn gen_random_kmer(&self) -> Vec<u8> {
            let mut rng = thread_rng();
            (0..self.k)
                .map(|i| {
                    let weights: Vec<usize> = self.counts.iter().map(|v| v[i]).collect();
                    let dist = WeightedIndex::new(&weights).unwrap();
                    Self::NUCS[dist.sample(&mut rng)]
                })
                .collect()
        }

        pub fn find_most_probable(&self, dna: &'a [u8]) -> &'a [u8] {
            let dna: &[u8] = dna.as_ref();
            assert!(dna.len() >= self.k);
            let mut most_probable_kmer: Option<&[u8]> = None;
            let mut max_p: f64 = -1.0;
            for kmer in dna.windows(self.k) {
                let p = self.probability_of(kmer);
                if p > max_p {
                    most_probable_kmer.replace(kmer);
                    max_p = p;
                }
            }
            assert!(max_p >= 0.0);
            most_probable_kmer.unwrap()
        }
    }

    pub fn find_best_profile_greedy<'a, Dna: AsRef<[&'a [u8]]>>(dna: Dna, k: usize) -> Profile<'a> {
        let dna = dna.as_ref();
        assert!(!dna.is_empty());
        let mut best_profile = Profile::new(dna.iter().map(|d| &d[..k]).collect::<Vec<_>>());
        let mut best_score: usize = best_profile.score();
        for kmer in dna[0].windows(k) {
            let mut profile = Profile::with_one(kmer);
            for seq in &dna[1..] {
                let motif = profile.find_most_probable(seq);
                profile.add(motif);
            }
            let score = profile.score();
            if score < best_score {
                best_profile = profile;
                best_score = score;
            }
        }
        best_profile
    }

    #[cfg(test)]
    mod tests {
        use super::{find_best_profile_greedy, Profile};

        #[test]
        fn test_profile_most_probable_kmer() {
            let dna: &[u8] = b"GGTACGGGGATTACCT";
            let motifs: Vec<&[u8]> = vec![
                b"TCGGGGGTTTTT",
                b"CCGGTGACTTAC",
                b"ACGGGGATTTTC",
                b"TTGGGGACTTTT",
                b"AAGGGGACTTCC",
                b"TTGGGGACTTCC",
                b"TCGGGGATTCAT",
                b"TCGGGGATTCCT",
                b"TAGGGGAACTAC",
                b"TCGGGTATAACC",
            ];
            let profile = Profile::new(&motifs);
            assert_eq!(profile.k, 12);
            assert_eq!(profile.score(), 30);
            assert_eq!(profile.consensus(), b"TCGGGGATTTCC");
            assert_eq!(profile.find_most_probable(dna), b"ACGGGGATTACC");
        }

        #[test]
        fn test_greedy_motif_search() {
            let dna: Vec<&[u8]> = vec![
                b"ATGACCGGGATACTGATAAAAAAAAGGGGGGGGGCGTACACATTAGATAAACGTATGAAGTACGTTAGACTCGGCGCCGCCG",
                b"ACCCCTATTTTTTGAGCAGATTTAGTGACCTGGAAAAAAAATTTGAGTACAAAACTTTTCCGAATACAATAAAACGGCGGGA",
                b"TGAGTATCCCTGGGATGACTTAAAATAATGGAGTGGTGCTCTCCCGATTTTTGAATATGTAGGATCATTCGCCAGGGTCCGA",
                b"GCTGAGAATTGGATGCAAAAAAAGGGATTGTCCACGCAATCGCGAACCAACGCGGACCCAAAGGCAAGACCGATAAAGGAGA",
                b"TCCCTTTTGCGGTAATGTGCCGGGAGGCTGGTTACGTAGGGAAGCCCTAACGGACTTAATATAATAAAGGAAGGGCTTATAG",
                b"GTCAATCATGTTCTTGTGAATGGATTTAACAATAAGGGCTGGGACCGCTTGGCGCACCCAAATTCAGTGTGGGCGAGCGCAA",
                b"CGGTTTTGGCCCTTGTTAGAGGCCCCCGTATAAACAAGGAGGGCCAATTATGAGAGAGCTAATCTATCGCGTGCGTGTTCAT",
                b"AACTTGAGTTAAAAAATAGGGAGCCCTGGGGCACATACAAGAGGAGTCTTCCTTATCAGTTAATGCTGTATGACACTATGTA",
                b"TTGGCCCATTGGCTAAAAGCCCAACTTGACAAATGGAAGATAGAATCCTTGCATACTAAAAAGGAGCGGACCGAAAGGGAAG",
                b"CTGGTGAGCAACGACAGATTCTTACGTGCATTAGCTCGCTTCCGGGGATCTAATAGCACGAAGCTTACTAAAAAGGAGCGGA"
            ];
            let best_profile = find_best_profile_greedy(dna, 15);
            assert_eq!(best_profile.score(), 58);
        }
    }
}
