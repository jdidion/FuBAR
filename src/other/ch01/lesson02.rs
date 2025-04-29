use pcre2::bytes::Regex;
use std::collections::HashMap;

/// Implementation of the PatternCount function from Lession 1.2
pub fn pattern_count<T: AsRef<[u8]>, P: AsRef<[u8]>>(text: T, pattern: P) -> usize {
    let text = text.as_ref();
    let pattern = pattern.as_ref();
    // return early if the text is shorter than the pattern
    if text.len() < pattern.len() {
        return 0;
    }
    // count the occurrences of `pattern` in `text`
    let mut count = 0;
    for i in 0..=(text.len() - pattern.len()) {
        if text[i..i + pattern.len()] == *pattern {
            count += 1;
        }
    }
    count
}

/// Implement PatternCount using an iterator over `text` (i.e., functional style)
pub fn pattern_count_func<T: AsRef<[u8]>, P: AsRef<[u8]>>(text: T, pattern: P) -> usize {
    let text = text.as_ref();
    let pattern = pattern.as_ref();
    text.windows(pattern.len()).filter(|value| *value == pattern).count()
}

/// Implement PatternCount using regular expression search. Note that we need to use `pcre2` rather
/// than `regex` because the latter does not support look-ahead.
pub fn pattern_count_re<T: AsRef<[u8]>, P: Into<String>>(text: T, pattern: P) -> usize {
    let text = text.as_ref();
    let pattern = pattern.into();
    let re = Regex::new(format!("(?={pattern})").as_str()).unwrap();
    re.find_iter(text).count()
}

/// Implementation of BetterFrequentWordS from Lesson 1.2.
pub fn frequent_words<'a, T: Into<&'a [u8]>>(text: T, k: usize) -> (Vec<&'a [u8]>, usize) {
    let frequency_table = frequency_table(text.into(), k);
    let max_count = match frequency_table.values().max() {
        None => 0,
        Some(max_count) => *max_count,
    };
    let most_frequent = match max_count {
        0 => Vec::new(),
        max_count => {
            let mut most_frequent: Vec<_> = frequency_table
                .into_iter()
                .filter_map(|(kmer, count)| if count == max_count { Some(kmer) } else { None })
                .collect();
            most_frequent.sort();
            most_frequent
        }
    };
    (most_frequent, max_count)
}

/// Helper function to compute the frequency table for kmers of size `k` within `text`.
fn frequency_table(text: &[u8], k: usize) -> HashMap<&[u8], usize> {
    match text.len() {
        n if n < k => HashMap::new(),
        n if n == k => HashMap::from([(text, 1)]),
        _ => {
            let mut frequency_table = HashMap::new();
            for kmer in text.windows(k) {
                frequency_table.entry(kmer).and_modify(|value| *value += 1).or_insert(1);
            }
            frequency_table
        }
    }
}

#[cfg(test)]
mod tests {
    use std::sync::LazyLock;

    #[test]
    fn test_pattern_count() {
        assert_eq!(super::pattern_count(String::from("ACGT"), String::from("ACGT")), 1);
    }

    #[test]
    fn test_pattern_count_func() {
        assert_eq!(super::pattern_count_func(String::from("ACGT"), String::from("ACGT")), 1);
    }

    #[test]
    fn test_pattern_count_re() {
        assert_eq!(super::pattern_count_re(String::from("ACGT"), String::from("ACGT")), 1);
    }

    // const TEXT_ORI: &[u8] =
    //     b"atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaacctgagtggatgacatcaagatagg\
    //         tcgttgtatctccttcctctcgtactctcatgaccacggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaata\
    //         cttgtgacttgtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggattacgaaagcatgatcat\
    //         ggctgttgttctgtttatcttgttttgactgagacttgttaggatagacggtttttcatcactgactagccaaagccttactctgc\
    //         ctgacatcgaccgtaaattgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaagatcttcaa\
    //         ttgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtttccttaaccctctattttttacggaagaatgatc\
    //         aagctgctgctcttgatcatcgtttc";

    // lazy_static! {
    //     static ref EXPECTED_ORI: [(usize, (Vec<&'static [u8]>, usize)); 7] = [
    //         (3, (vec![b"tga"], 25)),
    //         (4, (vec![b"atga"], 12)),
    //         (5, (vec![b"gatca", b"tgatc"], 8)),
    //         (6, (vec![b"tgatca"], 8)),
    //         (7, (vec![b"atgatca"], 5)),
    //         (8, (vec![b"atgatcaa"], 4)),
    //         (9, (vec![b"atgatcaag", b"ctcttgatc", b"cttgatcat", b"tcttgatca"], 3))
    //     ];
    // }

    #[test]
    fn test_frequent_words() {
        const TEXT1: &[u8] = b"ACAACTATGCATACTATCGGGAACTATCCT";
        const K1: usize = 5;
        const TEXT2: &[u8] = b"CGATATATCCATAG";
        const K2: usize = 3;
        const TEXT3: &[u8] = b"ACACACA";
        const K3: usize = 2;

        static EXPECTED1: LazyLock<(Vec<&'static [u8]>, usize)> =
            LazyLock::new(|| (vec![b"ACTAT"], 3));
        static EXPECTED2: LazyLock<(Vec<&'static [u8]>, usize)> =
            LazyLock::new(|| (vec![b"ATA"], 3));
        static EXPECTED3: LazyLock<(Vec<&'static [u8]>, usize)> =
            LazyLock::new(|| (vec![b"AC", b"CA"], 3));

        assert_eq!(super::frequent_words(TEXT1, K1), *EXPECTED1);
        assert_eq!(super::frequent_words(TEXT2, K2), *EXPECTED2);
        assert_eq!(super::frequent_words(TEXT3, K3), *EXPECTED3);
    }
}
