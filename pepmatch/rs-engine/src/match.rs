use memmap2::Mmap;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;

const PROTEIN_INDEX_MULTIPLIER: u64 = 100_000_000;
const HEADER_SIZE: usize = 8 + 1 + 4 + 8 + 8 + 8 + 8;

struct PepIndex {
    mmap: Mmap,
    k: usize,
    num_proteins: usize,
    seq_offset: usize,
    seq_len: usize,
    protein_offsets_offset: usize,
    metadata_offsets_offset: usize,
    metadata_offset: usize,
    kmer_offsets_offset: usize,
    kmer_data_offset: usize,
    num_kmers: usize,
}

unsafe impl Sync for PepIndex {}

impl PepIndex {
    fn open(path: &str) -> Self {
        let file = File::open(path).expect("Cannot open .pepidx file");
        let mmap = unsafe { Mmap::map(&file).expect("Cannot mmap file") };

        let mut pos = 8;
        let k = mmap[pos] as usize;
        pos += 1;
        let num_proteins = u32::from_le_bytes(mmap[pos..pos + 4].try_into().unwrap()) as usize;
        pos += 4;
        let seq_len = u64::from_le_bytes(mmap[pos..pos + 8].try_into().unwrap()) as usize;
        pos += 8;
        let num_kmers = u64::from_le_bytes(mmap[pos..pos + 8].try_into().unwrap()) as usize;
        pos += 8;
        let _total_positions = u64::from_le_bytes(mmap[pos..pos + 8].try_into().unwrap());
        pos += 8;
        let metadata_buf_len = u64::from_le_bytes(mmap[pos..pos + 8].try_into().unwrap()) as usize;

        let seq_offset = HEADER_SIZE;
        let protein_offsets_offset = seq_offset + seq_len;
        let metadata_offsets_offset = protein_offsets_offset + num_proteins * 8;
        let metadata_offset = metadata_offsets_offset + num_proteins * 8;
        let kmer_offsets_offset = metadata_offset + metadata_buf_len;
        let kmer_data_offset = kmer_offsets_offset + num_kmers * 8;

        PepIndex {
            mmap, k, num_proteins, seq_offset, seq_len,
            protein_offsets_offset, metadata_offsets_offset,
            metadata_offset, kmer_offsets_offset, kmer_data_offset, num_kmers,
        }
    }

    fn protein_offset(&self, protein_number: usize) -> u64 {
        let pos = self.protein_offsets_offset + (protein_number - 1) * 8;
        u64::from_le_bytes(self.mmap[pos..pos + 8].try_into().unwrap())
    }

    fn resolve(&self, encoded_idx: u64) -> Option<&[u8]> {
        let protein_number = (encoded_idx / PROTEIN_INDEX_MULTIPLIER) as usize;
        let position = (encoded_idx % PROTEIN_INDEX_MULTIPLIER) as usize;
        if protein_number == 0 || protein_number > self.num_proteins {
            return None;
        }
        let base = self.protein_offset(protein_number) as usize;
        let start = self.seq_offset + base + position;
        let end = start + self.k;
        if end > self.seq_offset + self.seq_len {
            return None;
        }
        Some(&self.mmap[start..end])
    }

    fn kmer_offset(&self, kmer_index: usize) -> usize {
        let pos = self.kmer_offsets_offset + kmer_index * 8;
        u64::from_le_bytes(self.mmap[pos..pos + 8].try_into().unwrap()) as usize
    }

    fn kmer_at(&self, kmer_index: usize) -> &[u8] {
        let offset = self.kmer_data_offset + self.kmer_offset(kmer_index);
        &self.mmap[offset..offset + self.k]
    }

    fn positions_at(&self, kmer_index: usize) -> &[u8] {
        let offset = self.kmer_data_offset + self.kmer_offset(kmer_index) + self.k;
        let count = u32::from_le_bytes(self.mmap[offset..offset + 4].try_into().unwrap()) as usize;
        let data_start = offset + 4;
        &self.mmap[data_start..data_start + count * 8]
    }

    fn lookup(&self, kmer: &[u8]) -> Option<Vec<u64>> {
        let mut lo = 0usize;
        let mut hi = self.num_kmers;

        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let mid_kmer = self.kmer_at(mid);
            match mid_kmer.cmp(kmer) {
                std::cmp::Ordering::Equal => {
                    let pos_bytes = self.positions_at(mid);
                    let count = pos_bytes.len() / 8;
                    let mut positions = Vec::with_capacity(count);
                    for i in 0..count {
                        let idx = u64::from_le_bytes(
                            pos_bytes[i * 8..(i + 1) * 8].try_into().unwrap(),
                        );
                        positions.push(idx);
                    }
                    return Some(positions);
                }
                std::cmp::Ordering::Less => lo = mid + 1,
                std::cmp::Ordering::Greater => hi = mid,
            }
        }
        None
    }

    fn read_metadata_str(&self, pos: usize) -> (&str, usize) {
        let len = u16::from_le_bytes(self.mmap[pos..pos + 2].try_into().unwrap()) as usize;
        let s = std::str::from_utf8(&self.mmap[pos + 2..pos + 2 + len]).unwrap_or("");
        (s, pos + 2 + len)
    }

    fn get_metadata(&self, protein_number: usize) -> [String; 9] {
        let offset_pos = self.metadata_offsets_offset + (protein_number - 1) * 8;
        let meta_rel = u64::from_le_bytes(
            self.mmap[offset_pos..offset_pos + 8].try_into().unwrap(),
        ) as usize;
        let mut pos = self.metadata_offset + meta_rel;

        let mut fields: [String; 9] = Default::default();
        for field in fields.iter_mut() {
            let (s, next) = self.read_metadata_str(pos);
            *field = s.to_string();
            pos = next;
        }
        fields
    }
}

fn hamming(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).filter(|(x, y)| x != y).count()
}

fn exact_match(peptide: &str, k: usize, index: &PepIndex) -> Vec<u64> {
    let pep_bytes = peptide.as_bytes();
    if pep_bytes.len() < k {
        return vec![];
    }

    let num_kmers = pep_bytes.len() - k + 1;
    let mut target_indices: Vec<usize> = (0..num_kmers).step_by(k).collect();
    let last = num_kmers - 1;
    if !target_indices.contains(&last) {
        target_indices.push(last);
    }
    target_indices.sort();
    target_indices.dedup();

    let mut hit_counts: HashMap<u64, usize> = HashMap::new();

    for &idx in &target_indices {
        let kmer = &pep_bytes[idx..idx + k];
        if let Some(positions) = index.lookup(kmer) {
            for db_index in positions {
                let candidate = db_index - idx as u64;
                *hit_counts.entry(candidate).or_insert(0) += 1;
            }
        }
    }

    let required = target_indices.len();
    hit_counts
        .into_iter()
        .filter(|&(_, count)| count == required)
        .map(|(hit, _)| hit)
        .collect()
}

fn check_left_neighbors(
    pep_bytes: &[u8], idx: usize, kmer_hit: u64, index: &PepIndex,
    k: usize, max_mismatches: usize, mut mismatches: usize,
) -> usize {
    let mut i = 0;
    while i < idx {
        if mismatches > max_mismatches { return 100; }
        let check_idx = (kmer_hit as i64 + i as i64 - idx as i64) as u64;
        match index.resolve(check_idx) {
            Some(kmer) => mismatches += hamming(kmer, &pep_bytes[i..i + k]),
            None => return 100,
        }
        i += k;
    }
    mismatches
}

fn check_right_neighbors(
    pep_bytes: &[u8], idx: usize, kmer_hit: u64, index: &PepIndex,
    k: usize, max_mismatches: usize, mut mismatches: usize,
) -> usize {
    let mut i = k + idx;
    while i + k <= pep_bytes.len() {
        if mismatches > max_mismatches { return 100; }
        let check_idx = (kmer_hit as i64 + i as i64 - idx as i64) as u64;
        match index.resolve(check_idx) {
            Some(kmer) => mismatches += hamming(kmer, &pep_bytes[i..i + k]),
            None => return 100,
        }
        i += k;
    }
    mismatches
}

fn check_left_residues(
    pep_bytes: &[u8], idx: usize, kmer_hit: u64, index: &PepIndex,
    max_mismatches: usize, mut mismatches: usize,
) -> usize {
    for i in 0..idx {
        if mismatches > max_mismatches { return 100; }
        let check_idx = (kmer_hit as i64 + i as i64 - idx as i64) as u64;
        match index.resolve(check_idx) {
            Some(kmer) => { if kmer[0] != pep_bytes[i] { mismatches += 1; } }
            None => return 100,
        }
    }
    mismatches
}

fn check_right_residues(
    pep_bytes: &[u8], idx: usize, kmer_hit: u64, index: &PepIndex,
    k: usize, max_mismatches: usize, mut mismatches: usize,
) -> usize {
    let num_kmers = pep_bytes.len() - k + 1;
    for i in (idx + 1)..num_kmers {
        if mismatches > max_mismatches { return 100; }
        let check_idx = (kmer_hit as i64 + i as i64 - idx as i64) as u64;
        match index.resolve(check_idx) {
            Some(kmer) => { if kmer[k - 1] != pep_bytes[i + k - 1] { mismatches += 1; } }
            None => return 100,
        }
    }
    mismatches
}

fn mismatch_match(peptide: &str, k: usize, max_mismatches: usize, index: &PepIndex) -> Vec<(String, usize, u64)> {
    let pep_bytes = peptide.as_bytes();
    if pep_bytes.len() < k { return vec![]; }

    let num_kmers = pep_bytes.len() - k + 1;
    let peptide_len = pep_bytes.len();
    let mut matches: Vec<(String, usize, u64)> = Vec::new();
    let mut seen: HashSet<u64> = HashSet::new();

    if peptide_len % k == 0 {
        let mut idx = 0;
        while idx < num_kmers {
            let kmer = &pep_bytes[idx..idx + k];
            if let Some(positions) = index.lookup(kmer) {
                for kmer_hit in positions {
                    let start = (kmer_hit as i64 - idx as i64) as u64;
                    if seen.contains(&start) { continue; }
                    let mismatches = check_left_neighbors(pep_bytes, idx, kmer_hit, index, k, max_mismatches, 0);
                    let mismatches = check_right_neighbors(pep_bytes, idx, kmer_hit, index, k, max_mismatches, mismatches);
                    if mismatches <= max_mismatches {
                        let mut matched = Vec::with_capacity(peptide_len);
                        let mut i = 0;
                        let mut valid = true;
                        while i < peptide_len {
                            match index.resolve(start + i as u64) {
                                Some(km) => matched.extend_from_slice(km),
                                None => { valid = false; break; }
                            }
                            i += k;
                        }
                        if valid {
                            seen.insert(start);
                            matches.push((String::from_utf8(matched).unwrap_or_default(), mismatches, start));
                        }
                    }
                }
            }
            idx += k;
        }
    } else {
        for idx in 0..num_kmers {
            let kmer = &pep_bytes[idx..idx + k];
            if let Some(positions) = index.lookup(kmer) {
                for kmer_hit in positions {
                    let start = (kmer_hit as i64 - idx as i64) as u64;
                    if seen.contains(&start) { continue; }
                    let mismatches = check_left_residues(pep_bytes, idx, kmer_hit, index, max_mismatches, 0);
                    let mismatches = check_right_residues(pep_bytes, idx, kmer_hit, index, k, max_mismatches, mismatches);
                    if mismatches <= max_mismatches {
                        let mut matched = Vec::with_capacity(peptide_len);
                        let mut i = 0;
                        let mut valid = true;
                        while i < peptide_len {
                            if i + k > peptide_len {
                                let remaining = peptide_len - i;
                                let back = k - remaining;
                                match index.resolve(start + i as u64 - back as u64) {
                                    Some(km) => matched.extend_from_slice(&km[back..]),
                                    None => { valid = false; break; }
                                }
                            } else {
                                match index.resolve(start + i as u64) {
                                    Some(km) => matched.extend_from_slice(km),
                                    None => { valid = false; break; }
                                }
                            }
                            i += k;
                        }
                        if valid {
                            seen.insert(start);
                            matches.push((String::from_utf8(matched).unwrap_or_default(), mismatches, start));
                        }
                    }
                }
            }
        }
    }

    matches
}

fn search_peptide(query_id: &str, peptide: &str, k: usize, max_mismatches: usize, index: &PepIndex) -> Vec<Vec<String>> {
    if max_mismatches == 0 {
        let hits = exact_match(peptide, k, index);
        if hits.is_empty() {
            return vec![make_empty_row(query_id, peptide)];
        }
        hits.iter().map(|&hit| make_hit_row(query_id, peptide, peptide, 0, hit, index)).collect()
    } else {
        let hits = mismatch_match(peptide, k, max_mismatches, index);
        if hits.is_empty() {
            return vec![make_empty_row(query_id, peptide)];
        }
        hits.iter().map(|(matched_seq, mm, start)| make_hit_row(query_id, peptide, matched_seq, *mm, *start, index)).collect()
    }
}

fn make_hit_row(query_id: &str, peptide: &str, matched_seq: &str, mismatches: usize, encoded_start: u64, index: &PepIndex) -> Vec<String> {
    let protein_number = (encoded_start / PROTEIN_INDEX_MULTIPLIER) as usize;
    let position = (encoded_start % PROTEIN_INDEX_MULTIPLIER) as usize;
    let meta = index.get_metadata(protein_number);

    let mutated: Vec<String> = peptide.as_bytes().iter()
        .zip(matched_seq.as_bytes())
        .enumerate()
        .filter(|(_, (a, b))| a != b)
        .map(|(i, _)| (i + 1).to_string())
        .collect();

    vec![
        query_id.to_string(),
        peptide.to_string(),
        matched_seq.to_string(),
        meta[0].clone(),
        meta[1].clone(),
        meta[2].clone(),
        meta[3].clone(),
        meta[4].clone(),
        mismatches.to_string(),
        format!("[{}]", mutated.join(", ")),
        (position + 1).to_string(),
        (position + peptide.len()).to_string(),
        meta[5].clone(),
        meta[6].clone(),
        meta[7].clone(),
        meta[8].clone(),
    ]
}

fn make_empty_row(query_id: &str, peptide: &str) -> Vec<String> {
    vec![
        query_id.to_string(),
        peptide.to_string(),
        String::new(), String::new(), String::new(), String::new(),
        String::new(), String::new(), String::new(), String::from("[]"),
        String::new(), String::new(), String::new(), String::new(),
        String::new(), String::new(),
    ]
}

pub(crate) fn run_discontinuous(
    pepidx_path: &str,
    epitopes: Vec<(String, Vec<(char, usize)>)>,
    max_mismatches: usize,
) -> Vec<Vec<String>> {
    let index = PepIndex::open(pepidx_path);
    let seq = &index.mmap[index.seq_offset..index.seq_offset + index.seq_len];

    let mut all_results: Vec<Vec<String>> = Vec::new();

    for (query_id, residues) in &epitopes {
        let mut found = false;
        let query_str: String = residues.iter()
            .map(|(r, p)| format!("{}{}", r, p))
            .collect::<Vec<_>>()
            .join(", ");

        for pn in 1..=index.num_proteins {
            let base = index.protein_offset(pn) as usize;
            let end = if pn < index.num_proteins {
                index.protein_offset(pn + 1) as usize
            } else {
                index.seq_len
            };
            let protein_len = end - base;

            let mut mismatches = 0;
            let mut valid = true;

            for &(residue, position) in residues {
                if position == 0 || position > protein_len {
                    valid = false;
                    break;
                }
                if seq[base + position - 1] as char != residue {
                    mismatches += 1;
                    if mismatches > max_mismatches {
                        valid = false;
                        break;
                    }
                }
            }

            if valid {
                found = true;
                let meta = index.get_metadata(pn);
                let matched_str: String = residues.iter()
                    .map(|&(_, p)| format!("{}{}", seq[base + p - 1] as char, p))
                    .collect::<Vec<_>>()
                    .join(", ");
                let mutated: Vec<String> = residues.iter()
                    .filter(|&&(r, p)| seq[base + p - 1] as char != r)
                    .map(|&(_, p)| p.to_string())
                    .collect();

                all_results.push(vec![
                    query_id.clone(),
                    query_str.clone(),
                    matched_str,
                    meta[0].clone(), meta[1].clone(), meta[2].clone(),
                    meta[3].clone(), meta[4].clone(),
                    mismatches.to_string(),
                    format!("[{}]", mutated.join(", ")),
                    residues.first().map(|&(_, p)| p.to_string()).unwrap_or_default(),
                    residues.last().map(|&(_, p)| p.to_string()).unwrap_or_default(),
                    meta[5].clone(), meta[6].clone(), meta[7].clone(), meta[8].clone(),
                ]);
            }
        }

        if !found {
            let mut row = vec![query_id.clone(), query_str.clone()];
            row.extend(std::iter::repeat(String::new()).take(6));
            row.push(String::new());
            row.push(String::from("[]"));
            row.extend(std::iter::repeat(String::new()).take(6));
            all_results.push(row);
        }
    }

    all_results
}

fn minimal_coverage_seeds(query: &[u8], k: usize) -> Vec<(&[u8], usize)> {
    let query_len = query.len();
    let mut seeds: Vec<(&[u8], usize)> = Vec::new();
    let mut j = 0;
    while j + k <= query_len {
        seeds.push((&query[j..j + k], j));
        j += k;
    }
    if query_len >= k {
        let last_start = query_len - k;
        if seeds.last().map(|(_, s)| *s) != Some(last_start) {
            seeds.push((&query[last_start..], last_start));
        }
    }
    seeds
}

fn is_terminal_deletion(q_idx: isize, query_len: usize, p_idx: isize, protein_len: usize) -> bool {
    if q_idx == 0 || q_idx == query_len as isize - 1 {
        return true;
    }
    if p_idx <= 0 || p_idx >= protein_len as isize - 1 {
        return true;
    }
    false
}

fn dfs(
    query: &[u8],
    q_idx: isize,
    protein: &[u8],
    p_idx: isize,
    indels_left: usize,
    direction: isize,
) -> Vec<usize> {
    if (direction == 1 && q_idx >= query.len() as isize)
        || (direction == -1 && q_idx < 0)
    {
        return vec![0];
    }

    let mut all_paths: Vec<usize> = Vec::new();
    let query_len = query.len();
    let protein_len = protein.len();

    // Match branch: both pointers advance, consumes 1 protein char.
    if p_idx >= 0
        && (p_idx as usize) < protein_len
        && query[q_idx as usize] == protein[p_idx as usize]
    {
        for consumed in dfs(query, q_idx + direction, protein, p_idx + direction, indels_left, direction) {
            all_paths.push(consumed + 1);
        }
    }

    if indels_left > 0 {
        // Deletion branch: query pointer advances, protein stays, 0 protein chars consumed.
        if !is_terminal_deletion(q_idx, query_len, p_idx, protein_len) {
            all_paths.extend(dfs(
                query, q_idx + direction, protein, p_idx, indels_left - 1, direction,
            ));
        }
        // Insertion branch: protein pointer advances, query stays, 1 protein char consumed.
        if p_idx >= 0 && (p_idx as usize) < protein_len {
            for consumed in dfs(
                query, q_idx, protein, p_idx + direction, indels_left - 1, direction,
            ) {
                all_paths.push(consumed + 1);
            }
        }
    }

    all_paths
}

fn extend_bidirectional(
    query: &[u8],
    q_seed_start: usize,
    p_hit_idx: usize,
    protein: &[u8],
    k: usize,
    max_indels: usize,
) -> Vec<(usize, Vec<u8>)> {
    let mut seen: HashSet<(usize, Vec<u8>)> = HashSet::new();
    let mut results: Vec<(usize, Vec<u8>)> = Vec::new();

    for r_budget in 0..=max_indels {
        let l_budget = max_indels - r_budget;

        let r_paths = dfs(
            query, (q_seed_start + k) as isize, protein, (p_hit_idx + k) as isize,
            r_budget, 1,
        );
        let l_paths = dfs(
            query, q_seed_start as isize - 1, protein, p_hit_idx as isize - 1,
            l_budget, -1,
        );

        for &r_consumed in &r_paths {
            for &l_consumed in &l_paths {
                if l_consumed > p_hit_idx {
                    continue;
                }
                let start = p_hit_idx - l_consumed;
                let end = p_hit_idx + k + r_consumed;
                if end > protein.len() {
                    continue;
                }
                let matched = protein[start..end].to_vec();
                let key = (start, matched.clone());
                if !seen.contains(&key) {
                    seen.insert(key);
                    results.push((start, matched));
                }
            }
        }
    }

    results
}

fn indel_search_peptide(
    query_id: &str,
    peptide: &str,
    indels_allowed: usize,
    index: &PepIndex,
) -> Vec<Vec<String>> {
    let pep_bytes = peptide.as_bytes();
    let peptide_len = pep_bytes.len();
    let k = index.k;

    if peptide_len < k {
        return vec![make_empty_row(query_id, peptide)];
    }

    let seeds = minimal_coverage_seeds(pep_bytes, k);
    let mut seen: HashSet<(usize, usize, Vec<u8>)> = HashSet::new();
    let mut matches: Vec<Vec<String>> = Vec::new();

    for (seed_bytes, q_seed_start) in &seeds {
        if let Some(hit_positions) = index.lookup(seed_bytes) {
            for encoded_hit in hit_positions {
                let prot_num = (encoded_hit / PROTEIN_INDEX_MULTIPLIER) as usize;
                if prot_num == 0 || prot_num > index.num_proteins {
                    continue;
                }
                let local_p_idx = (encoded_hit % PROTEIN_INDEX_MULTIPLIER) as usize;

                let prot_base = index.protein_offset(prot_num) as usize;
                let prot_end = if prot_num < index.num_proteins {
                    index.protein_offset(prot_num + 1) as usize
                } else {
                    index.seq_len
                };
                let protein = &index.mmap[index.seq_offset + prot_base..index.seq_offset + prot_end];

                for (true_start, matched) in
                    extend_bidirectional(pep_bytes, *q_seed_start, local_p_idx, protein, k, indels_allowed)
                {
                    let dedup_key = (prot_num, true_start, matched.clone());
                    if seen.contains(&dedup_key) {
                        continue;
                    }
                    seen.insert(dedup_key);

                    let matched_str = String::from_utf8_lossy(&matched).into_owned();
                    let n_indels = matched.len().abs_diff(peptide_len);
                    let meta = index.get_metadata(prot_num);

                    matches.push(vec![
                        query_id.to_string(),
                        peptide.to_string(),
                        matched_str,
                        meta[0].clone(), meta[1].clone(), meta[2].clone(),
                        meta[3].clone(), meta[4].clone(),
                        n_indels.to_string(),
                        String::from("[]"),
                        (true_start + 1).to_string(),
                        (true_start + matched.len()).to_string(),
                        meta[5].clone(), meta[6].clone(), meta[7].clone(), meta[8].clone(),
                    ]);
                }
            }
        }
    }

    if matches.is_empty() {
        vec![make_empty_row(query_id, peptide)]
    } else {
        matches
    }
}

pub(crate) fn run_indel(
    pepidx_path: &str,
    peptides: Vec<(String, String)>,
    indels_allowed: usize,
) -> Vec<Vec<String>> {
    let index = PepIndex::open(pepidx_path);

    let all_results: Vec<Vec<Vec<String>>> = peptides
        .par_iter()
        .map(|(query_id, peptide)| {
            indel_search_peptide(query_id, peptide, indels_allowed, &index)
        })
        .collect();

    all_results.into_iter().flatten().collect()
}

pub(crate) fn run(pepidx_path: &str, peptides: Vec<(String, String)>, k: usize, max_mismatches: usize) -> Vec<Vec<String>> {
    let index = PepIndex::open(pepidx_path);

    let all_results: Vec<Vec<Vec<String>>> = peptides
        .par_iter()
        .map(|(query_id, peptide)| {
            search_peptide(query_id, peptide, k, max_mismatches, &index)
        })
        .collect();

    all_results.into_iter().flatten().collect()
}
