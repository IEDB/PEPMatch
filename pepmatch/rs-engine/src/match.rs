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

struct HitRecord {
    query_id: String,
    query_seq: String,
    matched: Option<String>,
    protein_num: Option<u32>,
    mismatches: Option<i64>,
    mutated: String,
    idx_start: Option<i64>,
    idx_end: Option<i64>,
}

fn mutated_positions(peptide: &str, matched: &str) -> String {
    let v: Vec<String> = peptide.as_bytes().iter()
        .zip(matched.as_bytes())
        .enumerate()
        .filter(|(_, (a, b))| a != b)
        .map(|(i, _)| (i + 1).to_string())
        .collect();
    format!("[{}]", v.join(", "))
}

fn hit_record(query_id: &str, peptide: &str, matched_seq: &str, mismatches: usize, encoded_start: u64) -> HitRecord {
    let protein_number = (encoded_start / PROTEIN_INDEX_MULTIPLIER) as u32;
    let position = (encoded_start % PROTEIN_INDEX_MULTIPLIER) as usize;
    HitRecord {
        query_id: query_id.to_string(),
        query_seq: peptide.to_string(),
        matched: Some(matched_seq.to_string()),
        protein_num: Some(protein_number),
        mismatches: Some(mismatches as i64),
        mutated: mutated_positions(peptide, matched_seq),
        idx_start: Some((position + 1) as i64),
        idx_end: Some((position + peptide.len()) as i64),
    }
}

fn miss_record(query_id: &str, peptide: &str) -> HitRecord {
    HitRecord {
        query_id: query_id.to_string(),
        query_seq: peptide.to_string(),
        matched: None,
        protein_num: None,
        mismatches: None,
        mutated: String::from("[]"),
        idx_start: None,
        idx_end: None,
    }
}

fn search_peptide(query_id: &str, peptide: &str, k: usize, max_mismatches: usize, index: &PepIndex) -> Vec<HitRecord> {
    if max_mismatches == 0 {
        let hits = exact_match(peptide, k, index);
        if hits.is_empty() {
            return vec![miss_record(query_id, peptide)];
        }
        hits.iter().map(|&hit| hit_record(query_id, peptide, peptide, 0, hit)).collect()
    } else {
        let hits = mismatch_match(peptide, k, max_mismatches, index);
        if hits.is_empty() {
            return vec![miss_record(query_id, peptide)];
        }
        hits.iter().map(|(matched_seq, mm, start)| hit_record(query_id, peptide, matched_seq, *mm, *start)).collect()
    }
}

pub(crate) type Columns = (
    Vec<String>, Vec<String>, Vec<Option<String>>, Vec<Option<u32>>,
    Vec<Option<i64>>, Vec<String>, Vec<Option<i64>>, Vec<Option<i64>>,
);

pub(crate) type MetaColumns = (
    Vec<u32>, Vec<String>, Vec<String>, Vec<String>, Vec<String>,
    Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<String>,
);

fn unzip_records(records: Vec<HitRecord>) -> Columns {
    let n = records.len();
    let mut qid = Vec::with_capacity(n);
    let mut qseq = Vec::with_capacity(n);
    let mut matched = Vec::with_capacity(n);
    let mut pnum = Vec::with_capacity(n);
    let mut mm = Vec::with_capacity(n);
    let mut mutated = Vec::with_capacity(n);
    let mut s = Vec::with_capacity(n);
    let mut e = Vec::with_capacity(n);
    for r in records {
        qid.push(r.query_id);
        qseq.push(r.query_seq);
        matched.push(r.matched);
        pnum.push(r.protein_num);
        mm.push(r.mismatches);
        mutated.push(r.mutated);
        s.push(r.idx_start);
        e.push(r.idx_end);
    }
    (qid, qseq, matched, pnum, mm, mutated, s, e)
}

pub(crate) fn run_discontinuous(
    pepidx_path: &str,
    epitopes: Vec<(String, Vec<(char, usize)>)>,
    max_mismatches: usize,
) -> Columns {
    let index = PepIndex::open(pepidx_path);
    let seq = &index.mmap[index.seq_offset..index.seq_offset + index.seq_len];

    let mut records: Vec<HitRecord> = Vec::new();

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
                let matched_str: String = residues.iter()
                    .map(|&(_, p)| format!("{}{}", seq[base + p - 1] as char, p))
                    .collect::<Vec<_>>()
                    .join(", ");
                let mutated: Vec<String> = residues.iter()
                    .filter(|&&(r, p)| seq[base + p - 1] as char != r)
                    .map(|&(_, p)| p.to_string())
                    .collect();

                records.push(HitRecord {
                    query_id: query_id.clone(),
                    query_seq: query_str.clone(),
                    matched: Some(matched_str),
                    protein_num: Some(pn as u32),
                    mismatches: Some(mismatches as i64),
                    mutated: format!("[{}]", mutated.join(", ")),
                    idx_start: residues.first().map(|&(_, p)| p as i64),
                    idx_end: residues.last().map(|&(_, p)| p as i64),
                });
            }
        }

        if !found {
            records.push(HitRecord {
                query_id: query_id.clone(),
                query_seq: query_str.clone(),
                matched: None,
                protein_num: None,
                mismatches: None,
                mutated: String::from("[]"),
                idx_start: None,
                idx_end: None,
            });
        }
    }

    unzip_records(records)
}

fn metadata_columns(index: &PepIndex) -> MetaColumns {
    let n = index.num_proteins;
    let mut pnum = Vec::with_capacity(n);
    let mut f: [Vec<String>; 9] = Default::default();
    for v in f.iter_mut() {
        v.reserve(n);
    }
    for pn in 1..=n {
        pnum.push(pn as u32);
        let m = index.get_metadata(pn);
        for (i, field) in m.into_iter().enumerate() {
            f[i].push(field);
        }
    }
    let [a, b, c, d, e, g, h, i, j] = f;
    (pnum, a, b, c, d, e, g, h, i, j)
}

pub(crate) fn run_metadata(pepidx_path: &str) -> MetaColumns {
    let index = PepIndex::open(pepidx_path);
    metadata_columns(&index)
}

pub(crate) fn run(pepidx_path: &str, peptides: Vec<(String, String)>, k: usize, max_mismatches: usize) -> Columns {
    let index = PepIndex::open(pepidx_path);
    let records: Vec<HitRecord> = peptides
        .par_iter()
        .flat_map_iter(|(query_id, peptide)| {
            search_peptide(query_id, peptide, k, max_mismatches, &index).into_iter()
        })
        .collect();
    unzip_records(records)
}

// ── Counts-only path (aggregate; O(unique queries), no per-hit materialization) ──

pub(crate) type CountColumns = (Vec<String>, Vec<String>, Vec<i64>, Vec<u64>);

/// Tally accepted matches per mismatch level for one peptide, mirroring
/// mismatch_match's dedup + validity walk exactly but without building the
/// matched sequence, metadata, or any per-hit row.
fn mismatch_count(peptide: &str, k: usize, max_mismatches: usize, index: &PepIndex, counts: &mut [u64]) {
    let pep_bytes = peptide.as_bytes();
    if pep_bytes.len() < k { return; }

    let num_kmers = pep_bytes.len() - k + 1;
    let peptide_len = pep_bytes.len();
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
                        let mut i = 0;
                        let mut valid = true;
                        while i < peptide_len {
                            if index.resolve(start + i as u64).is_none() { valid = false; break; }
                            i += k;
                        }
                        if valid {
                            seen.insert(start);
                            counts[mismatches] += 1;
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
                        let mut i = 0;
                        let mut valid = true;
                        while i < peptide_len {
                            if i + k > peptide_len {
                                let remaining = peptide_len - i;
                                let back = k - remaining;
                                if index.resolve(start + i as u64 - back as u64).is_none() { valid = false; break; }
                            } else if index.resolve(start + i as u64).is_none() {
                                valid = false; break;
                            }
                            i += k;
                        }
                        if valid {
                            seen.insert(start);
                            counts[mismatches] += 1;
                        }
                    }
                }
            }
        }
    }
}

fn count_peptide(query_id: &str, peptide: &str, k: usize, max_mismatches: usize, index: &PepIndex) -> Vec<(String, String, i64, u64)> {
    let mut counts = vec![0u64; max_mismatches + 1];
    if max_mismatches == 0 {
        counts[0] = exact_match(peptide, k, index).len() as u64;
    } else {
        mismatch_count(peptide, k, max_mismatches, index, &mut counts);
    }
    let mut out = Vec::new();
    for (mm, &c) in counts.iter().enumerate() {
        if c > 0 {
            out.push((query_id.to_string(), peptide.to_string(), mm as i64, c));
        }
    }
    out
}

pub(crate) fn run_counts(pepidx_path: &str, peptides: Vec<(String, String)>, k: usize, max_mismatches: usize) -> CountColumns {
    let index = PepIndex::open(pepidx_path);
    let rows: Vec<(String, String, i64, u64)> = peptides
        .par_iter()
        .flat_map_iter(|(query_id, peptide)| {
            count_peptide(query_id, peptide, k, max_mismatches, &index).into_iter()
        })
        .collect();
    let n = rows.len();
    let mut qid = Vec::with_capacity(n);
    let mut qseq = Vec::with_capacity(n);
    let mut mm = Vec::with_capacity(n);
    let mut cnt = Vec::with_capacity(n);
    for (a, b, c, d) in rows {
        qid.push(a);
        qseq.push(b);
        mm.push(c);
        cnt.push(d);
    }
    (qid, qseq, mm, cnt)
}
