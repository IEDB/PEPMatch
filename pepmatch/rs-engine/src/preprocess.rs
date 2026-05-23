use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter, Read, Seek, SeekFrom};

const PROTEIN_INDEX_MULTIPLIER: u64 = 100_000_000;

struct Protein {
    header: String,
    sequence: String,
}

struct Metadata {
    protein_id: String,
    protein_name: String,
    species: String,
    taxon_id: String,
    gene: String,
    pe_level: String,
    sequence_version: String,
    gene_priority: String,
    swissprot: String,
}

fn parse_fasta(path: &str) -> Vec<Protein> {
    let file = File::open(path).expect("Cannot open FASTA file");
    let reader = BufReader::new(file);
    let mut proteins: Vec<Protein> = Vec::new();
    let mut current_header = String::new();
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line.expect("Cannot read line");
        if line.starts_with('>') {
            if !current_header.is_empty() {
                proteins.push(Protein {
                    header: current_header.clone(),
                    sequence: current_seq.clone(),
                });
                current_seq.clear();
            }
            current_header = line[1..].to_string();
        } else {
            current_seq.push_str(line.trim());
        }
    }
    if !current_header.is_empty() {
        proteins.push(Protein {
            header: current_header,
            sequence: current_seq,
        });
    }
    proteins
}

fn extract_between<'a>(header: &'a str, start_tag: &str, end_tag: &str) -> &'a str {
    if let Some(start) = header.find(start_tag) {
        let value_start = start + start_tag.len();
        let rest = &header[value_start..];
        let end = rest.find(end_tag).unwrap_or(rest.len());
        rest[..end].trim()
    } else {
        ""
    }
}

fn extract_field<'a>(header: &'a str, prefix: &str) -> &'a str {
    if let Some(start) = header.find(prefix) {
        let value_start = start + prefix.len();
        let rest = &header[value_start..];
        let end = rest.find(|c: char| c.is_whitespace()).unwrap_or(rest.len());
        rest[..end].trim()
    } else {
        ""
    }
}

fn extract_protein_id(header: &str) -> String {
    if let Some(first_pipe) = header.find('|') {
        let rest = &header[first_pipe + 1..];
        if let Some(second_pipe) = rest.find('|') {
            return rest[..second_pipe].to_string();
        }
    }
    header.split_whitespace().next().unwrap_or("").to_string()
}

fn extract_protein_name(header: &str) -> String {
    if let Some(space_pos) = header.find(' ') {
        let after_id = &header[space_pos + 1..];
        if let Some(os_pos) = after_id.find(" OS=") {
            return after_id[..os_pos].trim().to_string();
        }
        return after_id.trim().to_string();
    }
    String::new()
}

fn is_swissprot(header: &str) -> bool {
    header.starts_with("sp|")
}

fn has_isoform_dash(protein_id: &str) -> bool {
    protein_id.contains('-')
}

fn base_accession(protein_id: &str) -> &str {
    if let Some(dash_pos) = protein_id.find('-') {
        &protein_id[..dash_pos]
    } else {
        protein_id
    }
}

fn extract_all_metadata(proteins: &[Protein]) -> Vec<Metadata> {
    let mut canonical_pe: HashMap<String, String> = HashMap::new();

    for protein in proteins {
        let pid = extract_protein_id(&protein.header);
        if !has_isoform_dash(&pid) {
            let pe = extract_field(&protein.header, "PE=");
            if !pe.is_empty() {
                canonical_pe.insert(pid.clone(), pe.to_string());
            }
        }
    }

    let mut seen_genes: HashSet<String> = HashSet::new();
    let mut metadata_list: Vec<Metadata> = Vec::with_capacity(proteins.len());

    for protein in proteins {
        let h = &protein.header;
        let protein_id = extract_protein_id(h);
        let protein_name = extract_protein_name(h);
        let species = extract_between(h, "OS=", " OX=").to_string();
        let taxon_id = extract_field(h, "OX=").to_string();
        let gene = extract_field(h, "GN=").to_string();
        let sv = extract_field(h, "SV=");
        let sequence_version = if sv.is_empty() { "1".to_string() } else { sv.to_string() };
        let swissprot = if is_swissprot(h) { "1".to_string() } else { "0".to_string() };

        let mut pe_level = extract_field(h, "PE=").to_string();
        if has_isoform_dash(&protein_id) && (pe_level.is_empty() || pe_level == "0") {
            let base = base_accession(&protein_id);
            if let Some(canonical) = canonical_pe.get(base) {
                pe_level = canonical.clone();
            }
        }
        if pe_level.is_empty() {
            pe_level = "0".to_string();
        }

        let gene_priority = if is_swissprot(h)
            && !has_isoform_dash(&protein_id)
            && !gene.is_empty()
            && !seen_genes.contains(&gene)
        {
            seen_genes.insert(gene.clone());
            "1".to_string()
        } else {
            "0".to_string()
        };

        metadata_list.push(Metadata {
            protein_id,
            protein_name,
            species,
            taxon_id,
            gene,
            pe_level,
            sequence_version,
            gene_priority,
            swissprot,
        });
    }

    metadata_list
}

fn write_length_prefixed(buf: &mut Vec<u8>, s: &str) {
    let bytes = s.as_bytes();
    buf.extend_from_slice(&(bytes.len() as u16).to_le_bytes());
    buf.extend_from_slice(bytes);
}

fn write_kmer_chunk(kmer_map: HashMap<u128, Vec<u64>>, chunk_path: &str, k: usize) -> u64 {
    let mut sorted: Vec<(u128, Vec<u64>)> = kmer_map.into_iter().collect();
    sorted.sort_unstable_by_key(|a| a.0);

    let file = File::create(chunk_path).expect("Cannot create chunk file");
    let mut writer = BufWriter::new(file);
    let mut count: u64 = 0;

    for (kmer, positions) in &sorted {
        let kmer_bytes = &kmer.to_be_bytes()[16 - k..];
        writer.write_all(kmer_bytes).unwrap();
        writer.write_all(&(positions.len() as u32).to_le_bytes()).unwrap();
        for &pos in positions {
            writer.write_all(&pos.to_le_bytes()).unwrap();
        }
        count += 1;
    }

    writer.flush().unwrap();
    count
}

struct ChunkReader {
    reader: BufReader<File>,
    k: usize,
    current_kmer: Option<u128>,
    current_positions: Vec<u64>,
}

impl ChunkReader {
    fn new(path: &str, k: usize) -> Self {
        let file = File::open(path).expect("Cannot open chunk file");
        let mut cr = ChunkReader {
            reader: BufReader::new(file),
            k,
            current_kmer: None,
            current_positions: Vec::new(),
        };
        cr.advance();
        cr
    }

    fn advance(&mut self) {
        let mut kmer_bytes = vec![0u8; self.k];
        if self.reader.read_exact(&mut kmer_bytes).is_err() {
            self.current_kmer = None;
            self.current_positions.clear();
            return;
        }

        let mut kmer_val: u128 = 0;
        for &b in &kmer_bytes {
            kmer_val = (kmer_val << 8) | b as u128;
        }

        let mut count_bytes = [0u8; 4];
        self.reader.read_exact(&mut count_bytes).unwrap();
        let count = u32::from_le_bytes(count_bytes) as usize;

        self.current_positions.clear();
        for _ in 0..count {
            let mut pos_bytes = [0u8; 8];
            self.reader.read_exact(&mut pos_bytes).unwrap();
            self.current_positions.push(u64::from_le_bytes(pos_bytes));
        }

        self.current_kmer = Some(kmer_val);
    }
}

fn write_pepidx(path: &str, k: usize, proteins: &[Protein], metadata_list: &[Metadata]) {
    let file = File::create(path).expect("Cannot create output file");
    let mut writer = BufWriter::new(file);

    let mut concat_seq = Vec::new();
    let mut protein_offsets: Vec<u64> = Vec::with_capacity(proteins.len());
    for protein in proteins {
        protein_offsets.push(concat_seq.len() as u64);
        concat_seq.extend_from_slice(protein.sequence.as_bytes());
    }

    let chunk_size = 50_000;
    let mut chunk_paths: Vec<String> = Vec::new();
    let mut total_positions: u64 = 0;

    let mut chunk_start = 0;
    while chunk_start < proteins.len() {
        let chunk_end = std::cmp::min(chunk_start + chunk_size, proteins.len());
        let mut kmer_map: HashMap<u128, Vec<u64>> = HashMap::new();

        for protein_count in chunk_start..chunk_end {
            let seq = proteins[protein_count].sequence.as_bytes();
            if seq.len() < k {
                continue;
            }
            let protein_number = (protein_count + 1) as u64;
            for j in 0..=(seq.len() - k) {
                let mut kmer_val: u128 = 0;
                for &b in &seq[j..j + k] {
                    kmer_val = (kmer_val << 8) | b as u128;
                }
                let idx = protein_number * PROTEIN_INDEX_MULTIPLIER + j as u64;
                kmer_map.entry(kmer_val).or_default().push(idx);
                total_positions += 1;
            }
        }

        let chunk_path = format!("{}.chunk_{}", path, chunk_paths.len());
        write_kmer_chunk(kmer_map, &chunk_path, k);
        chunk_paths.push(chunk_path);
        chunk_start = chunk_end;
    }

    let mut metadata_buf: Vec<u8> = Vec::new();
    let mut metadata_offsets: Vec<u64> = Vec::with_capacity(proteins.len());
    for meta in metadata_list {
        metadata_offsets.push(metadata_buf.len() as u64);
        write_length_prefixed(&mut metadata_buf, &meta.protein_id);
        write_length_prefixed(&mut metadata_buf, &meta.protein_name);
        write_length_prefixed(&mut metadata_buf, &meta.species);
        write_length_prefixed(&mut metadata_buf, &meta.taxon_id);
        write_length_prefixed(&mut metadata_buf, &meta.gene);
        write_length_prefixed(&mut metadata_buf, &meta.pe_level);
        write_length_prefixed(&mut metadata_buf, &meta.sequence_version);
        write_length_prefixed(&mut metadata_buf, &meta.gene_priority);
        write_length_prefixed(&mut metadata_buf, &meta.swissprot);
    }

    let mut readers: Vec<ChunkReader> = chunk_paths.iter()
        .map(|p| ChunkReader::new(p, k))
        .collect();

    let merged_path = format!("{}.merged", path);
    let merged_file = File::create(&merged_path).expect("Cannot create merged file");
    let mut merged_writer = BufWriter::new(merged_file);
    let mut num_unique_kmers: u64 = 0;

    loop {
        let mut min_kmer: Option<u128> = None;
        for reader in &readers {
            if let Some(kmer) = reader.current_kmer {
                if min_kmer.is_none() || kmer < min_kmer.unwrap() {
                    min_kmer = Some(kmer);
                }
            }
        }

        let min_kmer = match min_kmer {
            Some(k) => k,
            None => break,
        };

        let mut all_positions: Vec<u64> = Vec::new();
        for reader in &mut readers {
            if reader.current_kmer == Some(min_kmer) {
                all_positions.extend_from_slice(&reader.current_positions);
                reader.advance();
            }
        }

        let kmer_bytes = &min_kmer.to_be_bytes()[16 - k..];
        merged_writer.write_all(kmer_bytes).unwrap();
        merged_writer.write_all(&(all_positions.len() as u32).to_le_bytes()).unwrap();
        for &pos in &all_positions {
            merged_writer.write_all(&pos.to_le_bytes()).unwrap();
        }

        num_unique_kmers += 1;
    }
    merged_writer.flush().unwrap();
    drop(merged_writer);

    for chunk_path in &chunk_paths {
        std::fs::remove_file(chunk_path).ok();
    }

    let mut merged_reader = ChunkReader::new(&merged_path, k);
    let mut kmer_offsets: Vec<u64> = Vec::with_capacity(num_unique_kmers as usize);
    let mut kmer_data: Vec<u8> = Vec::new();

    loop {
        match merged_reader.current_kmer {
            None => break,
            Some(kmer) => {
                kmer_offsets.push(kmer_data.len() as u64);
                let kmer_bytes = &kmer.to_be_bytes()[16 - k..];
                kmer_data.extend_from_slice(kmer_bytes);
                kmer_data.extend_from_slice(&(merged_reader.current_positions.len() as u32).to_le_bytes());
                for &pos in &merged_reader.current_positions {
                    kmer_data.extend_from_slice(&pos.to_le_bytes());
                }
                merged_reader.advance();
            }
        }
    }

    std::fs::remove_file(&merged_path).ok();
    writer.write_all(b"PEPIDX05").unwrap();
    writer.write_all(&(k as u8).to_le_bytes()).unwrap();
    writer.write_all(&(proteins.len() as u32).to_le_bytes()).unwrap();
    writer.write_all(&(concat_seq.len() as u64).to_le_bytes()).unwrap();
    writer.write_all(&num_unique_kmers.to_le_bytes()).unwrap();
    writer.write_all(&total_positions.to_le_bytes()).unwrap();
    writer.write_all(&(metadata_buf.len() as u64).to_le_bytes()).unwrap();
    writer.write_all(&concat_seq).unwrap();
    for &offset in &protein_offsets {
        writer.write_all(&offset.to_le_bytes()).unwrap();
    }
    for &offset in &metadata_offsets {
        writer.write_all(&offset.to_le_bytes()).unwrap();
    }
    writer.write_all(&metadata_buf).unwrap();
    for &offset in &kmer_offsets {
        writer.write_all(&offset.to_le_bytes()).unwrap();
    }
    writer.write_all(&kmer_data).unwrap();
    writer.flush().unwrap();
}

pub(crate) fn run(fasta_path: &str, k: usize, output_path: &str) {
    let proteins = parse_fasta(fasta_path);
    let metadata_list = extract_all_metadata(&proteins);
    write_pepidx(output_path, k, &proteins, &metadata_list);
}
