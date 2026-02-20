use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};

const PROTEIN_INDEX_MULTIPLIER: u64 = 10_000_000;

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

fn write_pepidx(path: &str, k: usize, proteins: &[Protein], metadata_list: &[Metadata]) {
    let file = File::create(path).expect("Cannot create output file");
    let mut writer = BufWriter::new(file);

    let mut concat_seq = Vec::new();
    let mut protein_offsets: Vec<u64> = Vec::with_capacity(proteins.len());

    for protein in proteins {
        protein_offsets.push(concat_seq.len() as u64);
        concat_seq.extend_from_slice(protein.sequence.as_bytes());
    }

    let mut kmer_map: HashMap<Vec<u8>, Vec<u64>> = HashMap::new();
    let mut total_positions: u64 = 0;

    for (protein_count, protein) in proteins.iter().enumerate() {
        let seq = protein.sequence.as_bytes();
        if seq.len() < k {
            continue;
        }
        let protein_number = (protein_count + 1) as u64;
        for j in 0..=(seq.len() - k) {
            let kmer = seq[j..j + k].to_vec();
            let idx = protein_number * PROTEIN_INDEX_MULTIPLIER + j as u64;
            kmer_map.entry(kmer).or_default().push(idx);
            total_positions += 1;
        }
    }

    let mut sorted_kmers: Vec<(Vec<u8>, Vec<u64>)> = kmer_map.into_iter().collect();
    sorted_kmers.sort_by(|a, b| a.0.cmp(&b.0));

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

    let mut kmer_offsets: Vec<u64> = Vec::with_capacity(sorted_kmers.len());
    let mut kmer_data = Vec::new();

    for (kmer, positions) in &sorted_kmers {
        kmer_offsets.push(kmer_data.len() as u64);
        kmer_data.extend_from_slice(kmer);
        kmer_data.extend_from_slice(&(positions.len() as u32).to_le_bytes());
        for &pos in positions {
            kmer_data.extend_from_slice(&pos.to_le_bytes());
        }
    }

    writer.write_all(b"PEPIDX05").unwrap();
    writer.write_all(&(k as u8).to_le_bytes()).unwrap();
    writer.write_all(&(proteins.len() as u32).to_le_bytes()).unwrap();
    writer.write_all(&(concat_seq.len() as u64).to_le_bytes()).unwrap();
    writer.write_all(&(sorted_kmers.len() as u64).to_le_bytes()).unwrap();
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
