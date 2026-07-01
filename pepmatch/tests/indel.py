def _terminal_deletion_blocked(gap_left, gap_right, protein_len):
  """A deletion at query position d found at protein offset `start` has its gap
  sitting between protein indices (start+d-1) and (start+d). Both neighbors are
  real, existing residues as long as they fall within [0, protein_len) — only
  an out-of-bounds reference means the gap is indistinguishable from the
  sequence simply ending there rather than a genuine biological indel.
  """
  return gap_left < 0 or gap_right >= protein_len


def brute_force_search(query, protein):
  """Enumerate every (start, matched_sequence) pair where query aligns to a
  window of protein using at most one single-residue insertion or deletion.

  Checks each candidate directly against the definition of a valid 1-indel
  match (no seeding, no recursive extension) so it can serve as an
  independent oracle for the Rust/seed-based implementation.
  """
  query_len = len(query)
  protein_len = len(protein)
  results = set()

  for start in range(protein_len - query_len + 1):
    if protein[start:start + query_len] == query:
      results.add((start, protein[start:start + query_len]))

  window_len = query_len - 1
  for d in range(query_len):
    shortened = query[:d] + query[d + 1:]
    for start in range(protein_len - window_len + 1):
      if protein[start:start + window_len] == shortened:
        if _terminal_deletion_blocked(start + d - 1, start + d, protein_len):
          continue
        results.add((start, shortened))

  # i in [1, query_len) keeps the inserted residue interior to the alignment;
  # an insertion at i=0 or i=query_len would just pad an exact match with an
  # arbitrary flanking residue and isn't a real edit.
  window_len = query_len + 1
  for i in range(1, query_len):
    for start in range(protein_len - window_len + 1):
      window = protein[start:start + window_len]
      if window[:i] == query[:i] and window[i + 1:] == query[i:]:
        results.add((start, window))

  return results
