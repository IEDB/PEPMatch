def _terminal_deletion_blocked(d, query_len, gap_left, gap_right, protein_len):
  """Blocks a deletion at the query's own boundary (no query-side context to
  confirm it's genuinely absent) or where its gap references an out-of-bounds
  protein index (indistinguishable from the sequence simply ending there).
  """
  if d == 0 or d == query_len - 1:
    return True
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
        if _terminal_deletion_blocked(d, query_len, start + d - 1, start + d, protein_len):
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
