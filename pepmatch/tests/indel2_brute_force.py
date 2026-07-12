def _out_of_bounds_deletion(p_idx, protein_len):
  """Blocks a deletion whose protein index is truly out of bounds — a database edge
  effect rather than an indel. p_idx 0 and protein_len-1 are real residues, allowed."""
  return p_idx < 0 or p_idx >= protein_len


def _align_deletions(query, q_idx, window, w_idx, dels_left, protein, w_offset):
  """query[q_idx:] aligns to window[w_idx:] using only deletions — no insertions, no
  mismatches. Keeping the two edit kinds in separate passes is what makes a match
  homogeneous: an alignment can never mix a deletion with an insertion."""
  if q_idx == len(query):
    return w_idx == len(window)

  if w_idx < len(window) and query[q_idx] == window[w_idx]:
    if _align_deletions(query, q_idx + 1, window, w_idx + 1, dels_left, protein, w_offset):
      return True

  if dels_left > 0:
    # A deletion at the query's own first/last residue has no query-side context to
    # confirm the residue is genuinely absent, so it is barred.
    query_terminal = q_idx == 0 or q_idx == len(query) - 1
    if not query_terminal and not _out_of_bounds_deletion(w_offset + w_idx, len(protein)):
      if _align_deletions(query, q_idx + 1, window, w_idx, dels_left - 1, protein, w_offset):
        return True

  return False


def _align_insertions(query, q_idx, window, w_idx, ins_left):
  """query[q_idx:] aligns to window[w_idx:] using only insertions."""
  if q_idx == len(query):
    return w_idx == len(window)

  if w_idx < len(window) and query[q_idx] == window[w_idx]:
    if _align_insertions(query, q_idx + 1, window, w_idx + 1, ins_left):
      return True

  # q_idx > 0 bars an insertion before the first query residue; one after the last is
  # unreachable because the base case returns as soon as the query is consumed. Both
  # would only pad an exact match with an arbitrary flanking residue.
  if ins_left > 0 and q_idx > 0 and w_idx < len(window):
    if _align_insertions(query, q_idx, window, w_idx + 1, ins_left - 1):
      return True

  return False


def brute_force_search(query, protein, max_indels=2):
  """Enumerate every (start, matched_sequence) where query aligns to a window of
  protein using at most max_indels homogeneous indels — insertions only or deletions
  only, never both.

  Checks each candidate window directly against the definition of a valid match (no
  seeding, no bidirectional extension), so it can serve as an independent oracle for
  the Rust engine.
  """
  query_len = len(query)
  protein_len = len(protein)
  results = set()

  for start in range(protein_len):
    for window_len in range(max(0, query_len - max_indels), query_len + max_indels + 1):
      end = start + window_len
      if end > protein_len:
        continue
      window = protein[start:end]
      if window_len <= query_len:
        if _align_deletions(
          query, 0, window, 0, query_len - window_len, protein, start
        ):
          results.add((start, window))
      elif _align_insertions(query, 0, window, 0, window_len - query_len):
        results.add((start, window))

  return results
