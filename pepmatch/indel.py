def _is_terminal_deletion(q_idx, query_len, p_idx, protein_len):
  """Return True if a deletion at this position is at a query or protein boundary.

  Terminal deletions are forbidden because a deletion at the first or last
  query position produces a match indistinguishable from a shorter peptide,
  and a deletion at a protein boundary arises from an edge effect rather
  than a true biological indel event.
  """
  if q_idx == 0 or q_idx == query_len - 1:
    return True
  if p_idx <= 0 or p_idx >= protein_len - 1:
    return True
  return False


def _dfs(query, q_idx, protein, p_idx, indels_left, direction):
  """Exhaustive depth-first search extending one direction from a seed hit.

  Explores three branches at each step: a match branch (both pointers
  advance), a deletion branch (query pointer advances, protein stays,
  consumes 0 protein chars), and an insertion branch (protein pointer
  advances, query stays, consumes 1 protein char). Returns a list of
  integers, each the number of protein characters consumed by one valid
  complete alignment path.

  Args:
    query: the full query peptide string.
    q_idx: current query position, advances by direction.
    protein: the full protein sequence string.
    p_idx: current protein position, advances by direction.
    indels_left: remaining indel budget for this direction.
    direction: +1 for right extension, -1 for left extension.
  """
  if (direction == 1 and q_idx >= len(query)) or \
     (direction == -1 and q_idx < 0):
    return [0]

  all_paths = []

  if 0 <= p_idx < len(protein) and query[q_idx] == protein[p_idx]:
    for consumed in _dfs(query, q_idx + direction, protein, p_idx + direction, indels_left, direction):
      all_paths.append(consumed + 1)

  if indels_left > 0:
    if not _is_terminal_deletion(q_idx, len(query), p_idx, len(protein)):
      all_paths.extend(
        _dfs(query, q_idx + direction, protein, p_idx, indels_left - 1, direction)
      )

    if 0 <= p_idx < len(protein):
      for consumed in _dfs(query, q_idx, protein, p_idx + direction, indels_left - 1, direction):
        all_paths.append(consumed + 1)

  return all_paths


def extend_bidirectional(query, q_seed_start, p_hit_idx, protein, k, max_indels=1):
  """Verify and enumerate all indel-consistent matches from a seed hit.

  For each allocation of the indel budget between left and right extensions,
  runs _dfs in both directions and combines results. Returns all unique
  (start_0based, matched_sequence) pairs found at this seed hit position.

  Args:
    query: the full query peptide string.
    q_seed_start: 0-based start of the seed k-mer within the query.
    p_hit_idx: 0-based start of the seed hit within protein.
    protein: the full protein sequence string.
    k: seed k-mer length.
    max_indels: maximum total indels allowed across both extensions.
  """
  seen = set()
  results = []

  for r_budget in range(max_indels + 1):
    l_budget = max_indels - r_budget

    r_paths = _dfs(query, q_seed_start + k, protein, p_hit_idx + k, r_budget, direction=1)
    l_paths = _dfs(query, q_seed_start - 1, protein, p_hit_idx - 1, l_budget, direction=-1)

    for r_consumed in r_paths:
      for l_consumed in l_paths:
        start = p_hit_idx - l_consumed
        end = p_hit_idx + k + r_consumed
        matched_seq = protein[start:end]
        key = (start, matched_seq)
        if key not in seen:
          seen.add(key)
          results.append(key)

  return results


def minimal_coverage_seeds(query, k):
  """Return non-overlapping k-mer seeds plus a final overlapping seed for
  full query coverage.

  The pigeonhole principle guarantees that for a query of length L with at
  most d indels, at least one of floor(L / (d+1)) non-overlapping seeds
  must hit exactly. Appends one final seed anchored at the query end if
  the last strided seed does not reach it.

  Args:
    query: the query peptide string.
    k: seed k-mer length.
  """
  seeds = []
  query_len = len(query)
  for j in range(0, query_len - k + 1, k):
    seeds.append((query[j:j + k], j))
  last_start = query_len - k
  if seeds and seeds[-1][1] != last_start:
    seeds.append((query[last_start:], last_start))
  return seeds