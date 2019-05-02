# A table mapping protein ids to names
CREATE TABLE protein_name (
	protein_id integer PRIMARY KEY,
	protein_name text,
);

# A breakdown of each protein postion and the nmer ID at that position
CREATE TABLE protein_peptide (
	protein_id integer REFERENCES protein_name(protein_id) ON UPDATE CASCADE ON DELETE CASCADE,
	position integer,
	nmer_id integer,
);

# indexing the nmer ID for fast lookups
# TODO: potentially make this a covering index with all relevant fields
CREATE INDEX nmer_idx ON protein_peptide(nmer_id);