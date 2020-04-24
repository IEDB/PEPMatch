-- NOTE that a given nmer can be at multiple positions in a protein,
-- so more than 1 row for a given nmer_id & protein_id may exist

-- A breakdown of each protein postion and the nmer ID at that position
DROP TABLE IF EXISTS protein_peptide;
CREATE TABLE protein_peptide (
	nmer_id integer,
	protein_id integer,
	position integer
);

.separator ,
.import [% protein_peptides_csv %] protein_peptide

-- indexing the nmer id & protein_id for fast lookups
CREATE INDEX nmer_protein_idx ON protein_peptide(nmer_id, protein_id);