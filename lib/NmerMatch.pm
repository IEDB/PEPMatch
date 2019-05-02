package NmerMatch;
require Exporter;
@ISA = qw(Exporter);
@EXPORT_OK = qw(read_fasta build_db);

use Bio::SeqIO;

our $VERSION = '0.01';

# given the name of a fasta file
# return a reference to a an array of hashrefs of sequences {id => 1, seq=>DFDSFSDF}
sub read_fasta {

	my $fasta = shift;

	my $seqin = Bio::SeqIO->new(-file => $fasta, '-format' => 'Fasta');

	my @seqs;
	my %seen_id;
	my $seq_num = 1;
	while (my $seq = $seqin->next_seq()) {
		my $seq_id = $seq->display_id();
		if ($seen{$seq_id}) {
			warn "Multiple sequences with same identifier: $seq_id";
			warn "Appending number to id";
			$seq_id .= ".$seq_num";
		}
		push @seqs, { id => $seq_id, seq => $seq->seq() };
		$seen_id{$seq_id} = 1;
		$seq_num++;
	}

	return \@seqs;

}

# given:
# 1. an array of hashrefs of protein sequences as returned by read_fasta
# 2. a peptide length
# do:
# find all unique nmer sequences
# catalog the proteins in the input file, recording the unique nmers and their positions
# return:
# a reference to the unqiue nmer hash
# a reference to the protein catalog (maybe this should be in an actual DB?)
sub build_db {

	my $db_seq_ref = shift;
	my $peptide_length = shift;

	my $nmer_id = 0;
	my %unique_nmer_id;
	my %db_catalog; # indexed first by nmer_id, then sequence id, an array
	                # of positions in that sequence in which the nmer occurs
	my @nmer_ids;
    
	foreach my $p (@$db_seq_ref) {

		my $protein_nmers_ref = break_protein($p->{seq}, $peptide_length);

		foreach my $start (0 .. $#$protein_nmers_ref) {
	    
			my $s = $protein_nmers_ref->[$start];

	    	if (!exists $unique_nmer_id{$s}) {
		  		$unique_nmer_id{$s} = ++$nmer_id;
		  	}
		  	# add the position of the nmer onto an array
		  	push @{$db_catalog{$nmer_id}{$p->{id}}}, $start;
	    }
	}

	return (\%unique_nmer_id, \%db_catalog);

}

# given a protein sequence and peptide length, return an
# array ref of overlapping peptides, indexed by position
sub break_protein {

	my $protein_seq = shift;
	my $peptide_length = shift;

	my $windows = length($protein_seq) - $peptide_length;
    
	my @nmers;
    
	foreach my $start (0 .. $windows) {
	    my $s = substr $protein_seq, $start, $peptide_length;
        push @nmers, $s;
	}

	return \@nmers;

}

# given:
# 1. the name of a database containing unqiue nmers
# do:
# load it into a hash
sub load_db {

}

## DO we need this?
# given
# 1. a pointer to a DB OR hash containing the db database
# 2. a fasta file of peptides or list of peptides
# do:
# identify all unique peptides in the list
# append them to the db if they do not already exist
sub append_db {


}

1;