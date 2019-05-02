package NmerMatch;
require Exporter;
@ISA = qw(Exporter);
@EXPORT_OK = qw(read_fasta build_catalog retrieve_catalog get_catalog_info read_query_file);

use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;
use Storable qw(nstore retrieve);
use JSON::XS;
use POSIX qw(strftime);
use Inline 'C';
use Inline (C => Config =>
            OPTIMIZE => '-O3',
            #CCFLAGSEX => '-march=x86-64',
            );

#TODO: test saving & loading of db_catalog & unique_nmer_id versus storing in DB & retrieval
#      I suspect the latter will be faster for the db_catalog, even for sqlite

our $VERSION = '0.01';

# the backend for the protein peptide catalog
# by default it is an in-memory hash, but can be changed
my $catalog_backend = 'memory';
my $catalog_source; # this will be set when the input fasta file is read so that it can be
                    # saved when the catalog is built
my $catalog_info;   # this will become a hashref of metadata about the catalog upon build/restore

# given the name of a query file containing a list (or fasta)
# of peptides, return an array reference to the peptides
sub read_query_file {

	my $input_file = shift;

	my @peptides;
	my $peptide_length;
	open my $INF, '<', $input_file or die "Can't open input file ($input_file) for reading: $!";
	while (<$INF>) {
		chomp;
		push @peptides, $_;
		my $cur_peptide_length = length($peptides[-1]);
		# set the peptide length to the length of the first peptide_length
		# all subsequent peptides are checked that they match in length
		if (!defined $peptide_length) {
			$peptide_length = $cur_peptide_length;
		}
		if ($cur_peptide_length ne $peptide_length) {
			die "Peptides of different lengths ($peptide_length and $cur_peptide_length) found in query file: $input_file";
		}
	}
	close $INF;

	if ($peptide_length ne $catalog_info->{peptide_length}) {
		die "Peptide length ($peptide_length) differs from catalog peptide length (" . $catalog_info->{peptide_length} . ')';
	}

	return \@peptides;
}

# given the name of a fasta file
# return a reference to a an array of hashrefs of sequences {id => 1, seq=>DFDSFSDF}
sub read_fasta {

	my $fasta = shift;

	# setting the catalog source globally
	$catalog_source = $fasta;

	my $seqin = Bio::SeqIO->new(-file => $fasta, '-format' => 'Fasta');

	my @seqs;
	my %seen_id;
	my $seq_num = 1;
	while (my $seq = $seqin->next_seq()) {
		my $seq_id = $seq->display_id();
		if ($seen_id{$seq_id}) {
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
# serialize the catalog according to the selected backend
# return:
# a reference to the unqiue nmer hash
# a reference to the protein catalog (maybe this should be in an actual DB?)
sub build_catalog {

	my $catalog_seq_ref = shift;
	my $peptide_length = shift;
	my $catalog_name = shift;

	my $nmer_id = 0;
	my %unique_nmer_id;
	my %nmer_catalog; # indexed first by nmer_id, then sequence id, an array
	                # of positions in that sequence in which the nmer occurs
	my @nmer_ids;
    
	foreach my $p (@$catalog_seq_ref) {

		my $protein_nmers_ref = break_protein($p->{seq}, $peptide_length);

		foreach my $start (0 .. $#$protein_nmers_ref) {
	    
			my $s = $protein_nmers_ref->[$start];

			# TODO: Here is where we branch depending on the catlog backend
			# if we find the current nmer in the unique nmers, use the ID, otherwise create a new id
			if (!exists $unique_nmer_id{$s}) {
				$unique_nmer_id{$s} = ++$nmer_id;
			}
			my $cur_nmer_id = $unique_nmer_id{$s};

		  	# add the position of the nmer onto an array
		  	push @{$nmer_catalog{$cur_nmer_id}{$p->{id}}}, $start;
	    }
	}

	# now serialize the uniqe nmers & nmer_catalog, and include a description
	# TODO: this should be refactored and dependent upon the backend
	my ($unique_nmers_file, $protein_peptides_file, $catalog_info_file) = get_catalog_file_names($catalog_name);

	if (!-d $catalog_name) {
		print "Creating directory for catalog: $catalog_name\n";
		mkdir $catalog_name or die "Failed creating directory for catalog: $!";
	}

	print "Serializing uniqe nmers to: $unique_nmers_file\n";
	nstore \%unique_nmer_id, $unique_nmers_file;
	print "Done\n";
	print "Serializing protein peptides to: $protein_peptides_file\n";
	nstore \%nmer_catalog, $protein_peptides_file;
	print "Done\n";
	print "Storing catalog info in: $catalog_info_file\n";
	$catalog_info =    { catalog_name   => $catalog_name,
		                 data_source    => $catalog_source,
	                     peptide_length => $peptide_length,
	                     build_date     => strftime "%F %R", localtime,
	                   };
	my $catalog_json = JSON::XS->new
	                           ->utf8
	                           ->pretty(1)
	                           ->canonical(1)
	                           ->encode($catalog_info);
	open my $CAT, '>', $catalog_info_file or die "Can't open catalog info file ($catalog_info_file): $!";
	print $CAT $catalog_json, "\n";
	close $CAT;

	return (\%unique_nmer_id, \%nmer_catalog);

}

# read in the catalog info file & return a hashref
sub get_catalog_info {

	my $catalog_name = shift;

	my ($unique_nmers_file, $protein_peptides_file, $catalog_info_file) = get_catalog_file_names($catalog_name);

	my $catalog_info_string;
	print "Reading catalog info from: $catalog_info_file\n";
	open my $CAT, '<', $catalog_info_file or die "Can't open catalog info file ($catalog_info_file): $!";
	while (<$CAT>) {
		$catalog_info_string .= $_;
	}
	close $CAT;
	$catalog_info = JSON::XS->new->utf8->decode($catalog_info_string);
	print "Catalog info\n";
	print Dumper $catalog_info;

	return $catalog_info;

}


# given a catalog name, return the names of the unique nmers & protein peptides
# files
sub get_catalog_file_names {

	my $catalog_name = shift;

	my $unique_nmers_file = $catalog_name . '/unique_nmers.store';
	my $protein_peptides_file = $catalog_name . '/protein_peptides.store';
	my $catalog_info_file = $catalog_name . '/catalog_info.json';

	return ($unique_nmers_file, $protein_peptides_file, $catalog_info_file);
}


# given:
# a catalog name
# do:
# retrieve the catalog and return references to the uniqe nmer ids and nmer_catalog
# depending upon the backend
sub retrieve_catalog {

	my $catalog_name = shift;

	# TODO: the functionlity will change depending on the backend
	my ($unique_nmers_file, $protein_peptides_file) = get_catalog_file_names($catalog_name);

	print "Restoring unique nmer reference from: $unique_nmers_file\n";
	my $unique_nmers_ref = retrieve($unique_nmers_file);
	print "Done\n";
	print "Restoring protein peptides reference from: $protein_peptides_file\n";
	my $protein_peptides_ref = retrieve($protein_peptides_file);
	print "Done\n";

	return ($unique_nmers_ref, $protein_peptides_ref);

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

__DATA__
__C__
  
  /* count mismatches in a string */
int count_mismatches(char* x, char* y, int max_mm) {                                           
    int i;
    int mismatches = 0;                                                                                                                                
    for(i=0; x[i] && y[i]; ++i) {                                               
      if(x[i] != y[i]) {                                                        
          mismatches++;
          if(mismatches > max_mm) {
            return max_mm + 1;
          }
      }                                                                         
    }                                                                           
    return mismatches;                                                         
 }