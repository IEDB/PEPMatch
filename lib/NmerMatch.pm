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

# Workflow for searching
# 1. get catalog info to set the CATALOG_NAME
# 2. read in the query sequences
# 3. read in the catalog (UNIQUE_NMER_ID)
# 4. add unique nmers from query seqs to UNIQUE_NMER_ID
# 5. encode each nmer as an integer ID (UNIQUE_NMER) pointing to the sequence
# 6. create list of catalog nmer IDs
# 7. create list of query nmer IDs
# 8. determine which ones need to be compared



#TODO: think about the functionality that should be available outside this package
#      versus inside and potentially wrap up some functions
#TODO: similarly, think about some of the variables that should be globally set rather
#      than being passed from function to function
#TODO: test saving & loading of db_catalog & unique_nmer_id versus storing in DB & retrieval
#      I suspect the latter will be faster for the db_catalog, even for sqlite

our $VERSION = '0.01';

# the backend for the protein peptide catalog
# by default it is an in-memory hash, but can be changed
my $CATALOG_BACKEND = 'memory';
my $CATALOG_NAME;
my $CATALOG_SOURCE; # this will be set when the input fasta file is read so that it can be
                    # saved when the catalog is built
my @CATALOG_SEQS;   # the sequences to be cataloguged
my @QUERY_PEPTIDES; # a list of peptides to query
my $CATALOG_INFO;   # this will become a hashref of metadata about the catalog upon build/restore/get_catalog_info
my $PEPTIDE_LENGTH;
my $MAX_MISMATCHES;
my %UNIUQE_NMER_ID; # nmer sequence => ID
my %NMER_CATALOG;   # nmer_ID => sequence ID => position

# given:
# a list of query peptides
# a catalog name
# max_mismatches
# do:
# identify the query peptides that have X or less mismatches in the catalog
# return:
# hashref indexed by query peptides, contianing arrays of hashrefs indicating
# the matches
# e.g,
# GHSDHFAFH => [{ nmer_id => 454, num_mm => 1}, {nmer_id => 652, num_mm => 0}]  
#{ }
#
}
sub query_vs_catalog {

	my $max_mismatches = shift;

	if (defined $MAX_MISMATCHES) {
		die "Max mismatches already set to $MAX_MISMATCHES; Cannot set to $max_mismatches";
	}
	$MAX_MISMATCHES = $max_mismatches;

	# read in the catalog; setting the %UNIUQE_NMER_ID and %NMER_CATALOG hashes
	retrieve_catalog();

	# update the unique_nmer_id hash to include any new nmers found in the query peptides
	update_unique_nmers(\@QUERY_PEPTIDES);

	# determine the offsets, given the number of mismatches
	my @offsets = determine_offsets();

	# determine the comparisons that need to be made by breaking the
	# query & catalog peptides into non-overlapping segments
	# the number of segments is num_mm + 1
	my %query_segments;
	my %catalog_segments;




}

# given a list of peptides, assign IDs to them
sub update_unique_nmers {

	my $peptide_list = shift;

	# get the maximum nmer id
	my $nmer_id = pop (sort {$UNIUQE_NMER_ID{$a} <=> $UNIUQE_NMER_ID{$b}});

	foreach my $p (@$peptide_list) {
		if (!defined $UNIUQE_NMER_ID{$p}) {
			$UNIUQE_NMER_ID{$p} = ++$nmer_id;
		}
	}

}

sub determine_comparisons {
        
    my $num_comparisons=0;
 
    my $t0 = Benchmark->new();
    
    my $offset_start = 0;
    
    foreach my $offset (@OFFSETS) {

    	$offset_start = $offset->{start};
    	my $word_length = $offset->{length};

        print "\tworking on offset ", $offset_start, " length: ", $word_length, "\n";
        my $num_checked =0;
        
        my $t1 = Benchmark->new;
        
        my %query_index = create_hash($word_length, $offset_start, \@QUERY_PEPTIDES);
        my %catalog_index = create_hash($word_length, $offset_start, \@peptides2);
        
        my $to_check = (keys %index1);
        print "Nmers to check for offset $offset_start: $to_check\n";
        my %tb_myco;
        foreach my $sub_peptide1 (keys %index1) {

            $num_checked++;
            if ($num_checked % 1000 == 0) {
                print "Checking $num_checked of $to_check\n";
            }

            # moving on if there are no myco peptides that also begin with this sub peptide
            # skip if complete

            $num_checked++;
            if ($num_checked % 1000 == 0) {
                print "Checking $num_checked of $to_check\n";
            }

            next unless ($index2{$sub_peptide1});
            
            # moving on if we've don this comparison already
            #my $done_key = $offset_start . $sub_peptide1;
            #next if ($done{$done_key});
            # loopeing through all peptides from list 1 that start ith this subpeptide
            foreach my $id1 (@{$index1{$sub_peptide1}}) {
                # this tb sub peptide must be compared to all myco peptides that also contain this sub ppeptide
                do_comparisons($id1, $index2{$sub_peptide1},$offset);
                $num_comparisons+= @{$index2{$sub_peptide1}};
                #do the comparison and print output...we can uniquify later

            }
            # mark this sub peptide comparison as complte
            #$done{$done_key} = 1;

        }
        
        my $t2 = Benchmark->new;
        
        my $td = timediff($t2, $t1);
        my $tds = timediff($t2, $t0);
        
        print "Time for loop: ", timestr($td), "\n";
        print "Time since start: ", timestr($tds), "\n";

        print "Total comparisons after offset $offset_start: ", format_number($num_comparisons,0,0), "\n";

    
    }
    
    return $num_comparisons;
    
}


# given:
# the hashref returned from query_vs_catalong
# the protein_peptides hashrefs
# a filename
# output:
# a TSV with the columns:
# query_peptide	matching_peptide num_mm matching_proteins matching_positions
sub output_matching_peptides {


}

# given the name of a query file containing a list (or fasta)
# of peptides, return an array reference to the peptides
sub read_query_file {

	my $input_file = shift;

	open my $INF, '<', $input_file or die "Can't open input file ($input_file) for reading: $!";
	while (<$INF>) {
		chomp;
		push @QUERY_PEPTIDES, $_;
		my $cur_peptide_length = length($QUERY_PEPTIDES[-1]);
		# set the peptide length to the length of the first peptide_length
		# all subsequent peptides are checked that they match in length
		if (!defined $PEPTIDE_LENGTH) {
			$PEPTIDE_LENGTH = $cur_peptide_length;
		}
		if ($cur_peptide_length ne $PEPTIDE_LENGTH) {
			die "Peptides of different lengths ($peptide_length and $cur_peptide_length) found in query file: $input_file";
		}
	}
	close $INF;

	if ($PEPTIDE_LENGTH ne $catalog_info->{peptide_length}) {
		die "Peptide length ($peptide_length) differs from catalog peptide length (" . $catalog_info->{peptide_length} . ')';
	}

	return \@QUERY_PEPTIDES;
}

# given the name of a fasta file
# import all sequences into the global @CATALOG_SEQS
# return a reference to this array hashrefs of sequences {id => 1, seq=>DFDSFSDF}
sub read_fasta {

	my $fasta = shift;

	# setting the catalog source globally
	# die if it's already set
	if (defined $CATALOG_SOURCE) {
		die "Catalog source already set to $CATALOG_SOURCE; Cannot set to $fasta";
	}
	$CATALOG_SOURCE = $fasta;

	my $seqin = Bio::SeqIO->new(-file => $fasta, '-format' => 'Fasta');

	my %seen_id;
	my $seq_num = 1;
	while (my $seq = $seqin->next_seq()) {
		my $seq_id = $seq->display_id();
		if ($seen_id{$seq_id}) {
			warn "Multiple sequences with same identifier: $seq_id";
			warn "Appending number to id";
			$seq_id .= ".$seq_num";
		}
		push @CATALOG_SEQS, { id => $seq_id, seq => $seq->seq() };
		$seen_id{$seq_id} = 1;
		$seq_num++;
	}

	return \@CATALOG_SEQS;

}

# given:
# 1. an array of hashrefs of protein sequences as returned by read_fasta
# 2. a peptide length
# do:
# find all unique nmer sequences
# catalog the proteins in the input file, recording the unique nmers and their positions
# serialize the catalog according to the selected backend
# return:
# a reference to the unqiue nmer hash (%UNIQUE_NMER_ID)
# a reference to the protein catalog (maybe this should be in an actual DB?) (%NMER_CATALOG)
sub build_catalog {

	my $peptide_length = shift;
	my $catalog_name = shift;

	# set the catalog name & peptide length globally
	# die if their already set
	if (defined $PEPTIDE_LENGTH) {
		die "Peptide length already set to $PEPTIDE_LENGTH; Cannot set to $peptide_length";
	}
	$PEPTIDE_LENGTH = $peptide_length;

	if (defined $CATALOG_NAME) {
		die "Catalog name already set to $CATALOG_NAME; Cannot set to $catalog_name";
	}
	$CATALOG_NAME = $catalog_name;

	# the first nmer ID will be 0 so we can store them efficiently in an array
	my $nmer_id = 0;

	my @nmer_ids;
    
	foreach my $p (@CATALOG_SEQS) {

		my $protein_nmers_ref = break_protein($p->{seq}, $peptide_length);

		foreach my $start (0 .. $#$protein_nmers_ref) {
	    
			my $s = $protein_nmers_ref->[$start];

			# TODO: Here is where we branch depending on the catlog backend
			# if we find the current nmer in the unique nmers, use the ID, otherwise create a new id
			if (!exists $UNIUQE_NMER_ID{$s}) {
				$UNIUQE_NMER_ID{$s} = $nmer_id++;
			}
			my $cur_nmer_id = $UNIUQE_NMER_ID{$s};

		  	# add the position of the nmer onto an array
		  	push @{$NMER_CATALOG{$cur_nmer_id}{$p->{id}}}, $start;
	    }
	}

	# now serialize the uniqe nmers & nmer_catalog, and include a description
	# TODO: this should be refactored and dependent upon the backend
	my ($unique_nmers_file, $protein_peptides_file, $catalog_info_file) = get_catalog_file_names($CATALOG_NAME);

	if (!-d $CATALOG_NAME) {
		print "Creating directory for catalog: $CATALOG_NAME\n";
		mkdir $CATALOG_NAME or die "Failed creating directory for catalog: $!";
	}

	print "Serializing uniqe nmers to: $unique_nmers_file\n";
	nstore \%UNIUQE_NMER_ID, $unique_nmers_file;
	print "Done\n";
	print "Serializing protein peptides to: $protein_peptides_file\n";
	nstore \%NMER_CATALOG, $protein_peptides_file;
	print "Done\n";
	print "Storing catalog info in: $catalog_info_file\n";
	$catalog_info =    { catalog_name   => $CATALOG_NAME,
		                 data_source    => $CATALOG_SOURCE,
	                     peptide_length => $PEPTIDE_LENGTH,
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

	return (\%UNIUQE_NMER_ID, \%NMER_CATALOG);

}

# read in the catalog info file & return a hashref
# set potentially set the CATALOG_NAME globally here
sub get_catalog_info {

	my $catalog_name = shift;

	if (defined $CATALOG_NAME) {
		die "Catalog name already set to $CATALOG_NAME; Cannot set to $catalog_name";
	}
	$CATALOG_NAME = $catalog_name;

	my ($unique_nmers_file, $protein_peptides_file, $catalog_info_file) = get_catalog_file_names($CATALOG_NAME);

	my $catalog_info_string;
	print "Reading catalog info from: $catalog_info_file\n";
	open my $CAT, '<', $catalog_info_file or die "Can't open catalog info file ($catalog_info_file): $!";
	while (<$CAT>) {
		$catalog_info_string .= $_;
	}
	close $CAT;
	$CATALOG_INFO = JSON::XS->new->utf8->decode($catalog_info_string);
	print "Catalog info\n";
	print Dumper $CATALOG_INFO;

	return $CATALOG_INFO;

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

	# TODO: the functionlity will change depending on the backend
	my ($unique_nmers_file, $protein_peptides_file) = get_catalog_file_names($CATALOG_NAME);

	print "Restoring unique nmer reference from: $unique_nmers_file\n";
	my $unique_nmers_ref = retrieve($unique_nmers_file);
	print "Done\n";

	#TODO: this can be done immediately befor the lookups are done for the output file
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

# this will determine the offsets based on the length of the nmers
# and the maximium number of mismathces
# TODO: enforce a minimum length for the blocks - we don't want to end up
# witih the terminal block being only 1 or 2 residues
sub determine_offsets {
	
	my $num_offsets = $MAX_MISMATCHES + 1;
	
	my $avg_block_length = $PEPTIDE_LENGTH / $num_offsets;
	my $large_block_length = ceil($avg_block_length);
	my $small_block_length = floor($avg_block_length);
	
	my @offsets;
	
	my $num_large_blocks = $num_offsets;
	my $num_small_blocks = 0;
	
	if ($large_block_length != $small_block_length) {
		
		$num_large_blocks = ($avg_block_length - $small_block_length) * $num_offsets;
		$num_small_blocks = $num_offsets - $num_large_blocks;
		
	}
	
	my $offset = 0;
	foreach my $i (1..$num_large_blocks) {
		push @offsets, {
			start => $offset,
			length => $large_block_length
		};
		$offset += $large_block_length;
	}
	
	foreach my $i (1..$num_small_blocks) {
		push @offsets, {
			start => $offset,
			length => $small_block_length
		};
		$offset += $small_block_length;
	}
	
	return @offsets;
	
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