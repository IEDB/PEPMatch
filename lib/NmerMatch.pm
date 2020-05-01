package NmerMatch;
require Exporter;
@ISA = qw(Exporter);
@EXPORT_OK = qw(read_fasta break_protein build_catalog retrieve_catalog get_catalog_info read_query_file query_vs_catalog output_matching_peptides);

use strict;
use warnings;
use Data::Dumper;
use FAST::Bio::SeqIO;
use FindBin;
use Storable qw(nstore retrieve);
use JSON::XS;
use POSIX qw(strftime);
use POSIX qw(floor ceil);
use Benchmark;
use DBI;
use Template;
use List::Util qw(min max);
use Inline 'C';
use Inline (C => Config =>
            OPTIMIZE => '-O3',
            #CCFLAGSEX => '-march=x86-64',
            );

# Workflow for searching
# 1. get catalog info to set the CATALOG_NAME
# 2. read in the query sequences
# 3. read in the catalog (UNIQUE_NMER_ID)
# ??4. add unique nmers from query seqs to UNIQUE_NMER_ID - is this necessary
# 5. create list of catalog nmers
# 6. create list of query nmers
# 7. determine which ones need to be compared & run the comparisons
# ??8. encode each nmer as an integer ID (UNIQUE_NMER) pointing to the sequence



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
my $MAX_MISMATCHES;
my $UNIQUE_NMER_ID; # nmer sequence => nmer ID
my $NMER_CATALOG;   # [nmer_ID] => {sequence ID} => [positions]
my $SEQ_NAMES_CATALOG; # [sequence ID] => sequence name
my %NUM_MISMATCHES;  # {query nmer ID}{catalog nmer ID} => number of mismatches
my %NMER_ID2SEQ;     # captures the nmer IDs that will need to be reported later
                     # and stores their sequence

my $t0 = Benchmark->new;


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

	# TODO: Track the comparisons that were done so we don't need to redo

	my $num_comparisons_done = 0;
	foreach my $o (@offsets) {

		my $t1 = Benchmark->new;

		print "Working on offset:\n";
		print Dumper $o;

		my ($query_word_hash_ref, $query_nmers_ref) = create_query_word_hash($o->{length}, $o->{start});
		my ($catalog_word_hash_ref, $catalog_nmers_ref) = create_catalog_word_hash($o->{length}, $o->{start}, $query_word_hash_ref);

		# now make the comparisons
		foreach my $sub_peptide (keys %$query_word_hash_ref) {
			#print "Sub-peptide: $sub_peptide\n";
			foreach my $query_nmer_id (@{$query_word_hash_ref->{$sub_peptide}}) {
				#print "Query:\n";
				#print Dumper @{$query_word_hash_ref->{$sub_peptide}};

				my $query_p = $query_nmers_ref->{$query_nmer_id};

				foreach my $catalog_nmer_id (@{$catalog_word_hash_ref->{$sub_peptide}}) {
					#print "Catalog:\n";
					#print Dumper @{$catalog_word_hash_ref->{$sub_peptide}};
					$num_comparisons_done++;
					#print "Comparison:\n";
					#print $query_nmers_ref->{$query_nmer_id}, "\n";
					#print $catalog_nmers_ref->{$catalog_nmer_id}, "\n"; 
					# if 0 mismatches, no need to compare
					my $num_mm;
					if ($query_nmer_id ne $catalog_nmer_id) {
						$num_mm = count_mismatches($query_p, 
							                       $catalog_nmers_ref->{$catalog_nmer_id},
							                       $MAX_MISMATCHES);
					} else {
						$num_mm = 0;
					}
					if ($num_mm <= $MAX_MISMATCHES) {
						$NUM_MISMATCHES{$query_nmer_id}{$catalog_nmer_id} = $num_mm;
						# now store the sequences by their ID in the NMER_ID2SEQ hash
						# so we can recall them later
						if (!defined $NMER_ID2SEQ{$query_nmer_id}) {
							$NMER_ID2SEQ{$query_nmer_id} = $query_p;
						}
						if (!defined $NMER_ID2SEQ{$catalog_nmer_id}) {
							$NMER_ID2SEQ{$catalog_nmer_id} = $catalog_nmers_ref->{$catalog_nmer_id};
						}
					} 
					#print "Mismatches: $NUM_MISMATCHES{$query_nmer_id}{$catalog_nmer_id}\n\n";
				}
			}
		}

        my $t2 = Benchmark->new;
        
        my $td = timediff($t2, $t1);
        my $tds = timediff($t2, $t0);
        
        print "Time for loop: ", timestr($td), "\n";
        print "Time since start: ", timestr($tds), "\n";

	}

	my $num_query_peptides = @QUERY_PEPTIDES;
	my $possible_comparisions = $CATALOG_INFO->{num_nmers} * $num_query_peptides;
	print "Possible comparisons: $possible_comparisions\n";
	print "Comparisons done: $num_comparisons_done\n";

	return \%NUM_MISMATCHES;

}

# given a word hash length, an offset, and a list of peptides
# return a hash keyed by substring, each value being an array of
# unique nmer ids that contain that substring at that postion
sub create_query_word_hash {
    
    my $length = shift;
    my $offset = shift;
    
    my %word_index;
    my %id_to_nmer;  # a map of nmer_id to nmer so we can look up later
    
    # pushing the IDs of the unique peptides starting with this substring onto
    # the index
    foreach my $query_p (@QUERY_PEPTIDES) {
    	# we only want to push the the peptide into the array once
    	next if (defined $id_to_nmer{$UNIQUE_NMER_ID->{$query_p}});
    	$id_to_nmer{$UNIQUE_NMER_ID->{$query_p}} = $query_p;
        push @{$word_index{substr $query_p, $offset, $length}}, $UNIQUE_NMER_ID->{$query_p};
    }
    
    return (\%word_index, \%id_to_nmer);
    
}

# similar to the query word hash, create a hash of the catalog words,
# skipping those that do not contain substrings observed in the query word hash
sub create_catalog_word_hash {

	my $length = shift;
	my $offset = shift;
	my $query_hash_ref = shift; # a reference to the query subpeptide hash

    my %word_index;
    my %id_to_nmer;
    
    my $catalog_max_nmer_id = $CATALOG_INFO->{max_nmer_id};

    # pushing the IDs of the unique peptides starting with this substring onto
    # the index
    foreach my $catalog_p (keys %$UNIQUE_NMER_ID) {
    	# skip if it's not part of the catalog nmers
    	next if $UNIQUE_NMER_ID->{$catalog_p} > $catalog_max_nmer_id;
    	my $sub_peptide = substr $catalog_p, $offset, $length;
    	# skip if the subpeptide is not included in the query
    	next if (!defined $query_hash_ref->{$sub_peptide});
    	# we only want to push the the peptide into the array once
    	next if (defined $id_to_nmer{$UNIQUE_NMER_ID->{$catalog_p}});
    	$id_to_nmer{$UNIQUE_NMER_ID->{$catalog_p}} = $catalog_p;
        push @{$word_index{$sub_peptide}}, $UNIQUE_NMER_ID->{$catalog_p};
    }
    
    return (\%word_index, \%id_to_nmer);

}

# given a list of peptides, assign IDs to them
sub update_unique_nmers {

	my $peptide_list = shift;

	# get the maximum nmer id
	my $nmer_id = $CATALOG_INFO->{max_nmer_id};

	foreach my $p (@$peptide_list) {
		if (!defined $UNIQUE_NMER_ID->{$p}) {
			$UNIQUE_NMER_ID->{$p} = ++$nmer_id;
		}
	}

}

# retrieve the nmer catalog from the protein_peptides sql db
# return only the peptides that were matched
sub retreive_nmer_catalog_sql {

	my ($unique_nmers_file, $protein_peptides_file,
		$protein_names_file, $catalog_info_file,
		$protein_peptides_csv) = get_catalog_file_names($CATALOG_NAME);

	print "Retrieving nmer catalog from: $protein_peptides_file\n";

	my $dbh = DBI->connect("dbi:SQLite:dbname=$protein_peptides_file","","");

	my %ids_to_retrieve;

	my $num_nmers_with_matches = (keys %NUM_MISMATCHES);
	my $num_retrieved = 0;

	# build a list of IDs that we need to retrieve
	foreach my $query_nmer_id (keys %NUM_MISMATCHES) {
		foreach my $match_nmer_id (keys %{$NUM_MISMATCHES{$query_nmer_id}}) {
			$ids_to_retrieve{$match_nmer_id} = 1;
		}
		#%ids_to_retrieve = (%ids_to_retrieve, map {$_ => 1} keys %{$NUM_MISMATCHES{$query_nmer_id}});
		if (++$num_retrieved % 1000 == 0) {
			print  "Retrieved $num_retrieved of $num_nmers_with_matches query peptides with matches\n";
		}
	}

	# let's split this up into batches of 1000 so as to limit the size of the queries
	my $batch_size = 1000;
	my @ids_to_retrieve = (keys %ids_to_retrieve);
	my $num_ids = @ids_to_retrieve;
	my $max_stop = $num_ids - 1;
	my $max_start = $max_stop - $batch_size;

	#print Dumper @ids_to_retrieve;
	print "IDs to retrieve: $num_ids\n";

	my $start = 0;
	my $stop = 0;

	while ($stop < $max_stop) {

		$stop = min($start+$batch_size-1, $max_stop);

		print "Fetching batch from $start to $stop\n";
		my @id_batch = @ids_to_retrieve[$start .. $stop];

		my $id_string = join ',', @id_batch;

		my $query = "select * from protein_peptide where nmer_id IN ($id_string)";

		#print "executing: $query\n";

		my $sth = $dbh->prepare($query);
		$sth->execute();
		while (my $row = $sth->fetch()) {
			push @{$NMER_CATALOG->{$row->[0]}{$row->[1]}}, $row->[2];
		}
		$sth->finish;

		$start += $batch_size;

	}
	$dbh->disconnect;

	print "Done retreiving nmer catalog!\n";

}

# given:
# the hashref returned from query_vs_catalong
# the protein_peptides hashrefs
# a filename
# output:
# a TSV with the columns:
# query_peptide	matching_peptide num_mm num_protein_matches matching_proteins_and_positions
sub output_matching_peptides {

	my $outfile = shift;
	my $long_outfile = shift;

	# if the NMER CATALOG is empty, retreive it from the db
	if (! keys %{$NMER_CATALOG}) {
		retreive_nmer_catalog_sql();
	}

	open my $OUTF, '>', $outfile or croak("Can't open ($outfile) for writing: $!");
	my @header_fields = qw(query_peptide	matching_peptide num_mm num_protein_matches matching_proteins_and_positions);
	print $OUTF join("\t", @header_fields), "\n";

	my $LNGF;
	if (defined $long_outfile) {
		open $LNGF, '>', $long_outfile or croak("Can't open ($long_outfile) for writing: $!");
		my @long_header_fields = qw(query_peptide	matching_peptide num_mm matching_protein pos_start);
		print $LNGF join("\t", @long_header_fields), "\n";
	}

	foreach my $query_nmer_id (keys %NUM_MISMATCHES) {
		my $query_p = $NMER_ID2SEQ{$query_nmer_id};
		foreach my $catalog_nmer_id (keys %{$NUM_MISMATCHES{$query_nmer_id}}) {
			my $catalog_p = $NMER_ID2SEQ{$catalog_nmer_id};
			my $catalog_matches = $NMER_CATALOG->{$catalog_nmer_id};
			my @matching_protein_ids = (keys %$catalog_matches);
			my $num_protein_matches = @matching_protein_ids;
			my @matching_proteins_and_positions;
			foreach my $seq_id (@matching_protein_ids) {
				foreach my $position (@{$catalog_matches->{$seq_id}}) {
					push @matching_proteins_and_positions, $SEQ_NAMES_CATALOG->[$seq_id] . ':' . $position;
				}
			}

			my $num_mm = $NUM_MISMATCHES{$query_nmer_id}{$catalog_nmer_id};
			my $matching_prots_and_pos_string = join '; ', @matching_proteins_and_positions;
			my @output_fields = ($query_p, $catalog_p, $num_mm, $num_protein_matches,
				$matching_prots_and_pos_string);

			print $OUTF join("\t", @output_fields), "\n";

			# print a row per match, if the long outfile is defined
			if (defined $long_outfile) {
				foreach my $seq_id (@matching_protein_ids) {
					foreach my $position (@{$catalog_matches->{$seq_id}}) {
						my @long_output_fields = ($query_p, $catalog_p, $num_mm, 
							$SEQ_NAMES_CATALOG->[$seq_id], $position);
						print $LNGF join("\t", @long_output_fields), "\n";
					}
				}
			}

		}
	}
	close $OUTF;
	if (defined $long_outfile) {
		close $LNGF;
	}
}

# given the name of a query file containing a list (or fasta)
# of peptides, return an array reference to the peptides
sub read_query_file {

	my $input_file = shift;
	my $expected_length = shift;

	open my $INF, '<', $input_file or die "Can't open input file ($input_file) for reading: $!";
	while (<$INF>) {
		chomp;
		push @QUERY_PEPTIDES, $_;
		my $cur_peptide_length = length($QUERY_PEPTIDES[-1]);
		if ($cur_peptide_length ne $expected_length) {
			die "Mismatch between catalog peptide length ($expected_length) and length of peptide ($cur_peptide_length - $_) on line $. of $input_file";
		}
	}
	close $INF;

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

	my $seqin = FAST::Bio::SeqIO->new(-file => $fasta, '-format' => 'Fasta');

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
# a reference to the unqiue nmer hash (%$UNIQUE_NMER_ID)
# a reference to the protein catalog (maybe this should be in an actual DB?) (%NMER_CATALOG)
sub build_catalog {

	my $peptide_length = shift;
	my $catalog_name = shift;

	# set the catalog name globally

	if (defined $CATALOG_NAME) {
		die "Catalog name already set to $CATALOG_NAME; Cannot set to $catalog_name";
	}
	$CATALOG_NAME = $catalog_name;

	if (!-d $CATALOG_NAME) {
		print "Creating directory for catalog: $CATALOG_NAME\n";
		mkdir $CATALOG_NAME or die "Failed creating directory for catalog: $!";
	}

	my ($unique_nmers_file, $protein_peptides_file,
		$protein_names_file, $catalog_info_file,
		$protein_peptides_csv) = get_catalog_file_names($CATALOG_NAME);

	# create a temporary CSV for the protein_peptide data, that will
	# eventually be loaded into a sqlite table
	open my $CSV, '>', $protein_peptides_csv or die "Can't open CSV for writing: $!";
	print $CSV join(',', qw(nmer_id protein_id position)), "\n";

	# the first nmer ID will be 0 so we can store them efficiently in an array
	my $nmer_id = 0;
	my $seq_id = 0;
    
	foreach my $p (@CATALOG_SEQS) {

		$SEQ_NAMES_CATALOG->[$seq_id] = $p->{id};

		my $protein_nmers_ref = break_protein($p->{seq}, $peptide_length);

		foreach my $start (0 .. $#$protein_nmers_ref) {
	    
			my $s = $protein_nmers_ref->[$start];

			# TODO: Here is where we branch depending on the catlog backend
			# if we find the current nmer in the unique nmers, use the ID, otherwise create a new id
			if (!exists $UNIQUE_NMER_ID->{$s}) {
				$UNIQUE_NMER_ID->{$s} = $nmer_id++;
			}
			my $cur_nmer_id = $UNIQUE_NMER_ID->{$s};

		  	# add the position of the nmer onto an array
		  	#push @{$NMER_CATALOG->[$cur_nmer_id]{$seq_id}}, $start;
		  	print $CSV join(',', ($cur_nmer_id, $seq_id, $start)), "\n";
	    }

	    $seq_id++;
	    if ($seq_id % 1000 == 0) {
	    	print "Catalogued $seq_id sequences\n";
	    }
	}
	close $CSV;

	# now serialize the uniqe nmers & nmer_catalog, and include a description
	# TODO: this should be refactored and dependent upon the backend

	print "Serializing unique nmers to: $unique_nmers_file\n";
	nstore \%$UNIQUE_NMER_ID, $unique_nmers_file;
	print "Done\n";

	print "Serializing protein peptides to: $protein_peptides_file\n";
    my $tt = Template->new( { ABSOLUTE => 1, } );
    my $sqlite_template = "$FindBin::Bin/../lib/sqlite.template.sql";
    my $db_build_file = "$CATALOG_NAME/build_db.sql";
    my $data = {protein_peptides_csv => $protein_peptides_csv};
    $tt->process( $sqlite_template,
                  $data,
                  $db_build_file)
        || die $tt->error();
    # now we need to send the commands to sqlite
    my $sqlite_CMD = qq(sqlite3 $protein_peptides_file < $db_build_file);
    print "Executing: $sqlite_CMD\n";
    system($sqlite_CMD) == 0 ||
    die "Error building sqlite database: $?";
	print "Done\n";

	print "Serializing protein names to: $protein_names_file\n";
	nstore \@$SEQ_NAMES_CATALOG, $protein_names_file;
	print "Done\n";

	print "Storing catalog info in: $catalog_info_file\n";
	$CATALOG_INFO =    { catalog_name   => $CATALOG_NAME,
		                 data_source    => $CATALOG_SOURCE,
	                     peptide_length => $peptide_length,
                         num_nmers      => $nmer_id,
                         max_nmer_id    => $nmer_id - 1,
	                     build_date     => strftime "%F %R", localtime,
	                   };
	my $catalog_json = JSON::XS->new
	                           ->utf8
	                           ->pretty(1)
	                           ->canonical(1)
	                           ->encode($CATALOG_INFO);
	open my $CAT, '>', $catalog_info_file or die "Can't open catalog info file ($catalog_info_file): $!";
	print $CAT $catalog_json, "\n";
	close $CAT;

	return($UNIQUE_NMER_ID, $NMER_CATALOG, $SEQ_NAMES_CATALOG, $CATALOG_INFO)

}

# read in the catalog info file & return a hashref
# set potentially set the CATALOG_NAME globally here
sub get_catalog_info {

	my $catalog_name = shift;

	if (defined $CATALOG_NAME) {
		die "Catalog name already set to $CATALOG_NAME; Cannot set to $catalog_name";
	}
	$CATALOG_NAME = $catalog_name;

	my ($unique_nmers_file, $protein_peptides_file, $protein_names_file, $catalog_info_file) = get_catalog_file_names($CATALOG_NAME);

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
# TODO: fix this so that it returns a hashref
sub get_catalog_file_names {

	my $catalog_name = shift;

	my $unique_nmers_file = $catalog_name . '/unique_nmers.store';
	my $protein_peptides_file = $catalog_name . '/protein_peptides.db';
	my $protein_names_file = $catalog_name . '/protein_names.store';
	my $catalog_info_file = $catalog_name . '/catalog_info.json';
	my $protein_peptides_csv = $catalog_name . '/protein_peptides.csv';


	return ($unique_nmers_file, 
		$protein_peptides_file, 
		$protein_names_file, 
		$catalog_info_file,
		$protein_peptides_csv);
}


# given:
# a catalog name
# do:
# retrieve the catalog and return references to the uniqe nmer ids and nmer_catalog
# depending upon the backend
sub retrieve_catalog {

	# TODO: the functionlity will change depending on the backend
	my ($unique_nmers_file, $protein_peptides_file, $protein_names_file) = get_catalog_file_names($CATALOG_NAME);

	print "Restoring unique nmer reference from: $unique_nmers_file\n";
	$UNIQUE_NMER_ID = retrieve($unique_nmers_file);
	print "Done\n";

	print "Restoring protein names reference from: $protein_names_file\n";
	$SEQ_NAMES_CATALOG = retrieve($protein_names_file);
	print "Done\n";

	return($UNIQUE_NMER_ID, $NMER_CATALOG, $SEQ_NAMES_CATALOG);

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
	
	my $peptide_length = $CATALOG_INFO->{peptide_length};

	my $avg_block_length = $peptide_length / $num_offsets;
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
	# NOTE, we should be able to do 1..$num_large_blocks, but
	# this fails in certain situations, so we must do
	# 0..$num_large_blocks-1
	foreach my $i (0..$num_large_blocks-1) {
		push @offsets, {
			start => $offset,
			length => $large_block_length
		};
		$offset += $large_block_length;
	}

	foreach my $i (0..$num_small_blocks-1) {
		push @offsets, {
			start => $offset,
			length => $small_block_length
		};
		$offset += $small_block_length;
	}


	# validate that each position of the sequence is covered once
	my %position_covered = map {$_ => 0} (0..$peptide_length-1);
	foreach my $o (@offsets) {
		foreach my $i ($o->{start}..$o->{start} + $o->{length}-1) {
			if ($position_covered{$i}) {
				warn "Error in offset calculations as position $i is covered by multiple offsets\n";
				print "Offsets:\n";
				print Dumper @offsets;
				die();
			}
			$position_covered{$i} = 1;
		}
	}

	foreach my $i (0..$peptide_length-1) {
		if (!defined $position_covered{$i}) {
			warn "Error in offset calculations as position $i is not covered by any offset\n";
			print "Offsets:\n";
			print Dumper @offsets;
		    die();
		}
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