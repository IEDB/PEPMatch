use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../local_lib/lib/perl5"; # install dependencies here
use POSIX qw(floor);
use Getopt::Long qw(GetOptions);
use NmerMatch qw(read_fasta build_catalog get_catalog_info
	             read_query_file query_vs_catalog output_matching_peptides
	             remove_matched_query_peptides retrieve_catalog filter_matches);
use Data::Dumper;
use Scalar::Util::Numeric qw(isint);

# actions:
#  1. run a search
#  2. build a catalog
# 
# run a search
# input: 
#   catalog name
#   nmer list
#
# build a catalog
# input:
#   catalog name
#   catalog fasta
#   peptide length

my $action;              # build or search
my $nmer_length;         # peptide length for the catalog
my $catalog_fasta;       # input fasta file to build the catalog
my $catalog_name;        # name of the catalog to build/use
                         # for now, this is a file, but it could be a reference
                         # to a schema in a database
my $query_input;         # fasta/list of peptides to search against the databse
my $max_mismatches;      # the maximum mismatches for a search
my $mm_start;            # the number of mismatches to start with for the search-deep option
my $output_file;         # short version of output file with 1 row summarizing all
                         # matches for a given peptide
my $long_output_file;    # a longer version of the output, with 1 row per match
my $force = 0;           # if --force specified, ignore warnings for large values of max-mismatches
my $help = 0;

GetOptions ("action|a=s"         => \$action,
	        "nmer-length|l=i"    => \$nmer_length,
            "catalog-fasta|s=s"  => \$catalog_fasta,
            "catalog-name|c=s"   => \$catalog_name,
            "nmer-query|q=s"     => \$query_input,
            "max-mismatches|m=i" => \$max_mismatches,
            "mm-start=i"         => \$mm_start,
            "output-file|o=s"    => \$output_file,
            "long-output-file|v=s" => \$long_output_file,
            "force"              => $force,
            "help|h"             => \$help,
            );

my $USAGE = qq(
	perl run_nmer_match.pl -a [build/search] -c [catalog name] [OPTIONS]
	
	  --action|-a             'build', 'search', 'search-deep' - required

	  --catalog-name|-c       name of the catalog for building/searching - required

	  --catalog-fasta|-s      fasta file for cataloging - required for 'build'

	  --nmer-length|l         length of nmers to catalog - required for 'build'

	  --nmer-query|-q         file containing list of nmers (one per line) for searching
	                          against the catalog - required for 'search/search-deep'

	  --max-mismatches|-m     maximum number of mismatches for the query sequences
	                          against the catalog - required for 'search'

	  --mm-start              for the 'search-deep' action, the number of mismatches to begin
	                          with (default = floor(nmer_length/2))

	  --output-file|-o        output TSV file containing one row per unique match and
	                          with a consolidated list of catalog sequences and postions
	                          containing the match - required for 'search/search-deep'

	  --long-output-file|-v   output TSV file with one row per query sequence and
	                          matching catalog sequence. This is the same content as
	                          output-file, but may be easier to parse

	  --force                 allow very large values of max-mismatches

	  --help                  print this help
);

if ($help) {
	print $USAGE, "\n";
	exit;
}

# check the parameters
# action

if (!defined $action) {
	die "'action' must be specified\n\n$USAGE";
}

my %valid_actions = map { $_ => 1} qw(build search search-deep);
if (!defined $valid_actions{$action}) {
	die "\nSpecified action ($action) is invalid\n\n$USAGE\n";
}

# catalog name
if (!defined $catalog_name) {
	die "'catalog-name' must be specified\n\n$USAGE";
}

if ($action eq 'build') {

	# check the parameters specific for build
	# catalog fasta
	if (!defined $catalog_fasta) {
		die "'catalog-fasta' must be specified when building a catalog\n\n$USAGE";
	}

	# nmer length
	if (!defined $nmer_length) {
		die "'nmer-length' must be specified when building a catalog\n\n$USAGE";
	}

	if (!isint $nmer_length) {
		die "'nmer-length' ($nmer_length) must be an integer\n\n$USAGE"
	}

	if ($nmer_length >= 50) {
		warn "Long 'nmer-length' of $nmer_length may result in increased memory usage and poor performance.";
	}

	# read in the db fasta file
	print "Reading fasta: $catalog_fasta\n";
	my $catalog_seq_ref = read_fasta($catalog_fasta);

	# catalog the input sequences by unique nmer
	print "Cataloging sequences\n";
	build_catalog($nmer_length, $catalog_name);

}
elsif (($action eq 'search') or ($action eq 'search-deep')) {

	# check the parameters specific for search
	# query_input
	if (!defined $query_input) {
		die "'nmer-query' must be specified when searching against a catalog\n\n$USAGE";
	}

	# max_mismatches
	if ($action eq 'search' && !defined $max_mismatches) {
		die "'max-mismatches' must be specified when searching against a catalog\n\n$USAGE";
	}

	if ($action eq 'search' && !isint $max_mismatches) {
		die "'max-mismatches' ($max_mismatches) must be an integer\n\n$USAGE"
	}

	# output file
	if (!defined $output_file) {
		die "'output-file' must be specified when searching against a catalog\n\n$USAGE";
	}

	# read in the catalog info only
	# NOTE: this sets the 'CATALOG_NAME' global in the package, so no need to pass this
	# later on
	my $catalog_info = get_catalog_info($catalog_name);

	# read in the nmer query file (list format)
	# NOTE: this sets the 'QUERY_PEPTIDES' in the package, so no need to pass it
	#       to the next function call
	my $catalog_peptide_length = $catalog_info->{peptide_length};

	# mm-start (search-deep only)
	if (!defined $mm_start) {
		$mm_start = floor($catalog_peptide_length/2);
	}


	if ($force) {
		if ($max_mismatches >= $catalog_peptide_length) {
			die "The max-mismatches of $max_mismatches is greater than or equal to the peptide length. Choose a value lower than the peptide length."
		}
	}
	else {
		if ($action eq 'search') {
			# if max mimsmatches is 90% or more of the peptide length, quit
			if ($max_mismatches >= 0.9 * $catalog_peptide_length) {
				die "The max-mismatches of $max_mismatches is too high. Choose a value lower than 90% of the peptide length."
			}

			# if max mismatches is more than 60% of the peptide length, throw a warning
			if ($max_mismatches >= 0.6 * $catalog_peptide_length) {
				print "The max-mismatches of $max_mismatches is very high. Note that performance may suffer at high values of max-mismatches relative to peptide length.  Are you sure you want to continue? (Y/N): ";
			    my $reply = <STDIN>;
			    chomp $reply;
			    if (uc($reply) ne 'Y') {
			    	exit;
			    }
			}
		}
	}

	my $query_peptides_ref = read_query_file($query_input, $catalog_peptide_length);

	my ($uid, $nmer_cat, $seq_names_cat) = retrieve_catalog();

	if ($action eq 'search') {

		# compare & output
		query_vs_catalog($max_mismatches);


	} elsif ($action eq 'search-deep') {

		$max_mismatches = $mm_start;
		while (@$query_peptides_ref) {
			my $num_query_peptides = scalar @$query_peptides_ref;
			print "Searching $num_query_peptides remaining query peptides for max mismatches: $max_mismatches\n";
			query_vs_catalog($max_mismatches);
			$query_peptides_ref = remove_matched_query_peptides();
			$max_mismatches++;
		}

		# now filter the matching peptide lists for each query peptide
		# to only include those matches with the lowest number of mm
		print "Filtering matches to include only the best hits for each peptide\n";
		filter_matches();

	}

	print "Writing output files\n";
	output_matching_peptides($output_file, $long_output_file);

}
else {
	die "Unknown action: $action";
}
