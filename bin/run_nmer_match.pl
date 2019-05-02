use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Getopt::Long qw(GetOptions);
use NmerMatch qw(read_fasta build_db);
use Data::Dumper;

# actions:
#  1. run a search
#  2. build a database
# 
# run a search
# input: 
#   database name
#   nmer fasta
#
# build a database
# input:
#   database name
#   db fasta
#   peptide length

my $action;              # build or search
my $nmer_length = 15;    # peptide length for the database
my $db_fasta;            # input fasta file to build the database
my $db_name;             # name of the database to build/use
                         # for now, this is a file, but it could be a reference
                         # to a schema in a database
my $nmer_fasta;          # fasta of peptides to search against the databse

GetOptions ("action|a=s"       => \$action,
	        "nmer-length|l=i"  => \$nmer_length,
            "db-fasta|d=s"     => \$db_fasta,
            "db-name|b=s"      => \$db_name,
            "nmer-fasta|n=s"   => \$nmer_fasta);

if ($action eq 'build') {

	# read in the db fasta file
	print "Reading fasta: $db_fasta\n";
	my $db_seq_ref = read_fasta($db_fasta);

	# catalog the input sequences by unique nmer
	print "Cataloging sequences\n";
	my ($unique_nmer_ref, $db_catalog_ref) = build_db($db_seq_ref, $nmer_length);

	print "Done cataloging - sleeping!\n";
	while (1) {

	};
	#print Dumper $unique_nmer_ref;
	#print Dumper $db_catalog_ref;

	# write the database

}
elsif ($action eq 'search') {

	# read in the nmer fasta file

	# read the database

	# compare & output

}
else {
	die "Unknown action: $action";
}
