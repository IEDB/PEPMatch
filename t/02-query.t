use 5.006;
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib";
use NmerMatch qw(get_catalog_info read_query_file query_vs_catalog output_matching_peptides);
use Test::More;
use Data::Dumper;
use Data::Compare;

plan tests => 2;

my $expected_results = get_expected_results();

my $catalog_name = "$FindBin::Bin/test_cat1";
my $query_input = "$FindBin::Bin/test_mm.lst";
my $max_mismatches = 1;

my $catalog_info = get_catalog_info($catalog_name);

my $query_peptides_ref = read_query_file($query_input);

ok(Compare($query_peptides_ref, $expected_results->{query_peptides_ref}), 'query peptide import');

my $num_mismatches_ref = query_vs_catalog($max_mismatches);

ok(Compare($num_mismatches_ref, $expected_results->{num_mismatches_ref}), 'results as expected');

# load the expected results to compare with what is generated
sub get_expected_results {
 
	my $query_pepties = ['MTMDKSELVQKAKLA','XTMDKSELVQKAKLA','MTMDKXELVQKAKLA','MTMDKSELVQKAKLX','GKEYREKIEAELQDI','XKEYREKXEAELQDI','XKEYREKIEAELQDX','GKEYREKIXAELXDI','MDDREDLVYQAKLAE','EQNKEALQDVEDENQ','EQNKEALQDVEDENX','EQNKXALQDVEDENQ','EQNKSDFASDADDSD'];
	my $num_mm_ref = {'0' => {'0' => 0},'456' => {'0' => 1},'460' => {'453' => 1},'232' => {'232' => 0},'454' => {'0' => 1},'455' => {'0' => 1},'453' => {'453' => 0},'461' => {'453' => 1},'80' => {'80' => 0}};

	my $expected_results = {
		query_peptides_ref => $query_pepties,
		num_mismatches_ref => $num_mm_ref,
	};

	return $expected_results;

};
