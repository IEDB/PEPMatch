#!perl -T
use 5.006;
use strict;
use warnings;
use Findbin::Bin;
use lib "$Findbin::Bin/../lib";
use Test::More;

plan tests => 1;

BEGIN {
    use_ok( 'NmerMatch' ) || print "Bail out!\n";
}

diag( "Testing NmerMatch $NmerMatch::VERSION, Perl $], $^X" );