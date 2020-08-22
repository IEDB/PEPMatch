use 5.006;
use strict;
use warnings;
use Module::Build;
Module::Build->VERSION('0.4004');

my $builder = Module::Build->new(
    module_name         => 'NmerMatch',
    license             => 'GPL_3',
    dist_author         => q{Jason Greenbaum <jgbaum@gmail.com>},
    dist_version_from   => 'lib/NmerMatch.pm',
    dist_abstract       => 'NmerMatch library and dependencies',
    release_status      => 'stable',
    configure_requires => {
        'Module::Build' => '0.4004',
    },
    test_requires => {
        'Test::More' => '0',
    },
    requires => {
        'Data::Compare' => '0',
        'DBD::SQLite'   => '0',
        'DBI'           => '0',
        'FAST'          => '0',
        'Inline::C'     => '0',
        'JSON::XS'      => '0',
        'List::Util'    => '0',
        'Scalar::Util::Numeric' => '0',
        'Template::Toolkit' => '0',
    },
    add_to_cleanup     => [ 'NmerMatch-*' ],
);

$builder->create_build_script();
