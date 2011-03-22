use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'MarkerDivergence',
    license             => 'perl',
    dist_author         => q{Jeffrey Barrick <jbarrick@msu.edu>},
    dist_version_from   => 'lib/MarkerDivergence.pm',
    build_requires => {
        'Test::More' => 0,
    },
    add_to_cleanup      => [ 'MarkerDivergence-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();
