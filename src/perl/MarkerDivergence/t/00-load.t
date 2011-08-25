#!perl -T

use Test::More tests => 1;

BEGIN {
    use_ok( 'MarkerDivergence' ) || print "Bail out!
";
}

diag( "Testing MarkerDivergence $MarkerDivergence::VERSION, Perl $], $^X" );
