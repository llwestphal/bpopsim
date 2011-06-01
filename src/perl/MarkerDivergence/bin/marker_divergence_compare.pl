#!/usr/bin/perl -w

###
# Pod Documentation
###

=head1 NAME

marker_divergence_compare.pl

=head1 SYNOPSIS

Usage: marker_divergence_compare.pl -1 file1.fit -2 file2.fit

Compares the alpha and tau values in two different fit files using the 2 dimensional KS test.

=head1 DESCRIPTION

=over

=item B<-1> <file path> 

First input fit file.

=item B<-o> <file path> 

Second input fit file.

=back

=head1 AUTHOR

Jeffrey Barrick <jbarrick@msu.edu>

=head1 COPYRIGHT

Copyright 2008-2009.  All rights reserved.

=cut

###
# End Pod Documentation
###


use strict;

use FindBin;
use lib $FindBin::Bin;
use Data::Dumper;

use KS_test_2d;

#Get options
use Getopt::Long;
use Pod::Usage;

my ($file_name_1, $file_name_2);

my ($help, $man);
GetOptions(
	'help|?' => \$help, 'man' => \$man,
	'file1|1=s' => \$file_name_1,
	'file2|2=s' => \$file_name_2,
);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
pod2usage(1) if (!defined $file_name_1 || !defined $file_name_2);


my $list_ref_1 = load_fit_file($file_name_1);
my $list_ref_2 = load_fit_file($file_name_2);

my $pr = bf_KS_test_2d($list_ref_1, $list_ref_2);
print "p-value that distributions are the same = $pr\n";


sub load_fit_file
{
	my ($file_name) = @_;
	
	open FIT, "<$file_name";
	my @lines = <FIT>;
	
	@lines = grep !/^\s*#/, @lines;
	@lines = grep !/^\s*$/, @lines;
	shift @lines;

	my @pts;
	foreach my $line (@lines)
	{
		my @split_line = split "\t", $line;
		my $pt;
		($pt->{x}, $pt->{y}) = ($split_line[4], $split_line[5]);
		push @pts, $pt;
	}
	
	return \@pts;
}

