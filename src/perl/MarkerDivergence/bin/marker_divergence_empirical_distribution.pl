#!/usr/bin/perl -w

###
# Pod Documentation
###

=head1 NAME

marker_divergence_significance.pl

=head1 SYNOPSIS

Usage: marker_divergence_significance.pl -i input.fit -d path_to_fits -o output.tab [-r -p output.pdf]

Perform KS test to compare alpha and tau distributions loaded from ".fit" files
in the indicated directory to empirical values loaded from the input file.

=head1 DESCRIPTION

=over

=item B<-i> <file path> 

Path to fit file containing experimental marker divergence data. REQUIRED

=item B<-f> <file path> 

Path to fit file containing simulated marker divergence data. May be provided multiple times
(for example if local average mode is being used) REQUIRED

=item B<-p> <file path> 

Path to output plot. REQUIRED


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

#Get options
use Getopt::Long;
use Pod::Usage;

use KS_test_2d;

my $ratio_mode = 0;
my $experimental_mode;

my ($help, $man);
my ($input_file, $plot_file, $matrix_file, $output_file, @fit_files, $recursive, $new);
our $detailed_colors = 0;
my $p_value_cutoff = 0.05;
my $local_average;
GetOptions(
	'help|?' => \$help, 'man' => \$man,
	'input|i=s' => \$input_file,
	'fit-files|f=s' => \@fit_files,
	'plot-file|p=s' => \$plot_file,
);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
pod2usage(1) if (!defined $input_file || scalar @fit_files == 0 || !defined $plot_file);
die "Input file does not exist" if (!-e $input_file);


###
# Create the R script file
###
print STDERR "Creating R script...\n";

open R_SCRIPT, ">$$.r_script";
print R_SCRIPT <<EOS;
experimental <- read.table("$input_file", header = TRUE)
theoretical <- c()
EOS

foreach my $fit_file (@fit_files)
{
	print R_SCRIPT <<EOS;
	theoretical_file <- read.table("$fit_file", header = TRUE)
	add <- cbind(theoretical_file\$alpha, theoretical_file\$tau)
	theoretical <- rbind(theoretical, add)
EOS

}


print R_SCRIPT <<EOS;

#	minx <- 0
#	maxx <- 1.5
	
#	miny <- 0
#	maxy <- 70

	maxx <- max(theoretical[,1], experimental\$alpha)
	minx <- min(theoretical[,1], experimental\$alpha)

	maxy <- max(theoretical[,2], experimental\$tau)
	miny <- min(theoretical[,2], experimental\$tau)
	
	pdf("$plot_file", height=6, width=6)
	
	
	plot(0, type="n", lty="solid", xlim=c(minx, maxx), ylim=c(miny, maxy), xaxs="i", yaxs="i", labels=T, lwd=1, axes=F, xlab="alpha", ylab="tau", las=1, cex.lab=1.2, cex.axis=1.2 )
	box()
	axis(2, cex.lab=1.2, las=1, cex.axis=1.2, yaxs="i", at = c(0:7)*10 )
	axis(1, cex.lab=1.2, las=1, cex.axis=1.2, xaxs="i", at = c(0:6)*0.25 )
	
	points(theoretical, pch=3, col=rgb(120,120,120, maxColorValue=255), cex=1)
	points(experimental\$alpha, experimental\$tau, pch=3, col="red", cex=1.2, lwd=2)
	dev.off()
EOS

close R_SCRIPT;

#Run R
print STDERR "Running R script: $$.r_script\n";
my $command = "R --vanilla < $$.r_script &> $$.r_output";
my $res = system $command;
die "Error running command: $command" if ($res);

#Read R output
print STDERR "Parsing R output: $$.r_output\n";
open R_OUTPUT, "$$.r_output";
my @lines = <R_OUTPUT>;
close R_OUTPUT;

unlink "$$.r_output";
unlink "$$.r_script";
