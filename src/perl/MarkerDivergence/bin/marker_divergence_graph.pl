#!/usr/bin/perl -w

###
# Pod Documentation
###

=head1 NAME

marker_divergence_graph.pl

=head1 SYNOPSIS

Usage: marker_divergence_graph.pl -i file.tab -f file.fit -o file.pdf

Use R to graph and fit curves.

=head1 DESCRIPTION

=over

=item B<-i> <file path> 

Input raw data file.

=item B<-m|--data-mode> <mode> 

Format of data. Either 'ratio', 'log_ratio', 'log10_ratio', or 'fraction'. Defaults to 'ratio'.

=item B<-f> <format> 

Fit file.

=item B<-o> <file path> 

Output graph file.

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

my ($help, $man, $verbose);

my $input_file;
my $output_file;
my $fit_file;
my $distribution_file;
my $data_mode = 'ratio';
my $transfer_multiplier = 8;
my $max_generations = 0;

GetOptions(
	'help|?' => \$help, 'man' => \$man,
	'input|i=s' => \$input_file,
	'output|o=s' => \$output_file,
	'fit|f=s' => \$fit_file,
	'distribution|d=s' => \$distribution_file,
	'data-mode|m=s' => \$data_mode,
	'max-generations|g=s' => \$max_generations,
	'transfer-multiplier|t=s' => \$transfer_multiplier,
);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
pod2usage(1) if (!defined $input_file || !defined $output_file);

print "$$\n";

## create R script
open R_SCRIPT, ">/tmp/$$.r_script";
print R_SCRIPT <<EOS;
RAW<-read.table("$input_file", header=T, sep="\t")
transfers<-as.matrix(RAW[1]*$transfer_multiplier);
RAW<-as.matrix(RAW)
RAW<-RAW[,-1]
is.na(RAW) <- RAW==0
EOS

if ($fit_file)
{
	print R_SCRIPT "FIT<-read.table(\"$fit_file\", header=T, sep=\"\t\")\n";
}

## convert according to data mode
if ($data_mode eq 'ratio')
{
	print R_SCRIPT "RAW = log10(RAW[i])\n";
}
elsif ($data_mode eq 'fraction')
{
	print R_SCRIPT "RAW = log10(RAW/(1-RAW))\n";
}
elsif ($data_mode eq 'log10_ratio')
{
}
elsif ( ($data_mode eq 'log_ratio') || ($data_mode eq 'ln_ratio') )
{
	print R_SCRIPT "RAW = RAW/log(10)\n";
}
elsif ($data_mode eq 'log2_ratio')
{
	print R_SCRIPT "RAW = RAW*(log2)/log(10)\n";
}
else
{
	print STDERR "Unrecognized data format $data_mode\n";
	pod2usage(-exitstatus => 0, -verbose => 2);
}

print R_SCRIPT <<EOS;

mycolors = c("black", "red", "blue", "green", "orange", "brown", "grey", "purple", "darkgreen", "cyan", "magenta", "tan");
mybgcolors = c("black", "red", "blue", "green", "orange", "brown", "grey", "purple", "darkgreen", "cyan", "magenta", "tan");

EOS

if ($fit_file && $distribution_file)
{
print R_SCRIPT <<EOS;

pdf("$output_file", height=6, width=6)
par(mar=c(3,3,2,2))
DIST<-read.table("$distribution_file", header=T, sep="\t")
	
plot(1, type="n", axes=F, xlim=c(min(FIT\$alpha, DIST\$alpha),max(FIT\$alpha, DIST\$alpha)), ylim=c(min(FIT\$tau, DIST\$tau),max(FIT\$tau, DIST\$tau)),  xlab="", ylab="")
box(which = "plot", lty = "solid")
axis(side=1, cex.axis=1.6)
axis(side=2,las = 1, cex.axis=1.6)
points(DIST\$alpha, DIST\$tau, pch="+", bg="white", col="grey", cex=1)
points(FIT\$alpha, FIT\$tau, pch="+", col="red", cex=1.5)

EOS


}
elsif ($fit_file)
{
	print R_SCRIPT <<EOS;
pdf("$output_file", height=6, width=10)
par(mar=c(3,3,2,2))
max_gen <- $max_generations	
if (max_gen == 0) 
{ 
	max_gen <- max(transfers)
}
plot(1, type="n", axes=F, xlim=c(min(transfers),max_gen), ylim=c(-2, 2),  xlab="", ylab="", xaxs = "i", yaxs = "i")
box(which = "plot", lty = "solid")
axis(side=1, at=c(min(transfers):(max_gen/80)*80), cex.axis=1.6)
axis(side=2, at=c((-2:2)*1), las = 1, cex.axis=1.6)

for (i in 1:ncol(RAW)) { 

	THEORETICAL_X = c(1:(transfers[nrow(transfers)]*10))/10
	THEORETICAL_X_GEN<-THEORETICAL_X / $transfer_multiplier
	THEORETICAL_Y = FIT\$log_baseline_ratio[i] + THEORETICAL_X_GEN * log(1+FIT\$initial_fitness_difference[i]) + FIT\$sign[i] * log(1+FIT\$winner_init_fract[i]*exp(FIT\$alpha[i] * (THEORETICAL_X_GEN - FIT\$tau[i])))

	print(paste(FIT\$log_baseline_ratio[i], FIT\$initial_fitness_difference[i], FIT\$sign[i], FIT\$winner_init_fract[i], FIT\$alpha[i], FIT\$tau[i], sep = ", "))
	THEORETICAL_Y = THEORETICAL_Y / log(10);
	print(THEORETICAL_Y)
	lines(THEORETICAL_X, THEORETICAL_Y, col=mycolors[(i-1)%%length(mycolors)+1], lwd=2)
}

for (i in 1:ncol(RAW)) { 

	num_fit_pts<-FIT\$fit_pts[i]
	print(num_fit_pts)
	FIT_TRANSFERS<-transfers[c(1:num_fit_pts)] 
	FIT_PTS<-RAW[c(1:num_fit_pts),i]
	points(FIT_TRANSFERS, FIT_PTS, pch=21+i%%5, bg=mybgcolors[(i-1)%%length(mybgcolors)+1], col=mycolors[(i-1)%%length(mycolors)+1], cex=1.4)
	
	if (num_fit_pts+1<=nrow(RAW)) { 
		NOT_FIT_TRANSFERS<-transfers[c((num_fit_pts+1):nrow(RAW))]
		NOT_FIT_PTS<-RAW[c((num_fit_pts+1):nrow(RAW)),i]
		points(NOT_FIT_TRANSFERS, NOT_FIT_PTS, pch=21+i%%5, col=mycolors[(i-1)%%length(mycolors)+1], cex=1.4, lwd=1.5)
	}

}

EOS
}
else
{
	print R_SCRIPT <<EOS
pdf("$output_file", height=6, width=10)
par(mar=c(3,3,2,2))
plot(1, type="n", axes=F, xlim=c(min(transfers),max(transfers)), ylim=c(-2.5, 2.5), xlab="", ylab="", xaxs = "i", yaxs = "i", cex.axis=1.6, mgp=c(2, 1, 0))
box(which = "plot", lty = "solid")
axis(side=1, at=c(min(transfers):(max(transfers)/80)*80), cex.axis=1.6)
axis(side=2, at=c((-2:2)*1), las = 1, cex.axis=1.6)
for (i in 1:ncol(RAW)) { 
	lines(transfers, RAW[,i],col=mycolors[(i-1)%%length(mycolors)+1], cex=1)
}
for (i in 1:ncol(RAW)) { 
		points(transfers, RAW[,i], pch=21+i%%5, bg=mybgcolors[(i-1)%%length(mybgcolors)+1], col="white", cex=1.5, lwd=2)

}
EOS

}

print R_SCRIPT <<EOS;
dev.off()
EOS

close R_SCRIPT;

## run R
my $command = "R --vanilla < /tmp/$$.r_script > /tmp/$$.r_output";
my $res = system $command;
die ($command) if ($res);	


