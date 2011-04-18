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

=item B<-d> <file path> 

Path to directory file containing fits of simulated marker divergence data. 
The experimental data is compared to all ".fit" files in this path. REQUIRED

=item B<-r>

Recursively search sub-directories of -d for ".fit" files. OPTIONAL

=item B<-o> <file path> 

Path to tab-delimitted output file that will be created showing the KS test 
results for each input fit file. REQUIRED

=item B<-p> <file path> 

Path to optional output file containing a plot showing regions of significant agreement
between experimental and simulated alpha and tau distributions and the maximum likelihood
values of s and mu. OPTIONAL

=item B<-m> <file path> 

Save the matrix file used to create the plot here. OPTIONAL

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

use IO::Handle;

use KS_test_2d;

my $ratio_mode = 0;
my $experimental_mode;

my ($help, $man);
my ($input_file, $plot_file, $matrix_file, $output_file, $fit_directory, $recursive, $new);
our $detailed_colors = 0;
my $p_value_cutoff = 0.05;
my $local_average;
GetOptions(
	'help|?' => \$help, 'man' => \$man,
	'input|i=s' => \$input_file,
	'fit-directory|d=s' => \$fit_directory,
	'recursive|r' => \$recursive,
	'output-file|o=s' => \$output_file,
	'plot-file|p=s' => \$plot_file,
	'matrix-file|m=s' => \$matrix_file,
	'experimental|x' => \$experimental_mode,
	'new|n' => \$new,
	'p-value|l=s' => \$p_value_cutoff,
	'full-colors|f' => \$detailed_colors,
	'local_average|a' => \$local_average,
);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
pod2usage(1) if (!defined $input_file || !defined $output_file || !defined $fit_directory);
die "Input file does not exist" if (!-e $input_file);


$matrix_file = "$$.matrix" if (!defined $matrix_file);

my @simulation_files;

if (!$recursive)
{
	opendir FITDIR, "$fit_directory";
	@simulation_files = readdir FITDIR;
	@simulation_files = grep /.fit$/, @simulation_files;
	@simulation_files = grep $_="$fit_directory/$_", @simulation_files;	
}
else
{
	@simulation_files = recursive_file_search($fit_directory);	
}


sub by_s_mu
{
 	$a =~ m/s(?:_|=)([0-9.]+)\_mu(?:_|=)([0-9.\-]+)\.fit/;
  	my ($a1, $a2) = ($1, $2);
  	 	
	$b =~ m/s(?:_|=)([0-9.]+)\_mu(?:_|=)([0-9.\-]+)\.fit/;	
	my ($b1, $b2) = ($1, $2);

	return (($a1 <=> $b1) || ($a2 <=> $b2));
}
@simulation_files = sort by_s_mu @simulation_files;

die "No .fit files found in specified location" if (scalar @simulation_files == 0);

my (%mu_hash, %s_hash, %data_hash);
if ($new)
{
	open OUT, ">$output_file";
	OUT->autoflush(1);
	print OUT +(join "\t", ('s', 'mu', 'p-value', 'alpha_sig', '95%_CI', 'file_name')) . "\n";
	
	my ($best, $best_mu, $best_s, $best_is_significant);
	my $max_like;

	### Run the new treatment
	my $experimental_pts_ref = load_fit_file($input_file);
#	print Dumper($experimental_pts_ref);
	
	## Extract s and mu values from the file name
	my %file_hash;
	foreach my $simulation_file (@simulation_files)
	{
		my $simulation_file_name = $simulation_file;
		$simulation_file_name =~ s/^.+\///;
		$simulation_file_name =~ s/\.fit$//;
		$simulation_file_name =~ m/s(_|=)(.+?)(_|$)/;
		my $on_s = sprintf("%.5f", $2);
		$simulation_file_name =~ m/mu(_|=)(.+?)(_|$)/;
		my $on_mu = sprintf("%.5f", $2);
		$file_hash{$on_s}->{$on_mu} = $simulation_file;
	}
	
	my $s_step;
	my $mu_step;
	if ($local_average)
	{
		my @s_list = sort {$a <=> $b} keys %file_hash;
		$s_step = $s_list[1]-$s_list[0];
		
		my @mu_list = sort {$a <=> $b} keys %{$file_hash{$s_list[0]}};
		$mu_step = $mu_list[1]-$mu_list[0];
	}
	
	foreach my $simulation_file (@simulation_files)
	{
		## what is the file name?
		my $simulation_file_name = $simulation_file;
		$simulation_file_name =~ s/^.+\///;
		$simulation_file_name =~ s/\.fit$//;
		
		$simulation_file_name =~ m/s(_|=)(.+?)(_|$)/;
		my $on_s = $2;
		$s_hash{$on_s}++;

		$simulation_file_name =~ m/mu(_|=)(.+?)(_|$)/;
		my $on_mu = $2;
		$mu_hash{$on_mu}++;
	
		my $simulation_pts_ref;
		if (!$local_average)
		{
			$simulation_pts_ref = load_fit_file($simulation_file);
		}
		else
		{		
			my @offsets = (-1, 0, +1);
			foreach my $o1 (@offsets)
			{
				my $test_s = sprintf("%.5f", $on_s + $o1 * $s_step); 
				foreach my $o2 (@offsets)
				{
					my $test_mu = sprintf("%.5f", $on_mu + $o2 * $mu_step); 
					if (defined $file_hash{$test_s}->{$test_mu})
					{
						my $add_pts_ref = load_fit_file($file_hash{$test_s}->{$test_mu});
						push @$simulation_pts_ref, @$add_pts_ref;
					}
				}
			}

			my @add_list		
		}
				
		my $p_value = bf_KS_test_2d($experimental_pts_ref, $simulation_pts_ref);
		
		print "$simulation_file $p_value\n";
		
		## is it significant?
		my $significance_code = '';
		if ($p_value >= $p_value_cutoff)
		{
			$significance_code = "***";
		}		
		$data_hash{"$on_mu\_$on_s"} = ($significance_code) ? 1 : 0;

		##new color scheme
		if ($detailed_colors)
		{
			$data_hash{"$on_mu\_$on_s"} = $p_value;
		}
		else ##simple color scheme
		{
			$data_hash{"$on_mu\_$on_s"} = ($significance_code) ? 1 : 0;
		}

		if ( (!defined $best) || ($significance_code && !$best_is_significant) || ($p_value > $best) )
		{
			$max_like = $simulation_file_name;
			$best = $p_value;
			$best_s = $on_s;
			$best_mu = $on_mu;
			$best_is_significant = $significance_code;
		}
		
		print OUT +join("\t", $on_s, $on_mu, $p_value, $significance_code, $simulation_file_name) . "\n";
	}

	my $sig_string = ($best_is_significant) ? "" : " NOT";
	print OUT "Maximum Likelihood: mu = $best_mu, s = $best_s\n";
	print OUT "Compound significance = $best, and is $sig_string within the 95% CI.\n";
	close OUT;
	
	##simple color scheme, best point is special
	$data_hash{"$best_mu\_$best_s"} = 2 if (!$detailed_colors && $best_is_significant);
}
else
{


	###
	# Create the R script file
	###
	print STDERR "Creating R script...\n";

	open R_SCRIPT, ">$$.r_script";
	print R_SCRIPT <<EOS;
experimental <- read.table("$input_file", header = TRUE)
EOS

	print R_SCRIPT "library(Matching)\n" if ($experimental_mode);

	foreach my $simulation_file (@simulation_files)
	{
		print R_SCRIPT	(!$experimental_mode ? <<END1 : <<END2);
theoretical <- read.table("$simulation_file", header = TRUE)
tt<-ks.test( experimental\$tau, theoretical\$tau, alternative = c("two.sided") )
print(tt\$p.value)
at<-ks.test( experimental\$alpha, theoretical\$alpha, alternative = c("two.sided") )
print(at\$p.value)
END1
theoretical <- read.table("$simulation_file", header = TRUE)
tt<-ks.boot( experimental\$tau, theoretical\$tau, alternative = c("two.sided") )
print(tt\$ks\$p.value)
at<-ks.boot( experimental\$alpha, theoretical\$alpha, alternative = c("two.sided") )
print(at\$ks\$p.value)
END2

	}
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

	#Parse the output for p-values
	chomp @lines;
	@lines = grep s/^\[1\] //, @lines;
	{
		open OUT, ">$output_file";
		print OUT +(join "\t", ('s', 'mu', 'tau_sig', 'alpha_sig', '95%_CI', 'file_name')) . "\n";

		my ($best, $best_mu, $best_s, $best_is_significant);
		my $max_like;

		foreach my $simulation_file (@simulation_files)
		{
			## get the next significance values
			my $tau_significance = shift @lines;
			my $alpha_significance = shift @lines;
			my $compound_significance = $alpha_significance * $tau_significance;

			## what is the file name?
			my $simulation_file_name = $simulation_file;
			$simulation_file_name =~ s/^.+\///;
			$simulation_file_name =~ s/\.fit$//;
			
			## Extract s and mu values from the file name
			$simulation_file_name =~ m/s(_|=)(.+?)(_|$)/;
			my $on_s = $2;
			$s_hash{$on_s}++;

			$simulation_file_name =~ m/mu(_|=)(.+?)(_|$)/;
			my $on_mu = $2;
			$mu_hash{$on_mu}++;
			
			## is it significant?
			my $significance_code = '';
			if ( ($tau_significance > 0.025) && ($alpha_significance > 0.025) )
			{
				$significance_code = "***";
			}

			$data_hash{"$on_mu\_$on_s"} = ($significance_code) ? 1 : 0;

			if ( (!defined $best) || ($significance_code && !$best_is_significant) || ($compound_significance > $best) )
			{
				$max_like = $simulation_file_name;
				$best = $compound_significance;
				$best_s = $on_s;
				$best_mu = $on_mu;
				$best_is_significant = $significance_code;
			}
			
			print OUT +join("\t", $on_s, $on_mu, $tau_significance, $alpha_significance, $significance_code, $simulation_file_name) . "\n";
		}
		my $sig_string = ($best_is_significant) ? "" : " NOT";
		print OUT "Maximum Likelihood: mu = $best_mu, s = $best_s\n";
		print OUT "Compound significance = $best, and is $sig_string within the 95% CI.\n";
		close OUT;
		$data_hash{"$best_mu\_$best_s"} = (($detailed_colors) ? 6 : 2) if ($best_is_significant);
	}

	#clean up 
	unlink "$$.r_script";
	unlink "$$.r_output";
}

#### Continue to plotting if requested
if ($plot_file)
{
	##
	# Create data matrix for R to read
	## 
	my $mu_string = join(",", sort {$a <=> $b} keys %mu_hash);
	my $s_string = join(",", sort {$a <=> $b} keys %s_hash);

	print "$mu_string\n";
	print "$s_string\n";
	open MATRIX, ">$matrix_file";
	foreach my $s (sort {$a <=> $b} keys %s_hash)
	{
		my @matrix_line;
		
		foreach my $mu (sort {$a <=> $b} keys %mu_hash)
		{
			my $item = $data_hash{"$mu\_$s"};
			
			if (!defined $item)
			{
				print STDERR "Missing matrix item: mu = $mu, s = $s\n";
				$item = 'NA'; #R code for missing value
			}
			
			push @matrix_line, $item;
		}
		print MATRIX +join("\t", @matrix_line) . "\n";
	}

	##
	# Plot
	## 
	print STDERR "Creating R plot script: $$.plot.r_script\n";
	open R_PLOT_SCRIPT, ">$$.plot.r_script";
	
	my $colors = "\"gray90\", \"blue\", \"black\"";
	
	##simple plot (3 colors)
	if (!$detailed_colors)
	{	
		print R_PLOT_SCRIPT <<EOS;
library(fields);
data <- read.table("$matrix_file", header = FALSE)
pdf("$plot_file", height=6, width=6)
#my_colors = c("gray90", "blue", "black"); to show missing data!
my_colors = c("white", "blue", "black");
image(c($s_string), c($mu_string), as.matrix(data), zlim=c(0,2), xlab="s", ylab="mu", las=1, col=my_colors)
dev.off()
EOS
	}
	else
	{
		print R_PLOT_SCRIPT <<EOS;
library(fields);
data <- read.table("$matrix_file", header = FALSE)
pdf("$plot_file", height=6, width=6)
my_colors = tim.colors(100);
my_colors[100] = "gray95";
my_colors[99] = "gray80";
my_colors[98] = "gray70";
my_colors[97] = "gray60";
my_colors[96] = "gray50";
my_colors = rev(my_colors)
br<-0:100;
brlab<-sapply(br,function(x) if (x%%5==0) x/100 else NA)
br<-br/100;
image.plot(c($s_string), c($mu_string), as.matrix(data), zlim=c(0,1), breaks=br, lab.breaks=brlab, xlab="s", ylab="mu", nlevel=100, las=1, col=my_colors)
dev.off()
EOS
	}
	
#for plotting CI legend rather than p-values	
#remove: my_colors = rev(my_colors)
#brlab<-sapply(br,function(x) if (x%%5==0) paste(x,"%", sep="") else NA)	
	close R_PLOT_SCRIPT;

	#Run R
	print STDERR "Running R plot script: $$.plot.r_script\n";
	my $command = "R --vanilla < $$.plot.r_script &> $$.plot.r_output";
	my $res = system $command;
	die "Error running command: $command" if ($res);
	
	#clean up
	unlink "$$.plot.r_script";
	unlink "$$.plot.r_output";
	unlink $matrix_file if ($matrix_file eq "$$.matrix");
}



#foreach my $simulation_file (@simulation_files)
#{
#	print R_SCRIPT	(!$experimental_mode ? <<END1 : <<END2);
#theoretical <- read.table("$simulation_file", header = TRUE)
#tt<-ks.test( experimental\$tau, theoretical\$tau, alternative = c("two.sided") )
#print(tt\$p.value)
#at<-ks.test( experimental\$alpha, theoretical\$alpha, alternative = c("two.sided") )
#print(at\$p.value)
#END1
#theoretical <- read.table("$simulation_file", header = TRUE)
#tt<-ks.boot( experimental\$tau, theoretical\$tau, alternative = c("two.sided") )
#print(tt\$ks\$p.value)
#at<-ks.boot( experimental\$alpha, theoretical\$alpha, alternative = c("two.sided") )
#print(at\$ks\$p.value)
#END2
#}
#close R_SCRIPT;


####
sub recursive_file_search
{
	my ($directory_name) = @_;
	my @simulation_files = ();

	#print STDERR "Examining $directory_name\n";
	opendir DIR, "$directory_name";
	my @dir_files = readdir DIR;
	@dir_files = grep !/^\./, @dir_files;
	
	#print STDERR "@dir_files\n";
	
	foreach my $dir_item (@dir_files)
	{
		my $full_name = "$directory_name/$dir_item";
		push @simulation_files, recursive_file_search($full_name) if (-d $full_name);
		if ($full_name =~ m/.fit$/)
		{
			push @simulation_files, $full_name;
		}
	}	
	return @simulation_files;

}

sub load_fit_file
{
	my ($file_name) = @_;
	
	open FIT, "<$file_name" or die "could not open $file_name\n";
	print STDERR "Loading $file_name\n";
	my @lines = <FIT>;
	
	@lines = grep !/^\s*#/, @lines;
	@lines = grep !/^\s*$/, @lines;
	
	my $header_line = shift @lines;
	my @header_items = split "\t", $header_line;
	my %header_hash;
	for (my $i=0; $i< scalar @header_items; $i++)
	{
		$header_hash{$header_items[$i]} = $i;
	}

	my @pts;
	foreach my $line (@lines)
	{
		my @split_line = split "\t", $line;
		my $pt;
		($pt->{x}, $pt->{y}) = ($split_line[$header_hash{alpha}], $split_line[$header_hash{tau}]);
		push @pts, $pt;
	}
	
	return \@pts;
}
