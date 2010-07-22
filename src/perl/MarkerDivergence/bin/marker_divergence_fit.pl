#!/usr/bin/perl -w

###
# Pod Documentation
###

=head1 NAME

marker_divergence_fit.pl

=head1 SYNOPSIS

Usage: marker_divergence_fit.pl -i file.tab -o file.fit [-p plot.pdf -m ratio|log10_ratio|log_ratio|fraction -b 5]

Fit alpha and tau estimates for each column of a data file
containing ratios of R/W colonies at generation intervals.

Input format is tab delimited. The first line must start with 'transfer' as the heading 
for the first column, as the unit of time. Other columns headed by the name of series.

Note that calculated alpha values are in units of per transfer and tau values are in units 
of transfers. To compare these to per generation values multiply and divide by the number 
of generations per transfer, respectively.

=head1 DESCRIPTION

=over

=item B<-i> <file path> REQUIRED

Input data file.

=item B<-o> <file path> REQUIRED

Output data file.

=item B<-p> <file path> OPTIONAL

Output file in PDF series with fits plotted to each data series.

=item B<-m|--data-mode> <mode> 

Format of data. Either 'ratio', 'log_ratio', 'log10_ratio', or 'fraction'. Defaults to 'ratio'.

=item B<-b|--fit-baseline> <int> 

Fit a baseline from this many initial points. (For when the initial red/white ratio is not equal to 1.0).
Tau and alpha values are corrected for the new baseline. If a comma delimitted list of values with one
for each series is given (no spaces), then a different number of baseline points can be specified for each series.

=item B<-w|--initial-fitness-difference> fit|float

Difference in relative fitness of the two strains. Note that it should be in units of 1/transfers if transfers
are the unit of time in the marker divergence data. Setting this value to the string "fit" causes the initial
fitness difference to be fit from the baseline points. DEFAULT 0.

=item B<-r|--residual-standard-error> double

Reject fits with a residual standard error greater than this value. DEFAULT 0.15.

=item B<-l|--lilliefors-p-value> double

Required p-value for the Lilliefors test. Fits that reject the null hypothesis that the residuals 
are normally distributed at the input level are rejected, unless none of the fits for a given marker
trajectory reach this threshold. Requires that R have the "nortest" module installed. DEFAULT 0.05. 
Setting to zero disables this test.

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

my $data_mode = 'ratio';
my $residual_error_accept = 0.15;
my $lillie_accept = 0.05; #0 = off

my $idealized_lillie;
my $calculate_durbin_watson;
my $initial_offset = 0;
my $plot;


##currently unused -- should allow changing the intervals
##of initial estimates for both parameters.
my $entered_alpha_estimate;
my $entered_tau_estimate;

my $alpha_start = 10;
my $alpha_factor = 10;
my $alpha_end = 0.01;

my $tau_start = 0;
my $tau_end = 0;
my $tau_add = 20;

my $initial_fitness_difference = 0;
my $gaussian_noise_stdev = 0;

my $baseline_pts = 0;
my $minimum_out_of_bounds_pts = 0;

my ($help, $man, $verbose);
my $input_file;
my $output_file;
GetOptions(
	'help|?' => \$help, 'man' => \$man,
	'input|i=s' => \$input_file,
	'output|o=s' => \$output_file,
	'baseline-pts|b=s' => \$baseline_pts,
	'verbose|v' => \$verbose,
	'data-mode|m=s' => \$data_mode,
	'lilliefors-p-value|l=s' => \$lillie_accept,
	'residual-standard-error|e=s' => \$residual_error_accept,
	'durbin-watson|d' => \$calculate_durbin_watson,
	'alpha-estimate|a=s' => \$entered_alpha_estimate,
	'tau-estimate|t=s' => \$entered_tau_estimate,
	'plot|p=s' => \$plot,
	'gaussian-noise-stdev|g=s' => \$gaussian_noise_stdev,
	'initial-fitness-difference|w=s' => \$initial_fitness_difference,
	'minimum-out-of-bounds-points|0=s' => \$minimum_out_of_bounds_pts,
);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
pod2usage(1) if (!defined $input_file || !defined $output_file);

##Minimum pts requires us to use at least this many pts if the residual std error criterion is not met
my @minimum_out_of_bounds_pts_list;
@minimum_out_of_bounds_pts_list = split /,/, $minimum_out_of_bounds_pts if ($minimum_out_of_bounds_pts =~ m/,/);

##if baseline had commas, use a different number for each one.
my @baseline_pts_list;
@baseline_pts_list = split /,/, $baseline_pts if ($baseline_pts =~ m/,/);

my $fit_initial_fitness_difference;
if (defined $initial_fitness_difference)
{
	if ($initial_fitness_difference eq "fit")
	{
		$fit_initial_fitness_difference = 1;
	}
	else
	{
		## convert from units to the change in log_ratio per transfer this will cause.
		$initial_fitness_difference = log(1+$initial_fitness_difference);
	}
}

## Create output file
open OUTPUT, ">$output_file" or die "Could not open $output_file";

## Load all of the series present in the input file
open INPUT, "<$input_file" or die "Could not open $input_file";
my @lines = <INPUT>;
chomp @lines;

print STDERR "PID = $$\n" if ($verbose);
my %done_hash;

my @header_list = split /\t/, shift @lines;

my %series;
my @generations;
foreach my $line (@lines)
{
	my @line_list = split /\t/, $line;
	push @generations, $line_list[0];
	for (my $i=1; $i<scalar @header_list; $i++)
	{
		#terminate collecting this column
		#if it fixed either red or white
		next if ($done_hash{$i});
		
		my $value = $line_list[$i];
		$value = '' if (!defined $value);
		
		if (!missing_value($value))
		{		
			if ($data_mode eq 'ratio')
			{
				($value > 0) or die "In data set $header_list[$i] found invalid value $value for data mode \"$data_mode\".\n";
				$value = log($value);
			}
			elsif ($data_mode eq 'fraction')
			{
				if ( ($line_list[$i] == 0) || ($line_list[$i] == 1) ) 
				{
					$done_hash{$i} = 1;
					next;
				}
					
				#keep blank intermediate values (missing measurements), they will be skipped later
				($value >= 0) or die "In data set $header_list[$i] found invalid value $value for data mode \"$data_mode\".\n";
				($value <= 1) or die "In data set $header_list[$i] found invalid value $value for data mode \"$data_mode\".\n";
				$value = log($value / (1-$value));
			}
			elsif ($data_mode eq 'log10_ratio')
			{
				$value = log(10) * $value;
				#print STDERR "$line_list[$i] $value\n";
			}
			elsif ( ($data_mode eq 'log_ratio') || ($data_mode eq 'ln_ratio') )
			{
				$value = $value;
				#print STDERR "$line_list[$i] $value\n";
			}
			elsif ($data_mode eq 'log2_ratio')
			{
				$value = log(2)*$value;
				#print STDERR "$line_list[$i] $value\n";
			}
			else
			{
				print STDERR "Unrecognized data format";
				pod2usage(-exitstatus => 0, -verbose => 2) if $man;
			}
			
			$value += gaussian($gaussian_noise_stdev, 0) if ($gaussian_noise_stdev);
		}
		push @{$series{$header_list[$i]}}, $value;
	}

}

#print STDERR Dumper(@generations,%series);

##For each series fit the data until we have a best estimate of tau and alpha
##along with an optional graph for quality control purposes

## Print output header
print OUTPUT "#SETTINGS:\n";
print OUTPUT "#input marker ratio format = $data_mode\n";
print OUTPUT "#residiual standard error (R-squared) accepted <= $residual_error_accept\n";
print OUTPUT "#Lilliefors p-value accepted > $lillie_accept\n";
print OUTPUT "#initial pts used to determine baseline = $baseline_pts\n";
print OUTPUT "#initial fitness difference (per transfer) $initial_fitness_difference\n";
print OUTPUT "#alpha estimates: from $alpha_start to $alpha_end by factors of $alpha_end\n";
print OUTPUT "#tau estimates: from $tau_start to $tau_end by addition of $tau_add\n";

print OUTPUT +(join "\t", ("line", "sign", "log_baseline_ratio", "winner_init_fract", "initial_fitness_difference", "tau", "alpha", "baseline_pts", "min_out_of_bounds_pts", "fit_pts", "total_pts", "residual_standard_error"));
print OUTPUT "\tlilliefors_p_value" if ($lillie_accept);
print OUTPUT "\tintercept_p\tslope_p" if ($fit_initial_fitness_difference);
print OUTPUT "\tdurbin_watson" if ($calculate_durbin_watson);
print OUTPUT "\n";

my @fits;
my $baseline_variance = 0;
my $total_baseline_pts = 0;

SERIES: foreach my $i (1..$#header_list)
{
	print STDERR "==Fitting $header_list[$i]==\n" if ($verbose);

	my @s = @{$series{$header_list[$i]}};

	### Create input file for R
	open R_INPUT, ">/tmp/$$.r_input.tab" or die "Could not open /tmp/$$.r_input.tab";
	print R_INPUT"$header_list[0]\tlog_ratio\n";
	PRINTED_POINTS: for (my $k=0; $k<scalar @s; $k++)
	{
		next PRINTED_POINTS if ( missing_value($s[$k]) ); #Skip missing values
		print R_INPUT +($generations[$k]) . "\t";
		print R_INPUT +($s[$k]) . "\n";
	}
	close R_INPUT;

	###Optional baseline
	my $baseline = 0.0;
	my $winner_initial_fraction = 0.5;
	my $this_initial_fitness_difference = 0;
		
	my $this_baseline_pts = $baseline_pts;	
	if (@baseline_pts_list)
	{
		die "Not enough values for baseline pts provided in comma-delimitted format.\n" if (!defined $baseline_pts_list[$i-1]);
		$this_baseline_pts = $baseline_pts_list[$i-1];
	}
	
	## Fit baseline pts to get initial fitness difference 
	my $intercept_p;
	my $slope_p;

	if ($fit_initial_fitness_difference)
	{
		## create R script
		open R_SCRIPT, ">/tmp/$$.r_script";
		print R_SCRIPT "X <- read.table(\"/tmp/$$.r_input.tab\", header = TRUE)\n";
		print R_SCRIPT "X_sub <- subset(X, $header_list[0] <= $generations[$this_baseline_pts-1])\n";
		print R_SCRIPT "lmfit <- lm(log_ratio ~ $header_list[0], X_sub)\n";
		print R_SCRIPT "summary(lmfit)\n";
		close R_SCRIPT;
		
		## run R
		my $command = "R --vanilla < /tmp/$$.r_script > /tmp/$$.r_output";
		my $res = system $command;
		die ($command) if ($res);			
					
		## track down the stats in the R output file...
		open R_OUTPUT, "</tmp/$$.r_output" or die ("Could not open /tmp/$$.r_output");
		my @output_lines = <R_OUTPUT>;
		close R_OUTPUT;
		chomp @output_lines;

		## parse R output
		my $intercept;
		my $slope;
		
		foreach my $line (@output_lines)
		{
			if ($line =~ m/^\(Intercept\)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/)
			{
				($intercept, $intercept_p) = ($1, $4);
			}
			elsif ($line =~ m/^$header_list[0]\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/)
			{
				($slope, $slope_p) = ($1, $4);
			}
		}		
		#print STDERR Dumper(@output_lines);
		$this_initial_fitness_difference = $slope;
		$baseline = $intercept;
		$winner_initial_fraction = exp($baseline)/(1+exp($baseline));
	}
	## Fit baseline pts as line with no slope or user supplied slope from initial fitness difference	
	else 
	{
		$this_initial_fitness_difference = $initial_fitness_difference;
		if ($this_baseline_pts)
		{
			my $baseline_pts = 0;
			
			for (my $b=0; $b<$this_baseline_pts; $b++)
			{
				next if ( missing_value($s[$b]) );

				die "In data set $header_list[$i]: could not find the requested number of baseline points." if ($b >= scalar @s);
				if ($s[$b] ne '')
				{
					$baseline_pts++;
					$baseline += $s[$b] - log(2) * $generations[$b] * $this_initial_fitness_difference;
				}
			}
			die "In data set $header_list[$i]: could not find any baseline points." if ($baseline_pts == 0);
			$baseline /= $baseline_pts;
			$winner_initial_fraction = exp($baseline)/(1+exp($baseline));
		}
	}
	
	#keep track of variance
	for (my $b=0; $b<$this_baseline_pts; $b++)
	{
		next if (missing_value($s[$b]));
		
		my $diff = $baseline + log(2) * $generations[$b] * $this_initial_fitness_difference - $s[$b];
		$baseline_variance += $diff * $diff;
		$total_baseline_pts++;
	}
		
	#$baseline = sprintf("%.5e", $baseline);
	my $total_pt = 0;

	### initialize alpha estimate
	my $converged = 0;
		
	### Create R script, needs to be repeated for each alpha estimate 
	open R_SCRIPT, ">/tmp/$$.r_script";
	print R_SCRIPT "library(nortest)\n" if ($lillie_accept);
	print R_SCRIPT "library(lawstat)\n" if ($calculate_durbin_watson);
	print R_SCRIPT "X <- read.table(\"/tmp/$$.r_input.tab\", header = TRUE)\n";
	
	my $this_minimum_out_of_bounds_pts = $minimum_out_of_bounds_pts;	
	if (@minimum_out_of_bounds_pts_list)
	{
		die "Not enough values for baseline pts provided in comma-delimitted format.\n" if (!defined $minimum_out_of_bounds_pts_list[$i-1]);
		$this_minimum_out_of_bounds_pts = $minimum_out_of_bounds_pts_list[$i-1];
	}
	
	PRINTED_POINTS: for (my $k=$this_baseline_pts+1; $k<scalar @s; $k++)
	{
		next PRINTED_POINTS if (missing_value($s[$k]));
		next PRINTED_POINTS if ($s[$k] == $baseline); #no point if no divergence at all!
		
		print R_SCRIPT "result<-0\n";
		
		my $on_alpha_estimate = $alpha_start;	
		ALPHA: while ($on_alpha_estimate >= $alpha_end)
		{
		
			my $on_tau_estimate = $tau_start;	
			TAU: while ($on_tau_estimate <= $tau_end)
			{
				$total_pt = $k+1;
				
				my $compare_baseline = 	$baseline + $this_initial_fitness_difference * $generations[$k];			
				my $sign = ($s[$k] > $compare_baseline) ? +1 : -1;
				my $use_baseline = $baseline; #don't flip b/c ratio is flipped in this case
				my $use_winner_initial_fraction = ($sign == +1) ? $winner_initial_fraction : 1-$winner_initial_fraction;
				my $use_fitness_difference = $this_initial_fitness_difference;
				
				### Try several values of alpha, which seems to guarantee convergence	
				my $tau_estimate = $on_tau_estimate;
				my $alpha_estimate = $on_alpha_estimate;
				
				print STDERR "Points Included: $total_pt // baseline = $baseline // winner_init_fr = $use_winner_initial_fraction // sign = $sign // alpha_est = $alpha_estimate // tau_est = $tau_estimate\n" if ($verbose);

				print R_SCRIPT "doit <- function()\n";
				print R_SCRIPT "{\n";
				print R_SCRIPT "X_sub <- subset(X, $header_list[0] <= $generations[$k])\n";
				print R_SCRIPT"print(X_sub)\n";
				
				my $to_fit = ($sign == +1) ? "ratio" : "ratio_reciprocal";
				
				print R_SCRIPT "nlsfit <- nls(log_ratio ~ $use_baseline + $use_fitness_difference * $header_list[0] + $sign * log(1 + $use_winner_initial_fraction * exp(alpha * ($header_list[0] - tau))), control = nls.control(maxiter = 100000, minFactor = 1/1024), algorithm=\"port\", data = X_sub, start=list(alpha=$alpha_estimate, tau = $tau_estimate), trace = FALSE)\n";
				##Note, the default algorithm sometimes doesn't return fits (no output at all...) When you run that series alone, it then works. Odd and annoying, but "port" seems to work around this.
				print R_SCRIPT "print(summary(nlsfit))\n";
				print R_SCRIPT "print($total_pt)\n";
				print R_SCRIPT "print($sign)\n";
				print R_SCRIPT "print($alpha_estimate)\n";
				print R_SCRIPT "print($tau_estimate)\n";

				##calculate optional goodness-of-fit statistics		
				if ($lillie_accept)
				{
					print R_SCRIPT "L <- lillie.test(as.vector(residuals(nlsfit)))\n";
					print R_SCRIPT "print(L\$p.value)\n";
				}		
						
				if ($calculate_durbin_watson)
				{ 		
					print R_SCRIPT "res<-as.vector(residuals(nlsfit))\n";

					print R_SCRIPT "crse<-abs(res[length(res)]/(log10($s[$k])-log10($use_baseline)))\n";
					print R_SCRIPT "print(crse)\n";

					#print R_SCRIPT "b<-bartels.test(res, alternative=\"positive.correlated\")\n";
					#print R_SCRIPT "print(b\$p.value)\n";	
					
					#	print R_SCRIPT "d<-durbin.watson(res, alternative=\"two.sided\", max.lag=1, method=\"resample\")\n";
					#	print R_SCRIPT "print(d)\n";
				}
				print R_SCRIPT "print(\"done\")\n"; #note: done only printed if all were successful
				print R_SCRIPT "return(1)\n"; #indicates successful return
				print R_SCRIPT "}\n";
				
				#only try fitting again if we didn't already converge for a different alpha and tau estimate
				#I have never seen different fit values result from changing the estimates
				print R_SCRIPT "if(!result)\n";
				print R_SCRIPT "{\n";
				print R_SCRIPT "try(result<-doit(), silent = TRUE )\n";
				print R_SCRIPT "}\n";
			} continue {		
				$on_tau_estimate += $tau_add;
			}
		} continue {		
			$on_alpha_estimate /= $alpha_factor;
		}
	}
	close R_SCRIPT;	

	my $command = "R --vanilla < /tmp/$$.r_script > /tmp/$$.r_output";
	my $res = system $command;
	die ($command) if ($res);			
				
	#Track down the stats in the R output file...
	open R_OUTPUT, "</tmp/$$.r_output" or die ("Could not open /tmp/$$.r_output");
	my @output_lines = <R_OUTPUT>;
	close R_OUTPUT;

	#print STDERR "@output_lines";

	chomp @output_lines;
	my ($current_alpha, $current_tau, $current_pt, $current_lilliefors, $current_durbin_watson, $current_sign, $current_residual_standard_error);
	my ($current_alpha_estimate, $current_tau_estimate);
	my ($best_alpha, $best_tau, $best_pt, $best_lilliefors, $best_durbin_watson, $best_sign, $best_residual_standard_error);
	$best_residual_standard_error = 1.0;
	$current_lilliefors = -1;
	$current_durbin_watson = -1;

	my $get_next = 0;
	OUTPUT: foreach $_ (@output_lines)
	{
	
		if ($_ =~ /^\[1\]\s+(\S+)/)
		{
			if ($1 eq "\"done\"")
			{			
				#$current_residual_standard_error = $current_residual_standard_error/$current_durbin_watson;
				print STDERR ">>alpha = $current_alpha_estimate, tau = $current_tau_estimate, alpha = $current_alpha, tau = $current_tau, pts = $current_pt, res std err = $current_residual_standard_error," if ($verbose);
				print STDERR "lilliefors = $current_lilliefors, durbin_watson = $current_durbin_watson\n" if ($verbose);
					
				## fitting is done if lilliefors EVER goes above the cutoff.	
#				if ($lillie_accept && ($current_lilliefors <= $lillie_accept))
#				{
#					print STDERR "Lilliefors exit: $current_lilliefors <= $lillie_accept\n" if ($verbose);
#					last OUTPUT if (defined $best_alpha && ($best_lilliefors > $lillie_accept));
#				}
				
				## fitting is done if residual std error EVER goes above the cutoff.	
#				if ($current_residual_standard_error > $residual_error_accept)
#				{
#					print STDERR "Res std error exit: $current_residual_standard_error > $residual_error_accept\n" if ($verbose);
#					last OUTPUT if (defined $best_alpha && ($best_residual_standard_error <= $residual_error_accept));
#				}
				
					
				if ($current_alpha > 0)
				{				
					if (($current_residual_standard_error <= $residual_error_accept) || ($current_residual_standard_error < $best_residual_standard_error))
					{
					
						if (!$lillie_accept || ($current_lilliefors > $lillie_accept) || !defined($best_lilliefors) || ($current_lilliefors > $best_lilliefors) )
						{	
							##require a minimum number of points if it does not meet criteria
							if ( ($current_pt >= $this_minimum_out_of_bounds_pts) || (($current_lilliefors > $lillie_accept) && ($current_residual_standard_error > $best_residual_standard_error)) )
							{
							
								#if (($current_durbin_watson <= 0.05) || ($current_durbin_watson < $best_durbin_watson))
								#{

								($best_alpha, $best_tau, $best_pt, $best_lilliefors, $best_durbin_watson, $best_sign, $best_residual_standard_error) 
									= ($current_alpha, $current_tau, $current_pt, $current_lilliefors, $current_durbin_watson, $current_sign, $current_residual_standard_error);
								print STDERR ">>NEW BEST\n" if ($verbose);
							}
								#}
						}
							
					}
					
					undef $current_pt;
					undef $current_sign;
					undef $current_alpha_estimate;
					undef $current_tau_estimate;
					$current_lilliefors = -1;
					$current_durbin_watson = -1;
				}
			}
			elsif (!defined $current_pt)
			{
				$current_pt = $1;
			}
			elsif (!defined $current_sign)
			{
				$current_sign = $1;
			}
			elsif (!defined $current_alpha_estimate)
			{
				$current_alpha_estimate = $1;
			}
			elsif (!defined $current_tau_estimate)
			{
				$current_tau_estimate = $1;
			}
			elsif ($lillie_accept && ($current_lilliefors == -1))
			{
				$current_lilliefors = $1;
			}
			elsif ($calculate_durbin_watson && ($current_durbin_watson == -1))
			{
				$current_durbin_watson = $1;
			}
		}
				
		if ($_ =~ /^Residual standard error:\s+(\S+)/)
		{
			$current_residual_standard_error = $1;
			
			undef $current_pt;
			undef $current_sign;
			undef $current_alpha_estimate;
			undef $current_tau_estimate;
			$current_lilliefors = -1;
			$current_durbin_watson = -1;
		}
	
		if ($_ =~ /^alpha\s+(\S+)/)
		{
			$current_alpha = $1;
		}
		
		if ($_ =~ /^tau\s*(\S+)/)
		{
			$current_tau = $1;
		}
	}

	### Success!!

	if (defined $best_alpha)
	{
		#warn if residual standard error was high
		if ($best_residual_standard_error > $residual_error_accept)
		{
			print STDERR "Warning: accepted R-squared of $best_residual_standard_error for data series $header_list[$i] exceeds threshold of $residual_error_accept\n";
		}

		my $print_initial_fitness_difference = exp($this_initial_fitness_difference) - 1;
		my $print_winner_initial_fraction = ($best_sign == +1) ? $winner_initial_fraction : 1-$winner_initial_fraction;
		
		print OUTPUT +(join "\t", $header_list[$i], $best_sign, $baseline, $print_winner_initial_fraction, $print_initial_fitness_difference, $best_tau, $best_alpha, $this_baseline_pts, $this_minimum_out_of_bounds_pts, $best_pt, $total_pt, $best_residual_standard_error);
		print OUTPUT "\t$best_lilliefors" if ($lillie_accept);
		print OUTPUT "\t$intercept_p\t$slope_p" if ($fit_initial_fitness_difference);
		print OUTPUT "\t$best_durbin_watson" if ($calculate_durbin_watson);
		print OUTPUT "\n";


		print STDERR "BEST: $best_alpha, $best_tau, $best_pt, $best_sign, $best_residual_standard_error, $best_lilliefors, $best_durbin_watson\n" if ($verbose);
		
		### Save information about the fit for plotting later
		my $fit_data;
		$fit_data->{alpha} = $best_alpha;
		$fit_data->{tau} = $best_tau;
		$fit_data->{sign} = $best_sign;
		$fit_data->{baseline} = $baseline;
		$fit_data->{winner_initial_fraction} = ($fit_data->{sign} == +1) ? $winner_initial_fraction : 1-$winner_initial_fraction;
		$fit_data->{pts} = $best_pt;
		$fit_data->{index} = $i;
		$fit_data->{initial_fitness_difference} = $this_initial_fitness_difference;


		push @fits, $fit_data;
	}
	else
	{
		print STDERR "No fits converged for data series $header_list[$i]\n";
		#die;
	}
}

if ($baseline_pts)
{
	$baseline_variance /= $total_baseline_pts;
	print OUTPUT "#baseline variance: $baseline_variance (n=$total_baseline_pts)\n";
}


#clean up
#unlink "/tmp/$$.r_input.tab";
#unlink "/tmp/$$.r_script";
#unlink "/tmp/$$.r_output";

##One PDF is created with all of the relevant fits.
if ($plot)
{
	print STDERR "Creating R Script for plots: /tmp/$$.R_plot_script\n";
	print STDERR "Experimental data input file: /tmp/$$.exp_plot_values\n";
	print STDERR "Creating fit curve input file: /tmp/$$.fit_plot_values\n";

	open R_PLOT_SCRIPT, ">/tmp/$$.r_plot_script";
	print R_PLOT_SCRIPT <<EOS;
pdf("$plot")
EOS

	foreach my $fit (@fits)
	{
		#write input file containing log10 experimental values
		my @s = @{$series{$header_list[$fit->{index}]}};

		open EXP, ">/tmp/$$.exp_plot_values_$fit->{index}";
		print EXP "transfer\tlog10_ratio\n";
		PRINTED_POINTS: for (my $k=0; $k<scalar @s; $k++)
		{
			#next PRINTED_POINTS if ( $s[$k] eq ''); #Skip missing values
			print EXP +join("\t", $generations[$k], (missing_value($s[$k]) ? 'NA' : ($s[$k]/log(10)))) . "\n";
		}
		close EXP.

		#write fit file with a higher time resolution.
		open FIT, ">/tmp/$$.fit_plot_values_$fit->{index}";
		print FIT "transfer\tlog10_ratio\n";
			
		my $use_baseline = $fit->{baseline}; #don't flip b/c ratio is flipped in this case
		my $use_winner_initial_fraction = $fit->{winner_initial_fraction};
		my $use_fitness_difference = $fit->{initial_fitness_difference};

		my $max_generations = $generations[-1];
		my $num_fit_points = 1000;
		for (my $i=0; $i<=$num_fit_points; $i++)
		{
			my $fit_time = $i * $max_generations / $num_fit_points;
			my $fit_value = $use_baseline + $use_fitness_difference * $fit_time + $fit->{sign} * log(1 + $use_winner_initial_fraction * exp($fit->{alpha}*($fit_time-$fit->{tau})));
			$fit_value = $fit_value / log(10);
			print FIT "$fit_time\t$fit_value\n";		
		}
		close FIT;
						
		print R_PLOT_SCRIPT <<EOS
exp_data<-read.delim("/tmp/$$.exp_plot_values_$fit->{index}", sep = "\t", header=TRUE)	
fit_data<-read.delim("/tmp/$$.fit_plot_values_$fit->{index}", sep = "\t", header=TRUE)	
plot(exp_data\$transfer,exp_data\$log10_ratio, xlab="$header_list[0]", ylab="log10 marker ratio")
title(main="Fit for series \#$header_list[$fit->{index}]", sub="alpha = $fit->{alpha}, tau = $fit->{tau}")
points(exp_data\$transfer[1:$fit->{pts}],exp_data\$log10_ratio[1:$fit->{pts}], col="blue", bg="blue", pch=19)
lines(fit_data\$transfer,fit_data\$log10_ratio)
EOS

	}
	
	## Close the PDF file so that it will be written
	print R_PLOT_SCRIPT "dev.off()\n";

	print STDERR "Running R plot script: /tmp/$$.r_plot_output\n";
	my $command = "R --vanilla < /tmp/$$.r_plot_script > /tmp/$$.r_plot_output";
	my $res = system $command;
	die ($command) if ($res);

	# clean up
	unlink "/tmp/$$.r_plot_script";
	unlink "/tmp/$$.r_plot_output";
	foreach my $fit (@fits)
	{
		unlink "/tmp/$$.exp_plot_values_$fit->{index}";
		unlink "/tmp/$$.fit_plot_values_$fit->{index}";
	}
}



sub missing_value
{
	my ($val) = @_;
	return (($val eq '') || ($val eq 'NA'))
}

