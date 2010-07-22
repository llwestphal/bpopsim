#!/usr/bin/perl -w

###
# Pod Documentation
###

=head1 NAME

marker_divergence_background_model.pl 

=head1 SYNOPSIS

Usage: marker_divergence_background_model.pl -T 6.64 -N 5E6 -s 0.04:0.1:0.01 -u -8:-6:0.5 -p pr_establishment_T_6.64_No_5E6.tab

Perform simulations of marker divergence data at specified s and mu pairs. Fit alpha and tau empirical parameters for each s and mu pair.

=head1 DESCRIPTION

=over

=item B<Marker Trajectory Simulation Options>

=over

=item B<-T> <double> 

Number of generations per growth phase. Required.

=item B<-N> <double> 

Initial population size at each growth cycle (after dilution). N sub zero. Required.

=item B<-u> <double> 

Log10 of the beneficial mutation rate. To perform runs with different 
parameters, use the form start:end:step. Required.

=item B<-s> <double> 

Beneficial mutation selection coefficient. To perform runs with different 
parameters, use the form start:end:step. Required.

=item B<-b|--beneficial-mutation-distribution> <name of distribution>

Probability distribution of beneficial mutation selection coefficients. 
Either 'dirac_delta' or 'exponential'. Default = dirac_delta

=item B<-r> <int> 

Number of replicates to simulate and fit for each parameter pair. Default = 100.

=item B<-i> <int> 

Time interval for which marker ratio data is available (in transfers). Default = 1.

=item B<-p> <file path> 

File containing probabilities of a new beneficial mutation not going extinct. 

=item B<-k|--skip-generations> <double>

Skip this many initial generations before printing data. Used to account for initial outgrowth. Default = 0.

=back

=item B<Marker Trajectory Fitting Options>

=over

=item B<-r|--residual-standard-error> double

Reject fits with a residual standard error greater than this value. DEFAULT 0.15.

=item B<-l|--lilliefors-p-value> double

Required p-value for the Lilliefors test. Fits that reject the null hypothesis that the residuals 
are normally distributed at the input level are rejected, unless none of the fits for a given marker
trajectory reach this threshold. Requires that R have the "nortest" module installed. DEFAULT 0.05. 
Setting to zero disables this test.

=back

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
$ENV{PATH} = "$ENV{PATH}:" . $FindBin::Bin;
use Data::Dumper;

#Get options
use Getopt::Long;
use Pod::Usage;

my $ratio_mode = 0;

my ($help, $man);
my ($s, $u, $N, $T, $stochastic);
my $replicates = 100;
my $interval = 1;
my $input_file;
my $skip_generations = 0;
my $initial_population_size;
my $single_mutation_only;
my $lillie_accept;
my $residual_error_accept;
my $gaussian_noise_stdev;
my $beneficial_mutation_distribution = 'dirac_delta';


my ($path, $lib_path, $probability_file, $verbose);
GetOptions(
	'help|?' => \$help, 'man' => \$man,
	'mutation-rate|u=s' => \$u,
	'average-selection-coefficient|s=s' => \$s,
	'beneficial-mutation-distribution|d=s' => \$beneficial_mutation_distribution,
	'population-size-after-dilution|N=s' => \$N,
	'growth-phase-generations|T=s' => \$T,
	'replicates|r=s' => \$replicates,
	'stochastic|a' => \$stochastic,
	'interval|i=s' => \$interval,
	'lilliefors-p-value|l=s' => \$lillie_accept,
	'residual-standard-error|e=s' => \$residual_error_accept,
	'probability-file|p=s' => \$probability_file,
	'skip-generations|k=s' => \$skip_generations,  ##stochastic, this is actually in transfers
	'verbose|v' => \$verbose, 
	'initial_population-size|0=s' => \$initial_population_size, ## stochastic only
	'single-mutation-only|1' => \$single_mutation_only,
	'gaussian-noise-stdev|g=s' => \$gaussian_noise_stdev,
);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
pod2usage(1) if (!defined $u || !defined $s);

$ENV{PATH} = "$ENV{PATH}:$path" if $path;

my ($s_start, $s_end, $s_interval);
if ($s =~ m/(.+)\:(.+)\:(.+)/)
{
	($s_start, $s_end, $s_interval) = ($1, $2, $3);
}
else
{
	($s_start, $s_end, $s_interval) = ($s, $s, $s);
}

my ($u_start, $u_end, $u_interval);
if ($u =~ m/(.+)\:(.+)\:(.+)/)
{
	($u_start, $u_end, $u_interval) = ($1, $2, $3);
}
else
{
	($u_start, $u_end, $u_interval) = ($u, $u, abs($u));
}

print STDERR "Mutation rates from: $u_start to $u_end by step size $u_interval\n";
print STDERR "Selection coefficients from: $s_start to $s_end by step size $s_interval\n";

for (my $s=$s_start; $s<=$s_end+$s_interval/2; $s+=$s_interval)
{
	for (my $u=$u_start; $u<=$u_end+$u_interval/2; $u+=$u_interval)
	{
		#avoid precision errors
		$s = sprintf("%.6g",$s);
		$u = sprintf("%.6g",$u);

		print STDERR "s = $s, u = $u\n";
	
		my $command;
		my $res;
		my $true_u = exp($u * log(10));

		my $prefix = (defined $stochastic) ? "stochastic" : "pop_gen";
		my $simulated_data_file_name = "$prefix\_s\=$s\_mu\=$u.tab";
		my $fit_file_name = "$prefix\_s\=$s\_mu\=$u.fit";

		if (defined $stochastic)
		{
			$command = "marker_divergence_stochastic_simulation.pl -b $beneficial_mutation_distribution -T $T -N $N -r $replicates -i $interval -s $s -u $true_u -k $skip_generations -o $simulated_data_file_name";
			$command .= " -0 $initial_population_size" if ($initial_population_size);
		}
		else
		{
			$command = "marker_divergence_pop_gen_simulation.pl -b $beneficial_mutation_distribution -T $T -N $N -r $replicates -i $interval -s $s -u $true_u -k $skip_generations -o $simulated_data_file_name";
			$command .= " -p $probability_file" if ($probability_file);
			$command .= " -1" if ($single_mutation_only);
		}
		
		#print STDERR "$command\n";
		
		#Check for existence of non-empty simulation file before running
		my $simulation_file_exists = 0;
		if (-e $simulated_data_file_name)
		{
			open TESTFILE, "<$simulated_data_file_name";
			my @test_lines = <TESTFILE>;
			$simulation_file_exists = (scalar @test_lines > 2);
		}
		
		if (!$simulation_file_exists)
		{
			print "$command\n" if ($verbose);
			$res = system $command;
			die "$command" if ($res);
		}
		
		#Check for existence of complete fit file before running
		my $fit_file_exists = 0;
		if (-e $fit_file_name)
		{
			open TESTFILE, "<$fit_file_name";
			my @test_lines = <TESTFILE>;
			$fit_file_exists = (scalar @test_lines == $replicates+1);
		}

		if (!$fit_file_exists)
		{
			$command = "marker_divergence_fit.pl -m log_ratio -i $simulated_data_file_name";
			$command .= " -l $lillie_accept" if (defined $lillie_accept);
			$command .= " -e $residual_error_accept" if (defined $residual_error_accept);
			$command .= " -g $gaussian_noise_stdev" if (defined $gaussian_noise_stdev);
			$command .= " -o $fit_file_name";
			print "$command\n" if ($verbose);
			$res = system $command;
			die "Err in command: $command" if ($res);
		}
	}
}
