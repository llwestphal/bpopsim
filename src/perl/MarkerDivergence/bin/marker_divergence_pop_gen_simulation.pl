#!/usr/bin/perl -w

###
# Pod Documentation
###

=head1 NAME

marker_divergence_pop_gen_simulation.pl

=head1 SYNOPSIS

Usage: marker_divergence_pop_gen_simulation.pl -T 6.64 -N 5E6 -s 0.08 -u 1E-8 -p pr_establishment_T_6.64_No_5E6.tab -o output.tab

Simulate marker divergence given a population size, transfer ratio, and beneficial 
mutation rate, and beneficial selection coefficient.

=head1 DESCRIPTION

=over

=item B<-o|--output-file> <file path> 

Name of output file with simulated marker divergence trajectories.

=item B<-T> <double> 

Number of generations per growth phase. Default = 6.64.

=item B<-N> <double> 

Initial population size at each growth cycle (after dilution). N sub zero.
Default = 5E6.

=item B<-u> <double>

Per generation rate of beneficial mutations. Default = 2E-7.

=item B<-s> <double>

Mean selection coefficient for beneficial mutations. Default = 0.05.

=item B<-b|--beneficial-mutation-distribution> <name of distribution>

Probability distribution of beneficial mutation selection coefficients. 
Either 'dirac_delta' or 'exponential'. Default = dirac_delta

=item B<-i> <int>

Time interval in growth cycles to print marker ratio. Default = 1.

=item B<-m> <int>

Stop simulating when the marker ratio diverges by more than this factor. Default = 100. 
When the --single-mutations-only option is active, the run ends when the ancestral population
makes up 1/m of the total population.

=item B<-r> <int>

Replicates to perform. Default = 10.

=item B<-p> <path_to_file> 

Path to a table that gives corrected probability that a beneficial mutant will survive dilution,
given its current selective advantage over the population average.
If this parameter is not supplied the classic approximation of 2s is used. 
This approximation is not very accurate for s > 0.05.

=item B<-k|--skip-generations> <double>

Skip this many initial generations before printing data. Used to account for initial outgrowth.

=item B<-z|--minimum-data-points> <int>

Require at least this many data points past time = 0 to be output before ending simulation.

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

#Improved random number generator
use Math::Random::MT::Auto qw(srand rand exponential);

#Get options
use Getopt::Long;
use Pod::Usage;
my $input_file;
my ($low_abs, $high_abs, $blank_time);
my ($help, $man);

##Model parameters
##Defaults are set to Woods experiment done under long-term conditions
my $growth_phase_generations = 6.64385619;     #T
my $pop_size_after_dilution = 5E6;        #N sub 0
my $ancestral_growth_rate = log 2;        #r
my $mutation_rate_per_division = 2E-7; #mu
my $average_mutation_s = 0.05; #s
my $single_mutation_only = 0;

my $multiplicative_selection_coefficients = 0;
#normally the growth rate is increased additively s = s1 + s2 + s3,... starting at 0
#setting this flag makes it increase as s = (1+s1)*(1+s2)*(1+s3)*... starting at 1

my $beneficial_mutation_distribution = 'dirac_delta';

##Calculated
my $effective_pop_size; # N sub e

##Output parameters
my $output_file;
my $transfer_interval_to_print = 1; #at what transfers
my $total_transfers = 1000000000;
my $max_divergence_factor = 100;
my $minimum_data_points = 10; #output must have at least this many data points (+1 for time zero)
my $maximum_data_points;
our $drop_frequency;
my $replicates = 10;

##Establishment probability file information
our $probability_file;
our @probability_table;
our $probability_precision;

my @fitnesses;
my @allele_counts;
our $model = 'kishony';
my $minimum_generations;
my $skip_generations = 0;
my $rand_seed;
my $input_initial_w;
my $verbose = 0;

##Keep track of mutations and population variation
my $detailed;

GetOptions(
	'help|?' => \$help, 'man' => \$man,
	'growth-phase-generations|T=s' => \$growth_phase_generations,
	'population-size-after-dilution|N=s' => \$pop_size_after_dilution,
	'mutation-rate-per-generation|u=s' => \$mutation_rate_per_division,
	'average_selection_coefficient|s=s' => \$average_mutation_s,
	'time-interval|i=s' => \$transfer_interval_to_print,
	'replicates|r=s' => \$replicates,
	'marker-divergence|m=s' => \$max_divergence_factor,
	'probability-file|p=s' => \$probability_file,
	'fitnesses|f=s' => \@fitnesses,
	'detailed' => \$detailed,
	'minimum-data-points|z=s' => \$minimum_data_points,
	'maximum-data-points|x=s' => \$maximum_data_points,
	'skip-generations|k=s' => \$skip_generations,
	'input_initial_w|0=s' => \$input_initial_w,
	'drop_frequency|2=s' => \$drop_frequency,
	'multiplicative' => \$multiplicative_selection_coefficients,
	'single-mutation-only|1' => \$single_mutation_only,	
	'beneficial-mutation-distribution|d=s' => \$beneficial_mutation_distribution,
	'output-file|o=s' => \$output_file,
	'seed=s' => \$rand_seed,
	'verbose|v' => \$verbose,
);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
pod2usage(1) if (!defined $growth_phase_generations || !defined $pop_size_after_dilution || !defined $mutation_rate_per_division || !defined $average_mutation_s);
pod2usage(1) if (!defined $output_file);
die "Mutation rate (-u) must be > 0" if ($mutation_rate_per_division < 0);
my $single_mutation_ancestor_fraction = 10000;

srand($rand_seed) if ($rand_seed);

my %beneficial_mutation_distribution_translation_hash = (
	'dirac_delta' => 0,
	'exponential' => 1,
);

my $beneficial_mutation_distribution_code = $beneficial_mutation_distribution_translation_hash{$beneficial_mutation_distribution};
defined $beneficial_mutation_distribution_code or die("Beneficial mutation distribution (-d) must be one of: " . (join ", ", sort keys %beneficial_mutation_distribution_translation_hash) . "\n");

$model = 'classic' if (!defined $probability_file);
$effective_pop_size = $pop_size_after_dilution * $ancestral_growth_rate * $growth_phase_generations;

if ($detailed)
{
	$drop_frequency = 0.1 if (!defined $drop_frequency);
	$drop_frequency = $drop_frequency / $effective_pop_size;
}

my $require_log_10_divergence = 0.1;

##ADD: Print out settings to output file as comments
print STDERR "Transfer interval (generations): T = $growth_phase_generations\n";
print STDERR "Population size after transfer: No = $pop_size_after_dilution\n";
print STDERR "Mutation rate per cell division: u = $mutation_rate_per_division\n";
print STDERR "Mean beneficial mutation selection coefficient: s = $average_mutation_s\n";
print STDERR "Beneficial mutation distribution: $beneficial_mutation_distribution\n";
print STDERR "Effective population size: Ne = $effective_pop_size\n";
print STDERR "Method for determining probability of establishment: $model\n";

# Effective population after conversion from batch to continuous growth model

if ($model eq 'classic')
{
	$effective_pop_size = $pop_size_after_dilution * $ancestral_growth_rate * $growth_phase_generations;
	
	# For the simplified chance of loss of a beneficial mutation to drift calculation
	# (p sub e of s = 2 * s), this product should be << 1.
	my $r_s_T = $ancestral_growth_rate * $average_mutation_s * $growth_phase_generations;
	print STDERR "Using classic approximation of the probability of a beneficial mutant surviving (2s).\n";
	print STDERR "Assumption for p sub e of
	 s calculations: (r*s*T) ($r_s_T ) << 1\n";

}
elsif ($model eq 'kishony')
{
	$effective_pop_size = $pop_size_after_dilution * $ancestral_growth_rate * $growth_phase_generations;
	
	#load the file if necessary
	open PR, "<$probability_file" or die "Could not open $probability_file";
	my @pr_lines = <PR>;
	chomp @pr_lines;
	foreach my $l (@pr_lines)
	{
		my @spl = split "\t", $l;
		push @probability_table, \@spl;
	}
	
	#first entry is zero, next gives the precision!
	$probability_precision = 1 / $probability_table[1]->[0];
	print STDERR "Probability file: $probability_file.\nPrecision: 1/$probability_precision\n";
}
else
{
	die "Unknown model.";
}

# Exponential constant for mutational Poisson process lamba is 
# mutation rate per cell division times the number of cells (N sub e)
# after the conversion to a continuous model

my $lambda = $mutation_rate_per_division * $effective_pop_size;
print STDERR "Lambda = $lambda\n" if ($verbose);

##
## Main Simulation Loop
##
my @runs;
for (my $on_run=0; $on_run < $replicates; $on_run++)
{ 
	my $num_mutations;
	my $unique_id_counter = 0;
	my %mutations_different;

	print STDERR "Replicate: " . ($on_run+1) . "\n";

	my @populations;

	my $initial_w = (scalar @fitnesses > 0) ? $fitnesses[0] : ($multiplicative_selection_coefficients ? 1 : 0);
	$initial_w = $input_initial_w if (defined $input_initial_w);
	
	push @populations, { 'n' => 0.5, 'color' => 0, 'w' => $initial_w, 'muts' => 0, id => $unique_id_counter++ };
	push @populations, { 'n' => 0.5, 'color' => 1, 'w' => $initial_w, 'muts' => 0, id => $unique_id_counter++ };
	
	#extra mutation tracking
	if ($detailed)
	{
		$mutations_different{'1-0'} = 0;
		$populations[0]->{mut_hash} = {};
		$populations[1]->{mut_hash} = {};
	}
	
	my	$on_generation = -$skip_generations;
	my	$last_printed_transfer = -$transfer_interval_to_print;
	my	$next_transfer_to_print = 0;
	
	my $total_mutations = 0;
	my $total_mutations_lost = 0;
	my $max_w = 1;
	my $avg_w = 1;
	
	MUTATION: while ($last_printed_transfer <= $total_transfers)
	{
		## Move time forward until another mutation occurs.
		#Waiting time in generations until next mutation occurs.
		my $mutation_wait_time = exponential(1/$lambda); 
		
		$num_mutations++;
		
		#Generation when mutation will happen
		my $mutation_at_generation = $on_generation + $mutation_wait_time; 
		
		## Calculate points to output between last time and time of next mutation
		## First time through the loop, a partial time interval
		## may be calculated to get back on $print_interval
		
		my $update_generations;
		TRANSFER: while ($mutation_at_generation >= $next_transfer_to_print * $growth_phase_generations)
		{
			$last_printed_transfer = $next_transfer_to_print;
			$next_transfer_to_print += $transfer_interval_to_print;
			$update_generations = $next_transfer_to_print * $growth_phase_generations - $on_generation;
			
			#print STDERR Dumper(@populations);
			my $spliced = 0;
			##drop extremely low frequencey subpopulations
			if ($drop_frequency)
			{
				for (my $i=0; $i<scalar @populations; $i++)
				{
					if ($populations[$i]->{n} < $drop_frequency)
					{				
						my $id = $populations[$i]->{id};
						splice @populations, $i, 1;
						$i--;
					
						if ($detailed)
						{
							foreach my $key (keys %mutations_different)
							{
								my ($f, $s) = split '-', $key;
								if (($f == $id) || ($s == $id))
								{
									#print STDERR "$key $f $s $id\n";
									delete $mutations_different{$key};
								}
							}
						}
						
					}
				}
			}
						
			## update subpopulation frequencies
			($on_generation, $avg_w) = update_subpopulations(\@populations, $on_generation, $update_generations);
			
			## calculate new marker ratio
			my @by_color;
			$by_color[0] = 0;
			$by_color[1] = 0;
			foreach my $p (@populations)
			{
				$by_color[$p->{color}] += $p->{n};
			}
			my $ratio = ($by_color[1] == 0) ? "INF" : $by_color[0] / $by_color[1];
									
			
			print STDERR "Transfers  = $last_printed_transfer\n" if ($verbose);
			
			## Save the current ratio and end if we are diverged
			if (!$detailed)
			{
				push @{$runs[$on_run]}, $ratio;
			}
			else
			{
				my $info;
				$info->{fixed_mutations} = 0;
				##numbers have been normalized at this point
				my $het = 0;
				my $max = 0;
				my $dominant_mut_hash = {};
				my $dominant_n = 0;
				for (my $p=0; $p < scalar @populations; $p++)
				{
					my $lineage_1 = $populations[$p];
					$info->{fixed_mutations} += $lineage_1->{muts} * $lineage_1->{n};				
					for (my $q=$p+1; $q < scalar @populations; $q++)
					{
						my $lineage_2 = $populations[$q];
						#print STDERR "$lineage_2->{id}" . "-" . $lineage_1->{id} . "\n" if (!defined $mutations_different{$lineage_2->{id} . '-' . $lineage_1->{id}});
						my $diff = $mutations_different{$lineage_2->{id} . '-' . $lineage_1->{id}};
						$het += $diff * $lineage_1->{n} * $lineage_2->{n}; 						

						if (($lineage_1->{n} > 0.01) && ($lineage_2->{n} > 0.01))
						{
							$max = $diff if ($diff > $max);
						}
					}
					if ($lineage_1->{n} > $dominant_n)
					{
						$dominant_n = $lineage_1->{n};
						$dominant_mut_hash = $lineage_1->{mut_hash};
					}
				}				
				$info->{max_mut_diff} = $max;
				$info->{average_mut_diff} = $het;
				$info->{additional_mutations} = 0;
				$info->{average_fitness} = $avg_w;
				$info->{dominant_mut_hash} = $dominant_mut_hash;
				
				push @{$runs[$on_run]}, $info;
				
				## print all populations if verbose
				print STDERR "Run: $on_run Gen: $on_generation Ratio: $ratio Avg W: $avg_w\n";
			}
				
			## logic to decide when to end	
			if ( (defined $maximum_data_points) && ($last_printed_transfer / $transfer_interval_to_print >= $maximum_data_points) )
			{
				last MUTATION;
			}	
				
			if (($ratio > $max_divergence_factor) || ($ratio < 1/$max_divergence_factor))
			{
				last MUTATION if (!defined $minimum_data_points);
				last MUTATION if ($last_printed_transfer / $transfer_interval_to_print >= $minimum_data_points);
			}
			
			if ($single_mutation_only && ($populations[0]->{n} < 1/$single_mutation_ancestor_fraction))
			{
				last MUTATION if (!defined $minimum_data_points);
				last MUTATION if ($last_printed_transfer / $transfer_interval_to_print >= $minimum_data_points);
			}
		}
		
		#my $num_populations = scalar @$pop_list_ref;
		#print STDERR "  Lineages: $num_populations  Total mutations: $total_mutations  Lost to Drift: $total_mutations_lost\n";
		#print STDERR "  Maximum Fitness: $max_w  Averge Fitness: $avg_w\n";
		
		## Update subpopulations for anything past the last printed time and the time of mutation
		$update_generations = $mutation_at_generation - $on_generation;
		($on_generation, $avg_w) = update_subpopulations(\@populations, $on_generation, $update_generations);


		print STDERR "Mutation occurs at generation: $mutation_at_generation\n" if ($verbose);
		
		##In which lineage did the mutation happen?
		my $mutated_lineage;
		my $rand = rand 1.0;
		for (my $i=0; $i<scalar @populations; $i++)
		{
			my $p = $populations[$i];
			$rand -= $p->{n};
			if ($rand <= 0)
			{
				$mutated_lineage = $i;
				print STDERR "Mutated lineage $i\n" if ($verbose);
				last;
			}
		}
		die "Undefined lineage mutated?" if (!defined $mutated_lineage);
		my $ancestor = $populations[$mutated_lineage];

		next MUTATION if ($single_mutation_only && $ancestor->{'muts'} > 0);

		##What is the new lineage's fitness? (For now Dirac delta function for fitness distribution)
		my $new_w;
		if (scalar @fitnesses == 0)
		{
			if ($beneficial_mutation_distribution_code == 0) #dirac delta
			{
				$new_w = ($multiplicative_selection_coefficients) ? $ancestor->{w} * (1+$average_mutation_s) : $ancestor->{w} + $average_mutation_s;
			}
			else #exponential
			{
				my $this_mutation_s = exponential($average_mutation_s);
				$new_w = ($multiplicative_selection_coefficients) ? $ancestor->{w} * (1+$this_mutation_s) : $ancestor->{w} + $this_mutation_s;
			}
		}
		else
		{
			$new_w = $ancestor->{w};
			if ($ancestor->{muts}+1 < scalar @fitnesses)
			{
				$new_w = $fitnesses[$ancestor->{'muts'}+1];
			}
		}
		$max_w = $new_w if ($new_w > $max_w);
		
		##Does the new lineage survive drift?
		my $diff_w = $new_w - $avg_w;
		my $probability_of_surviving_drift = probability_of_beneficial_mutant_surviving($diff_w);
		print STDERR "Probability of fitness $new_w advantage $diff_w surviving drift = $probability_of_surviving_drift\n" if ($verbose);
		
		$rand = rand 1.0;
		$total_mutations++;
		if ($rand > $probability_of_surviving_drift)
		{
			$total_mutations_lost++;
			next MUTATION;
		}
		
		##Create and add the new lineage
		my $new_mut_hash = {};
		if ($detailed)
		{
			$new_mut_hash->{$unique_id_counter}++;
			foreach my $key (%{$ancestor->{mut_hash}})
			{
				$new_mut_hash->{$key}++;
			}
		}
				
		my $new_lineage = 	{ 
								'n'=> 1/($effective_pop_size*$probability_of_surviving_drift), 
								'color' => $ancestor->{'color'}, 
								'w' => $new_w, 
								'muts' => $ancestor->{'muts'} + 1,
								id => $unique_id_counter++,	
								mut_hash => $new_mut_hash,
							};
				
		#print STDERR "id = $unique_id_counter\n";		
									
		## calculate all differences vs currently living ids
		if ($detailed)
		{			
			foreach my $p (@populations)
			{
				my $diff = 0;
				foreach my $key (keys %$new_mut_hash)
				{
					$diff++ if (!defined $p->{mut_hash}->{$key});
				}
				
				$mutations_different{$new_lineage->{id} . '-' . $p->{id}} = $diff;
				#print STDERR ">>$new_lineage->{id}-$p->{id}\n";
			}
		}
		
		## $probability_of_surviving_drift takes care up super-exponential growth of ones lucky enough to make it through sampling
		push @populations, $new_lineage;

		
		if ($verbose)
		{
			print STDERR "ADDED! $new_lineage->{color}:$new_lineage->{w}:$new_lineage->{n}\n" 
		}
	}
	
	#Check to see if the current run diverged by enough to fit
	if ($single_mutation_only && ( abs(log($runs[$on_run]->[$#{$runs[$on_run]}])/log(10)) < $require_log_10_divergence ) )
	{
		print STDERR "Encountered marker divergence curve that did not diverge by log10($require_log_10_divergence).\n";
		@{$runs[$on_run]} = ();
		$on_run--;
	}

	print STDERR "Average mutations per generation = " . ($num_mutations/($on_generation*$effective_pop_size)) . "\n"  if ($verbose);

}


#Print out all results
open OUT, ">$output_file";
if (!$detailed)
{

	my $line = "transfer";
	for (my $on_run=0; $on_run < $replicates; $on_run++)
	{
		$line .= "\t$on_run";
	}

	my $still_going = 1;
	my $transfer = 0;
	while ($still_going)
	{
		print OUT "$line\n";
		
		$line = $transfer * $transfer_interval_to_print;
		$still_going = 0;
		for (my $on_run=0; $on_run < $replicates; $on_run++)
		{
			$line .= "\t";
			if (defined $runs[$on_run]->[$transfer])
			{
				$line .= sprintf("%.4f",log($runs[$on_run]->[$transfer]));
				$still_going = 1;
			}
		}
		$transfer++;
	}
}
else
{

	my $line = "transfer";
	for (my $on_run=0; $on_run < $replicates; $on_run++)
	{
		$line .= "\tmutations-$on_run";
	}
	for (my $on_run=0; $on_run < $replicates; $on_run++)
	{
		$line .= "\tfitness-$on_run";
	}
	
	for (my $on_run=0; $on_run < $replicates; $on_run++)
	{
		$line .= "\tavg_diff_mutations-$on_run";
	}
	for (my $on_run=0; $on_run < $replicates; $on_run++)
	{
		$line .= "\tmax_diff_mutations-$on_run";
	}
	
	my $still_going = 1;
	my $transfer = 0;
	while ($still_going)
	{
		print OUT "$line\n";
		
		$line = $transfer * $transfer_interval_to_print;
		$still_going = 0;
		for (my $on_run=0; $on_run < $replicates; $on_run++)
		{
			$line .= "\t";
			if (defined $runs[$on_run]->[$transfer])
			{
				my $info = $runs[$on_run]->[$transfer];
				$line .= join "\t", ($info->{fixed_mutations});
				$still_going = 1;
			}
		}
		for (my $on_run=0; $on_run < $replicates; $on_run++)
		{
			$line .= "\t";
			if (defined $runs[$on_run]->[$transfer])
			{
				my $info = $runs[$on_run]->[$transfer];
				$line .= join "\t", ($info->{average_fitness});
				$still_going = 1;
			}
		}
		
		for (my $on_run=0; $on_run < $replicates; $on_run++)
		{
			$line .= "\t";
			if (defined $runs[$on_run]->[$transfer])
			{
				my $info = $runs[$on_run]->[$transfer];
				$line .= join "\t", ($info->{average_mut_diff});
				$still_going = 1;
			}
		}
		
		for (my $on_run=0; $on_run < $replicates; $on_run++)
		{
			$line .= "\t";
			if (defined $runs[$on_run]->[$transfer])
			{
				my $info = $runs[$on_run]->[$transfer];
				$line .= join "\t", ($info->{max_mut_diff});
				$still_going = 1;
			}
		}
		$transfer++;
	}
	
	#print the final dominant for each run...
	open FD, ">final_dominant.dat";
	
	print Dumper($runs[0]);
	for (my $on_run=0; $on_run < $replicates; $on_run++)
	{
		print FD +join(",", keys %{$runs[$on_run]->[-1]->{dominant_mut_hash}}) . "\n";
	}
}


#####
# probability_of_beneficial_mutant_surviving
######

sub probability_of_beneficial_mutant_surviving
{
	my ($net_s) = @_;

	#not really a beneficial mutation, zero probability of surviving
	return 0 if ($net_s <= 0);

	#classic approximation
	return 2 * $net_s if (!$probability_file);

	my $i = int($net_s*$probability_precision);
		
	if ($i > scalar @probability_table)
	{
		my $i = $#probability_table;
		print STDERR "Warning: beneficial mutation with fitness advantage $net_s exceeds probability table size."; 
		print STDERR "Using last probability entry, index $i, probability $probability_table[$i]\n";
	}
	
	return $probability_table[$i]->[1];
}

#####
# update_subpopulations
#
# Run time forward by the specified number of generations.
#####

sub update_subpopulations
{
	my ($pop_list_ref, $on_generation, $update_generations) = @_;
			
	## update frequencies of each subpopulation
	my $renormalize_n = 0;
	foreach my $p (@$pop_list_ref)
	{
		$p->{n} *= exp($ancestral_growth_rate * $p->{w} * $update_generations);
		$renormalize_n += $p->{n};
	}
	
	##renormalize to a total population size of 1, update the average fitness
	my $avg_w = 0;
	foreach my $p (@$pop_list_ref)
	{
		$p->{n} /= $renormalize_n; 
		$avg_w += $p->{n} * $p->{w};
	}	
	
	$on_generation += $update_generations;
	
	###
	print STDERR "== Updating time by $update_generations generations, now at $on_generation\n" if ($verbose);
	if ($verbose)
	{
		foreach my $p (@$pop_list_ref)
		{
	    	print STDERR "\t$p->{color}:$p->{w}:$p->{n}\n";
		}
	}
	
	return ($on_generation, $avg_w);
}