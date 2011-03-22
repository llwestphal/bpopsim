#!/usr/bin/perl -w

###
# Pod Documentation
###

=head1 NAME

marker_divergence_stochastic_simulation.pl

=head1 SYNOPSIS

Usage: marker_divergence_stochastic_simulation.pl

Experimental. No documentation yet.

=head1 DESCRIPTION

=over

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
use POSIX;

#Improved random number generator
use Math::Random::MT::Auto qw(rand exponential);

#Get options
use Getopt::Long;
use Pod::Usage;
my $input_file;
my ($low_abs, $high_abs, $blank_time);
my ($help, $man);


##Model parameters
##Defaults are set to Woods experiment done under long-term conditions
my $pop_size_after_dilution = 5E6;        #N sub 0
my $mutation_rate_per_division = 1E-8; #mu
my $average_mutation_s = 0.05; #s
my $beneficial_mutation_distribution = 'dirac_delta';
my $growth_phase_generations = 6.64;

##Output parameters
my $total_transfers = 200000;
my $max_divergence_factor = 100;
my $initial_population_size = 2;
my $verbose = 0;

my $replicates = 10;
my $transfer_interval_to_print = 1;
my $skip_transfers = 0;
my $minimum_printed = 8;

GetOptions(
	'help|?' => \$help, 'man' => \$man,
	'population-size-after-dilution|N=s' => \$pop_size_after_dilution,
	'mutation-rate-per-division|u=s' => \$mutation_rate_per_division,
	'average_selection_coefficient|s=s' => \$average_mutation_s,
	'beneficial-mutation-distribution|d=s' => \$beneficial_mutation_distribution,
	'transfer-interval|i=s' => \$transfer_interval_to_print,
	'growth-phase-generations|T=s' => \$growth_phase_generations,
	'replicates|r=s' => \$replicates,
	'initial_population_size|0' => $initial_population_size,
	'verbose|v' => \$verbose,
	'skip_transfers|k' => \$skip_transfers,
);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
pod2usage(1) if (!defined $pop_size_after_dilution || !defined $mutation_rate_per_division || !defined $average_mutation_s);
die "Mutation rate (-u) must be > 0" if ($mutation_rate_per_division < 0);

my %beneficial_mutation_distribution_translation_hash = (
	'dirac_delta' => 0,
	'exponential' => 1,
);

my $beneficial_mutation_distribution_code = $beneficial_mutation_distribution_translation_hash{$beneficial_mutation_distribution};
defined $beneficial_mutation_distribution_code or die("Beneficial mutation distribution (-d) must be one of: " . (join ", ", sort keys %beneficial_mutation_distribution_translation_hash) . "\n");

my $dilution_factor = exp(log(2) * $growth_phase_generations);
my $pop_size_before_dilution = $pop_size_after_dilution * $dilution_factor;

print STDERR "u = $mutation_rate_per_division\n" if ($verbose);
print STDERR "s = $average_mutation_s\n" if ($verbose);
print STDERR "N = $pop_size_after_dilution\n" if ($verbose);
print STDERR "T = $growth_phase_generations\n" if ($verbose);
print STDERR "dil = $dilution_factor\n" if ($verbose);

##Print out settings to output file

#####ADD


# Exponential constant for mutational Poisson process lamba is 
# mutation rate per cell division times the number of cells (N sub e)
# after the conversion to a continuous model

my $lambda = $mutation_rate_per_division;


#my $lambda = $mutation_rate_per_division * 2;

	## Update: no factor of two is needed. we typically count in terms of new genomes.
	## Note that the factor of two is necessary to account for the fact that both
	## of the new cells produced by a division have this chance of happening.
	## Correction is necessary since we count in terms of new cells only.

print STDERR "Lambda = $lambda\n" if ($verbose);

my @runs;
for (my $on_run=0; $on_run < $replicates; $on_run++)
{ 
	print STDERR "Replicate " . ($on_run+1) . "\n";

	##
	## Main simulation Loop
	##

	my @populations;

	push @populations, { 'n' => $initial_population_size/2, 'color' => 0, 'w' => 1 };
	push @populations, { 'n' => $initial_population_size/2, 'color' => 1, 'w' => 1 };
	my $max_w = 1;
	my $total_pop_size = $initial_population_size;
	my $total_mutations = 0;
	my $total_mutations_lost = 0;

	#only include the initial value if we are not skipping any transfers
	if ($skip_transfers == -1)
	{
		#print "transfers\tratio\n";
		push @{$runs[$on_run]}, 1;
	}
	
	my $transfers = -$skip_transfers;
	my $divisions_until_mutation = 0;
	TRANSFER: while ($transfers < $total_transfers)
	{
		## Move time forward until another mutation occurs.
		$divisions_until_mutation += exponential(1/$lambda);
		print STDERR "Divisions before next mutation: $divisions_until_mutation\n" if ($verbose);
				
		## Calculate points to output between last time and current time
		## First time through the loop, a partial time interval
		## may be calculated to get back on $print_interval

		## Move forward by a large chunk of time which assumes
		## all populations have the ancestor fitness (=1.0)
		## This will, at worst, overestimate how long.
		## We can then backtrack by single divisions to find the division where the mutation occurs
		
		my @divided_lineages = ();
		while ($divisions_until_mutation > 0)
		{	
			##Which lineage is replicating next?
			##Keep track of time in a perfectly integrated fashion
			
			my $desired_divisions = $divisions_until_mutation;
			$desired_divisions = $pop_size_before_dilution - $total_pop_size if ($desired_divisions + $total_pop_size > $pop_size_before_dilution);
			
			## Note: we underestimate by a few divisions so that we can step forward by single division increments
			## as we get close to the one where the mutation happened (or right before a transfer).			
			$desired_divisions = 1 if ($desired_divisions < 1);
			
			print STDERR "Total pop size: $total_pop_size\n" if ($verbose);
			print STDERR "Desired divisions: $desired_divisions\n" if ($verbose);
			
			## How much time would we like to pass to achieve the desired number of divisions?
			## (assuming the entire population has the maximum fitness)
			my $est_update_time = log( ($desired_divisions+$total_pop_size) / $total_pop_size ) / ($max_w * log(2));
		
			##What is the minimum time required to get a single division?
			my $min_update_time;
			for (my $i=0; $i< scalar @populations; $i++) 	
			{
				my $p = $populations[$i];
			
				#what is the time to get to the next whole number of cells?
				my $current = $p->{n};
				my $whole_cells = POSIX::floor($p->{n}) + 1;
				
				# WC = N * exp(growth_rate * t) 
				my $time = log($whole_cells / $current) / ($p->{w});
				
				if ( !defined $min_update_time || ($time <= $min_update_time) )
				{
					@divided_lineages = () if ( (defined $min_update_time) && ($time < $min_update_time) );
					$min_update_time = $time;
					push @divided_lineages, $i;  #a list, because there can be ties
				}
			}
			
			my $update_time = ($min_update_time > $est_update_time) ? $min_update_time : $est_update_time;
					
			print STDERR "Update time: $update_time\n" if ($verbose);
					
			##Now update all lineages by the time that actually passed
			my $new_pop_size = 0;
			for (my $i=0; $i< scalar @populations; $i++) 	
			{
				my $p = $populations[$i];
				$p->{n} = $p->{n} * exp(log(2) * $update_time * $p->{w});
				$new_pop_size += POSIX::floor($p->{n});
			}				
			
			my $completed_divisions = $new_pop_size - $total_pop_size;
			$divisions_until_mutation -= $completed_divisions;
			$total_pop_size = $new_pop_size;
			
			#print "$update_time $new_pop_size $total_pop_size $completed_divisions\n";	
						
			
			#print "Divisions until mutation: $divisions_until_mutation\n";
					
			if ($divisions_until_mutation <= 0)
			{
				$total_mutations++;
				#print "Mutating!\n";
			
				##Mutation happened in the one that just divided.
				##Break ties randomly here.
				die "Undefined lineage divided?" if (scalar @divided_lineages == 0);
				my $ancestor = $populations[$divided_lineages[int rand scalar @divided_lineages]];
				
				##What is the new lineage's fitness? (For now Dirac delta function for fitness distribution)
				my $new_w;
	
				if ($beneficial_mutation_distribution_code == 0) #dirac delta
				{
					$new_w = $ancestor->{w} + $average_mutation_s;
				}
				else #exponential
				{
					my $this_mutation_s = exponential($average_mutation_s);
					$new_w = $ancestor->{w} + $this_mutation_s;
				}
				
				##Create and add the new lineage
				my $new_lineage = { 'n'=> 1, 'color' => $ancestor->{'color'}, 'w'=> $new_w};
				push @populations, $new_lineage;
				
				##Update maximum fitness
				$max_w = $new_w if ($new_w > $max_w);
				
				##One ancestor organism was converted to the new lineage
				$ancestor->{n}--;
			}
									
			##When it is time for a transfer, resample population
			if ($total_pop_size >= $pop_size_before_dilution)
			{
				my @by_color;
				
				my $new_pop_size = 0;
				POPULATION: for (my $i=0; $i< scalar @populations; $i++) 	
				{
					my $p = $populations[$i];
					
					## Low N, treat stochastically
					if ($p->{n} <= 1000)
					{
						my $new_n = 0;
						for (my $n=0; $n<$p->{n}; $n++)
						{
							$new_n++ if (rand() < 1/$dilution_factor);
						}
						$p->{n} = $new_n;
						
						##Extinction!!!
						if ($new_n == 0)
						{
							$total_mutations_lost++;
							splice @populations, $i, 1;
							$i--;
							next POPULATION;
						}
					}
					## High N, just divide
					else
					{
						$p->{n} = $p->{n} / $dilution_factor;
					}
					
					
					$new_pop_size += POSIX::floor($p->{n});
					$by_color[$p->{color}] += $p->{n};
				}
				##One color was lost, bail	
				last TRANSFER if (!$by_color[0] || !$by_color[1]);
				my $ratio = $by_color[0] / $by_color[1];
				$transfers++;
			
				if ( ($transfers >= 0) && ($transfers % $transfer_interval_to_print == 0) )
				{
					print STDERR "Transfer $transfers: $total_pop_size => $new_pop_size [R/W Ratio $ratio]\n" if ($verbose);
					print STDERR "  Total mutations: $total_mutations  Lost to Drift: $total_mutations_lost  Maximum Fitness: $max_w\n" if ($verbose);
					
					push @{$runs[$on_run]}, $ratio;
				}
				
				$total_pop_size = $new_pop_size;
				last TRANSFER if ((scalar @{$runs[$on_run]} >= $minimum_printed) && (($ratio > $max_divergence_factor) || ($ratio < 1/$max_divergence_factor)));

			}	
		}
	}
}


##Print everything out

my $line = "transfer";
for (my $on_run=0; $on_run < $replicates; $on_run++)
{
	$line .= "\t$on_run";
}

my $still_going = 1;
my $transfer = 0;
while ($still_going)
{
	print "$line\n";
	
	$line = $transfer * $transfer_interval_to_print;
	$still_going = 0;
	for (my $on_run=0; $on_run < $replicates; $on_run++)
	{
		$line .= "\t";
		if (defined $runs[$on_run]->[$transfer])
		{
			$line .= $runs[$on_run]->[$transfer];
			$still_going = 1;
		}
	}
	
	$transfer++;
}