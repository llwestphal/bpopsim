#include "cSubpopulation.h"
#include "cPopulation.h"

#include <iostream>
#include <string>

// Boost
#include <boost/program_options.hpp>

// GSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace boost::program_options;
using namespace std;

// setup and parse configuration options:
void get_cmdline_options(variables_map &options, int argc, char* argv[]) {

  options_description cmdline_options("Allowed options");
  cmdline_options.add_options()
  ("help,h", "produce this help message")
  ("generations-per-transfer,T", value<double>(), "Number of generations per transfer")
  ("population-size-after-transfer,N", value<uint64_t>(), "Population size after transfer")
  ("output-file,o", value<string>(), "Output file")
  ;

/* Need to add these as options...

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
	'seed=s' => \$rand_seed,
	'verbose|v' => \$verbose,
*/
  store(parse_command_line(argc, argv, cmdline_options), options);
  notify(options);

  // check here for required options
  if(options.count("help")
    || !options.count("output-file")) {
    cerr << "Usage: bpopsim -o output" << endl;
    cerr << cmdline_options << endl;
    exit(0);
  }
}

int main(int argc, char* argv[])
{
  //set up command line options
  variables_map cmdline_options;
  get_cmdline_options(cmdline_options, argc, argv);
  string output_file = cmdline_options["output-file"].as<string>();
  
	//create generator and seed
	const gsl_rng_type * T;
  gsl_rng * randgen;
	gsl_rng_env_setup();
  T = gsl_rng_mt19937;
  randgen = gsl_rng_alloc (T);



	cPopulation population;
	population.SetParameters(cmdline_options);
	population.DisplayParameters();

	for (int on_run=0; on_run < population.GetReplicates(); on_run++)
	{	
		population.ClearRuns();
    std::cout << "Replicate " << on_run+1 << std::endl;    

		population.SeedSubpopulations();
		population.ResetRunStats();		

		while ( (population.GetTransfers() < population.GetTotalTransfers()) && population.GetKeepTransferring() )
		{
			
			population.SetDivisionsUntilMutation(population.GetDivisionsUntilMutation() + gsl_ran_exponential(randgen, population.GetLambda()));
			
			if (population.GetVerbose()) std::cout << "  New divisions before next mutation: " << population.GetDivisionsUntilMutation() << std::endl;
		      
			while ( (population.GetDivisionsUntilMutation() > 0) && (population.GetTransfers() < population.GetTotalTransfers()) && 					population.GetKeepTransferring()) 
			{
        population.CalculateDivisions();
				population.Mutate(randgen);
				population.Resample(randgen);			
			}
		}
    population.RunSummary();
		population.PushBackRuns();
	}
	
	population.PrintOut(output_file);
}


