#include "common.h"
#include "cPopulation.h"

#include <iostream>

// Global variable for keeping track of verbosity
bool g_verbose = false;

// setup and parse configuration options:
void get_cmdline_options(variables_map &options, uint16_t argc, char* argv[]) {

  options_description cmdline_options("Allowed options");
  cmdline_options.add_options()
  ("help,h", "produce this help message")
  ("generations-per-transfer,T", value<double>(), "Generations per transfer")
  ("population-size-after-transfer,N", value<uint32_t>(), "Population size after transfer")
  ("number-of-transfers,n", value<uint16_t>(), "Max number of transfer to replicate")
  ("output-file,o", value<std::string>(), "Output file")
  ("mutation-rate-per-division,u", value<double>(), "Mutation rate per division")
  ("average-selection-coefficient,s", value<double>(), "Average selection coefficient")
  ("time-interval,i", value<uint16_t>(), "Time interval")
  ("replicates,r", value<uint16_t>(), "Replicates")
  ("marker-divergence,m", value<uint16_t>(), "Max divergence factor")
  ("type-of-mutations,f", value<char>(), "Type of mutations")
  ("verbose,v", "Verbose")
  ("lineage-tree,l", value<uint16_t>(), "Lineage Tree")
  ("seed,d", value<uint16_t>(), "Seed for random number generator")
  ("red-white,m", value<bool>(), "Only care about red/white lineages. For marker divergence.")
  ("log-approximation,a", value<char>(), "Less precise/faster approximation")
  ("log-approximation-value", value<int>(), "Precise-ness of log approximation")
  ;

/* Need to add these as options...
  'fitnesses|f=s' => \@fitnesses,
  'detailed' => \$detailed,
  'minimum-data-points|z=s' => \$minimum_data_points,
  'maximum-data-points|x=s' => \$maximum_data_points,
  'skip-generations|k=s' => \$skip_generations,
  'input_initial_w|0=s' => \$input_initial_w,
  'multiplicative' => \$multiplicative_selection_coefficients
*/

  store(parse_command_line(argc, argv, cmdline_options), options);
  notify(options);

  // check here for required options
  if(options.count("help")
      || !options.count("output-file")) {
      std::cerr << "Usage: GSL_RNG_SEED=123 bpopsim -o output" << std::endl << std::endl;
      std::cerr << cmdline_options << std::endl;
      exit(0);
  }
  
  if (options.count("verbose")) g_verbose = true;
}

int main(int argc, char* argv[])
{
  tree<cGenotype>::iterator_base loc;
	
  //set up command line options
  variables_map cmdline_options;
  get_cmdline_options(cmdline_options, argc, argv);
  std::string output_file = cmdline_options["output-file"].as<std::string>();
  
  //Initialize Population object
  cPopulation population;
  
  //Build lookup table for logs 
  //Currently it is the the 15th, I should take it as a command line option
  population.ConstructLookUpTable();
  
  //Set cli options
  population.SetParameters(cmdline_options);
  population.DisplayParameters();
	
  //Create Random Number generator and Seed
  //@agm Program defaults to system time seed if not specified at cli
  //     Changed generator to taus2 because it's a little faster and still "simulation quality"
  
  // @JEB: use one main RNG (so we don't reset for each replicate)
  // Let cPopulation store a copy of it.
  const gsl_rng_type *T;
  gsl_rng * randgen;
  gsl_rng_env_setup();
  T = gsl_rng_taus2;
  randgen = gsl_rng_alloc(T);
  
  uint16_t seed = 0;
  if (cmdline_options.count("seed")) {
    seed = cmdline_options["seed"].as<uint16_t>();
  } else {
    seed = time(NULL) * getpid();
  }  
  gsl_rng_set(randgen, seed);
  population.SetRNG(randgen);
  
  std::vector< std::vector<cGenotypeFrequency> > frequencies, subpops;
  
  for (int on_run=1; on_run <= population.GetReplicates(); on_run++)
  {
    // Re-initialize the population for a new run 
    // (should really clean up at end of loop, not beginning @jeb)
    population.ClearRuns();
    population.ResetRunStats();
    
    uint32_t count(0);
    std::cout << "Replicate " << on_run << std::endl;   
		 
    //Initialize the population
    if (population.GetRedWhiteOnly()) {
      population.SeedSubpopulationForRedWhite(); 
    } else {
      population.SeedPopulationWithOneColony();
    }
    
    // Print the initial tree
    //if (g_verbose) population.PrintTree();
    
    //std::cout << node_id << std::endl;
		 
    while( (population.GetTransfers() < population.GetTotalTransfers()) && population.GetKeepTransferring() ) {
				
      // Calculate the number of divisions until the next mutation 
      population.SetDivisionsUntilMutation(population.GetDivisionsUntilMutation() + round(gsl_ran_exponential(randgen, population.GetLambda())));
				
      if (g_verbose) { 
        std::cout << "  New divisions before next mutation: " << population.GetDivisionsUntilMutation() << std::endl; 
      }
         
      while( population.GetDivisionsUntilMutation() > 0 && population.GetKeepTransferring() ) 
      {
        population.CalculateDivisions();
					 
        if( population.GetDivisionsUntilMutation() <= 0) { 
          population.Mutate(); 
        }
					 
        if( population.GetPopulationSize() >= population.GetPopSizeBeforeDilution()) {
          population.FrequenciesPerTransferPerNode(&frequencies);
          population.Resample(); 
          count++;
          std::cout << std::endl << "Passing.... " << count << std::endl;
        }
        
        //if (g_verbose) population.PrintTree();
      }
    }
    std::cout << std::endl << std::endl;
    population.RunSummary();
    population.PushBackRuns();

    if (g_verbose) {
      population.PrintTree();
    }
     
    //Cout << Endl << Endl << "Printing to screen.... " << Endl;
    //population.PrintFrequenciesToScreen(&frequencies);
    //std::cout << std::endl << std::endl << "Printing to file.... " << std::endl;
    //population.PrintOut(output_file, &frequencies);
    std::cout << std::endl << std::endl << "Printing max difference of relevant mutations.... " << std::endl;
    population.CalculateSimilarity(&frequencies);
    std::cout << std::endl << "Generating Muller Matrix.... " << std::endl;
    std::vector< std::vector<int> > muller_matrix;
    population.DrawMullerMatrix(output_file, muller_matrix, &frequencies);
  }
	 
   //population.PrintOut(output_file, frequencies);
}
