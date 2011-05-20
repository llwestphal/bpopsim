#include "common.h"
#include "cPopulation.h"

#include <iostream>

// Global variable for keeping track of verbosity
bool g_verbose = false;
bool g_ro_only = false;

// setup and parse configuration options:
void get_cmdline_options(variables_map &options, uint16_t argc, char* argv[]) {

  options_description cmdline_options("Allowed options");
  cmdline_options.add_options()
  ("help,h", "produce this help message")
  ("generations-per-transfer,T", value<double>(), "Generations per transfer")
  ("population-size-after-transfer,N", value<uint32_t>(), "Population size after transfer")
  ("number-of-transfers,n", value<uint16_t>(), "Max number of transfer to replicate")
  ("output-folder,o", value<std::string>(), "Output folder")
  ("mutation-rate-per-division,u", value<double>(), "Mutation rate per division")
  ("initial-population-size,i", value<uint32_t>(), "Initial Population Size")
  ("replicates,r", value<uint16_t>(), "Replicates")
  ("marker-divergence,m", value<uint16_t>(), "Max divergence factor")
  ("type-of-mutations,f", value<char>(), "Type of mutations")
  ("verbose,v", "Verbose")
  ("lineage-tree,l", value<uint16_t>(), "Lineage Tree")
  ("seed,d", value<uint16_t>(), "Seed for random number generator")
  ("red-white,k", "Only care about red/white lineages. For marker divergence.")
  ("transfer-interval-to-print,t", value<uint16_t>(), "Red/White Printing intervals.")
  ("imv,x", value< std::vector<double> >(), "Initial Mutational Values.")
  ("coarse-graining,c", value<uint16_t>(), "Amount to coarse-grain output.")
  ;

  store(parse_command_line(argc, argv, cmdline_options), options);
  notify(options);

  // check here for required options
  if(options.count("help")
      || !options.count("output-folder")) {
      std::cerr << "Usage: GSL_RNG_SEED=123 bpopsim -o output" << std::endl << std::endl;
      std::cerr << cmdline_options << std::endl;
      exit(0);
  }
  
  if ( options.count("verbose") ) g_verbose = true;
  if ( options.count("red-white") ) g_ro_only = true;
}

int main(int argc, char* argv[])
{
	
  //set up command line options
  variables_map cmdline_options;
  get_cmdline_options(cmdline_options, argc, argv);
  
  std::string output_folder = cmdline_options["output-folder"].as<std::string>();
  uint16_t num_replicates = cmdline_options["replicates"].as<uint16_t>();
  
  uint16_t transfer_interval_to_print(1);
  if ( cmdline_options.count("transfer-interval-to-print") ) {
      transfer_interval_to_print = cmdline_options["transfer-interval-to-print"].as<uint16_t>();
  }
  
  std::vector< std::vector<double> > red_white_ratios;
	
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
  if ( cmdline_options.count("seed") ) {
    seed = cmdline_options["seed"].as<uint16_t>();
  } else {
    seed = time(NULL) * getpid();
  }
  
  gsl_rng_set(randgen, seed);
  
  if (g_ro_only) std::cout << std::endl << "You chose to use red and white only." << std::endl;
  
  //std::vector< std::vector<uint16_t> > number_unique_genotypes_in_all_replicates;
  
  for (uint32_t on_run=0; on_run < num_replicates; on_run++)
  {
    std::vector<double> current_ro_ratio;
    
    //Initialize Population object
    cPopulation population;
    
    //Build lookup table for logs 
    //Currently it is the the 15th, I should take it as a command line option
    //population.ConstructLookUpTable();
    
    //Set cli options
    population.SetParameters(cmdline_options);
    population.DisplayParameters();
    
    tree<cGenotype>::iterator_base loc;
    std::vector< std::vector<cGenotypeFrequency> > frequencies, subpops;
    
    population.SetRNG(randgen);
    
    uint32_t count(0);
    std::cout << "Replicate " << on_run << std::endl;   
		 
    //Initialize the population
    if (g_ro_only) {
      population.SeedSubpopulationForRedWhite(); 
    } else {
      population.SeedPopulationWithOneColony();
    }
    
    //Get an initial time points
    population.CalculateAverageFitness();
    population.FrequenciesPerTransferPerNode(&frequencies);
    
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
          
          if( on_run == 0 ) 
            population.CalculateAverageFitness();
          
          population.Resample();
          population.CullPopulations();
          
          count++;
          std::cout << std::endl << "Passing.... " << count << std::endl;
          
          if ( population.GetTransfers() %  transfer_interval_to_print == 0 ) {  
            current_ro_ratio.push_back(population.GetRatio());
            
            if (g_verbose == 1)
            {
              std::cout << "Transfer " << population.GetTransfers() << " : " << 
              "=>" << population.GetPopulationSize() << "  R/W Ratio: " << population.GetRatio() << std::endl;  
              std::cout << "Total mutations: " << population.GetTotalMutations() << " Maximum Fitness: " << population.GetMaxW() << std::endl;
              std::cout << "Size = " << current_ro_ratio.size() << std::endl;
            }
          }
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
    
    if (g_ro_only) {
      red_white_ratios.push_back(current_ro_ratio);
    }
    else {
      //number_unique_genotypes_in_all_replicates.push_back(population.CurrentUniqueGenotypes());
      
      //std::cout << std::endl << std::endl << "Printing to screen.... " << std::endl;
      //population.PrintFrequenciesToScreen(output_folder, &frequencies);
    
      std::cout << std::endl << std::endl << "Printing max difference of relevant mutations.... " << std::endl;
      population.CalculateSimilarity(output_folder, &frequencies);
  
      std::cout << std::endl << std::endl << "Printing time to sweep.... " << std::endl;
      population.TimeToSweep(output_folder, &frequencies);
      
      if( on_run == 0 ) {
        std::cout << std::endl << std::endl << "Printing to file.... " << std::endl;
        population.PrintOut(output_folder, &frequencies);
        
        std::cout << std::endl << std::endl << "Printing average fitness.... " << std::endl;
        population.PrintFitness(output_folder);
    
        std::cout << std::endl << "Generating Muller Matrix.... " << std::endl;
        std::vector< std::vector<int> > muller_matrix;
        population.DrawMullerMatrix(output_folder, muller_matrix, &frequencies);
      }
    }
  }
  /*
  if (g_ro_only) {
    //Initialize Population object
    cPopulation population;
    std::cout << std::endl << "Printing to r/w ratio file.... " << std::endl;
    population.PrintOut_RedWhiteOnly(output_folder, &red_white_ratios, transfer_interval_to_print);
  }
  else {
    //Initialize Population object
    cPopulation population;
    std::cout << std::endl << "Printing unique genotypes to file.... " << std::endl;
    population.PrintUniqueGenotypes(output_folder, &number_unique_genotypes_in_all_replicates);
  }*/
}
