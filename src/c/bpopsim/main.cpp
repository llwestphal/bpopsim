#include "cPopulation.h"

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
  ("verbose,v", value<uint16_t>(), "Verbose")
  ("lineage-tree,l", value<uint16_t>(), "Lineage Tree")
  ("seed,d", value<uint16_t>(), "Seed for random number generator")
  ("redwhite-only,w", value<char>(), "Only care about red/white lineages")
  ("log-approximation,a", value<char>(), "Less precise/faster approximation")
  ("log-approximation-value,v", value<int>(), "Precise-ness of log approximation")
  ;

/* Need to add these as options...
  'fitnesses|f=s' => \@fitnesses,
  'detailed' => \$detailed,
  'minimum-data-points|z=s' => \$minimum_data_points,
  'maximum-data-points|x=s' => \$maximum_data_points,
  'skip-generations|k=s' => \$skip_generations,
  'input_initial_w|0=s' => \$input_initial_w,
  'drop_frequency|2=s' => \$drop_frequency,
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
}

int main(int argc, char* argv[])
{
	 tree<cGenotype>::iterator_base loc;
	 uint32_t node_id;
	 uint16_t seed;
	
   //set up command line options
   variables_map cmdline_options;
   get_cmdline_options(cmdline_options, argc, argv);
   std::string output_file = cmdline_options["output-file"].as<std::string>();
  
   //Initialize Tree object 
   cLineageTree newtree;
  
   //Initialize Population object
	 cPopulation population;
  
   //Build lookup table for logs 
   //Currently it is the the 15th, I should take it as a command line option
   population.ConstructLookUpTable();
  
   //Set cli options
   population.SetParameters(cmdline_options);
   population.DisplayParameters();
	
	 //create generator and seed
	 //@agm defaults to system time seed if not specified at cli
	 const gsl_rng_type *T;
	 gsl_rng * randgen;
	 T = gsl_rng_mt19937;
	 randgen = gsl_rng_alloc (T);
	 if ( population.GetSeed() == 0 ) { seed = time (NULL) * getpid(); }
	 else { seed = population.GetSeed(); }
	 gsl_rng_set(randgen, seed);
  
	 std::vector< std::vector<cGenotypeFrequency> > frequencies;
   //Cout << Endl << population.ReturnLog(5) << Endl;
  
   for (int on_run=0; on_run < population.GetReplicates(); on_run++)
   {
      population.ClearRuns(&newtree);
		 
      uint32_t count(0);
      std::cout << "Replicate " << on_run+1;   
		 
      //Check to see if the user wants the red/white lineages included
      if (population.GetRedWhiteOnly() == 't') population.SeedSubpopulationForRedWhite(&newtree, node_id);
      if (population.GetRedWhiteOnly() == 'f') population.SeedPopulationWithOneColony(&newtree, node_id);
		  //std::cout << node_id << std::endl;
		 
      population.ResetRunStats();
		 
      while( (population.GetTransfers() < population.GetTotalTransfers()) && population.GetKeepTransferring()) 
			{
				
         // Calculate the number of divisions until the next mutation 
				 population.SetDivisionsUntilMutation(population.GetDivisionsUntilMutation() + round(gsl_ran_exponential(randgen, population.GetLambda())));
				
				 // 
				if (population.GetVerbose()) { 
					std::cout << "  New divisions before next mutation: " << population.GetDivisionsUntilMutation() << std::endl; }
         
			 	 while( population.GetDivisionsUntilMutation() > 0 && population.GetKeepTransferring()) 
				 {
					 population.CalculateDivisions();
					 
					 if( population.GetDivisionsUntilMutation() <= 0) { population.Mutate(randgen, &newtree, node_id); }
					 
					 if( population.GetPopulationSize() >= population.GetPopSizeBeforeDilution()) {
						 population.FrequenciesPerTransferPerNode(&newtree, &frequencies);
						 population.Resample(randgen); 
						 count++;
						 Cout << Endl << "Passing.... " << count << Endl;
           }
           //If the user does not want the red/white lineages included then KeepTransferring should always be set to true
           if( population.GetRedWhiteOnly() == 'f' ) population.SetKeepTransferring(true);
         }
       }
       Cout << Endl << Endl;
       population.RunSummary();
       population.PushBackRuns();
       //kptree::print_tree_bracketed(newtree);
       Cout << Endl << Endl << "Printing to screen.... " << Endl;
       population.PrintFrequenciesToScreen(frequencies);
       Cout << Endl << Endl << "Printing to screen.... " << Endl;
       population.PrintOut(output_file, frequencies);
     }
	 
   //population.PrintOut(output_file, frequencies);
}
