#include "cPopulation.h"
#include "tree_util.hh"
#include "tree.hh"

// setup and parse configuration options:
void get_cmdline_options(variables_map &options, int argc, char* argv[]) {

  options_description cmdline_options("Allowed options");
  cmdline_options.add_options()
  ("help,h", "produce this help message")
  ("generations-per-transfer,T", value<double>(), "Generations per transfer")
  ("population-size-after-transfer,N", value<uint64_t>(), "Population size after transfer")
  ("output-file,o", value<std::string>(), "Output file")
  ("mutation-rate-per-division,u", value<double>(), "Mutation rate per division")
  ("average-selection-coefficient,s", value<double>(), "Average selection coefficient")
  ("time-interval,i", value<int>(), "Time interval")
  ("replicates,r", value<int>(), "Replicates")
  ("marker-divergence,m", value<int>(), "Max divergence factor")
  ("type-of-mutations,f", value<char>(), "Type of mutations")
  ("verbose,v", value<int>(), "Verbose")
  ("lineage-tree,l", value<int>(), "Lineage Tree")
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

void print_tree(const tree<cGenotype>& tr, tree<cGenotype>::pre_order_iterator it, tree<cGenotype>::pre_order_iterator end)
{
	if(!tr.is_valid(it)) return;
	int rootdepth=tr.depth(it);
	std::cout << "-----" << std::endl;
	while(it!=end) {
		for(int i=0; i<tr.depth(it)-rootdepth; ++i) 
			std::cout << "  ";
		  std::cout << (*it).unique_node_id << "," << (*it).fitness << std::endl << std::flush;
		++it;
	}
	std::cout << "-----" << std::endl;
}

int main(int argc, char* argv[])
{
	 tree<cGenotype>::iterator_base loc;
   //set up command line options
   variables_map cmdline_options;
   get_cmdline_options(cmdline_options, argc, argv);
   std::string output_file = cmdline_options["output-file"].as<std::string>();
  
   //create generator and seed
   const gsl_rng_type * T;
   gsl_rng * randgen;
   gsl_rng_env_setup();
   T = gsl_rng_mt19937;
   randgen = gsl_rng_alloc (T);
	
   cPopulation population;
   population.SetParameters(cmdline_options);
   population.DisplayParameters();
	
	 cLineageTree newtree;
	
   for (int on_run=0; on_run < population.GetReplicates(); on_run++)
   {
		  int node_id = 3;
		  population.ClearRuns(newtree);
		 
      std::cout << "Replicate " << on_run+1 << std::endl;   
		 
      //population.SeedSubpopulations(tree);
		 population.NewSeedSubpopulation(newtree);
		 
		 population.ResetRunStats();
		 
      while ( (population.GetTransfers() < population.GetTotalTransfers()) && population.GetKeepTransferring() ){
				
         population.SetDivisionsUntilMutation(population.GetDivisionsUntilMutation() + gsl_ran_exponential(randgen, population.GetLambda()));
				 if (population.GetVerbose()) std::cout << "  New divisions before next mutation: " << population.GetDivisionsUntilMutation() << std::endl;
         while ( (population.GetDivisionsUntilMutation() > 0) && (population.GetTransfers() < population.GetTotalTransfers()) && population.GetKeepTransferring()) {
					 population.CalculateDivisions();
					 if (population.GetDivisionsUntilMutation() <= 0) { node_id++; population.NewMutate(randgen, newtree, node_id); };
					 if (population.GetTotalPopSize() >= population.GetPopSizeBeforeDilution()) { population.Resample(randgen); };  
					 
         }
				
      }
		  //print_tree(newtree, newtree.begin(), newtree.end());
      population.RunSummary();
      population.PushBackRuns();
		  kptree::print_tree_bracketed(newtree);
		  std::cout << std::endl;
   }
	 
   //population.PrintOut(output_file);
}
