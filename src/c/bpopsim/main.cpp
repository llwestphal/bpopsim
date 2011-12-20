#include "common.h"
#include "cPopulation.h"

#include <iostream>

using namespace bpopsim;
using namespace std;

// Global variable for keeping track of verbosity
bool g_verbose = false;
bool g_ro_only = false;

int main(int argc, char* argv[])
{
  AnyOption options("Usage: bpopsim etc");
  options
  ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
  ("generations-per-transfer,T", "Generations per transfer", 6.64)
  ("population-size-after-transfer,N", "Population size after transfer", double(5E6))
  ("number-of-transfers,n", "Max number of transfer to replicate", 50)
  ("output-folder,o", "Output folder")
  ("input_file,i", "Input file for converting external tree")
  ("mutation-rate-per-division,u", "Mutation rate per division", 0)
  ("initial-population-size,p", "Initial Population Size", double(5E6))
  ("replicates,r", "Replicates")
  ("marker-divergence,m", "Max divergence factor", 100)
  ("type-of-mutations,f", "Type of mutations", 'u')
  ("verbose,v", "Verbose", TAKES_NO_ARGUMENT)
  ("lineage-tree,l", "Lineage Tree")
  ("seed,d", "Seed for random number generator")
  ("red-white,k", "Only care about red/white lineages. For marker divergence.", TAKES_NO_ARGUMENT)
  ("transfer-interval-to-print,t", "Red/White Printing intervals.", 1)
  ("imv,x", "Fitness gain per mutation.")
  ("coarse-graining,c", "Amount to coarse-grain output.", 1)
  ("initial-fitness,z", "Initial Fitness.", 1.0)
  ("mut_num,y", "Specify an exact number of mutations per transfer.", 1)
  ("muller_res,w", "Muller plot vertical resolution.", 2500)
  
  //Here are the output options
  //They are all set to false by default
  ("frequencies", "Print Frequencies", TAKES_NO_ARGUMENT)
  ("muller", "Print Muller Matrix", TAKES_NO_ARGUMENT)
  ("average_fit", "Print Average Fitness", TAKES_NO_ARGUMENT)
  ("print_screen", "Print Frequencies to Screen", TAKES_NO_ARGUMENT)
  ("time_sweep", "Print time to sweep for each mutation", TAKES_NO_ARGUMENT)
  ("max_diff", "Print max difference of sweeping mutations", TAKES_NO_ARGUMENT)
  ("single_fit", "Print the fitness of the single cell.", TAKES_NO_ARGUMENT)
  ("convert_tree", "Convert an external phylogenetic tree.", TAKES_NO_ARGUMENT)
  ("sweeping_descent_fitness", "Output the average fitness of sweeping descent per time.", TAKES_NO_ARGUMENT)
  ("new_simulator", "Use the new simulator type.", TAKES_NO_ARGUMENT)
  .processCommandArgs(argc, argv);
  
  if(options.count("help")
     && !options.count("frequencies")
     && !options.count("muller")
     && !options.count("average_fit")
     && !options.count("print_screen")
     && !options.count("time_sweep")
     && !options.count("max_diff")
     && !options.count("single_fit")
     && !options.count("sweeping_descent_fitness")
     ) {
    cout << "You have not chosen any output options." << endl;
		options.printUsage();
		return -1;
	}                       
  
	try {
    
    bool print_freq(false), print_muller(false), print_average_fit(false),
    print_screen(false), print_max_diff(false), print_time_to_sweep(false),
    print_single_fit(false), use_mute_num(false), new_simulator(false),
    print_sweeping_descent_fitness(false);
    
    if( options.count("verbose") ) g_verbose = true;
    if( options.count("frequencies") ) print_freq = true;
    if( options.count("muller") ) print_muller = true;
    if( options.count("average_fit") ) print_average_fit = true;
    if( options.count("print_screen") ) print_screen = true;
    if( options.count("time_sweep") ) print_time_to_sweep = true;
    if( options.count("max_diff") ) print_max_diff = true;
    if( options.count("single_fit") ) print_single_fit = true;
    if( options.count("mut_num") ) use_mute_num = true;
    if( options.count("new_simulator") ) new_simulator = true;
    if( options.count("sweeping_descent_fitness") ) print_sweeping_descent_fitness = true;
    
    if( options.count("convert_tree") ) {
      cPopulation access_to_functions;
      
      access_to_functions.ConvertExternalData(options["input_file"]);
        
      if( print_muller ) {
        std::cout << std::endl << "Generating Muller Matrix.... " << std::endl;
        std::vector< std::vector<int> > muller_matrix;
        access_to_functions.DrawMullerMatrix(options["output-folder"], muller_matrix);
      }
    
    }
    
    //@agm The new simulator
    else if( new_simulator ) {
      // @JEB: use one main RNG (so we don't reset for each replicate)
      // Let cPopulation store a copy of it.
      const gsl_rng_type *T;
      gsl_rng * randgen;
      gsl_rng_env_setup();
      T = gsl_rng_taus2;
      randgen = gsl_rng_alloc(T);
      
      uint16_t seed = 0;
      if ( options.count("seed") ) {
        seed = from_string<uint16_t>(options["seed"]);
      } else {
        seed = time(NULL) * getpid();
      }
      
      gsl_rng_set(randgen, seed);
      
      for (uint32_t on_run=0; on_run < from_string<uint32_t>(options["replicates"]); on_run++)
      {
        
        std::vector<double> current_ro_ratio;
        
        //Initialize Population object
        cPopulation population;
        
        //Set cli options
        population.SetParameters(options);
        population.DisplayParameters();
        
        tree<cGenotype>::iterator_base loc;
        std::vector< std::vector<cGenotypeFrequency> > subpops;
        
        population.SetRNG(randgen);
        
        std::cout << "Replicate " << on_run << std::endl;   
        
        //Initialize the population
        if (g_ro_only) {
          population.SeedSubpopulationForRedWhite(); 
        } else {
          population.SeedPopulationWithOneColony();
        }
        
        //Get an initial time points
        population.CalculateAverageFitness();
        population.FrequenciesPerTransferPerNode();
        
        // Print the initial tree
        //if (g_verbose) population.PrintTree();
        
        //std::cout << node_id << std::endl;
        
        uint32_t count(0);

        while( population.GetTransfers() < population.GetTotalTransfers() ) {
          
          population.SetDivisionsUntilMutation(round(gsl_ran_exponential(randgen, population.GetLambda())));
          
          while( population.GetDivisionsUntilMutation() > 0 && (population.GetTransfers() < population.GetTotalTransfers())) {
            population.CalculateDivisionsNew();
            
            if( population.GetDivisionsUntilMutation() <= 0 ) population.MutateNew();
            
            if( population.GetPopulationSize() >= population.GetPopSizeBeforeDilution() ) {
              population.FrequenciesPerTransferPerNode();
              
              if( on_run == 0 ) 
                population.CalculateAverageFitness();
              
              if( print_single_fit )
                population.Deterministic_Resample();
              else {
                population.Resample();
                population.CullPopulations();
              }
              
              count++;
              std::cout << "Passing.... " << count << std::endl;
              
              if( print_single_fit && !use_mute_num) {
                population.PrintSingleFitness(options["output-folder"]);
                //std::cout << "Population size: " << population.GetPopulationSize() << std::endl;
              }
            }
          }
        }
        
        if( print_freq ) {
          std::cout << "Printing to file.... \n";
          population.PrintOut(options["output-folder"], on_run);
          cout << endl;
        }
        
      }
      
    }
    else {
      
      //for testing!!!!
      vector<double> mutations_until_division(0);
      
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
      if ( options.count("seed") ) {
        seed = from_string<uint16_t>(options["seed"]);
      } else {
        seed = time(NULL) * getpid();
      }
      
      gsl_rng_set(randgen, seed);
      
      if (g_ro_only) std::cout << std::endl << "You chose to use red and white only." << std::endl;
      
      //std::vector< std::vector<uint16_t> > number_unique_genotypes_in_all_replicates;
      
      uint32_t replicates(0);
      
      if( use_mute_num && from_string<uint32_t>(options["mut_num"]) == 1) replicates = pow(2, from_string<double>(options["generations-per-transfer"]));
      else replicates = from_string<uint32_t>(options["replicates"]);
      
      uint32_t start_num(0);
      if( use_mute_num ) start_num = 1;
      
      for (uint32_t on_run=start_num; on_run < replicates; on_run++)
      {
        std::vector<double> current_ro_ratio;
        
        //Initialize Population object
        cPopulation population;
        
        //Build lookup table for logs 
        //Currently it is to the 15th, I should take it as a command line option
        //population.ConstructLookUpTable();
        
        //Set cli options
        population.SetParameters(options);
        population.DisplayParameters();
        
        tree<cGenotype>::iterator_base loc;
        std::vector< std::vector<cGenotypeFrequency> > subpops;
        
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
        population.FrequenciesPerTransferPerNode();
        
        // Print the initial tree
        //if (g_verbose) population.PrintTree();
        
        //std::cout << node_id << std::endl;
        
        uint16_t number_of_mutations(0);
        vector<uint32_t> mutation_division(1, 0);

        double division_of_mutation;
        
        if( use_mute_num ) {
          for(uint32_t i=0; i<from_string<uint32_t>(options["mut_num"]); i++) {
            division_of_mutation = (double) gsl_rng_uniform_int(randgen, pow(2, from_string<double>(options["generations-per-transfer"])));
            mutation_division.push_back( division_of_mutation );
            mutations_until_division.push_back( division_of_mutation );
          }
          
          /*uint32_t sum(0);
          for(uint32_t i=0; i<from_string<uint32_t>(options["mut_num"]); i++) {
            sum+=mutation_division[i];
            //cout << mutation_division[i] << " " << pow(2,from_string<double>(options["generations-per-transfer"])) << endl;
          }*/
          
          //cout << sum << " " << pow(2,from_string<double>(options["generations-per-transfer"])-1) << endl;
          
          sort(mutation_division.begin(), mutation_division.end());
        }
        
        while( population.GetTransfers() < population.GetTotalTransfers() ) {
          
          if (g_verbose) {
            for(vector<cSubpopulation>::iterator this_time = population.GetPopulation().begin(); this_time != population.GetPopulation().end(); this_time++) {
              cout << "Genotype3: " << this_time->GetNode_id() << " Frequency3: " << this_time->GetNumber() << endl;
            }
          }
          
          // Calculate the number of divisions until the next mutation 
          if( use_mute_num && number_of_mutations < from_string<uint32_t>(options["mut_num"]) && from_string<uint32_t>(options["mut_num"]) == 1) {
            population.SetDivisionsUntilMutation( on_run );
          }
          
          else if( use_mute_num && number_of_mutations >= from_string<uint32_t>(options["mut_num"]) && from_string<uint32_t>(options["mut_num"]) == 1 ) {
            population.SetDivisionsUntilMutation( pow(2, from_string<double>(options["generations-per-transfer"])) );
          }
          
          else if ( use_mute_num && from_string<uint32_t>(options["mut_num"]) != 1 ) {
            population.SetDivisionsUntilMutation( mutation_division[number_of_mutations+1] - mutation_division[number_of_mutations] );
            
            //population.SetDivisionsUntilMutation( pow(2, from_string<double>(options["generations-per-transfer"]))/pow(2, from_string<double>(options["generations-per-transfer"])/2) );
            //cout << population.GetPopulationSize() << " " << population.GetDivisionsUntilMutation() << " " << mutation_division.size() << " " << mutation_division[number_of_mutations] << " " << pow(2,from_string<double>(options["generations-per-transfer"])) << endl;
          }
          
          else if ( !use_mute_num ) {
            population.SetDivisionsUntilMutation(population.GetDivisionsUntilMutation() + round(gsl_ran_exponential(randgen, population.GetLambda())));
          }
          
          if (g_verbose) {
            for(vector<cSubpopulation>::iterator this_time = population.GetPopulation().begin(); this_time != population.GetPopulation().end(); this_time++) {
              cout << "Genotype4: " << this_time->GetNode_id() << " Frequency4: " << this_time->GetNumber() << endl;
            }
          }
          
          while( population.GetDivisionsUntilMutation() > 0 && (population.GetTransfers() < population.GetTotalTransfers())) 
          {
            
            if (g_verbose) {
              for(vector<cSubpopulation>::iterator this_time = population.GetPopulation().begin(); this_time != population.GetPopulation().end(); this_time++) {
                cout << "Genotype5: " << this_time->GetNode_id() << " Frequency5: " << this_time->GetNumber() << endl;
              }
            }
            
            //std::cout << population.GetTransfers() << " " << population.GetTotalTransfers() << " " << population.GetDivisionsUntilMutation() << std::endl;
            
            population.CalculateDivisions();
            
            if (g_verbose) {
              for(vector<cSubpopulation>::iterator this_time = population.GetPopulation().begin(); this_time != population.GetPopulation().end(); this_time++) {
                cout << "Genotype7: " << this_time->GetNode_id() << " Frequency7: " << this_time->GetNumber() << endl;
              }
            }
            
            if( population.GetDivisionsUntilMutation() <= 0 ) { 
              if( !use_mute_num ) population.Mutate();
              
              else if( use_mute_num && number_of_mutations <  from_string<uint32_t>(options["mut_num"]) ) {
                population.Mutate();
                number_of_mutations++;
                cout << "Number of mutations: " << number_of_mutations << endl;
              }
              
            }
            
            if( population.GetPopulationSize() >= population.GetPopSizeBeforeDilution()) {
              
              population.FrequenciesPerTransferPerNode();
              population.CalculateAverageFitness();
              
              if( print_single_fit )
                population.Deterministic_Resample();
              else {
                population.Resample();
                population.CullPopulations();
              }
              
              if( print_single_fit && !use_mute_num) {
                population.PrintSingleFitness(options["output-folder"]);
                //std::cout << "Population size: " << population.GetPopulationSize() << std::endl;
              }
              
              if( g_verbose ) std::cout << "Total pop size3: " << population.GetPopulationSize() <<std::endl;
              
              count++;
              std::cout << "Passing.... " << count << std::endl;

              if ( population.GetTransfers() % from_string<uint16_t>(options["transfer-interval-to-print"]) == 0 ) {  
                current_ro_ratio.push_back(population.GetRatio());
                
                if (g_verbose)
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
        
        //std::cout << std::endl << std::endl;
        //population.RunSummary();
        population.PushBackRuns();

        if (g_verbose) {
          population.PrintTree();
        }
        
        if (g_ro_only) {
          red_white_ratios.push_back(current_ro_ratio);
        }
        else {
          //number_unique_genotypes_in_all_replicates.push_back(population.CurrentUniqueGenotypes());
          
          if( use_mute_num ) {
            std::cout << "Printing expectation value for set number of mutations.... \n";
            population.PrintExpectationValue(options["output-folder"]);
            cout << endl;
          }
          if( print_screen ) {
            std::cout << "Printing to screen.... \n";
            population.PrintFrequenciesToScreen(options["output-folder"]);
            cout << endl;
          }
          
          if( print_max_diff ) {
            std::cout << "Printing max difference of relevant mutations.... \n";
            population.CalculateSimilarity(options["output-folder"]);
            cout << endl;
          }
      
          if( print_time_to_sweep ) {
            std::cout << "Printing time to sweep.... \n";
            population.TimeToSweep(options["output-folder"]);
            cout << endl;
          }
          
          if( print_freq ) {
            std::cout << "Printing to file.... \n";
            population.PrintOut(options["output-folder"], on_run);
            cout << endl;
          }
          
          if( print_sweeping_descent_fitness ) {
            cout << "Printing fitness of the winning line of descent.... \n";
            population.PrintWinningFitness(options["output-folder"], on_run);
            cout << endl;
          }
          
          if( print_average_fit ) {
            std::cout << "Printing average fitness.... \n";
            population.PrintFitness(options["output-folder"], on_run);
            cout << endl;
          }
          
          if( on_run == 0 ) {
        
            if( print_muller ) {
              std::cout << std::endl << "Generating Muller Matrix.... \n";
              std::vector< std::vector<int> > muller_matrix;
              population.DrawMullerMatrix(options["output-folder"], muller_matrix);
              cout << endl;
            }
          }
        }
        
        //cout << endl << endl;
        //population.PrintFreqsQuick();
        //population.PrintTree();
        
      }
      
      if (g_ro_only) {
        //Initialize Population object
        cPopulation population;
        std::cout << std::endl << "Printing to r/w ratio file.... " << std::endl;
        population.PrintOut_RedWhiteOnly(options["output-folder"], &red_white_ratios, from_string<uint16_t>(options["transfer-interval-to-print"]));
      }
      /*else {
        //Initialize Population object
        cPopulation population;
        std::cout << std::endl << "Printing unique genotypes to file.... " << std::endl;
        population.PrintUniqueGenotypes(output_folder, &number_unique_genotypes_in_all_replicates);
      }*/
      
      /*double sum(0);
      for(vector<double>::iterator i = mutations_until_division.begin(); i != mutations_until_division.end(); i++) {
        sum+=log(*i);
      }
      cout << "Here is the average of time of mutation: " << sum/mutations_until_division.size();*/
    }
  } catch(...) {
    // failed; 
    return -1;
  }
    
  return 0;
}