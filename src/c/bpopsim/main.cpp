#include "cSubpopulation.h"
#include "cPopulation.h"

#include <iostream>
#include <string>

// Boost
#include <boost/program_options.hpp>

// GSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// Tree
#include "tree.hh"

using namespace boost::program_options;
using namespace std;

// setup and parse configuration options:
void get_cmdline_options(variables_map &options, int argc, char* argv[]) {

  options_description cmdline_options("Allowed options");
  cmdline_options.add_options()
  ("help,h", "produce this help message")
  ("generations-per-transfer,T", value<double>(), "Generations per transfer")
  ("population-size-after-transfer,N", value<uint64_t>(), "Population size after transfer")
  ("output-file,o", value<string>(), "Output file")
  ("mutation-rate-per-generation,u", value<double>(), "Mutation rate per generation")
  ("average-selection-coefficient,s", value<double>(), "Average selection coefficient")
  ("time-interval,i", value<int>(), "Time interval")
  ("replicates,r", value<int>(), "Replicates")
  ("marker-divergence,m", value<int>(), "Max divergence factor")
  ("verbose,v", value<int>(), "Verbose")
  ;

/* Need to add these as options...
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


  //TREE

  cSubpopulation testone;

  testone.SetNumber(1);
  cSubpopulation testtwo;
  testtwo.SetNumber(2);
  cSubpopulation testthree;
  testthree.SetNumber(3);
  cSubpopulation testfour;
  testfour.SetNumber(4);

  cSubpopulation testfive;
  testfive.SetNumber(5);
  cSubpopulation testsix;
  testsix.SetNumber(6);
  cSubpopulation testseven;
  testseven.SetNumber(7);
  cSubpopulation testeight;
  testeight.SetNumber(8);


  tree<cSubpopulation*> trtr;
  tree<cSubpopulation*>::iterator toptop,oneone,twotwo,threethree,locloc;
  
  toptop=trtr.begin();
  
  cSubpopulation* testpointer = &testone;
  cSubpopulation* testpointertwo = &testtwo;
  cSubpopulation* testpointerthree = &testthree;
  cSubpopulation* testpointerfour = &testfour;

  cSubpopulation* testpointerfive = &testfive;
  cSubpopulation* testpointersix = &testsix;
  cSubpopulation* testpointerseven = &testseven;
  cSubpopulation* testpointereight = &testeight;
  

  //cout << testpointer->GetNumber() << endl; 
  oneone=trtr.insert(toptop,testpointer);
  twotwo=trtr.append_child(oneone,testpointertwo);
  
  trtr.append_child(twotwo,testpointerthree);
  threethree=trtr.append_child(twotwo,testpointerfour);
  trtr.append_child(threethree,testpointerfive);
  trtr.append_child(twotwo,testpointersix);
  trtr.append_child(oneone,testpointerseven);
  trtr.append_child(oneone,testpointereight);
  locloc=trtr.begin();
  tree<cSubpopulation*>::iterator sibsib=trtr.begin();
  testone.SetNumber(2000);
  

  while(sibsib!=trtr.end()) {
  cout << (*sibsib)->GetNumber() << endl;
  ++sibsib;
  }

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
          
      while ( (population.GetDivisionsUntilMutation() > 0) && (population.GetTransfers() < population.GetTotalTransfers()) &&           population.GetKeepTransferring()) 
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


