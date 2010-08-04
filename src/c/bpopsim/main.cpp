#include "cSubpopulation.h"
#include "commandline.h"
#include "cPopulation.h"
#include "cPopulation.cc"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{

	
	//create generator and seed
	const gsl_rng_type * T;
        gsl_rng * randgen;
	gsl_rng_env_setup();
        T = gsl_rng_mt19937;
        randgen = gsl_rng_alloc (T);

	cline(argc, argv);
	if (argc!=1) 
	{
		return 0;
	}

	cPopulation population;

	population.SetParameters();

	std::string output_file_name("output.txt");
	
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
	


	population.PrintOut();

	
}


