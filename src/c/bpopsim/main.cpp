#include "cSubpopulation.h"
#include "commandline.h"
#include "cPopulation.h"
#include "cPopulation.cc"

const int size_of_by_color = 2;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <iostream>

#include <cmath>
#include <string>
#include <vector>

using namespace std;

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

	// Simulation parameters that should be arguments
	uint64_t initial_population_size = 2;
	uint64_t pop_size_after_dilution = int(5E6);             // N sub 0 --int is to get rid of warning
	double mutation_rate_per_division = 1E-8;         // mu
	double average_mutation_s = 0.05;                 // s
  	double growth_phase_generations = 6.64;
	//char beneficial_mutation_distribution_code = 'e'; //e = exponential, u = uniform


	population.SetBinomialSamplingThreshold(1000);
	std::string output_file_name("output.txt");
  
	

	

	//Output parameters that should be arguments
	
	population.SetTransferIntervalToPrint(1);
	population.SetVerbose(1);
	
	population.SetTotalTransfers(200000);
	population.SetMaxDivergenceFactor(100);
	population.SetReplicates(1000);
	population.SetMinimumPrinted(8);

	
	
	// Simulation parameters that are pre-calculated
		
	population.SetDilutionFactor(exp(log(2) * growth_phase_generations));
	population.SetTransferBinomialSamplingP(1/population.GetDilutionFactor());
	population.SetPopSizeBeforeDilution(pop_size_after_dilution * population.GetDilutionFactor());
	population.SetLambda(1/mutation_rate_per_division);
	//std::vector< std::vector<double> > runs;
	//std::vector<double> this_run;

	//Create population of subpopulations
	
	if (population.GetVerbose()==1) 
	{
		cout << "u = " << mutation_rate_per_division << endl;
		cout << "s = " << average_mutation_s << endl;
		cout << "N = " << pop_size_after_dilution << endl;
		cout << "dil = " << population.GetDilutionFactor() << endl;
	}

	for (int on_run=0; on_run < population.GetReplicates(); on_run++)
	{	
		population.ClearRuns();
    		cout << "Replicate " << on_run+1 << endl;    

		//This is the vector for the ratio information for each run.
		//vector<double> this_run;	 
                //Push all of the ratio info for one run into this vector, then push vector into runs	

		//vector<cSubpopulation> populations; //Vector containing all of the populations;
		//Seed a red population


		
		//Create red population
		cSubpopulation r;
		//Set parameters
		r.SetNumber(initial_population_size/2);
    		r.SetFitness(1);
    		r.SetMarker('r');
		//populations.push_back(r);
		//Add subpopulation to population
		population.AddSubpopulation(r);
		
		//Seed a white population
		cSubpopulation w;
		w.SetNumber(initial_population_size/2);
    		w.SetFitness(1);
    		w.SetMarker('w');	
    		population.AddSubpopulation(w);
		population.SetTotalPopSize(initial_population_size);
		population.SetTotalMutations(0);		
		population.SetTotalSubpopulationsLost(0);
  
		population.SetTransfers(1);		
		population.SetDivisionsUntilMutation(0);
		population.SetKeepTransferring(true);

		while ( (population.GetTransfers() < population.GetTotalTransfers()) && population.GetKeepTransferring() )
		{
			// Move time forward until another mutation occurs.

			//divisions_until_mutation += floor(returnExp(lambda));
			population.SetDivisionsUntilMutation(population.GetDivisionsUntilMutation() + gsl_ran_exponential(randgen, population.GetLambda()));
			if (population.GetVerbose()) cout << "  New divisions before next mutation: " << population.GetDivisionsUntilMutation() << endl;
			
			// Calculate points to output between last time and current time
			// First time through the loop, a partial time interval
			// may be calculated to get back on $print_interval

      			// Move forward by a large chunk of time which assumes
     		        // all populations have the maximum fitness in the population
		        // This will, at worst, underestimate how long.
		        // We can then move forward by single divisions to find the exact division where the mutation occurs
		      
			while ( (population.GetDivisionsUntilMutation() > 0) && (population.GetTransfers() < population.GetTotalTransfers()) && 					population.GetKeepTransferring()) 
			{
        		//Keep track of time in a perfectly integrated fashion
        
       			 	population.SetDesiredDivisions(population.GetDivisionsUntilMutation());
       			 	if (population.GetDesiredDivisions() + population.GetTotalPopSize() > population.GetPopSizeBeforeDilution())
				{
       			   		population.SetDesiredDivisions(population.GetPopSizeBeforeDilution() - population.GetTotalPopSize());
        			}
        			if(population.GetVerbose() == 1) 
				{
			        	cout << "Divisions before next mutation: " << population.GetDivisionsUntilMutation() << endl;
			        }
			        // Note: we underestimate by a few divisions so that we can step forward by single division increments
			        // as we get close to the one where the mutation happened (or right before a transfer).			
			        if (population.GetDesiredDivisions() < 1)
				{
        				population.SetDesiredDivisions(1);
        			}
        			if (population.GetVerbose() == 1) 
				{
			        	cout << "Total pop size: " << population.GetTotalPopSize() << endl;
			        	cout << "Desired divisions " << population.GetDesiredDivisions() << endl;
			        }
			
			        // How much time would we like to pass to achieve the desired number of divisions?
			        // (assuming the entire population has the maximum fitness)
			        population.SetUpdateTime(log((population.GetDesiredDivisions()+(double)population.GetTotalPopSize()) / population.GetTotalPopSize()) / (population.GetMaxW() * log(2)));
        
    				//What is the minimum time required to get a single division?
				population.SetTimeToNextWholeCell(0);

				population.DetermineDivisionTime();
        			// At a minumum, we want to make sure that one cell division took place
        			if (population.GetTimeToNextWholeCell() > population.GetUpdateTime()) 
				{
        				if (population.GetVerbose()) cout << "Time to next whole cell greater than update time: " << population.GetTimeToNextWholeCell() << 										     " < " << population.GetUpdateTime() << endl;		
        				population.SetUpdateTime(population.GetTimeToNextWholeCell());
        			}

        			if (population.GetVerbose() == 1) 
				{
 				       	cout << "Update time: " << population.GetUpdateTime() << endl;		
        			}
            
        			//Now update all lineages by the time that actually passed
	
				population.SetNewPopSize(0);

				population.UpdateLineages();
				cout << population.GetNewPopSize() << " " << population.GetTotalPopSize() << endl;
			        population.SetCompletedDivisions(population.GetNewPopSize() - population.GetTotalPopSize());
			        
				if (population.GetVerbose()) cout << "Completed divisions: " << population.GetCompletedDivisions() << endl;
			        population.SetDivisionsUntilMutation(population.GetDivisionsUntilMutation() - population.GetCompletedDivisions());
			        population.SetTotalPopSize(population.GetNewPopSize());
	
				population.Mutate(randgen);
      
				population.Resample(randgen);				

					        
			}
		
		}
		
    	
		cout << "Total mutations: " << population.GetTotalMutations() << endl;
	    	cout << "Total subpopulations lost: " << population.GetTotalSubpopulationsLost() << endl;
	    	cout << "Transfers: " << population.GetTransfers() << endl;
	    	cout << "Maximum Fitness: " << population.GetMaxW() << endl;		
		population.PushBackRuns();
	}
	


	population.PrintOut();

	
}


