#include "cSubpopulation.h"
#include "commandline.h"

const int size_of_by_color = 2;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <iostream>
#include <fstream>
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

	double update_time;

	// Simulation parameters that should be arguments
	uint64_t initial_population_size = 2;
	uint64_t pop_size_after_dilution = int(5E6);             // N sub 0 --int is to get rid of warning
	double mutation_rate_per_division = 1E-8;         // mu
	double average_mutation_s = 0.05;                 // s
  	double growth_phase_generations = 6.64;
	char beneficial_mutation_distribution_code = 'e'; //e = exponential, u = uniform
	double binomial_sampling_threshold = 1000;
	std::string output_file_name("output.txt");
  
	//Output parameters that should be arguments
	int total_transfers = 200000;
	double max_divergence_factor = 100;
	double verbose = 1;  
	int replicates = 1000;
	int transfer_interval_to_print = 1;
	int minimum_printed = 8;

	// Simulation parameters that are pre-calculated
	const double dilution_factor = exp(log(2) * growth_phase_generations);
	const double  transfer_binomial_sampling_p = 1/dilution_factor;
	double pop_size_before_dilution = pop_size_after_dilution * dilution_factor;
	double lambda = mutation_rate_per_division;
	vector< vector<double> > runs;

	if (verbose==1) {
		cout << "u = " << mutation_rate_per_division << endl;
		cout << "s = " << average_mutation_s << endl;
		cout << "N = " << pop_size_after_dilution << endl;
		cout << "T = " << growth_phase_generations << endl;
		cout << "dil = " << dilution_factor << endl;
	}

	for (int on_run=0; on_run < replicates; on_run++) {
    		cout << "Replicate " << on_run+1 << endl;
    
		vector<double> this_run; //This is the vector for the ratio information for each run. 
                             //Push all of the ratio info for one run into this vector, then push vector into runs	

		vector<cSubpopulation> populations; //Vector containing all of the populations;

		//Seed a red population
		cSubpopulation r;
		r.SetNumber(initial_population_size/2);
    r.SetFitness(1);
    r.SetMarker('r');
		populations.push_back(r);

		//Seed a white population
		cSubpopulation w;
		w.SetNumber(initial_population_size/2);
    w.SetFitness(1);
    w.SetMarker('w');
		populations.push_back(w);		
    
		double max_w = 1;
		double total_pop_size = initial_population_size;
		int total_mutations = 0;
		int total_subpopulations_lost = 0;
	
  
		int transfers = 1;
		double divisions_until_mutation = 0;
		bool keep_transferring = true;

		while ( (transfers < total_transfers) && keep_transferring )
		{
			// Move time forward until another mutation occurs.

			//divisions_until_mutation += floor(returnExp(lambda));
			divisions_until_mutation += gsl_ran_exponential(randgen, 1/lambda);
	if (verbose) cout << "  New divisions before next mutation: " << divisions_until_mutation << endl;
			
			// Calculate points to output between last time and current time
			// First time through the loop, a partial time interval
			// may be calculated to get back on $print_interval

      // Move forward by a large chunk of time which assumes
      // all populations have the maximum fitness in the population
      // This will, at worst, underestimate how long.
      // We can then move forward by single divisions to find the exact division where the mutation occurs
      
		while ( (divisions_until_mutation > 0) && (transfers < total_transfers) && keep_transferring ) {
        	//Keep track of time in a perfectly integrated fashion
        
        	double desired_divisions = divisions_until_mutation;
        	if (desired_divisions + total_pop_size > pop_size_before_dilution) {
          		desired_divisions = pop_size_before_dilution - total_pop_size;
        	}
        if(verbose == 1) {
        	cout << "Divisions before next mutation: " << divisions_until_mutation << endl;
        }
        // Note: we underestimate by a few divisions so that we can step forward by single division increments
        // as we get close to the one where the mutation happened (or right before a transfer).			
        if (desired_divisions < 1) {
        	desired_divisions = 1;
        }
        if (verbose == 1) {
        	cout << "Total pop size: " << total_pop_size << endl;
        	cout << "Desired divisions " << desired_divisions << endl;
        }
			
        // How much time would we like to pass to achieve the desired number of divisions?
        // (assuming the entire population has the maximum fitness)
        update_time = log( (desired_divisions+(double)total_pop_size) / total_pop_size ) / (max_w * log(2));
        
        //What is the minimum time required to get a single division?
        double time_to_next_whole_cell = 0;
        
        vector<int> divided_lineages;	
        for (int i=0; i< int(populations.size()); i++) {
        	cSubpopulation &sp = populations[i];
        	if (sp.GetNumber() == 0) continue;
        
          //what is the time to get to the next whole number of cells?
        double current_cells = sp.GetNumber();
        double whole_cells = floor(sp.GetNumber())+1;
        // WC = N * exp(growth_rate * t) 
          
        double this_time_to_next_whole_cell = log(whole_cells / current_cells) / (sp.GetFitness());
          
        if ( time_to_next_whole_cell == 0 || (this_time_to_next_whole_cell < time_to_next_whole_cell) ) {
          	divided_lineages.clear();
          	time_to_next_whole_cell = this_time_to_next_whole_cell;
          	divided_lineages.push_back(i); //a list, because there can be ties
          }
          else if (this_time_to_next_whole_cell == time_to_next_whole_cell)
          {
          	divided_lineages.push_back(i); //a list, because there can be ties
          }
        }
	
        // At a minumum, we want to make sure that one cell division took place
        if (time_to_next_whole_cell > update_time) {
        	if (verbose) cout << "Time to next whole cell greater than update time: " << time_to_next_whole_cell << " > " << update_time << endl;		
        	update_time = time_to_next_whole_cell;
        }

        if (verbose == 1) {
        	cout << "Update time: " << update_time << endl;		
        }
            
        //Now update all lineages by the time that actually passed
        long double new_pop_size = 0;
        
        for (vector<cSubpopulation>::iterator it = populations.begin(); it!=populations.end(); ++it) {
        	if (it->GetNumber() == 0) continue;
        	it->SetNumber(it->GetNumber() * exp(log(2) * update_time * it->GetFitness()));
          	new_pop_size += it->GetNumber();
        }				
                
        long double completed_divisions = new_pop_size - total_pop_size;
        if (verbose) cout << "Completed divisions: " << completed_divisions << endl;
        divisions_until_mutation = divisions_until_mutation - completed_divisions;
        total_pop_size = new_pop_size;
		
        //cout << "Update_time\t" << "New_pop_size\t" << "Total_pop_size\t" << "Completed_divisions" << endl;
        //cout << update_time << "\t" << new_pop_size << "\t" << total_pop_size << "\t" << completed_divisions << endl;
        //cout << "Divisions until mutation: " << divisions_until_mutation << endl;				
        if (divisions_until_mutation <= 0) {
         	total_mutations++;
          	if (verbose) cout << "* Mutating!" << endl;
        
          	//Mutation happened in the one that just divided.
          	// ==> This seems to not be entirely accurate...
          
          	//Break ties randomly here.
          	cSubpopulation& ancestor = populations[divided_lineages[rand() % divided_lineages.size()]];
          


          	//What is the new lineage's fitness?
          	//Option #1: dirac delta
          	


		//double new_w = 0;
          	

		//if (beneficial_mutation_distribution_code == 0) {
            	//new_w = ancestor.GetFitness() + average_mutation_s;
        	//}
        	//Option #2: exponential
        	//else {
		//double this_mutation_s = gsl_ran_exponential(randgen,1/average_mutation_s);
         	//double this_mutation_s = returnExp(average_mutation_s);
         	//new_w = ancestor.GetFitness() + this_mutation_s;
        	//}
          
        	
          
          //Create and add the new lineage

        cSubpopulation new_lineage = ancestor.CreateDescendant(randgen);

	if (verbose) cout << "  Color: " << ancestor.GetMarker() << endl;
        if (verbose) cout << "  New Fitness: " << new_lineage.GetFitness() << endl;


        populations.push_back(new_lineage);
          
        //Update maximum fitness
        if(new_lineage.GetFitness() > max_w) {
        	max_w = new_lineage.GetFitness();
        }
        
        }
      
        //When it is time for a transfer, resample population
        if (total_pop_size >= pop_size_before_dilution) {
        
        	if (verbose) cout << ">> Transfer!" << endl;
        	//Is there an exists() in C++?

        	double by_color[size_of_by_color];
        	by_color[0] = 0;
        	by_color[1] = 0;
          
        	long double new_pop_size = 0;
        	//POPULATION: for (my $i=0; $i< scalar @populations; $i++) 	
          	for (vector<cSubpopulation>::iterator it = populations.begin(); it!=populations.end(); ++it) {
            		if (it->GetNumber() == 0) continue;
            
            // Perform accurate binomial sampling only if below a certain population size
            	if (it->GetNumber() < binomial_sampling_threshold) {
              		if (verbose) cout << "binomial" << it->GetNumber() << endl;
				
            			it->Transfer(transfer_binomial_sampling_p, randgen);
				//it->SetNumber(gsl_ran_binomial(randgen, transfer_binomial_sampling_p, uint64_t(it->GetNumber())));
				//it->n = returnBin(it->n, transfer_binomial_sampling_p);
            		}
            // Otherwise, treat as deterministic and take expectation...
            else {
            	it->SetNumber(it->GetNumber() * transfer_binomial_sampling_p);
            }
            
            // Keep track of lineages we lost
            if (it->GetNumber() == 0) {
            	total_subpopulations_lost++;
            }
            

            new_pop_size += floor(it->GetNumber());
            //There is probably a better way to do this, I just don't know the syntax ??by_color[p.color] += p.n;??
            if (verbose) cout << it->GetMarker() << endl;
            if (verbose) cout << it->GetNumber() << endl;
            if (verbose) cout << it->GetFitness() << endl;
            if (it->GetMarker() == 'r') {
            	by_color[0] += it->GetNumber();
            }
            else {
            	by_color[1] += it->GetNumber();
            }
          }
          //One color was lost, bail	
          if ( (by_color[0] == 0) || (by_color[1] == 0) ) {
          	keep_transferring = false;
          }
            
          if (verbose) cout << "Colors: " << by_color[0] << " / " << by_color[1] << endl;
          long double ratio = by_color[0] / by_color[1];
          
          transfers++;

          if ( (transfers >= 0) && (transfers % transfer_interval_to_print == 0) ) {
            this_run.push_back(ratio);
            if (verbose == 1) {
                cout << "Transfer " << transfers << " : " << total_pop_size << "=>" << new_pop_size << "  R/W Ratio: " << ratio << endl;	
                cout << "Total mutations: " << total_mutations << " Maximum Fitness: " << max_w << endl;
                cout << "Size = " << this_run.size() << endl;
            }
          }
          
          total_pop_size = new_pop_size;

          if ( (int(this_run.size()) >= minimum_printed) && ((ratio > max_divergence_factor) || (ratio < 1/max_divergence_factor)) ) {
            if (verbose) cout << "DIVERGENCE CONDITION MET" << endl;
            keep_transferring = false;
          }	
        }
      }
    }
    
    cout << "Total mutations: " << total_mutations << endl;
    cout << "Total subpopulations lost: " << total_subpopulations_lost << endl;
    cout << "Transfers: " << transfers << endl;
    cout << "Maximum Fitness: " << max_w << endl;

    runs.push_back(this_run);
  }
	
  //Print everything out
  ofstream output_file;
  output_file.open ("output.txt");
  output_file << "transfer";
  for (int on_run=0; on_run < replicates; on_run++) {
    output_file << "\t" << on_run;
  }
  output_file << endl;

  bool still_going = true;
  int i = 0;
  while (still_going) {    
    int on_transfer = i * transfer_interval_to_print;
    output_file << on_transfer;
    still_going = false;
    for (int on_run=0; on_run < replicates; on_run++) {
      if(i<int(runs[on_run].size())) {
        
        output_file << "\t" << runs[on_run][i];
        still_going = true;
      }
      else {
        output_file << "\t";
      }
    }
    output_file << endl;
    i++;
  }
}


