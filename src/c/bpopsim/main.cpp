#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#if !defined(__SUNPRO_CC) || (__SUNPRO_CC > 0x530)
#include <boost/generator_iterator.hpp>
#endif

#include <boost/math/special_functions/factorials.hpp>

#include <boost/random/linear_congruential.hpp>
#include <boost/random.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>


using namespace std;
using namespace boost;
//methods
//--------
//transfer
//resample

struct Individual {
  double w;
  double n;
  char color;
};

// Global random number generator
mt19937 rng(time(NULL));

const int size_of_by_color = 2;

double returnExp(double a) 
{
    exponential_distribution<> exp_dist(a);
    variate_generator<mt19937&, exponential_distribution <> > next_value(rng, exp_dist);
    return next_value(); //Not sure why I need to do this...
}

uint64_t returnBin(const uint64_t n, const double p)
{
    binomial_distribution<> my_binomial(n,p);
    variate_generator<mt19937&, binomial_distribution<> > next_value(rng, my_binomial);
    return next_value();
}

// NOTE: The BOOST Poisson random number generator DOES NOT WORK for large means!!!
double returnPoisson(const double n, const double p)
{
    poisson_distribution<> my_poisson(n * p);
    variate_generator<mt19937&, poisson_distribution<> > next_value(rng, my_poisson);
    cout << "Poisson with mean " << n << " * " << p << " = " << next_value() << endl;
    return next_value();
}

int main()
{
	double update_time;

  // Simulation parameters that should be arguments
  double initial_population_size = 2;
	double pop_size_after_dilution = 5E6;             // N sub 0
	double mutation_rate_per_division = 1E-8;         // mu
	double average_mutation_s = 0.05;                 // s
  double growth_phase_generations = 6.64;
	double beneficial_mutation_distribution_code = 0; // 0 = uniform, 1 = exponential
  double binomial_sampling_threshold = 1000;
  string output_file_name("output.txt");
  
	//Output parameters that should be arguments
	int total_transfers = 200000;
	double max_divergence_factor = 100;
	double verbose = 0;  
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

    vector<Individual> populations; //Vector containing all of the populations;

		//Seed a red population
		Individual r;
		r.n = initial_population_size/2;
		r.w = 1;
		r.color = 'r';
    populations.push_back(r);

		//Seed a white population
		Individual w;
		w.n = initial_population_size/2;
    w.w = 1;
		w.color = 'w';
		populations.push_back(w);		
    
		double max_w = 1;
		uint64_t total_pop_size = initial_population_size;
		int total_mutations = 0;
		int total_subpopulations_lost = 0;
	
  
		int transfers = 1;
		double divisions_until_mutation = 0;
		bool keep_transferring = true;

		while ( (transfers < total_transfers) && keep_transferring )
		{
			// Move time forward until another mutation occurs.

			divisions_until_mutation += floor(returnExp(lambda));
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
        for (int i=0; i< populations.size(); i++) {
          Individual &this_subpop = populations[i];
          if (this_subpop.n == 0) continue;
        
          //what is the time to get to the next whole number of cells?
          double current_cells = this_subpop.n;
          double whole_cells = floor(this_subpop.n)+1;
          // WC = N * exp(growth_rate * t) 
          
          double this_time_to_next_whole_cell = log(whole_cells / current_cells) / (this_subpop.w);
          
          if ( time_to_next_whole_cell == 0 || (this_time_to_next_whole_cell < time_to_next_whole_cell) ) 
          {
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
        uint64_t new_pop_size = 0;
        
        for (vector<Individual>::iterator it = populations.begin(); it!=populations.end(); ++it) {
          if (it->n == 0) continue;
          it->n = it->n * exp(log(2) * update_time * it->w);
          new_pop_size += it->n;
        }				
                
        uint64_t completed_divisions = new_pop_size - total_pop_size;
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
          Individual ancestor = populations[divided_lineages[rand() % divided_lineages.size()]];
          
          //What is the new lineage's fitness?
          //Option #1: dirac delta
          double new_w = 0;
          if (beneficial_mutation_distribution_code == 0) {
            new_w = ancestor.w + average_mutation_s;
          }
          //Option #2: exponential
          else {
            double this_mutation_s = returnExp(average_mutation_s);
            new_w = ancestor.w + this_mutation_s;
          }
          
          if (verbose) cout << "  Color: " << ancestor.color << endl;
          if (verbose) cout << "  New Fitness: " << new_w << endl;
          
          //Create and add the new lineage

          Individual new_lineage;
          new_lineage.n = 1;
          new_lineage.w = new_w;
          new_lineage.color = ancestor.color;

          populations.push_back(new_lineage);
          
          //Update maximum fitness
          if(new_w > max_w) {
            max_w = new_w;
          }
          
          //One ancestor organism was converted to the new lineage
          ancestor.n--;
        }
      
        //When it is time for a transfer, resample population
        if (total_pop_size >= pop_size_before_dilution) {
        
          if (verbose) cout << ">> Transfer!" << endl;
          //Is there an exists() in C++?

          double by_color[size_of_by_color];
          by_color[0] = 0;
          by_color[1] = 0;
          
          double new_pop_size = 0;
          //POPULATION: for (my $i=0; $i< scalar @populations; $i++) 	
          for (vector<Individual>::iterator it = populations.begin(); it!=populations.end(); ++it) {
            if (it->n == 0) continue;
            
            // Perform accurate binomial sampling only if below a certain population size
            if (it->n < binomial_sampling_threshold) {
              if (verbose) cout << "binomial" << it->n << endl;
              it->n = returnBin(it->n, transfer_binomial_sampling_p);
            }
            // Otherwise, treat as deterministic and take expectation...
            else {
              it->n *= transfer_binomial_sampling_p;
            }
            
            // Keep track of lineages we lost
            if (it->n == 0) {
              total_subpopulations_lost++;
            }
            

            new_pop_size += floor(it->n);
            //There is probably a better way to do this, I just don't know the syntax ??by_color[p.color] += p.n;??
            if (verbose) cout << it->color << endl;
            if (verbose) cout << it->n << endl;
            if (verbose) cout << it->w << endl;
            if (it->color == 'r') {
              by_color[0] += it->n;
            }
            else {
              by_color[1] += it->n;
            }
          }
          //One color was lost, bail	
          if ( (by_color[0] == 0) || (by_color[1] == 0) ) {
            keep_transferring = false;
          }
            
          if (verbose) cout << "Colors: " << by_color[0] << " / " << by_color[1] << endl;
          double ratio = by_color[0] / by_color[1];
          
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

          if ( (this_run.size() >= minimum_printed) && ((ratio > max_divergence_factor) || (ratio < 1/max_divergence_factor)) ) {
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
      if(i<runs[on_run].size()) {
        
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


