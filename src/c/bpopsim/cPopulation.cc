#include "cPopulation.h"

void cPopulation::SetParameters(const variables_map &options)
{
	
  SetGrowthPhaseGenerations(
		options.count("generations-per-transfer") ?
		options["generations-per-transfer"].as<double>() : 6.64
		);
  SetPopSizeAfterDilution(
		options.count("population-size-after-transfer") ?
		options["population-size-after-transfer"].as<uint64_t>() : u_int64_t(5E6)
		); 
  SetInitialPopulationSize(
		options.count("initial-population-size") ?
		options["initial-population-size"].as<uint64_t>() : u_int64_t(2)
		);  
  SetMutationRatePerDivision(
		options.count("mutation-rate-per-division") ?
		options["mutation-rate-per-division"].as<double>() : 5E-8
		);  
  SetAverageMutationS(
		options.count("average-selection-coefficient") ?
		options["average-selection-coefficient"].as<double>() : 0.05
		);  
  SetTransferIntervalToPrint(
		options.count("transfer-interval-to-print") ?
	  options["transfer-interval-to-print"].as<int>() : 1
		);  
  SetVerbose(
		options.count("verbose") ?
		options["verbose"].as<int>() : 0
		);  
  SetTotalTransfers(
		options.count("number-of-transfers") ?
		options["number-of-transfers"].as<u_int64_t>() : 50
		);  
  SetMaxDivergenceFactor(
		options.count("marker-divergence") ?
		options["marker-divergence"].as<u_int64_t>() : 100
		);  
  SetReplicates(
		options.count("replicates") ?
		options["replicates"].as<int>() : 10
		);  
  SetMinimumPrinted(
		options.count("minimum-printed") ?
		options["minimum-printed"].as<u_int64_t>() : 8
		);
  SetBeneficialMutationDistribution(
		options.count("type-of-mutations") ?
		options["type-of-mutations"].as<char>() : 'u'
		);
  SetLineageTree(
		options.count("lineage-tree") ?
		options["lineage-tree"].as<int>() : 1
		);
	SetSeedParams(
		options.count("seed") ?
		options["seed"].as<long>() : 0
		);
	
  // Simulation parameters that are pre-calculated
  SetDilutionFactor(exp(log(2)*GetGrowthPhaseGenerations()));
  SetTransferBinomialSamplingP(1/GetDilutionFactor());
  SetPopSizeBeforeDilution(GetPopSizeAfterDilution() * GetDilutionFactor());
  SetLambda(1/GetMutationRatePerDivision());
  SetBinomialSamplingThreshold(1000);
}

void cPopulation::UpdateSubpopulations(long double update_time) 
{
  // @JEB note that m_divided_lineages is only valid when we assume our 
  // chunking is good such that each subpop can divide only once when mutation is happening
  m_divided_lineages.clear();
  
  m_population_size_stale = true; // we are changind the size of the population here.
  
  u_int64_t i=-1;
  for (std::vector<cSubpopulation>::iterator it = m_populations.begin(); it!=m_populations.end(); ++it) {
    i++; // must advance iterator before continue statement
    if (it->GetNumber() == 0) continue;

    // N = No * exp(log(2)*growth_rate * t) 
    double new_number = it->GetNumber() * exp(log(2) * update_time * it->GetFitness());     
    if (static_cast<uint32_t>(new_number) - static_cast<uint32_t>(it->GetNumber()) >= 1) {
      m_divided_lineages.push_back(i);
    }
    it->SetNumber(new_number);
  }
}

// @JEB for efficiency, we only recalculate this if it has been marked stale
const uint64_t cPopulation::GetPopulationSize() 
{
  if (!m_population_size_stale) return m_population_size;

  m_population_size = 0;
  for (std::vector<cSubpopulation>::iterator it = m_populations.begin(); it!=m_populations.end(); ++it) {
     m_population_size += static_cast<u_int64_t>(it->GetNumber());
  }  
  m_population_size_stale = false;
  return m_population_size;
}

long double cPopulation::TimeToNextWholeCell() 
{
   long double time_to_next_whole_cell = -1;
   for (std::vector<cSubpopulation>::iterator it = m_populations.begin(); it!=m_populations.end(); ++it) { 
   
      if(it->GetNumber() == 0) continue;

      //what is the time to get to the next whole number of cells?     
      long double current_cells = it->GetNumber();
      long double next_whole_cells = static_cast<u_int64_t>(it->GetNumber());
      
      // N = No * exp(log(2) * growth_rate * t) 
      long double this_time_to_next_whole_cell = (log(next_whole_cells / current_cells) / (it->GetFitness())) / log(2);   

      if ( time_to_next_whole_cell == -1 || (this_time_to_next_whole_cell < time_to_next_whole_cell) ) {
        time_to_next_whole_cell = this_time_to_next_whole_cell;
      }
  }
  
  return time_to_next_whole_cell;
}

//@agm Here I try to iterate through the populations in m_populations, request the size of 
//     each node, and divide the size of each node (subpopulation) by the total population size
//     This should give the relative frequence of a given unique_node_id in the population.

//@agm Now the information is stored in a vector and passed back to the main function for final printing.

void cPopulation::FrequenciesPerTransferPerNode(tree<cGenotype> newtree, 
                                                std::vector< std::vector<cGenotypeFrequency> >& frequencies)
{
	tree<cGenotype>::iterator update_location;
	double total_freqs(0);
	std::vector<cGenotypeFrequency> freq_per_node;
	
	u_int64_t total_cells(0);
	
	for(std::vector<cSubpopulation>::iterator it = m_populations.begin(); it!=m_populations.end(); ++it) {
		total_cells += it -> GetNumber();
	}
	
	freq_per_node.resize(newtree.size());
	std::vector<u_int64_t> number_per_subpop (newtree.size(),0);
	
	for(std::vector<cSubpopulation>::iterator it = m_populations.begin(); it!=m_populations.end(); ++it) {
		update_location = it -> GetGenotypeIter();
			
		while(update_location != NULL) {
			number_per_subpop[(*update_location).unique_node_id] += it -> GetNumber();
			update_location = newtree.parent(update_location);
		}
	}
	
	for(std::vector<cSubpopulation>::iterator it = m_populations.begin(); it!=m_populations.end(); ++it) {
		cGenotypeFrequency this_node;
		update_location = it -> GetGenotypeIter();
		this_node.unique_node_id = (*update_location).unique_node_id;
		this_node.frequency = (long double) number_per_subpop[this_node.unique_node_id]/total_cells;
		
		freq_per_node[this_node.unique_node_id] = this_node;
		total_freqs += this_node.frequency;
	}
	
	//Cout << Endl << "There are " << total_cells << " cells, in " << newtree.size() << " nodes."<< Endl;

	//@agm Printing sum of frequencies and building the doubly deep vector
	//Cout << Endl << "Sum of all freqs: " << total_freqs << Endl;
	frequencies.push_back(freq_per_node);
}


void cPopulation::Resample(gsl_rng * randgen) 
{
  //When it is time for a transfer, resample population
  
  uint32_t population_size_before_transfer = m_population_size;
  m_population_size_stale = true; // we are changind the size of the population here.
 
	 if (GetVerbose()) std::cout << ">> Transfer!" << std::endl;
	 //Is there an exists() in C++?
	 m_by_color[RED] = 0;
	 m_by_color[WHITE] = 0;
	
	 for (std::vector<cSubpopulation>::iterator it = m_populations.begin(); it!=m_populations.end(); ++it) {
			// Perform accurate binomial sampling only if below a certain population size
			if (it->GetNumber() < GetBinomialSamplingThreshold()) {
				 if (GetVerbose()) std::cout << "binomial " << it->GetNumber() << std::endl;
				 it->Transfer(GetTransferBinomialSamplingP(), randgen);
			}
			// Otherwise, treat as deterministic and take expectation...
			else {
				 it->SetNumber(it->GetNumber() * GetTransferBinomialSamplingP());
			}
					 
		  if (GetVerbose()) { 
				std::cout << it->GetMarker() << std::endl;
			  std::cout << it->GetNumber() << std::endl;
			  std::cout << it->GetFitness() << std::endl;
			}	
			if (it->GetMarker() == 'r') m_by_color[RED] += it->GetNumber();
			else m_by_color[WHITE] += it->GetNumber();
		  if (it->GetNumber() == 0) {
				SetTotalSubpopulationsLost(GetTotalSubpopulationsLost()+1);
				it = m_populations.erase(it);
				it--;
			}
	 }
	 //One color was lost, bail  
	 if ( (m_by_color[RED] == 0) || (m_by_color[WHITE] == 0) ) {
			SetKeepTransferring(false);
	 }     
	 if (GetVerbose()) std::cout << "Colors: " << m_by_color[RED] << " / " << m_by_color[WHITE] << std::endl;
	 SetRatio(m_by_color[RED] / m_by_color[WHITE]);
	 SetTransfers(GetTransfers()+1);

	 if ( /*(GetTransfers() >= 0) && */(GetTransfers() % GetTransferIntervalToPrint() == 0) ) {  
			m_this_run.push_back(GetRatio());
		 
			if (GetVerbose() == 1) { 
				 std::cout << "Transfer " << GetTransfers() << " : " << population_size_before_transfer << 
				 "=>" << GetPopulationSize() << "  R/W Ratio: " << GetRatio() << std::endl;  
				 std::cout << "Total mutations: " << GetTotalMutations() << " Maximum Fitness: " << GetMaxW() << std::endl;
				 std::cout << "Size = " << m_this_run.size() << std::endl;
			}
	 }  

	 if ( (u_int64_t(m_this_run.size()) >= GetMinimumPrinted()) && 
			 ((GetRatio() > GetMaxDivergenceFactor()) || 
				(GetRatio() < 1/GetMaxDivergenceFactor())) )   
	 {
			if (GetVerbose()) std::cout << "DIVERGENCE CONDITION MET" << std::endl;
			SetKeepTransferring(false);
	 }  
  
}

void cPopulation::PushBackRuns()
{
   m_runs.push_back(m_this_run);
}

void cPopulation::ClearRuns(cLineageTree& newtree)
{
   m_this_run.clear();
   m_populations.clear();
   newtree.clear();
}

void cPopulation::RunSummary()
{
   std::cout << "Total mutations: " << GetTotalMutations() << std::endl;
   std::cout << "Total subpopulations lost: " << GetTotalSubpopulationsLost() << std::endl;
   std::cout << "Transfers: " << GetTransfers() << std::endl;
   std::cout << "Maximum Fitness: " << GetMaxW() << std::endl;
}

void cPopulation::ResetRunStats()
{

   SetTotalMutations(0);    
   SetTotalSubpopulationsLost(0);
   SetTransfers(1);    
   SetDivisionsUntilMutation(0);
   SetKeepTransferring(true);

}

void cPopulation::DisplayParameters()
{
  if (GetVerbose()==1) {
    std::cout << "u = " << GetMutationRatePerDivision() << std::endl;
    std::cout << "s = " << GetAverageMutationS() << std::endl;
    std::cout << "N = " << GetPopSizeAfterDilution() << std::endl;
    std::cout << "dil = " << GetDilutionFactor() << std::endl;
  }
}

void cPopulation::CalculateDivisions()
{
  // Move time forward until another mutation occurs.
  // Calculate points to output between last time and current time
  // First time through the loop, a partial time interval
  // may be calculated to get back on $print_interval

  // Move forward by a large chunk of time which assumes
  // all populations have the maximum fitness in the population
  // This will, at worst, underestimate how long.
  // We can then move forward by single divisions to find the exact division where the mutation occurs

  // We would like to move forward by as many divisions as it takes to
  // get to either (1) the next mutation OR (2) the next transfer.
  long double desired_divisions = GetDivisionsUntilMutation();
  if (desired_divisions + GetPopulationSize() > GetPopSizeBeforeDilution()) {
     desired_divisions = GetPopSizeBeforeDilution() - GetPopulationSize();
  }
  
  if(GetVerbose() == 1) {
    std::cout << "Divisions before next mutation: " << GetDivisionsUntilMutation() <<std::endl;
  }
  
  // Note: we underestimate by a few divisions so that we can step forward by single division increments
  // as we get close to the one where the mutation happened (or right before a transfer).      
  
  if (desired_divisions < 1) {
     desired_divisions = 1;
  }
  
  if (GetVerbose() == 1) {
     std::cout << "Total pop size: " << GetPopulationSize() <<std::endl;
     std::cout << "Desired divisions " << desired_divisions <<std::endl;
  }
      
  // How much time would we like to pass to achieve the desired number of divisions?
  // (Assuming the entire population has the maximum fitness, makes us underestimate by a few)
  long double update_time = log((desired_divisions+(double)GetPopulationSize()) / (double)GetPopulationSize()) / (GetMaxW());

  // At a minumum, we want to make sure that one cell division took place

  // What is the minimum time required to get a single division?
  long double time_to_next_whole_cell = TimeToNextWholeCell();

  if (time_to_next_whole_cell > update_time) {
     if (GetVerbose())std::cout << "Time to next whole cell greater than update time: " << 
			 time_to_next_whole_cell << " < " << update_time <<std::endl;    
     update_time = time_to_next_whole_cell;
  }

  if (GetVerbose() == 1) {
     std::cout << "Update time: " << update_time <<std::endl;    
  }
            
  //Now update all lineages by the time that actually passed
 
  uint32_t previous_population_size = m_population_size;
  UpdateSubpopulations(update_time);
  SetCompletedDivisions(GetPopulationSize() - previous_population_size);
              
  if (GetVerbose())std::cout << "Completed divisions: " << GetCompletedDivisions() <<std::endl;
  SetDivisionsUntilMutation(GetDivisionsUntilMutation() - GetCompletedDivisions());
}

/*@agm The functions below should build a new tree using the tree.h header
       I thought the best way to do this was to comment out all of the previous code in both this file
       and the associated header so it would be clear where stuff was changed. */

void cPopulation::NewSeedSubpopulation(cLineageTree& newtree, 
                                       u_int64_t& node_id) 
{	
	cGenotype r, w;
	tree<cGenotype>::iterator top, red_side, white_side;
	long double starting_fitness = 1.0;
	
	//initialize object of cSubpopulation type
	cSubpopulation red, white;
	node_id = 0;
	
	/*@agm The unique code and fitness is set. */
	
	/*@agm This function leads to a really odd problem when incrementing from the AddSubpopulation function.
           Notice the red and white objects must by set before their the function can be called which increments
	       node_id... ahh but if the function call cannot be made without first assigning a node_id... irritating*/

	r.fitness = starting_fitness;
	r.unique_node_id = node_id;
    w.fitness = starting_fitness;
	w.unique_node_id = (node_id+1);
	
	//Start building the tree
	red_side = newtree.insert(newtree.begin(), r);
	white_side = newtree.insert(newtree.begin(), w);
	
	red.SetNumber(GetInitialPopulationSize()/2);
	red.SetGenotype(red_side);
	red.SetMarker('r');
	
	white.SetNumber(GetInitialPopulationSize()/2);
	white.SetGenotype(white_side);
	white.SetMarker('w');
	
	AddSubpopulation(red, node_id);
	AddSubpopulation(white, node_id);	
}

void cPopulation::AddSubpopulation(cSubpopulation& subpop, 
                                   u_int64_t& node_id) 
{
  m_population_size_stale = true; // we have just changed the population size
	m_populations.push_back(subpop);
	SetNumberOfSubpopulations(GetNumberOfSubpopulations()+1);
	node_id++;
	//std::cout << subpop.GetNode_id() << " " << subpop.GetFitness() << std::endl;
}

void cPopulation::NewMutate(gsl_rng * randgen, 
							cLineageTree& newtree, 
							u_int64_t& node_id) 
{	
	m_total_mutations++;
	
	if (m_verbose) std::cout << "* Mutating!" << std::endl;
	
	//Mutation happened in the one that just divided
	//Break ties randomly here.
	cSubpopulation& ancestor = m_populations[m_divided_lineages[rand() % m_divided_lineages.size()]];          
	cSubpopulation new_subpop;
  
  //std::cout << "Divided has number: " << ancestor.GetNumber() << std::endl;
  // There must be at least two cells for a mutation to have occurred...
  assert(ancestor.GetNumber() >= 2);
	
  new_subpop.NewCreateDescendant(randgen, 
                                 ancestor, 
								 GetAverageMutationS(), 
								 GetBeneficialMutationDistribution(), 
							     newtree, 
								 node_id);
	
	if (GetVerbose()) std::cout << "  Color: " << new_subpop.GetMarker() << std::endl;
	if (GetVerbose()) std::cout << "  New Fitness: " << new_subpop.GetFitness() << std::endl;
	
	AddSubpopulation(new_subpop, node_id);
	
	if(new_subpop.GetFitness() > GetMaxW()) 
	{
		SetMaxW(new_subpop.GetFitness());
	}
	
}

//@agm As the function name implies, this prints the frequnecies above some threshold to screen at whatever
//     time it is called and passed the frequencies vector

void cPopulation::PrintFrequenciesToScreen(std::vector< std::vector<cGenotypeFrequency> > frequencies) {
	Cout << "Done with round... Here's the Output:" << Endl << Endl;
	
	for (u_int64_t i = 0; i<frequencies.size(); i++) {
		double total_freqs = 0;
		for (u_int64_t j = 0; j<frequencies[i].size(); j++) {
			//@agm set up a minimum frequency to report the print out the number so it isn't overwhelming.
			if (frequencies[i][j].frequency > 0.001) {
				Cout << "Frequency of mutation # " << (frequencies[i][j]).unique_node_id << " at time " << i << " is: " << (frequencies[i][j]).frequency << Endl;
			}
			total_freqs += frequencies[i][j].frequency;
		}
		Cout << Endl << "Round # " << i << " sum of frequencies is: " << total_freqs << Endl << Endl;
	}
	Cout << "Number of cells before dilution: " << GetPopSizeBeforeDilution();
	Cout << Endl << "Number of cells after dilution: " << GetPopulationSize() << Endl << Endl;
	Cout << "------------------------------------------------" << Endl << Endl;
} 

//@agm I comandeered this function to print stuff out in the manner I see fit
void cPopulation::PrintOut(const std::string& output_file_name, 
                           std::vector< std::vector<cGenotypeFrequency> > frequencies)
{  
	//Print everything out
	std::ofstream output_file;
	output_file.open(output_file_name.c_str(),std::ios_base::app);
	
	//int last_time(frequencies.size()-1);

	/*for (int i = 0; i<frequencies[last_time].size(); i++) { 
		if (frequencies[last_time][i].frequency > .0001) output_file << std::setw(15) << frequencies[last_time][i].unique_node_id;
	}*/
	//int width = 20;
    
	for (u_int64_t i = 0; i<frequencies.size(); ++i) {
        //output_file.setf(std::ios::left);
        output_file << i;
		for (u_int64_t j = 0; j<frequencies[i].size(); ++j) {
			if (j == 0 || j == 1) {
                output_file << " " << frequencies[i][j].frequency;
            }
		}
        output_file << "\n";
	}
}

