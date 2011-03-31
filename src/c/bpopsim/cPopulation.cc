#include "cPopulation.h"

void cPopulation::SetParameters(const variables_map &options)
{
	
  SetGrowthPhaseGenerations(
		options.count("generations-per-transfer") ?
		options["generations-per-transfer"].as<double>() : 6.64
		);
  SetPopSizeAfterDilution(
		options.count("population-size-after-transfer") ?
		options["population-size-after-transfer"].as<uint32_t>() : uint32_t(5E6)
		); 
  SetInitialPopulationSize(
		options.count("initial-population-size") ?
		options["initial-population-size"].as<uint32_t>() : uint32_t(2)
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
	  options["transfer-interval-to-print"].as<uint16_t>() : 1
		);  
  SetVerbose(
		options.count("verbose") ?
		options["verbose"].as<uint16_t>() : 0
		);  
  SetTotalTransfers(
		options.count("number-of-transfers") ?
		options["number-of-transfers"].as<uint16_t>() : 50
		);  
  SetMaxDivergenceFactor(
		options.count("marker-divergence") ?
		options["marker-divergence"].as<uint16_t>() : 100
		);  
  SetReplicates(
		options.count("replicates") ?
		options["replicates"].as<uint16_t>() : 10
		);  
  SetMinimumPrinted(
		options.count("minimum-printed") ?
		options["minimum-printed"].as<uint16_t>() : 8
		);
  SetBeneficialMutationDistribution(
		options.count("type-of-mutations") ?
		options["type-of-mutations"].as<char>() : 'u'
		);
  SetLineageTree(
		options.count("lineage-tree") ?
		options["lineage-tree"].as<uint16_t>() : 1
		);
	SetSeedParams(
		options.count("seed") ?
		options["seed"].as<uint16_t>() : 0
		);
  SetRedWhiteOnly(
    options.count("redwhite-only") ?
    options["redwhite-only"].as<char>() : 'f'
    );
  SetLogApproximation(
    options.count("log-approximation") ?
    options["log-approximation"].as<char>() : 'f'
    );
  SetLogApproximationValue(
    options.count("log-approximation-value") ?
    options["log-approximation-value"].as<int>() : 14
    );
  
  // Simulation parameters that are pre-calculated
  SetDilutionFactor(exp(log(2)*GetGrowthPhaseGenerations()));
  SetTransferBinomialSamplingP(1/GetDilutionFactor());
  SetPopSizeBeforeDilution(GetPopSizeAfterDilution() * GetDilutionFactor());
  SetLambda(1/GetMutationRatePerDivision());
  SetBinomialSamplingThreshold(1000);
}

void cPopulation::UpdateSubpopulations(double update_time) 
{
  // @JEB note that m_divided_lineages is only valid when we assume our 
  // chunking is good such that each subpop can divide only once when mutation is happening
  m_divided_lineages.clear();
  
  m_population_size_stale = true; // we are changind the size of the population here.
  
  uint32_t i=-1;
  for (std::vector<cSubpopulation>::iterator it = m_populations.begin(); it!=m_populations.end(); ++it) {
    i++; // must advance iterator before continue statement
    if (it->GetNumber() == 0) continue;

    // N = No * exp(log(2)*growth_rate * t) 
    //@agm This is the slow use of the log and exp functions
    double new_number = it->GetNumber() * exp(ReturnLog(2) * update_time * it->GetFitness());     
    if (static_cast<uint32_t>(new_number) - static_cast<uint32_t>(it->GetNumber()) >= 1) {
      m_divided_lineages.push_back(i);
    }
    it->SetNumber(new_number);
  }
}

// @JEB for efficiency, we only recalculate this if it has been marked stale
const uint32_t cPopulation::GetPopulationSize() 
{
  if (!m_population_size_stale) return m_population_size;

  m_population_size = 0;
  for (std::vector<cSubpopulation>::iterator it = m_populations.begin(); it!=m_populations.end(); ++it) {
     m_population_size += static_cast<uint32_t>(it->GetNumber());
  }  
  m_population_size_stale = false;
  return m_population_size;
}

double cPopulation::TimeToNextWholeCell() 
{
   double time_to_next_whole_cell = -1;
   for (std::vector<cSubpopulation>::iterator it = m_populations.begin(); it!=m_populations.end(); ++it) { 
   
      if(it->GetNumber() == 0) continue;

      //what is the time to get to the next whole number of cells?     
      double current_cells = it->GetNumber();
      double next_whole_cells = static_cast<uint32_t>(it->GetNumber());
      
      // N = No * exp(log(2) * growth_rate * t) 
      double this_time_to_next_whole_cell = (ReturnLog(next_whole_cells / current_cells) / (it->GetFitness())) / ReturnLog(2);   

      if ( time_to_next_whole_cell == -1 || (this_time_to_next_whole_cell < time_to_next_whole_cell) ) {
        time_to_next_whole_cell = this_time_to_next_whole_cell;
      }
  }
  
  return time_to_next_whole_cell;
}

//@agm Here I iterate through the populations in m_populations, build a vector to store the values 
//     for the sizes of each node by iterating up, and divide the size of each node (subpopulation) by the total 
//     population size. This should give the relative frequency of a given unique_node_id in the population.

//@agm Now the information is stored in a vector and passed back to the main function for later use.

void cPopulation::FrequenciesPerTransferPerNode(tree<cGenotype> * newtree, 
                                                std::vector< std::vector<cGenotypeFrequency> > * frequencies)
{
	tree<cGenotype>::iterator update_location;
	uint32_t total_cells(0);
  double total_freqs;
	
	std::vector<cGenotypeFrequency> freq_per_node(newtree->size());
	std::vector<uint32_t> number_per_subpop(newtree->size(),0), node_ids(newtree->size(),0);
  std::vector<int> subpops(newtree->size(),0);
  
  for (std::vector<cSubpopulation>::iterator it = m_populations.begin(); it!=m_populations.end(); ++it) {
    update_location = it -> GetGenotypeIter();
    
    subpops[(*update_location).unique_node_id] = it->GetNumber();
    node_ids[(*update_location).unique_node_id] = (*update_location).unique_node_id;
    
    total_cells += it->GetNumber();
  
		while(update_location != NULL) {
			number_per_subpop[(*update_location).unique_node_id] += it -> GetNumber();
      update_location = newtree->parent(update_location);
		}
	}
  
  m_total_cells.push_back(total_cells);
  m_subpops.push_back(subpops);
  
  /*for (std::vector<uint32_t>::iterator it = number_per_subpop.begin(); it!=number_per_subpop.end(); ++it) {
    Cout << *it << " ";
  }*/
  
	
  //@agm It is critical to iterate over the number_per_subpop vector here rather than the m_populations
  //     iterating over m_populations will cause you to lose populations after 350 or so transfers
  //     This is necessary because eventually there are no cells left who ONLY have red or white causing
  //     them to be dropped.
  //for (std::vector<cGenotypeFrequency>::iterator it = freq_per_node.begin(); it!=freq_per_node.end(); ++it) total_cells += (*it).subpop_size;
  //std::cout << total_cells << " " << number_per_subpop[0] << std::endl;
  int count(0);
	for (std::vector<uint32_t>::iterator it = node_ids.begin(); it!=node_ids.end(); ++it) {
		cGenotypeFrequency this_node;
    
    this_node.unique_node_id = (*it);
    this_node.frequency = (double) number_per_subpop[(*it)]/number_per_subpop[0];
    freq_per_node[(*it)] = this_node;
		
		total_freqs += this_node.frequency;
    //Cout << Endl << this_node.unique_node_id << "  " << this_node.frequency;
    count++;
	}

	//Cout << Endl << "There are " << total_cells << " cells, in " << newtree.size() << " nodes."<< Endl;

	//@agm Printing sum of frequencies and building the doubly deep vector
	//Cout << Endl << "Sum of all freqs: " << total_freqs << Endl;
	frequencies->push_back(freq_per_node);
}

//@agm Here I want to calculate the frequencies of subpopulations rather than mutations
//     It's not yet clear to me how to do that.

void cPopulation::AssignChildFreq(tree<cGenotype>::sibling_iterator child_node,
                                  std::vector<cGenotypeFrequency> * frequencies,
                                  tree<cGenotype> * newtree,
                                  double parent_low,
                                  double parent_high,
                                  std::vector<cChildFrequency> * child_freqs) {
  
  //kptree::print_tree_bracketed(*newtree);
  double child_low(parent_low), child_high;
  
  child_high = child_low + (*frequencies)[(*child_node).unique_node_id].frequency;
  
  (*child_freqs)[(*child_node).unique_node_id].child_high = child_high;
  (*child_freqs)[(*child_node).unique_node_id].child_low = child_low;
  
  //std::cout << (*child_node).unique_node_id << " " << child_high << " " << child_low << std::endl;
  
  for (tree<cGenotype>::sibling_iterator node = (*newtree).begin(child_node); node!=(*newtree).end(child_node); ++node) {
    AssignChildFreq(node, frequencies, newtree, child_low, child_high, child_freqs);
    child_low = child_high;
  }
  //std::cout << (*child_node).unique_node_id << " " << child_high << " " << child_low << std::endl;
  
}

void cPopulation::DrawMullerMatrix(tree<cGenotype> * newtree, 
                                   std::vector< std::vector<int> > muller_matrix, 
                                   std::vector< std::vector<cGenotypeFrequency> > * frequencies){
  
  std::vector< tree<cGenotype>::iterator > where(newtree->size());
  
  //double threshold(.025);
  //std::vector<bool> relevant_columns(newtree->size(), true);
  
  //Build vector called where to store all iterator in tree for parental recal later
  for (tree<cGenotype>::iterator node_loc = (*newtree).begin(); node_loc != (*newtree).end(); node_loc++)
    where[(*node_loc).unique_node_id] = node_loc;
  
  std::string filename = "/Users/austin/Desktop/test_muller.out";
  std::ofstream output_handle(filename.c_str());
  
  //step through simulation time
  for (uint32_t time=0; time<(*frequencies).size(); time++) {
    //int total_cells_per_time(m_total_cells[time]);
    //std::vector<int> this_time_point(1/threshold, 0);
    
    //int count(0), mutation_counter(0);
    
    std::vector<cChildFrequency> child_freqs((*newtree).size());
    
    tree<cGenotype>::sibling_iterator location;
    location = (*newtree).begin();
    
    AssignChildFreq(location, &((*frequencies)[time]), newtree, 0, (*frequencies)[time][0].frequency, &child_freqs);
    
    //Iterator and pointer horror!!!!!!!!!
    /*for (std::vector<cGenotypeFrequency>::iterator it = (*frequencies)[time].begin(); it!=(*frequencies)[time].end(); ++it) {
      //skip to end of this time step if most of this time point will be zeros
      if( relevant_columns[(*it).unique_node_id] == true ) {
        //std::cout << relevant_columns[(*it).unique_node_id] << std::endl;
        if( (*it).unique_node_id == 0 ) {
          double portion; 
          portion = (double) m_subpops[time][count] / total_cells_per_time;
          if( portion == 1 ) goto skip_zeros;
        }
        else {
          double portion; 
          portion = (double) m_subpops[time][count] / total_cells_per_time;
          int number_of_spots (int(portion*this_time_point.size()));
    
          //find parent of it
          cGenotype node_parent;
          node_parent = ( *(newtree->parent(where[(*it).unique_node_id])) );
        
          //std::cout << (*it).unique_node_id << " " << portion << " " << (*it).subpop_size << " " << total_cells_per_time << " " << node_parent.unique_node_id << std::endl;
          
          for (int j=0; j<number_of_spots; j++) {
            for (int i=0; i<this_time_point.size(); i++) {
              if( node_parent.unique_node_id == this_time_point[i]) {
                this_time_point[i] = mutation_counter;
                break;
              }
            }
          }
          if( number_of_spots > 0 ) mutation_counter++;;
        }
      }
      count++;
    }*/
    
    skip_zeros:
    std::cout << time << std::endl;
    for(int i=0; i<child_freqs.size(); i++) {
      if( ((child_freqs[i].child_high) - (child_freqs[i].child_low)) > .01 ) output_handle << std::left << std::setw(8) << child_freqs[i].child_low << std::setw(8) << child_freqs[i].child_high << std::endl;
    }
    output_handle << std::endl;
    //muller_matrix.push_back(this_time_point);
  }
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

	 if ( (uint32_t(m_this_run.size()) >= GetMinimumPrinted()) && 
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
  double desired_divisions = GetDivisionsUntilMutation();
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
  // Can't change this log to ReturnLog with log table or nothing works
  double update_time = log((desired_divisions+(double)GetPopulationSize()) / (double)GetPopulationSize()) / (GetMaxW());

  // At a minumum, we want to make sure that one cell division took place

  // What is the minimum time required to get a single division?
  double time_to_next_whole_cell = TimeToNextWholeCell();

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

/*@agm The functions below should build a new tree using the tree.h header */

void cPopulation::SeedSubpopulationForRedWhite(tree<cGenotype>* newtree, 
                                               uint32_t &node_id) 
{	
	cGenotype r, w;
	tree<cGenotype>::iterator top, red_side, white_side;
	double starting_fitness = 1.0;
	
	//initialize object of cSubpopulation type
	cSubpopulation red, white;
	node_id = 0;
	
	/*@agm The unique code and fitness is set. */

	r.fitness = starting_fitness;
	r.unique_node_id = node_id;
  w.fitness = starting_fitness;
	w.unique_node_id = (node_id+1);
	
	//Start building the tree
	red_side = newtree->insert(newtree->begin(), r);
	white_side = newtree->insert(newtree->begin(), w);
	
	red.SetNumber(GetInitialPopulationSize()/2);
	red.SetGenotype(red_side);
	red.SetMarker('r');
	
	white.SetNumber(GetInitialPopulationSize()/2);
	white.SetGenotype(white_side);
	white.SetMarker('w');
	
	AddSubpopulation(red, r.unique_node_id);
	AddSubpopulation(white, w.unique_node_id);	
}

//@agm This function seeds the population with only one colony to avoid the red/white problem
//     when possible. For this to work properly, I need to get rid of victory conditions.

//update: Victory conditions are now a command line option.  If no set to false the simulation
//        will run to the max number of iterations.

void cPopulation::SeedPopulationWithOneColony(tree<cGenotype>* newtree, 
                                              uint32_t &node_id) {
  cGenotype neutral;
  tree<cGenotype>::iterator start_position;
  double starting_fitness(1.0);
  
  cSubpopulation begin_here;
  node_id = 0;
  
  neutral.fitness = starting_fitness;
  neutral.unique_node_id = node_id;
  
  start_position = newtree->insert(newtree->begin(), neutral);
  
  begin_here.SetNumber(GetInitialPopulationSize());
  begin_here.SetGenotype(start_position);
  begin_here.SetMarker('n');
  
  AddSubpopulation(begin_here, neutral.unique_node_id);
  
}

void cPopulation::AddSubpopulation(cSubpopulation& subpop, 
                                   uint32_t& node_id) 
{
  m_population_size_stale = true; // we have just changed the population size
	m_populations.push_back(subpop);
	SetNumberOfSubpopulations(GetNumberOfSubpopulations()+1);
	node_id++;
	//std::cout << subpop.GetNode_id() << " " << subpop.GetFitness() << std::endl;
}

//Generates new mutant and adds it to the tree
void cPopulation::Mutate(gsl_rng * randgen, 
                            tree<cGenotype> * newtree, 
                            uint32_t & node_id) 
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

//Utilities Section

//@agm I comandeered this function to print stuff out in the manner I see fit
//@agm I basically do a non-human readable raw dump because it's easier for R to deal with
//     If you want human readable use the PrintToScreen function

/***** Important for post-hoc use with R ******
    The output file will have a header and it
    will have the largest number of columns in file.
    Thus, in R, you can simply set header to 
    true, and fill NA spaces as the read.table
    function allows. */

void cPopulation::PrintOut(const std::string& output_file_name, 
                           std::vector< std::vector<cGenotypeFrequency> > * frequencies)
{  
  std::vector<bool> ColsToPrint = MutationAboveThreshold(&(*frequencies), .1);
  
	//Print everything out
	std::ofstream output_file;
	output_file.open(output_file_name.c_str(),std::ios_base::app);
  
  uint32_t last_time((*frequencies).size()-1);
  
	for (uint32_t i = 0; i<(*frequencies)[last_time].size(); i++) { 
    if ( ColsToPrint[i] == true ) output_file << i << " ";
	}
	//int width = 20;
  output_file << "\n";
  
	for (uint32_t i = 0; i<(*frequencies).size(); i++) {
    uint32_t count(0);
    for ( std::vector<cGenotypeFrequency>::iterator it = (*frequencies)[i].begin(); it!=(*frequencies)[i].end(); ++it) {
      if ( ColsToPrint[count] == true ) output_file << (*it).frequency << " ";
      count++;
		}
    output_file << "\n";
	}
}

void cPopulation::ClearRuns(tree<cGenotype>* newtree)
{
  m_this_run.clear();
  m_populations.clear();
  newtree->clear();
}

//@agm As the function name implies, this prints the frequencies above some threshold to screen at whatever
//     time it is called and passed the frequencies vector

void cPopulation::PrintFrequenciesToScreen(std::vector< std::vector<cGenotypeFrequency> > * frequencies) {
	Cout << "Done with round... Here's the Output:" << Endl << Endl;
  
  int count(0);
  
	for (uint32_t i = 0; i<(*frequencies).size(); i++) {
		double total_freqs = 0;
    for ( std::vector<cGenotypeFrequency>::iterator it = (*frequencies)[i].begin(); it!=(*frequencies)[i].end(); ++it) {
    //@agm set up a minimum frequency to report the print out the number so it isn't overwhelming.
      if ((*it).frequency > 0.01) {
        Cout << "Frequency of mutation # " << std::right << std::setw(6) << (*it).unique_node_id << " at time " << std::right << std::setw(4) << i << " is: " << std::left << std::setw(10) << (*it).frequency << Endl;
      }
      total_freqs += (*it).frequency;
      count++;
    }
		Cout << Endl << "Round # " << i << " sum of frequencies is: " << total_freqs << Endl << Endl;
	}
	Cout << "Number of cells before dilution: " << GetPopSizeBeforeDilution();
	Cout << Endl << "Number of cells after dilution: " << GetPopulationSize() << Endl << Endl;
	Cout << "------------------------------------------------" << Endl << Endl;
} 

//@agm This function determines the maximum difference in genotype frequency between a mutation
//     and its predecessor if they got above the threshold passed to the threshold function below

void cPopulation::CalculateSimilarity(std::vector< std::vector<cGenotypeFrequency> > * frequencies) {
  std::vector<bool> relevant_mutations (MutationAboveThreshold(frequencies, .1));
  std::vector< std::vector<cGenotypeFrequency> > only_relevant_mutations;
  int counter(0);
  
  for (uint32_t i = 0; i<(*frequencies).size(); i++) {
    std::vector<cGenotypeFrequency> relevant_mutations_per_time(relevant_mutations.size());
    uint32_t count(0);
    counter = 0;
    for ( std::vector<cGenotypeFrequency>::iterator it = (*frequencies)[i].begin(); it!=(*frequencies)[i].end(); ++it) {
      if ( relevant_mutations[count] == true ) {
        relevant_mutations_per_time[counter] = (*it);
        counter++;
      }
      count++;
		}
    only_relevant_mutations.push_back(relevant_mutations_per_time);
	}
  
  std::vector<float> max_diff(counter, 0);
  float current_diff;
  
  for (int i = 0 ; i < counter; i++) {
    for (int time = 0; time < only_relevant_mutations.size(); time++) {
      current_diff = fabs(only_relevant_mutations[time][i].frequency - only_relevant_mutations[time][i+1].frequency);
      if( max_diff[i] < current_diff ) max_diff[i] = current_diff;
    }
  }
  for (int i = 0; i<max_diff.size(); i++) std::cout << std::endl << i << " " << max_diff[i] << std::endl;
}

//@agm This function takes the frequency pointer and returns a boolean vector
//     The boolean vector contains a true in the mutations that got above the passed threshold

std::vector<bool> cPopulation::MutationAboveThreshold(std::vector< std::vector<cGenotypeFrequency> > * frequencies, float threshold) {
  uint32_t last_time((*frequencies).size()-1);
  
  std::vector<bool> Fixed((*frequencies)[last_time].size(),false);
  
  for (uint16_t i = 0; i<(*frequencies).size(); i++) {
    uint32_t count(0);
    for (std::vector<cGenotypeFrequency>::iterator it = (*frequencies)[i].begin(); it!=(*frequencies)[i].end(); ++it) {
      if ( (*it).frequency >= threshold ) Fixed[count] = true;
      count++;
    }
  }
  return Fixed;
}

//@agm I wrote this based on the taylor series expansion... aren't we impressed
//     unfortunately it is significantly slower than the built in log function
//     I kept it around to test the effect of various precisions on the simulation
//     To get better precision change the iterations to whatever number you like

float cPopulation::Logarithm(float mantissa) {
  float value(0), num, topower;
  uint8_t iterations(2);
  
  for (uint8_t i = 0; i<iterations; i++) {
    topower = 1;
    num = (2*i)+1;
    for(uint8_t j = 0; j<num; j++) topower *= ((mantissa-1)/(mantissa+1));
    value += (1/num)*topower;
  }
  return 2*value;
}

//Building lookup table for use with icsilog
void cPopulation::ConstructLookUpTable() {
  m_lookuptable = (float*) malloc(((int) std::pow(exp(1),m_N))*sizeof(float));
  fill_icsi_log_table2(m_N, m_lookuptable);
}

/* ICSIlog V 2.0 */
void cPopulation::fill_icsi_log_table2(const unsigned precision, float* const   pTable)
{
  /* step along table elements and x-axis positions
   (start with extra half increment, so the steps intersect at their midpoints.) */
  float oneToTwo = 1.0f + (1.0f / (float)( 1 <<(precision + 1) ));
  int i;
  for(i = 0;  i < (1 << precision);  ++i )
  {
    // make y-axis value for table element
    pTable[i] = logf(oneToTwo) / log(2);
    
    oneToTwo += 1.0f / (float)( 1 << precision );
  }
}

double cPopulation::ReturnLog(double num) {
  if( m_approx_bool == 't' || m_approx_bool == 'T') return icsi_log_v2(num, m_lookuptable, m_N);
  else return log(num);
}