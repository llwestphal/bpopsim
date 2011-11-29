#include "common.h"
#include "cPopulation.h"

using namespace bpopsim;
using namespace std;

void cPopulation::SetParameters(AnyOption &options)
{
	
  SetGrowthPhaseGenerations(
		from_string<double>(options["generations-per-transfer"])
		);
  SetPopSizeAfterDilution(
		from_string<double>(options["population-size-after-transfer"])
		); 
  SetInitialPopulationSize(
    from_string<double>(options["initial-population-size"])
		);  
  SetMutationRatePerDivision(
		from_string<double>(options["mutation-rate-per-division"])
		);
  SetTotalTransfers(
		from_string<uint32_t>(options["number-of-transfers"])
		);  
  SetMaxDivergenceFactor(
		from_string<double>(options["marker-divergence"])
		);  
  SetBeneficialMutationDistribution(
		from_string<char>(options["type-of-mutations"])
		);
  SetInitialMutVals(
    from_string<vector< double > >(options["imv"])
    );
  SetCoarseGraining(
    from_string<uint16_t>(options["coarse-graining"])
    );
  SetInitialFitness(
    from_string<double>(options["initial-fitness"])
    );
  SetMullerRez(
    from_string<uint64_t>(options["muller_res"])
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
  
  uint32_t i=-1;
  m_population_size = 0; // Update the population size.
  for (std::vector<cSubpopulation>::iterator it = m_current_subpopulations.begin(); it!=m_current_subpopulations.end(); ++it) {
    i++; // must advance iterator before continue statement
    if (it->GetNumber() == 0) continue;

    // N = No * exp(log(2)*growth_rate * t) 
    //@agm This is the slow use of the log and exp functions
    double new_number = it->GetNumber() * exp(log(2) * update_time * it->GetFitness());     
    if (floor(new_number) - floor(it->GetNumber()) >= 1) {
      m_divided_lineages.push_back(i);
    }
    it->SetNumber(new_number);
    m_population_size += floor(new_number);
  }
}

// @JEB Calculates the total population size by iterating over subpops. 
// This is for checking the code only. Normally, use GetPopulationSize().
const double cPopulation::CalculatePopulationSize() 
{
  double calculated_population_size = 0;
  for (std::vector<cSubpopulation>::iterator it = m_current_subpopulations.begin(); it!=m_current_subpopulations.end(); ++it) {
     calculated_population_size += static_cast<double>(it->GetNumber());
  }  
  return calculated_population_size;
}

double cPopulation::TimeToNextWholeCell() 
{
   double time_to_next_whole_cell = -1;
   for (std::vector<cSubpopulation>::iterator it = m_current_subpopulations.begin(); it!=m_current_subpopulations.end(); ++it) { 
   
      if(it->GetNumber() == 0) continue;

      //what is the time to get to the next whole number of cells?     
      double current_cells(it->GetNumber());
      double next_whole_cells( floor(it->GetNumber())+1 );
      
      //cout << "Current cells: " << current_cells << " Next whole cells: " << next_whole_cells << endl;
      // N = No * exp(log(2) * growth_rate * t) 
      double this_time_to_next_whole_cell = (log(next_whole_cells / current_cells) / (it->GetFitness())) / log(2);   

      if ( time_to_next_whole_cell == -1 || (this_time_to_next_whole_cell < time_to_next_whole_cell) ) {
        time_to_next_whole_cell = this_time_to_next_whole_cell;
      }
  }
  
  return time_to_next_whole_cell;
}

//@agm Here I iterate through the populations in m_current_subpopulations, build a vector to store the values 
//     for the sizes of each node by iterating up, and divide the size of each node (subpopulation) by the total 
//     population size. This should give the relative frequency of a given unique_node_id in the population.

//@agm Now the information is stored in a vector and passed back to the main function for later use.

void cPopulation::FrequenciesPerTransferPerNode()
{	
  //  Debug: this checks to see if our population size was correct.
  //  Please leave commented out and don't remove.
  //uint32_t current_population_size = GetPopulationSize();
  //uint32_t calculated_population_size = CalculatePopulationSize();
  //assert(current_population_size == calculated_population_size);
  
  // There should be one node in the tree for each assigned genotype id.
  // Since these are the same size, use the one that doesn't require the
  // tree to calculate its size for speed (m_genotype_count). @JEB
  
  //@agm Now we cull nodes
  //assert(m_tree.size() == m_genotype_count);
  
  // Genotype frequencies (mutations at each node of tree)
	vector<cGenotypeFrequency> freq_per_node(m_genotype_count);
  
  // Genotype number (mutations at each node of tree)
	vector<double> number_with_genotype( m_genotype_count, 0.0 );
  
  m_cell_equivalents = 0;

	tree<cGenotype>::iterator update_location;
  for (uint32_t it = 0; it<m_current_subpopulations.size(); ++it) {
    
    cSubpopulation& this_pop(m_current_subpopulations[it]);
    
    update_location = this_pop.GetGenotypeIter();
      
    // traverse up the tree to the first ancestor, adding the number in
    // this subpopulation to every mutational node along the way
		while(update_location != NULL) {
			number_with_genotype[update_location->unique_node_id] += this_pop.GetNumber()*this_pop.GetTimer();
      update_location = m_tree.parent(update_location);
		}
    
    m_cell_equivalents += this_pop.GetNumber()*this_pop.GetTimer();
	}
  
  //cout << m_population_size << endl;
  m_total_cells.push_back(m_cell_equivalents); // shouldn't need to save this, remove later @JEB
	
	for (uint32_t i=0; i < m_genotype_count; i++) {
    
    if( number_with_genotype[i] > 0 ) {
      cGenotypeFrequency this_node;
      
      this_node.unique_node_id = i;
      this_node.frequency = number_with_genotype[i]/m_cell_equivalents;
      
      freq_per_node[i] = this_node;
    }
    else {
      cGenotypeFrequency this_node;
      
      this_node.unique_node_id = i;
      this_node.frequency = 0.0;
      
      freq_per_node[i]=this_node;
    }

	}
	m_frequencies.push_back(freq_per_node);
}

void cPopulation::CalculateAverageFitness() {
  double total_fitness(0);
  uint16_t pop_counter(0);
  
  for (std::vector<cSubpopulation>::iterator it = m_current_subpopulations.begin(); it!=m_current_subpopulations.end(); ++it) {    
    total_fitness += it->GetFitness();
    pop_counter++;
  }
  m_average_fitness.push_back(total_fitness/pop_counter);
}

void cPopulation::PrintFitness(std::string output_folder) {
  std::string output_file;
  std::ofstream output_handle;
  
  output_file = output_folder + "/AverageFitness.dat";
  
	output_handle.open(output_file.c_str(),std::ios_base::app);
  
  for (uint16_t time = 0; time<m_average_fitness.size(); time++) {
    if( time%m_coarse_graining == 0 ) {
      output_handle << m_average_fitness[time] << "\t";
      std::cout << m_average_fitness[time] << "\t";
    }
  }
  output_handle << "\n";
  output_handle.close();
}

void cPopulation::PrintSingleFitness(std::string output_folder) {
  std::string output_file;
  std::ofstream output_handle;
  
  float single_fitness;
  
  output_file = output_folder + "/SingleFitness.dat";
  output_handle.open(output_file.c_str(),std::ios_base::app);
  
  if(m_current_subpopulations.size() == 1) {
    for (std::vector<cSubpopulation>::iterator it = m_current_subpopulations.begin(); it!=m_current_subpopulations.end(); ++it) { 
      single_fitness = it->GetFitness();
    }
  }
  else {
    std::cerr << "There is more than one population remaining... there is a bug." << std::endl;
    abort();
  }
  
  output_handle << single_fitness << "\n";
  //if( single_fitness != 1 )
    //std::cout << "Fitness: " << single_fitness << "\n";
  output_handle.close();
}

void cPopulation::ConvertExternalData(const string &input_file) {
  
  fstream input_handle;
  input_handle.open(input_file.c_str(), fstream::in);
  assert(input_handle.is_open());
  
  char char_line[1000];
  
  uint16_t num_nodes = 1;
  
  while(!input_handle.getline(char_line,1000).eof()) {
    string line(char_line);
    vector<string> line_pieces(split(line, " ")), preformat_tree(split(line_pieces[0], ","));
    
    tree<cGenotype>::iterator top, next;
    
    if( num_nodes == 1 ) { 
      cGenotype origin;
      origin.unique_node_id = num_nodes;
      origin.name = preformat_tree[0];
      origin.fitness = 1;
      origin.mut_num = 0;
      
      top = m_tree.begin();
      m_tree.insert(top, origin);
      
      m_frequencies.resize(line_pieces.size()-1);
      
      for(uint16_t this_time = 0; this_time < m_frequencies.size(); ++this_time) {
        if( from_string<double>(line_pieces[this_time+1]) != 0 ) {
          cGenotypeFrequency temp_time;
          temp_time.unique_node_id = num_nodes;
          temp_time.name = *(preformat_tree.end()-1);
          temp_time.frequency = from_string<double>(line_pieces[this_time+1]);
          m_frequencies[this_time].push_back(temp_time);
        }
      }
    }
    else {
      for(tree<cGenotype>::iterator this_node = m_tree.begin(); this_node != m_tree.end(); ++this_node) {
        if( *(preformat_tree.end()-2) == (*this_node).name ) {
          cSubpopulation new_child;
          cGenotype child;
          child.name = *(preformat_tree.end()-1);
          child.unique_node_id = num_nodes;
          child.fitness = 1;
          child.mut_num = num_nodes;
          new_child.AddToTree(m_tree, this_node, child);
          
          for(uint16_t this_time = 0; this_time < m_frequencies.size(); ++this_time) {
            if( from_string<double>(line_pieces[this_time+1]) != 0 ) {
              cGenotypeFrequency temp_time;
              temp_time.unique_node_id = num_nodes;
              temp_time.name = *(preformat_tree.end()-1);
              temp_time.frequency = from_string<double>(line_pieces[this_time+1]);
              m_frequencies[this_time].push_back(temp_time);
            }
          }
        }
      }
    }
    
    num_nodes++;
  }
  
  input_handle.close();
  
  /*for(tree<cGenotype>::iterator this_node = m_tree.begin(); this_node != m_tree.end(); ++this_node) {
    cout << this_node->unique_node_id << " " << this_node->name << endl;
  }*/
  
  //PrintFreqsQuick();
  //PrintTree();
}

void cPopulation::PrintFreqsQuick() {
  uint32_t time_keeper(0);
  for(vector< vector<cGenotypeFrequency> >::iterator this_time=m_frequencies.begin(); this_time!=m_frequencies.end(); ++this_time) {
    time_keeper++;
    cout << "Time: " << time_keeper << endl;
    for(vector<cGenotypeFrequency>::iterator this_genotype=this_time->begin(); this_genotype!=this_time->end(); ++this_genotype) {
      cout << this_genotype->unique_node_id << " " << this_genotype->name << " " << this_genotype->frequency << endl;
    }
    cout << endl << endl;
  }
}

std::vector<cGenotypeFrequency>::iterator cPopulation::Find_Node_in_Freq(std::vector<cGenotypeFrequency> &frequencies, tree<cGenotype>::sibling_iterator this_node) {
  
  cGenotypeFrequency return_node;
  return_node.unique_node_id = this_node->unique_node_id;
  
  for( std::vector<cGenotypeFrequency>::iterator a_node = frequencies.begin(); a_node < frequencies.end(); ++a_node ){
    if( a_node->unique_node_id == this_node->unique_node_id ) {
      return a_node;
    }
    else
      return_node.frequency = 0;
  }
  
  std::vector<cGenotypeFrequency> return_vector;
  return_vector.push_back(return_node);
  
  return return_vector.begin();
}

double cPopulation::Find_Node_in_Freq_By_NodeID(std::vector<cGenotypeFrequency> &frequencies,
                                                uint32_t this_node) {
  for( std::vector<cGenotypeFrequency>::iterator a_node = frequencies.begin(); a_node < frequencies.end(); ++a_node ) {
    if( a_node->unique_node_id == this_node )
      return a_node->frequency;
  }
  
  return 0;
}

cSubpopulation* cPopulation::Find_Node_in_Populations_By_NodeID(uint32_t this_node) {
  for( uint32_t it = 0; it < m_current_subpopulations.size(); ++it ) {
    if( m_current_subpopulations[it].GetNode_id() == this_node ) {
      return &m_current_subpopulations[it];
    }
  }
  assert(0==1);
}


// @JEB function returns the high frequency of the tallest child
//      (to tell the parent where to put its low frequency)
double cPopulation::AssignChildFreq(tree<cGenotype>::sibling_iterator this_node,
                                  double in_low,
                                  double in_high,
                                  std::vector<cFrequencySlice> * child_freqs,
                                  std::vector<cGenotypeFrequency> &frequencies, 
                                  int depth)  // current depth in tree, defaults to zero
{
  
  //kptree::print_tree_bracketed(*newtree);
  
  //The low for this mutation should be the low of the input interval
  double this_low = in_low;
  double this_high = this_low + Find_Node_in_Freq_By_NodeID(frequencies, this_node->unique_node_id);
  double size_depth1_children(0), half_size_parent_swath((this_high-this_low)/2);
  
  for (tree<cGenotype>::sibling_iterator it_node = m_tree.begin(this_node); it_node!=m_tree.end(this_node); ++it_node) {
    size_depth1_children += Find_Node_in_Freq_By_NodeID(frequencies, it_node->unique_node_id);
  }
  /*if(size_depth1_children != 0)
    std::cout << size_depth1_children << std::endl;*/
  
  double this_bottom_high(this_low + half_size_parent_swath - (size_depth1_children/2));
  double this_top_low(this_bottom_high);
  
  // The swath for this mutation may shrink on the bottom
  // due to its children taking a bite out of it.  
  double last_assigned_child_high = this_top_low;
  for (tree<cGenotype>::sibling_iterator it_node = m_tree.begin(this_node); it_node!=m_tree.end(this_node); ++it_node) {
    
    // is the frequency > 0? (It may have gone extinct, in which case it is a waste to keep going down the tree!)
    if( Find_Node_in_Freq_By_NodeID(frequencies, it_node->unique_node_id) > 0 ) {
      last_assigned_child_high = AssignChildFreq(it_node, this_top_low, this_high, child_freqs, frequencies, depth+1);
    }
    
    // Our new low is the last high assigned to a child
    this_top_low = last_assigned_child_high;
  }
  
  // At this point we know the top and bottom of this node...
  // This is where we would want to paint into a bitmap between them!!
  // Draw between this_low and this_high ... rounding to only paint whole numbers ...
  
  (*child_freqs).push_back(cFrequencySlice(this_node->unique_node_id, this_low, this_bottom_high));
  (*child_freqs).push_back(cFrequencySlice(this_node->unique_node_id, this_top_low, this_high));
  
  // print indented version
  /*for (int i=0; i<depth; i++) {
    std::cout << " ";
  }*/
  //std::cout << " ID:" << this_node->unique_node_id << " Freq:" << (*frequencies)[this_node->unique_node_id].frequency;
  //std::cout << " [" << this_low << "," << this_high << "]" << std::endl;
                              
  return this_high;  
}

void cPopulation::DrawMullerMatrix(std::string output_folder,
                                   std::vector< std::vector<int> > muller_matrix){
  
  //std::vector< tree<cGenotype>::iterator > where(m_tree.size());
  std::map<uint32_t, uint32_t> renumber;
  uint32_t renumber_value(0);
  
  //double threshold(.025);
  //std::vector<bool> relevant_columns(newtree->size(), true);
  
  //Build vector to store all iterators in tree for parental recall later
  //@JEB should we store this iterator in cGenotypeFrequency?
  //@AGM storing transient information like iterators in an object won't work
  
  /*for (tree<cGenotype>::iterator node_loc = m_tree.begin(); node_loc != m_tree.end(); node_loc++)
    where[(*node_loc).unique_node_id] = node_loc;*/
  
  std::string output_file;
  output_file.append(output_folder);
  output_file.append("/MullerMatrix.dat");
  std::ofstream output_handle(output_file.c_str());
  
  uint32_t time = 0;
  
  //step through simulation time
  for (std::vector< std::vector<cGenotypeFrequency> >::iterator this_time_freq = m_frequencies.begin(); this_time_freq < m_frequencies.end(); ++this_time_freq) {
    std::cout << time << std::endl;
    time++;
    
    std::vector<cFrequencySlice> child_freqs;
    
    tree<cGenotype>::sibling_iterator location;
    location = m_tree.begin();
    
    AssignChildFreq(location, 0, 1, &child_freqs, *this_time_freq);
    sort(child_freqs.begin(), child_freqs.end(), cSortByLow());
    
    /*cout << "Time " << time << " Child freqs size: " << child_freqs.size() << endl;
    for (uint32_t j=0; j<child_freqs.size(); j++) {
      cout << child_freqs[j].low << " " << child_freqs[j].high << endl;
    }
    PrintTree();*/
    
    uint32_t resolution(m_muller_rez), last_node_meeting_span;
    double pixel_step, min_step;
    min_step = (double) 1/resolution;
    
    //@agm Here I first iterate through the number of pixels
    for (uint32_t i=1; i<=resolution; i++) {
      
      //Determine the position of the current pixel_step
      pixel_step = (double) i/resolution;
      
      for (uint32_t j=0; j<child_freqs.size(); j++) {
        
        if( child_freqs[j].high >= pixel_step ) {
          
          //if( span > min_step ) 
          {
            
            //Add new significant mutations to a map for renumbering
            if( renumber.count(child_freqs[j].unique_node_id) == 0 ) {
              renumber.insert(std::make_pair(child_freqs[j].unique_node_id, renumber_value));
              renumber_value++;
              renumber_value %= resolution;
            }
            //Return the renumbered-number for the unique_node_id from the built map
            output_handle << std::left << std::setw(6) << renumber.find(child_freqs[j].unique_node_id)->second;
            last_node_meeting_span = renumber.find(child_freqs[j].unique_node_id)->second;
          }

          break;
        }
      }
    }
    output_handle << std::endl;
  }
}

void cPopulation::Resample() 
{
  //When it is time for a transfer, resample population
  assert(m_rng);
  
  SetTransfers(GetTransfers()+1);
 
  std::cout << ">> Transfer!" << std::endl;
  //Is there an exists() in C++?
  m_by_color[RED] = 0;
  m_by_color[WHITE] = 0;
	
  m_population_size = 0; // recalculate population size
  for (std::vector<cSubpopulation>::iterator it = m_current_subpopulations.begin(); it!=m_current_subpopulations.end(); ++it) {
    // Perform accurate binomial sampling only if below a certain population size
    if (it->GetNumber() < GetBinomialSamplingThreshold()) {
      if (g_verbose) std::cout << "binomial " << it->GetNumber() << std::endl;
      it->Transfer(GetTransferBinomialSamplingP(), m_rng);
    }
    // Otherwise, treat as deterministic and take expectation...
    else {
      it->SetNumber(it->GetNumber() * GetTransferBinomialSamplingP());
    }
         
    if (g_verbose) { 
      std::cout << it->GetMarker() << std::endl;
      std::cout << it->GetNumber() << std::endl;
      std::cout << it->GetFitness() << std::endl;
    }	
    if( it->GetMarker() == 'r' ) m_by_color[RED] += it->GetNumber();
    else if( it->GetMarker() == 'w' ) m_by_color[WHITE] += it->GetNumber();
    
    //@agm This section is to delete new mutations from the tree that do not get passed, 
    //     it also deletes subpopulations from the list that have zero population
    
    m_population_size += floor(it->GetNumber());
    
    if( it->GetNumber() == 0 ) {
      it = m_current_subpopulations.erase(it);
      --it;
    
      SetTotalSubpopulationsLost(GetTotalSubpopulationsLost()+1);
    }
  }
    
  if (g_verbose) std::cout << "Colors: " << m_by_color[RED] << " / " << m_by_color[WHITE] << std::endl;
  SetRatio( (double) m_by_color[RED]/m_by_color[WHITE] );

  // Checks for stopping early if in marker divergence mode
  if (g_ro_only) {
    
    // One color was lost -- bail
    if ( (m_by_color[RED] == 0) || (m_by_color[WHITE] == 0) ) {
      m_keep_transferring = false;
    }
    
    // We have the minimum number of transfers to be printed and have diverged sufficiently -- bail
    else {
      if  ( (GetRatio() > GetMaxDivergenceFactor()) || 
            (GetRatio() < 1/GetMaxDivergenceFactor()) )   
      {
        if (g_verbose) std::cout << "DIVERGENCE CONDITION MET" << std::endl;
        m_keep_transferring = false;
      }  
    }
  }
  
  //cout << "Pop size: " << m_population_size << endl;
}

//Currently this only works for picking one cell after resample
void cPopulation::Deterministic_Resample() {
  long int starting_size(0), random_number;
  //std::vector<cSubpopulation> temp_populations;
  
  assert(m_rng);
  
  SetTransfers(GetTransfers()+1);
  m_population_size = 0; // recalculate population size
  
  for (std::vector<cSubpopulation>::iterator it = m_current_subpopulations.begin(); it!=m_current_subpopulations.end(); ++it) {
    starting_size += it->GetNumber();
  }
  
  //std::cout << "Population size before dilution: " << starting_size << std::endl;
  
  //for (uint32_t num_cells = 0; num_cells < m_pop_size_after_dilution; num_cells++) {
    random_number = gsl_rng_uniform_int(m_rng, starting_size);
    
    for (std::vector<cSubpopulation>::iterator it = m_current_subpopulations.begin(); it!=m_current_subpopulations.end(); ++it) {
      if( random_number < it->GetNumber() ) {
        it->SetNumber(1);
        it++;
        while(it != m_current_subpopulations.end()) {
          m_current_subpopulations.erase(it);
        }
        break;
      }
      else {
        random_number -= it->GetNumber();
        m_current_subpopulations.erase(it);
        --it;
      }
    }
  //}
  
  /*for (std::vector<cSubpopulation>::iterator it = m_current_subpopulations.begin(); it!=m_current_subpopulations.end(); ++it) {
    std::cout << it->GetNumber() << " " << it->GetNode_id() << " " << it->GetFitness() << std::endl;
  }*/
  
  m_population_size = 1;
}



//@agm For some reason (that I don't fully appreciate), the way we were deleting subpopulations 
//     (based on the GetNumber()==0) Therefore, I am now deleting subpopulations Only when they 
//     arose between transfers And do not get passed to the next round just like those that are
//     removed from the tree.  This is a much more conservative pruning so I should try to figure
//     out why the other does not work.
     
void cPopulation::CullPopulations() {
  
  std::cout << std::endl;
  
  std::vector<cSubpopulation>::iterator this_mutation(m_mutations_since_last_transfer.begin());
  bool test(false);
  while( this_mutation != m_mutations_since_last_transfer.end() ) {
    for (std::vector<cSubpopulation>::iterator it = m_current_subpopulations.begin(); it!=m_current_subpopulations.end(); ++it) {
      if( this_mutation->GetNode_id() == it->GetNode_id() ) {
        this_mutation = m_mutations_since_last_transfer.erase(this_mutation);
        test = true;
        //std::cout << it->GetNode_id() << std::endl;
        break;
      }
    }
    
    if( test == false )
      this_mutation++;
    
    test = false;
    //std::cout << "I'm Here: 1 " << this_mutation->GetNode_id() << std::endl;
  }
  
  //std::cout << "I'm Here: 2 " << this_mutation->GetNode_id() << std::endl;
  
  this_mutation = m_mutations_since_last_transfer.begin();
  
  for (std::vector<cSubpopulation>::reverse_iterator rit = m_mutations_since_last_transfer.rbegin(); rit!=m_mutations_since_last_transfer.rend(); ++rit) {
    
    //std::cout << "I'm Here: 2 " << rit->GetNode_id() << std::endl;

    m_tree.erase(rit->GetGenotypeIter());
  }
  m_mutations_since_last_transfer.clear();
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

void cPopulation::DisplayParameters()
{
  if (g_verbose) {
    std::cout << "u = " << GetMutationRatePerDivision() << std::endl;
    std::cout << "s = " << GetAverageMutationS() << std::endl;
    std::cout << "N = " << GetPopSizeAfterDilution() << std::endl;
    std::cout << "dil = " << GetDilutionFactor() << std::endl;
  }
}

void cPopulation::CalculateDivisions()
{
  
  if (g_verbose) {
    std::cout << "Completed divisions: " << GetCompletedDivisions() << std::endl;
    for(vector<cSubpopulation>::iterator this_time = m_current_subpopulations.begin(); this_time != m_current_subpopulations.end(); this_time++) {
      cout << "Genotype: " << this_time->GetNode_id() << " Frequency: " << this_time->GetNumber() << endl;
    }
  }
  
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
  
  // Note: we underestimate by a few divisions so that we can step forward by single division increments
  // as we get close to the one where the mutation happened (or right before a transfer).      
  
  if (desired_divisions < 1) {
     desired_divisions = 1;
  }
  
  if (g_verbose) {
     std::cout << "Total pop size: " << GetPopulationSize() <<std::endl;
     std::cout << "Desired divisions: " << desired_divisions <<std::endl;
  }
      
  // How much time would we like to pass to achieve the desired number of divisions?
  // (Assuming the entire population has the maximum fitness, makes us underestimate by a few)
  // Can't change this log to ReturnLog with log table or nothing works
  
  double update_time(0);

  update_time = log((desired_divisions+(double)GetPopulationSize()) / (double)GetPopulationSize()) / (GetMaxW());

  // At a minumum, we want to make sure that one cell division took place

  // What is the minimum time required to get a single division?
  double time_to_next_whole_cell = TimeToNextWholeCell();

  if (time_to_next_whole_cell > update_time) {
     if (g_verbose) std::cout << "Time to next whole cell greater than update time: " << 
			 time_to_next_whole_cell << " < " << update_time <<std::endl;    
     update_time = time_to_next_whole_cell;
  }

  if (g_verbose) {
     std::cout << "Update time: " << update_time <<std::endl;    
  }
            
  //Now update all lineages by the time that actually passed
 
  //FIX THIS!!!!
  //Population size shouldn't have to be recaculated here
  //I think this is fixed... the problem was in the AddSubpopulation() method
  //m_population_size = CalculatePopulationSize();
  
  double previous_population_size(m_population_size);
  
  UpdateSubpopulations(update_time);
  
  SetCompletedDivisions(GetPopulationSize() - previous_population_size);
              
  if (g_verbose) {
    std::cout << "Completed divisions: " << GetCompletedDivisions() << std::endl;
    for(vector<cSubpopulation>::iterator this_time = m_current_subpopulations.begin(); this_time != m_current_subpopulations.end(); this_time++) {
      cout << "Genotype: " << this_time->GetNode_id() << " Frequency: " << this_time->GetNumber() << endl;
    }
  }
  
  //std::cout << "Divisions until mutation: " << GetDivisionsUntilMutation() << " Completed divisions: " << GetCompletedDivisions() << " Population Size: " << GetPopulationSize() << " Previous Population Size: " << previous_population_size << std::endl;
  SetDivisionsUntilMutation(GetDivisionsUntilMutation() - GetCompletedDivisions());

}

void cPopulation::CalculateDivisionsNew() {
  //Calculate time to next division by having the initial time be 1
  //The time remains 1 until a cell mutates
  //at that point each population is iterated through and the time variable
  //in each subpopulation is set by time = old_time - (elapsed_time * Fitness)
  //Thus, every step is going to entail cycling through to find the smallest remaining time
  //doubling that pop, updating the time in pops, then finding the smallest pop again.
  
  //This will happen until the mutation time is breached at which point the next population to 
  //divide will get a mutation.
  
  double gets_to_double(numeric_limits<double>::max());
  double number_of_resources_here(0), mean_weighted_fitness(0);
  double number_of_resources_to_transfer(m_pop_size_before_dilution - m_population_size);
  
  for(uint32_t time = 0; time < m_current_subpopulations.size(); time++) {
    cSubpopulation&  this_time(m_current_subpopulations[time]);
    this_time.SetDivided(false);
    
    double a_time( ( 2.0 - this_time.GetTimer() ) / this_time.GetFitness() );
    this_time.SetTimer( 2.0 - a_time );
    
    if( a_time < gets_to_double ) {
      gets_to_double = a_time;
    }
  }
  
  m_elapsed_time = gets_to_double;
  
  for(uint32_t time = 0; time < m_current_subpopulations.size(); time++) {
    cSubpopulation&  this_time(m_current_subpopulations[time]);
    
    double a_time( ( 2.0 - this_time.GetTimer() ) / this_time.GetFitness() );
    
    mean_weighted_fitness += this_time.GetNumber()*this_time.GetFitness();
    
    if( a_time == m_elapsed_time ) {
      number_of_resources_here += this_time.GetNumber();
      this_time.SetDivided(true);
      //cout << "Node id of next doubling: " << this_time.GetNode_id() << " " << m_elapsed_time << endl;
    }
  }
  
  mean_weighted_fitness /= m_population_size;
  
  //cout << " Doubled! " << endl;
  
  if( number_of_resources_here < number_of_resources_to_transfer ) {
  
    for(uint32_t time = 0; time < m_current_subpopulations.size(); time++) {
      cSubpopulation&  this_time(m_current_subpopulations[time]);
      //cout << this_time.GetNode_id() << " " << this_time.GetNumber() << " " << this_time.GetFitness() << endl;
      
      this_time.SetTimer( this_time.GetTimer() + m_elapsed_time );
      assert(this_time.GetTimer() <= 2);
      
      if( this_time.GetDivided()  ) {
        
        SetCompletedDivisions( GetCompletedDivisions() + this_time.GetNumber() );
        SetDivisionsUntilMutation( GetDivisionsUntilMutation() - this_time.GetNumber() );
        m_population_size += this_time.GetNumber();
        
        this_time.SetNumber( this_time.GetNumber() * this_time.GetTimer() );
        
        this_time.SetTimer(1.0);
      }
    }
  }
  
  //@agm To fix the stair-steppy thing going on in the frequency
  //     In the frequency vector, simply multiply the Timer and
  //     number of cells to get the number of cell equivalents
  
  else if( number_of_resources_to_transfer < number_of_resources_here ) {
    double population_size_now = m_population_size;
    
    for(uint32_t time = 0; time < m_current_subpopulations.size(); time++) {
      cSubpopulation&  this_time(m_current_subpopulations[time]);
      //cout << this_time.GetNode_id() << " " << this_time.GetNumber() << endl;
      
      double number_cells_added( this_time.GetNumber() * number_of_resources_to_transfer / population_size_now );
      
      SetCompletedDivisions( GetCompletedDivisions() + number_cells_added );
      SetDivisionsUntilMutation( GetDivisionsUntilMutation() - number_cells_added );
      m_population_size += number_cells_added;
      
      this_time.SetNumber( this_time.GetNumber() + number_cells_added );
      this_time.SetTimer( this_time.GetTimer() - number_cells_added / this_time.GetNumber() );
      assert(this_time.GetTimer() <= 2);

      this_time.SetDivided(false);
    }
    
  }
  
  //cout << "Size of population: " << m_population_size << " " << CalculatePopulationSize() << endl;
  
}

/*@agm The functions below should build a new tree using the tree.h header */

void cPopulation::SeedSubpopulationForRedWhite() 
{	
	cGenotype r, w;
	tree<cGenotype>::iterator top, red_side, white_side;
	double starting_fitness = GetInitialFitness();
	
	//initialize object of cSubpopulation type
	cSubpopulation red, white;
	
	/*@agm The unique code and fitness is set. */
  
  /*@jeb note that we increment genotype count here.
   It should really happen in AddSubpopulation?
   */
  
	r.fitness = starting_fitness;
	r.unique_node_id = m_genotype_count++;
  w.fitness = starting_fitness;
	w.unique_node_id = m_genotype_count++;
	
	//Start building the tree
	red_side = m_tree.insert(m_tree.begin(), r);
	white_side = m_tree.insert(m_tree.begin(), w);
	
	red.SetNumber(GetInitialPopulationSize()/2);
	red.SetGenotype(red_side);
	red.SetMarker('r');
	
	white.SetNumber(GetInitialPopulationSize()/2);
	white.SetGenotype(white_side);
	white.SetMarker('w');
	
	AddSubpopulation(red);
	AddSubpopulation(white);	
}

//@agm This function seeds the population with only one colony to avoid the red/white problem
//     when possible. For this to work properly, I need to get rid of victory conditions.

//update: Victory conditions are now a command line option.  If not set to false the simulation
//        will run to the max number of iterations.

void cPopulation::SeedPopulationWithOneColony() {
  tree<cGenotype>::iterator start_position;
  double starting_fitness(GetInitialFitness());
  
  cGenotype neutral;
  neutral.fitness = starting_fitness;
  neutral.unique_node_id = m_genotype_count++;
  neutral.mut_num = 0;
  start_position = m_tree.insert(m_tree.begin(), neutral);
  
  cSubpopulation begin_here;
  begin_here.SetNumber(GetInitialPopulationSize());
  begin_here.SetGenotype(start_position);
  begin_here.SetMarker('n');
  AddSubpopulation(begin_here);
  
  if (g_verbose) {
    std::cout << "Completed divisions: " << GetCompletedDivisions() << std::endl;
    for(vector<cSubpopulation>::iterator this_time = m_current_subpopulations.begin(); this_time != m_current_subpopulations.end(); this_time++) {
      cout << "Genotype1: " << this_time->GetNode_id() << " Frequency1: " << this_time->GetNumber() << endl;
    }
  }
  
}

void cPopulation::AddSubpopulation(cSubpopulation& subpop) 
{
  m_population_size += subpop.GetNumber(); // we have just changed the population size
	m_current_subpopulations.push_back(subpop);
  if(subpop.GetNode_id() > 0 && !g_ro_only)
    m_mutations_since_last_transfer.push_back(subpop); // Keep track of all mutations since last transfer

	//std::cout << subpop.GetNode_id() << " " << subpop.GetFitness() << std::endl;
}

//Generates new mutant and adds it to the tree
void cPopulation::Mutate() 
{	
  // we better have a random number generator
  assert(m_rng);
	
	if (g_verbose) std::cout << "* Mutating!" << std::endl;
  if (g_verbose) std::cout << "Total population: " << CalculatePopulationSize() << std::endl;
	
	//Mutation happened in the one that just divided
	//Break ties randomly here.
	cSubpopulation& ancestor = m_current_subpopulations[m_divided_lineages[rand() % m_divided_lineages.size()]];          
	cSubpopulation new_subpop;
  
  //std::cout << "Divided has number: " << ancestor.GetNumber() << std::endl;
  // There must be at least two cells for a mutation to have occurred...
  assert(ancestor.GetNumber() >= 2);
	
  //This is necessary so the first few mutation have a much larger fitness advantage
  
  if( m_first_mutational_vals.size() > ancestor.GetMutNum() ) {
    new_subpop.CreateDescendant(m_rng, 
                                ancestor, 
                                m_first_mutational_vals[ancestor.GetMutNum()], 
                                GetBeneficialMutationDistribution(),
                                m_tree,
                                m_genotype_count++);
  }
  else {
    new_subpop.CreateDescendant(m_rng, 
                                ancestor, 
                                m_first_mutational_vals[m_first_mutational_vals.size()-1], 
                                GetBeneficialMutationDistribution(),
                                m_tree,
                                m_genotype_count++);
  }
  
	if (g_verbose) std::cout << "  Color: " << new_subpop.GetMarker() << std::endl;
	if (g_verbose) std::cout << "  New Fitness: " << new_subpop.GetFitness() << std::endl;
  
	AddSubpopulation(new_subpop);
  
  //@agm Since an existent cell is picked to mutate,
  //     rather than doubling a cell and picking its
  //     progeny to mutate... we should not add a new
  //     cell to the population, but AddSubpopulation does.
  m_population_size-=1;
	
	if(new_subpop.GetFitness() > GetMaxW()) 
	{
		SetMaxW(new_subpop.GetFitness());
	}
	
  //ancestor.SetNumber(ancestor.GetNumber()-1);
  
  if (g_verbose) std::cout << "  New Number: " << new_subpop.GetNumber() << std::endl;
  if (g_verbose) std::cout << "  Old Number: " << ancestor.GetNumber() << std::endl;
  
  m_total_mutations++;
}

void cPopulation::MutateNew() {
  
  //@agm The code gets really finicky here
  //     if you try to change the ancestor to a pointer
  //     or a reference you will get memory leaks
  //     I have not figured out why, but it does not leak if you pass the full object
	
	if (g_verbose) std::cout << "* Mutating!" << std::endl;
  if (g_verbose) std::cout << "Total population: " << CalculatePopulationSize() << std::endl;
  
  vector< uint32_t > populations_to_mutate;
  
  //find the subpopulation or subpopulations that just divided and mutate them
  for(vector< cSubpopulation >::iterator pos = m_current_subpopulations.begin(); pos != m_current_subpopulations.end(); pos++) {
    //cSubpopulation& ancestor = m_current_subpopulations[pos];
    //cout << pos->GetNode_id() << " " << pos->GetNumber() << endl;
    
    if( pos->GetDivided() ) {
      populations_to_mutate.push_back(pos->GetNode_id());
      //cout << "First loop: " << pos->GetNode_id() << " " << pos->GetDivided() << " " << pos->GetNumber() << endl;
    }
  }
  
  //find the subpopulation or subpopulations that just divided and mutate them
  for(uint32_t pos = 0; pos < populations_to_mutate.size(); pos++) {

    cSubpopulation ancestor;
    ancestor = (* Find_Node_in_Populations_By_NodeID( populations_to_mutate[pos] ) );
    double number_of_cells( ancestor.GetNumber() );
    //cout << ancestor.GetNode_id() << " " << ancestor.GetNumber() << endl;

    cSubpopulation new_subpop;
      
    //std::cout << "Divided has number: " << ancestor.GetNumber() << std::endl;
    // There must be at least two cells for a mutation to have occurred...
    //cout << "Second loop: " << ancestor->GetNode_id() << " " << ancestor->GetDivided() << " " << ancestor->GetNumber() << endl;
    assert(number_of_cells >= 2);
    
    //This is necessary so the first few mutation have a much larger fitness advantage
      
    if( m_first_mutational_vals.size() > ancestor.GetMutNum() ) {
      new_subpop.CreateDescendant(m_rng, 
                                  ancestor, 
                                  m_first_mutational_vals[ancestor.GetMutNum()], 
                                  GetBeneficialMutationDistribution(),
                                  m_tree,
                                  m_genotype_count++);
    }
    else {
      new_subpop.CreateDescendant(m_rng, 
                                  ancestor, 
                                  m_first_mutational_vals[m_first_mutational_vals.size()-1], 
                                  GetBeneficialMutationDistribution(),
                                  m_tree,
                                  m_genotype_count++);
    }
    
    if (g_verbose) std::cout << "  Color: " << new_subpop.GetMarker() << std::endl;
    if (g_verbose) std::cout << "  New Fitness: " << new_subpop.GetFitness() << std::endl;
    
    AddSubpopulation(new_subpop);
    
    //@agm Since an existent cell is picked to mutate,
    //     rather than doubling a cell and picking its
    //     progeny to mutate... we should not add a new
    //     cell to the population, but AddSubpopulation does.
    m_population_size-=1;
    
    if(new_subpop.GetFitness() > GetMaxW()) SetMaxW(new_subpop.GetFitness());
    
    number_of_cells--;
    
    ancestor.SetNumber(number_of_cells);
    
    if (g_verbose) std::cout << "  New Number: " << new_subpop.GetNumber() << std::endl;
    if (g_verbose) std::cout << "  Old Number: " << ancestor.GetNumber() << std::endl;
    
    m_total_mutations++;
  }
  
}

//Utilities Section

std::vector<uint16_t> cPopulation::CurrentUniqueGenotypes() {
  std::vector<uint16_t> number_of_unique_genotypes;
  
  for (uint32_t time = 0; time<m_all_subpopulations_at_all_times.size(); time++) {
    uint16_t current_number(0);
    for (uint32_t this_one = 0; this_one < sizeof(m_all_subpopulations_at_all_times[time]); this_one++) {
      double& this_pop = m_all_subpopulations_at_all_times[time][this_one];
      if( (double) this_pop / m_total_cells[time] > .1 ) current_number++;
    }
    number_of_unique_genotypes.push_back(current_number);
  }
  return number_of_unique_genotypes;
}

void cPopulation::PrintUniqueGenotypes(const std::string& output_folder,
                                       std::vector< std::vector<uint16_t> > * number_of_unique_genotypes) {
  //Will print out only red and white
	std::ofstream output_handle;
  std::string output_file;
  
  output_file = output_folder + "/Number_Unique_Genotypes.dat";
  
	output_handle.open(output_file.c_str(),std::ios_base::app);
  
  uint16_t largest_replicate(0);
  for (uint16_t replicate=0; replicate<(*number_of_unique_genotypes).size(); replicate++) {
    if( (*number_of_unique_genotypes)[replicate].size() > (*number_of_unique_genotypes)[largest_replicate].size() ) largest_replicate = replicate;
  }
  
  output_handle << "transfer";
  for (uint16_t replicate = 0; replicate<(*number_of_unique_genotypes).size(); replicate++) {
    output_handle << "\t" << replicate ;
  }
  
  output_handle << std::endl;
  
  for (uint16_t time = 0; time<(*number_of_unique_genotypes)[largest_replicate].size(); time++) {
    output_handle << time;
    for (uint16_t replicate = 0; replicate<(*number_of_unique_genotypes).size(); replicate++) {
      if( time >= (*number_of_unique_genotypes)[replicate].size() )
        output_handle << "\t";
      else
        output_handle << "\t" << (*number_of_unique_genotypes)[replicate][time];
    }
    output_handle << std::endl;
	}
  
}

//@agm I comandeered this function to print stuff out in the manner I see fit
//@agm I basically do a non-human readable raw dump because it's easier for R to deal with
//     If you want human readable use the PrintToScreen function

/***** Important for post-hoc use with R ******
 The output file will have a header and it
 will have the largest number of columns in file.
 Thus, in R, you can simply set header to 
 true, and fill NA spaces as the read.table
 function allows. */

void cPopulation::PrintOut(const std::string& output_folder, uint32_t on_run)
{  
  
  std::vector<uint32_t> all_relevant_nodes( MutationAboveThreshold(1.0) );
  
  //Print everything out
	std::ofstream output_handle;
  std::string output_file;
  
  output_file.append(output_folder);
  output_file.append("/Genotype_Frequencies_");
  output_file.append(to_string(on_run));
  output_file.append(".dat");
	output_handle.open(output_file.c_str(),std::ios_base::app);
  
  for (uint32_t a_node=0; a_node<all_relevant_nodes.size(); a_node++) {
    output_handle << "Genotype_" << all_relevant_nodes[a_node] << " ";
  }
  
  output_handle << endl;
  
  for (uint32_t i = 0; i < m_frequencies.size(); ++i) {
    vector< cGenotypeFrequency >& this_time( m_frequencies[i] );
    for (uint32_t a_node=0; a_node<all_relevant_nodes.size(); a_node++) {
      //cout << a_node << endl;
      double frequency( Find_Node_in_Freq_By_NodeID(this_time, all_relevant_nodes[a_node]) );
      output_handle << frequency << " ";
    }
    output_handle << endl;
  }
  
  output_handle.close();

}

void cPopulation::PrintExpectationValue(const std::string& output_folder) {
  //Print everything out
	std::ofstream output_handle;
  std::string output_file;
  
  output_file.append(output_folder);
  output_file.append("/ExpectationVal.dat");
	output_handle.open(output_file.c_str(),std::ios_base::app);
  
  double calculated_mutant_population_freq(0);

  for ( std::vector<cGenotypeFrequency>::iterator it = m_frequencies[m_frequencies.size()-1].begin(); it!=m_frequencies[m_frequencies.size()-1].end(); ++it) {
    if( it->unique_node_id != 0 ) {
      calculated_mutant_population_freq += it->frequency;
    }
  }
    
  cout << calculated_mutant_population_freq << endl;
  
}

void cPopulation::PrintOut_RedWhiteOnly(const std::string& output_folder, 
                                        std::vector< std::vector<double> > * red_white_ratios,
                                        uint16_t transfer_interval_to_print) 
{  
	//Will print out only red and white
	std::ofstream output_handle;
  std::string output_file;
  
  output_file = output_folder + "_Gfreqs_RO.dat";
  
	output_handle.open(output_file.c_str(),std::ios_base::app);
  
  uint16_t largest_replicate(0);
  //find largest replicate
  for (uint16_t replicate = 0; replicate<(*red_white_ratios).size(); replicate++) {
    if( (*red_white_ratios)[replicate].size() > largest_replicate ) largest_replicate = replicate;
  }
  
  output_handle << "transfer";
  for (uint16_t replicate = 0; replicate<(*red_white_ratios).size(); replicate++) {
    output_handle << "\t" << replicate ;
  }

  output_handle << std::endl;
  
  for (uint16_t time = 0; time<(*red_white_ratios)[largest_replicate].size(); time++) {
    output_handle << time*transfer_interval_to_print;
    for (uint16_t replicate = 0; replicate<(*red_white_ratios).size(); replicate++) {
      if( time >= (*red_white_ratios)[replicate].size() )
        output_handle << "\t";
      else
        output_handle << "\t" << log((*red_white_ratios)[replicate][time]);
    }
    output_handle << std::endl;
	}
}

//@agm As the function name implies, this prints the frequencies above some threshold to screen at whatever
//     time it is called and passed the frequencies vector

void cPopulation::PrintFrequenciesToScreen(std::string output_folder) {
  int count(0);
  std::ofstream output_handle;
  std::string output_file;
  
  //@agm this should output exactly the same thing to file as it does to screen
  //     the virtue of this format is that it is human readable
  output_file.append(output_folder);
  output_file.append("/");
  output_file.append("HumanReadable_GenotypeFrequencies.dat");
	output_handle.open(output_file.c_str(),std::ios_base::app);
  
	for (uint32_t i = 0; i<m_frequencies.size(); i++) {
		double total_freqs = 0;
    for ( std::vector<cGenotypeFrequency>::iterator it = m_frequencies[i].begin(); it!=m_frequencies[i].end(); ++it) {
    //@agm set up a minimum frequency to report the print out the number so it isn't overwhelming.
      if ((*it).frequency > 0.01) {
        std::cout << "Frequency of mutation # " << std::right << std::setw(6) << (*it).unique_node_id << " at time ";
        output_handle << "Frequency of mutation # " << std::right << std::setw(6) << (*it).unique_node_id << " at time ";
        std::cout << std::right << std::setw(4) << i << " is: " << std::left << std::setw(10) << (*it).frequency << std::endl;
        output_handle << std::right << std::setw(4) << i << " is: " << std::left << std::setw(10) << (*it).frequency << std::endl;
      }
      total_freqs += (*it).frequency;
      count++;
    }
		std::cout << "Round # " << i << " sum of frequencies is: " << total_freqs << std::endl << std::endl;
    output_handle << "Round # " << i << " sum of frequencies is: " << total_freqs << std::endl << std::endl;
	}
} 

void cPopulation::PrintFrequenciesToScreen_RedWhiteOnly(std::string output_folder) {;}

//@agm This function determines the maximum difference in genotype frequency between a mutation
//     and its predecessor if they got above the threshold passed to the threshold function below

double cPopulation::CalculateSimilarity(std::string output_folder) {
  
  uint32_t youngest_sweep(MutationAboveThreshold_2(.98));
  
  std::vector<uint32_t> all_sweep_ids;
  
  for (tree<cGenotype>::iterator it = m_tree.begin(); it!=m_tree.end(); it++) {
    if( it->unique_node_id == youngest_sweep ) {
      while (it != NULL) {
        all_sweep_ids.push_back(it->unique_node_id);
        it = m_tree.parent(it);
      }
      break;
    }
  }
  
  std::sort(all_sweep_ids.begin(), all_sweep_ids.end());
  
  std::ofstream output_handle;
  std::string output_file;
  
  std::vector<float> max_diff(all_sweep_ids.size(), 0);
  float current_diff;
  double num_below_threshold(0);
  
  uint32_t i=0;
  
  for (uint16_t a_node=0; a_node<all_sweep_ids.size(); a_node++) {
    uint32_t time = 0;
    for (std::vector< std::vector<cGenotypeFrequency> >::iterator this_time = m_frequencies.begin(); this_time < m_frequencies.end(); ++this_time) {
      
      //This conditional is to coarse-grain the sampling to every 75 transfers (~500 generations)
      //and to make sure there is a value in the vector at that position
      //another way to think about it is this conditional is checking to make sure that
      //a particular mutation has arisen in time before comparing it to something.
      
      if ( time%m_coarse_graining == 0 ) {
        current_diff = fabs(Find_Node_in_Freq_By_NodeID(*this_time, all_sweep_ids[a_node]) - Find_Node_in_Freq_By_NodeID(*this_time, all_sweep_ids[a_node+1]));
        if( max_diff[i] < current_diff ) max_diff[i] = current_diff;
      }
      time++;
    }
    i++;
  }
  
  output_file = output_folder + "/" + "SweepClumpiness.dat";
  output_handle.open(output_file.c_str(), std::ios_base::app);
  
  uint16_t num_simultaneous(1);

  for (int i = 0; i<max_diff.size()-1; i++) {
    if( max_diff[i] <= .1 ) num_simultaneous++; 
    else {
      output_handle << num_simultaneous << std::endl;
      num_simultaneous = 1;
    }
  }
  
  output_handle.close();
  
  output_file = output_folder + "/" + "SignificantParallelMutations.dat";
	output_handle.open(output_file.c_str(), std::ios_base::app);
  
  for (int i = 0; i<max_diff.size()-1; i++) {
    std::cout << std::endl << i << " " << max_diff[i] << std::endl;
    output_handle << max_diff[i] << std::endl;
    if( max_diff[i] <= .1 ) num_below_threshold++;
  }
  
  output_handle.close();
  return num_below_threshold;
}

void cPopulation::TimeToSweep(std::string output_folder) {
  
  uint32_t youngest_sweep(MutationAboveThreshold_2(.98));
  
  std::vector<uint32_t> all_sweep_ids;
  
  for (tree<cGenotype>::iterator it = m_tree.begin(); it!=m_tree.end(); it++) {
    if( it->unique_node_id == youngest_sweep ) {
      while (it != NULL) {
        all_sweep_ids.push_back(it->unique_node_id);
        it = m_tree.parent(it);
      }
      break;
    }
  }
  
  std::sort(all_sweep_ids.begin(), all_sweep_ids.end());
  
  std::ofstream output_handle;
  std::string output_file;
  
  std::vector<uint16_t> time_to_sweep(all_sweep_ids.size(), 0);
  
  uint32_t i=0;
  
  for (uint16_t a_node=0; a_node<all_sweep_ids.size(); a_node++) {
    uint32_t time = 0;
    for (std::vector< std::vector<cGenotypeFrequency> >::iterator this_time = m_frequencies.begin(); this_time < m_frequencies.end(); ++this_time) {
      
      //This conditional is to coarse-grain the sampling to every 75 transfers (~500 generations)
      //and to make sure there is a value in the vector at that position
      //another way to think about it is this conditional is checking to make sure that
      //a particular mutation has arisen in time before comparing it to something.
      
      double frequency( Find_Node_in_Freq_By_NodeID(*this_time, all_sweep_ids[a_node]) );
      
      if ( time%m_coarse_graining == 0 && frequency != 1 && frequency > .01) {
        time_to_sweep[a_node]+=3;
      }
      time++;
    }
    i++;
  }
  
  //@agm this should output exactly the same thing to file as it does to screen
  //     the virtue of this format is that it is human readable
  output_file.append(output_folder);
  output_file.append("/");
  output_file.append("TimeToSweep.dat");
	output_handle.open(output_file.c_str(),std::ios_base::app);
  
  for (std::vector<uint16_t>::iterator it_node = time_to_sweep.begin(); it_node != time_to_sweep.end(); it_node++) {
    if( (*it_node) != 0 ) {
      std::cout << (*it_node) << std::endl;
      output_handle << (*it_node) << std::endl;
    }
  }
}

uint32_t cPopulation::MutationAboveThreshold_2(float threshold) {
  
  for( std::vector< std::vector<cGenotypeFrequency> >::reverse_iterator time = m_frequencies.rbegin(); time < m_frequencies.rend(); ++time ) {
    for( std::vector<cGenotypeFrequency>::reverse_iterator node = time->rbegin(); node < time->rend(); ++node) {
      if( node->frequency >= threshold )
        return node->unique_node_id;
    }
  }
  
  std::cout << "Nothing met threshold." << std::endl;
  
  return 0;
}

//@agm This function takes the frequency pointer and returns a boolean vector
//     The boolean vector contains a true in the mutations that got above the passed threshold

std::vector<uint32_t> cPopulation::MutationAboveThreshold(float threshold) {
  
  vector< uint32_t > Fixed;
  
  for (vector< vector< cGenotypeFrequency > >::iterator it_time = m_frequencies.begin(); it_time!=m_frequencies.end(); ++it_time) {
    for( vector<cGenotypeFrequency>::iterator it_node = it_time->begin(); it_node != it_time->end(); ++it_node) {
      if( it_node->frequency >= threshold && !binary_search(Fixed.begin(), Fixed.end(), it_node->unique_node_id)) {
        Fixed.push_back(it_node->unique_node_id);
        sort(Fixed.begin(), Fixed.end());
        cout << it_node->unique_node_id << " ";
      }
    }
  }
  cout << endl;
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

void cPopulation::PrintTree() {
  kptree::print_tree_bracketed(m_tree);
  std::cout << std::endl;
}