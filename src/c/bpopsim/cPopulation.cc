#include "common.h"
#include "cPopulation.h"

using namespace bpopsim;
using namespace std;


cPopulation::cPopulation(AnyOption& options, gsl_rng* in_rng, uint32_t in_replicate)
: simulation_parameters(options, in_rng)
, output_parameters(options)
, marker_divergence(options)
, replicate(in_replicate)
, unique_genotype_count(0)
, existing_genotype_count(0)
, maximum_subpopulation_fitness(simulation_parameters.initial_fitness)
, current_population_size(0)
, num_completed_transfers(0)
, run_end_condition_met(false)
, num_divisions_until_next_mutation(0)
, num_completed_divisions(0)
, rng(in_rng)
, log_2(log(2))
{
  debug = options.count("debug");
}

void cPopulation::DisplayParameters()
{
  simulation_parameters.Print();
}

void cPopulation::RunSimulation(bool new_simulation_model)
{
  //Initialize the population
  if (marker_divergence.marker_states) {
    SeedPopulationWithMarkedGenotypes(marker_divergence.marker_states); 
  } else {
    SeedPopulationWithOneGenotype();
  }
  
  //Get an initial time points
  RecordStatisticsAtTransfer();
  
  // Print the initial tree
  //if (g_verbose) population.PrintTree();

  while( num_completed_transfers < simulation_parameters.maximum_number_of_transfers ) {
    
    num_divisions_until_next_mutation = CalculateDivisionsUntilNextBeneficialMutation();
    
    while( num_divisions_until_next_mutation > 0 && (num_completed_transfers < simulation_parameters.maximum_number_of_transfers)) {
      
      if (new_simulation_model)
        ProcessCellDivisionTimeStepNew();
      else
        ProcessCellDivisionTimeStep();
      
      // Debug check on mutate function
      int64_t cells_before_mutate;
      if (debug) cells_before_mutate = CalculatePopulationSize();
      
      if( num_divisions_until_next_mutation <= 0 ) {
        if (new_simulation_model)
          MutateNew();
        else
          Mutate();
      }
      
      if (debug) {
        int64_t cells_after_mutate = CalculatePopulationSize();
        assert(cells_before_mutate == cells_after_mutate);
      }
        
      if( current_population_size >= simulation_parameters.final_population_size_at_transfer ) {
        
        cerr << "    Replicate: " << setw(3) << replicate << "  Transfer: " << setw(4) << num_completed_transfers;
        cerr << "  Current Genotypes: " << setw(6) << existing_genotype_count; 
        cerr << endl;
        
        // Record statistics
        RecordStatisticsAtTransfer();

        // @JEB: Would need to re-implement mutation accumulation version
        /*
          if( print_single_fit )
            TransferResampleExactlyOne();
          else {
        */ 
        TransferResampleDilution();
        CullSubpopulationsThatDidNotEstablish();
        //  }
          
        /* @JEB: Would need to re-implement...
         
          if( print_single_fit && !use_mute_num) {
            population.PrintSingleFitness(options["output-folder"]);
            //cout << "Population size: " << population.GetPopulationSize() << endl;
          }
         */
      }
    }
  }
}

void cPopulation::DisplaySimulationSummary()
{
  cerr << "Total mutations: " << replicate_statistics.total_mutations << endl;
  cerr << "Total subpopulations lost: " << replicate_statistics.total_subpopulations_lost << endl;
  cerr << "Transfers: " << num_completed_transfers << endl;
  cerr << "Maximum Fitness: " << maximum_subpopulation_fitness   << endl;
}

//@agm I comandeered this function to print stuff out in the manner I see fit

/***** Important for post-hoc use with R ******
 The output file will have a header and it
 will have the largest number of columns in file.
 Thus, in R, you can simply set header to 
 true, and fill NA spaces as the read.table
 function allows. */

void cPopulation::OutputCladeFrequencies(double frequency_threshold)
{  
  //vector<uint32_t> all_sweep_ids = GenotypesFromAncestorToFinalDominant(frequency_threshold);
  
  vector<uint32_t> all_sweep_ids = CladesAboveThreshold(frequency_threshold);

  
  //Print each sweeping genotype and its frequency per time
	ofstream output_file;
  string output_file_name = output_parameters.output_directory_name + "/clade_frequencies_" 
    + to_string(frequency_threshold) + "_" + to_string(replicate) + ".dat";
	output_file.open(output_file_name.c_str(),ios::out);
    
  cerr << "Output: " << output_file_name << endl;
  
  output_file << "transfer";
  for (uint32_t a_node=0; a_node<all_sweep_ids.size(); a_node++) {
    output_file << "\t" << "clade_" << all_sweep_ids[a_node];
  }
  output_file << endl;
  
  for (uint32_t transfer=0; transfer<replicate_statistics.clade_frequencies.size(); transfer++) {
    output_file << transfer;
    vector<cGenotypeFrequency>& transfer_clade_frequencies = replicate_statistics.clade_frequencies[transfer];
    if( transfer % output_parameters.coarse_graining == 0 ) {
      for (uint32_t a_node=0; a_node<all_sweep_ids.size(); a_node++) {
        double frequency( GenotypeFrequency( transfer_clade_frequencies, all_sweep_ids[a_node]) );
        output_file << "\t" << frequency;
      }
      output_file << endl;
    }
  }
  
  output_file.close();
}

void cPopulation::OutputGenotypeFrequencies(double frequency_threshold)
{  
  vector<uint32_t> all_sweep_ids = GenotypesAboveThreshold(frequency_threshold);
  
  //Print each sweeping genotype and its frequency per time
	ofstream output_file;
  string output_file_name = output_parameters.output_directory_name + "/genotype_frequencies_" 
  + to_string(frequency_threshold) + "_" + to_string(replicate) + ".dat";
	output_file.open(output_file_name.c_str(),ios::out);
  
  cerr << "Output: " << output_file_name << endl;
  
  output_file << "transfer";
  for (uint32_t a_node=0; a_node<all_sweep_ids.size(); a_node++) {
    output_file << "\t" << "genotype_" << all_sweep_ids[a_node];
  }
  output_file << endl;
  
  
  for (uint32_t transfer=0; transfer<replicate_statistics.genotype_frequencies.size(); transfer++) {
    output_file << transfer;
    vector<cGenotypeFrequency>& transfer_genotype_frequencies = replicate_statistics.genotype_frequencies[transfer];
    if( transfer % output_parameters.coarse_graining == 0 ) {
      for (uint32_t a_node=0; a_node<all_sweep_ids.size(); a_node++) {
        double frequency( GenotypeFrequency( transfer_genotype_frequencies, all_sweep_ids[a_node]) );
        output_file << "\t" << frequency;
      }
      output_file << endl;
    }
  }
  
  output_file.close();
}

// Prints the frequency of every genotype at the current time above frequency_threshold

/* This version is NOT doing what it says
void cPopulation::OutputGenotypeFrequencies(double frequency_threshold) 
{  
  ofstream output_file;
  string output_file_name = output_parameters.output_directory_name + "/clade_frequencies_" 
    + to_string(frequency_threshold) + "_" + to_string(replicate) + ".dat";
	output_file.open(output_file_name.c_str(),ios_base::out);
  
  cerr << "Output: " << output_file_name << endl;
  
	for (uint32_t i = 0; i<replicate_statistics.clade_frequencies.size(); i++) {
		double total_freqs = 0;
    for ( vector<cGenotypeFrequency>::iterator it = replicate_statistics.genotype_frequencies[i].begin(); 
         it!=replicate_statistics.genotype_frequencies[i].end(); ++it) {
      
      //@agm set up a minimum frequency to report the print out the number so it isn't overwhelming.
      if ((*it).m_frequency > 0.01) {
       // cout << "Frequency of mutation # " << right << setw(6) << (*it).unique_node_id << " at time ";
        output_file << "Frequency of mutation # " << right << setw(6) << (*it).unique_node_id << " at time ";
      //  cout << right << setw(4) << i << " is: " << left << setw(10) << (*it).m_frequency << endl;
        output_file << right << setw(4) << i << " is: " << left << setw(10) << (*it).m_frequency << endl;
      }
      total_freqs += (*it).m_frequency;
    }
		//cout << "Round # " << i << " sum of frequencies is: " << total_freqs << endl << endl;
    output_file << "Round # " << i << " sum of frequencies is: " << total_freqs << endl << endl;
	}
} 
*/

void cPopulation::OutputMullerMatrix(uint32_t frequency_resolution)
{
  ofstream output_file;
  string output_file_name = output_parameters.output_directory_name + "/muller_matrix_" 
    + to_string(frequency_resolution) + "_" + to_string(replicate) + ".dat";
	output_file.open(output_file_name.c_str(),ios_base::out);
  
  cerr << "Output: " << output_file_name << endl;
  
  //vector< tree<cGenotype>::iterator > where(m_tree.size());
  map<uint32_t, uint32_t> renumber;
  uint32_t renumber_value(0);
  
  //double threshold(.025);
  //vector<bool> relevant_columns(newtree->size(), true);
  
  //Build vector to store all iterators in tree for parental recall later
  //@JEB should we store this iterator in cGenotypeFrequency?
  //@AGM storing transient information like iterators in an object won't work
  
  /*for (tree<cGenotype>::iterator node_loc = m_tree.begin(); node_loc != m_tree.end(); node_loc++)
   where[(*node_loc).unique_node_id] = node_loc;*/
  
  uint32_t time = 0;
  
  //step through simulation time
  for (vector< vector<cGenotypeFrequency> >::iterator this_time_freq = replicate_statistics.clade_frequencies.begin(); 
       this_time_freq < replicate_statistics.clade_frequencies.end(); ++this_time_freq) {
    cerr << " Writing Transfer: " << setw(4) << time << endl;
    time++;
    vector<cFrequencySlice> child_freqs;
    
    tree<cGenotype>::sibling_iterator location;
    location = genotype_tree.begin();
    
    AssignChildFrequency(location, 0, 1, &child_freqs, *this_time_freq);
    sort(child_freqs.begin(), child_freqs.end(), cSortByLow());
    
    /*cout << "Time " << time << " Child freqs size: " << child_freqs.size() << endl;
     for (uint32_t j=0; j<child_freqs.size(); j++) {
     cout << child_freqs[j].low << " " << child_freqs[j].high << endl;
     }
     PrintTree();*/
    
    uint32_t last_node_meeting_span;
    
    double pixel_step, min_step;
    min_step = (double) 1/frequency_resolution;
    
    //@agm Here I first iterate through the number of pixels
    for (uint32_t i=1; i<=frequency_resolution; i++) {
      
      //Determine the position of the current pixel_step
      pixel_step = (double) i/frequency_resolution;
      
      for (uint32_t j=0; j<child_freqs.size(); j++) {
        
        if( child_freqs[j].high >= pixel_step ) {
          
          //if( span > min_step ) 
          {
            
            //Add new significant mutations to a map for renumbering
            if( renumber.count(child_freqs[j].unique_node_id) == 0 ) {
              renumber.insert(make_pair(child_freqs[j].unique_node_id, renumber_value));
              renumber_value++;
              renumber_value %= frequency_resolution;
            }
            //Return the renumbered-number for the unique_node_id from the built map
            output_file<< left << setw(6) << renumber.find(child_freqs[j].unique_node_id)->second;
            last_node_meeting_span = renumber.find(child_freqs[j].unique_node_id)->second;
          }
          
          break;
        }
      }
    }
    output_file << endl;
  }
}

void cStatistics::OutputAveragePopulationFitness() {
  
	ofstream output_file;
  string output_file_name = output_directory_name + "/average_fitness.tab";
	output_file.open(output_file_name.c_str(),ios::out);
  
  cerr << "Output: " << output_file_name << endl;
  
  int32_t max_transfers_printed = 0;
  for (uint32_t replicate = 0; replicate < this->size(); ++replicate) {
    cReplicateStatistics& replicate_statistics = (*this)[replicate];
    if (replicate_statistics.average_population_fitness.size() > max_transfers_printed)
      max_transfers_printed = replicate_statistics.average_population_fitness.size();
  }
  
  output_file << "replicate";
  for (uint32_t transfer = 0; transfer < max_transfers_printed; ++transfer) {
    if( transfer % coarse_graining == 0 )
      output_file << "\t" << transfer;
  }
  output_file << endl;
  
  for (uint32_t replicate = 0; replicate < this->size(); ++replicate) {
    
    cReplicateStatistics& replicate_statistics = (*this)[replicate];
    output_file << (replicate+1);
    
    for (uint32_t transfer = 0; transfer < max_transfers_printed; ++transfer) {
      
      if( transfer % coarse_graining == 0 ) {
        
        if (transfer < replicate_statistics.average_population_fitness.size() )
          output_file << "\t" << replicate_statistics.average_population_fitness[transfer];
        else
          output_file << "\t" << "NA" << endl;

      }
    }
    output_file << endl;
  }
  output_file.close();
}



// Move time forward by the requested amount, growing all populations
// 

void cPopulation::UpdateSubpopulationsForGrowth(double update_time) 
{
  // Lineages that divided during this time step
  just_divided_lineages.clear();
  
  uint32_t i=-1;
  current_population_size = 0; // Update the population size.
  for (vector<cSubpopulation>::iterator it = current_subpopulations.begin(); it!=current_subpopulations.end(); ++it) {
    i++; // must advance iterator before continue statement
    
    // @JEB: It doesn't seem like this is possible...
    assert(it->GetNumber() != 0);

    // N = No * exp(log(2) * relative_growth_rate * t) 
    double new_number = it->GetNumber() * exp(log_2 * it->GetFitness() * update_time);     
    if (floor(new_number) - floor(it->GetNumber()) >= 1) {
      just_divided_lineages.push_back(i);
    }
    it->SetNumber(new_number);
    current_population_size += floor(new_number);
  }
  
  // Check how we have updated the population size
  if (debug) {
    int64_t calculated_population_size = CalculatePopulationSize();
    assert(calculated_population_size == current_population_size);
  }
}

double cPopulation::TimeToNextWholeCell() 
{
   double time_to_next_whole_cell = -1;
   for (vector<cSubpopulation>::iterator it = current_subpopulations.begin(); it!=current_subpopulations.end(); ++it) { 
   
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

 
//
// void cPopulation::ProcessCellDivisionTimeStep()
//
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

void cPopulation::ProcessCellDivisionTimeStep()
{
  
  if (g_verbose) {
    cout << "Completed divisions: " << num_completed_divisions << endl;
    for(vector<cSubpopulation>::iterator this_time = current_subpopulations.begin(); this_time != current_subpopulations.end(); this_time++) {
      cout << "Genotype: " << this_time->GetNode_id() << " Fitness: " << this_time->GetFitness() << " Number: " << this_time->GetNumber() << endl;
    }
  }
  
  double desired_divisions = num_divisions_until_next_mutation;
  
  // Adjust divisions to be to the end of the transfer if it is larger
  if (desired_divisions + current_population_size > simulation_parameters.final_population_size_at_transfer) {
    desired_divisions = simulation_parameters.final_population_size_at_transfer - current_population_size;
  }
  
  // Note: we underestimate by a few divisions so that we can step forward by single division increments
  // as we get close to the one where the mutation happened (or right before a transfer).      
  
  if (desired_divisions < 1) {
    desired_divisions = 1;
  }
  
  if (g_verbose) {
    cout << "Total pop size: " << current_population_size <<endl;
    cout << "Desired divisions: " << desired_divisions <<endl;
  }
  
  // How much time would we like to pass to achieve the desired number of divisions?
  // (Assuming the entire population has the maximum fitness, makes us underestimate by a few)
  // Can't change this log to ReturnLog with log table or nothing works
  
  double update_time(0);
  
  update_time = log((desired_divisions+current_population_size) / current_population_size) / maximum_subpopulation_fitness;
  
  // At a minumum, we want to make sure that one cell division took place
  
  // What is the minimum time required to get a single division?
  double time_to_next_whole_cell = TimeToNextWholeCell();
  
  if (time_to_next_whole_cell > update_time) {
    if (g_verbose) cout << "Time to next whole cell greater than update time: " << 
      time_to_next_whole_cell << " < " << update_time <<endl;    
    update_time = time_to_next_whole_cell;
  }
  
  if (g_verbose) {
    cout << "Update time: " << update_time <<endl;    
  }
  
  //Now update all lineages by the time that actually passed
  
  //FIX THIS!!!!
  //Population size shouldn't have to be recaculated here
  //I think this is fixed... the problem was in the AddSubpopulation() method
  //m_population_size = CalculatePopulationSize();
  
  double previous_population_size = current_population_size;
  
  UpdateSubpopulationsForGrowth(update_time);
  
  double completed_divisions = current_population_size - previous_population_size;
  
  if (g_verbose) {
    cout << "Completed divisions: " << completed_divisions << endl;
    for(vector<cSubpopulation>::iterator this_time = current_subpopulations.begin(); this_time != current_subpopulations.end(); this_time++) {
      //if( this_time->GetNumber() == 0 )
      cout << "Genotype: " << this_time->GetNode_id() << " Frequency: " << this_time->GetNumber() << endl;
    }
  }
  
  //cout << "Divisions until mutation: " << GetDivisionsUntilMutation() << " Completed divisions: " << GetCompletedDivisions() << " Population Size: " << GetPopulationSize() << " Previous Population Size: " << previous_population_size << endl;
  num_divisions_until_next_mutation -= completed_divisions;
}

//@JEB: Not sure what "new" is doing differently

void cPopulation::ProcessCellDivisionTimeStepNew() 
{
  
  /*
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
  double number_of_resources_to_transfer(m_final_population_size_per_transfer - m_population_size);
  
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
  */
}

void cPopulation::TransferResampleDilution() 
{
  //When it is time for a transfer, resample population
  double by_color[0];
  
  // @JEB: Move to main simulation loop
  num_completed_transfers++;
  
  if (g_verbose) cerr << ">> Transfer!" << endl;
	
  current_population_size = 0; // recalculate population size
  for (vector<cSubpopulation>::iterator it = current_subpopulations.begin(); it!=current_subpopulations.end(); ++it) {
    // Perform accurate binomial sampling only if below a certain population size
    if (it->GetNumber() < simulation_parameters.binomial_sampling_threshold) {
      if (g_verbose) cerr << "binomial " << it->GetNumber() << endl;
      it->Transfer(simulation_parameters.binomial_sampling_transfer_probability, rng);
    }
    // Otherwise, treat as deterministic and take expectation...
    else {
      it->SetNumber(it->GetNumber() * simulation_parameters.binomial_sampling_transfer_probability);
    }
    
    if (g_verbose) { 
      cerr << it->GetMarker() << endl;
      cerr << it->GetNumber() << endl;
      cerr << it->GetFitness() << endl;
    }	
    if( it->GetMarker() == 'r' ) by_color[RED] += it->GetNumber();
    else if( it->GetMarker() == 'w' ) by_color[WHITE] += it->GetNumber();
    
    //@agm This section is to delete new mutations from the tree that do not get passed, 
    //     it also deletes subpopulations from the list that have zero population
    
    current_population_size += floor(it->GetNumber());
    
    if( it->GetNumber() == 0 ) {
      it = current_subpopulations.erase(it);
      --it;
      replicate_statistics.total_subpopulations_lost++;
      existing_genotype_count--;
    }
  }
  
  // Check our calculated population size
  // Check how we have updated the population size
  if (debug) {
    int64_t calculated_population_size = CalculatePopulationSize();
    assert(calculated_population_size == current_population_size);
  }
  
  /*
   if (g_verbose) cout << "Colors: " << by_color[RED] << " / " << by_color[WHITE] << endl;
   SetRatio( (double) by_color[RED]/by_color[WHITE] );
   
   // Checks for stopping early if in marker divergence mode
   if (g_ro_only) {
   
   // One color was lost -- bail
   if ( (by_color[RED] == 0) || (by_color[WHITE] == 0) ) {
   run_end_condition_met = true;
   }
   
   // We have the minimum number of transfers to be printed and have diverged sufficiently -- bail
   else {
   if  ( (GetRatio() > GetMaxDivergenceFactor()) || 
   (GetRatio() < 1/GetMaxDivergenceFactor()) )   
   {
   if (g_verbose) cout << "DIVERGENCE CONDITION MET" << endl;
   run_end_condition_met = true;
   }  
   }
   }
   */
  
  //cout << "Pop size: " << m_population_size << endl;
}

//Currently this only works for picking one cell after resample
void cPopulation::TransferResampleExactlyOne() {
  long int starting_size(0), random_number;
  //vector<cSubpopulation> temp_populations;
  
  assert(rng);
  
  num_completed_transfers++; // @JEB should be moved to main loop
  current_population_size = 0; // recalculate population size
  
  for (vector<cSubpopulation>::iterator it = current_subpopulations.begin(); it!=current_subpopulations.end(); ++it) {
    starting_size += it->GetNumber();
  }
  
  //cout << "Population size before dilution: " << starting_size << endl;
  
  //for (uint32_t num_cells = 0; num_cells < m_initial_population_size_per_transfer; num_cells++) {
  random_number = gsl_rng_uniform_int(rng, starting_size);
  
  vector<cSubpopulation>::iterator it;
  for (it = current_subpopulations.begin(); it!=current_subpopulations.end(); ++it) {
    if( random_number < it->GetNumber() ) {
      it->SetNumber(1);
      
      // Erase empty subpopulations
      ++it;
      while(it != current_subpopulations.end()) {
        current_subpopulations.erase(it);
      }
      break;
    }
    else {
      random_number -= it->GetNumber();
      current_subpopulations.erase(it);
      --it;
    }
  }
  
  //}
  
  /*for (vector<cSubpopulation>::iterator it = m_current_subpopulations.begin(); it!=m_current_subpopulations.end(); ++it) {
   cout << it->GetNumber() << " " << it->GetNode_id() << " " << it->GetFitness() << endl;
   }*/
  
  current_population_size = 1;
}


void cPopulation::SeedPopulationWithMarkedGenotypes(uint32_t genetic_marker_states) 
{	
  // For now assumes...
  assert(genetic_marker_states == 2);
  
	cGenotype r, w;
	tree<cGenotype>::iterator top, red_side, white_side;
	double starting_fitness = simulation_parameters.initial_fitness;
	
	//initialize object of cSubpopulation type
	cSubpopulation red, white;
	
	/*@agm The unique code and fitness is set. */
  
  /*@jeb note that we increment genotype count here.
   It should really happen in AddSubpopulation?
   */
  
	r.fitness = starting_fitness;
	r.unique_node_id = unique_genotype_count++;
  w.fitness = starting_fitness;
	w.unique_node_id = unique_genotype_count++;
  
	//Start building the tree
	red_side = genotype_tree.insert(genotype_tree.begin(), r);
	white_side = genotype_tree.insert(genotype_tree.begin(), w);
	
	red.SetNumber(simulation_parameters.initial_population_size / genetic_marker_states);
	red.SetGenotype(red_side);
	red.SetMarker('r');
	
	white.SetNumber(simulation_parameters.initial_population_size / genetic_marker_states);
	white.SetGenotype(white_side);
	white.SetMarker('w');
	
	AddSubpopulation(red);
	AddSubpopulation(white);	
  
  current_population_size = simulation_parameters.initial_population_size;

}

//@agm This function seeds the population with only one colony to avoid the red/white problem
//     when possible. For this to work properly, I need to get rid of victory conditions.

//update: Victory conditions are now a command line option.  If not set to false the simulation
//        will run to the max number of iterations.

void cPopulation::SeedPopulationWithOneGenotype() {
  tree<cGenotype>::iterator start_position;
  double starting_fitness(simulation_parameters.initial_fitness);
  
  cGenotype neutral;
  neutral.fitness = starting_fitness;
  neutral.unique_node_id = unique_genotype_count++;
  neutral.mut_num = 0;
  start_position = genotype_tree.insert(genotype_tree.begin(), neutral);
  
  cSubpopulation begin_here;
  begin_here.SetNumber(simulation_parameters.initial_population_size);
  begin_here.SetGenotype(start_position);
  begin_here.SetMarker('n');
  AddSubpopulation(begin_here);
  
  current_population_size = simulation_parameters.initial_population_size;
  
  if (g_verbose) {
    cout << "Completed divisions: " << num_completed_divisions << endl;
    for(vector<cSubpopulation>::iterator this_time = current_subpopulations.begin(); this_time != current_subpopulations.end(); this_time++) {
      cout << "Genotype1: " << this_time->GetNode_id() << " Frequency1: " << this_time->GetNumber() << endl;
    }
  }
  
}

// Note that if you call add subpopulation, you may need to update the population size!
void cPopulation::AddSubpopulation(cSubpopulation& subpop) 
{
	current_subpopulations.push_back(subpop);
  existing_genotype_count++;
  
  mutations_since_last_transfer.push_back(subpop); // Keep track of all mutations since last transfer
	//cout << subpop.GetNode_id() << " " << subpop.GetFitness() << endl;
}

//@agm For some reason (that I don't fully appreciate), the way we were deleting subpopulations 
//     (based on the GetNumber()==0) Therefore, I am now deleting subpopulations Only when they 
//     arose between transfers And do not get passed to the next round just like those that are
//     removed from the tree.  This is a much more conservative pruning so I should try to figure
//     out why the other does not work.

void cPopulation::CullSubpopulationsThatDidNotEstablish() {
  
  uint32_t num_deleted(0);
  
  vector<cSubpopulation>::iterator this_mutation(mutations_since_last_transfer.begin());
  bool test(false);
  while( this_mutation != mutations_since_last_transfer.end() ) {
    for (vector<cSubpopulation>::iterator it = current_subpopulations.begin(); it!=current_subpopulations.end(); ++it) {
      if( this_mutation->GetNode_id() == it->GetNode_id() ) {
        this_mutation = mutations_since_last_transfer.erase(this_mutation);
        test = true;
        num_deleted++;
        break;
      }
    }
    
    if( test == false )
      this_mutation++;
    
    test = false;
  }
  
  this_mutation = mutations_since_last_transfer.begin();
  
  // Delete in reverse order from tree (since later can never be ancestor of earlier) this prevents orphaning subtrees
  for (vector<cSubpopulation>::reverse_iterator rit = mutations_since_last_transfer.rbegin(); rit!=mutations_since_last_transfer.rend(); ++rit) {
    genotype_tree.erase(rit->GetGenotypeIter());
  }
  //cout << "Mutations culled: " << num_mutant_subpopulations_culled << " Out of: " << m_total_mutations << endl;
  replicate_statistics.num_subpopulations_culled += num_deleted;
  mutations_since_last_transfer.clear();
  
}

//Generates new mutant and adds it to the tree
void cPopulation::Mutate() 
{	
  // we better have a random number generator
  assert(rng);
	
	if (g_verbose) cout << "* Mutating!" << endl;
  if (g_verbose) cout << "Total population: " << CalculatePopulationSize() << endl;
	
	//Mutation happened in the one that just divided
	//Break ties randomly here.
	cSubpopulation& ancestor = current_subpopulations[just_divided_lineages[rand() % just_divided_lineages.size()]];          
	cSubpopulation new_subpop;
  
  //cout << "Divided has number: " << ancestor.GetNumber() << endl;
  // There must be at least two cells for a mutation to have occurred...
  assert(ancestor.GetNumber() >= 2);
	
  //This is necessary so the first few mutation have a much larger fitness advantage
  
  if( simulation_parameters.first_beneficial_mutation_effects.size() > ancestor.GetMutNum() ) {
    new_subpop.CreateDescendant(rng, 
                                ancestor, 
                                simulation_parameters.first_beneficial_mutation_effects[ancestor.GetMutNum()], 
                                simulation_parameters.beneficial_mutation_effect_model,
                                genotype_tree,
                                unique_genotype_count++);
  }
  else {
    new_subpop.CreateDescendant(rng, 
                                ancestor, 
                                simulation_parameters.beneficial_mutation_effect, 
                                simulation_parameters.beneficial_mutation_effect_model,
                                genotype_tree,
                                unique_genotype_count++);
  }
  
	if (g_verbose) cout << "  Color: " << new_subpop.GetMarker() << endl;
	if (g_verbose) cout << "  New Fitness: " << new_subpop.GetFitness() << endl;
  
	AddSubpopulation(new_subpop);
  
  //@agm Since an existent cell is picked to mutate,
  //     rather than doubling a cell and picking its
  //     progeny to mutate, the population size does not change.
	
  
  // Update the maximum subpopulation fitness
  maximum_subpopulation_fitness = max(maximum_subpopulation_fitness, new_subpop.GetFitness());
  
  if (g_verbose) cout << "  New Number: " << new_subpop.GetNumber() << endl;
  if (g_verbose) cout << "  Old Number: " << ancestor.GetNumber() << endl;
  
  replicate_statistics.total_mutations++;
}


void cPopulation::MutateNew() {
  /*
  //@agm The code gets really finicky here
  //     if you try to change the ancestor to a pointer
  //     or a reference you will get memory leaks
  //     I have not figured out why, but it does not leak if you pass the full object
	
	if (g_verbose) cout << "* Mutating!" << endl;
  if (g_verbose) cout << "Total population: " << CalculatePopulationSize() << endl;
  
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
    
    //cout << "Divided has number: " << ancestor.GetNumber() << endl;
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
    
    if (g_verbose) cout << "  Color: " << new_subpop.GetMarker() << endl;
    if (g_verbose) cout << "  New Fitness: " << new_subpop.GetFitness() << endl;
    
    AddSubpopulation(new_subpop);
    
    //@agm Since an existent cell is picked to mutate,
    //     rather than doubling a cell and picking its
    //     progeny to mutate... we should not add a new
    //     cell to the population, but AddSubpopulation does.
    m_population_size-=1;
    
    if(new_subpop.GetFitness() > GetMaxW()) SetMaxW(new_subpop.GetFitness());
    
    number_of_cells--;
    
    ancestor.SetNumber(number_of_cells);
    
    if (g_verbose) cout << "  New Number: " << new_subpop.GetNumber() << endl;
    if (g_verbose) cout << "  Old Number: " << ancestor.GetNumber() << endl;
    
    m_total_mutations++;
  }
  */
}

//////////  Methods for recording statistics //////////

// Pushes the current frequencies of all genotypes to the run_statistics.clade_frequencies vector
//
//@agm Here I iterate through the populations in m_current_subpopulations, build a vector to store the values 
//     for the sizes of each node by iterating up, and divide the size of each node (subpopulation) by the total 
//     population size. This should give the relative frequency of a given unique_node_id in the population.

void cPopulation::RecordStatisticsAtTransfer()
{	
  //  Debug: this checks to see if our population size was correct.
  if (debug) {
    uint64_t calculated_population_size = CalculatePopulationSize();
    assert(current_population_size == calculated_population_size);
  }
  
  /////////////////////////////////////////
  /// Statistic: genotype_frequencies, 
  /// Statistic: clade_frequencies
  /////////////////////////////////////////
  
  // Initialize all entries
  map<uint32_t,cGenotypeFrequency> current_genotype_frequency_map; // only of exactly this genotype 
  map<uint32_t,cGenotypeFrequency> current_clade_frequency_map;    // AND including all descendants

	tree<cGenotype>::iterator update_location, genotype;
  for (uint32_t it = 0; it<current_subpopulations.size(); ++it) {
    
    cSubpopulation& this_pop = current_subpopulations[it];
    
    update_location = this_pop.GetGenotypeIter();
    assert(update_location != NULL);
    
    
    cGenotypeFrequency& current_genotype_frequency = current_genotype_frequency_map[update_location->unique_node_id];
    current_genotype_frequency.unique_node_id = update_location->unique_node_id;
    current_genotype_frequency.m_frequency = floor(this_pop.GetNumber()) / current_population_size;
    
    // Add frequency to all parents to get clade numbers
		while(update_location != NULL) {
      
      cGenotypeFrequency& current_clade_frequency = current_clade_frequency_map[update_location->unique_node_id];
      current_clade_frequency.unique_node_id = update_location->unique_node_id;
      current_clade_frequency.m_frequency += floor(this_pop.GetNumber());
      
      //current_clade_frequencies[update_location->unique_node_id].unique_node_id = update_location->unique_node_id;
      //current_clade_frequencies[update_location->unique_node_id].m_frequency += floor(this_pop.GetNumber());
      
      update_location = genotype_tree.parent(update_location);
		}    
    
  }
  
  //
  // Convert maps to vectors (but really, why not leave them as maps...)
  //
  
  vector<cGenotypeFrequency> current_genotype_frequencies; 
	vector<cGenotypeFrequency> current_clade_frequencies;   
  
  for( map<uint32_t,cGenotypeFrequency>::iterator it = current_genotype_frequency_map.begin(); it != current_genotype_frequency_map.end(); ++it ) {
    current_genotype_frequencies.push_back( it->second );
  }
  
  for( map<uint32_t,cGenotypeFrequency>::iterator it = current_clade_frequency_map.begin(); it != current_clade_frequency_map.end(); ++it ) {
    it->second.m_frequency = it->second.m_frequency / current_population_size;
    current_clade_frequencies.push_back( it->second );
  }
  
  
  
  // Push to record
	replicate_statistics.genotype_frequencies.push_back(current_genotype_frequencies);
	replicate_statistics.clade_frequencies.push_back(current_clade_frequencies);
  
   
  /////////////////////////////////////////
  /// Statistic: average_population_fitness
  /////////////////////////////////////////
  
  double total_fitness(0);
  for (vector<cSubpopulation>::iterator it = current_subpopulations.begin(); it!=current_subpopulations.end(); ++it) {    
    total_fitness += it->GetFitness()*floor(it->GetNumber());
  }
  replicate_statistics.average_population_fitness.push_back(total_fitness/current_population_size);
}


//////////  Utility Functions //////////


// Calculates the total population size by iterating over subpops. 
// This is really for debug checking of code only. Normally, use current_population_size.
uint64_t cPopulation::CalculatePopulationSize() 
{
  uint64_t calculated_population_size = 0;
  for (vector<cSubpopulation>::iterator it = current_subpopulations.begin(); it!=current_subpopulations.end(); ++it) {
    calculated_population_size += floor(it->GetNumber());
  }  
  return calculated_population_size;
}

// Calculates the maximum population fitness by iterating over subpops. 
// This is really for debug checking of code only. Normally, use maximum_subpopulation_fitness.

double cPopulation::CalculateMaximumSubpopulationFitness() {
  double max_fitness = 0.0;
  for(vector<cSubpopulation>::iterator it=current_subpopulations.begin(); it!=current_subpopulations.begin(); ++it) {
    max_fitness = max(max_fitness, it->GetFitness());
  }
  return max_fitness;
}



// Helper function for producing Muller plot matrix
// Assigns bounds to child genotype centered (and equally dividing) parent genotype slice.
double cPopulation::AssignChildFrequency(tree<cGenotype>::sibling_iterator this_node,
                                              double in_low,
                                              double in_high,
                                              vector<cFrequencySlice> * child_freqs,
                                              vector<cGenotypeFrequency> &frequencies, 
                                              int depth)  // current depth in tree, defaults to zero
{
  
  //kptree::print_tree_bracketed(*newtree);
  
  //The low for this mutation should be the low of the input interval
  double this_low = in_low;
  double this_high = this_low + GenotypeFrequency(frequencies, this_node->unique_node_id);
  double size_depth1_children(0), half_size_parent_swath((this_high-this_low)/2);
  
  for (tree<cGenotype>::sibling_iterator it_node = genotype_tree.begin(this_node); it_node!=genotype_tree.end(this_node); ++it_node) {
    size_depth1_children += GenotypeFrequency(frequencies, it_node->unique_node_id);
  }
  /*if(size_depth1_children != 0)
   cout << size_depth1_children << endl;*/
  
  double this_bottom_high(this_low + half_size_parent_swath - (size_depth1_children/2));
  double this_top_low(this_bottom_high);
  
  // The swath for this mutation may shrink on the bottom
  // due to its children taking a bite out of it.  
  double last_assigned_child_high = this_top_low;
  for (tree<cGenotype>::sibling_iterator it_node = genotype_tree.begin(this_node); it_node!=genotype_tree.end(this_node); ++it_node) {
    
    // is the frequency > 0? (It may have gone extinct, in which case it is a waste to keep going down the tree!)
    if( GenotypeFrequency(frequencies, it_node->unique_node_id) > 0 ) {
      last_assigned_child_high = AssignChildFrequency(it_node, this_top_low, this_high, child_freqs, frequencies, depth+1);
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
   cout << " ";
   }*/
  //cout << " ID:" << this_node->unique_node_id << " Freq:" << (*frequencies)[this_node->unique_node_id].frequency;
  //cout << " [" << this_low << "," << this_high << "]" << endl;
  
  return this_high;  
}





// JEB: This is leftover from Austin's MutationsAboveThreshold. 
// Not exactly sure what it is doing...
// It appears to collect all genotypes leading to the final dominant and return their id's

vector<uint32_t> cPopulation::GenotypesFromAncestorToFinalDominant(float threshold) {
  
  uint32_t youngest_sweep(Last_Sweep(threshold));
  
  vector<uint32_t> all_sweep_ids;
  
  for (tree<cGenotype>::iterator it = genotype_tree.begin(); it!=genotype_tree.end(); it++) {
    if( it->unique_node_id == youngest_sweep ) {
      while (it != NULL) {
        all_sweep_ids.push_back(it->unique_node_id);
        it = genotype_tree.parent(it);
      }
      break;
    }
  }
  
  sort(all_sweep_ids.begin(), all_sweep_ids.end());
  return all_sweep_ids;
}

// Returns a list of all clades that were ever above the
// requested threshold at any transfer during the simulation.

vector<uint32_t> cPopulation::CladesAboveThreshold(float threshold) {
  
  set<uint32_t> genotype_set;
  
  for (uint32_t transfer=0; transfer < replicate_statistics.clade_frequencies.size(); transfer++)
  {
    vector<cGenotypeFrequency>& transfer_clade_frequencies = replicate_statistics.clade_frequencies[transfer];
    
    for (uint32_t i = 0; i<transfer_clade_frequencies.size(); i++) {
      cGenotypeFrequency& gf = transfer_clade_frequencies[i];
      if (gf.m_frequency >= threshold)
        genotype_set.insert(gf.unique_node_id);
    }
  }
  
  vector<uint32_t> genotype_list(genotype_set.begin(), genotype_set.end());
  sort(genotype_list.begin(), genotype_list.end());
  return genotype_list;
}

// Returns a list of all genotypes that were ever above the
// requested threshold at any transfer during the simulation.

vector<uint32_t> cPopulation::GenotypesAboveThreshold(float threshold) {
  
  set<uint32_t> genotype_set;
  
  for (uint32_t transfer=0; transfer < replicate_statistics.genotype_frequencies.size(); transfer++)
  {
    vector<cGenotypeFrequency>& transfer_genotype_frequencies = replicate_statistics.genotype_frequencies[transfer];
    
    for (uint32_t i = 0; i<transfer_genotype_frequencies.size(); i++) {
      cGenotypeFrequency& gf = transfer_genotype_frequencies[i];
      if (gf.m_frequency >= threshold)
        genotype_set.insert(gf.unique_node_id);
    }
  }
  
  vector<uint32_t> genotype_list(genotype_set.begin(), genotype_set.end());
  sort(genotype_list.begin(), genotype_list.end());
  return genotype_list;
}


uint32_t cPopulation::Last_Sweep(float threshold) {
  
  for( vector< vector<cGenotypeFrequency> >::reverse_iterator time = replicate_statistics.clade_frequencies.rbegin(); 
      time < replicate_statistics.clade_frequencies.rend(); ++time ) {
    for( std::vector<cGenotypeFrequency>::reverse_iterator node = time->rbegin(); node < time->rend(); ++node) {
      if( node->m_frequency >= threshold )
        return node->unique_node_id;
    }
  }
  
  cout << "Nothing met threshold." << std::endl;
  
  return 0;
}


//Used for MA?

void cPopulation::PrintSingleFitness(string output_folder) {
  string output_file;
  ofstream output_handle;
  
  float single_fitness;
  
  output_file = output_folder + "/SingleFitness.dat";
  output_handle.open(output_file.c_str(),ios_base::app);
  
  if(current_subpopulations.size() == 1) {
    for (vector<cSubpopulation>::iterator it = current_subpopulations.begin(); it!=current_subpopulations.end(); ++it) { 
      single_fitness = it->GetFitness();
    }
  }
  else {
    cerr << "There is more than one population remaining... there is a bug." << endl;
    abort();
  }
  
  output_handle << single_fitness << "\n";
  //if( single_fitness != 1 )
    //cout << "Fitness: " << single_fitness << "\n";
  output_handle.close();
}

void cPopulation::ConvertExternalData(const string &input_file) {
  
  fstream input_handle;
  input_handle.open(input_file.c_str(), fstream::in);
  assert(input_handle.is_open());
  
  char char_line[1000];
  
  uint32_t num_nodes = 1;
  
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
      
      top = genotype_tree.begin();
      genotype_tree.insert(top, origin);
      
      replicate_statistics.clade_frequencies.resize(line_pieces.size()-1);
      
      for(uint32_t this_time = 0; this_time < replicate_statistics.clade_frequencies.size(); ++this_time) {
        if( from_string<double>(line_pieces[this_time+1]) != 0 ) {
          cGenotypeFrequency temp_time;
          temp_time.unique_node_id = num_nodes;
          temp_time.name = *(preformat_tree.end()-1);
          temp_time.m_frequency = from_string<double>(line_pieces[this_time+1]);
          replicate_statistics.clade_frequencies[this_time].push_back(temp_time);
        }
      }
    }
    else {
      for(tree<cGenotype>::iterator this_node = genotype_tree.begin(); this_node != genotype_tree.end(); ++this_node) {
        if( *(preformat_tree.end()-2) == (*this_node).name ) {
          cSubpopulation new_child;
          cGenotype child;
          child.name = *(preformat_tree.end()-1);
          child.unique_node_id = num_nodes;
          child.fitness = 1;
          child.mut_num = num_nodes;
          new_child.AddToTree(genotype_tree, this_node, child);
          
          for(uint32_t this_time = 0; this_time < replicate_statistics.clade_frequencies.size(); ++this_time) {
            if( from_string<double>(line_pieces[this_time+1]) != 0 ) {
              cGenotypeFrequency temp_time;
              temp_time.unique_node_id = num_nodes;
              temp_time.name = *(preformat_tree.end()-1);
              temp_time.m_frequency = from_string<double>(line_pieces[this_time+1]);
              replicate_statistics.clade_frequencies[this_time].push_back(temp_time);
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



vector<cGenotypeFrequency>::iterator cPopulation::Find_Node_in_Freq(vector<cGenotypeFrequency> &frequencies, tree<cGenotype>::sibling_iterator this_node) {
  
  cGenotypeFrequency return_node;
  return_node.unique_node_id = this_node->unique_node_id;
  
  for( vector<cGenotypeFrequency>::iterator a_node = frequencies.begin(); a_node < frequencies.end(); ++a_node ){
    if( a_node->unique_node_id == this_node->unique_node_id ) {
      return a_node;
    }
    else
      return_node.m_frequency = 0;
  }
  
  vector<cGenotypeFrequency> return_vector;
  return_vector.push_back(return_node);
  
  return return_vector.begin();
}

double cPopulation::GenotypeFrequency(vector<cGenotypeFrequency> &frequencies,
                                                uint32_t this_node) {
  for( vector<cGenotypeFrequency>::iterator a_node = frequencies.begin(); a_node < frequencies.end(); ++a_node ) {
    if( a_node->unique_node_id == this_node )
      return a_node->m_frequency;
  }
  
  return 0;
}

cSubpopulation* cPopulation::Find_Node_in_Populations_By_NodeID(uint32_t this_node) {
  for( uint32_t it = 0; it < current_subpopulations.size(); ++it ) {
    if( current_subpopulations[it].GetNode_id() == this_node ) {
      return &current_subpopulations[it];
    }
  }
  assert(0==1);
}





//Utilities Section

vector<uint32_t> cPopulation::CurrentUniqueGenotypes() {
  
  /*
  vector<uint32_t> number_of_unique_genotypes;
  
  for (uint32_t time = 0; time<m_all_subpopulations_at_all_times.size(); time++) {
    uint32_t current_number(0);
    for (uint32_t this_one = 0; this_one < sizeof(m_all_subpopulations_at_all_times[time]); this_one++) {
      double& this_pop = m_all_subpopulations_at_all_times[time][this_one];
      if( (double) this_pop / m_total_cells[time] > .1 ) current_number++;
    }
    number_of_unique_genotypes.push_back(current_number);
  }
  return number_of_unique_genotypes;
}

void cPopulation::PrintUniqueGenotypes(const string& output_folder,
                                       vector< vector<uint32_t> > * number_of_unique_genotypes) {
  //Will print out only red and white
	ofstream output_handle;
  string output_file;
  
  output_file = output_folder + "/Number_Unique_Genotypes.dat";
  
	output_handle.open(output_file.c_str(),ios_base::app);
  
  uint32_t largest_replicate(0);
  for (uint32_t replicate=0; replicate<(*number_of_unique_genotypes).size(); replicate++) {
    if( (*number_of_unique_genotypes)[replicate].size() > (*number_of_unique_genotypes)[largest_replicate].size() ) largest_replicate = replicate;
  }
  
  output_handle << "transfer";
  for (uint32_t replicate = 0; replicate<(*number_of_unique_genotypes).size(); replicate++) {
    output_handle << "\t" << replicate ;
  }
  
  output_handle << endl;
  
  for (uint32_t time = 0; time<(*number_of_unique_genotypes)[largest_replicate].size(); time++) {
    output_handle << time;
    for (uint32_t replicate = 0; replicate<(*number_of_unique_genotypes).size(); replicate++) {
      if( time >= (*number_of_unique_genotypes)[replicate].size() )
        output_handle << "\t";
      else
        output_handle << "\t" << (*number_of_unique_genotypes)[replicate][time];
    }
    output_handle << endl;
	}
 */ 
}

//The goal here is the find the fitness of the winning line of descent at every time point
//This allows comparison eventualy winners and losers style

void cPopulation::PrintWinningFitness(string output_folder, uint32_t run_number) {
  
/*
  vector<uint32_t> all_sweep_ids(GenotypesFromAncestorToFinalDominant(1.0));
  
  ofstream output_handle;
  string output_file;
  
  output_file.append(output_folder);
  output_file.append("/Winning_Fitness_");
  output_file.append(to_string(run_number));
  output_file.append(".dat");
	output_handle.open(output_file.c_str(),ios_base::app);
  
  uint32_t count(0);
  
  //First I iterate over the frequency vector per time
  for (vector< vector<cGenotypeFrequency> >::iterator this_time = m_frequencies.begin(); this_time < m_frequencies.end(); ++this_time) {
    if( count % m_coarse_graining == 0 ) {
      
      for (uint32_t a_node=0; a_node<all_sweep_ids.size(); a_node++) {
        
        //Here we grab the frequency of each of the nodes that are known to have swept by the END of the simulation
        double frequency( GenotypeFrequency( *this_time, all_sweep_ids[a_node]) );
        
        //On the first node that has zero frequency,
        //Find the fitness in the tree of the one that came immediately before it
        if( frequency == 0 ) {
          for (tree<cGenotype>::iterator it = m_tree.begin(); it!=m_tree.end(); it++) {
            if( it->unique_node_id == all_sweep_ids[a_node-1] )
               output_handle << it->fitness << "\t";
               
          }
          //After finding the youngest one at this time point in the winning line...
          //Go to the next time point
          break;
        }
      }
    }
    count++;
  }
  
  output_handle.close();
 */
}

void cPopulation::PrintExpectationValue(const string& output_folder) {
  /*
  //Print everything out
	ofstream output_handle;
  string output_file;
  
  output_file.append(output_folder);
  output_file.append("/ExpectationVal.dat");
	output_handle.open(output_file.c_str(),ios_base::app);
  
  double calculated_mutant_population_freq(0);

  for ( vector<cGenotypeFrequency>::iterator it = m_frequencies[m_frequencies.size()-1].begin(); it!=m_frequencies[m_frequencies.size()-1].end(); ++it) {
    if( it->unique_node_id != 0 ) {
      calculated_mutant_population_freq += it->m_frequency;
    }
  }
    
  cout << calculated_mutant_population_freq << endl;
  */
}

void cPopulation::PrintOut_RedWhiteOnly(const string& output_folder, 
                                        vector< vector<double> > * red_white_ratios,
                                        uint32_t transfer_interval_to_print) 
{  
	//Will print out only red and white
	ofstream output_handle;
  string output_file;
  
  output_file = output_folder + "_Gfreqs_RO.dat";
  
	output_handle.open(output_file.c_str(),ios_base::app);
  
  uint32_t largest_replicate(0);
  //find largest replicate
  for (uint32_t replicate = 0; replicate<(*red_white_ratios).size(); replicate++) {
    if( (*red_white_ratios)[replicate].size() > largest_replicate ) largest_replicate = replicate;
  }
  
  output_handle << "transfer";
  for (uint32_t replicate = 0; replicate<(*red_white_ratios).size(); replicate++) {
    output_handle << "\t" << replicate ;
  }

  output_handle << endl;
  
  for (uint32_t time = 0; time<(*red_white_ratios)[largest_replicate].size(); time++) {
    output_handle << time*transfer_interval_to_print;
    for (uint32_t replicate = 0; replicate<(*red_white_ratios).size(); replicate++) {
      if( time >= (*red_white_ratios)[replicate].size() )
        output_handle << "\t";
      else
        output_handle << "\t" << log((*red_white_ratios)[replicate][time]);
    }
    output_handle << endl;
	}
}

void cPopulation::PrintFrequenciesToScreen_RedWhiteOnly(string output_folder) {;}

//@agm This function determines the maximum difference in genotype frequency between a mutation
//     and its predecessor if they got above the threshold passed to the threshold function below

double cPopulation::CalculateSimilarity(string output_folder) {
  
  /*
  vector<uint32_t> all_sweep_ids(GenotypesFromAncestorToFinalDominant(.99));
  
  ofstream output_handle;
  string output_file;
  
  vector<float> max_diff(all_sweep_ids.size(), 0);
  float current_diff;
  double num_below_threshold(0);
  
  uint32_t i=0;
  
  for (uint32_t a_node=0; a_node<all_sweep_ids.size(); a_node++) {
    uint32_t time = 0;
    for ( vector<cGenotypeFrequency>::iterator it = run_statistics.clade_frequencies[i].begin(); it!=run_statistics.clade_frequencies[i].end(); ++it) {
      
      //This conditional is to coarse-grain the sampling to every 75 transfers (~500 generations)
      //and to make sure there is a value in the vector at that position
      //another way to think about it is this conditional is checking to make sure that
      //a particular mutation has arisen in time before comparing it to something.
      
      if ( time% output_parameters.coarse_graining == 0 ) {
        current_diff = fabs(GenotypeFrequency(*this_time, all_sweep_ids[a_node]) - GenotypeFrequency(*this_time, all_sweep_ids[a_node+1]));
        if( max_diff[i] < current_diff ) max_diff[i] = current_diff;
      }
      time++;
    }
    i++;
  }
  
  output_file = output_folder + "/" + "SweepClumpiness.dat";
  output_handle.open(output_file.c_str(), ios_base::app);
  
  uint32_t num_simultaneous(1);

  for (int i = 0; i<max_diff.size()-1; i++) {
    if( max_diff[i] <= .1 ) num_simultaneous++; 
    else {
      output_handle << num_simultaneous << endl;
      num_simultaneous = 1;
    }
  }
  
  output_handle.close();
  
  output_file = output_folder + "/" + "SignificantParallelMutations.dat";
	output_handle.open(output_file.c_str(), ios_base::app);
  
  for (int i = 0; i<max_diff.size()-1; i++) {
    cout << endl << i << " " << max_diff[i] << endl;
    output_handle << max_diff[i] << endl;
    if( max_diff[i] <= .1 ) num_below_threshold++;
  }
  
  output_handle.close();
  return num_below_threshold;
   */
}

void cPopulation::TimeToSweep(string output_folder) {
  
  /*
  vector<uint32_t> all_sweep_ids(GenotypesFromAncestorToFinalDominant(.98));
  
  ofstream output_handle;
  string output_file;
  
  vector<uint32_t> time_to_sweep(all_sweep_ids.size(), 0);
  
  uint32_t i=0;
  
  for (uint32_t a_node=0; a_node<all_sweep_ids.size(); a_node++) {
    uint32_t time = 0;
    for (vector< vector<cGenotypeFrequency> >::iterator this_time = m_frequencies.begin(); this_time < m_frequencies.end(); ++this_time) {
      
      //This conditional is to coarse-grain the sampling to every 75 transfers (~500 generations)
      //and to make sure there is a value in the vector at that position
      //another way to think about it is this conditional is checking to make sure that
      //a particular mutation has arisen in time before comparing it to something.
      
      double frequency( GenotypeFrequency(*this_time, all_sweep_ids[a_node]) );
      
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
	output_handle.open(output_file.c_str(),ios_base::app);
  
  for (vector<uint32_t>::iterator it_node = time_to_sweep.begin(); it_node != time_to_sweep.end(); it_node++) {
    if( (*it_node) != 0 ) {
      cout << (*it_node) << endl;
      output_handle << (*it_node) << endl;
    }
  }
}

vector<uint32_t> cPopulation::GenotypesFromAncestorToFinalDominant(float threshold) {
  
  uint32_t youngest_sweep(Last_Sweep(threshold));
  
  vector<uint32_t> all_sweep_ids;
  
  for (tree<cGenotype>::iterator it = m_tree.begin(); it!=m_tree.end(); it++) {
    if( it->unique_node_id == youngest_sweep ) {
      while (it != NULL) {
        all_sweep_ids.push_back(it->unique_node_id);
        it = m_tree.parent(it);
      }
      break;
    }
  }
  
  sort(all_sweep_ids.begin(), all_sweep_ids.end());
  
  return all_sweep_ids;
}

uint32_t cPopulation::Last_Sweep(float threshold) {
  
  for( vector< vector<cGenotypeFrequency> >::reverse_iterator time = m_frequencies.rbegin(); time < m_frequencies.rend(); ++time ) {
    for( vector<cGenotypeFrequency>::reverse_iterator node = time->rbegin(); node < time->rend(); ++node) {
      if( node->m_frequency >= threshold )
        return node->unique_node_id;
    }
  }
  
  cout << "Nothing met threshold." << endl;
  
  return 0;
   */
}




////////// THINGS WE DONT REALLY NEED BELOW


//@agm I wrote this based on the taylor series expansion... aren't we impressed
//     unfortunately it is significantly slower than the built in log function
//     I kept it around to test the effect of various precisions on the simulation
//     To get better precision change the iterations to whatever number you like
void cPopulation::PrintFreqsQuick() {
  uint32_t time_keeper(0);
  for(vector< vector<cGenotypeFrequency> >::iterator this_time=replicate_statistics.clade_frequencies.begin(); 
      this_time!=replicate_statistics.clade_frequencies.end(); ++this_time) {
    time_keeper++;
    cout << "Time: " << time_keeper << endl;
    for(vector<cGenotypeFrequency>::iterator this_genotype=this_time->begin(); this_genotype!=this_time->end(); ++this_genotype) {
      cout << this_genotype->unique_node_id << " " << this_genotype->name << " " << this_genotype->m_frequency << endl;
    }
    cout << endl << endl;
  }
}

float cPopulation::Logarithm(float mantissa) {
  float value(0), num, topower;
  uint32_t iterations(2);
  
  for (uint32_t i = 0; i<iterations; i++) {
    topower = 1;
    num = (2*i)+1;
    for(uint32_t j = 0; j<num; j++) topower *= ((mantissa-1)/(mantissa+1));
    value += (1/num)*topower;
  }
  return 2*value;
}