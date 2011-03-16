#include "cPopulation.h"
#include "cLineageTree.h"
#include "tree_util.hh"

void cPopulation::UpdateLineages() 
{
  for (std::vector<cSubpopulation>::iterator it = m_populations.begin(); it!=m_populations.end(); ++it) {
     if (it->GetNumber() == 0) continue;
     it->SetNumber(it->GetNumber() * exp(log(2) * GetUpdateTime() * it->GetFitness()));
     SetNewPopSize(GetNewPopSize() + it->GetNumber());
  }
}

void cPopulation::DetermineDivisionTime() 
{
   int itCount = 0;
	
   for (std::vector<cSubpopulation>::iterator it = m_populations.begin(); it!=m_populations.end(); ++it) { 
      if(it->GetNumber() == 0) continue;
      //what is the time to get to the next whole number of cells?
     
      SetCurrentCells(it->GetNumber());

      SetWholeCells(floor(it->GetNumber())+1);
      // WC = N * exp(growth_rate * t) 

      SetThisTimeToNextWholeCell(log(GetWholeCells() / GetCurrentCells()) / (it->GetFitness()));   

      if ( GetTimeToNextWholeCell() == 0 || (GetThisTimeToNextWholeCell() < GetTimeToNextWholeCell()) ) {
        m_divided_lineages.clear();
        SetTimeToNextWholeCell(GetThisTimeToNextWholeCell());
        m_divided_lineages.push_back(itCount); //a list, because there can be ties 
      }
      else if (GetThisTimeToNextWholeCell() == GetTimeToNextWholeCell()) {
        m_divided_lineages.push_back(itCount); //a list, because there can be ties
      }
    itCount++;
  }
}

//Here I try to iterate through the populations in m_population twice
//The first time, I simply sum the number extant populations
//The second time, I iterate through and count all of the children of 
//each node, then I divide that number of children by the total population size
//This should give the relative frequence of a given unique_node_id in the 
//population.

//Later I plan to put this information into a vector of vectors

void cPopulation::FrequenciesPerTransferPerNode(tree<cGenotype> newtree, 
																								std::vector< std::vector<double> >& frequencies)
{
	tree<cGenotype>::iterator loc;
	int counter;
	
	for (std::vector<cSubpopulation>::iterator it = m_populations.begin(); it!=m_populations.end(); ++it) {
		counter = 0;
		loc = it -> GetGenotypeIter();
		tree<cGenotype>::iterator node_children = newtree.begin(loc);
		
		//std::cout << (*loc).unique_node_id << ", ";
		
		if (node_children != newtree.end(loc)) std::cout << std::endl << "Children of " << (*loc).unique_node_id << ":" << std::endl;
	  while(node_children != newtree.end(loc)) {
			std::cout << (*node_children).unique_node_id << ", ";
			node_children++;
			counter++;
		}
		
		tree<cGenotype>::sibling_iterator start_loc2 = newtree.begin(loc);
		if (start_loc2 != newtree.end(loc)) std::cout << std::endl << "Number of Children Above: " << counter << std::endl;
	}
}


void cPopulation::Resample(gsl_rng * randgen) 
{
  //When it is time for a transfer, resample population
		
	 if (GetVerbose()) std::cout << ">> Transfer!" << std::endl;
	 //Is there an exists() in C++?
	 m_by_color[RED] = 0;
	 m_by_color[WHITE] = 0;
	 SetNewPopSize(0);        
	 for (std::vector<cSubpopulation>::iterator it = m_populations.begin(); it!=m_populations.end(); ++it) {
			if (it->GetNumber() == 0) continue;       
			// Perform accurate binomial sampling only if below a certain population size
			if (it->GetNumber() < GetBinomialSamplingThreshold()) {
				 if (GetVerbose()) std::cout << "binomial " << it->GetNumber() << std::endl;
				 it->Transfer(GetTransferBinomialSamplingP(), randgen);
			}
			// Otherwise, treat as deterministic and take expectation...
			else {
				 it->SetNumber(it->GetNumber() * GetTransferBinomialSamplingP());
			}
					
			// Keep track of lineages we lost
			if (it->GetNumber() == 0) {  
				 SetTotalSubpopulationsLost(GetTotalSubpopulationsLost()+1);
			}
			
			SetNewPopSize(GetNewPopSize() + floor(it->GetNumber()));
			if (GetVerbose()) std::cout << it->GetMarker() << std::endl;
			if (GetVerbose()) std::cout << it->GetNumber() << std::endl;
			if (GetVerbose()) std::cout << it->GetFitness() << std::endl;
			if (it->GetMarker() == 'r') m_by_color[RED] += it->GetNumber();
			else m_by_color[WHITE] += it->GetNumber();
	 }
	 //One color was lost, bail  
	 if ( (m_by_color[RED] == 0) || (m_by_color[WHITE] == 0) ) {
			SetKeepTransferring(false);
	 }     
	 if (GetVerbose()) std::cout << "Colors: " << m_by_color[RED] << " / " << m_by_color[WHITE] << std::endl;
	 SetRatio(m_by_color[RED] / m_by_color[WHITE]);
	 SetTransfers(GetTransfers()+1);

	 if ( (GetTransfers() >= 0) && (GetTransfers() % GetTransferIntervalToPrint() == 0) ) {  
			m_this_run.push_back(GetRatio());
		 
			if (GetVerbose() == 1) { 
				 std::cout << "Transfer " << GetTransfers() << " : " << GetTotalPopSize() << 
				 "=>" << GetNewPopSize() << "  R/W Ratio: " << GetRatio() << std::endl;  
				 std::cout << "Total mutations: " << GetTotalMutations() << " Maximum Fitness: " << GetMaxW() << std::endl;
				 std::cout << "Size = " << m_this_run.size() << std::endl;
			}
	 }  
	 SetTotalPopSize(GetNewPopSize());
	 if ( (int(m_this_run.size()) >= GetMinimumPrinted()) && ((GetRatio() > GetMaxDivergenceFactor()) || (GetRatio() < 1/GetMaxDivergenceFactor())) )   
	 {
			if (GetVerbose()) std::cout << "DIVERGENCE CONDITION MET" << std::endl;
			SetKeepTransferring(false);
	 }  
  
}

void cPopulation::PushBackRuns()
{
   m_runs.push_back(m_this_run);
}

void cPopulation::PrintOut(const std::string& output_file_name)
{  

   //Print everything out
   std::ofstream output_file;
   output_file.open(output_file_name.c_str(),std::ios_base::app);
   output_file << "transfer";
   for (int on_run=0; on_run < GetReplicates(); on_run++) 
   {
     output_file << "\t" << on_run;
   }
   output_file << std::endl;

   bool still_going = true;
   int i = 0;
   while (still_going) 
   {    
      int on_transfer = i * GetTransferIntervalToPrint();
      output_file << on_transfer;
      still_going = false;
      for (int on_run=0; on_run < GetReplicates(); on_run++) 
      {      
         if(i<int(m_runs[on_run].size())) 
         {
            output_file << "\t" << m_runs[on_run][i];
            still_going = true;
         }
         else 
         {
            output_file << "\t";
         }
      }
      output_file << std::endl;
      i++;
   }
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

   SetTotalPopSize(GetInitialPopulationSize());
   SetTotalMutations(0);    
   SetTotalSubpopulationsLost(0);
   SetTransfers(1);    
   SetDivisionsUntilMutation(0);
   SetKeepTransferring(true);

}

void cPopulation::SetParameters(const variables_map &options)
{

  SetGrowthPhaseGenerations(
    options.count("generations-per-transfer") ?
    options["generations-per-transfer"].as<double>() : 6.64
  );

  SetPopSizeAfterDilution(
    options.count("population-size-after-transfer") ?
    options["population-size-after-transfer"].as<uint64_t>() : int(5E6)
  );  
  SetInitialPopulationSize(
    options.count("initial-population-size") ?
    options["initial-population-size"].as<uint64_t>() : int(2)
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
    options.count("total-transfers") ?
    options["total-transfers"].as<int>() : 200
  );  
  SetMaxDivergenceFactor(
    options.count("marker-divergence") ?
    options["marker-divergence"].as<int>() : 100
  );  
  SetReplicates(
    options.count("replicates") ?
    options["replicates"].as<int>() : 10
  );  
  SetMinimumPrinted(
    options.count("minimum-printed") ?
    options["minimum-printed"].as<int>() : 8
  );
  SetBeneficialMutationDistribution(
    options.count("type-of-mutations") ?
    options["type-of-mutations"].as<char>() : 'u'
  );
  SetLineageTree(
    options.count("lineage-tree") ?
    options["lineage-tree"].as<int>() : 1
  );

  // Simulation parameters that are pre-calculated
  SetDilutionFactor(exp(log(2) * GetGrowthPhaseGenerations()));
  SetTransferBinomialSamplingP(1/GetDilutionFactor());
  SetPopSizeBeforeDilution(GetPopSizeAfterDilution() * GetDilutionFactor());
  SetLambda(1/GetMutationRatePerDivision());
  SetBinomialSamplingThreshold(1000);
}

void cPopulation::DisplayParameters()
{
  if (GetVerbose()==1) 
  {
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

  SetDesiredDivisions(GetDivisionsUntilMutation());
  if (GetDesiredDivisions() + GetTotalPopSize() > GetPopSizeBeforeDilution())
  {
     SetDesiredDivisions(GetPopSizeBeforeDilution() - GetTotalPopSize());
  }
  
  if(GetVerbose() == 1) 
      {
         std::cout << "Divisions before next mutation: " << GetDivisionsUntilMutation() <<std::endl;
      }
  
  // Note: we underestimate by a few divisions so that we can step forward by single division increments
  // as we get close to the one where the mutation happened (or right before a transfer).      
  
  if (GetDesiredDivisions() < 1)
  {
     SetDesiredDivisions(1);
  }
  
  if (GetVerbose() == 1) 
  {
     std::cout << "Total pop size: " << GetTotalPopSize() <<std::endl;
     std::cout << "Desired divisions " << GetDesiredDivisions() <<std::endl;
  }
      
   // How much time would we like to pass to achieve the desired number of divisions?
  // (assuming the entire population has the maximum fitness)
  SetUpdateTime(log((GetDesiredDivisions()+(double)GetTotalPopSize()) / GetTotalPopSize()) / (GetMaxW() * log(2)));

  //What is the minimum time required to get a single division?
  SetTimeToNextWholeCell(0);
  DetermineDivisionTime();

  // At a minumum, we want to make sure that one cell division took place
  if (GetTimeToNextWholeCell() > GetUpdateTime()) 
  {
     if (GetVerbose())std::cout << "Time to next whole cell greater than update time: " << GetTimeToNextWholeCell() << " < " << GetUpdateTime() <<std::endl;    
     SetUpdateTime(GetTimeToNextWholeCell());
  }

  if (GetVerbose() == 1) 
  {
     std::cout << "Update time: " << GetUpdateTime() <<std::endl;    
  }
            
  //Now update all lineages by the time that actually passed
 
  SetNewPopSize(0);
  UpdateLineages();

  SetCompletedDivisions(GetNewPopSize() - GetTotalPopSize());
              
  if (GetVerbose())std::cout << "Completed divisions: " << GetCompletedDivisions() <<std::endl;
  SetDivisionsUntilMutation(GetDivisionsUntilMutation() - GetCompletedDivisions());
  SetTotalPopSize(GetNewPopSize());
}

/*@agm The functions below should build a new tree using the tree.h header
       I thought the best way to do this was to comment out all of the previous code in both this file
       and the associated header so it would be clear where stuff was changed. */

void cPopulation::NewSeedSubpopulation(cLineageTree& newtree, 
																			 unsigned int& node_id) 
{	
	cGenotype r, w;
	tree<cGenotype>::iterator top, red_side, white_side;
	long double starting_fitness = 1.0;
	
	//initialize object of cSubpopulation type
	cSubpopulation red, white;
	node_id = 0;
	
	/*@agm The head is set, though I don't think it is necessary.  The unique code and fitness is set. */
	
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
																	 unsigned int& node_id) 
{
	m_populations.push_back(subpop);
	SetNumberOfSubpopulations(GetNumberOfSubpopulations()+1);
	node_id++;
	//std::cout << subpop.GetNode_id() << " " << subpop.GetFitness() << std::endl;
}

void cPopulation::NewMutate(gsl_rng * randgen, 
														cLineageTree& newtree, 
														unsigned int& node_id) 
{	
	SetTotalMutations(GetTotalMutations()+1);
	
	if (m_verbose) std::cout << "* Mutating!" << std::endl;
	
	//Mutation happened in the one that just divided
	
	//Break ties randomly here.
	cSubpopulation& ancestor = m_populations[m_divided_lineages[rand() % m_divided_lineages.size()]];          
	
	cSubpopulation new_subpop;
	
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


