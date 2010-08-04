#include "cPopulation.h"
/* */

void cPopulation::AddSubpopulation(cSubpopulation& subpop)
{

	m_populations.push_back(subpop);	

	SetNumberOfSubpopulations(GetNumberOfSubpopulations()+1);
}

void cPopulation::UpdateLineages()
{
	for (std::vector<cSubpopulation>::iterator it = m_populations.begin(); it!=m_populations.end(); ++it)
	{
		if (it->GetNumber() == 0) continue;
        	it->SetNumber(it->GetNumber() * exp(log(2) * GetUpdateTime() * it->GetFitness()));

		SetNewPopSize(GetNewPopSize() + it->GetNumber());
		
	}
}

void cPopulation::DetermineDivisionTime()
{
	for (int i=0; i < int(GetNumberOfSubpopulations()); i++) 
	{
		if (m_populations[i].GetNumber() == 0) continue;
       		//what is the time to get to the next whole number of cells?
		SetCurrentCells(m_populations[i].GetNumber());
		SetWholeCells(floor(m_populations[i].GetNumber())+1);
       		// WC = N * exp(growth_rate * t) 
	
      		SetThisTimeToNextWholeCell(log(GetWholeCells() / GetCurrentCells()) / (m_populations[i].GetFitness()));       
        	if ( GetTimeToNextWholeCell() == 0 || (GetThisTimeToNextWholeCell() < GetTimeToNextWholeCell()) ) 
		{
			m_divided_lineages.clear();
          		SetTimeToNextWholeCell(GetThisTimeToNextWholeCell());
       			m_divided_lineages.push_back(i); //a list, because there can be ties
			
   		}
          	else if (GetThisTimeToNextWholeCell() == GetTimeToNextWholeCell())
          	{
          		m_divided_lineages.push_back(i); //a list, because there can be ties
			
          	}
        }
	
}

void cPopulation::Mutate(gsl_rng * randgen)
{
	if (GetDivisionsUntilMutation() <= 0) 
	{
		SetTotalMutations(GetTotalMutations()+1);
		if (m_verbose) std::cout << "* Mutating!" << std::endl;
        
		//Mutation happened in the one that just divided

		//Break ties randomly here.
		cSubpopulation& ancestor = m_populations[m_divided_lineages[rand() % m_divided_lineages.size()]];        	
			          	
		cSubpopulation new_lineage = ancestor.CreateDescendant(randgen);

		if (GetVerbose()) std::cout << "  Color: " << new_lineage.GetMarker() << std::endl;
		if (GetVerbose()) std::cout << "  New Fitness: " << new_lineage.GetFitness() << std::endl;
		AddSubpopulation(new_lineage);
		//Update maximum fitness
		if(new_lineage.GetFitness() > GetMaxW()) 
		{
			SetMaxW(new_lineage.GetFitness());
		}
        
	}
}

void cPopulation::Resample(gsl_rng * randgen)
{
	//When it is time for a transfer, resample population
	if (GetTotalPopSize() >= GetPopSizeBeforeDilution()) 
	{
        
		if (GetVerbose()) std::cout << ">> Transfer!" << std::endl;
        	//Is there an exists() in C++?

        	m_by_color[RED] = 0;
        	m_by_color[WHITE] = 0;
		
		SetNewPopSize(0);        
           	for (std::vector<cSubpopulation>::iterator it = m_populations.begin(); it!=m_populations.end(); ++it)
		{
			if (it->GetNumber() == 0) continue;
            
               		// Perform accurate binomial sampling only if below a certain population size
            		if (it->GetNumber() < GetBinomialSamplingThreshold()) 
			{
              			if (GetVerbose()) std::cout << "binomial" << it->GetNumber() << std::endl;
				
            			it->Transfer(GetTransferBinomialSamplingP(), randgen);
            		}
            		// Otherwise, treat as deterministic and take expectation...
            		else 
			{
			       	it->SetNumber(it->GetNumber() * GetTransferBinomialSamplingP());
			}
            
            		// Keep track of lineages we lost
            		if (it->GetNumber() == 0) 
			{
			       	SetTotalSubpopulationsLost(GetTotalSubpopulationsLost()+1);
			}
            

            		SetNewPopSize(GetNewPopSize() + floor(it->GetNumber()));
			if (GetVerbose()) std::cout << it->GetMarker() << std::endl;
			if (GetVerbose()) std::cout << it->GetNumber() << std::endl;
			if (GetVerbose()) std::cout << it->GetFitness() << std::endl;
			if (it->GetMarker() == 'r') 
			{
		            	m_by_color[RED] += it->GetNumber();
		        }
 			else 
			{
				m_by_color[WHITE] += it->GetNumber();
		        }
		}
		//One color was lost, bail	
	        if ( (m_by_color[RED] == 0) || (m_by_color[WHITE] == 0) ) 
		{
			SetKeepTransferring(false);
		}
            
       		if (GetVerbose()) std::cout << "Colors: " << m_by_color[RED] << " / " << m_by_color[WHITE] << std::endl;
       		SetRatio(m_by_color[RED] / m_by_color[WHITE]);
          
       		SetTransfers(GetTransfers()+1);

       		if ( (GetTransfers() >= 0) && (GetTransfers() % GetTransferIntervalToPrint() == 0) ) 
		{
	        	m_this_run.push_back(GetRatio());
     	    		if (GetVerbose() == 1) 
			{
		                std::cout << "Transfer " << GetTransfers() << " : " << GetTotalPopSize() << 
				"=>" << GetNewPopSize() << "  R/W Ratio: " << GetRatio() << std::endl;	
		                std::cout << "Total mutations: " << GetTotalMutations() << " Maximum Fitness: " << GetMaxW() << std::endl;
		                std::cout << "Size = " << m_this_run.size() << std::endl;
		        }
       	  	}
          
  		SetTotalPopSize(GetNewPopSize());

  		if ( (int(m_this_run.size()) >= GetMinimumPrinted()) && ((GetRatio() > GetMaxDivergenceFactor()) || 									(GetRatio() < 1/GetMaxDivergenceFactor())) ) 	
		{
			if (GetVerbose()) std::cout << "DIVERGENCE CONDITION MET" << std::endl;
			SetKeepTransferring(false);
		}	
	}
}

void cPopulation::PushBackRuns()
{
	m_runs.push_back(m_this_run);
}
void cPopulation::PrintOut()
{

	//Print everything out
	std::ofstream output_file;
	output_file.open ("output.txt");
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

void cPopulation::ClearRuns()
{
	m_this_run.clear();
	m_populations.clear();
}



