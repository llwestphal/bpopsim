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

void cPopulation::SetParameters()
{

	// Simulation parameters that should be arguments
	SetInitialPopulationSize(2);
	SetPopSizeAfterDilution(int(5E6));             // N sub 0 --int is to get rid of warning
	SetMutationRatePerDivision(1E-8);         // mu
	SetAverageMutationS(0.05);                 // s
  	SetGrowthPhaseGenerations(6.64);


	SetTransferIntervalToPrint(1);
	SetVerbose(1);
	SetTotalTransfers(200000);
	SetMaxDivergenceFactor(100);
	SetReplicates(1000);
	SetMinimumPrinted(8);
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
        	if (GetVerbose())std::cout << "Time to next whole cell greater than update time: " << GetTimeToNextWholeCell() << 										     " < " << GetUpdateTime() <<std::endl;		
        	SetUpdateTime(GetTimeToNextWholeCell());
        }

        if (GetVerbose() == 1) 
	{
		std::cout << "Update time: " << GetUpdateTime() <<std::endl;		
        }
            
        //Now update all lineages by the time that actually passed
	
	SetNewPopSize(0);

	UpdateLineages();
	std::cout << GetNewPopSize() << " " << GetTotalPopSize() <<std::endl;
	SetCompletedDivisions(GetNewPopSize() - GetTotalPopSize());
			        
	if (GetVerbose())std::cout << "Completed divisions: " << GetCompletedDivisions() <<std::endl;
	SetDivisionsUntilMutation(GetDivisionsUntilMutation() - GetCompletedDivisions());
	SetTotalPopSize(GetNewPopSize());
}

void cPopulation::SeedSubpopulations()
{
		//Create red population
		cSubpopulation r;
		//Set parameters
		r.SetNumber(GetInitialPopulationSize()/2);
    		r.SetFitness(1);
    		r.SetMarker('r');
		
		//Add subpopulation to population
		AddSubpopulation(r);
		
		//Seed a white population
		cSubpopulation w;
		w.SetNumber(GetInitialPopulationSize()/2);
    		w.SetFitness(1);
    		w.SetMarker('w');	
    		AddSubpopulation(w);
}
