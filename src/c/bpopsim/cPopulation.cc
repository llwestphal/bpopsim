#include "cPopulation.h"
/* */

using namespace std;


/* Copy constructor */

//cPopulation::cPopulation(const cPopulation& in)
//{
//  	m_total_pop_size = in.m_total_pop_size;
//	m_total_mutations = in.m_total_mutations;
//	m_total_subpopulations_lost = in.m_total_subpopulations;
 // 	m_transfers = in.m_transfers;
//}



void cPopulation::AddSubpopulation(cSubpopulation& subpop)
{

	//GetSubpopulations().push_back(subpop);
	m_populations.push_back(subpop);	

	SetNumberOfSubpopulations(GetNumberOfSubpopulations()+1);
}

	
