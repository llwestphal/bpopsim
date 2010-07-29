#ifndef cPopulation_h
#define cPopulation_h

#include "cSubpopulation.h"

class cPopulation {

private:
	long double m_total_pop_size;
	int m_total_mutations;
	int m_total_subpopulations_lost;
  	int m_transfers;
	
public:
	cPopulation();
  
	virtual ~Population() { ; }; 
  
	const long double GetTotalPopSize() { return m_total_pop_size; }
	const int GetTotalMutations() { return m_total_mutations; }
	const int GetTotalSubpopulationsLost() { return m_total_subpopulations_lost; }
	const int GetTransfers() { return m_transfers; }
//  virtual MutationList& GetMutations(const char in_marker) {};

	void SetTotalPopSize() { m_total_pop_size = in_total_pop_size; }
	void SetTotalMutations() { m_total_mutations = in_total_mutations; }
	void SetTotalSubpopulationsLost() { m_total_subpopulations_lost=in_total_subpopulations_lost; }
	void SetTransfers() { m_transfers = in_transfers; }


	virtual cPopulation& CreatePopulation(cSubpopulation& subpop);
};


#endif
