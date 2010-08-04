#ifndef cPopulation_h
#define cPopulation_h

#include "cSubpopulation.h"

#include <vector>

class cPopulation {

private:
	long double m_total_pop_size;
	int m_total_mutations;
	int m_total_subpopulations_lost;
  	int m_transfers;
	std::vector<cSubpopulation> m_populations;
	double m_divisions_until_mutation;
	double m_desired_divisions;
	int m_number_of_subpopulations;
	long double m_new_pop_size;
	double m_completed_divisions;
	double m_update_time;
public:

	cPopulation()
	{
		m_total_pop_size = 2;
		m_total_mutations = 0;
		m_total_subpopulations_lost = 0;
		m_divisions_until_mutation = 0;
		m_desired_divisions = 0;
  		m_transfers = 0;
	}
	//cPopulation(const cPopulation& in); // Copy constructor
	virtual ~cPopulation() { ; };
	//cPopulation(const cPopulation& in);
	
	const long double GetTotalPopSize() { return m_total_pop_size; }
	const int GetTotalMutations() { return m_total_mutations; }
	const int GetTotalSubpopulationsLost() { return m_total_subpopulations_lost; }
	const int GetTransfers() { return m_transfers; }
	std::vector<cSubpopulation> GetSubpopulations() { return m_populations; }
	const double GetDivisionsUntilMutation() { return m_divisions_until_mutation; }
	const double GetDesiredDivisions() { return m_desired_divisions; }
	const int GetNumberOfSubpopulations() { return m_number_of_subpopulations; }
	const int GetNewPopSize() { return m_new_pop_size; }
	const double GetCompletedDivisions() { return m_completed_divisions; }
	const double GetUpdateTime() { return m_update_time; }



//  virtual MutationList& GetMutations(const char in_marker) {};

	void SetTotalPopSize(long double in_total_pop_size) { m_total_pop_size = in_total_pop_size; }
	void SetTotalMutations(int in_total_mutations) { m_total_mutations = in_total_mutations; }
	void SetTotalSubpopulationsLost(int in_total_subpopulations_lost) { m_total_subpopulations_lost=in_total_subpopulations_lost; }
	void SetTransfers(int in_transfers) { m_transfers = in_transfers; }
	void SetDivisionsUntilMutation(double in_divisions_until_mutation){ m_divisions_until_mutation = in_divisions_until_mutation; }	
	void SetDesiredDivisions(double in_desired_divisions){ m_desired_divisions = in_desired_divisions; }
	void SetNumberOfSubpopulations(int in_number_of_subpopulations){ m_number_of_subpopulations = in_number_of_subpopulations; }
	void SetNewPopSize(int in_new_pop_size) {m_new_pop_size = in_new_pop_size; }
	void SetCompletedDivisions(int in_completed_divisions) {m_completed_divisions = in_completed_divisions; }
	void SetUpdateTime(double in_update_time) {m_update_time = in_update_time; }

	void AddSubpopulation(cSubpopulation& subpop);
};


#endif
