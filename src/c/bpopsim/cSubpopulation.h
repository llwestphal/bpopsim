#ifndef cSubpopulation_h
#define cSubpopulation_h


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "stdint.h"

//typedef std::vector<int> MutationList;

class cSubpopulation {

private:
	long double m_fitness;
	long double m_number;
	char m_marker;
  
public:
	cSubpopulation();
	cSubpopulation(const cSubpopulation& in); // Copy constructor
	virtual ~cSubpopulation() { ; }; 
	  
	const long double GetFitness() { return m_fitness; }
	const long double GetNumber() { return m_number; }
	const char GetMarker() { return m_marker; }
	//  virtual MutationList& GetMutations(const char in_marker) {};
	 
	void SetFitness(const long double in_fitness) { m_fitness = in_fitness; }
	void SetNumber(const long double in_number) { m_number = in_number; }
	void SetMarker(const char in_marker) { m_marker = in_marker; }
	void Transfer(long double success_prob, gsl_rng * randomgenerator);
	long double MutantFitness(long double in_fitness, double in_average_mutation_s, char in_type_of_mutation, gsl_rng * randomgenerator);  

	virtual cSubpopulation& CreateDescendant(gsl_rng * randomgenerator);


};


#endif
