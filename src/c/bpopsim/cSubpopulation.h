#ifndef cSubpopulation_h
#define cSubpopulation_h

#include <gsl/gsl_randist.h>
#include <stdint.h>
#include "cLineageTree.h"

class cSubpopulation {

private:
  //long double m_fitness;
  long double m_number;
  char m_marker;
  //int m_lineage;
	tree<cGenotype>::tree::iterator m_genotype;

public:
  cSubpopulation();
  cSubpopulation(const cSubpopulation& in); // Copy constructor
  virtual ~cSubpopulation() { ; }; 

  const long double GetFitness() { return *m_genotype; }
  const long double GetNumber() { return m_number; }
  const char GetMarker() { return m_marker; }
  //const int GetLineage() { return m_lineage; }
  //virtual MutationList& GetMutations(const char in_marker) {};
  	 
  //void SetLineage(const int in_lineage) {m_lineage = in_lineage; }
  //void SetFitness(const long double in_fitness) { m_fitness = in_fitness; }
	void SetGenotype(tree<cGenotype>::tree::iterator location) { m_genotype = location; }
  void SetNumber(const long double in_number) { m_number = in_number; }
  void SetMarker(const char in_marker) { m_marker = in_marker; } 
  void Transfer(long double success_prob, gsl_rng * randomgenerator);
  long double MutantFitness(long double in_fitness, double in_average_mutation_s, char in_type_of_mutation, gsl_rng * randomgenerator);  

  //virtual void CreateDescendant(gsl_rng * randomgenerator, cSubpopulation& ancestor, double averageselectioncoefficient, char distributiontype, int sizeoftree);
	virtual void NewCreateDescendant(gsl_rng * randomgenerator, cSubpopulation &ancestor, double averageselectioncoefficient, char beneficialmutationdistribution, tree<cGenotype> in_tree);

};
#endif
