#ifndef cSubpopulation_h
#define cSubpopulation_h

#include <gsl/gsl_randist.h>
#include "stdint.h"

class cSubpopulation {

private:
  long double m_fitness;
  long double m_number;
  char m_marker;
  int m_pointer;
  
public:
  cSubpopulation();
  cSubpopulation(const cSubpopulation& in); // Copy constructor
  virtual ~cSubpopulation() { ; }; 

  const long double GetFitness() { return m_fitness; }
  const long double GetNumber() { return m_number; }
  const char GetMarker() { return m_marker; }
  const int GetPointer() { return m_pointer; }
  //virtual MutationList& GetMutations(const char in_marker) {};
  	 
  void SetPointer(const int in_pointer) {m_pointer = in_pointer; }
  void SetFitness(const long double in_fitness) { m_fitness = in_fitness; }
  void SetNumber(const long double in_number) { m_number = in_number; }
  void SetMarker(const char in_marker) { m_marker = in_marker; } 
  void Transfer(long double success_prob, gsl_rng * randomgenerator);
  long double MutantFitness(long double in_fitness, double in_average_mutation_s, char in_type_of_mutation, gsl_rng * randomgenerator);  

  virtual void CreateDescendant(gsl_rng * randomgenerator, cSubpopulation& ancestor, double averageselectioncoefficient, char distributiontype);

};
#endif
