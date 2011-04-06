#ifndef cSubpopulation_h
#define cSubpopulation_h

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <iomanip> 
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <stdint.h>
#include <time.h>
#include "tree.hh"
#include "tree_util.hh"
#include "icsilog.h"
#include <functional>
#include <list>
#include <utility>

#include <boost/program_options.hpp>


/*@agm Rather than create a new class, I created a struct... it seemed simpler
       since cGenotype really ought to be entirely public with no methods. */

struct cGenotype{
 uint32_t unique_node_id;
 double fitness;
 };

struct cGenotypeFrequency{
	uint32_t unique_node_id;
	double frequency;
};

struct cFrequencySlice {
  cFrequencySlice(uint32_t node_id, double in_low, double in_high) : unique_node_id(node_id), low(in_low), high(in_high) {}
  uint32_t unique_node_id;
  double low;
  double high;
};

struct cSortByLow {
  bool operator () (const cFrequencySlice & lhs , const cFrequencySlice & rhs) const {
    return lhs.low < rhs.low;
  }
};

class cSubpopulation {

private:
  double m_number;
  char m_marker;
  tree<cGenotype>::iterator m_genotype;

public:
  cSubpopulation();
  cSubpopulation(const cSubpopulation& in); // Copy constructor
  virtual ~cSubpopulation() { ; }; 

  const double GetFitness() { return (*m_genotype).fitness; }
  tree<cGenotype>::iterator GetGenotypeIter() { return m_genotype; }
  const uint32_t GetNode_id() { return (*m_genotype).unique_node_id; }
  const double GetNumber() { return m_number; }
  const char GetMarker() { return m_marker; }
  
  void SetGenotype(tree<cGenotype>::iterator location) { m_genotype = location; }
  void SetNumber(const double in_number) { m_number = in_number; }
  void SetMarker(const char in_marker) { m_marker = in_marker; } 
  void Transfer(double success_prob, gsl_rng * randomgenerator);
  double MutantFitness(double in_fitness, 
                            double in_average_mutation_s, 
                            char in_type_of_mutation, 
                            gsl_rng * randomgenerator);  
  virtual void CreateDescendant(gsl_rng * randomgenerator, 
                                   cSubpopulation &ancestor, 
                                   double averageselectioncoefficient, 
                                   char beneficialmutationdistribution, 
                                   tree<cGenotype>& in_tree, 
                                   uint32_t node_id);

};
#endif
