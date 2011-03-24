#ifndef cSubpopulation_h
#define cSubpopulation_h

//I added this macros section to avoid scoping for common function calls
//Anything added here will be available throughout the program
#define Endl std::endl
#define Cout std::cout


#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <iomanip> 
#include <gsl/gsl_randist.h>
#include <stdint.h>
#include <time.h>
#include "tree.hh"
#include "tree_util.hh"

#include <boost/program_options.hpp>


/*@agm Rather than create a new class, I created a struct... it seemed simpler
       since the cGenotype type will not have any methods either way. */

struct cGenotype{
 uint64_t unique_node_id;
 long double fitness;
 };

//@agm I put this in class form as well, just in case we need to add methods one day.
//     For now, the struct works better.

/*class cGenotype{
  public:
    uint64 unique_node_id;
    long double fitness;
 };*/

struct cGenotypeFrequency{
	uint64_t unique_node_id;
	long double frequency;
};

class cSubpopulation {

private:
  long double m_number;
  char m_marker;
  tree<cGenotype>::iterator m_genotype;

public:
  cSubpopulation();
  cSubpopulation(const cSubpopulation& in); // Copy constructor
  virtual ~cSubpopulation() { ; }; 

  const long double GetFitness() { return (*m_genotype).fitness; }
  tree<cGenotype>::iterator GetGenotypeIter() { return m_genotype; }
  const uint64_t GetNode_id() { return (*m_genotype).unique_node_id; }
  const long double GetNumber() { return m_number; }
  const char GetMarker() { return m_marker; }
  
  void SetGenotype(tree<cGenotype>::iterator location) { m_genotype = location; }
  void SetNumber(const long double in_number) { m_number = in_number; }
  void SetMarker(const char in_marker) { m_marker = in_marker; } 
  void Transfer(long double success_prob, gsl_rng * randomgenerator);
  long double MutantFitness(long double in_fitness, 
                            double in_average_mutation_s, 
                            char in_type_of_mutation, 
                            gsl_rng * randomgenerator);  
  virtual void NewCreateDescendant(gsl_rng * randomgenerator, 
                                   cSubpopulation &ancestor, 
                                   double averageselectioncoefficient, 
                                   char beneficialmutationdistribution, 
                                   tree<cGenotype> in_tree, 
                                   uint64_t node_id);

};
#endif
