#ifndef cSubpopulation_h
#define cSubpopulation_h

#include "common.h"

namespace bpopsim {

  struct cGenotype{
    uint16_t unique_node_id;
    string   name;
    double   fitness;
    uint8_t  mut_num;
    double   m_timer;
    bool     m_divided;
    
    cGenotype() : 
      unique_node_id(0),
      name(""),
      fitness(1.0),
      mut_num(0),
      m_timer(1.0),
      m_divided(false) { };
   };

  struct cGenotypeFrequency{
    uint32_t unique_node_id;
    string name;
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
    char   m_marker;
    tree<cGenotype>::iterator m_genotype;

  public:
    cSubpopulation() :   //
    m_number(0),
    m_marker('n') { };
    
    cSubpopulation(const cSubpopulation& in); // Copy constructor
    
    virtual ~cSubpopulation() { ; }; 

    const double GetFitness() { return m_genotype->fitness; }
    const string GetName() { return m_genotype->name; }
    tree<cGenotype>::iterator GetGenotypeIter() { return m_genotype; }
    const uint32_t GetNode_id() { return m_genotype->unique_node_id; }
    const double GetNumber() { return m_number; }
    const char GetMarker() { return m_marker; }
    const uint8_t GetMutNum() { return m_genotype->mut_num; }
    const double GetTimer() { return m_genotype->m_timer; }
    const bool GetDivided() { return m_genotype->m_divided; }
    
    void SetGenotype(tree<cGenotype>::iterator location) { m_genotype = location; }
    void SetNumber(const double in_number) { m_number = in_number; }
    void SetMarker(const char in_marker) { m_marker = in_marker; } 
    void SetTimer(const double in_timer) { m_genotype->m_timer = in_timer; }
    void SetDivided(const double in_divided) { m_genotype->m_divided = in_divided; }
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
    
    virtual void AddToTree(tree<cGenotype>& in_tree,
                           tree<cGenotype>::iterator parent,
                           cGenotype child);
  };
  
}
  
#endif
