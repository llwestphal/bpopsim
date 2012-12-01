#include "common.h"
#include "cSubpopulation.h"
#include "cPopulation.h"

using namespace std;

namespace bpopsim {

cGenotype::cGenotype(
          SimulationParameters& in_simulation_parameters, 
          uint32_t in_node_id
          ) 
: unique_node_id(in_node_id)
, name("")
, fitness(1.0)
, this_mutation_fitness_effect(0.0)
, total_mutation_count(0)
, marked_for_deletion(false)
{ 
  // set mutation counts vector to be as large as the number of categories
  mutation_counts.resize(in_simulation_parameters.mutation_rates_per_division.size(), 0); 
};
  
cGenotype::cGenotype(
                     SimulationParameters& in_simulation_parameters, 
                     uint32_t in_node_id,
                     cGenotype& in_ancestor
                     ) 
: unique_node_id(in_node_id)
, name("")
, fitness(in_ancestor.fitness)
, this_mutation_fitness_effect(0.0)
, total_mutation_count(in_ancestor.total_mutation_count)
, mutation_counts(in_ancestor.mutation_counts)
, marked_for_deletion(false)
{ 
  // set mutation counts vector to be as large as the number of categories
  mutation_counts.resize(in_simulation_parameters.mutation_rates_per_division.size(), 0); 
};
  
void cGenotype::AddOneMutation(
                               SimulationParameters& simulation_parameters,
                               gsl_rng * rng
                               )
{
  // @JEB: more of this could be moved to SimulationParameters object
  uint32_t this_mutation_category = simulation_parameters.DetermineMutationCategory(rng);
  double average_mutation_fitness_effect = simulation_parameters.mutation_fitness_effects[this_mutation_category];
  
  if(simulation_parameters.mutation_fitness_effect_model=="e")
  {
    this->this_mutation_fitness_effect = gsl_ran_exponential(rng,average_mutation_fitness_effect);	
  }
    
  if (simulation_parameters.mutation_fitness_effect_model=="u")
  {
    this->this_mutation_fitness_effect = average_mutation_fitness_effect;
  }
  
  // Update our information
  this->total_mutation_count++;
  this->mutation_counts[this_mutation_category]++;
  this->fitness += this->this_mutation_fitness_effect;
}
  
/* Copy constructor */
cSubpopulation::cSubpopulation(const cSubpopulation& in) :
  m_number(0),
  m_marker(0),
  m_genotype(0) { 
  m_number = in.m_number;
  m_marker = in.m_marker;
  m_genotype = in.m_genotype;
};

void cSubpopulation::CreateDescendant(  gsl_rng * randomgenerator, 
                                        cSubpopulation &ancestor, 
                                        SimulationParameters& simulation_parameters,
                                        tree<cGenotype>& in_tree, 
                                        uint32_t node_id
                                      ) 
{	
  if (g_verbose) cout << "Creating descendant with node_id: " << node_id << endl;
  
	tree<cGenotype>::iterator_base new_geno_it;
	// There is only one new one...
	SetNumber(1.0);

	// ...taken from the ancestor.
	ancestor.SetNumber(ancestor.GetNumber()-1.0);
	  
	SetMarker(ancestor.GetMarker());
	
	cGenotype new_genotype(simulation_parameters, node_id, ancestor.GetGenotype());
  new_genotype.AddOneMutation(simulation_parameters, randomgenerator);
  new_geno_it = in_tree.append_child(ancestor.m_genotype, new_genotype); 
	SetGenotype(new_geno_it);
	
}

void cSubpopulation::AddToTree(tree<cGenotype> &in_tree, 
                               tree<cGenotype>::iterator parent, 
                               cGenotype child) {
  in_tree.append_child(parent, child);
}

void cSubpopulation::Transfer(double success_prob, 
                              gsl_rng * randomgenerator) {
		
	uint32_t random_gsl_int = gsl_ran_binomial(randomgenerator, 
                                             success_prob, 
                                             GetIntegralNumber());
   
	m_number = random_gsl_int;
}

}
