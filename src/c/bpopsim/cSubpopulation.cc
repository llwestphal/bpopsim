#include "cSubpopulation.h"

/* */
cSubpopulation::cSubpopulation()
{
   m_number = 0;
   m_marker = 0;
}

/* Copy constructor */
cSubpopulation::cSubpopulation(const cSubpopulation& in)
{
   m_number = in.m_number;
   m_marker = in.m_marker;
   m_genotype = in.m_genotype;
}

long double cSubpopulation::MutantFitness(long double in_fitness, double in_average_mutation_s, char in_type_of_mutation, gsl_rng * randomgenerator)
{
   if(in_type_of_mutation=='e')
   {
     return  gsl_ran_exponential(randomgenerator,in_average_mutation_s) + in_fitness;	
   }

   if (in_type_of_mutation=='u')
   {
     return in_fitness + in_average_mutation_s;
   }
  
   return 0;

}

void cSubpopulation::NewCreateDescendant(gsl_rng * randomgenerator, cSubpopulation &ancestor, double averageselectioncoefficient, char beneficialmutationdistribution, tree<cGenotype> in_tree, int node_id) {
	tree<cGenotype>::iterator_base new_geno_it;
	// There is only one new one...
	SetNumber(1);
	// ...taken from the ancestor.
	ancestor.SetNumber(ancestor.GetNumber()-1);
	
	SetMarker(ancestor.GetMarker());
	
	cGenotype new_genotype;
	new_genotype.fitness = ancestor.MutantFitness(ancestor.GetFitness(), averageselectioncoefficient, beneficialmutationdistribution, randomgenerator);
	new_genotype.unique_node_id = node_id;
	
	/* @agm I'm checking to see if the new fitness is the same as the old fitness.
		 Ostensibly this ought not be possible, but you never know.
	   Either way, the new genotype will not be added to the list if it has the 
		 same fitness as its parent... at least that's what I think it's doing. */
	
	if (new_genotype.fitness != (*ancestor.m_genotype).fitness) { 
		new_geno_it = in_tree.append_child(ancestor.m_genotype, new_genotype); 
		//std::cout << new_genotype.unique_node_id << " " << new_genotype.fitness << std::endl;
	
	}
	
	
	
	SetGenotype(new_geno_it);
	
}

void cSubpopulation::Transfer(long double success_prob, gsl_rng * randomgenerator)
{
   SetNumber(gsl_ran_binomial(randomgenerator, success_prob, uint64_t(GetNumber())));
}

