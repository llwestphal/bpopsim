#include <iostream>
#include <cmath>

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

long double cSubpopulation::MutantFitness(long double in_fitness, 
																					double in_average_mutation_s, 
																					char in_type_of_mutation, 
																					gsl_rng * randomgenerator) 
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

void cSubpopulation::NewCreateDescendant(gsl_rng * randomgenerator, 
																				 cSubpopulation &ancestor, 
																				 double averageselectioncoefficient, 
																				 char beneficialmutationdistribution, 
																				 tree<cGenotype> in_tree, 
																				 unsigned int node_id) 
{	
	tree<cGenotype>::iterator_base new_geno_it;
	// There is only one new one...
	SetNumber(1);
	
	//**********************************************************************************
	//For some reason uncommenting this will cause everything to get screwed up
	//For example, the GetNumber() function will start returning a negative number
	//This will cause the population number to not properly decrement after dilution
	// etc, etc, etc
	//It's being left for posterity
	//**********************************************************************************

	// ...taken from the ancestor.
	ancestor.SetNumber(ancestor.GetNumber()-1);
	  
	SetMarker(ancestor.GetMarker());
	
	cGenotype new_genotype;
	new_genotype.fitness = ancestor.MutantFitness(ancestor.GetFitness(), 
																								averageselectioncoefficient, 
																								beneficialmutationdistribution, 
																								randomgenerator);
	new_genotype.unique_node_id = node_id;
	
	/* @agm I'm checking to see if the new fitness is the same as the old fitness.
					Ostensibly this ought not be possible, but you never know.
	        Either way, the new genotype will not be added to the list if it has the 
		      same fitness as its parent... at least that's what I think it's doing. */
	
	if (new_genotype.fitness != (*ancestor.m_genotype).fitness) { 
		new_geno_it = in_tree.append_child(ancestor.m_genotype, new_genotype); 
	}

	SetGenotype(new_geno_it);
	
}

void cSubpopulation::Transfer(long double success_prob, 
															gsl_rng * randomgenerator) {
		
	int random_gsl_int = gsl_ran_binomial(randomgenerator, 
													 success_prob, 
													 u_int64_t(GetNumber()));
	

	//if ( random_gsl_int > 10 ) std::cout << std::endl << "This is gsl out: " << " " << 
		//randomgenerator << " " << success_prob << " " << random_gsl_int << " " << random_gsl_int << std::endl;
   
	SetNumber(random_gsl_int);
}

