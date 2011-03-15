#include "cSubpopulation.h"

/* */
cSubpopulation::cSubpopulation()
{
   //m_fitness = 0;
   m_number = 0;
   m_marker = 0;
   //m_lineage = 0;

}

/* Copy constructor */
cSubpopulation::cSubpopulation(const cSubpopulation& in)
{
   //m_fitness = in.m_fitness;
   m_number = in.m_number;
   m_marker = in.m_marker;
   //m_lineage = in.m_lineage;

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
/*
void cSubpopulation::CreateDescendant(gsl_rng * randomgenerator, cSubpopulation &ancestor, double averageselectioncoefficient, char beneficialmutationdistribution, int sizeoftree)
{
   // There is only one new one...
   SetNumber(1);
   // ...taken from the ancestor.
   ancestor.SetNumber(ancestor.GetNumber()-1);

   SetMarker(ancestor.GetMarker());
  
   //and give new fitness
   SetFitness(ancestor.MutantFitness(ancestor.GetFitness(), averageselectioncoefficient, beneficialmutationdistribution, randomgenerator));
 
   //Set the subpopulation's internal pointer to itself, i.e. the position of itself in array

   SetLineage(sizeoftree);

}*/

void cSubpopulation::NewCreateDescendant(gsl_rng * randomgenerator, cSubpopulation &ancestor, double averageselectioncoefficient, char beneficialmutationdistribution, tree<cGenotype> in_tree) {
	tree<cGenotype>::iterator new_geno_it;
	// There is only one new one...
	SetNumber(1);
	// ...taken from the ancestor.
	ancestor.SetNumber(ancestor.GetNumber()-1);
	
	SetMarker(ancestor.GetMarker());
  
	//and give new fitness
	//SetFitness(ancestor.MutantFitness(ancestor.GetFitness(), averageselectioncoefficient, beneficialmutationdistribution, randomgenerator));
	
	
	//Need fixing
	
	cGenotype new_genotype;
	new_genotype = ancestor.MutantFitness(ancestor.GetFitness(), averageselectioncoefficient, beneficialmutationdistribution, randomgenerator);
	
	new_geno_it = in_tree.append_child(ancestor.m_genotype, new_genotype);
	
	SetGenotype(new_geno_it);
	
}

void cSubpopulation::Transfer(long double success_prob, gsl_rng * randomgenerator)
{
   SetNumber(gsl_ran_binomial(randomgenerator, success_prob, uint64_t(GetNumber())));
}

