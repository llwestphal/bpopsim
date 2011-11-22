#include "common.h"
#include "cSubpopulation.h"

using namespace std;
using namespace bpopsim;

/* Copy constructor */
cSubpopulation::cSubpopulation(const cSubpopulation& in) :
  m_number(0),
  m_marker(0),
  m_genotype(0) { 
  m_number = in.m_number;
  m_marker = in.m_marker;
  m_genotype = in.m_genotype;
};

double cSubpopulation::MutantFitness(double in_fitness, 
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

void cSubpopulation::CreateDescendant(  gsl_rng * randomgenerator, 
                                        cSubpopulation &ancestor, 
                                        double averageselectioncoefficient, 
                                        char beneficialmutationdistribution, 
                                        tree<cGenotype>& in_tree, 
                                        uint32_t node_id
                                      ) 
{	
  if (g_verbose) cout << "Creating descendant with node_id: " << node_id << endl;
  
	tree<cGenotype>::iterator_base new_geno_it;
	// There is only one new one...
	SetNumber(1);

	// ...taken from the ancestor.
	ancestor.SetNumber(ancestor.GetNumber()-1);
	  
	SetMarker(ancestor.GetMarker());
	
	cGenotype new_genotype;
	new_genotype.fitness = ancestor.MutantFitness(ancestor.GetFitness(), 
                                                averageselectioncoefficient, 
                                                beneficialmutationdistribution, 
                                                randomgenerator);
	new_genotype.unique_node_id = node_id;
  new_genotype.mut_num = ancestor.GetMutNum()+1;
  new_geno_it = in_tree.append_child(ancestor.m_genotype, new_genotype); 
	SetGenotype(new_geno_it);
	
}

void cSubpopulation::AddToTree(tree<bpopsim::cGenotype> &in_tree, 
                               tree<cGenotype>::iterator parent, 
                               bpopsim::cGenotype child) {
  in_tree.append_child(parent, child);
}

void cSubpopulation::Transfer(double success_prob, 
                              gsl_rng * randomgenerator) {
		
	uint32_t random_gsl_int = gsl_ran_binomial(randomgenerator, 
                                             success_prob, 
                                             double(GetNumber()));
	

	//if ( random_gsl_int > 10 ) std::cout << std::endl << "This is gsl out: " << " " << 
		//randomgenerator << " " << success_prob << " " << random_gsl_int << " " << random_gsl_int << std::endl;
   
	SetNumber(random_gsl_int);
}

