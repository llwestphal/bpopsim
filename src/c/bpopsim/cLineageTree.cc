#include "cLineageTree.h"
#include "cPopulation.h"

void cLineageTree::PrintTree(const std::string& output_file_name)
{
   //Print everything out
   std::ofstream output_file;
   output_file.open(output_file_name.c_str(),std::ios_base::app);
   output_file << "List of Mutations for Each Lineage (Each line is a lineage) " << std::endl;     
   int current_position;    

   //Start at end of array. The first position was just a dummy node.
   for (int on_branch=GetSizeOfTree()-1; on_branch > 2; on_branch--) 
   {
     current_position = on_branch;
      //Iterate until we trace the lineage back to the ancestor, who points at zero
      while (m_tree[current_position]!=0)
      {  
         //Put a tab between mutations for a particular lineage
         output_file << m_tree[current_position] << "	";
         current_position = m_lineages[current_position];            
      }
      //Put a return key after each lineage
      output_file << std::endl;
   }
   output_file.close();
}

void cLineageTree::CalculateFrequencies(cPopulation& in, const std::string& output_file_name)
{

   std::ofstream output_file;

   //Append to the file
   output_file.open(output_file_name.c_str(),std::ios_base::app);
   output_file << "pop-size	size-of-last-mutation " << std::endl;

   long double* frequency = NULL;
   int size = in.GetNumberOfSubpopulations()+1;
   frequency = new long double[size];  
  
   long double sum = 0;


   std::cout << std::setprecision(10);

	 /*
   //Iterate through the sub populations
   for (std::vector<cSubpopulation>::iterator it = in.GetPopulation().end(); it!=in.GetPopulation().begin(); --it)
   {
      if(m_tree[it->GetLineage()]!=0)
      {
         //Figure out the ancestor
         //std::cout << it->GetLineage() << std::endl;
         //In the frequency array, add the frequency of a particular sublineage to its position in the array
         frequency[it->GetLineage()]+=(it->GetNumber())/(in.GetTotalPopSize());
         //Write the total population size of that lineage, and the size of the mutation it took to get there
         output_file << it->GetNumber() << "	" << it->GetFitness() << std::endl;
         sum+= (it->GetNumber())/(in.GetTotalPopSize());
      }
   }*/
   output_file << std::endl;
   //Verify that the sum of all the frequencies of all the lineages is 1
   std::cout << sum << std::endl;
   delete [] frequency;
   frequency = NULL;


}       
     
/*void cLineageTree::ClearRuns()
{
   m_tree.clear();
   m_lineages.clear();
}*/

void cLineageTree::ClearRuns()
{
	this -> clear();
	//m_lineages.clear();
}
  
