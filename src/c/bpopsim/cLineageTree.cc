#include "cLineageTree.h"
#include "cPopulation.h"

void cLineageTree::PrintTree(const std::string& output_file_name)
{
   //Print everything out
   std::ofstream output_file;
   output_file.open(output_file_name.c_str(),std::ios_base::app);
   output_file << "tree ";     
   int current_position;    
   for (int on_branch=GetSizeOfTree()-1; on_branch > 0; on_branch--) 
   {
      current_position = on_branch;
      while (m_tree[current_position]!=0)
      {
         output_file << m_tree[current_position] << " ";
         current_position = m_pointers[current_position];            
      }
      output_file << std::endl;
   }
   output_file.close();
}

void cLineageTree::CalculateFrequencies(cPopulation& in, const std::string& output_file_name)
{

   std::ofstream output_file;
   output_file.open(output_file_name.c_str(),std::ios_base::app);
   output_file << "freq ";

   long double* frequency = NULL;
   int size = in.GetNumberOfSubpopulations()+1;
   frequency = new long double[size];  
  
   long double sum = 0;


   std::cout << std::setprecision(10);


   for (std::vector<cSubpopulation>::iterator it = in.GetPopulation().end(); it!=in.GetPopulation().begin(); --it)
   {
      if(m_tree[it->GetPointer()]!=0)
      {
         std::cout << it->GetPointer() << std::endl;
         std::cout << size << std::endl;
         frequency[it->GetPointer()]+=(it->GetNumber())/(in.GetTotalPopSize());
         output_file << it->GetNumber() << " " << it->GetFitness() << std::endl;
         sum+= (it->GetNumber())/(in.GetTotalPopSize());
      }
   }
   std::cout << sum << std::endl;
   delete [] frequency;
   frequency = NULL;


}       
     
void cLineageTree::ClearRuns()
{
   m_tree.clear();
   m_pointers.clear();
}
  
