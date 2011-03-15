#ifndef cLineageTree_h
#define cLineageTree_h
typedef long double cGenotype;

//#include "cPopulation.h"

class cPopulation;
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip> 
#include "tree.hh"

class cLineageTree : public tree<cGenotype> {
	private:
     std::vector<long double> m_tree;
     std::vector<int> m_lineages;
     std::vector<long double> m_frequencies;
   
   public:
     cLineageTree() {;};
     
   virtual ~cLineageTree() {;};
	
	//tree<long double>::iterator BeginTree(){ return m_newtree.begin(); }
	//void AddChild(tree<long double>::iterator parent, long double child){ m_newtree.append_child(parent, child); }
	//void AddSibling(tree<long double>::iterator old_parent, long double sibling){ m_newtree.insert(old_parent, sibling); }
	
   void AddNode(long double in_mutation){ m_tree.push_back(in_mutation);}
   void SetLineage(int in_previous){m_lineages.push_back(in_previous);}
   void ClearRuns();
   void PrintTree(const std::string& output_file_name);
   void CalculateFrequencies(cPopulation& in, const std::string& output_file_name);     
      

   std::vector<long double> GetRuns() { return m_tree;}
   std::vector<int> GetLineages() { return m_lineages;}
   int GetSizeOfTree() {return m_tree.size();}
   
   
}; 
#endif
