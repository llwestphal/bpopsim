#ifndef lineageTree_h
#define lineageTree_h

//#include "cPopulation.h"

class cPopulation;
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip> 

class lineageTree
{
   private:
     std::vector<long double> m_tree;
     std::vector<int> m_pointers;
     std::vector<long double> m_frequencies;
   
   public:
     lineageTree() {;};
     
   virtual ~lineageTree() {;};

   void AddNode(long double in_mutation){ m_tree.push_back(in_mutation);}
   void SetPointer(int in_previous){m_pointers.push_back(in_previous);}
   void ClearRuns();
   void PrintTree(const std::string& output_file_name);
   void CalculateFrequencies(cPopulation& in, const std::string& output_file_name);     
      

   std::vector<long double> GetRuns() { return m_tree;}
   std::vector<int> GetPointers() { return m_pointers;}
   int GetSizeOfTree() {return m_tree.size();}
   
   
}; 
#endif
