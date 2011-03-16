#ifndef cLineageTree_h
#define cLineageTree_h

/*@agm Rather than create a new class, I created a struct... it seemed simpler
       this is not resistance, I'm just not sure what @jeb wanted as a class with static variables. */

typedef struct {
	unsigned int unique_node_id;
	long double fitness;
} cGenotype;

class cPopulation;
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip> 
#include "tree.hh"

/*@agm I chopped out all of the function declarations here and simply told cLineageTree 
       to inherit everything from tree.h with the cGenotype struct as the variable type.*/

class cLineageTree : public tree<cGenotype> {
	private:
     std::vector<long double> m_tree;
     std::vector<int> m_lineages;
     std::vector<long double> m_frequencies;
}; 
#endif
