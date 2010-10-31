#ifndef lineageTree_h
#define lineageTree_h

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "cSubpopulation.h"
#include "tree.hh"
class lineageTree 
{

  private:
    tree<cSubpopulation*>::iterator m_start_node, m_first_branch;
    tree<cSubpopulation*> m_base;  
 

  public:
      lineageTree()
      {
         m_start_node = m_base.begin();
      }


      virtual ~lineageTree() {};

  void EstablishRoots(cSubpopulation in_pop){ m_first_branch = m_base.insert(m_start_node, &in_pop); }
  
  tree<cSubpopulation*>::iterator GetFirstBranch() {return m_first_branch; }
  tree<cSubpopulation*>::iterator GetStartNode() {return m_start_node; }
  tree<cSubpopulation*> GetBase() {return m_base; }
	
}; 

#endif
