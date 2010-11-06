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
    tree<cSubpopulation*>::iterator m_start_node, m_first_branch, m_current_branch;
    tree<cSubpopulation*> m_base;  
 

  public:
      lineageTree()
      {
         tree<cSubpopulation*>::iterator m_start_node = m_base.begin();

      }


      virtual ~lineageTree() {};

  void EstablishRoots(cSubpopulation* in_pop){ m_first_branch = m_base.insert(m_start_node, in_pop); }
  
  void SetFirstBranch() {m_current_branch = m_first_branch;}
  
  void AppendProgeny(cSubpopulation* in_pop, tree<cSubpopulation*>::iterator& in_it) {in_it = m_base.append_child(m_current_branch,in_pop);}
 
  void SetCurrentBranch(tree<cSubpopulation*>::iterator in_current_branch) {m_current_branch = in_current_branch; }

  tree<cSubpopulation*>::iterator GetFirstBranch() {return m_first_branch; }
  tree<cSubpopulation*>::iterator GetStartNode() {return m_start_node; }
  tree<cSubpopulation*> GetBase() {return m_base; }
  tree<cSubpopulation*>::iterator GetCurrentBranch() {return m_current_branch; }
}; 

#endif
