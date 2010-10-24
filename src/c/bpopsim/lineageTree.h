#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "lList.h"

class lineageTree : public cSubpopulation
{
public:
  lineageTree() {

  } 
  lineageTree(cSubpopulation in)
  {
    m_poppointer(in);
  } 

  ~lineageTree() {}

  void SetMutation(const long double in_mutation) { m_mutation = in_mutation; }
  //void SetLineageSize(cSubpopulation in) {*m_lineagesize = *(&in.m_number); }
  
  //long double* GetMutation() { return *m_lineagesize; }
  long double GetLineageSize() { return m_poppointer.GetNumber(); }

	
private:
	long double m_mutation;
	long double m_lineagesize;
	cSubpopulation m_poppointer;
}; 
