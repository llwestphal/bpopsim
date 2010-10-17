#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "lList.h"

class lineageTree
{
public:
  lineageTree() {

  }

  ~lineageTree() {}

  void SetMutation(const long double in_mutation) { m_mutation = in_mutation; }
  void SetLineageSize(cSubpopulation& in) { * m_lineagesize = in.GetNumber(); }
  long double* GetMutation() { return *m_lineagesize; }

	
private:
	long double m_mutation;
	long double * m_lineagesize;
}; 
