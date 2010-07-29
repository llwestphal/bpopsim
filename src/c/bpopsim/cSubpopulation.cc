#include "cSubpopulation.h"

/* */
cSubpopulation::cSubpopulation()
{
  m_fitness = 0;
  m_number = 0;
  m_marker = 0;
}

/* Copy constructor */
cSubpopulation::cSubpopulation(const cSubpopulation& in)
{
  m_fitness = in.m_fitness;
  m_number = in.m_number;
  m_marker = in.m_marker;
}

cSubpopulation& cSubpopulation::CreateDescendant(long double in_fitness)
{
  //clone ourself
  cSubpopulation& new_sp = *(new cSubpopulation(*this));
  
  //and give new fitness
  new_sp.SetFitness(in_fitness);
  
  //and there is only one new one
  new_sp.SetNumber(1);
  //taken from the ancestor
  SetNumber(GetNumber()-1);

  return new_sp;
}