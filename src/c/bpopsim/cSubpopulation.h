#ifndef cSubpopulation_h
#define cSubpopulation_h

#include <vector>

typedef std::vector<int> MutationList;

class cSubpopulation {

private:
	long double m_fitness;
	long double m_number;
	char m_marker;
  
public:
  cSubpopulation();
  cSubpopulation(const cSubpopulation& in); // Copy constructor
  virtual ~cSubpopulation() { ; }; 
  
  const long double GetFitness() { return m_fitness; }
  const long double GetNumber() { return m_number; }
  const char GetMarker() { return m_marker; }
//  virtual MutationList& GetMutations(const char in_marker) {};
 
  void SetFitness(const long double in_fitness) { m_fitness = in_fitness; }
  void SetNumber(const long double in_number) { m_number = in_number; }
  void SetMarker(const char in_marker) { m_marker = in_marker; }
  
  virtual cSubpopulation& CreateDescendant(long double in_fitness);
};


#endif
