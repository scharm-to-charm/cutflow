#ifndef CUT_COUNTER_HH
#define CUT_COUNTER_HH

#include <map> 
#include <string> 
#include <vector> 

class CutCounter
{
public: 
  double& operator[](std::string key); 
  std::vector< std::pair<std::string, double> > get_ordered_cuts() const; 
private: 
  std::map<std::string, double> m_counts; 
  std::vector<std::string> m_cuts; 
}; 

#endif
