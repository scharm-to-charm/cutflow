#include "CutCounter.hh"

double& CutCounter::operator[](std::string key)
{
  bool new_cut = (m_counts.find(key) == m_counts.end());
  if (new_cut) {
    m_cuts.push_back(key);
  }
  return m_counts[key];
}

std::vector< std::pair<std::string, double> > CutCounter::get_ordered_cuts()
  const
{
  typedef std::vector<std::string>::const_iterator IdxItr;
  std::vector< std::pair<std::string, double> > ordered_cuts;
  for (IdxItr itr = m_cuts.begin(); itr != m_cuts.end(); itr++) {
    ordered_cuts.push_back(*m_counts.find(*itr));
  }
  return ordered_cuts;
}
