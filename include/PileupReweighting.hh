#ifndef PILEUPREWEIGHTING_HH
#define PILEUPREWEIGHTING_HH

class SusyBuffer;
class RunInfo;
namespace Root {
  class TPileupReweighting;
}

#include <string>

class PileupReweighting
{
public:
  PileupReweighting(const std::string& pu_config,
		    const std::string& pu_lumicalc);
  ~PileupReweighting();
  float get_pileup_weight(const SusyBuffer&);
  unsigned random_run_number(const SusyBuffer&) const;
private:
  Root::TPileupReweighting* m_prw;
};


#endif
