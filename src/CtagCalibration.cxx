#include "CtagCalibration.hh"
#include "CalibrationDataInterface/CalibrationDataInterfaceROOT.h"
#include "TEnv.h"
#include <stdexcept> 
#include <cassert> 
#include <fstream>

// translator from flavor truth label to enum
ctag::Flavor get_flavor(int flavor_truth_label) { 
  switch (flavor_truth_label) { 
  case -1: return ctag::DATA; 
  case 0: return ctag::U; 
  case 4: return ctag::C; 
  case 5: return ctag::B; 
  case 15: return ctag::T; 
  default: throw std::domain_error(
    "got weird flavor truth label in " __FILE__); 
  }
}


// fill from the CalibrationDataInterfaceROOT output
JetTagSF::JetTagSF(const std::pair<double, double>& cal_result): 
  nominal(cal_result.first), 
  up(cal_result.first + cal_result.second), 
  down(cal_result.first - cal_result.second)
{}

// fill with NaN
JetTagSF::JetTagSF(): 
  nominal(0.0/0.0), 
  up(0.0/0.0), 
  down(0.0/0.0)
{}

JetTagSF operator*(const JetTagSF& left, const JetTagSF& right) { 
  JetTagSF out; 
  out.nominal = left.nominal * right.nominal; 
  out.up = left.up * right.up; 
  out.down = left.down * right.down; 
  return out; 
}
JetTagSF& operator*=(JetTagSF& left, const JetTagSF& right) { 
  left.nominal *= right.nominal; 
  left.up *= right.up; 
  left.down *= right.down; 
  return left; 
}

// ---- forward declare some utilities ----
namespace { 
  // shorthand for variations
  double up(const std::pair<double, double>& res) { 
    return res.first + res.second;
  } 
  double down(const std::pair<double, double>& res) { 
    return res.first - res.second;
  } 

  // get function (workaround for lacking map::at() in c++98)
  template<typename T, typename K> 
  const T& get(const std::map<K,T>&, const K&); 

  // check if files exist
  bool exists(std::string file_name); 
}

// ----- c-tag calibration class

CtagCalibration::CtagCalibration(std::string cdi_file, std::string env_file): 
  m_cdi(0),
  m_jet_author("AntiKt4TopoLCJVF")
{
  // create the env file and add the root file to it
  // (this is a terrible hack to get around the werid interface to the cdi)
  TEnv env(env_file.c_str()); 
  env.SetValue("File", cdi_file.c_str()); 
  int error = env.WriteFile(env_file.c_str()); 
  if (error) throw std::runtime_error("couldn't save " + env_file); 
  if (!exists(env_file)) throw std::runtime_error(
    "no file " + env_file + " found"); 
  if (!exists(cdi_file)) throw std::runtime_error(
    "no file " + cdi_file + " found"); 
    
  m_cdi = new Analysis::CalibrationDataInterfaceROOT(
    "JetFitterCOMBCharm", env_file, ""); 

  using namespace ctag; 
  // WARNING: these are hacks until we get a better CDI
  m_op_string = "-1_0_-0_82"; 

  check_cdi(); 

  // note that the cuts used here aren't consistent with the ones
  // listed above. 
  m_anti_u_cut = 0.95; 
  m_anti_b_cut = -0.9;

  // flavors start with B and end with DATA, so we can hack a for loop
  for (int flavor = ctag::B; flavor < ctag::DATA; flavor++) { 
    set_indices(static_cast<ctag::Flavor>(flavor)); 
  }

}

CtagCalibration::~CtagCalibration() { 
  delete m_cdi; 
  m_cdi = 0; 
}

JetTagSF CtagCalibration::scale_factor(
  const JetTagFactorInputs& tf_inputs) const { 

  Analysis::CalibrationDataVariables vars = get_vars(
    tf_inputs.pt, tf_inputs.eta); 
  Analysis::Uncertainty unct = Analysis::Total; 

  unsigned sf_index = get(m_flav_op_sf_index, tf_inputs.flavor); 
  unsigned eff_index = get(m_flav_op_eff_index, tf_inputs.flavor); 

  if (pass(tf_inputs) ) { 
    return JetTagSF(m_cdi->getScaleFactor(vars, sf_index, unct));
  }
  return JetTagSF(m_cdi->getInefficiencyScaleFactor(
		    vars, sf_index, eff_index, unct));
}

bool CtagCalibration::pass(const JetTagFactorInputs& tf_inputs) 
  const { 
  bool pass_u = tf_inputs.anti_u > m_anti_u_cut; 
  bool pass_b = tf_inputs.anti_b > m_anti_b_cut; 
  return pass_b && pass_u; 
}


// ----- private stuff ---------
void CtagCalibration::check_cdi() const { 
  if (! m_cdi->getBinnedScaleFactors(
	m_jet_author, 
	get_label(ctag::B), 
	m_op_string)) { 
    throw std::runtime_error("ctag calibration information not found"); 
  }
}

Analysis::CalibrationDataVariables CtagCalibration::get_vars(double pt, 
							     double eta) 
  const { 
  Analysis::CalibrationDataVariables vars; 
  vars.jetAuthor = m_jet_author; 
  vars.jetPt = pt; 
  vars.jetEta = eta; 
  return vars; 
}

std::string CtagCalibration::get_label(ctag::Flavor flavor) const { 
  switch (flavor) { 
  case ctag::B: return std::string("B"); 
  case ctag::C: return std::string("C"); 
  case ctag::U: return std::string("Light"); 
  case ctag::T: return std::string("T"); 
  case ctag::DATA: throw std::domain_error(
    "shouldn't be asking for a data truth label in " __FILE__); 
  default: 
    assert(false); 
  }
}

void CtagCalibration::set_indices(ctag::Flavor flav) 
{ 
  const std::string& label = get_label(flav); 
  bool ok_sf = m_cdi->retrieveCalibrationIndex(
    label, 
    m_op_string, 
    m_jet_author, 
    true, 			// is sf
    m_flav_op_sf_index[flav]); 

  bool ok_eff = m_cdi->retrieveCalibrationIndex(
    label, 
    m_op_string, 
    m_jet_author, 
    false, 			// is not sf
    m_flav_op_eff_index[flav]); 
  if (!ok_eff || !ok_sf) { 
    std::string problem = "problem setting op " + m_op_string + 
      " " + get_label(flav) + " there may be something wrong with your CDI"; 
    throw std::runtime_error(problem); 
  }
}

// disabled assignment operator 
CtagCalibration& CtagCalibration::operator=(const CtagCalibration&){ 
  assert(false); 
  return *this; 
}

namespace { 
  template<typename T, typename K> 
  const T& get(const std::map<K,T>& map, const K& index) { 
    typename std::map<K, T>::const_iterator pos = map.find(index);
    if (pos == map.end()) throw std::domain_error(
      "asked for unknown index in " __FILE__); 
    return pos->second; 
  }

  bool exists(std::string file_name) { 
    std::ifstream file(file_name.c_str(), std::ios::binary); 
    if (!file) { 
      file.close(); 
      return false; 
    }
    file.close(); 
    return true; 
  }
}

