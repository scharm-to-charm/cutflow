#ifndef CTAG_CALIBRATION_HH
#define CTAG_CALIBRATION_HH

#include <string> 
#include <map>
#include "ctag_defs.hh"

class BaselineJet; 

namespace Analysis { 
  class CalibrationDataInterfaceROOT; 
  class CalibrationDataVariables; 
}

// input structure for both cuts and SF
struct JetTagFactorInputs { 
  double pt; 			// in MeV
  double eta; 			
  double anti_b; 		// log(pc/pb)
  double anti_u; 		// log(pc/pu)
  // these are the flavor_truth_label values, translated to enums
  ctag::Flavor flavor; 		// B, C, U, T, or DATA
};

// translator from the int value in D3PDs to enums used here
ctag::Flavor get_flavor(int flavor_truth_label); 


// the output structure. 'up' and 'down' variations may be symmetric in 
// many cases, but we keep both for flexibility. 
struct JetTagSF { 
  JetTagSF(const std::pair<double, double>&); 
  JetTagSF(); 
  double nominal; 
  double up; 
  double down; 
}; 
// Inevatably we'll multiply the SF
JetTagSF operator*(const JetTagSF&, const JetTagSF&); 
JetTagSF& operator*=(JetTagSF&, const JetTagSF&); 


// The class responsible for the actual calibration 
class CtagCalibration
{
public: 
  CtagCalibration(std::string calibration_file, 
		  std::string env_file = "/tmp/btag.env"); 
  ~CtagCalibration(); 
  JetTagSF scale_factor(const JetTagFactorInputs& jet_tf_inputs) const; 
  bool pass(const JetTagFactorInputs& jet_tf_inputs) const; 
private: 
  // don't allow copying
  CtagCalibration(const CtagCalibration&){} 
  CtagCalibration& operator=(const CtagCalibration&); 

  typedef std::pair<double, double> CalResult; 
  void check_cdi() const; 
  Analysis::CalibrationDataVariables get_vars(double pt, double eta) const; 
  std::string get_label(ctag::Flavor) const; 
  void set_indices(ctag::Flavor); 

  Analysis::CalibrationDataInterfaceROOT* m_cdi; 
  std::string m_op_string; 
  std::string m_jet_author; 
  double m_anti_u_cut; 
  double m_anti_b_cut; 
  std::map<ctag::Flavor, unsigned> m_flav_op_eff_index; 
  std::map<ctag::Flavor, unsigned> m_flav_op_sf_index; 
}; 


#endif 
