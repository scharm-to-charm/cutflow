//////////////////////// -*-c++-*-//////////////////////////
// This class has been automatically generated on
// Wed Jun 20 17:06:57 2012 by ROOT version 5.30/00
// from TTree susy/susy
// found on file: file1.root
//////////////////////////////////////////////////////////

#ifndef susy_h
#define susy_h

#include <vector> 
#include <string> 


class SmartChain; 

class SusyBuffer {
public :
  unsigned int RunNumber; 
  unsigned int EventNumber; 
  unsigned lbn; 
  bool          xe80_tclcw_tight; 
  bool          xe80T_tclcw_loose; 
  bool          xe80_tclcw_loose; 
  float         averageIntPerXing;
  unsigned        larError;
  unsigned        tileError;
  unsigned        coreFlags; 

  //  --- truth info ---
  unsigned int mc_channel_number; 
  int hfor_type; 
  int spart1_pdgid; 
  int spart2_pdgid; 
  // float mc_event_weight; 

  // hack for sherpa bug
  std::vector<std::vector<double> >* mcevt_weight; 

  // needed for boson filter
  int mc_n; 
  std::vector<float>* mc_pt; 
  std::vector<float>* mc_eta; 
  std::vector<float>* mc_phi; 
  std::vector<float>* mc_m; 
  std::vector<int>* mc_status; 
  std::vector<int>* mc_pdgId; 

  // for overlap removal in stop samples
  float MET_Truth_NonInt_etx; 
  float MET_Truth_NonInt_ety; 

  // --- Random Event ---
  float Eventshape_rhoKt4LC; 

  //  --- MET things ---
  std::vector<std::vector<float> >*    jet_MET_Egamma10NoTau_wet;
  std::vector<std::vector<float> >*    jet_MET_Egamma10NoTau_wpx;
  std::vector<std::vector<float> >*    jet_MET_Egamma10NoTau_wpy;
  std::vector<std::vector<unsigned> >* jet_MET_Egamma10NoTau_statusWord;

  std::vector<std::vector<float> >*    el_MET_Egamma10NoTau_wet;
  std::vector<std::vector<float> >*    el_MET_Egamma10NoTau_wpx;
  std::vector<std::vector<float> >*    el_MET_Egamma10NoTau_wpy;
  std::vector<std::vector<unsigned> >* el_MET_Egamma10NoTau_statusWord;

  float MET_Egamma10NoTau_CellOut_etx; //CellOut
  float MET_Egamma10NoTau_CellOut_ety; //CellOut
  float MET_Egamma10NoTau_CellOut_sumet; //CellOut
  float MET_CellOut_Eflow_STVF_etx;
  float MET_CellOut_Eflow_STVF_ety; 
  float MET_CellOut_Eflow_STVF_sumet;
  float MET_Egamma10NoTau_RefGamma_etx;
  float MET_Egamma10NoTau_RefGamma_ety;
  float MET_Egamma10NoTau_RefGamma_sumet;

  std::vector<float>   *jet_jvtxf; 
  int             el_n;
  std::vector<float>   *el_eta; 
  std::vector<float>   *el_phi; 
  std::vector<int>     *el_author;
  std::vector<unsigned int> *el_OQ;
  std::vector<int>     *el_mediumPP;
  std::vector<int>     *el_tightPP;
  std::vector<float>     *el_ptcone20;
  // std::vector<float>     *el_trackd0pv;
  // std::vector<float>     *el_trackz0pv;
  std::vector<float>     *el_charge;
  std::vector<float>   *el_cl_E;
  std::vector<float>   *el_cl_eta;
  std::vector<float>   *el_cl_phi;
  std::vector<float>   *el_trackphi;
  std::vector<float>   *el_tracketa;
  std::vector<int>     *el_nPixHits;
  std::vector<int>     *el_nSCTHits;
  int           mu_staco_n;
  std::vector<float>   *mu_staco_pt;
  std::vector<float>   *mu_staco_eta;
  std::vector<float>   *mu_staco_phi;
  std::vector<float>   *mu_staco_ptcone20; 
  std::vector<float>   *mu_staco_charge;
  std::vector<int>     *mu_staco_isCombinedMuon;
  std::vector<int>     *mu_staco_isSegmentTaggedMuon;
  std::vector<int>     *mu_staco_loose;
  std::vector<float>   *mu_staco_id_theta_exPV;
  std::vector<float>   *mu_staco_id_qoverp_exPV;
  std::vector<float>   *mu_staco_me_theta_exPV;
  std::vector<float>   *mu_staco_me_qoverp_exPV;
  std::vector<float>   *mu_staco_ms_phi;
  std::vector<float>   *mu_staco_ms_theta;
  std::vector<float>   *mu_staco_ms_qoverp;
  std::vector<float>   *mu_staco_id_theta;
  std::vector<int>     *mu_staco_nPixHits;
  std::vector<int>     *mu_staco_nSCTHits;
  std::vector<int>     *mu_staco_nTRTHits;
  std::vector<int>     *mu_staco_nPixHoles;
  std::vector<int>     *mu_staco_nSCTHoles;
  std::vector<int>     *mu_staco_nTRTOutliers;
  std::vector<int>     *mu_staco_nPixelDeadSensors;
  std::vector<int>     *mu_staco_nSCTDeadSensors;
  std::vector<float>   *mu_staco_energyLossPar; 

  // muon quality cuts
  std::vector<float>   *mu_staco_qoverp_exPV; 
  std::vector<float>   *mu_staco_cov_qoverp_exPV; 
  std::vector<float>   *mu_staco_d0_exPV; 
  std::vector<float>   *mu_staco_z0_exPV; 

  // Jet Stuff
  int           jet_n;
  std::vector<float>   *jet_pt;
  std::vector<float>   *jet_eta;
  std::vector<float>   *jet_phi;
  std::vector<float>   *jet_E;
  std::vector<float>   *jet_constscale_eta;
  std::vector<float>   *jet_constscale_phi;
  std::vector<float>   *jet_constscale_E;
  std::vector<float>   *jet_constscale_m;
  std::vector<float>   *jet_ActiveAreaPx; 
  std::vector<float>   *jet_ActiveAreaPy; 
  std::vector<float>   *jet_ActiveAreaPz; 
  std::vector<float>   *jet_ActiveAreaE; 
  std::vector<float>   *jet_BCH_CORR_JET; 
  // Event Shape is defined above
  std::vector<float>   *jet_emfrac;
  std::vector<float>   *jet_hecf;
  std::vector<float>   *jet_LArQuality;
  std::vector<float>   *jet_HECQuality;
  std::vector<float>   *jet_AverageLArQF;
  std::vector<float>   *jet_Timing;
  std::vector<float>   *jet_sumPtTrk;
  std::vector<float>   *jet_fracSamplingMax;
  std::vector<int>     *jet_SamplingMax;
  std::vector<float>   *jet_NegativeE;
  std::vector<int>     *jet_flavor_truth_label;
  // mu = averageIntPerXing defined elsewhere 
  std::vector<int>     *vx_nTracks;
  // run number is somewhere else 

  std::vector<float>   *jet_flavor_weight_JetFitterCOMBNN;
  std::vector<float>     *jet_flavor_component_jfitcomb_pu;
  std::vector<float>     *jet_flavor_component_jfitcomb_pb;
  std::vector<float>     *jet_flavor_component_jfitcomb_pc;
  std::vector<float>     *jet_flavor_component_jfitc_pu;
  std::vector<float>     *jet_flavor_component_jfitc_pb;
  std::vector<float>     *jet_flavor_component_jfitc_pc;

  // track branches for single track veto
  std::vector<float>* trk_pt; 
  std::vector<float>* trk_eta; 
  std::vector<float>* trk_phi_wrtPV; 
  std::vector<float>* trk_d0_wrtPV; 
  std::vector<float>* trk_z0_wrtPV; 
  std::vector<int>* trk_ndof; 
  std::vector<float>* trk_chi2; 
  std::vector<int>* trk_nPixHits; 
  std::vector<int>* trk_nSCTHits; 

  SusyBuffer(SmartChain *tree); 
  bool is_data() const; 

private: 
  void set_mc_branches(SmartChain* chain, std::string jc); 
  bool m_is_data; 
};

#endif

