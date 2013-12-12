#include <cstdio>
#include "SmartChain.hh"
#include "SusyBuffer.h"
#include "CutCounter.hh"
#include "ctag_defs.hh"		// for JFC_MEDIUM_* cuts
#include "CtagCalibration.hh"
#include "mctlib.h"

#include "SUSYTools/SUSYObjDef.h"
#include "TLorentzVector.h"
#include "GoodRunsLists/TGRLCollection.h"
#include "GoodRunsLists/TGoodRunsListReader.h"

#include <fstream>
#include <stdexcept> 
#include <cstdlib>
#include <vector> 
#include <algorithm>
#include <iostream> 

#define xxx printf("line %i\n", __LINE__); 

// ============= external files ============
const std::string g_grl_file = "grl.xml"; 
const std::string g_btag_file = "cdi.root"; 
// =========================================

// minimal class to keep track of particles
class IdLorentzVector : public TLorentzVector
{
public: 
  int index; 
  bool pass; 
};

// fill a common object container so it can be passed to each cutflow
struct SelectionObjects 
{
  std::vector<IdLorentzVector> veto_jets; 
  std::vector<IdLorentzVector> veto_electrons; 
  std::vector<IdLorentzVector> veto_muons; 
  std::vector<IdLorentzVector> control_electrons; 
  std::vector<IdLorentzVector> control_muons; 
  std::vector<IdLorentzVector> signal_jets; 
  TVector2 met; 

  TVector2 mu_met; 
  IdLorentzVector electron_jet; 
  bool is_data; 
  bool pass_grl; 

}; 

// multiple object selections are called within the event loop
void signal_selection(const SelectionObjects&, SUSYObjDef* def, 
		      const SusyBuffer& buffer, CutCounter& counter, 
		      double weight = 1.0); 
void el_cr_selection(const SelectionObjects&, SUSYObjDef* def, 
		     const SusyBuffer& buffer, CutCounter& counter); 
void mu_cr_selection(const SelectionObjects&, SUSYObjDef* def, 
		     const SusyBuffer& buffer, CutCounter& counter); 

// these functions check to see if the object passed the SUSYObjDef cuts
std::vector<IdLorentzVector> filter_pass(const std::vector<IdLorentzVector>&); 
std::vector<IdLorentzVector> filter_fail(const std::vector<IdLorentzVector>&); 

// check for c-tags
bool has_medium_tag(int jet_index, const SusyBuffer& buffer); 
// for sorting
bool has_higher_pt(const TLorentzVector& v1, const TLorentzVector& v2); 
// scalar sum for first n jets
double scalar_sum_pt(const std::vector<IdLorentzVector>& obj, size_t num); 
// m_ct function 
double get_m_ct(const IdLorentzVector& v1, const IdLorentzVector& v2); 
// ctag sf function (wrapper for CtagCalibration)
double get_ctag_sf(const IdLorentzVector& jet, const SusyBuffer& buffer, 
		   const CtagCalibration& ctag_cal); 
double get_mctcorr(const TLorentzVector& v1, const TLorentzVector& v2, const TVector2& vmet);

// IO functions
void dump_counts(const CutCounter&, std::string); 
bool exists(std::string file_name); 
std::string red(std::string); 

template<typename M, typename A>
A remove_overlaping(const M& mask, A altered, const float delta_r); 

// =================================================================
// =============== main event loop starts here =====================
// =================================================================

int main (int narg, const char* argv[]) { 

  SmartChain* chain = new SmartChain("susy"); 
  for (int iii = 1; iii < narg; iii++) { 
    printf("file: %s, %i of %i\n", argv[iii], iii, narg - 1); 
    chain->add(argv[iii]); 
  }
  SusyBuffer buffer(chain); 

  const bool is_data = buffer.is_data(); 
  printf("running on %s\n", is_data ? "data":"MC"); 

  // ------ initialize susytools here -----------------
  SUSYObjDef* def = new SUSYObjDef; 
  def->initialize(is_data, true); // not data, atlfast
  printf("initalized\n"); 
  Root::TGRLCollection* grl = 0; 
  if (is_data) { 
    std::string grl_name = g_grl_file; 
    if (!exists(grl_name)) throw std::runtime_error(grl_name + " not found");
    Root::TGoodRunsListReader reader(grl_name.c_str(),true); 
    reader.Interpret(); 
    grl = new Root::TGRLCollection(reader.GetMergedGRLCollection()); 
  }

  CtagCalibration* ctag_cal = 0; 
  if (!is_data) { 
    try { 
      ctag_cal = new CtagCalibration(g_btag_file); 
    } catch (std::runtime_error& err) { 
      printf(red("disabled tagging SF: %s\n").c_str(), err.what()); 
    }
  }

  CutCounter counter; 
  CutCounter signal_counter; 
  CutCounter signal_counter_mc_wt; 
  CutCounter signal_counter_ctag_wt; 
  CutCounter el_cr_counter; 
  CutCounter mu_cr_counter; 

  // const long long max_entries = 10LL; 
  const long long max_entries = 100000000LL; 
  const unsigned n_entries = std::min(chain->GetEntries(),max_entries); 
  printf("looping over %i entries\n", n_entries); 
  for (unsigned nnn = 0; nnn < n_entries; nnn++) {
    if (nnn % 1000 == 0) { 
      printf("%ik entries processed\n", nnn/1000); 
    }
    def->Reset(); 
    chain->GetEntry(nnn); 
    
    // ---- start filling objects here -----
    std::vector<IdLorentzVector> all_jets; 
    for (int jeti = 0; jeti < buffer.jet_n; jeti++) { 
      def->FillJet(
	jeti, 
	buffer.jet_pt                 ->at(jeti), 
	buffer.jet_eta                ->at(jeti), 
	buffer.jet_phi                ->at(jeti),
	buffer.jet_E                  ->at(jeti), 
	buffer.jet_constscale_eta        ->at(jeti), 
	buffer.jet_constscale_phi        ->at(jeti), 
	buffer.jet_constscale_E        ->at(jeti), 
	buffer.jet_constscale_m        ->at(jeti),
	buffer.jet_ActiveAreaPx->at(jeti), 
	buffer.jet_ActiveAreaPy->at(jeti), 
	buffer.jet_ActiveAreaPz->at(jeti), 
	buffer.jet_ActiveAreaE->at(jeti), 
	buffer.Eventshape_rhoKt4LC, 
	buffer.averageIntPerXing,
	buffer.vx_nTracks); 

      if (!is_data) def->ApplyJetSystematics(
	jeti, 
	buffer.jet_constscale_eta        ->at(jeti), 
	buffer.jet_flavor_truth_label ->at(jeti), 
	buffer.averageIntPerXing,
	buffer.vx_nTracks, 
	SystErr::NONE); 

      bool good_jet = def->IsGoodJet(
	jeti, 
	buffer.jet_constscale_eta        ->at(jeti), 
	buffer.jet_emfrac             ->at(jeti), 
	buffer.jet_hecf               ->at(jeti),
	buffer.jet_LArQuality         ->at(jeti), 
	buffer.jet_HECQuality         ->at(jeti), 
	buffer.jet_AverageLArQF       ->at(jeti), 
	buffer.jet_Timing             ->at(jeti), 
	buffer.jet_sumPtTrk           ->at(jeti), 
	buffer.jet_fracSamplingMax    ->at(jeti),
	buffer.jet_SamplingMax        ->at(jeti), 
	buffer.jet_NegativeE          ->at(jeti), 
	buffer.RunNumber, 
	20e3, 	//pt cut 20 GeV
	10,	//eta cut (sort of nonexistent except for extreme cases)
	JetID::VeryLooseBad);
      TLorentzVector jet_tlv = def->GetJetTLV(); 
      IdLorentzVector jet; 
      jet.SetPxPyPzE(jet_tlv.Px(), jet_tlv.Py(), jet_tlv.Pz(), jet_tlv.E()); 
      jet.index = jeti; 
      jet.pass = good_jet; 
      all_jets.push_back(jet); 
    } // end jet filling loop
    std::sort(all_jets.begin(),all_jets.end(),has_higher_pt); 

    std::vector<IdLorentzVector> all_electrons; 
    for (int eli = 0; eli < buffer.el_n; eli++) { 
      bool good_el = def->FillElectron(
	eli,
	buffer.el_eta                   ->at(eli), 
	buffer.el_phi                   ->at(eli), 
	buffer.el_cl_eta                ->at(eli),
	buffer.el_cl_phi                ->at(eli),
	buffer.el_cl_E                  ->at(eli),
	buffer.el_tracketa              ->at(eli),
	buffer.el_trackphi              ->at(eli),
	buffer.el_author                ->at(eli),
	buffer.el_mediumPP              ->at(eli),
	buffer.el_OQ                    ->at(eli),
	buffer.el_nPixHits              ->at(eli),
	buffer.el_nSCTHits              ->at(eli),
	buffer.el_MET_Egamma10NoTau_wet->at(eli).at(0), 
	10e3,			// et cut
	2.47);
      TLorentzVector el_tlv = def->GetElecTLV(); 
      IdLorentzVector electron; 
      electron.SetPxPyPzE(el_tlv.Px(), el_tlv.Py(), el_tlv.Pz(), el_tlv.E()); 
      electron.index = eli; 
      electron.pass = good_el; 
      all_electrons.push_back(electron); 
    } // end electron filling loop

    std::vector<IdLorentzVector> all_muons; 
    for (int mui = 0; mui < buffer.mu_staco_n; mui++) { 
      bool good_muon = def->FillMuon
	(mui,
	 buffer.mu_staco_pt                           ->at(mui),
	 buffer.mu_staco_eta                          ->at(mui),
	 buffer.mu_staco_phi                          ->at(mui),
	 buffer.mu_staco_me_qoverp_exPV               ->at(mui),
	 buffer.mu_staco_id_qoverp_exPV               ->at(mui),
	 buffer.mu_staco_me_theta_exPV                ->at(mui),
	 buffer.mu_staco_id_theta_exPV                ->at(mui),
	 buffer.mu_staco_id_theta                     ->at(mui),
	 buffer.mu_staco_charge                       ->at(mui), 
	 buffer.mu_staco_isCombinedMuon               ->at(mui),
	 buffer.mu_staco_isSegmentTaggedMuon          ->at(mui),
	 buffer.mu_staco_loose                        ->at(mui),
	 buffer.mu_staco_nPixHits                     ->at(mui),
	 buffer.mu_staco_nPixelDeadSensors            ->at(mui),
	 buffer.mu_staco_nPixHoles                    ->at(mui),
	 buffer.mu_staco_nSCTHits                     ->at(mui),
	 buffer.mu_staco_nSCTDeadSensors              ->at(mui),
	 buffer.mu_staco_nSCTHoles                    ->at(mui),
	 buffer.mu_staco_nTRTHits                     ->at(mui),
	 buffer.mu_staco_nTRTOutliers                 ->at(mui),
	 10e3, 
	 2.4);
      TLorentzVector muon_tlv = def->GetMuonTLV(mui); 
      IdLorentzVector muon; 
      muon.SetPxPyPzE(
	muon_tlv.Px(), muon_tlv.Py(), muon_tlv.Pz(), muon_tlv.E()); 
      muon.index = mui; 
      muon.pass = good_muon; 
      all_muons.push_back(muon); 
    } // end muon filling loop

    //  ----- object preselection ------
    std::vector<IdLorentzVector> preselected_el = filter_pass(all_electrons); 
    std::vector<IdLorentzVector> preselected_mu = filter_pass(all_muons); 
    std::vector<IdLorentzVector> preselected_jets; 
    for (std::vector<IdLorentzVector>::const_iterator jitr = all_jets.begin(); 
	 jitr != all_jets.end(); jitr++) { 
      bool is_good_pt = jitr->Pt() > 20e3; 
      bool is_good_eta = std::abs(jitr->Eta()) < 2.8; 
      if (is_good_eta && is_good_pt) { 
	preselected_jets.push_back(*jitr); 
      }
    }

    counter["preselected_jets"] += preselected_jets.size(); 
    counter["preselected_el"] += preselected_el.size(); 
    counter["preselected_mu"] += preselected_mu.size(); 


    // ---- overlap removal ------
    std::vector<IdLorentzVector> after_overlap_jets = remove_overlaping(
      preselected_el, preselected_jets, 0.2); 
    std::vector<IdLorentzVector> after_overlap_el = remove_overlaping(
      after_overlap_jets, preselected_el, 0.4); 
    std::vector<IdLorentzVector> after_overlap_mu = remove_overlaping(
      after_overlap_jets, preselected_mu, 0.4); 
    counter["after_overlap_jets"] += after_overlap_jets.size(); 
    counter["after_overlap_el"]   += after_overlap_el.size(); 
    counter["after_overlap_mu"]   += after_overlap_mu.size(); 

    // ---- veto object selection -----
    // selection objects are used in various cutflows
    SelectionObjects so; 
    so.is_data = is_data; 
    so.veto_jets = filter_fail(after_overlap_jets); 
    for (std::vector<IdLorentzVector>::const_iterator 
	   itr = after_overlap_el.begin(); 
	 itr != after_overlap_el.end(); itr++) { 
      float rel_isolation = (
	buffer.el_ptcone20->at(itr->index) / itr->Pt()); 
      bool isolated_el = rel_isolation < 0.1; 
      if (isolated_el) {
	so.veto_electrons.push_back(*itr); 
      }
    }
    for (std::vector<IdLorentzVector>::const_iterator
	   itr = after_overlap_mu.begin(); 
	 itr != after_overlap_mu.end(); itr++) { 
      float isolation = buffer.mu_staco_ptcone20->at(itr->index); 
      bool isolated_mu = isolation < 1.8e3; 
      if (isolated_mu) { 
	so.veto_muons.push_back(*itr); 
      }
    }
    counter["veto_jets"] += so.veto_jets.size(); 
    counter["veto_electrons"] += so.veto_electrons.size(); 
    counter["veto_muons"] += so.veto_muons.size(); 

    // ---- signal object selection -----
    std::vector<IdLorentzVector> good_jets = filter_pass(after_overlap_jets);
    for (std::vector<IdLorentzVector>::const_iterator
	   itr = good_jets.begin(); 
	 itr != good_jets.end(); itr++) { 
      bool signal_pt = itr->Pt() > 30e3; 
      bool tag_eta = std::abs(itr->Eta()) < 2.5; 
      float jet_jvf = buffer.jet_jvtxf->at(itr->index); 
      bool ok_jvf = (jet_jvf > 0.5) || (itr->Pt() > 50e3); 
      if (signal_pt && tag_eta && ok_jvf) { 
	so.signal_jets.push_back(*itr); 
      }
    }
    for (std::vector<IdLorentzVector>::const_iterator
	   itr = so.veto_electrons.begin(); 
	 itr != so.veto_electrons.end(); itr++) { 
      bool control_pt = itr->Pt() > 20e3; 
      if (control_pt) { 
	so.control_electrons.push_back(*itr); 
      }
    }
    for (std::vector<IdLorentzVector>::const_iterator
	   itr = so.veto_muons.begin(); 
	 itr != so.veto_muons.end(); itr++) { 
      bool control_pt = itr->Pt() > 10e3; 
      if (control_pt) { 
	so.control_muons.push_back(*itr); 
      }
    }
    counter["good_jets"] += good_jets.size(); 
    counter["signal_jets"] += so.signal_jets.size(); 
    counter["control_electrons"] += so.control_electrons.size(); 
    counter["control_muons"] += so.control_muons.size(); 

    // --- calculate met ----
    std::vector<int> preselected_mu_idx; 
    for (std::vector<IdLorentzVector>::const_iterator 
	   itr = preselected_mu.begin(); itr != preselected_mu.end(); 
	 itr++){ 
      preselected_mu_idx.push_back(itr->index); 
    }
    std::vector<int> el_index; 
    for (int eln = 0; eln < buffer.el_n; eln++) { 
      float wet = buffer.el_MET_Egamma10NoTau_wet->at(eln).at(0); 
      if (wet != 0.0) el_index.push_back(eln); 
    }
    so.met = def->GetMET(
      buffer.jet_MET_Egamma10NoTau_wet,
      buffer.jet_MET_Egamma10NoTau_wpx,
      buffer.jet_MET_Egamma10NoTau_wpy,
      buffer.jet_MET_Egamma10NoTau_statusWord,
      el_index,
      buffer.el_MET_Egamma10NoTau_wet,
      buffer.el_MET_Egamma10NoTau_wpx,
      buffer.el_MET_Egamma10NoTau_wpy,
      buffer.el_MET_Egamma10NoTau_statusWord,
      buffer.MET_Egamma10NoTau_CellOut_etx, //CellOut
      buffer.MET_Egamma10NoTau_CellOut_ety, //CellOut
      buffer.MET_Egamma10NoTau_CellOut_sumet, //CellOut
      buffer.MET_CellOut_Eflow_STVF_etx, 
      buffer.MET_CellOut_Eflow_STVF_ety,
      buffer.MET_CellOut_Eflow_STVF_sumet,		  
      buffer.MET_Egamma10NoTau_RefGamma_etx,
      buffer.MET_Egamma10NoTau_RefGamma_ety,
      buffer.MET_Egamma10NoTau_RefGamma_sumet,
      preselected_mu_idx, 
      buffer.mu_staco_ms_qoverp, 
      buffer.mu_staco_ms_theta, 
      buffer.mu_staco_ms_phi, 
      buffer.mu_staco_charge, 
      buffer.mu_staco_energyLossPar,
      buffer.averageIntPerXing); 
    so.mu_met = so.met; 
    if (so.control_muons.size() == 1) { 
      TLorentzVector met_4vec; 
      met_4vec.SetPtEtaPhiE(so.met.Mod(), 0, so.met.Phi(), so.met.Mod()); 
      met_4vec += so.control_muons.at(0); 
      so.mu_met.Set(met_4vec.Px(), met_4vec.Py()); 
    }
    if (so.control_electrons.size() == 1) { 
      float min_delta_r = 10; 
      for (std::vector<IdLorentzVector>::const_iterator 
	     itr = all_jets.begin(); itr != all_jets.end(); itr++) {
	float delta_r = so.control_electrons.at(0).DeltaR(*itr); 
	if (delta_r < min_delta_r) { 
	  so.electron_jet = *itr; 
	  min_delta_r = delta_r; 
	}
      }
    }

    // ---- grl -----
    if (grl) { 
      so.pass_grl = grl->HasRunLumiBlock(buffer.RunNumber, buffer.lbn);
    } else { 
      so.pass_grl = true; 
    }

    // ---- event weights ----
    double ctag_wt = 1; 
    if (ctag_cal) {
      if (so.signal_jets.size() > 0) { 
	ctag_wt *= get_ctag_sf(so.signal_jets.at(0), buffer, *ctag_cal); 
      }
      if (so.signal_jets.size() > 1) { 
	ctag_wt *= get_ctag_sf(so.signal_jets.at(1), buffer, *ctag_cal); 
      }
    }

    // ---- run the event-wise selections -----
    signal_selection(so, def, buffer, signal_counter); 
    signal_selection(so, def, buffer, signal_counter_mc_wt, 
		     buffer.mcevt_weight->at(0).at(0)); 
    signal_selection(so, def, buffer, signal_counter_ctag_wt, ctag_wt); 
    el_cr_selection(so, def, buffer, el_cr_counter); 
    mu_cr_selection(so, def, buffer, mu_cr_counter); 

  } // end of event loop

  // ------ dump results ------
  dump_counts(counter, "objects"); 
  dump_counts(signal_counter, "signal region"); 
  dump_counts(signal_counter_mc_wt, "signal region evt wt"); 
  dump_counts(signal_counter_ctag_wt, "signal region tag wt"); 
  // dump_counts(el_cr_counter, "el control region"); 
  // dump_counts(mu_cr_counter, "mu control region"); 
}

// ----- common preselection runs before the other event wise selections
// will return false if a cut fails
bool common_preselection(const SelectionObjects& so, SUSYObjDef* def, 
			 const SusyBuffer& buffer, CutCounter& counter, 
			 double weight) { 
  counter["total"] += weight; 

  if (!so.pass_grl) return false; 
  counter["grl"] += weight; 

  if (!buffer.trigger) return false; 
  counter["trigger"] += weight; 

  bool pass_vxp = def->IsGoodVertex(buffer.vx_nTracks); 
  if (!pass_vxp) return false; 
  counter["primary_vertex"] += weight; 

  if (so.veto_jets.size()) return false; 
  counter["bad_jet_veto"] += weight; 

  bool has_tile_trip = def->IsTileTrip(buffer.RunNumber, buffer.lbn, 
               buffer.EventNumber); 
  if (has_tile_trip) return false; 
  counter["tile_trip"] += weight; 

  if (buffer.tileError == 2 && so.is_data) return false; 
  counter["tile_error"] += weight; 

  bool has_lar_error = (buffer.larError == 2); 
  if(has_lar_error && so.is_data) return false; 
  counter["lar_error"] += weight; 
    
  if (buffer.coreFlags & 0x40000 && so.is_data) return false; 
  counter["core_flags"] += weight; 
  
  return true; 
}

// ===============================================
// ========= eventwise signal selection ==========
// ===============================================

void signal_selection(const SelectionObjects& so, SUSYObjDef* def, 
		      const SusyBuffer& buffer, CutCounter& counter, 
		      double weight){
    
  bool pass_preselection = common_preselection(
    so, def, buffer, counter, weight); 
  if (!pass_preselection) return; 

  if (so.veto_electrons.size()) return; 
  counter["electron_veto"] += weight; 

  if (so.veto_muons.size()) return; 
  counter["muon_veto"] += weight; 

  if (so.veto_jets.size()) return; 
  counter["bad_jet_veto"] += weight; 
    
  const size_t n_jets = 2; 
  if (so.signal_jets.size() < n_jets) return; 
  counter["n_jet"] += weight; 
    
  if (so.met.Mod() < 150e3) return; 
  counter["met_150"] += weight; 

  const size_t n_jets = 2; 
  if (so.signal_jets.size() < n_jets) return; 
  counter["n_jet"] += weight; //Will's label: Minimum jet multiplicity
    
  if (so.signal_jets.at(0).Pt() < 130e3) return; 
  counter["leading_jet_130"] += weight; 

  if (so.signal_jets.at(1).Pt() < 50e3) return; 
  counter["second_jet_50"] += weight; 

  if (so.signal_jets.size() > 2) {
    if (so.signal_jets.at(2).Pt() > 50e3) return;
  }
  counter["third_jet_veto50"] += weight; 

  std::vector<size_t> jet_indices
  for (std::vector<IdLorentzVector>::const_iterator 
	 itr = so.signal_jets.begin(); itr != so.signal_jets.end(); itr++) { 
    jet_indices.push_back(itr->index); 
  }
  bool clean_for_chf = ChfCheck(jet_indices, buffer, *def); 
  if (clean_for_chf) return; 
  counter["pass_chf"] += weight; 

  bool medium_first = has_medium_tag(so.signal_jets.at(0).index, buffer); 
  bool medium_second = has_medium_tag(so.signal_jets.at(1).index, buffer); 
  if (! (medium_first && medium_second) ) return; 
  counter["two_ctag"] += weight; 

  TLorentzVector met_4vec; 
  met_4vec.SetPtEtaPhiE(1, 0, so.met.Phi(), 1); 
  float min_dphi = 1000; 
  for (std::vector<IdLorentzVector>::const_iterator 
	 itr = so.signal_jets.begin(); itr < so.signal_jets.begin() + n_jets; 
       itr++){
    float deltaphi = std::abs(met_4vec.DeltaPhi(*itr)); 
    min_dphi = std::min(deltaphi, min_dphi); 
  }

  if (min_dphi < 0.4) return; 
  counter["dphi_jetmet_min"] += weight; 

  double mass_eff = so.met.Mod() + scalar_sum_pt(so.signal_jets, 2); 
  if (so.met.Mod() / mass_eff < 0.25) return; 
  counter["met_eff"] += weight; 

  double mass_ct = get_mctcorr(so.signal_jets.at(0), so.signal_jets.at(1), so.met); 
  if (mass_ct < 150e3) return; 
  counter["m_ct_150"] += weight; 
  
  double mass_bb = (so.signal_jets.at(0) + so.signal_jets.at(1)).M(); 
  if (mass_bb < 200e3) return; 
  counter["m_bb"] += weight; 

  double mass_ct = get_m_ct(so.signal_jets.at(0), so.signal_jets.at(1)); 
  if (mass_ct < 150e3) return; 
  counter["m_ct_150"] += weight; 
  
} // end of signal region cutflow

void el_cr_selection(const SelectionObjects& so, SUSYObjDef* def, 
		     const SusyBuffer& buffer, CutCounter& counter){
  bool pass_preselection = common_preselection(so, def, buffer, counter, 1.0); 
  if (!pass_preselection) return; 
    
  if (so.veto_muons.size()) return; 
  counter["muon_veto"]++; 

  if (so.veto_jets.size()) return; 
  counter["bad_jet_veto"]++; 

  if (so.control_electrons.size() != 1) return; 
  counter["electron_jet"]++; 
  std::vector<IdLorentzVector> jets = so.signal_jets; 
  jets.push_back(so.electron_jet); 
  std::sort(jets.begin(), jets.end(), has_higher_pt); 
    
  const size_t n_jets = 3; 
  if (jets.size() < n_jets) return; 
  counter["n_jet"]++; 
    
  TLorentzVector met_4vec; 
  met_4vec.SetPtEtaPhiE(1, 0, so.met.Phi(), 1); 
  float min_dphi = 1000; 
  for (std::vector<IdLorentzVector>::const_iterator 
	 itr = jets.begin(); itr < jets.begin() + n_jets; 
       itr++) { 
    float deltaphi = std::abs(met_4vec.DeltaPhi(*itr)); 
    min_dphi = std::min(deltaphi, min_dphi); 
  }
  if (min_dphi < 0.4) return; 
  counter["dphi_jetmet_min"]++; 
    
  if (so.met.Mod() < 150e3) return; 
  counter["met_150"]++; 
  
  if (jets.at(0).Pt() < 150e3) return; 
  counter["leading_jet_150"]++; 

} // end of el cr cutflow

void mu_cr_selection(const SelectionObjects& so, SUSYObjDef* def, 
		     const SusyBuffer& buffer, CutCounter& counter){

  bool pass_preselection = common_preselection(so, def, buffer, counter, 1.0); 
  if (!pass_preselection) return; 
    
  if (so.veto_electrons.size()) return; 
  counter["electron_veto"]++; 

  if (so.veto_muons.size() != 1) return; 
  counter["muon_requirement"]++; 

  if (so.veto_jets.size()) return; 
  counter["bad_jet_veto"]++; 
    
  const size_t n_jets = 3; 
  if (so.signal_jets.size() < n_jets) return; 
  counter["n_jet"]++; 
    
  TLorentzVector met_4vec; 
  met_4vec.SetPtEtaPhiE(1, 0, so.met.Phi(), 1); 
  float min_dphi = 1000; 
  for (std::vector<IdLorentzVector>::const_iterator 
	 itr = so.signal_jets.begin(); itr < so.signal_jets.begin() + n_jets; 
       itr++) { 
    float deltaphi = std::abs(met_4vec.DeltaPhi(*itr)); 
    min_dphi = std::min(deltaphi, min_dphi); 
  }
  if (min_dphi < 0.4) return; 
  counter["dphi_jetmet_min"]++; 
  
  if (so.mu_met.Mod() < 150e3) return; 
  counter["mu_met_150"]++; 
    
  if (so.signal_jets.at(0).Pt() < 150e3) return; 
  counter["leading_jet_150"]++; 

} // end of mu control region cutflow


// ===== selection functions =======

std::vector<IdLorentzVector> filter_pass(
  const std::vector<IdLorentzVector>& in) { 
  std::vector<IdLorentzVector> out; 
  for (std::vector<IdLorentzVector>::const_iterator itr = in.begin(); 
       itr != in.end(); itr++) { 
    if (itr->pass) { 
      out.push_back(*itr); 
    }
  }
  return out; 
}

std::vector<IdLorentzVector> filter_fail(
  const std::vector<IdLorentzVector>& in) { 
  std::vector<IdLorentzVector> out; 
  for (std::vector<IdLorentzVector>::const_iterator itr = in.begin(); 
       itr != in.end(); itr++) { 
    if (!itr->pass) { 
      out.push_back(*itr); 
    }
  }
  return out; 
}

template<typename M, typename A>
A remove_overlaping(const M& mask, A altered, const float delta_r) { 
  for (typename M::const_iterator mitr = mask.begin(); 
       mitr != mask.end(); mitr++) { 
    A new_container; 
    for (typename M::const_iterator vic = altered.begin(); 
	 vic != altered.end(); vic++) { 
      assert(mitr->Pt() > 0); 
      double delr = mitr->DeltaR(*vic); 
      if (delr > delta_r) { 
	new_container.push_back(*vic); 
      }
    }
    altered = new_container; 
  }
  return altered; 
} 

// ================= calc functions ==================

bool has_medium_tag(int jet_index, const SusyBuffer& buffer) { 
  double pb = buffer.jet_flavor_component_jfitc_pb->at(jet_index); 
  double pc = buffer.jet_flavor_component_jfitc_pc->at(jet_index); 
  double pu = buffer.jet_flavor_component_jfitc_pu->at(jet_index); 

  // medium tag values are defined in ctag_defs.hh
  if (log(pc / pu) < JFC_MEDIUM_ANTI_U_CUT) return false; 
  if (log(pb / pu) < JFC_MEDIUM_ANTI_B_CUT) return false; 
  return true; 
}

bool has_higher_pt(const TLorentzVector& v1, const TLorentzVector& v2) { 
  return v1.Pt() > v2.Pt(); 
}

double scalar_sum_pt(const std::vector<IdLorentzVector>& obj, size_t num) { 
  int n_jets = std::min(num, obj.size()); 
  double sum = 0.0; 
  for (std::vector<IdLorentzVector>::const_iterator itr = obj.begin(); 
       itr < obj.begin() + n_jets; itr++) { 
    sum += itr->Pt(); 
  }
  return sum; 
}

double get_m_ct(const IdLorentzVector& v1, const IdLorentzVector& v2) { 
  double et1 = v1.Pt(); 
  double et2 = v2.Pt(); 
  TVector2 diff_pt = v1.Vect().XYvector() - v2.Vect().XYvector(); 
  double mct2 = std::pow(et1 + et2, 2) - std::pow(diff_pt.Mod(), 2); 
  return std::sqrt(mct2); 
}

double get_mctcorr(const TLorentzVector& tv1, const TLorentzVector& tv2, const TVector2& vmet)
{
  mctlib mct_object;
  
  double v1[4] = {tv1.E(), tv1.Px(), tv1.Py(), tv1.Pz()};
  double v2[4] = {tv2.E(), tv2.Px(), tv2.Py(), tv2.Pz()};
  double vds[4] = {0.0, 0.0, 0.0, 0.0};
  double ptm[2] = {vmet.X(), vmet.Y()};
  return mct_object.mctcorr(v1, v2, vds, ptm, 8000000.0, 0.0);
}

// ================= reweighting functions ============
void dump(const JetTagFactorInputs& inputs) { 
  printf("pt %f, eta %f, anti_b %f, anti_u %f, flavor %i\n", 
	 inputs.pt, inputs.eta, inputs.anti_b, inputs.anti_u, inputs.flavor);
}

double get_ctag_sf(const IdLorentzVector& jet, const SusyBuffer& buffer, 
		   const CtagCalibration& ctag_cal){ 
  double pb = buffer.jet_flavor_component_jfitc_pb->at(jet.index); 
  double pc = buffer.jet_flavor_component_jfitc_pc->at(jet.index); 
  double pu = buffer.jet_flavor_component_jfitc_pu->at(jet.index); 
  int ftl = buffer.jet_flavor_truth_label->at(jet.index); 
  
  JetTagFactorInputs inputs; 
  inputs.pt = jet.Pt(); 
  inputs.eta = jet.Eta(); 
  inputs.anti_b = log(pc / pb); 
  inputs.anti_u = log(pc / pu); 
  inputs.flavor = get_flavor(ftl); 
  // dump(inputs); 
  return ctag_cal.scale_factor(inputs).nominal; 
}


// ================= IO functions ==================

void dump_counts(const CutCounter& counter, std::string name) { 
  printf(" ========== %s ========== \n",name.c_str()); 
  typedef std::vector<std::pair<std::string, double> > OrdCuts; 
  OrdCuts ordered_cuts = counter.get_ordered_cuts(); 
  for (OrdCuts::const_iterator itr = ordered_cuts.begin(); 
       itr != ordered_cuts.end(); itr++) { 
    size_t pad_size = 20; 
    std::string name = itr->first; 
    if (name.size() < pad_size) { 
      name.insert(0, pad_size - name.size(), ' '); 
    }
    printf("%s: %.1f\n", name.c_str(), itr->second); 
  }
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

std::string red(std::string st) { 
  return "\033[31;1m" + st + "\033[m"; 
}
