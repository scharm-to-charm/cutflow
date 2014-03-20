#include <cstdio>
#include "SmartChain.hh"
#include "SusyBuffer.h"
#include "CutCounter.hh"
#include "ctag_defs.hh"		// for JFC_MEDIUM_* cuts
#include "CtagCalibration.hh"
#include "mctlib.hh"
#include "sbottom_functions.hh"

#include "SUSYTools/SUSYObjDef.h"
#include "TLorentzVector.h"
#include "GoodRunsLists/TGRLCollection.h"
#include "GoodRunsLists/TGoodRunsListReader.h"
#include "egammaAnalysisUtils/egammaTriggerMatching.h"

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
// ================= output ================
// list of set branches (set to "" for std::cout)
const std::string g_set_branches = "branches.txt"; 
// const std::string g_set_branches = "/dev/null"; 
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
  std::vector<IdLorentzVector> after_overlap_el;
  std::vector<IdLorentzVector> after_overlap_mu;
  std::vector<IdLorentzVector> signal_electrons; 
  std::vector<IdLorentzVector> signal_muons; 
  std::vector<IdLorentzVector> signal_jets; 
  TVector2 met; 

  TVector2 mu_met; 
  IdLorentzVector electron_jet; 
  bool is_data; 

  bool pass_grl; 
  bool has_cosmic_muon; 
  bool has_bad_muon; 
  bool has_bad_tile;
  bool has_trigger_matched_muon; 
  bool has_trigger_matched_electron; 
  bool has_dilep_trigger_matched_muon; 
  bool has_dilep_trigger_matched_electron; 
  double energy_weighted_time; 
}; 

// multiple object selections are called within the event loop
void signal_selection(const SelectionObjects&, SUSYObjDef* def, 
		      const SusyBuffer& buffer, CutCounter& counter, 
		      double weight = 1.0); 
void cra_1l_selection(const SelectionObjects&, SUSYObjDef* def, 
		      const SusyBuffer& buffer, CutCounter& counter, 
		      double weight = 1.0); 
void cra_sf_selection(const SelectionObjects&, SUSYObjDef* def, 
		      const SusyBuffer& buffer, CutCounter& counter, 
		      double weight = 1.0); 
void cra_of_selection(const SelectionObjects&, SUSYObjDef* def, 
		      const SusyBuffer& buffer, CutCounter& counter, 
		      double weight = 1.0); 

// these functions check to see if the object passed the SUSYObjDef cuts
std::vector<IdLorentzVector> filter_pass(const std::vector<IdLorentzVector>&); 
std::vector<IdLorentzVector> filter_fail(const std::vector<IdLorentzVector>&); 

// function to get indices (functions inherited from sbottom want them)
std::vector<size_t> get_indices(const std::vector<IdLorentzVector>&); 

// check for c-tags
bool has_medium_tag(const IdLorentzVector& jet, const SusyBuffer& buffer); 
// for sorting
bool has_higher_pt(const TLorentzVector& v1, const TLorentzVector& v2); 
// scalar sum for first n jets
double scalar_sum_pt(const std::vector<IdLorentzVector>& obj, size_t num); 
// opposite sign same flavor selection (for control regions)
bool has_os_sf_pair(const std::vector<IdLorentzVector>& electrons, 
		    const std::vector<IdLorentzVector>& muons, 
		    const SusyBuffer& buffer); 
// opposite sign different flavor selection (for control regions)
bool has_os_of_pair(const std::vector<IdLorentzVector>& electrons, 
		    const std::vector<IdLorentzVector>& muons, 
		    const SusyBuffer& buffer); 
// m_ct function 
double get_m_ct(const IdLorentzVector& v1, const IdLorentzVector& v2); 
double get_mctcorr(const TLorentzVector& v1, const TLorentzVector& v2, 
		   const TLorentzVector& vds, const TVector2& vmet);
double get_mt(const TLorentzVector& lep, const TVector2& met); 

// delta phi between jets and met
double get_min_dphi(const std::vector<IdLorentzVector>& jets, 
		    const TVector2& met); 

// ctag sf function (wrapper for CtagCalibration)
double get_ctag_sf(const IdLorentzVector& jet, const SusyBuffer& buffer, 
		   const CtagCalibration& ctag_cal); 
// bad tile veto 
bool has_bad_tile(const std::vector<IdLorentzVector>& jets, 
		  const TVector2& met, const SusyBuffer& buffer); 

// energy weighted time
double energy_weighted_time(const std::vector<IdLorentzVector>& jets,
			    int njets, const SusyBuffer& buffer);


// trigger match checks 
bool has_trigger_matched_electron(const std::vector<IdLorentzVector>& el, 
				  const SusyBuffer& buffer); 
bool has_trigger_matched_muon(const std::vector<IdLorentzVector>& mu, 
			      const SusyBuffer& buffer, 
			      SUSYObjDef& def); 
bool has_dilep_trigger_matched_electron(const std::vector<IdLorentzVector>&, 
					const SusyBuffer& buffer); 
bool has_dilep_trigger_matched_muon(const std::vector<IdLorentzVector>& mu, 
				    const SusyBuffer& buffer, 
				    SUSYObjDef& def); 

// IO functions
void dump_counts(const CutCounter&, std::string); 
bool exists(std::string file_name); 
std::string red(std::string); 
void dump_branches(const std::vector<std::string>& branch_names, 
		   const std::string file_name = ""); 

template<typename M, typename A>
A remove_overlaping(const M& mask, A altered, const float delta_r); 

// =================================================================
// =============== main event loop starts here =====================
// =================================================================

int main (int narg, const char* argv[]) { 
  if (narg == 1) throw std::runtime_error("no files given"); 

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
  CutCounter cra_1l_counter; 
  CutCounter cra_sf_counter; 
  CutCounter cra_of_counter; 

  // const long long max_entries = 10LL; 
  const long long max_entries = 100000LL; 
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
	7e3,			// et cut
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
	 6e3, 
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
    // selection objects are used in various cutflows
    SelectionObjects so; 
    std::vector<IdLorentzVector> after_overlap_jets = remove_overlaping(
      preselected_el, preselected_jets, 0.2); 
    so.after_overlap_el = remove_overlaping(
      after_overlap_jets, preselected_el, 0.4); 
    so.after_overlap_mu = remove_overlaping(
      after_overlap_jets, preselected_mu, 0.4); 
    counter["after_overlap_jets"] += after_overlap_jets.size(); 
    counter["after_overlap_el"]   += so.after_overlap_el.size(); 
    counter["after_overlap_mu"]   += so.after_overlap_mu.size(); 

    // ---- veto object selection -----
    so.is_data = is_data; 
    so.veto_jets = filter_fail(after_overlap_jets); 

    counter["veto_jets"] += so.veto_jets.size(); 

    // ---- signal object selection -----
    std::vector<IdLorentzVector> good_jets = filter_pass(after_overlap_jets);
    for (std::vector<IdLorentzVector>::const_iterator
	   itr = good_jets.begin(); 
	 itr != good_jets.end(); itr++) { 
      bool signal_pt = itr->Pt() > 20e3; // was 30, for mindphi(jet-MET)

      // the eta requirement for tagging used to be made here, now it's 
      // done when we check for tags. 
      bool ok_eta = std::abs(itr->Eta()) < 2.8; // was 2.5 (for tagging)

      float jet_jvf = buffer.jet_jvtxf->at(itr->index); 
      bool ok_jvf_frac = jet_jvf > 0.5; 
      bool ok_jvf = ( ok_jvf_frac || (itr->Pt() > 50e3) || 
		      (std::abs(itr->Eta()) > 2.4) );
      bool no_tracks = jet_jvf < -0.5; // no-track jets should have -1
      if (signal_pt && ok_eta && ( ok_jvf || no_tracks) )  { 
	so.signal_jets.push_back(*itr); 
      }
    }
    
    for (std::vector<IdLorentzVector>::const_iterator
	   itr = so.after_overlap_el.begin(); 
	 itr != so.after_overlap_el.end(); itr++) { 
      bool control_pt = itr->Pt() > 20e3 || true; 
      bool tight_pp = buffer.el_tightPP->at(itr->index);
      bool rel_iso = buffer.el_ptcone20->at(itr->index) / itr->Pt() < 0.1;
      if (control_pt && tight_pp && rel_iso) { 
	so.signal_electrons.push_back(*itr); 
      }
    }
    for (std::vector<IdLorentzVector>::const_iterator
	   itr = so.after_overlap_mu.begin(); 
	 itr != so.after_overlap_mu.end(); itr++) { 
      bool control_pt = itr->Pt() > 20e3 || true; 
      bool iso = buffer.mu_staco_ptcone20->at(itr->index) < 1.8e3;
      if (control_pt && iso) { 
	so.signal_muons.push_back(*itr); 
      }
    }
    counter["good_jets"] += good_jets.size(); 
    counter["signal_jets"] += so.signal_jets.size(); 
    counter["signal_electrons"] += so.signal_electrons.size(); 
    counter["signal_muons"] += so.signal_muons.size(); 

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
    if (so.signal_muons.size() == 1) { 
      TLorentzVector met_4vec; 
      met_4vec.SetPtEtaPhiE(so.met.Mod(), 0, so.met.Phi(), so.met.Mod()); 
      met_4vec += so.signal_muons.at(0); 
      so.mu_met.Set(met_4vec.Px(), met_4vec.Py()); 
    }
    if (so.signal_electrons.size() == 1) { 
      float min_delta_r = 10; 
      for (std::vector<IdLorentzVector>::const_iterator 
	     itr = all_jets.begin(); itr != all_jets.end(); itr++) {
	float delta_r = so.signal_electrons.at(0).DeltaR(*itr); 
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

    // --- muon quality cuts ---
    so.has_cosmic_muon = false; 
    for (std::vector<IdLorentzVector>::const_iterator 
	   itr = so.after_overlap_mu.begin(); 
	 itr != so.after_overlap_mu.end(); itr++){ 
      float mu_z0_exPV = buffer.mu_staco_z0_exPV->at(itr->index); 
      float mu_d0_exPV = buffer.mu_staco_d0_exPV->at(itr->index); 
      if (def->IsCosmicMuon(mu_z0_exPV, mu_d0_exPV)) { 
	so.has_cosmic_muon = true; 
	break; 
      }
    }
    
    so.has_bad_muon = false; 
    for (std::vector<IdLorentzVector>::const_iterator 
	   itr = preselected_mu.begin(); 
	 itr != preselected_mu.end(); itr++){ 
      float mu_qoverp = buffer.mu_staco_qoverp_exPV->at(itr->index); 
      float mu_cov_qoverp = buffer.mu_staco_cov_qoverp_exPV->at(itr->index); 
      if (def->IsBadMuon(mu_qoverp, mu_cov_qoverp)) { 
	so.has_bad_muon = true; 
	break; 
      }
    }
    
    // ---- bad tile cut ----
    so.has_bad_tile = has_bad_tile(preselected_jets, so.met, buffer); 
    
    // ---- energy weighted time ----
    // 2 jets in SRA-type regions
    so.energy_weighted_time = energy_weighted_time(so.signal_jets, 2, buffer);
    //preselected, after_overlap, good, signal

    // ---- trigger matching ----
    so.has_trigger_matched_muon = has_trigger_matched_muon(
      so.signal_muons, buffer, *def); 
    so.has_dilep_trigger_matched_muon = has_dilep_trigger_matched_muon(
      so.signal_muons, buffer, *def);
    so.has_trigger_matched_electron = has_trigger_matched_electron(
      so.signal_electrons, buffer); 
    so.has_dilep_trigger_matched_electron = has_dilep_trigger_matched_electron(
      so.signal_electrons, buffer);

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
		     buffer.get_mcevt_weight()); 
    signal_selection(so, def, buffer, signal_counter_ctag_wt, ctag_wt); 
    cra_1l_selection(so, def, buffer, cra_1l_counter); 
    cra_sf_selection(so, def, buffer, cra_sf_counter); 
    cra_of_selection(so, def, buffer, cra_of_counter); 

  } // end of event loop

  // ------ dump results ------
  dump_counts(counter, "objects"); 
  dump_counts(signal_counter, "signal region"); 
  dump_counts(signal_counter_mc_wt, "signal region evt wt"); 
  dump_counts(signal_counter_ctag_wt, "signal region tag wt"); 
  dump_counts(cra_1l_counter, "CRA 1L"); 
  dump_counts(cra_sf_counter, "CRA SF"); 
  dump_counts(cra_of_counter, "CRA DF"); 
  dump_branches(chain->get_all_branch_names(), g_set_branches); 
}

// ----- common preselection runs before the other event wise selections
// will return false if a cut fails
bool common_preselection(const SelectionObjects& so, SUSYObjDef* def, 
			 const SusyBuffer& buffer, CutCounter& counter, 
			 double weight) { 

  if (!so.pass_grl) return false; 
  counter["grl"] += weight; 

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

  if (so.has_bad_tile) return false; 
  counter["bad_tile_veto"] += weight; 

  bool has_lar_error = (buffer.larError == 2); 
  if(has_lar_error && so.is_data) return false; 
  counter["lar_error"] += weight; 
    
  // 0L has >4, sbottoms had >5
  if (std::abs(so.energy_weighted_time) > 5) return false; 
  counter["energy_weighted_time"] += weight;

  if (buffer.coreFlags & 0x40000 && so.is_data) return false; 
  counter["core_flags"] += weight; 

  if (so.has_cosmic_muon) return false; 
  counter["cosmic_muon_veto"] += weight; 

  if (so.has_bad_muon) return false; 
  counter["bad_muon_veto"] += weight; 
  
  return true; 
}

bool common_lepton_trigger_matching(const SelectionObjects& so, 
				    const SusyBuffer& buffer, 
				    CutCounter& counter, 
				    double weight) {
 
  bool mu_trig = buffer.EF_mu24i_tight || buffer.EF_mu36_tight; 
  bool el_trig = buffer.EF_e24vhi_medium1 || buffer.EF_e60_medium1; 
  if (! (mu_trig || el_trig) ) return false; 
  counter["pass_lepton_trigger"] += weight; 

  if (! ( (so.has_trigger_matched_electron && el_trig) || 
  	  (so.has_trigger_matched_muon && mu_trig ) ) ) return false; 
  counter["pass_lepton_trigger_match"] += weight; 
  
  return true; 
}

// ===============================================
// ========= eventwise signal selection ==========
// ===============================================

void signal_selection(const SelectionObjects& so, SUSYObjDef* def, 
		      const SusyBuffer& buffer, CutCounter& counter, 
		      double weight){
  counter["total"] += weight; 
    
  bool met_trig = buffer.xe80_tclcw_tight || buffer.xe80T_tclcw_loose ||
    buffer.xe80T_tclcw_loose; 
  if (!met_trig) return; 
  counter["met_trigger"] += weight; 

  bool pass_preselection = common_preselection(
    so, def, buffer, counter, weight); 
  if (!pass_preselection) return; 

  if (so.after_overlap_mu.size()) return; 
  counter["muon_veto"] += weight; 

  if (so.after_overlap_el.size()) return; 
  counter["electron_veto"] += weight; 

  bool clean_for_chf = ChfCheck(get_indices(so.signal_jets), buffer, *def); 
  if (clean_for_chf) return; 
  counter["pass_chf"] += weight; 
    
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

  double min_dphi = get_min_dphi(so.signal_jets, so.met); 
  if (min_dphi < 0.4) return; 
  counter["dphi_jetmet_min"] += weight; 

  double mass_eff = so.met.Mod() + scalar_sum_pt(so.signal_jets, 2); 
  if (so.met.Mod() / mass_eff < 0.25) return; 
  counter["met_eff"] += weight; 

  bool medium_first =  has_medium_tag(so.signal_jets.at(0), buffer);
  bool medium_second = has_medium_tag(so.signal_jets.at(1), buffer);

  if (! (medium_first || medium_second) ) return; 
  counter["at_least_one_ctag"] += weight; 
  if (! (medium_first && medium_second) ) return; 
  counter["two_ctag"] += weight; 


  double mass_ct = get_mctcorr(so.signal_jets.at(0), 
			       so.signal_jets.at(1), 
			       TLorentzVector(),
			       so.met); 
  if (mass_ct < 150e3) return; 
  counter["m_ct_150"] += weight; 

  double mass_cc = (so.signal_jets.at(0) + so.signal_jets.at(1)).M(); 
  if (mass_cc < 200e3) return; 
  counter["m_cc"] += weight; 
  
} // end of signal region cutflow

void cra_1l_selection(const SelectionObjects& so, SUSYObjDef* def, 
		      const SusyBuffer& buffer, CutCounter& counter, 
		      double weight){
  counter["total"] += weight; 

  bool pass_lepton_trigger = common_lepton_trigger_matching(
    so, buffer, counter, weight); 
  if (!pass_lepton_trigger) return; 

  bool pass_preselection = common_preselection(
    so, def, buffer, counter, weight); 
  if (!pass_preselection) return; 

  int total_leptons = so.signal_electrons.size() + so.signal_muons.size();
  if (total_leptons != 1) return; 
  counter["pass_1l"] += weight; 
  
  // control leptons are a subset of the veto leptons, so we shouldn't
  // have any additional veto leptons. 
  int total_veto_leptons = so.after_overlap_el.size() + so.after_overlap_mu.size(); 
  if (total_veto_leptons != 1) return; 
  counter["pass_lepton_veto"] += weight; 

  bool clean_for_chf = ChfCheck(get_indices(so.signal_jets), buffer, *def); 
  if (clean_for_chf) return; 
  counter["pass_chf"] += weight; 
    
  if (so.met.Mod() < 100e3) return; 
  counter["met_100"] += weight; 

  const size_t n_jets = 2; 
  if (so.signal_jets.size() < n_jets) return; 
  counter["n_jet"] += weight; //Will's label: Minimum jet multiplicity
  
  if (so.signal_jets.at(0).Pt() < 130e3) return; 
  counter["leading_jet_130"] += weight; 

  if (so.signal_jets.at(1).Pt() < 50e3) return; 
  counter["second_jet_50"] += weight; 

  bool medium_first =  has_medium_tag(so.signal_jets.at(0), buffer);
  bool medium_second = has_medium_tag(so.signal_jets.at(1), buffer);

  if (! (medium_first || medium_second) ) return; 
  counter["at_least_one_ctag"] += weight; 
  if (! (medium_first && medium_second) ) return; 
  counter["two_ctag"] += weight; 

  TLorentzVector lep = so.signal_muons.size() == 1 ? 
    so.signal_muons.at(0) : so.signal_electrons.at(0); 
  double mass_ct = get_mctcorr(so.signal_jets.at(0), 
			       so.signal_jets.at(1), lep, so.met); 
  if (mass_ct < 150e3) return; 
  counter["m_ct_150"] += weight; 

  double mt = get_mt(lep, so.met); 
  if (!(40e3 < mt && mt < 100e3)) return; 
  counter["mt"] += weight; 
  
} // end of cra_1l_selection

void cra_sf_selection(const SelectionObjects& so, SUSYObjDef* def, 
		      const SusyBuffer& buffer, CutCounter& counter, 
		      double weight){
  counter["total"] += weight; 

  bool pass_mu = buffer.EF_mu18_tight_mu8_EFFS;
  bool pass_el = buffer.EF_2e12Tvh_loose1;
  bool pass_lepton_trigger = pass_el || pass_mu;
  if (!pass_lepton_trigger) return; 
  counter["pass_two_lepton_trig"] += weight;

  bool trig_matched = so.has_dilep_trigger_matched_muon ||
    so.has_dilep_trigger_matched_electron;
  if (!trig_matched) return;
  counter["pass_trigger_match"] += weight;
    
  bool pass_preselection = common_preselection(
    so, def, buffer, counter, weight); 
  if (!pass_preselection) return; 

  int n_el = so.signal_electrons.size();
  int n_mu = so.signal_muons.size();
  int total_leptons = n_mu + n_el; 
  bool ossf_pair = has_os_sf_pair(
    so.signal_electrons, so.signal_muons, buffer); 
  if (!ossf_pair) return; 
  counter["pass_ossf"] += weight; 
  
  IdLorentzVector lep1 = n_mu == 2 ? 
    so.signal_muons.at(0) : so.signal_electrons.at(0); 
  IdLorentzVector lep2 = n_mu == 2 ? 
    so.signal_muons.at(1) : so.signal_electrons.at(1); 

  if (total_leptons != 2) return;
  counter["pass_2_lepton"] += weight;

  // signal leptons are a subset of the veto leptons
  int total_veto_leptons = so.after_overlap_el.size() + so.after_overlap_mu.size(); 
  if (total_veto_leptons != 2) return; 
  counter["pass_lepton_veto"] += weight; 

  bool clean_for_chf = ChfCheck(get_indices(so.signal_jets), buffer, *def); 
  if (clean_for_chf) return; 
  counter["pass_chf"] += weight; 
    
  TVector2 lept_pxy = lep1.Vect().XYvector() + lep2.Vect().XYvector();
  TVector2 lept_met = so.met + lept_pxy;
  if (lept_met.Mod() < 100e3) return; 
  counter["met_100"] += weight; 

  const size_t n_jets = 2; 
  if (so.signal_jets.size() < n_jets) return; 
  counter["n_jet"] += weight; //Will's label: Minimum jet multiplicity
  
  if (so.signal_jets.at(0).Pt() < 50e3) return; 
  counter["leading_jet_50"] += weight; 

  if (so.signal_jets.at(1).Pt() < 50e3) return; 
  counter["second_jet_50"] += weight; 

  bool medium_first =  has_medium_tag(so.signal_jets.at(0), buffer);
  bool medium_second = has_medium_tag(so.signal_jets.at(1), buffer);

  if (! (medium_first || medium_second) ) return; 
  counter["at_least_one_ctag"] += weight; 
  if (! (medium_first && medium_second) ) return; 
  counter["two_ctag"] += weight; 

  double m_ll = (lep1 + lep2).M(); 
  if (! (75e3 < m_ll && m_ll < 105e3)) return; 
  counter["mll_zpeak"] += weight; 

  double lepton_pt = std::max(lep1.Pt(), lep2.Pt()); 
  if (! lepton_pt > 70e3) return; 
  counter["lepton_pt_70"] += weight; 

  double mass_cc = (so.signal_jets.at(0) + so.signal_jets.at(1)).M(); 
  if (mass_cc < 200e3) return; 
  counter["m_cc"] += weight; 
  
} // end of cra_sf_selection

void cra_of_selection(const SelectionObjects& so, SUSYObjDef* def, 
		      const SusyBuffer& buffer, CutCounter& counter, 
		      double weight){
  counter["total"] += weight; 

  bool pass_lepton_trigger = common_lepton_trigger_matching(
    so, buffer, counter, weight); 
  if (!pass_lepton_trigger) return; 
    
  bool pass_preselection = common_preselection(
    so, def, buffer, counter, weight); 
  if (!pass_preselection) return; 

  // control leptons are a subset of the veto leptons, so we shouldn't
  // have any additional veto leptons. 
  int extra_leptons = 
    (so.after_overlap_el.size() + so.after_overlap_mu.size() ) - 
    (so.signal_electrons.size() + so.signal_muons.size());
  if (extra_leptons != 0) return; 
  counter["pass_lepton_veto"] += weight; 

  bool ofos_pair = has_os_of_pair(
    so.signal_electrons, so.signal_muons, buffer); 
  if (!ofos_pair) return; 
  counter["pass_ofos"] += weight; 

  bool clean_for_chf = ChfCheck(get_indices(so.signal_jets), buffer, *def); 
  if (clean_for_chf) return; 
  counter["pass_chf"] += weight; 
    
  if (so.met.Mod() < 50e3) return; 
  counter["met_50"] += weight; 

  const size_t n_jets = 2; 
  if (so.signal_jets.size() < n_jets) return; 
  counter["n_jet"] += weight; //Will's label: Minimum jet multiplicity
  
  if (so.signal_jets.at(0).Pt() < 50e3) return; 
  counter["leading_jet_50"] += weight; 

  if (so.signal_jets.at(1).Pt() < 50e3) return; 
  counter["second_jet_50"] += weight; 

  bool medium_first =  has_medium_tag(so.signal_jets.at(0), buffer);
  bool medium_second = has_medium_tag(so.signal_jets.at(1), buffer);

  if (! (medium_first || medium_second) ) return; 
  counter["at_least_one_ctag"] += weight; 
  if (! (medium_first && medium_second) ) return; 
  counter["two_ctag"] += weight; 

  TLorentzVector sum_lept = so.signal_electrons.at(0) + so.signal_muons.at(0);
  double mll = sum_lept.M(); 
  if (mll < 50e3) return; 
  counter["mll"] += weight; 
  
} // end of cra_of_selection


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

bool has_os_sf_pair(const std::vector<IdLorentzVector>& electrons, 
		    const std::vector<IdLorentzVector>& muons, 
		    const SusyBuffer& buffer) { 
  int n_el = electrons.size();
  int n_mu = muons.size();
  int total_leptons = n_mu + n_el; 
  if (total_leptons != 2) return false; 
  if (n_mu != 2 && n_el != 2) return false; 
  
  if (n_el == 2) { 
    float l1_charge = buffer.el_charge->at(electrons.at(0).index); 
    float l2_charge = buffer.el_charge->at(electrons.at(1).index); 
    return l1_charge * l2_charge < 0.0; 
  } else if (n_mu == 2) { 
    float l1_charge = buffer.mu_staco_charge->at(muons.at(0).index); 
    float l2_charge = buffer.mu_staco_charge->at(muons.at(1).index); 
    return l1_charge * l2_charge < 0.0; 
  }
  assert(false); 		// shouldn't get here

}
bool has_os_of_pair(const std::vector<IdLorentzVector>& electrons, 
		    const std::vector<IdLorentzVector>& muons, 
		    const SusyBuffer& buffer) { 
  int n_el = electrons.size();
  int n_mu = muons.size();
  if (n_el != 1) return false; 
  if (n_mu != 1) return false; 

  float el_charge = buffer.el_charge->at(electrons.at(0).index); 
  float mu_charge = buffer.mu_staco_charge->at(muons.at(0).index); 
  return el_charge * mu_charge < 0.0; 
}


// ================= calc functions ==================

bool has_medium_tag(const IdLorentzVector& jet, const SusyBuffer& buffer) { 
  int jet_index = jet.index; 
  double pb = buffer.jet_flavor_component_jfitc_pb->at(jet_index); 
  double pc = buffer.jet_flavor_component_jfitc_pc->at(jet_index); 
  double pu = buffer.jet_flavor_component_jfitc_pu->at(jet_index); 

  if (std::abs(jet.Eta()) > 2.5) return false; 

  // medium tag values are defined in ctag_defs.hh
  if (log(pc / pu) < JFC_MEDIUM_ANTI_U_CUT) return false; 
  if (log(pc / pb) < JFC_MEDIUM_ANTI_B_CUT) return false; 
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

double get_mctcorr(const TLorentzVector& tv1, const TLorentzVector& tv2, 
		   const TLorentzVector& tvd, 
		   const TVector2& vmet)
{
  mctlib mct_object;
  
  double v1[4] =  {tv1.E(), tv1.Px(), tv1.Py(), tv1.Pz()};
  double v2[4] =  {tv2.E(), tv2.Px(), tv2.Py(), tv2.Pz()};
  double vds[4] = {tvd.E(), tvd.Px(), tvd.Py(), tvd.Pz()};
  double ptm[2] = {vmet.X(), vmet.Y()};
  return mct_object.mctcorr(v1, v2, vds, ptm, 8e6, 0.0);
}

double get_mt(const TLorentzVector& lep, const TVector2& met) { 
  TVector2 lep2vec = lep.Vect().XYvector(); 
  return std::sqrt(2*lep2vec.Mod()*met.Mod() - 2*lep2vec*met);
}

double get_min_dphi(const std::vector<IdLorentzVector>& jets, 
		    const TVector2& met) { 
  assert(jets.size() >= 2); 

  TLorentzVector met_4vec; 
  met_4vec.SetPtEtaPhiE(1, 0, met.Phi(), 1); 
  float min_dphi = 1000; 
  for (std::vector<IdLorentzVector>::const_iterator 
	 itr = jets.begin(); itr < jets.begin() + 2;
       itr++){
    float deltaphi = std::abs(met_4vec.DeltaPhi(*itr)); 
    min_dphi = std::min(deltaphi, min_dphi); 
  }
  if (jets.size() > 2) {
    float deltaphi = std::abs(met_4vec.DeltaPhi(jets.at(2)));
    min_dphi = std::min(deltaphi, min_dphi);
  }
  return min_dphi; 
}

std::vector<size_t> get_indices(const std::vector<IdLorentzVector>& vecs) { 
  std::vector<size_t> indices; 
  for (std::vector<IdLorentzVector>::const_iterator itr = vecs.begin(); 
       itr != vecs.end(); itr++) { 
    indices.push_back(itr->index); 
  }
  return indices; 
}

// ================= trigger matching =================
// trigger match checks 
bool has_trigger_matched_electron(const std::vector<IdLorentzVector>& el, 
                                  const SusyBuffer& buffer) { 
  for (std::vector<IdLorentzVector>::const_iterator itr = el.begin(); 
       itr != el.end(); itr++) { 
    if (itr->Pt() < 25e3) continue; // trigger threshold from sbottom
    int nothing; 
    if (PassedTriggerEF(
	  itr->Eta(), itr->Phi(), buffer.trig_EF_el_EF_e24vhi_medium1, 
	  nothing, buffer.trig_EF_el_eta->size(), 
	  buffer.trig_EF_el_eta, buffer.trig_EF_el_phi)) return true; 
    if (PassedTriggerEF(
	  itr->Eta(), itr->Phi(), buffer.trig_EF_el_EF_e60_medium1, 
	  nothing, buffer.trig_EF_el_eta->size(),  
	  buffer.trig_EF_el_eta, buffer.trig_EF_el_phi)) return true; 
  }
  return false; 
}
bool has_dilep_trigger_matched_electron(
  const std::vector<IdLorentzVector>& el, const SusyBuffer& buffer) {
  int n_passing = 0;
  for (std::vector<IdLorentzVector>::const_iterator itr = el.begin(); 
       itr != el.end(); itr++) { 
    if (itr->Pt() < 25e3) continue; // trigger threshold from sbottom
    int nothing; 
    if (PassedTriggerEF(
	  itr->Eta(), itr->Phi(), buffer.trig_EF_el_EF_2e12Tvh_loose1, 
	  nothing, buffer.trig_EF_el_eta->size(), 
	  buffer.trig_EF_el_eta, buffer.trig_EF_el_phi)) n_passing++; 
  }
  if (n_passing >= 2) return true;
  return false; 
}
bool has_trigger_matched_muon(const std::vector<IdLorentzVector>& mu, 
                              const SusyBuffer& buffer, 
                              SUSYObjDef& def) { 
  for (std::vector<IdLorentzVector>::const_iterator itr = mu.begin(); 
       itr != mu.end(); itr++) { 
    if (itr->Pt() < 25e3) continue; // trigger threshold from sbottom
    int nothing; 
    if (def.MuonHasTriggerMatch(
          itr->Eta(), itr->Phi(), 
          buffer.trig_EF_trigmuonef_EF_mu24i_tight, 
          nothing, nothing, 
          buffer.trig_EF_trigmuonef_track_CB_eta->size(), 
          buffer.trig_EF_trigmuonef_track_CB_eta, 
          buffer.trig_EF_trigmuonef_track_CB_phi, 
          buffer.trig_EF_trigmuonef_track_CB_hasCB)) return true; 
    if (def.MuonHasTriggerMatch(
          itr->Eta(), itr->Phi(), 
          buffer.trig_EF_trigmuonef_EF_mu36_tight, 
          nothing, nothing, 
          buffer.trig_EF_trigmuonef_track_CB_eta->size(), 
          buffer.trig_EF_trigmuonef_track_CB_eta, 
          buffer.trig_EF_trigmuonef_track_CB_phi, 
          buffer.trig_EF_trigmuonef_track_CB_hasCB)) return true; 
  }
  return false; 
}
bool has_dilep_trigger_matched_muon(const std::vector<IdLorentzVector>& mu, 
                                    const SusyBuffer& buffer, 
                                    SUSYObjDef& def) { 
  int n_passing = 0;
  for (std::vector<IdLorentzVector>::const_iterator itr = mu.begin(); 
       itr != mu.end(); itr++) { 
    if (itr->Pt() < 20e3) continue; // trigger threshold from sbottom
    int nothing; 
    if (def.MuonHasTriggerMatch(
          itr->Eta(), itr->Phi(), 
          buffer.trig_EF_trigmuonef_EF_mu18_tight_mu8_EFFS, 
          nothing, nothing, 
          buffer.trig_EF_trigmuonef_track_CB_eta->size(), 
          buffer.trig_EF_trigmuonef_track_CB_eta, 
          buffer.trig_EF_trigmuonef_track_CB_phi, 
          buffer.trig_EF_trigmuonef_track_CB_hasCB)) n_passing++; 
  }
  if (n_passing >= 2) return true;
  return false; 
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

// ================= object / event quality ===========
bool has_bad_tile(const std::vector<IdLorentzVector>& jets, 
		  const TVector2& met, 
		  const SusyBuffer& buffer) { 
  for (std::vector<IdLorentzVector>::const_iterator itr = jets.begin(); 
       itr != jets.end(); itr++) { 
    float BCH_CORR_JET = buffer.jet_BCH_CORR_JET->at(itr->index); 
    float dphi = met.DeltaPhi(itr->Vect().XYvector()); 
    if (itr->Pt() > 40e3 && BCH_CORR_JET > 0.05 && std::abs(dphi) < 0.3) { 
      return true; 
    }
  }
  return false; 
}

double energy_weighted_time(const std::vector<IdLorentzVector>& jets, /////NB: different function when include
			    int njets, const SusyBuffer& buffer) {    /////    leptons in control regions
  double denom = 0.;
  double num = 0.;
  for (std::vector<IdLorentzVector>::const_iterator itr = jets.begin(); 
       itr != jets.end(); itr++) {
    njets --;
    if (njets < 0) break;
    denom += itr->E();
    num += (itr->E())*(buffer.jet_Timing->at(itr->index));
  }
  if (jets.size()==0)
    return 0.;
  else
    return num/denom;
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

void dump_branches(const std::vector<std::string>& branch_names, 
		   const std::string file_name) { 
  std::ofstream out; 
  if (file_name.size()) { 
    out.open(file_name.c_str(), std::ofstream::trunc); 
  } else { 
    std::cout << "======= all set branches ======\n"; 
  }
  for (std::vector<std::string>::const_iterator itr = branch_names.begin(); 
       itr != branch_names.end(); itr++) { 
    out << *itr << std::endl; 
    if (file_name.size() == 0) { 
      std::cout << *itr << std::endl;
    }
  }
}
