#include "PileupReweighting.hh"
#include "SusyBuffer.h"
#include "PileupReweighting/TPileupReweighting.h"

#include <fstream>
#include <stdexcept>

namespace {
  bool exists(std::string file_name);
}

PileupReweighting::PileupReweighting(const std::string& pu_config,
				     const std::string& pu_lumicalc):
  m_prw(0)
{
  if (!exists(pu_config)) throw std::runtime_error("missing " + pu_config);
  if (!exists(pu_lumicalc)) throw std::runtime_error(
    "missing " + pu_lumicalc);
  m_prw = new Root::TPileupReweighting("PileupReweighting");
  m_prw->SetDefaultChannel(0); // this is what Brett does
  m_prw->AddConfigFile(pu_config);

  // this is taken from:
  // twiki.cern.ch/twiki/bin/view/AtlasProtected/ExtendedPileupReweighting#Recipe_A_MC12a_Pileup_Reweightin
  // and confirmed by Alex.
  m_prw->SetDataScaleFactors(1/1.09);

  m_prw->AddLumiCalcFile(pu_lumicalc);
  m_prw->MergeMCRunNumbers(195847,195848);
  m_prw->SetUnrepresentedDataAction(2);
  m_prw->Initialize();

}

PileupReweighting::~PileupReweighting() {
  delete m_prw;
  m_prw = 0;
}

unsigned PileupReweighting::random_run_number(const SusyBuffer& buf) const {
  m_prw->SetRandomSeed(
    314159 + buf.mc_channel_number * 2718 + buf.EventNumber);
  return m_prw->GetRandomRunNumber(buf.RunNumber);
}

float PileupReweighting::get_pileup_weight(const SusyBuffer& buffer) {
  float avx = buffer.averageIntPerXing;
  // taken from the same twiki as above
  if (buffer.lbn==1 && int(avx+0.5)==1) avx = 0;
  return m_prw->GetCombinedWeight(
    buffer.RunNumber, buffer.mc_channel_number,avx);
}

namespace {
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
