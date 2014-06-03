#include "PileupReweighting.hh"
#include "SusyBuffer.h"
#include "PileupReweighting/TPileupReweighting.h"

PileupReweighting::PileupReweighting(const std::string& pu_config,
				     const std::string& pu_lumicalc):
  m_prw(0)
{
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
  return m_prw->GetCombinedWeight(
    buffer.RunNumber, buffer.mc_channel_number,avx);
}

