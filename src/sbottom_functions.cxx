#include <vector>
#include "SUSYBuffer.h"
#include "sbottom_functions.hh"
#include "SUSYTools/SUSYObjDef.h"

bool ChfCheck(const std::vector<size_t>& jet_indices, const SusyBuffer &buf, SUSYObjDef &def)
{
  size_t max=jet_indices.size();
  if(jet_indices.size()>=2) max=2;
  bool shouldbecleaned=false;
  for (size_t ii=0;ii<max;ii++)
    {
      size_t ijet = jet_indices[ii];
      TLorentzVector jetTLV = def.GetJetTLV(static_cast<int>(ijet));
      double chf=0;
      chf=(buf.jet_sumPtTrk)->at(ijet)/jetTLV.Pt();
      
      if(jetTLV.Pt()>100000.&&chf<0.02&&std::abs(jetTLV.Eta())<2) shouldbecleaned=true;
      if(jetTLV.Pt()>100000.&&chf<0.05&&std::abs(jetTLV.Eta())<2
         &&(buf.jet_emfrac->at(ijet)>0.9)) shouldbecleaned=true;
    }
  return(shouldbecleaned);
}