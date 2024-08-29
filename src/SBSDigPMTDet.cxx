#include "SBSDigPMTDet.h"
#include "TMath.h"

using namespace std;

SBSDigPMTDet::SBSDigPMTDet()
{
}

SBSDigPMTDet::SBSDigPMTDet(UShort_t uniqueid, UInt_t nchan, std::vector<double> NpeChargeConv):
  fUniqueID(uniqueid), fNChan(nchan)
{
  //for(int i = 0; i<fNChan; i++)PMTmap[i] = PMTSignal();
  for(int i = 0; i<fNChan; i++)PMTmap.push_back(PMTSignal(NpeChargeConv[i]));
}

SBSDigPMTDet::SBSDigPMTDet(UShort_t uniqueid, UInt_t nchan, std::vector<double> NpeChargeConv, double sigmapulse, double gatewidth):
  fUniqueID(uniqueid), fNChan(nchan)
{
  //for(int i = 0; i<fNChan; i++)PMTmap[i] = PMTSignal(NpeChargeConv);
  for(int i = 0; i<fNChan; i++)PMTmap.push_back(PMTSignal(NpeChargeConv[i]));
  fRefPulse = new SPEModel(fUniqueID, sigmapulse, 0, -gatewidth/2., gatewidth/2.);
}

SBSDigPMTDet::~SBSDigPMTDet()
{
  
}

void SBSDigPMTDet::Digitize(g4sbs_tree* T, TRandom3* R)
{
  for(int i = 0; i<fNChan; i++)PMTmap[i].Digitize(i, fUniqueID, T, R, fPedestal, fPedSigma, fADCconv, fADCbits, fTDCconv, fTDCbits, int(fThreshold));
}
  
void SBSDigPMTDet::SetSamples(double sampsize)
{
  for(int i = 0; i<fNChan; i++)PMTmap[i].SetSamples(-fGateWidth/2+30.0, fGateWidth/2+30.0, sampsize);
}

void SBSDigPMTDet::Clear(bool dosamples)
{
  for(int i = 0; i<fNChan; i++)PMTmap[i].Clear(dosamples);
}
