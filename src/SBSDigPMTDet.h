#ifndef SBSDIGPMTDET_H
#define SBSDIGPMTDET_H

#include "SBSDigPMTSignal.h"
#include "g4sbs_tree.h"
#include <iostream>
#include <map>
#include <vector>

//________________________________
class SBSDigPMTDet {
public:
  SBSDigPMTDet();
  SBSDigPMTDet(UShort_t uniqueid, UInt_t nchan, std::vector<double> NpeChargeConv);
  SBSDigPMTDet(UShort_t uniqueid, UInt_t nchan, std::vector<double> NpeChargeConv, double sigmapulse, double gatewidth);
  virtual ~SBSDigPMTDet();
  void Digitize(g4sbs_tree *T, TRandom3 *R);
  void Clear(bool dosamples = false);
  void SetSamples(double sampsize);

  // private:
  UInt_t fNChan;
  UShort_t fUniqueID;
  // std::map<int, PMTSignal> PMTmap;
  std::vector<PMTSignal> PMTmap;
  SPEModel *fRefPulse;

  std::vector<double> fGain;
  double fPedestal;
  double fPedSigma;
  double fTrigOffset;
  double fThreshold;
  double fGateWidth;
  double fADCconv;
  double fADCbits;
  double fTDCconv;
  double fTDCbits;
  double fSigmaPulse;

  std::vector<Double_t> fTimeOffset = std::vector<Double_t>(fNChan, 0.0);
};

#endif
