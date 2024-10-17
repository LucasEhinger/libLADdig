#ifndef SBSDIGPMTSIGNAL_H
#define SBSDIGPMTSIGNAL_H

#include <iostream>
#include <vector>
#include <map>
#include <TROOT.h>
#include "TF1.h"
#include "TH1D.h"
#include "TRandom3.h"
//#include "gmn_tree.h"
#include "g4sbs_tree.h"

//
// classes for signal digitization
//
//________________________________
class SPEModel {
 public:
  SPEModel();
  SPEModel(UShort_t uniqueid, double sigma, double t0 = 0, double tmin = -50, double tmax = +50);
  virtual ~SPEModel();
  double Eval(double t){return fPulseHisto->Interpolate(t);};
  double Integral(int binmin, int binmax);
  bool   PulseOverThr(double charge, double thr);
  bool   FindLeadTrailTime(double charge, double thr, double &t_lead, double &t_trail);
  bool   FindPeakTimeAmp(double charge, double thr, double &amp_peak, double &t_peak);
  
 private:
  TH1D *fPulseHisto;//At least we'll have to setup and use one per detector
  // Eric: I start to suspect that even thoough we only initialize them once a run, 
  // using histograms slow us down... will have to investigate better.
  // NB: it is better than functions, but still...
  double GetHistoX(double y, double x1, double x2);
  
  //ClassDef(SPEModel,1);//it actually doesn't like classdef... (not sure why...)
};

//_________________________________
class PMTSignal {
 public:
  PMTSignal();
  PMTSignal(double npechargeconv);
  void Fill(SPEModel *model, int npe, double thr, double evttime, int signal);
  void Fill_FADCmode1(int npe, double thr, double evttime, double sigmatime, int signal);
  void Fill_FADCmode7(SPEModel *model, int npe, double thr, double evttime, int signal);
  void Digitize(int chan, int detid, g4sbs_tree* T, //gmn_tree* T, 
		TRandom3* R, double ped, double ped_noise, double ADCconv, double ADCbits, double TDCconv, double TDCbits, int thr_adc);
  void Clear(bool dosamples = false);
  ~PMTSignal(){Clear();};
  
  void AddSumEdep(double edep){
    fSumEdep+= edep;
  };
  void SetNpeChargeConv(double npechargeconv){fNpeChargeConv = npechargeconv;};
  void SetSamples(double tmin, double tmax, double sampsize);
    
  double SumEdep(){return fSumEdep;};
  UInt_t Npe(){return fNpe;};
  double Charge(){return fNpe*fNpeChargeConv;};
  UInt_t ADC(){return fADC;};

  double EventTime(){return fEventTime;};
  UInt_t LeadTimesSize(){return fLeadTimes.size();};
  double LeadTime(int i){return fLeadTimes.at(i);};
  UInt_t TrailTimesSize(){return fTrailTimes.size();};
  double TrailTime(int i){return fTrailTimes.at(i);};
  UInt_t TDCSize(){return fTDC_l.size();};
  UInt_t TDC_l(int i){return fTDC_l.at(i);};
  UInt_t TDC_t(int i){return fTDC_t.at(i);};
  //SimEncoder::tdc_data TDCData() { return fTDCData; }

  UInt_t ADCSamples(int i){return fADCSamples[i];};
  
  double Eval(double t){
    //printf("tau = %f, result = %f \n", ftau, TMath::Max(0., fNorm*((t-ft0+ftau*0.4)/(ftau*ftau*0.16))*TMath::Exp(-(t-ft0+ftau*0.4)/(ftau*0.4))) );
    return( TMath::Max(0., fNorm*((t-ft0+ftau*0.4)/(ftau*ftau*0.16))*TMath::Exp(-(t-ft0+ftau*0.4)/(ftau*0.4))) );
  }//can't be worse than TF1::Eval... can it?
  
  void SetPulseParam(double norm, double t0, double tau){
    fNorm = norm;
    ft0 = t0;
    ftau = tau;
  }
  
 private:
  //summing variables for dig...
  double fSumEdep;//Not forced to use it for everything
  UInt_t fNpe;
  double fNpeChargeConv;
  double fADC;// One unique ADC value ?
  std::vector<double> fADC_time;
  double time_offset;

  double fEventTime;
  //TDCs: multiple values possible.
  std::vector<double> fLeadTimes;
  std::vector<double> fTrailTimes;
  std::vector<double> fTDC_l;
  std::vector<double> fTDC_t;
  //SimEncoder::tdc_data fTDCData;
  //TRndmManager* fRN;
  std::vector<double> fPeakAmps;
  //std::vector<Int_t> fAmps;

  //let's try something for HCal
  //TF1* f1;
  //let's try something else for HCal 11
  double fNorm;
  double ft0;
  double ftau;
  int fNADCSamps;
  int fNSamps;
  double fSampSize;
  double fADCSampSize;
  double fTmin;
  double* fSamples;
  double* fADCSamples;
};

#endif
