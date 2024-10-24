#include "SBSDigPMTSignal.h"
#include "TFormula.h"
#include "TMath.h"
#include "g4sbs_types.h"

using namespace std;

//
// Class SPEModel
//
SPEModel::SPEModel() { fPulseHisto = new TH1D("fPulseHisto", "", 1000, -50, 50); }

SPEModel::SPEModel(UShort_t uniqueid, double sigma, double t0, double tmin, double tmax) {
  // TF1 fFunc("fFunc", "landaun", tmin, tmax);//garbage (sorry)
  //  power law x exp decay...
  TF1 fFunc("fFunc", "TMath::Max(0., [0]*((x-[1]+[2]*0.4)/([2]*[2]*0.16))*TMath::Exp(-(x-[1]+[2]*0.4)/([2]*0.4)) )",
            tmin, tmax);

  fFunc.SetParameters(1., t0, sigma);
  const int NbinsTotal = int(tmax - tmin) * 20; // 20 bins/ns should do... since we will extrapolate after...
  // let's try 20...
  fPulseHisto = new TH1D(Form("fPulseHisto_%d", uniqueid), "", NbinsTotal, tmin, tmax);
  double t_i;  //, t_j;
  double ps_i; //, g_j;
  for (int i = 1; i <= NbinsTotal; i++) {
    t_i  = fPulseHisto->GetBinCenter(i + 1);
    ps_i = fFunc.Eval(t_i);
    fPulseHisto->Fill(t_i, ps_i);
    // cout << fPulseHisto->GetBinContent(i) << " ";
  }
  // cout << endl;
}

SPEModel::~SPEModel() { fPulseHisto->Delete(); }

double SPEModel::Integral(int binmin, int binmax) {
  binmax = max(binmax, 0);
  return fPulseHisto->Integral(binmin, binmax, "width");
}

bool SPEModel::PulseOverThr(double charge, double thr) {
  if (fPulseHisto->GetMaximum() < thr / charge) {
    return false;
  } else {
    return true;
  }
};

bool SPEModel::FindLeadTrailTime(double charge, double thr, double &t_lead, double &t_trail) {
  // cout << charge << " " << thr << endl;
  if (!PulseOverThr(charge, thr)) {
    t_lead  = 1.0e38;
    t_trail = 1.0e38;
    return false;
  } else {
    double xmax = fPulseHisto->GetBinCenter(fPulseHisto->GetMaximumBin());
    // cout << fPulseHisto->GetBinContent(1)<< " " << thr/charge << " " <<
    // fPulseHisto->GetBinContent(fPulseHisto->GetNbinsX())<< endl;
    if (fPulseHisto->GetBinContent(1) < thr / charge) {
      t_lead = GetHistoX(thr / charge, fPulseHisto->GetBinLowEdge(1), xmax);
    } else {
      t_lead = 1.0e38;
      return false;
    }
    if (fPulseHisto->GetBinContent(fPulseHisto->GetNbinsX()) < thr / charge) {
      t_trail = GetHistoX(thr / charge, xmax, fPulseHisto->GetBinLowEdge(fPulseHisto->GetNbinsX() + 1));
    } else {
      t_trail = 1.0e38;
      return false;
    }
    // why does this function seem to oblitarate fMCHit containers when t_trail is calculated to be 1.e38... GetHistoX?
    if (t_trail > 1.e30)
      cout << "thr/charge" << thr / charge << " >? " << fPulseHisto->GetBinContent(fPulseHisto->GetNbinsX()) << " or "
           << fPulseHisto->GetBinContent(fPulseHisto->GetNbinsX() + 1) << endl;
    return true;
  }
}

bool SPEModel::FindPeakTimeAmp(double charge, double thr, double &amp_peak, double &t_peak) {
  if (!PulseOverThr(charge, thr)) {
    t_peak   = 1.0e38;
    amp_peak = 1.0e38;
    return false;
  } else {
#include <cmath>

    double result = std::pow(10, 123);
    amp_peak      = fPulseHisto->GetMaximum() * charge * (std::pow(10, 12));
    // Really necessary??? the time peak is at zero by definition of the reference histogram...
    /*
    int binmax = fPulseHisto->GetMaximumBin();
    t_peak = fPulseHisto->GetBinCenter(binmax-1)*fPulseHisto->GetBinContent(binmax-1)+
      amp_peak*fPulseHisto->GetBinCenter(binmax)+
      fPulseHisto->GetBinCenter(binmax+1)*fPulseHisto->GetBinContent(binmax+1);
    t_peak/= (fPulseHisto->GetBinContent(binmax-1)+
              amp_peak+
              fPulseHisto->GetBinContent(binmax+1));
    */
    t_peak = 0.0; //...
    return true;
  }
}

double SPEModel::GetHistoX(double y, double x1, double x2) {
  double splineslope;
  for (int k = fPulseHisto->FindBin(x1); k <= fPulseHisto->FindBin(x2); k++) {
    if (((fPulseHisto->GetBinContent(k + 1) - y) * (fPulseHisto->GetBinContent(k) - y)) < 0) {
      // threshold crossed if diff(TH1::GetBinContent-y) changes sign
      splineslope = (fPulseHisto->GetBinContent(k + 1) - fPulseHisto->GetBinContent(k)) /
                    (fPulseHisto->GetBinCenter(k + 1) - fPulseHisto->GetBinCenter(k));

      return fPulseHisto->GetBinCenter(k) + (y - fPulseHisto->GetBinContent(k)) / splineslope;
    }
  }
  return 1.0e38;
}

//
// Class PMTSignal
//
PMTSignal::PMTSignal()
    : fSumEdep(0), fNpe(0), fNpeChargeConv(1.0), fADC(0), fEventTime(0), fNorm(0), ft0(0), ftau(0), fNADCSamps(0),
      fNSamps(0), time_offset(100) {
  fLeadTimes.clear();
  fTrailTimes.clear();
  fTDC_l.clear();
  fTDC_t.clear();

  fPeakAmps.clear();
  // f1 = 0;
  // R = TRndmManager::GetInstance();
}

PMTSignal::PMTSignal(double npechargeconv)
    : fSumEdep(0), fNpe(0), fNpeChargeConv(npechargeconv), fADC(0), fEventTime(0), fNorm(0), ft0(0), ftau(0),
      fNADCSamps(0), fNSamps(0), time_offset(100) {
  fLeadTimes.clear();
  fTrailTimes.clear();
  fTDC_l.clear();
  fTDC_t.clear();

  fPeakAmps.clear();
  // f1 = 0;
}

void PMTSignal::Fill(SPEModel *model, int npe, double thr, double evttime, int signal) {
  if (signal == 0)
    fEventTime = evttime;
  fNpe += npe;
  // if(model->PulseOverThr(fCharge, thr))fNpe//fADC = model->GetCharge()*model->GetADCconversion();

  // cout << evttime << endl;

  // determine lead and trail times
  double t_lead, t_trail;
  // find the lead and trail time for *this* pulse, not the total pulse
  // the following line slows the digitization in the case of full background.
  // Needs improvement ASAP!
  bool goodtime = model->FindLeadTrailTime(npe * fNpeChargeConv, thr, t_lead, t_trail);

  // if(goodtime)cout << evttime << " " << t_lead << " " << t_trail << endl;

  if (goodtime) {
    // TODO: Isn't comparable with multiple peaks
    double amp_peak, t_peak;

    if (model->FindPeakTimeAmp(npe * fNpeChargeConv, thr, amp_peak, t_peak)) {
      fPeakAmps.push_back(amp_peak);
    }
  }

  t_lead += fEventTime;
  t_trail += fEventTime;

  if (goodtime) {
    // Filter here the lead and trail times
    if (fLeadTimes.size() > 0) {
      // Check if the lead and trail times are inside an existing lead time- trail time pair
      // we assume here that fLeadTimes and fTrailTimes are same size
      // *if we do things correctly, that should be the case*
      // we shall keep lead-trail times pair in timing order
      // we neglect pulse overlaps ftm.
      if (fLeadTimes.size() != fTrailTimes.size()) {
        cout << " B - Warning: size of lead times container: " << fLeadTimes.size()
             << " != size of trail times container: " << fTrailTimes.size() << endl;
      }

      for (size_t i = 0; i < fLeadTimes.size(); i++) {
        // possibility of the current pair straddling with others.... :/
        // treat those separately to simplify...
        // tL < tT_i-1 < tL_i < tT
        if (i > 0) {
          if (t_lead < fTrailTimes.at(i - 1) && fLeadTimes.at(i) < t_trail) {
            // do necessary substitutions first, then erase
            //  tL < tL_i-1 => tL *replaces* tL_i-1
            if (t_lead < fLeadTimes.at(i - 1)) {
              fLeadTimes.erase(fLeadTimes.begin() + i - 1);
              fLeadTimes.insert(fLeadTimes.begin() + i - 1, t_lead);
            }
            // tT_i < tT => tT *replaces* tT_i
            if (fTrailTimes.at(i) < t_trail) {
              fTrailTimes.erase(fTrailTimes.begin() + i);
              fTrailTimes.insert(fTrailTimes.begin() + i, t_trail);
            }
            fLeadTimes.erase(fLeadTimes.begin() + i);
            fTrailTimes.erase(fTrailTimes.begin() + i - 1);
            break;
          }
        } // end if(i>0)
        // tL < tT_i < tL_i+1 < tT
        if (i < fLeadTimes.size() - 1) {
          if (t_lead < fTrailTimes.at(i) && fLeadTimes.at(i + 1) < t_trail) {
            // do necessary substitutions first, then erase
            //  tL < tL_i => tL *replaces* tL_i
            if (t_lead < fLeadTimes.at(i)) {
              fLeadTimes.erase(fLeadTimes.begin() + i);
              fLeadTimes.insert(fLeadTimes.begin() + i, t_lead);
            }
            // tT_i+1 < tT => tT *replaces* tT_i+1
            if (fTrailTimes.at(i + 1) < t_trail) {
              fTrailTimes.erase(fTrailTimes.begin() + i + 1);
              fTrailTimes.insert(fTrailTimes.begin() + i + 1, t_trail);
            }
            fLeadTimes.erase(fLeadTimes.begin() + i + 1);
            fTrailTimes.erase(fTrailTimes.begin() + i);
            break;
          }
        } // end if(i<fLeadTimes.size()-1)

        // if not, 6 cases to consider:
        // tL < tT < tL_i < tT_i => both inserted *before* existing pair
        if (t_trail < fLeadTimes.at(i)) {
          fLeadTimes.insert(fLeadTimes.begin() + i, t_lead);
          fTrailTimes.insert(fTrailTimes.begin() + i, t_trail);
          break;
        }
        // tL < tL_i < tT < tT_i => tL *replaces* tL_i
        if (t_lead < fLeadTimes.at(i) && fLeadTimes.at(i) < t_trail && t_trail < fTrailTimes.at(i)) {
          fLeadTimes.erase(fLeadTimes.begin() + i);
          fLeadTimes.insert(fLeadTimes.begin() + i, t_lead);
          break;
        }
        // tL_i < tL < tT < tT_i => tL *replaces* tL_i AND tT *replaces* tT_i
        if (t_lead < fLeadTimes.at(i) && fTrailTimes.at(i) < t_trail) {
          fLeadTimes.erase(fLeadTimes.begin() + i);
          fLeadTimes.insert(fLeadTimes.begin() + i, t_lead);
          fTrailTimes.erase(fTrailTimes.begin() + i);
          fTrailTimes.insert(fTrailTimes.begin() + i, t_trail);
          break;
        }
        if (fLeadTimes.at(i) < t_lead && t_trail < fTrailTimes.at(i)) {
          break;
        }
        if (fLeadTimes.at(i) < t_lead && t_lead < fTrailTimes.at(i) && fTrailTimes.at(i) < t_trail) {
          fTrailTimes.erase(fTrailTimes.begin() + i);
          fTrailTimes.insert(fTrailTimes.begin() + i, t_trail);
        }
        if (fTrailTimes.at(i) < t_lead) {
          fLeadTimes.insert(fLeadTimes.begin() + i + 1, t_lead);
          fTrailTimes.insert(fTrailTimes.begin() + i + 1, t_trail);
          break;
        }
        if (fLeadTimes.size() != fTrailTimes.size()) {
          cout << " A - Warning: size of lead times container: " << fLeadTimes.size()
               << " != size of trail times container: " << fTrailTimes.size() << endl;
        }
        if (i >= fLeadTimes.size())
          cout << "Warning: i = " << i << " >= size of containers = " << fLeadTimes.size() << endl;
      }
    } else {
      // of course, if initial size was 0, just psuh it back
      // hopefully it will be the case most of the time
      fLeadTimes.push_back(t_lead);
      fTrailTimes.push_back(t_trail);
    }
    if (fLeadTimes.size() != fTrailTimes.size()) {
      cout << " A - Warning: size of lead times container: " << fLeadTimes.size()
           << " != size of trail times container: " << fTrailTimes.size() << endl;
    }
  } // end if(t_lead && t_trail<30)
}

void PMTSignal::Fill_FADCmode1(int npe, double thr, double evttime, double sigmatime, int signal) {
  if (signal == 0)
    fEventTime = evttime;
  fNpe += npe;
  // if(model->PulseOverThr(fCharge, thr))fNpe//fADC = model->GetCharge()*model->GetADCconversion();

  // cout << evttime << " <? " << fTmin << "+" << fNSamps << "*" << fSampSize << " = " << fTmin+fNSamps*fSampSize <<
  // endl;

  if (evttime >= fTmin + fNSamps * fSampSize)
    return;

  // cout << "fillfadcmode1 " << sigmatime << endl;

  SetPulseParam(npe * fNpeChargeConv, evttime, sigmatime);

  // cout << npe*fNpeChargeConv << " " << evttime << " " << sigmatime << endl;

  // f1->SetParameters(npe*fNpeChargeConv, evttime, sigmatime);
  // determine lead and trail times
  double t_lead, t_trail;
  bool goodtime = false; // model->FindLeadTrailTime(npe*fNpeChargeConv, thr, t_lead, t_trail);
  // the following block slows the digitization in the case of full background.
  // Needs improvement ASAP!
  if (fNSamps) {
    fSamples[0] += Eval(fTmin + (0.5) * fSampSize); // f1->Eval(fTmin+(0.5)*fSampSize);//*fSampSize;
    // cout << Eval(fTmin+(0.5)*fSampSize) << endl;
    // Evaluate this function might be a bit of a time drain!
    // cout << fSamples[0] << " ";
    for (int i = 1; i < fNSamps; i++) {
      fSamples[i] += Eval(fTmin + (i + 0.5) * fSampSize); // f1->Eval(fTmin+(i+0.5)*fSampSize);//*fSampSize;
      // if(i>0){
      // cout << Eval(fTmin+(i+0.5)*fSampSize) << endl;
      // cout << fSamples[i] << " ";
      if (fSamples[i - 1] <= thr && thr < fSamples[i]) {
        t_lead = fTmin + (i - 0.5) * fSampSize + fSampSize * (thr - fSamples[i - 1]) / (fSamples[i] - fSamples[i - 1]);
        goodtime = true;
      }
      if (fSamples[i - 1] >= thr && thr > fSamples[i]) {
        t_trail = fTmin + (i - 0.5) * fSampSize + fSampSize * (thr - fSamples[i - 1]) / (fSamples[i] - fSamples[i - 1]);
        goodtime = true;
      }
    }
    // cout << endl;
  }

  // find the lead and trail time for *this* pulse, not the total pulse
  // t_lead+=evttime;
  // t_trail+=evttime;
  // if(goodtime)cout << evttime << " " << t_lead << " " << t_trail << endl;

  if (goodtime && signal == 0) {
    // Filter here the lead and trail times
    if (fLeadTimes.size() > 0) {
      // Check if the lead and trail times are inside an existing lead time- trail time pair
      // we assume here that fLeadTimes and fTrailTimes are same size
      // *if we do things correctly, that should be the case*
      // we shall keep lead-trail times pair in timing order
      // we neglect pulse overlaps ftm.
      if (fLeadTimes.size() != fTrailTimes.size()) {
        cout << " B - Warning: size of lead times container: " << fLeadTimes.size()
             << " != size of trail times container: " << fTrailTimes.size() << endl;
      }

      for (size_t i = 0; i < fLeadTimes.size(); i++) {
        // possibility of the current pair straddling with others.... :/
        // treat those separately to simplify...
        // tL < tT_i-1 < tL_i < tT
        if (i > 0) {
          if (t_lead < fTrailTimes.at(i - 1) && fLeadTimes.at(i) < t_trail) {
            // do necessary substitutions first, then erase
            //  tL < tL_i-1 => tL *replaces* tL_i-1
            if (t_lead < fLeadTimes.at(i - 1)) {
              fLeadTimes.erase(fLeadTimes.begin() + i - 1);
              fLeadTimes.insert(fLeadTimes.begin() + i - 1, t_lead);
            }
            // tT_i < tT => tT *replaces* tT_i
            if (fTrailTimes.at(i) < t_trail) {
              fTrailTimes.erase(fTrailTimes.begin() + i);
              fTrailTimes.insert(fTrailTimes.begin() + i, t_trail);
            }
            fLeadTimes.erase(fLeadTimes.begin() + i);
            fTrailTimes.erase(fTrailTimes.begin() + i - 1);
            break;
          }
        } // end if(i>0)
        // tL < tT_i < tL_i+1 < tT
        if (i < fLeadTimes.size() - 1) {
          if (t_lead < fTrailTimes.at(i) && fLeadTimes.at(i + 1) < t_trail) {
            // do necessary substitutions first, then erase
            //  tL < tL_i => tL *replaces* tL_i
            if (t_lead < fLeadTimes.at(i)) {
              fLeadTimes.erase(fLeadTimes.begin() + i);
              fLeadTimes.insert(fLeadTimes.begin() + i, t_lead);
            }
            // tT_i+1 < tT => tT *replaces* tT_i+1
            if (fTrailTimes.at(i + 1) < t_trail) {
              fTrailTimes.erase(fTrailTimes.begin() + i + 1);
              fTrailTimes.insert(fTrailTimes.begin() + i + 1, t_trail);
            }
            fLeadTimes.erase(fLeadTimes.begin() + i + 1);
            fTrailTimes.erase(fTrailTimes.begin() + i);
            break;
          }
        } // end if(i<fLeadTimes.size()-1)

        // if not, 6 cases to consider:
        // tL < tT < tL_i < tT_i => both inserted *before* existing pair
        if (t_trail < fLeadTimes.at(i)) {
          fLeadTimes.insert(fLeadTimes.begin() + i, t_lead);
          fTrailTimes.insert(fTrailTimes.begin() + i, t_trail);
          break;
        }
        // tL < tL_i < tT < tT_i => tL *replaces* tL_i
        if (t_lead < fLeadTimes.at(i) && fLeadTimes.at(i) < t_trail && t_trail < fTrailTimes.at(i)) {
          fLeadTimes.erase(fLeadTimes.begin() + i);
          fLeadTimes.insert(fLeadTimes.begin() + i, t_lead);
          break;
        }
        // tL_i < tL < tT < tT_i => tL *replaces* tL_i AND tT *replaces* tT_i
        if (t_lead < fLeadTimes.at(i) && fTrailTimes.at(i) < t_trail) {
          fLeadTimes.erase(fLeadTimes.begin() + i);
          fLeadTimes.insert(fLeadTimes.begin() + i, t_lead);
          fTrailTimes.erase(fTrailTimes.begin() + i);
          fTrailTimes.insert(fTrailTimes.begin() + i, t_trail);
          break;
        }
        if (fLeadTimes.at(i) < t_lead && t_trail < fTrailTimes.at(i)) {
          break;
        }
        if (fLeadTimes.at(i) < t_lead && t_lead < fTrailTimes.at(i) && fTrailTimes.at(i) < t_trail) {
          fTrailTimes.erase(fTrailTimes.begin() + i);
          fTrailTimes.insert(fTrailTimes.begin() + i, t_trail);
        }
        if (fTrailTimes.at(i) < t_lead) {
          fLeadTimes.insert(fLeadTimes.begin() + i + 1, t_lead);
          fTrailTimes.insert(fTrailTimes.begin() + i + 1, t_trail);
          break;
        }
        if (fLeadTimes.size() != fTrailTimes.size()) {
          cout << " A - Warning: size of lead times container: " << fLeadTimes.size()
               << " != size of trail times container: " << fTrailTimes.size() << endl;
        }
        if (i >= fLeadTimes.size())
          cout << "Warning: i = " << i << " >= size of containers = " << fLeadTimes.size() << endl;
      }
    } else {
      // of course, if initial size was 0, just psuh it back
      // hopefully it will be the case most of the time
      fLeadTimes.push_back(t_lead);
      fTrailTimes.push_back(t_trail);
    }
    if (fLeadTimes.size() != fTrailTimes.size()) {
      cout << " A - Warning: size of lead times container: " << fLeadTimes.size()
           << " != size of trail times container: " << fTrailTimes.size() << endl;
    }
  } // end if(t_lead && t_trail<30)
}

void PMTSignal::Fill_FADCmode7(SPEModel *model, int npe, double thr, double evttime, int signal) {
  if (signal == 0)
    fEventTime = evttime;
  fNpe += npe;

  // Mode7 (?): pedestal, amplitude, integral, peak time (I assume there is a threshold?)
  int evttime_offset = int(evttime / 80.);
  for (int k = 0; k < fNADCSamps; k++) {
    // if evttime = 0, no offset
    //  fPulseHisto bears 20 bins/ns, each FADC sample is 4ns
    //  => 1 fADC sample = 80 ns
    fADCSamples[k] =
        0; // model->Integral(evttime_offset+k*80,
           // evttime_offset+(k+1)*80-1);//fPulseHisto->Integral(evttime_offset+k*80, evttime_offset+(k+1)*80, "width");
  }

  double amp_peak, t_peak;

  if (model->FindPeakTimeAmp(npe * fNpeChargeConv, thr, amp_peak, t_peak)) {
    fLeadTimes.push_back(t_peak + evttime);
    fPeakAmps.push_back(amp_peak);
  }
}

void PMTSignal::Digitize(int chan, int detid, g4sbs_tree *T, // gmn_tree* T,
                         TRandom3 *R, double ped, double ped_noise, double ADCconv, double ADCbits, double TDCconv,
                         double TDCbits, int thr_adc) {
  if (fNpe <= 0) {
    // fADC = R->Gaus(ped, ped_noise);
    return;
  }

  fADC = (Charge() * 1.0e15 / ADCconv + R->Gaus(ped, ped_noise)); // No longer casting to Int

  if (fADC > UInt_t((TMath::Power(2, ADCbits)))) { // No longer casting to Int
    fADC = (TMath::Power(2, ADCbits));             // No longer casting to Int
  }

  double tdc_value;
  if (fLeadTimes.size()) {
    for (size_t i = 0; i < fLeadTimes.size(); i++) {
      // cout << "detid " << detid << " fLeadTimes.at(" << i << ") " << fLeadTimes.at(i) << " fTrailTimes.at(" << i <<
      // ") " << fTrailTimes.at(i) << endl;
      //  trim "all" bits that are above the number of TDC bits - a couple to speed it up
      //  (since TDC have a revolving clock, as far as I understand)
      //  let's use an arbitrary reference time offset of 1000 (?) TDC chans before the trigger
      tdc_value = ((fLeadTimes.at(i)) / TDCconv) + time_offset/TDCconv; // No longer casting to Int
      fTDC_l.push_back(tdc_value);                              // they're already sorted in order, presumably
      // also mark the traling time with setting bin 31 to 1 // need to reconvert then
      if (fTrailTimes.size()) {
        tdc_value = ((fTrailTimes.at(i)) / TDCconv) + time_offset/TDCconv; // No longer casting to Int
        fTDC_t.push_back(tdc_value);
      }
    }
  }
  // cout << "detid " << detid << " TDC size " << fTDCs.size() << endl;

  int i_tc     = -1;
  int i_max    = -1;
  double vpeak = thr_adc; // to ensure we only look for the max if the threshold is crossed
  double vmin  = 0;

  if (fNSamps) {
    fADC        = 0;
    Int_t Nconv = fNSamps / fNADCSamps;
    for (int i = 0; i < fNADCSamps; i++) {
      for (int j = 0; j < Nconv; j++) {
        fADCSamples[i] += fSamples[i * Nconv + j] * fSampSize; // renormalize the sample for the integration;
        // if(detid==ACTIVEANA_UNIQUE_DETID && TMath::IsNaN(fADCSamples[i]) )cout << " i " << i << " j " << j << ", idx
        // " << i*Nconv+j << " size " << sizeof(fSamples) << " samp " << fSamples[i*Nconv+j] << "; ";
      }
      fADCSamples[i] *= 1.0e15 / ADCconv;
      fADCSamples[i] += R->Gaus(ped, ped_noise);

      if (fADCSamples[i] > 4095)
        fADCSamples[i] = 4095;

      fADC += fADCSamples[i];
      // /*
      // if(i_tc<0 && fADCSamples[i]>thr_adc)i_tc = i;
      // if(fADCSamples[i]>vpeak){
      //   vpeak = fADCSamples[i];
      //   i_max = i;
      // }
      // if(i<4)vmin+= fADCSamples[i]/4;
      // */
    }
    // if(detid==ACTIVEANA_UNIQUE_DETID)cout << endl;
  }
  double peak_amp = 0.0;
  if (!fPeakAmps.empty()) {
    peak_amp = *std::max_element(fPeakAmps.begin(), fPeakAmps.end());
  }
  // fSumEdep*=1.0e9;// store in eV.

  // Fill in directly (hoping it takes less time...)

  if (detid == HODO_UNIQUE_DETID) {
    // T->Earm_BBHodo_dighit_nchan++;
    // T->Earm_BBHodo_dighit_chan->push_back(chan);
    // T->Earm_BBHodo_dighit_adc->push_back(fADC);
    T->Earm_BBHodo_Dig.nchan++;
    T->Earm_BBHodo_Dig.chan->push_back(chan);
    T->Earm_BBHodo_Dig.adc->push_back(fADC);
    T->Earm_BBHodo_Dig.amp->push_back(peak_amp);
    // T->Earm_BBHodo_Dig.amp->push_back(850);
    T->Earm_BBHodo_Dig.adc_time->push_back(fEventTime + time_offset);
    if (fTDC_l.size()) {
      for (int j = 0; j < fTDC_l.size(); j++) {
        T->Earm_BBHodo_Dig.tdc_t->push_back(fTDC_t[j]);
        T->Earm_BBHodo_Dig.tdc_l->push_back(fTDC_l[j]);
        // T->Earm_BBHodo_Dig.tdc_t->push_back(835);// temp hack. fixme
        // T->Earm_BBHodo_Dig.tdc_l->push_back(830);
      }
      // equalize the hits:
      int max_size = max(T->Earm_BBHodo_Dig.tdc_l->size(), T->Earm_BBHodo_Dig.tdc_t->size());
      max_size     = max(max_size, T->Earm_BBHodo_Dig.nchan);
      while (T->Earm_BBHodo_Dig.nchan < max_size) {
        T->Earm_BBHodo_Dig.nchan++;
        T->Earm_BBHodo_Dig.chan->push_back(chan);
        T->Earm_BBHodo_Dig.adc->push_back(-1000000);
        T->Earm_BBHodo_Dig.amp->push_back(-1000000);
      }
      while (T->Earm_BBHodo_Dig.tdc_l->size() < max_size) {
        T->Earm_BBHodo_Dig.tdc_l->push_back(-1000000);
      }
      while (T->Earm_BBHodo_Dig.tdc_t->size() < max_size) {
        T->Earm_BBHodo_Dig.tdc_t->push_back(-1000000);
      }

    } else {
      T->Earm_BBHodo_Dig.tdc_l->push_back(-1000000);
      T->Earm_BBHodo_Dig.tdc_t->push_back(-1000000);
    }
  }
} //

void PMTSignal::SetSamples(double tmin, double tmax, double sampsize) {
  fTmin        = tmin;
  fADCSampSize = sampsize;
  // fSampSize = sampsize/10;//the bin size is too large for the tdc size
  fSampSize   = sampsize / 32; // bin size is similar to TDC bit size
  fNADCSamps  = round((tmax - tmin) / fADCSampSize);
  fNSamps     = round((tmax - tmin) / fSampSize);
  fADCSamples = new double[fNADCSamps];
  fSamples    = new double[fNSamps];

  // cout << fNADCSamps << " " << fNSamps << " " << fADCSampsSize << " " << fSampSize << endl;

  memset(fSamples, 0, fNSamps * sizeof(double));
  memset(fADCSamples, 0, fNADCSamps * sizeof(double));
  // f1 = new TF1("f1", "landaun", tmin, tmax);
  // f1 = new TF1("fFunc",
  //"TMath::Max(0., [0]*((x-[1]+[2]*0.4)/([2]*[2]*0.16))*TMath::Exp(-(x-[1]+[2]*0.4)/([2]*0.4)) )",
  // tmin, tmax);
}

void PMTSignal::Clear(bool dosamples) {
  // cout << " PMTSignal::Clear() " << endl;

  fSumEdep = 0;
  fNpe     = 0;
  fADC     = 0;

  fEventTime = 0;
  fLeadTimes.clear();
  fTrailTimes.clear();
  fTDC_t.clear();
  fTDC_l.clear();

  fPeakAmps.clear();

  if (dosamples) {
    memset(fSamples, 0, fNSamps * sizeof(double));
    memset(fADCSamples, 0, fNADCSamps * sizeof(double));
  }
}

// ClassImp(SPEModel)
