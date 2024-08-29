#include "SBSDigAuxi.h"
#include "TMath.h"
#include "g4sbs_types.h"
#define DEBUG 0

using namespace std;

bool UnfoldData(g4sbs_tree *T, double theta_sbs, double d_hcal, TRandom3 *R, std::vector<SBSDigPMTDet *> pmtdets,
                std::vector<int> detmap, std::vector<SBSDigGEMDet *> gemdets, std::vector<int> gemmap, double tzero,
                int signal) {
  bool has_data = false;

  int Npe;
  double t;

  double x_ref = -d_hcal * sin(theta_sbs);
  double z_ref = d_hcal * cos(theta_sbs);

  double z_hit, Npe_Edep_ratio, sigma_tgen;

  // double z_det, pz, E,
  double beta, sin2thetaC;

  int chan;
  int mod;

  int idet = 0;
  if (!detmap.empty()) {
    // ordering by increasing unique det ID
    while (idet < (int)detmap.size()) {
      if (detmap[idet] != HCAL_UNIQUE_DETID) {
        idet++;
      } else {
        break;
      }
    }
    if (idet >= detmap.size())
      idet = -1;
    // cout << " " << idet;
    if (idet >= 0) { // && T->Harm_HCalScint.sumedep) {
      for (size_t k = 0; k < T->Harm_HCalScint.sumedep->size(); k++) {
        chan = T->Harm_HCalScint.cell->at(k);
        // T->Harm_HCalScint_hit_sumedep->at(k);

        z_hit = -(T->Harm_HCalScint.xhitg->at(k) - x_ref) * sin(theta_sbs) +
                (T->Harm_HCalScint.zhitg->at(k) - z_ref) * cos(theta_sbs);

        // Evaluation of number of photoelectrons from energy deposit documented at:
        // https://sbs.jlab.org/DocDB/0000/000043/002/Harm_HCal_Digi_EdepOnly_2.pdf
        // TODO: put that stuff in DB...
        Npe_Edep_ratio = 5.242 + 11.39 * z_hit + 10.41 * pow(z_hit, 2);
        Npe            = R->Poisson(Npe_Edep_ratio * T->Harm_HCalScint.sumedep->at(k) * 1.0e3);
        t              = tzero + R->Gaus(T->Harm_HCalScint.tavg->at(k) + 10.11, 1.912) - pmtdets[idet]->fTrigOffset;

        sigma_tgen = pmtdets[idet]->fSigmaPulse * (0.4244 + 11380 / pow(Npe + 153.4, 2));
        // cout << sigma_tgen << endl;
        // Generate here,...
        // cout << " HCal : t = " << t << ", t_zero = " << tzero << ", t_avg = " << T->Harm_HCalScint.tavg->at(k) << ",
        // -t_offset = " << -pmtdets[idet]->fTrigOffset << ", Npe " << Npe << endl;
        pmtdets[idet]->PMTmap[chan].Fill_FADCmode1(Npe, pmtdets[idet]->fThreshold, t, sigma_tgen, signal);
      }
      has_data = true;
    }

    idet = 0;
    while (idet < (int)detmap.size()) {
      // if(idet<0)idet++;
      if (detmap[idet] != BBPS_UNIQUE_DETID) {
        idet++;
      } else {
        break;
      }
    }
    if (idet >= detmap.size())
      idet = -1;
    // cout << " " << idet;
    //  Process BB PS data
    if (idet >= 0) { // && T->Earm_BBPSTF1.nhits){
      for (int i = 0; i < T->Earm_BBPSTF1.nhits; i++) {
        // Evaluation of number of photoelectrons and time from energy deposit documented at:
        if (T->Earm_BBPSTF1.sumedep->at(i) < 1.e-4)
          continue;
        // check probability to generate p.e. yield
        bool genpeyield = true;
        if (T->Earm_BBPSTF1.sumedep->at(i) < 1.e-2)
          genpeyield = R->Uniform(0, 1) <= 1 - exp(0.29 - 950. * T->Earm_BBPSTF1.sumedep->at(i));
        // if we're go, let's generate
        if (genpeyield) {
          beta =
              sqrt(pow(m_e + T->Earm_BBPSTF1.sumedep->at(i), 2) - m_e * m_e) / (m_e + T->Earm_BBPSTF1.sumedep->at(i));
          sin2thetaC = TMath::Max(1. - 1. / pow(n_leadglass * beta, 2), 0.);
          // 1500. Used to be 454.: just wrong
          Npe =
              R->Poisson(300.0 * T->Earm_BBPSTF1.sumedep->at(i) * sin2thetaC / (1. - 1. / (n_leadglass * n_leadglass)));
          t = tzero + T->Earm_BBPSTF1.tavg->at(i) +
              R->Gaus(3.2 - 5.805 * T->Earm_BBPSTF1.zhit->at(i) - 17.77 * pow(T->Earm_BBPSTF1.zhit->at(i), 2), 0.5) -
              pmtdets[idet]->fTrigOffset;
          chan = T->Earm_BBPSTF1.cell->at(i);
          // T->Earm_BBPSTF1_hit_sumedep->at(i);

          // if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
          // pmtdets[idet]->PMTmap[chan].Fill(pmtdets[idet]->fRefPulse, Npe, 0, t, signal);
          // pmtdets[idet]->PMTmap[chan].Fill_FADCmode7(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t,
          // signal); cout << " BBPS : t = " << t << ", t_zero = " << tzero << ", t_avg = " <<
          // T->Earm_BBPSTF1.tavg->at(i) << ", -t_offset = " << -pmtdets[idet]->fTrigOffset << ", Npe " << Npe << endl;
          pmtdets[idet]->PMTmap[chan].Fill_FADCmode1(Npe, pmtdets[idet]->fThreshold, t, pmtdets[idet]->fSigmaPulse,
                                                     signal);
        }
      }
      has_data = true;
    }

    idet = 0;
    while (idet < (int)detmap.size()) {
      // if(idet<0)idet++;
      if (detmap[idet] != BBSH_UNIQUE_DETID) {
        idet++;
      } else {
        break;
      }
    }
    if (idet >= detmap.size())
      idet = -1;
    // cout << " " << idet;
    if (idet >= 0) { // && T->Earm_BBSHTF1.nhits){
      for (int i = 0; i < T->Earm_BBSHTF1.nhits; i++) {
        // Evaluation of number of photoelectrons and time from energy deposit documented at:
        //
        if (T->Earm_BBSHTF1.sumedep->at(i) < 1.e-4)
          continue;
        // check probability to generate p.e. yield
        bool genpeyield = true;
        if (T->Earm_BBSHTF1.sumedep->at(i) < 1.e-2)
          genpeyield = R->Uniform(0, 1) <= 1 - exp(0.29 - 950. * T->Earm_BBSHTF1.sumedep->at(i));
        // if we're go, let's generate
        if (genpeyield) {
          beta =
              sqrt(pow(m_e + T->Earm_BBSHTF1.sumedep->at(i), 2) - m_e * m_e) / (m_e + T->Earm_BBSHTF1.sumedep->at(i));
          sin2thetaC = TMath::Max(1. - 1. / pow(n_leadglass * beta, 2), 0.);
          // 1800. Used to be 932.: just wrong
          Npe =
              R->Poisson(360.0 * T->Earm_BBSHTF1.sumedep->at(i) * sin2thetaC / (1. - 1. / (n_leadglass * n_leadglass)));
          t = tzero + T->Earm_BBSHTF1.tavg->at(i) +
              R->Gaus(2.216 - 8.601 * T->Earm_BBSHTF1.zhit->at(i) - 7.469 * pow(T->Earm_BBSHTF1.zhit->at(i), 2), 0.8) -
              pmtdets[idet]->fTrigOffset;
          chan = T->Earm_BBSHTF1.cell->at(i);
          // T->Earm_BBSHTF1_hit_sumedep->at(i);

          // if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
          // pmtdets[idet]->PMTmap[chan].Fill(pmtdets[idet]->fRefPulse, Npe, 0, t, signal);
          // pmtdets[idet]->PMTmap[chan].Fill_FADCmode7(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t,
          // signal); cout << " BBSH : t = " << t << ", t_zero = " << tzero << ", t_avg = " <<
          // T->Earm_BBSHTF1.tavg->at(i) << ", -t_offset = " << -pmtdets[idet]->fTrigOffset << ", Npe " << Npe << endl;
          pmtdets[idet]->PMTmap[chan].Fill_FADCmode1(Npe, pmtdets[idet]->fThreshold, t, pmtdets[idet]->fSigmaPulse,
                                                     signal);
        }
      }
      has_data = true;
    }

    // ECal
    idet = 0;
    while (idet < (int)detmap.size()) {
      // if(idet<0)idet++;
      if (detmap[idet] != ECAL_UNIQUE_DETID) {
        idet++;
      } else {
        break;
      }
    }
    if (idet >= detmap.size())
      idet = -1;
    // cout << " " << idet;
    if (idet >= 0) { // && T->Earm_ECalTF1.nhits){
      for (int i = 0; i < T->Earm_ECalTF1.nhits; i++) {
        // Evaluation of number of photoelectrons and time from energy deposit documented at:
        //
        if (T->Earm_ECalTF1.sumedep->at(i) < 1.e-4)
          continue;
        // check probability to generate p.e. yield
        bool genpeyield = true;
        if (T->Earm_ECalTF1.sumedep->at(i) < 1.e-2)
          genpeyield = R->Uniform(0, 1) <= 1 - exp(0.29 - 950. * T->Earm_ECalTF1.sumedep->at(i));
        // if we're go, let's generate
        if (genpeyield) {
          beta =
              sqrt(pow(m_e + T->Earm_ECalTF1.sumedep->at(i), 2) - m_e * m_e) / (m_e + T->Earm_ECalTF1.sumedep->at(i));
          sin2thetaC = TMath::Max(1. - 1. / pow(n_leadglass * beta, 2), 0.);
          // 1800. Used to be 932.: just wrong
          Npe = R->Poisson(
              360.0 * T->Earm_ECalTF1.sumedep->at(i) * sin2thetaC /
              (1. - 1. / (n_leadglass * n_leadglass))); // ECal is a lead glass detector, so this could be correct
          // rough, first step: t of photons at PMT = tavg_TF1_hit + (0.172-z_TF1_hit)*n/C
          t = tzero + T->Earm_ECalTF1.tavg->at(i) + (0.172 - T->Earm_ECalTF1.zhit->at(i) * n_leadglass / 3.e-1) -
              pmtdets[idet]->fTrigOffset;
          chan = T->Earm_ECalTF1.cell->at(i);
          // T->Earm_ECalTF1_hit_sumedep->at(i);

          // if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
          // pmtdets[idet]->PMTmap[chan].Fill(pmtdets[idet]->fRefPulse, Npe, 0, t, signal);
          // ECal switched to FADC :)
          // cout << T->Earm_ECalTF1.sumedep->at(i) << " " << Npe << " " << t << " " << pmtdets[idet]->fSigmaPulse <<
          // endl;
          pmtdets[idet]->PMTmap[chan].Fill_FADCmode1(Npe, pmtdets[idet]->fThreshold, t, pmtdets[idet]->fSigmaPulse,
                                                     signal);
        }
      }
      has_data = true;
    }

    idet = 0;
    while (idet < (int)detmap.size()) {
      // if(idet<0)idet++;
      if (detmap[idet] != GRINCH_UNIQUE_DETID) {
        idet++;
      } else {
        break;
      }
    }
    if (idet >= detmap.size())
      idet = -1;
    // cout << " " << idet;
    //  Process GRINCH data
    if (idet >= 0) { // && T->Earm_GRINCH.nhits){
      for (int i = 0; i < T->Earm_GRINCH.nhits; i++) {
        chan = int(T->Earm_GRINCH.PMT->at(i) / 5) - 1;
        t    = tzero + T->Earm_GRINCH.Time_avg->at(i) + pmtdets[idet]->fTrigOffset;
        Npe  = T->Earm_GRINCH.NumPhotoelectrons->at(i);

        // if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
        pmtdets[idet]->PMTmap[chan].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, signal);
      }
      has_data = true;
    }

    idet = 0;
    while (idet < (int)detmap.size()) {
      // if(idet<0)idet++;
      if (detmap[idet] != HODO_UNIQUE_DETID) {
        idet++;
      } else {
        break;
      }
    }
    if (idet >= detmap.size())
      idet = -1;
    // cout << " " << idet;
    //  Process hodoscope data
    if (idet >= 0) { // && T->Earm_BBHodoScint.nhits){
      for (int i = 0; i < T->Earm_BBHodoScint.nhits; i++) {
        for (int j = 0; j < 2; j++) { // j = 0: close beam PMT, j = 1: far beam PMT //LHE: Probably the left & right??
          // Evaluation of number of photoelectrons and time from energy deposit documented at:
          // https://hallaweb.jlab.org/dvcslog/SBS/170711_172759/BB_hodoscope_restudy_update_20170711.pdf
          Npe = R->Poisson(1.0e7 * T->Earm_BBHodoScint.sumedep->at(i) * 0.113187 *
                           exp(-(0.3 + pow(-1, j) * T->Earm_BBHodoScint.xhit->at(i)) / 1.03533) * 0.24);
          t   = tzero + T->Earm_BBHodoScint.tavg->at(i) + (0.55 + pow(-1, j) * T->Earm_BBHodoScint.xhit->at(i)) / 0.15 -
              pmtdets[idet]->fTrigOffset;
          chan = T->Earm_BBHodoScint.cell->at(i) * 2 + j;
          // T->Earm_BBHodoScint_hit_sumedep->at(i);
          // if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
          pmtdets[idet]->PMTmap[chan].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, signal);
          // cout << "BBhodo: chan "<< chan << " Npe " << Npe << " t " << t << ", t_zero = " << tzero << ", t_avg = " <<
          // T->Earm_BBHodoScint.tavg->at(i) << ", -t_offset = " << -pmtdets[idet]->fTrigOffset << endl;
        }
      }
      has_data = true;
    }

    // CDet
    idet = 0;
    while (idet < (int)detmap.size()) {
      // if(idet<0)idet++;
      // cout << " " << idet << ":" << detmap[idet];
      if (detmap[idet] != CDET_UNIQUE_DETID) {
        idet++;
      } else {
        break;
      }
    }
    if (idet >= detmap.size())
      idet = -1;

    if (idet >= 0) {
      for (int i = 0; i < T->CDET_Scint.nhits; i++) {
        // Evaluation of number of photoelectrons and time from energy deposit documented at:
        //
        chan = T->CDET_Scint.cell->at(i);
        Npe  = R->Poisson(T->CDET_Scint.sumedep->at(i) * 5.634e3);
        t    = tzero + T->CDET_Scint.tavg->at(i) + 5.75 + T->CDET_Scint.xhit->at(i) / 0.16 - pmtdets[idet]->fTrigOffset;

        pmtdets[idet]->PMTmap[chan].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, signal);
      }
      has_data = true;
    }

    // ** How to add a new subsystem **
    // Unfold here the data for the new detector
    // genrp detectors beam side
    idet = 0;
    while (idet < (int)detmap.size()) {
      // if(idet<0)idet++;
      if (detmap[idet] != PRPOLBS_SCINT_UNIQUE_DETID) {
        idet++;
      } else {
        break;
      }
    }
    if (idet >= detmap.size())
      idet = -1;
    // Process hodoscope data
    if (idet >= 0) { // && T->Harm_PRPolScintBeamSide.nhits){
      for (int i = 0; i < T->Harm_PRPolScintBeamSide.nhits; i++) {
        for (int j = 0; j < 2; j++) { // j = 0: close beam PMT, j = 1: far beam PMT
          Npe = R->Poisson(1.0e7 * T->Harm_PRPolScintBeamSide.sumedep->at(i) * 0.113187 *
                           exp(-(0.25 + pow(-1, j) * T->Harm_PRPolScintBeamSide.xhit->at(i)) / 1.03533) * 0.24);
          t   = tzero + T->Harm_PRPolScintBeamSide.tavg->at(i) +
              (0.5 + pow(-1, j) * T->Harm_PRPolScintBeamSide.xhit->at(i)) / 0.15 - pmtdets[idet]->fTrigOffset;
          chan = T->Harm_PRPolScintBeamSide.cell->at(i) * 2 + j;
          pmtdets[idet]->PMTmap[chan].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, signal);
        }
      }
      has_data = true;
    }
    // genrp detectors far side scint
    idet = 0;
    while (idet < (int)detmap.size()) {
      // if(idet<0)idet++;
      // cout << " " << idet << ":" << detmap[idet];
      if (detmap[idet] != PRPOLFS_SCINT_UNIQUE_DETID) {
        idet++;
      } else {
        // cout << "detmap [idet] = " << detmap[idet] << " =? " << PRPOLFS_SCINT_UNIQUE_DETID << " => PR pol scint
        // unique ID FOUND! " << endl;
        break;
      }
    }
    if (idet >= detmap.size())
      idet = -1;
    if (idet >= 0) { // && T->Harm_PRPolScintFarSide.nhits){
      for (int i = 0; i < T->Harm_PRPolScintFarSide.nhits; i++) {
        for (int j = 0; j < 2; j++) { // j = 0: close beam PMT, j = 1: far beam PMT
          // EPAF: I don't know what scintillator the side hodoscope uses. so I will use the same photon yield and light
          // speed by default. section is ~4x larger, so attenuation should be better, but without better information, I
          // will keep the same attenuation as the hodoscope. From the 3D drawings shown in the presentations, I will
          // assume that each light guide is about 1/2 of the scintillator length (which is 50 cm).

          Npe = R->Poisson(1.0e7 * T->Harm_PRPolScintFarSide.sumedep->at(i) * 0.113187 *
                           exp(-(0.25 + pow(-1, j) * T->Harm_PRPolScintFarSide.xhit->at(i)) / 1.03533) * 0.24);
          t   = tzero + T->Harm_PRPolScintFarSide.tavg->at(i) +
              (0.50 + pow(-1, j) * T->Harm_PRPolScintFarSide.xhit->at(i)) / 0.15 - pmtdets[idet]->fTrigOffset;
          chan = T->Harm_PRPolScintFarSide.cell->at(i) * 2 + j;
          pmtdets[idet]->PMTmap[chan].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, signal);
        }
      }
      has_data = true;
    }

    // genrp detectors Activeana
    idet = 0;
    while (idet < (int)detmap.size()) {
      // if(idet<0)idet++;
      // cout << " " << idet << ":" << detmap[idet];
      if (detmap[idet] != ACTIVEANA_UNIQUE_DETID) {
        idet++;
      } else {
        // cout << "detmap [idet] = " << detmap[idet] << " =? " << ACTIVEANA_UNIQUE_DETID << " => Active analyzer FOUND!
        // " << endl;
        break;
      }
    }
    if (idet >= detmap.size())
      idet = -1;
    if (idet >= 0) { // && T->Harm_ActAnScint.nhits){
      // cout << T->Harm_ActAnScint.nhits << endl;
      for (int i = 0; i < T->Harm_ActAnScint.nhits; i++) {
        // Npe = R->Poisson(Npe_edep_unit*T->Harm_ActAnScint.sumedep->at(i)); //Find Npe_edep_unit: Ave. amount of light
        // produced per energy deposit
        //  The number of photoelectrons for each PMT is the product of the raw number of photoelectrons produced
        //  which depends on the energy deposit sumedep)times the light attenuation
        //  which depends on the distance between the light production and the PMT)
        //  => Npe = (Npe_edep_unit*sumedep)*exp(-(|x_hit-x_PMT|)/Lambda)
        //  EPAF: active analyzer uses same scintillator as hodoscope, hence we can use same photon yield and light
        //  speed. section is 2.56x larger, so attenuation should be better, but without better information, I will keep
        //  the same attenuation as the hodoscope. Also, I assume there is no light guide (or that its length is
        //  negligible...

        Npe = R->Poisson(1.0e7 * T->Harm_ActAnScint.sumedep->at(i) * 0.113187 *
                         exp((T->Harm_ActAnScint.yhit->at(i) + 0.125) / 1.03533) * 0.24);
        t   = tzero + T->Harm_ActAnScint.tavg->at(i) + (0.125 + T->Harm_ActAnScint.yhit->at(i)) / 0.15 -
            pmtdets[idet]->fTrigOffset;
        chan = T->Harm_ActAnScint.cell->at(i);
        // pmtdets[idet]->PMTmap[chan].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, signal);
        pmtdets[idet]->PMTmap[chan].Fill_FADCmode1(Npe, pmtdets[idet]->fThreshold, t, pmtdets[idet]->fSigmaPulse,
                                                   signal);
        // cout << "ActiveAna: chan "<< chan << " Npe " << Npe << " t " << t << ", t_zero = " << tzero << ", t_avg = "
        // << T->Harm_ActAnScint.tavg->at(i) << ", -t_offset = " << -pmtdets[idet]->fTrigOffset << endl;
      }
      has_data = true;
    }
  } // end if(!detmap.empty())

  // GEMs
  if (!gemmap.empty()) {
    idet = 0;
    while (idet < (int)gemmap.size()) {
      if (gemmap[idet] != BBGEM_UNIQUE_DETID) {
        idet++;
      } else {
        break;
      }
    }
    // while(gemmap[idet]!=BBGEM_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
    if (idet >= gemmap.size())
      idet = -1;
    // Now process the GEM data
    if (idet >= 0) { // && T->Earm_BBGEM.nhits){
      for (int k = 0; k < T->Earm_BBGEM.nhits; k++) {
        if (T->Earm_BBGEM.edep->at(k) > 0) {
          SBSDigGEMDet::gemhit hit;
          hit.source = signal;
          // Here... that's one source of errors when we get out of the 4 INFN GEMs patter
          mod = 0;
          // cout << gemdets[idet]->fNPlanes/2 << endl;
          while (mod < gemdets[idet]->fNPlanes / 2) {
            // cout << mod << " " << T->Earm_BBGEM.plane->at(k) << " == ? " << gemdets[idet]->GEMPlanes[mod*2].Layer()
            // << " ; " << (gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) << " <=
            // ? " << T->Earm_BBGEM.xin->at(k) << " <= ? " <<
            // (gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) << endl;
            if ((gemdets[idet]->GEMPlanes[mod * 2].Xoffset() - gemdets[idet]->GEMPlanes[mod * 2].dX() * 0.5) <=
                    T->Earm_BBGEM.xin->at(k) &&
                T->Earm_BBGEM.xin->at(k) <=
                    (gemdets[idet]->GEMPlanes[mod * 2].Xoffset() + gemdets[idet]->GEMPlanes[mod * 2].dX() * 0.5) &&
                T->Earm_BBGEM.plane->at(k) == gemdets[idet]->GEMPlanes[mod * 2].Layer())
              break;
            mod++;
          } // that does the job, but maybe can be optimized???
          if (mod == gemdets[idet]->fNPlanes / 2)
            continue;
          /*
            if(T->Earm_BBGEM.plane->at(k)==5){
              if(fabs(T->Earm_BBGEM.xin->at(k))>=1.024)continue;
              mod = 12 + floor((T->Earm_BBGEM.xin->at(k)+1.024)/0.512);
            }else{
              if(fabs(T->Earm_BBGEM.xin->at(k))>=0.768)continue;
              mod = (T->Earm_BBGEM.plane->at(k)-1)*3 + floor((T->Earm_BBGEM.xin->at(k)+0.768)/0.512);
            }
          */
          // if(mod<2)cout << mod << " " << T->Earm_BBGEM.plane->at(k) << " " << T->Earm_BBGEM.xin->at(k) << endl;
          hit.module = mod;
          hit.edep   = T->Earm_BBGEM.edep->at(k) * 1.0e9; // eV! not MeV!!!!
          // hit.tmin = T->Earm_BBGEM_hit_tmin->at(k);
          // hit.tmax = T->Earm_BBGEM_hit_tmax->at(k);
          hit.t = tzero + T->Earm_BBGEM.t->at(k);
          // cout << mod*2 << " " << gemdets[idet]->GEMPlanes[mod*2].Xoffset() << endl;
          hit.xin  = T->Earm_BBGEM.xin->at(k) - gemdets[idet]->GEMPlanes[mod * 2].Xoffset();
          hit.yin  = T->Earm_BBGEM.yin->at(k) - gemdets[idet]->GEMPlanes[mod * 2 + 1].Xoffset();
          hit.zin  = T->Earm_BBGEM.zin->at(k) - gemdets[idet]->fZLayer[T->Earm_BBGEM.plane->at(k) - 1]; //+0.8031825;
          hit.xout = T->Earm_BBGEM.xout->at(k) - gemdets[idet]->GEMPlanes[mod * 2].Xoffset();
          hit.yout = T->Earm_BBGEM.yout->at(k) - gemdets[idet]->GEMPlanes[mod * 2 + 1].Xoffset();
          hit.zout = T->Earm_BBGEM.zout->at(k) - gemdets[idet]->fZLayer[T->Earm_BBGEM.plane->at(k) - 1]; //+0.8031825;
          // cout << mod << " " << hit.zin << " " << hit.zout << endl;
          gemdets[idet]->fGEMhits.push_back(hit);
        } // end if(sumedep>0)
      }
      has_data = true;
    }

    // SBS GEM detectors
    idet = 0;
    while (idet < (int)gemmap.size()) {
      if (gemmap[idet] != SBSGEM_UNIQUE_DETID) {
        idet++;
      } else {
        break;
      }
    }
    // while(gemmap[idet]!=SBSGEM_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
    if (idet >= gemmap.size())
      idet = -1;
    // Now process the GEM data
    if (idet >= 0) { // && T->Harm_SBSGEM.nhits){
      for (int k = 0; k < T->Harm_SBSGEM.nhits; k++) {
        if (T->Harm_SBSGEM.edep->at(k) > 0) {
          SBSDigGEMDet::gemhit hit;
          hit.source = signal;
          // Here... that's one source of errors when we get out of the 4 INFN GEMs patter
          mod = 0;
          // cout << gemdets[idet]->fNPlanes/2 << endl;
          while (mod < gemdets[idet]->fNPlanes / 2) {
            // cout << mod << " " << T->Harm_SBSGEM.plane->at(k) << " == ? " << gemdets[idet]->GEMPlanes[mod*2].Layer()
            // << " ; " << (gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) << " <=
            // ? " << T->Harm_SBSGEM.xin->at(k) << " <= ? " <<
            // (gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) << endl;
            if ((gemdets[idet]->GEMPlanes[mod * 2].Xoffset() - gemdets[idet]->GEMPlanes[mod * 2].dX() * 0.5) <=
                    T->Harm_SBSGEM.xin->at(k) &&
                T->Harm_SBSGEM.xin->at(k) <=
                    (gemdets[idet]->GEMPlanes[mod * 2].Xoffset() + gemdets[idet]->GEMPlanes[mod * 2].dX() * 0.5) &&
                T->Harm_SBSGEM.plane->at(k) == gemdets[idet]->GEMPlanes[mod * 2].Layer())
              break;
            mod++;
          } // that does the job, but maybe can be optimized???
          if (mod == gemdets[idet]->fNPlanes / 2)
            continue;
          /*
            if(T->Harm_SBSGEM.plane->at(k)==5){
              if(fabs(T->Harm_SBSGEM.xin->at(k))>=1.024)continue;
              mod = 12 + floor((T->Harm_SBSGEM.xin->at(k)+1.024)/0.512);
            }else{
              if(fabs(T->Harm_SBSGEM.xin->at(k))>=0.768)continue;
              mod = (T->Harm_SBSGEM.plane->at(k)-1)*3 + floor((T->Harm_SBSGEM.xin->at(k)+0.768)/0.512);
            }
          */
          // if(mod<2)cout << mod << " " << T->Harm_SBSGEM.plane->at(k) << " " << T->Harm_SBSGEM.xin->at(k) << endl;
          hit.module = mod;
          hit.edep   = T->Harm_SBSGEM.edep->at(k) * 1.0e9; // eV! not MeV!!!!
          // hit.tmin = T->Harm_SBSGEM_hit_tmin->at(k);
          // hit.tmax = T->Harm_SBSGEM_hit_tmax->at(k);
          hit.t = tzero + T->Harm_SBSGEM.t->at(k);
          // cout << mod*2 << " " << gemdets[idet]->GEMPlanes[mod*2].Xoffset() << endl;
          hit.xin  = T->Harm_SBSGEM.xin->at(k) - gemdets[idet]->GEMPlanes[mod * 2].Xoffset();
          hit.yin  = T->Harm_SBSGEM.yin->at(k) - gemdets[idet]->GEMPlanes[mod * 2 + 1].Xoffset();
          hit.zin  = T->Harm_SBSGEM.zin->at(k) - gemdets[idet]->fZLayer[T->Harm_SBSGEM.plane->at(k) - 1]; //+0.8031825;
          hit.xout = T->Harm_SBSGEM.xout->at(k) - gemdets[idet]->GEMPlanes[mod * 2].Xoffset();
          hit.yout = T->Harm_SBSGEM.yout->at(k) - gemdets[idet]->GEMPlanes[mod * 2 + 1].Xoffset();
          hit.zout = T->Harm_SBSGEM.zout->at(k) - gemdets[idet]->fZLayer[T->Harm_SBSGEM.plane->at(k) - 1]; //+0.8031825;
          // cout << mod << " " << hit.zin << " " << hit.zout << endl;
          gemdets[idet]->fGEMhits.push_back(hit);
        } // end if(sumedep>0)
      }
      has_data = true;
    }

    // GEp GEM detectors
    idet = 0;
    while (idet < (int)gemmap.size()) {
      // if(idet<0)idet++;
      if (gemmap[idet] != FT_UNIQUE_DETID) {
        idet++;
      } else {
        break;
      }
    }
    if (idet >= gemmap.size())
      idet = -1;
    if (idet >= 0) {
      for (int k = 0; k < T->Harm_FT.nhits; k++) {
        if (T->Harm_FT.edep->at(k) > 0) {
          SBSDigGEMDet::gemhit hit;
          hit.source = signal;
          // Here... that's one source of errors when we get out of the 4 INFN GEMs patter
          mod = 0;
          // cout << gemdets[idet]->fNPlanes/2 << endl;
          while (mod < gemdets[idet]->fNPlanes / 2) {
            // cout << mod << " " << T->Harm_FT.plane->at(k) << " == ? " << gemdets[idet]->GEMPlanes[mod*2].Layer() << "
            // ; " << (gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) << " <= ? "
            // << T->Harm_FT.xin->at(k) << " <= ? " <<
            // (gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) << endl;
            if ((gemdets[idet]->GEMPlanes[mod * 2].Xoffset() - gemdets[idet]->GEMPlanes[mod * 2].dX() * 0.5) <=
                    T->Harm_FT.xin->at(k) &&
                T->Harm_FT.xin->at(k) <=
                    (gemdets[idet]->GEMPlanes[mod * 2].Xoffset() + gemdets[idet]->GEMPlanes[mod * 2].dX() * 0.5) &&
                T->Harm_FT.plane->at(k) == gemdets[idet]->GEMPlanes[mod * 2].Layer())
              break;
            mod++;
          } // that does the job, but maybe can be optimized???
          if (mod == gemdets[idet]->fNPlanes / 2)
            continue; // cout << mod << endl;

          // if(mod<2)cout << mod << " " << T->Harm_FT.plane->at(k) << " " << T->Harm_FT.xin->at(k) << endl;
          hit.module = mod;
          hit.edep   = T->Harm_FT.edep->at(k) * 1.0e9; // eV! not MeV!!!!
          // hit.tmin = T->Harm_FT_hit_tmin->at(k);
          // hit.tmax = T->Harm_FT_hit_tmax->at(k);
          hit.t = tzero + T->Harm_FT.t->at(k);
          // cout << mod*2 << " " << gemdets[idet]->GEMPlanes[mod*2].Xoffset() << endl;
          hit.xin  = T->Harm_FT.xin->at(k) - gemdets[idet]->GEMPlanes[mod * 2].Xoffset();
          hit.yin  = T->Harm_FT.yin->at(k) - gemdets[idet]->GEMPlanes[mod * 2 + 1].Xoffset();
          hit.zin  = T->Harm_FT.zin->at(k) - gemdets[idet]->fZLayer[T->Harm_FT.plane->at(k) - 1]; //+1.7886925;
          hit.xout = T->Harm_FT.xout->at(k) - gemdets[idet]->GEMPlanes[mod * 2].Xoffset();
          hit.yout = T->Harm_FT.yout->at(k) - gemdets[idet]->GEMPlanes[mod * 2 + 1].Xoffset();
          hit.zout = T->Harm_FT.zout->at(k) - gemdets[idet]->fZLayer[T->Harm_FT.plane->at(k) - 1]; //+1.7886925;
          // cout << T->Harm_FT.plane->at(k) << " " << T->Harm_FT.zin->at(k) << " " <<
          // gemdets[idet]->fZLayer[T->Harm_FT.plane->at(k)-1] << " " << mod << " " << hit.zin << " " << hit.zout <<
          // endl;
          gemdets[idet]->fGEMhits.push_back(hit);
        } // end if(sumedep>0)
      }
      has_data = true;
    }

    idet = 0;
    while (idet < (int)gemmap.size()) {
      // if(idet<0)idet++;
      if (gemmap[idet] != FPP1_UNIQUE_DETID) {
        idet++;
      } else {
        break;
      }
    }
    if (idet >= gemmap.size())
      idet = -1;
    if (idet >= 0) {
      for (int k = 0; k < T->Harm_FPP1.nhits; k++) {
        if (T->Harm_FPP1.edep->at(k) > 0) {
          SBSDigGEMDet::gemhit hit;
          hit.source = signal;
          // Here... that's one source of errors when we get out of the 4 INFN GEMs patter
          mod = 0;
          // cout << gemdets[idet]->fNPlanes/2 << endl;
          while (mod < gemdets[idet]->fNPlanes / 2) {
            // cout << mod << " " << T->Harm_FPP1.plane->at(k) << " == ? " << gemdets[idet]->GEMPlanes[mod*2].Layer() <<
            // " ; " << (gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) << " <= ? "
            // << T->Harm_FPP1.xin->at(k) << " <= ? " <<
            // (gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) << endl;
            if ((gemdets[idet]->GEMPlanes[mod * 2].Xoffset() - gemdets[idet]->GEMPlanes[mod * 2].dX() * 0.5) <=
                    T->Harm_FPP1.xin->at(k) &&
                T->Harm_FPP1.xin->at(k) <=
                    (gemdets[idet]->GEMPlanes[mod * 2].Xoffset() + gemdets[idet]->GEMPlanes[mod * 2].dX() * 0.5) &&
                T->Harm_FPP1.plane->at(k) == gemdets[idet]->GEMPlanes[mod * 2].Layer())
              break;
            mod++;
          } // that does the job, but maybe can be optimized???
          if (mod == gemdets[idet]->fNPlanes / 2)
            continue; // cout << mod << endl;

          // if(mod<2)cout << mod << " " << T->Harm_FPP1.plane->at(k) << " " << T->Harm_FPP1.xin->at(k) << endl;
          hit.module = mod;
          hit.edep   = T->Harm_FPP1.edep->at(k) * 1.0e9; // eV! not MeV!!!!
          // hit.tmin = T->Harm_FPP1_hit_tmin->at(k);
          // hit.tmax = T->Harm_FPP1_hit_tmax->at(k);
          hit.t = tzero + T->Harm_FPP1.t->at(k);
          // cout << mod*2 << " " << gemdets[idet]->GEMPlanes[mod*2].Xoffset() << endl;
          hit.xin  = T->Harm_FPP1.xin->at(k) - gemdets[idet]->GEMPlanes[mod * 2].Xoffset();
          hit.yin  = T->Harm_FPP1.yin->at(k) - gemdets[idet]->GEMPlanes[mod * 2 + 1].Xoffset();
          hit.zin  = T->Harm_FPP1.zin->at(k) - gemdets[idet]->fZLayer[T->Harm_FPP1.plane->at(k) - 1]; //+1.7886925;
          hit.xout = T->Harm_FPP1.xout->at(k) - gemdets[idet]->GEMPlanes[mod * 2].Xoffset();
          hit.yout = T->Harm_FPP1.yout->at(k) - gemdets[idet]->GEMPlanes[mod * 2 + 1].Xoffset();
          hit.zout = T->Harm_FPP1.zout->at(k) - gemdets[idet]->fZLayer[T->Harm_FPP1.plane->at(k) - 1]; //+1.7886925;
          // cout << mod << " " << hit.zin << " " << hit.zout << endl;
          gemdets[idet]->fGEMhits.push_back(hit);
        } // end if(sumedep>0)
      }
      has_data = true;
    }

    idet = 0;
    while (idet < (int)gemmap.size()) {
      // if(idet<0)idet++;
      if (gemmap[idet] != FPP2_UNIQUE_DETID) {
        idet++;
      } else {
        break;
      }
    }
    if (idet >= gemmap.size())
      idet = -1;
    if (idet >= 0) {
      for (int k = 0; k < T->Harm_FPP2.nhits; k++) {
        if (T->Harm_FPP2.edep->at(k) > 0) {
          SBSDigGEMDet::gemhit hit;
          hit.source = signal;
          // Here... that's one source of errors when we get out of the 4 INFN GEMs patter
          mod = 0;
          // cout << gemdets[idet]->fNPlanes/2 << endl;
          while (mod < gemdets[idet]->fNPlanes / 2) {
            // cout << mod << " " << T->Harm_FPP2.plane->at(k) << " == ? " << gemdets[idet]->GEMPlanes[mod*2].Layer() <<
            // " ; " << (gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) << " <= ? "
            // << T->Harm_FPP2.xin->at(k) << " <= ? " <<
            // (gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5) << endl;
            if ((gemdets[idet]->GEMPlanes[mod * 2].Xoffset() - gemdets[idet]->GEMPlanes[mod * 2].dX() * 0.5) <=
                    T->Harm_FPP2.xin->at(k) &&
                T->Harm_FPP2.xin->at(k) <=
                    (gemdets[idet]->GEMPlanes[mod * 2].Xoffset() + gemdets[idet]->GEMPlanes[mod * 2].dX() * 0.5) &&
                T->Harm_FPP2.plane->at(k) == gemdets[idet]->GEMPlanes[mod * 2].Layer())
              break;
            mod++;
          } // that does the job, but maybe can be optimized???
          // cout << "module " << mod << " " << T->Harm_FPP2.plane->at(k) << " " <<
          // gemdets[idet]->GEMPlanes[mod*2].Layer() << " " <<
          // gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5 << " < " <<
          // T->Harm_FPP2.xin->at(k) << " < " <<
          // gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5 << endl;
          if (mod == gemdets[idet]->fNPlanes / 2)
            continue; // cout << mod << endl;

          // if(mod<2)cout << mod << " " << T->Harm_FPP2.plane->at(k) << " " << T->Harm_FPP2.xin->at(k) << endl;
          hit.module = mod;
          hit.edep   = T->Harm_FPP2.edep->at(k) * 1.0e9; // eV! not MeV!!!!
          // hit.tmin = T->Harm_FPP2_hit_tmin->at(k);
          // hit.tmax = T->Harm_FPP2_hit_tmax->at(k);
          hit.t = tzero + T->Harm_FPP2.t->at(k);
          // cout << mod*2 << " " << gemdets[idet]->GEMPlanes[mod*2].Xoffset() << endl;
          hit.xin  = T->Harm_FPP2.xin->at(k) - gemdets[idet]->GEMPlanes[mod * 2].Xoffset();
          hit.yin  = T->Harm_FPP2.yin->at(k) - gemdets[idet]->GEMPlanes[mod * 2 + 1].Xoffset();
          hit.zin  = T->Harm_FPP2.zin->at(k) - gemdets[idet]->fZLayer[T->Harm_FPP2.plane->at(k) - 1]; //+1.7886925;
          hit.xout = T->Harm_FPP2.xout->at(k) - gemdets[idet]->GEMPlanes[mod * 2].Xoffset();
          hit.yout = T->Harm_FPP2.yout->at(k) - gemdets[idet]->GEMPlanes[mod * 2 + 1].Xoffset();
          hit.zout = T->Harm_FPP2.zout->at(k) - gemdets[idet]->fZLayer[T->Harm_FPP2.plane->at(k) - 1]; //+1.7886925;
          // cout << mod << " " << hit.zin << " " << hit.zout << endl;

          gemdets[idet]->fGEMhits.push_back(hit);
        } // end if(sumedep>0)
      }
      has_data = true;
    }

    // GEn-rp GEMs: CEPOL_Front
    idet = 0;
    while (idet < (int)gemmap.size()) {
      if (gemmap[idet] != CEPOL_GEMFRONT_UNIQUE_DETID) {
        idet++;
      } else {
        break;
      }
    }
    if (idet >= gemmap.size())
      idet = -1;
    if (idet >= 0) { // && T->Earm_BBGEM.nhits){
      for (int k = 0; k < T->Harm_CEPolFront.nhits; k++) {
        if (T->Harm_CEPolFront.edep->at(k) > 0) {
          SBSDigGEMDet::gemhit hit;
          hit.source = signal;
          mod        = 0;
          while (mod < gemdets[idet]->fNPlanes / 2) {
            if ((gemdets[idet]->GEMPlanes[mod * 2].Xoffset() - gemdets[idet]->GEMPlanes[mod * 2].dX() * 0.5) <=
                    T->Harm_CEPolFront.xin->at(k) &&
                T->Harm_CEPolFront.xin->at(k) <=
                    (gemdets[idet]->GEMPlanes[mod * 2].Xoffset() + gemdets[idet]->GEMPlanes[mod * 2].dX() * 0.5) &&
                T->Harm_CEPolFront.plane->at(k) == gemdets[idet]->GEMPlanes[mod * 2].Layer())
              break;
            mod++;
          } // that does the job, but maybe can be optimized???
          if (mod == gemdets[idet]->fNPlanes / 2)
            continue;
          hit.module = mod;
          hit.edep   = T->Harm_CEPolFront.edep->at(k) * 1.0e9; // eV! not MeV!!!!
          hit.t      = tzero + T->Harm_CEPolFront.t->at(k);
          hit.xin    = T->Harm_CEPolFront.xin->at(k) - gemdets[idet]->GEMPlanes[mod * 2].Xoffset();
          hit.yin    = T->Harm_CEPolFront.yin->at(k) - gemdets[idet]->GEMPlanes[mod * 2 + 1].Xoffset();
          hit.zin =
              T->Harm_CEPolFront.zin->at(k) - gemdets[idet]->fZLayer[T->Harm_CEPolFront.plane->at(k) - 1]; //+0.8031825;
          hit.xout = T->Harm_CEPolFront.xout->at(k) - gemdets[idet]->GEMPlanes[mod * 2].Xoffset();
          hit.yout = T->Harm_CEPolFront.yout->at(k) - gemdets[idet]->GEMPlanes[mod * 2 + 1].Xoffset();
          hit.zout = T->Harm_CEPolFront.zout->at(k) -
                     gemdets[idet]->fZLayer[T->Harm_CEPolFront.plane->at(k) - 1]; //+0.8031825;
          gemdets[idet]->fGEMhits.push_back(hit);
          //	    cout<<" Harm_CEPolFront  "<<"  zin  "<<hit.zin<<"  zout  "<<hit.zout<<" plane
          //"<<T->Harm_CEPolFront.plane->at(k)<<endl;
        } // end if(sumedep>0)
      }
      has_data = true;
    }

    // GEn-rp GEMs: CEPOL_Rear
    idet = 0;
    while (idet < (int)gemmap.size()) {
      if (gemmap[idet] != CEPOL_GEMREAR_UNIQUE_DETID) {
        idet++;

      } else {
        break;
      }
    }
    if (idet >= gemmap.size())
      idet = -1;
    if (idet >= 0) { // && T->Earm_BBGEM.nhits){

      for (int k = 0; k < T->Harm_CEPolRear.nhits; k++) {
        if (T->Harm_CEPolRear.edep->at(k) > 0) {
          SBSDigGEMDet::gemhit hit;
          hit.source = signal;
          mod        = 0;
          while (mod < gemdets[idet]->fNPlanes / 2) {
            if ((gemdets[idet]->GEMPlanes[mod * 2].Xoffset() - gemdets[idet]->GEMPlanes[mod * 2].dX() * 0.5) <=
                    T->Harm_CEPolRear.xin->at(k) &&
                T->Harm_CEPolRear.xin->at(k) <=
                    (gemdets[idet]->GEMPlanes[mod * 2].Xoffset() + gemdets[idet]->GEMPlanes[mod * 2].dX() * 0.5) &&
                T->Harm_CEPolRear.plane->at(k) == gemdets[idet]->GEMPlanes[mod * 2].Layer())
              break;
            mod++;
          } // that does the job, but maybe can be optimized???
          if (mod == gemdets[idet]->fNPlanes / 2)
            continue;
          hit.module = mod;
          hit.edep   = T->Harm_CEPolRear.edep->at(k) * 1.0e9; // eV! not MeV!!!!
          hit.t      = tzero + T->Harm_CEPolRear.t->at(k);
          hit.xin    = T->Harm_CEPolRear.xin->at(k) - gemdets[idet]->GEMPlanes[mod * 2].Xoffset();
          hit.yin    = T->Harm_CEPolRear.yin->at(k) - gemdets[idet]->GEMPlanes[mod * 2 + 1].Xoffset();
          hit.zin =
              T->Harm_CEPolRear.zin->at(k) - gemdets[idet]->fZLayer[T->Harm_CEPolRear.plane->at(k) - 1]; //+0.8031825;
          hit.xout = T->Harm_CEPolRear.xout->at(k) - gemdets[idet]->GEMPlanes[mod * 2].Xoffset();
          hit.yout = T->Harm_CEPolRear.yout->at(k) - gemdets[idet]->GEMPlanes[mod * 2 + 1].Xoffset();
          hit.zout =
              T->Harm_CEPolRear.zout->at(k) - gemdets[idet]->fZLayer[T->Harm_CEPolRear.plane->at(k) - 1]; //+0.8031825;
          gemdets[idet]->fGEMhits.push_back(hit);
          // cout<<" Harm_CEPolRear  "<<"  zin  "<<hit.zin<<"  zout  "<<hit.zout<<endl;

        } // end if(sumedep>0)
      }
      has_data = true;
    }
    // GEn-rp GEMs: prpolbs_gem
    idet = 0;
    while (idet < (int)gemmap.size()) {
      if (gemmap[idet] != PRPOLBS_GEM_UNIQUE_DETID) {
        idet++;
      } else {
        break;
      }
    }
    if (idet >= gemmap.size())
      idet = -1;
    if (idet >= 0) { // && T->Earm_BBGEM.nhits){
      for (int k = 0; k < T->Harm_PrPolGEMBeamSide.nhits; k++) {
        // cout<<" Nhits_prbol print:   "<<T->Harm_PrPolGEMBeamSide.nhits<<endl;
        if (T->Harm_PrPolGEMBeamSide.edep->at(k) > 0) {
          SBSDigGEMDet::gemhit hit;
          hit.source = signal;
          mod        = 0;
          while (mod < gemdets[idet]->fNPlanes / 2) {
            if ((gemdets[idet]->GEMPlanes[mod * 2].Xoffset() - gemdets[idet]->GEMPlanes[mod * 2].dX() * 0.5) <=
                    T->Harm_PrPolGEMBeamSide.xin->at(k) &&
                T->Harm_PrPolGEMBeamSide.xin->at(k) <=
                    (gemdets[idet]->GEMPlanes[mod * 2].Xoffset() + gemdets[idet]->GEMPlanes[mod * 2].dX() * 0.5) &&
                T->Harm_PrPolGEMBeamSide.plane->at(k) == gemdets[idet]->GEMPlanes[mod * 2].Layer())
              break;
            mod++;
          } // that does the job, but maybe can be optimized???
          if (mod == gemdets[idet]->fNPlanes / 2)
            continue;
          hit.module = mod;
          hit.edep   = T->Harm_PrPolGEMBeamSide.edep->at(k) * 1.0e9; // eV! not MeV!!!!
          hit.t      = tzero + T->Harm_PrPolGEMBeamSide.t->at(k);
          hit.xin    = T->Harm_PrPolGEMBeamSide.xin->at(k) - gemdets[idet]->GEMPlanes[mod * 2].Xoffset();
          hit.yin    = T->Harm_PrPolGEMBeamSide.zin->at(k) - gemdets[idet]->GEMPlanes[mod * 2 + 1].Xoffset() - 1.36;
          hit.zin    = T->Harm_PrPolGEMBeamSide.yin->at(k) -
                    gemdets[idet]->fZLayer[T->Harm_PrPolGEMBeamSide.plane->at(k) - 1]; //+0.8031825;
          hit.xout = T->Harm_PrPolGEMBeamSide.xout->at(k) - gemdets[idet]->GEMPlanes[mod * 2].Xoffset();
          hit.yout = T->Harm_PrPolGEMBeamSide.zout->at(k) - gemdets[idet]->GEMPlanes[mod * 2 + 1].Xoffset() - 1.36;
          hit.zout = T->Harm_PrPolGEMBeamSide.yout->at(k) -
                     gemdets[idet]->fZLayer[T->Harm_PrPolGEMBeamSide.plane->at(k) - 1]; //+0.8031825;
          gemdets[idet]->fGEMhits.push_back(hit);
          // cout<<" Harm_PrPolGEMBeamSide  "<<"  zin  "<<hit.zin<<"  zout  "<<hit.zout<<" plane
          // "<<T->Harm_PrPolGEMBeamSide.plane->at(k)<<endl;
        } // end if(sumedep>0)
      }
      has_data = true;
    }

    // GEn-rp GEMs: prpolfs_gem
    idet = 0;
    while (idet < (int)gemmap.size()) {
      if (gemmap[idet] != PRPOLFS_GEM_UNIQUE_DETID) {
        idet++;

      } else {
        break;
      }
    }
    if (idet >= gemmap.size())
      idet = -1;
    if (idet >= 0) { // && T->Earm_BBGEM.nhits){
      for (int k = 0; k < T->Harm_PrPolGEMFarSide.nhits; k++) {
        if (T->Harm_PrPolGEMFarSide.edep->at(k) > 0) {
          SBSDigGEMDet::gemhit hit;
          hit.source = signal;
          mod        = 0;
          while (mod < gemdets[idet]->fNPlanes / 2) {
            if ((gemdets[idet]->GEMPlanes[mod * 2].Xoffset() - gemdets[idet]->GEMPlanes[mod * 2].dX() * 0.5) <=
                    T->Harm_PrPolGEMFarSide.xin->at(k) &&
                T->Harm_PrPolGEMFarSide.xin->at(k) <=
                    (gemdets[idet]->GEMPlanes[mod * 2].Xoffset() + gemdets[idet]->GEMPlanes[mod * 2].dX() * 0.5) &&
                T->Harm_PrPolGEMFarSide.plane->at(k) == gemdets[idet]->GEMPlanes[mod * 2].Layer())
              break;
            // cout << gemdets[idet]->GEMPlanes[mod*2].Xoffset()-gemdets[idet]->GEMPlanes[mod*2].dX()*0.5 << " <? " <<
            // T->Harm_PrPolGEMFarSide.xin->at(k) << " <? " <<
            // gemdets[idet]->GEMPlanes[mod*2].Xoffset()+gemdets[idet]->GEMPlanes[mod*2].dX()*0.5 << endl; cout <<
            // T->Harm_PrPolGEMFarSide.plane->at(k) << " =? " << gemdets[idet]->GEMPlanes[mod*2].Layer() << endl;
            mod++;
          } // that does the job, but maybe can be optimized???
          if (mod == gemdets[idet]->fNPlanes / 2)
            continue;
          hit.module = mod;
          hit.edep   = T->Harm_PrPolGEMFarSide.edep->at(k) * 1.0e9; // eV! not MeV!!!!
          hit.t      = tzero + T->Harm_PrPolGEMFarSide.t->at(k);
          hit.xin    = T->Harm_PrPolGEMFarSide.xin->at(k) - gemdets[idet]->GEMPlanes[mod * 2].Xoffset();
          hit.yin    = T->Harm_PrPolGEMFarSide.zin->at(k) - gemdets[idet]->GEMPlanes[mod * 2 + 1].Xoffset();
          hit.zin =
              -(T->Harm_PrPolGEMFarSide.yin->at(k) - gemdets[idet]->fZLayer[T->Harm_PrPolGEMFarSide.plane->at(k) - 1]);
          hit.xout = T->Harm_PrPolGEMFarSide.xout->at(k) - gemdets[idet]->GEMPlanes[mod * 2].Xoffset();
          hit.yout = T->Harm_PrPolGEMFarSide.zout->at(k) - gemdets[idet]->GEMPlanes[mod * 2 + 1].Xoffset();
          hit.zout =
              -(T->Harm_PrPolGEMFarSide.yout->at(k) - gemdets[idet]->fZLayer[T->Harm_PrPolGEMFarSide.plane->at(k) - 1]);
          gemdets[idet]->fGEMhits.push_back(hit);
          // cout<<" Harm_PrPolGEMFarSide  "<<"  zin  "<< hit.zin <<"  zout   "<< hit.zout <<endl;
        } // end if(sumedep>0)
      }
      has_data = true;
    }
  }
  return has_data;
}
