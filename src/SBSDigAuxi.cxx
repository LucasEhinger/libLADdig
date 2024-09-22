#include "SBSDigAuxi.h"
#include "TMath.h"
#include "g4sbs_types.h"

using namespace std;

bool UnfoldData(g4sbs_tree *T, double theta_sbs, double d_hcal, TRandom3 *R, std::vector<SBSDigPMTDet *> pmtdets,
                std::vector<int> detmap, std::vector<SBSDigGEMDet *> gemdets, std::vector<int> gemmap, double tzero,
                int signal) {
  bool has_data = false;

  int Npe;
  double t;

  // These vars are unused, so commenting them out for now (LHE)
  // double x_ref = -d_hcal * sin(theta_sbs);
  // double z_ref = d_hcal * cos(theta_sbs);

  // double z_hit, Npe_Edep_ratio, sigma_tgen;

  // double z_det, pz, E,
  // double beta, sin2thetaC;

  int chan;
  int mod;

  int idet = 0;
  if (!detmap.empty()) {
    // ordering by increasing unique det ID

    // Hodoscope
    idet = 0;
    while (idet < (int)detmap.size()) {
      if (detmap[idet] != HODO_UNIQUE_DETID) {
        idet++;
      } else {
        break;
      }
    }
    if (idet >= detmap.size())
      idet = -1;

    double detector_z   = 2.0;
    double index_of_ref = 2.0;
    if (idet >= 0) { // && T->Earm_BBHodoScint.nhits){
      for (int i = 0; i < T->Earm_BBHodoScint.nhits; i++) {
        for (int j = 0; j < 2;
             j++) { // j = 0: close beam PMT, j = 1: far beam PMT //LHE: Probably the left & right (or t/b)??
          // Evaluation of number of photoelectrons and time from energy deposit documented at:
          // https://hallaweb.jlab.org/dvcslog/SBS/170711_172759/BB_hodoscope_restudy_update_20170711.pdf
          Npe = R->Poisson(1.0e7 * T->Earm_BBHodoScint.sumedep->at(i) * 0.113187 *
                           exp(-(detector_z / 2 + pow(-1, j) * T->Earm_BBHodoScint.zhit->at(i)) / 1.03533) * 0.24);
          t   = tzero + T->Earm_BBHodoScint.tavg->at(i) +
              (0.4 + detector_z / 2 + pow(-1, j) * T->Earm_BBHodoScint.zhit->at(i)) / (0.3 / index_of_ref) -
              pmtdets[idet]->fTrigOffset; // Assume 0.4 is some sort of distance offset for the cables + pmt's.
          t=0; //LHE: This is a hack to view the time walk for now
          chan = 2 * (T->Earm_BBHodoScint.plane->at(i) * 11 + T->Earm_BBHodoScint.paddle->at(i)) + j;
          t += pmtdets[idet]->fTimeOffset[chan];
          // T->Earm_BBHodoScint_hit_sumedep->at(i);
          // if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
          pmtdets[idet]->PMTmap[chan].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, signal);
          // cout << "BBhodo: chan "<< chan << " Npe " << Npe << " t " << t << ", t_zero = " << tzero << ", t_avg = " <<
          // T->Earm_BBHodoScint.tavg->at(i) << ", -t_offset = " << -pmtdets[idet]->fTrigOffset << endl;
        }
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
  }
  return has_data;
}
