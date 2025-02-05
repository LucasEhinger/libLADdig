// includes: standard
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

// includes: root
#include "TChain.h"
#include "TChainElement.h"
#include "TCut.h"
#include "TEventList.h"
#include "TMath.h"
#include "TObjString.h"
#include "TRandom3.h"
#include "TString.h"
#include <TROOT.h>

// includes: specific
#include "SBSDigAuxi.h"
#include "SBSDigBkgdGen.h"
#include "SBSDigGEMDet.h"
#include "SBSDigGEMSimDig.h"
#include "SBSDigPMTDet.h"
#include "g4sbs_tree.h"
#include "g4sbs_types.h"

#ifdef __APPLE__
#include "unistd.h"
#endif

#include <chrono>

using namespace std;
//____________________________________________________
int main(int argc, char **argv) {
  // Step 0: read out arguments
  string db_file, inputsigfile, inputbkgdfile = ""; // sources of files
  ULong64_t Nentries = -1;                          // number of events to process
  // UShort_t Nbkgd = 0;//number of background files to add to each event
  double BkgdTimeWindow = 0, LumiFrac = 0;
  bool pmtbkgddig = false;

  if (argc < 3 || argc > 4) {
    cout << "*** Inadequate number of arguments! ***" << endl
         << " Arguments: database (mandatory); " << endl
         << "           list_of_sig_input_files (str, mandatory); " << endl
         << "          nb_of_sig_evts_to_process (int, def=-1); " << endl;
    //<< "         bkgd_histo_input_file (str, def=''); " << endl
    // << "        bkgd_lumi_frac (double, def=0); " << endl;
    return (-1);
  }

  db_file = argv[1];
  cout << " database file " << db_file << endl;
  inputsigfile = argv[2];
  cout << " Signal input files from: " << inputsigfile << endl;
  if (argc > 3)
    Nentries = atoi(argv[3]);
  cout << " Number of (signal) events to process = " << Nentries << endl;
  /*
  if(argc>5){
    inputbkgdfile = argv[4];
    cout << " Background histgrams from: " << inputbkgdfile << endl;
    LumiFrac = max(0., atof(argv[5]));
    cout << " Fraction of background to superimpose to signal = " << LumiFrac << endl;
  }
  */

  // ------------------- // dev notes // ------------------- //
  // First, we want to extend the input tree (for signal only!!!)
  // I guess in order to avoid adding extra layers of code,
  // the tree extension might have to be coded in the custom tree class

  std::vector<SBSDigPMTDet *> PMTdetectors;
  std::vector<int> detmap;
  std::vector<SBSDigGEMDet *> GEMdetectors;
  std::vector<SBSDigGEMSimDig *> GEMsimDig;
  std::vector<int> gemdetmap;

  // Variable parameters.
  // Can be configured with the database, but are provided with defaults.
  Int_t Rseed            = 0;
  Double_t TriggerJitter = 3.0;

  std::vector<TString> detectors_list;

  // const int nparam_pmtdet_adc  = 12;
  const int nparam_pmtdet_fadc = 13;
  const int nparam_gemdet      = 12;

  int nparam_LADhodo_read    = 0;
  Int_t NChan_LADhodo        = 180;
  Double_t gatewidth_LADhodo = 100.;
  std::vector<Double_t> gain_LADhodo;
  Double_t ped_LADhodo        = 0.;
  Double_t pedsigma_LADhodo   = 0.;
  Double_t trigoffset_LADhodo = 18.6;
  Double_t threshold_LADhodo  = 3.e3;
  Double_t ADCconv_LADhodo    = 100.;
  Int_t ADCbits_LADhodo       = 12;
  Double_t TDCconv_LADhodo    = 0.1;
  Int_t TDCbits_LADhodo       = 19;
  Double_t sigmapulse_LADhodo = 1.6;
  std::vector<Double_t> timeoffset_LADhodo(110, 0);

  int nparam_LADgem_read    = 0;
  Int_t NPlanes_LADgem      = 32; // number of planes/modules/readout
  Double_t gatewidth_LADgem = 400.;
  Double_t ZsupThr_LADgem   = 240.;
  Int_t Nlayers_LADgem      = 5;
  std::vector<Double_t> LADgem_layer_z;
  Int_t *layer_LADgem;
  Int_t *nstrips_LADgem;
  Double_t *offset_LADgem;
  Double_t *RO_angle_LADgem;
  Double_t *triggeroffset_LADgem;
  Double_t *gain_LADgem; // one gain per module
  Double_t *commonmode_array_LADgem;
  UShort_t nAPV_LADgem = 0;

  //Unused for LAD
  int nparam_sbsgem_read    = 0;
  Int_t NPlanes_sbsgem      = 32; // number of planes/modules/readout
  Double_t gatewidth_sbsgem = 400.;
  Double_t ZsupThr_sbsgem   = 240.;
  Int_t Nlayers_sbsgem      = 5;
  std::vector<Double_t> sbsgem_layer_z;
  Int_t *layer_sbsgem;
  Int_t *nstrips_sbsgem;
  Double_t *offset_sbsgem;
  Double_t *RO_angle_sbsgem;
  Double_t *triggeroffset_sbsgem;
  Double_t *gain_sbsgem; // one gain per module
  Double_t *commonmode_array_sbsgem;
  UShort_t nAPV_sbsgem = 0;

  //-----------------------------
  //  Read database
  //-----------------------------
  cout << "read database: " << db_file.c_str() << endl;
  ifstream in_db(db_file.c_str());
  if (!in_db.is_open()) {
    cout << "database " << db_file.c_str() << " does not exist!!!" << endl;
    exit(-1);
  }

  TString currentline;
  while (currentline.ReadLine(in_db) && !currentline.BeginsWith("endconfig")) {
    if (!currentline.BeginsWith("#")) {
      Int_t ntokens = 0;
      std::unique_ptr<TObjArray> tokens(currentline.Tokenize(", \t"));
      if (!tokens->IsEmpty()) {
        ntokens = tokens->GetLast() + 1;
      }
      // TObjArray *tokens = currentline.Tokenize(" ");//vg: def lost => versions prior to 6.06; should be fixed! ???
      // int ntokens = tokens->GetEntries();

      if (ntokens >= 2) {
        TString skey = ((TObjString *)(*tokens)[0])->GetString();

        if (skey == "Rseed") {
          TString stemp = ((TObjString *)(*tokens)[1])->GetString();
          Rseed         = stemp.Atoi();
        }

        if (skey == "TriggerJitter") {
          TString stemp = ((TObjString *)(*tokens)[1])->GetString();
          TriggerJitter = stemp.Atof();
        }

        if (skey == "detectors_list") {
          for (int k = 1; k < ntokens; k++) {
            TString sdet = ((TObjString *)(*tokens)[k])->GetString();
            detectors_list.push_back(sdet);
          }
        }

        // BBHODO
        if (skey == "NChan_LADhodo") {
          TString stemp = ((TObjString *)(*tokens)[1])->GetString();
          NChan_LADhodo  = stemp.Atoi();
          nparam_LADhodo_read++;
        }

        if (skey == "gatewidth_LADhodo") {
          TString stemp    = ((TObjString *)(*tokens)[1])->GetString();
          gatewidth_LADhodo = stemp.Atof();
          nparam_LADhodo_read++;
        }

        if (skey == "gain_LADhodo") {
          gain_LADhodo.resize(NChan_LADhodo);
          if (ntokens == NChan_LADhodo + 1) {
            for (int k = 0; k < NChan_LADhodo; k++) {
              TString stemp  = ((TObjString *)(*tokens)[k + 1])->GetString();
              gain_LADhodo[k] = stemp.Atof() * qe;
            }
          } else {
            cout << ntokens - 1 << " entries for " << skey << " dont match number of channels = " << NChan_LADhodo
                 << endl
                 << " applying first value on all planes " << endl;
            TString stemp = ((TObjString *)(*tokens)[1])->GetString();
            for (int k = 0; k < NChan_LADhodo; k++) {
              gain_LADhodo[k] = stemp.Atof() * qe;
            }
          }
          nparam_LADhodo_read++;
        }

        if (skey == "ped_LADhodo") {
          TString stemp = ((TObjString *)(*tokens)[1])->GetString();
          ped_LADhodo    = stemp.Atof();
          nparam_LADhodo_read++;
        }

        if (skey == "pedsigma_LADhodo") {
          TString stemp   = ((TObjString *)(*tokens)[1])->GetString();
          pedsigma_LADhodo = stemp.Atof();
          nparam_LADhodo_read++;
        }

        if (skey == "trigoffset_LADhodo") {
          TString stemp     = ((TObjString *)(*tokens)[1])->GetString();
          trigoffset_LADhodo = stemp.Atof();
          nparam_LADhodo_read++;
        }

        if (skey == "threshold_LADhodo") {
          TString stemp    = ((TObjString *)(*tokens)[1])->GetString();
          threshold_LADhodo = stemp.Atof();
          nparam_LADhodo_read++;
        }

        if (skey == "ADCconv_LADhodo") {
          TString stemp  = ((TObjString *)(*tokens)[1])->GetString();
          ADCconv_LADhodo = stemp.Atof();
          nparam_LADhodo_read++;
        }

        if (skey == "ADCbits_LADhodo") {
          TString stemp  = ((TObjString *)(*tokens)[1])->GetString();
          ADCbits_LADhodo = stemp.Atoi();
          nparam_LADhodo_read++;
        }

        if (skey == "TDCconv_LADhodo") {
          TString stemp  = ((TObjString *)(*tokens)[1])->GetString();
          TDCconv_LADhodo = stemp.Atof();
          nparam_LADhodo_read++;
        }

        if (skey == "TDCbits_LADhodo") {
          TString stemp  = ((TObjString *)(*tokens)[1])->GetString();
          TDCbits_LADhodo = stemp.Atoi();
          nparam_LADhodo_read++;
        }

        if (skey == "sigmapulse_LADhodo") {
          TString stemp     = ((TObjString *)(*tokens)[1])->GetString();
          sigmapulse_LADhodo = stemp.Atof();
          nparam_LADhodo_read++;
        }

        if (skey == "timeoffset_LADhodo") {
          timeoffset_LADhodo.resize(NChan_LADhodo);
          if (ntokens == NChan_LADhodo + 1) {
            for (int k = 0; k < NChan_LADhodo; k++) {
              TString stemp        = ((TObjString *)(*tokens)[k + 1])->GetString();
              timeoffset_LADhodo[k] = stemp.Atof();
            }
          } else {
            cout << ntokens - 1 << " entries for " << skey << " dont match number of channels = " << NChan_LADhodo
                 << endl
                 << " applying first value on all planes " << endl;
            TString stemp = ((TObjString *)(*tokens)[1])->GetString();
            for (int k = 0; k < NChan_LADhodo; k++) {
              timeoffset_LADhodo[k] = stemp.Atof();
            }
          }
          nparam_LADhodo_read++;
        }

        // GEMs
        if (skey == "NPlanes_LADgem") {
          cout << "reading " << skey.Data() << endl;
          TString stemp = ((TObjString *)(*tokens)[1])->GetString();
          NPlanes_LADgem = stemp.Atoi();

          layer_LADgem         = new Int_t[NPlanes_LADgem];
          nstrips_LADgem       = new Int_t[NPlanes_LADgem];
          offset_LADgem        = new Double_t[NPlanes_LADgem];
          RO_angle_LADgem      = new Double_t[NPlanes_LADgem];
          triggeroffset_LADgem = new Double_t[NPlanes_LADgem / 2];
          gain_LADgem          = new Double_t[NPlanes_LADgem / 2];
          nparam_LADgem_read++;
        }

        if (skey == "gatewidth_LADgem") {
          cout << "reading " << skey.Data() << endl;
          TString stemp   = ((TObjString *)(*tokens)[1])->GetString();
          gatewidth_LADgem = stemp.Atof();
          nparam_LADgem_read++;
        }

        if (skey == "ZsupThr_LADgem") {
          cout << "reading " << skey.Data() << endl;
          TString stemp = ((TObjString *)(*tokens)[1])->GetString();
          ZsupThr_LADgem = stemp.Atof();
          nparam_LADgem_read++;
        }

        if (skey == "nlayers_LADgem") {
          cout << "reading " << skey.Data() << endl;
          TString stemp = ((TObjString *)(*tokens)[1])->GetString();
          Nlayers_LADgem = stemp.Atof();
          nparam_LADgem_read++;
        }

        if (skey == "LADgem_layer_z") {
          cout << "reading " << skey.Data() << endl;
          TString stemp = ((TObjString *)(*tokens)[1])->GetString();
          if (ntokens == Nlayers_LADgem + 1) {
            for (int k = 1; k < ntokens; k++) {
              TString stemp = ((TObjString *)(*tokens)[k])->GetString();
              LADgem_layer_z.push_back(stemp.Atof());
            }
          } else {
            cout << "number of entries for LADgem_layer_z = " << ntokens - 1
                 << " don't match nlayers = " << Nlayers_LADgem << endl;
            cout << "fix your db " << endl;
          }
          nparam_LADgem_read++;
        }

        if (skey == "layer_LADgem") {
          cout << "reading " << skey.Data() << endl;
          if (ntokens == NPlanes_LADgem + 1) {
            for (int k = 1; k < ntokens; k++) {
              TString stemp      = ((TObjString *)(*tokens)[k])->GetString();
              layer_LADgem[k - 1] = stemp.Atoi();
            }
          } else {
            cout << "number of entries for layer_LADgem = " << ntokens - 1 << " don't match Nplanes = " << NPlanes_LADgem
                 << endl;
            cout << "fix your db " << endl;
            exit(-1);
          }
          nparam_LADgem_read++;
        }

        if (skey == "nstrips_LADgem") {
          cout << "reading " << skey.Data() << endl;
          if (ntokens == NPlanes_LADgem + 1) {
            for (int k = 1; k < ntokens; k++) {
              TString stemp        = ((TObjString *)(*tokens)[k])->GetString();
              nstrips_LADgem[k - 1] = stemp.Atoi();
            }
          } else {
            cout << "number of entries for nstrips_LADgem = " << ntokens - 1
                 << " don't match Nplanes = " << NPlanes_LADgem << endl;
            cout << "fix your db " << endl;
            exit(-1);
          }
          nparam_LADgem_read++;
        }

        if (skey == "offset_LADgem") {
          cout << "reading " << skey.Data() << endl;
          if (ntokens == NPlanes_LADgem + 1) {
            for (int k = 1; k < ntokens; k++) {
              TString stemp       = ((TObjString *)(*tokens)[k])->GetString();
              offset_LADgem[k - 1] = stemp.Atof();
            }
          } else {
            cout << "number of entries for offset_LADgem = " << ntokens - 1 << " don't match Nplanes = " << NPlanes_LADgem
                 << endl;
            cout << "fix your db " << endl;
            exit(-1);
          }
          nparam_LADgem_read++;
        }

        if (skey == "RO_angle_LADgem") {
          cout << "reading " << skey.Data() << endl;
          if (ntokens == NPlanes_LADgem + 1) {
            for (int k = 1; k < ntokens; k++) {
              TString stemp         = ((TObjString *)(*tokens)[k])->GetString();
              RO_angle_LADgem[k - 1] = stemp.Atof() * TMath::DegToRad();
            }
          } else {
            cout << "number of entries for RO_angle_LADgem = " << ntokens - 1
                 << " don't match Nplanes = " << NPlanes_LADgem << endl;
            cout << "fix your db " << endl;
            exit(-1);
          }
          nparam_LADgem_read++;
        }

        if (skey == "triggeroffset_LADgem") {
          cout << "reading " << skey.Data() << endl;
          if (ntokens == NPlanes_LADgem / 2 + 1) {
            for (int k = 1; k < ntokens; k++) {
              TString stemp              = ((TObjString *)(*tokens)[k])->GetString();
              triggeroffset_LADgem[k - 1] = stemp.Atof();
            }
          } else {
            cout << "number of entries for triggeroffset_LADgem = " << ntokens - 1
                 << " don't match Nplanes = " << NPlanes_LADgem << endl;
            cout << "fix your db " << endl;
            exit(-1);
          }
          nparam_LADgem_read++;
        }

        if (skey == "gain_LADgem") {
          cout << "reading " << skey.Data() << endl;
          if (ntokens == NPlanes_LADgem / 2 + 1) {
            for (int k = 1; k < ntokens; k++) {
              TString stemp     = ((TObjString *)(*tokens)[k])->GetString();
              gain_LADgem[k - 1] = stemp.Atof();
            }
          } else {
            cout << "number of entries for gain_LADgem = " << ntokens - 1
                 << " don't match Nplanes/2 = " << NPlanes_LADgem / 2 << endl;
            if (ntokens >= 2) {
              cout << "applying first value on all planes " << endl;
              TString stemp = ((TObjString *)(*tokens)[1])->GetString();
              for (int k = 0; k < NPlanes_LADgem / 2; k++) {
                gain_LADgem[k] = stemp.Atof();
              }
            } else {
              cout << "fix your db " << endl;
              exit(-1);
            }
          }
          nparam_LADgem_read++;
        }

        if (skey == "commonmode_array_LADgem") {
          cout << "reading " << skey.Data() << endl;
          commonmode_array_LADgem = new Double_t[ntokens - 1];
          for (int k = 1; k < ntokens; k++) {
            nAPV_LADgem++;
            TString stemp                 = ((TObjString *)(*tokens)[k])->GetString();
            commonmode_array_LADgem[k - 1] = stemp.Atof();
          }
          nparam_LADgem_read++;
        }

        // SBSGEM
        if (skey == "NPlanes_sbsgem") {
          cout << "reading " << skey.Data() << endl;
          TString stemp  = ((TObjString *)(*tokens)[1])->GetString();
          NPlanes_sbsgem = stemp.Atoi();

          layer_sbsgem         = new Int_t[NPlanes_sbsgem];
          nstrips_sbsgem       = new Int_t[NPlanes_sbsgem];
          offset_sbsgem        = new Double_t[NPlanes_sbsgem];
          RO_angle_sbsgem      = new Double_t[NPlanes_sbsgem];
          triggeroffset_sbsgem = new Double_t[NPlanes_sbsgem / 2];
          gain_sbsgem          = new Double_t[NPlanes_sbsgem / 2];
          nparam_sbsgem_read++;
        }

        if (skey == "gatewidth_sbsgem") {
          cout << "reading " << skey.Data() << endl;
          TString stemp    = ((TObjString *)(*tokens)[1])->GetString();
          gatewidth_sbsgem = stemp.Atof();
          nparam_sbsgem_read++;
        }

        if (skey == "ZsupThr_sbsgem") {
          cout << "reading " << skey.Data() << endl;
          TString stemp  = ((TObjString *)(*tokens)[1])->GetString();
          ZsupThr_sbsgem = stemp.Atof();
          nparam_sbsgem_read++;
        }

        if (skey == "nlayers_sbsgem") {
          cout << "reading " << skey.Data() << endl;
          TString stemp  = ((TObjString *)(*tokens)[1])->GetString();
          Nlayers_sbsgem = stemp.Atof();
          nparam_sbsgem_read++;
        }

        if (skey == "sbsgem_layer_z") {
          cout << "reading " << skey.Data() << endl;
          TString stemp = ((TObjString *)(*tokens)[1])->GetString();
          if (ntokens == Nlayers_sbsgem + 1) {
            for (int k = 1; k < ntokens; k++) {
              TString stemp = ((TObjString *)(*tokens)[k])->GetString();
              sbsgem_layer_z.push_back(stemp.Atof());
            }
          } else {
            cout << "number of entries for sbsgem_layer_z = " << ntokens - 1
                 << " don't match nlayers = " << Nlayers_sbsgem << endl;
            cout << "fix your db " << endl;
          }
          nparam_sbsgem_read++;
        }

        if (skey == "layer_sbsgem") {
          cout << "reading " << skey.Data() << endl;
          if (ntokens == NPlanes_sbsgem + 1) {
            for (int k = 1; k < ntokens; k++) {
              TString stemp       = ((TObjString *)(*tokens)[k])->GetString();
              layer_sbsgem[k - 1] = stemp.Atoi();
            }
          } else {
            cout << "number of entries for layer_sbsgem = " << ntokens - 1
                 << " don't match Nplanes = " << NPlanes_sbsgem << endl;
            cout << "fix your db " << endl;
            exit(-1);
          }
          nparam_sbsgem_read++;
        }

        if (skey == "nstrips_sbsgem") {
          cout << "reading " << skey.Data() << endl;
          if (ntokens == NPlanes_sbsgem + 1) {
            for (int k = 1; k < ntokens; k++) {
              TString stemp         = ((TObjString *)(*tokens)[k])->GetString();
              nstrips_sbsgem[k - 1] = stemp.Atoi();
            }
          } else {
            cout << "number of entries for nstrips_sbsgem = " << ntokens - 1
                 << " don't match Nplanes = " << NPlanes_sbsgem << endl;
            cout << "fix your db " << endl;
            exit(-1);
          }
          nparam_sbsgem_read++;
        }

        if (skey == "offset_sbsgem") {
          cout << "reading " << skey.Data() << endl;
          if (ntokens == NPlanes_sbsgem + 1) {
            for (int k = 1; k < ntokens; k++) {
              TString stemp        = ((TObjString *)(*tokens)[k])->GetString();
              offset_sbsgem[k - 1] = stemp.Atof();
            }
          } else {
            cout << "number of entries for offset_sbsgem = " << ntokens - 1
                 << " don't match Nplanes = " << NPlanes_sbsgem << endl;
            cout << "fix your db " << endl;
            exit(-1);
          }
          nparam_sbsgem_read++;
        }

        if (skey == "RO_angle_sbsgem") {
          cout << "reading " << skey.Data() << endl;
          if (ntokens == NPlanes_sbsgem + 1) {
            for (int k = 1; k < ntokens; k++) {
              TString stemp          = ((TObjString *)(*tokens)[k])->GetString();
              RO_angle_sbsgem[k - 1] = stemp.Atof() * TMath::DegToRad();
            }
          } else {
            cout << "number of entries for RO_angle_sbsgem = " << ntokens - 1
                 << " don't match Nplanes = " << NPlanes_sbsgem << endl;
            cout << "fix your db " << endl;
            exit(-1);
          }
          nparam_sbsgem_read++;
        }

        if (skey == "triggeroffset_sbsgem") {
          cout << "reading " << skey.Data() << endl;
          if (ntokens == NPlanes_sbsgem / 2 + 1) {
            for (int k = 1; k < ntokens; k++) {
              TString stemp               = ((TObjString *)(*tokens)[k])->GetString();
              triggeroffset_sbsgem[k - 1] = stemp.Atof();
            }
          } else {
            cout << "number of entries for triggeroffset_sbsgem = " << ntokens - 1
                 << " don't match Nplanes = " << NPlanes_sbsgem << endl;
            cout << "fix your db " << endl;
            exit(-1);
          }
          nparam_sbsgem_read++;
        }

        if (skey == "gain_sbsgem") {
          cout << "reading " << skey.Data() << endl;
          if (ntokens == NPlanes_sbsgem / 2 + 1) {
            for (int k = 1; k < ntokens; k++) {
              TString stemp      = ((TObjString *)(*tokens)[k])->GetString();
              gain_sbsgem[k - 1] = stemp.Atof();
            }
          } else {
            cout << "number of entries for gain_sbsgem = " << ntokens - 1
                 << " don't match Nplanes/2 = " << NPlanes_sbsgem / 2 << endl;
            if (ntokens >= 2) {
              cout << "applying first value on all planes " << endl;
              TString stemp = ((TObjString *)(*tokens)[1])->GetString();
              for (int k = 0; k < NPlanes_sbsgem / 2; k++) {
                gain_sbsgem[k] = stemp.Atof();
              }
            } else {
              cout << "fix your db " << endl;
              exit(-1);
            }
          }
          nparam_sbsgem_read++;
        }

        if (skey == "commonmode_array_sbsgem") {
          cout << "reading " << skey.Data() << endl;
          commonmode_array_sbsgem = new Double_t[ntokens - 1];
          for (int k = 1; k < ntokens; k++) {
            nAPV_sbsgem++;
            TString stemp                  = ((TObjString *)(*tokens)[k])->GetString();
            commonmode_array_sbsgem[k - 1] = stemp.Atof();
          }
          nparam_sbsgem_read++;
        }

      } // end if( ntokens >= 2 )
      tokens->~TObjArray(); // ineffective... :(
    } // end if( !currentline.BeginsWith("#"))
  } // end while

  //-----------------------------
  //  Declare detectors
  //-----------------------------
  cout << " declaring detectors " << endl;
  for (int k = 0; k < detectors_list.size(); k++) {
    cout << "detector: " << detectors_list[k].Data() << "... " << endl;
    if (detectors_list[k] == "LADgem") {
      if (nparam_LADgem_read != nparam_gemdet) {
        cout << detectors_list[k] << " does not have the right number of parameters!!! " << endl
             << " fix database and retry! " << endl;
        exit(-1);
      }

      SBSDigGEMDet *LADgem     = new SBSDigGEMDet(BBGEM_UNIQUE_DETID, NPlanes_LADgem, layer_LADgem, nstrips_LADgem,
                                                 offset_LADgem, RO_angle_LADgem, 6, ZsupThr_LADgem);
      SBSDigGEMSimDig *gemdig = new SBSDigGEMSimDig(NPlanes_LADgem / 2, triggeroffset_LADgem, ZsupThr_LADgem, nAPV_LADgem,
                                                    commonmode_array_LADgem);
      for (int m = 0; m < Nlayers_LADgem; m++) {
        LADgem->fZLayer.push_back(LADgem_layer_z[m]);
      }
      LADgem->fGateWidth = gatewidth_LADgem;

      GEMdetectors.push_back(LADgem);
      gemdetmap.push_back(BBGEM_UNIQUE_DETID);
      GEMsimDig.push_back(gemdig);
      cout << " set up! " << endl;
    }
    if (detectors_list[k] == "sbsgem") {
      if (nparam_sbsgem_read != nparam_gemdet) {
        cout << detectors_list[k] << " does not have the right number of parameters!!! " << endl
             << " fix database and retry! " << endl;
        exit(-1);
      }

      SBSDigGEMDet *sbsgem    = new SBSDigGEMDet(SBSGEM_UNIQUE_DETID, NPlanes_sbsgem, layer_sbsgem, nstrips_sbsgem,
                                                 offset_sbsgem, RO_angle_sbsgem, 6, ZsupThr_sbsgem);
      SBSDigGEMSimDig *gemdig = new SBSDigGEMSimDig(NPlanes_sbsgem / 2, triggeroffset_sbsgem, ZsupThr_sbsgem,
                                                    nAPV_sbsgem, commonmode_array_sbsgem);
      for (int m = 0; m < Nlayers_sbsgem; m++) {
        sbsgem->fZLayer.push_back(sbsgem_layer_z[m]);
      }
      sbsgem->fGateWidth = gatewidth_sbsgem;

      GEMdetectors.push_back(sbsgem);
      gemdetmap.push_back(SBSGEM_UNIQUE_DETID);
      GEMsimDig.push_back(gemdig);
      cout << " set up! " << endl;
    }
    if (detectors_list[k] == "LADhodo") {
      if (nparam_LADhodo_read != nparam_pmtdet_fadc) {
        cout << detectors_list[k] << " does not have the right number of parameters!!! " << endl
             << " fix database and retry! " << endl;
        exit(-1);
      }
      SBSDigPMTDet *LADhodo =
          new SBSDigPMTDet(HODO_UNIQUE_DETID, NChan_LADhodo, gain_LADhodo, sigmapulse_LADhodo, gatewidth_LADhodo);

      LADhodo->fGain       = gain_LADhodo;
      LADhodo->fPedestal   = ped_LADhodo;
      LADhodo->fPedSigma   = pedsigma_LADhodo;
      LADhodo->fTrigOffset = trigoffset_LADhodo;
      LADhodo->fThreshold  = threshold_LADhodo * spe_unit / ROimpedance;
      LADhodo->fGateWidth  = gatewidth_LADhodo;
      LADhodo->fADCconv    = ADCconv_LADhodo;
      LADhodo->fADCbits    = ADCbits_LADhodo;
      LADhodo->fTDCconv    = TDCconv_LADhodo;
      LADhodo->fTDCbits    = TDCbits_LADhodo;
      LADhodo->fTimeOffset = timeoffset_LADhodo;

      PMTdetectors.push_back(LADhodo);
      detmap.push_back(HODO_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
  }

  TRandom3 *R = new TRandom3(Rseed);

  // Step 1: read input files build the input chains
  // build signal chain
  ifstream sig_inputfile(inputsigfile);
  TChain *C_s = new TChain("T");
  while (currentline.ReadLine(sig_inputfile) && !currentline.BeginsWith("endlist")) {
    if (!currentline.BeginsWith("#")) {
      C_s->Add(currentline.Data());
    }
  }
  TObjArray *fileElements_s = C_s->GetListOfFiles();
  TIter next_s(fileElements_s);
  TChainElement *chEl_s = 0;

  // restart reading to get background info
  while (currentline.ReadLine(sig_inputfile) && !currentline.BeginsWith("end_bkgdinfo")) {
    if (!currentline.BeginsWith("#")) {
      Int_t ntokens = 0;
      std::unique_ptr<TObjArray> tokens(currentline.Tokenize(", \t"));
      if (!tokens->IsEmpty()) {
        ntokens = tokens->GetLast() + 1;
      }
      if (ntokens >= 3) {
        inputbkgdfile = ((TObjString *)(*tokens)[0])->GetString();
        TString stemp = ((TObjString *)(*tokens)[1])->GetString();
        // EPAF: this "bkgd time window" notion is confusing.
        // I will replace it with more straightforward stuff:
        // number of events total and current:
        // then I will convert this to a to make it transparent to the code
        double nevtotal = stemp.Atof();
        stemp           = ((TObjString *)(*tokens)[2])->GetString();
        double Ibeam    = stemp.Atof();                     // in uA
        BkgdTimeWindow  = nevtotal * qe / Ibeam / spe_unit; // in ns!
        LumiFrac        = 1.0;
        if (Ibeam <= 0)
          LumiFrac = 0.0;
        if (ntokens == 4) {
          stemp      = ((TObjString *)(*tokens)[3])->GetString();
          pmtbkgddig = stemp.Atoi();
        }
      }
    }
  }

  TFile *f_bkgd;
  SBSDigBkgdGen *BkgdGenerator;
  if (LumiFrac > 0) {
    f_bkgd = TFile::Open(inputbkgdfile.c_str());
    if (f_bkgd->IsZombie()) {
      LumiFrac = 0;
    } else {
      BkgdGenerator = new SBSDigBkgdGen(f_bkgd, detectors_list, BkgdTimeWindow, pmtbkgddig);
      cout << "Includes background from file: " << inputbkgdfile.c_str() << " (integrated on " << BkgdTimeWindow
           << " ns time window);" << endl
           << " assuming " << LumiFrac * 100 << "% luminosity." << endl;
    }
  }

  double Theta_SBS, D_HCal;

  ULong64_t Nev_fs;
  ULong64_t ev_s;

  ULong64_t NEventsTotal = 0;

  int i_fs = 0;
  bool has_data;

  double timeZero;

  std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

  while ((chEl_s = (TChainElement *)next_s())) {
    if (NEventsTotal >= Nentries) {
      break;
    }
    TFile f_s(chEl_s->GetTitle(), "UPDATE");
    if (f_s.IsZombie())
      cout << "File " << chEl_s->GetTitle() << " cannot be found. Please check the path of your file." << endl;
    // run_data = (G4SBSRunData*)f_s.Get("run_data");
    // TODO: fixme. Uncommented rundata to avoid having to import it.
    // G4SBSRunData* run_data = (G4SBSRunData*)f_s.Get("run_data");
    // if(run_data==0){
    //   cout << "File does not have run data available!!! skip!" << endl;
    //   continue;
    // }
    // Theta_SBS = run_data->fSBStheta;
    // D_HCal = run_data->fHCALdist;
    Theta_SBS       = 0.5;
    D_HCal          = 0.5;
    C_s             = (TChain *)f_s.Get("T");
    g4sbs_tree *T_s = new g4sbs_tree(C_s, detectors_list, bool(LumiFrac));

    Nev_fs = C_s->GetEntries();

    for (ev_s = 0; ev_s < Nev_fs; ev_s++, NEventsTotal++) {
      if (NEventsTotal >= Nentries)
        break;
      if (NEventsTotal % 100 == 0)
        cout << NEventsTotal << "/" << Nentries << endl;

      timeZero = R->Gaus(0.0, TriggerJitter);

      // Clear detectors
      for (int k = 0; k < PMTdetectors.size(); k++) {
        if (detmap[k] == HCAL_UNIQUE_DETID || detmap[k] == ECAL_UNIQUE_DETID || detmap[k] == BBPS_UNIQUE_DETID ||
            detmap[k] == BBSH_UNIQUE_DETID || detmap[k] == ACTIVEANA_UNIQUE_DETID) {
          PMTdetectors[k]->Clear(true);
        } else {
          PMTdetectors[k]->Clear();
        }
      }
      for (int k = 0; k < GEMdetectors.size(); k++) {
        GEMdetectors[k]->Clear();
      }

      has_data = false;

      T_s->ClearDigBranches();
      T_s->GetEntry(ev_s);

      // unfold the thing then... but where???
      has_data = UnfoldData(T_s, Theta_SBS, D_HCal, R, PMTdetectors, detmap, GEMdetectors, gemdetmap, timeZero, 0);
      if (!has_data)
        continue;

      // if we want to add background, add background
      if (LumiFrac > 0) {
        // first digitize signal only...
        //	for(int k = 0; k<GEMdetectors.size(); k++){
        //   GEMsimDig[k]->Digitize(GEMdetectors[k], R);
        //   GEMsimDig[k]->CheckOut(GEMdetectors[k], gemdetmap[k], R, T_s, bool(LumiFrac>0));
        // }
        BkgdGenerator->GenerateBkgd(R, PMTdetectors, detmap, GEMdetectors, gemdetmap, LumiFrac);
      }

      // Digitize: PMT detectors
      for (int k = 0; k < PMTdetectors.size(); k++) {
        PMTdetectors[k]->Digitize(T_s, R);
      }
      // digitize: GEMs
      for (int k = 0; k < GEMdetectors.size(); k++) {
        GEMsimDig[k]->Digitize(GEMdetectors[k], R, bool(LumiFrac > 0));
        GEMsimDig[k]->CheckOut(GEMdetectors[k], gemdetmap[k], R, T_s);
      }
      T_s->FillDigBranches();
    } // end loop on signal events
    // if there are debugging histos to write, write them...
    for (int k = 0; k < GEMdetectors.size(); k++) {
      cout << "GEM det ID: " << GEMdetectors[k]->fUniqueID << endl;
      GEMsimDig[k]->write_histos();
      GEMsimDig[k]->print_time_execution();
    }
    if (LumiFrac > 0)
      BkgdGenerator->WriteXCHistos();
    // write expanded tree
    T_s->fChain->Write("", TObject::kOverwrite);
    f_s.Write();
    f_s.Close();
    i_fs++;
  } // end loop on signal files

  std::chrono::time_point<std::chrono::steady_clock> end = std::chrono::steady_clock::now();

  std::chrono::duration<double> diff = end - start;
  cout << " Total time " << std::setprecision(9) << diff.count() << " s " << endl;

  exit(0);
}
