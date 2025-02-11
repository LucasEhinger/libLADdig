#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

double hodo_minEdep = 1. / 1000;
double hodo_maxEdep = 1. / 100;
double gem_minEdep  = 1e-7;
double gem_maxEdep  = 5e-7;

void g4_to_dig(string energy = "400") {
  // Open the ROOT file
  string filepath     = "/volatile/hallc/c-lad/ehingerl/G4_LAD/carlos_proton/";
  string infile_name  = filepath + "raw/" + "ScanLAD_proton_" + energy + "MeV_10k_20240205.root";
  string outfile_name = filepath + "dig/" + "ScanLAD_proton_" + energy + "MeV_10k_20240205_dig.root";

  TFile *file = TFile::Open(infile_name.c_str(), "READ");
  if (!file || file->IsZombie()) {
    cerr << "Error opening file" << endl;
    return;
  }

  // Get the hodoposition tree
  TTree *hodoposition = (TTree *)file->Get("hodoposition");
  // Get the hodoenergy tree
  TTree *hodoenergy = (TTree *)file->Get("hodoenergy");
  // Get the gemana tree
  TTree *gemana = (TTree *)file->Get("gemana");
  if (!hodoposition || !hodoenergy || !gemana) {
    cerr << "Error getting trees" << endl;
    file->Close();
    return;
  }

  // Set branch addresses for hodoposition
  vector<double> *vXbar = nullptr;
  vector<double> *vYbar = nullptr;
  vector<double> *vZbar = nullptr;
  vector<double> *vTbar = nullptr;
  vector<int> *vPaddle  = nullptr;
  vector<int> *vTrackID = nullptr;

  hodoposition->SetBranchAddress("vXbar", &vXbar);
  hodoposition->SetBranchAddress("vYbar", &vYbar);
  hodoposition->SetBranchAddress("vZbar", &vZbar);
  hodoposition->SetBranchAddress("vTbar", &vTbar);
  hodoposition->SetBranchAddress("vPaddle", &vPaddle);
  hodoposition->SetBranchAddress("vTrackID", &vTrackID);

  // Set branch addresses for hodoenergy
  int EventID;
  vector<int> *vPadCopy   = nullptr;
  vector<double> *vEneDep = nullptr;
  vector<int> *vPDG       = nullptr;
  vector<int> *vLevel     = nullptr;

  hodoenergy->SetBranchAddress("EventID", &EventID);
  hodoenergy->SetBranchAddress("vPadCopy", &vPadCopy);
  hodoenergy->SetBranchAddress("vEneDep", &vEneDep);
  hodoenergy->SetBranchAddress("vPDG", &vPDG);
  hodoenergy->SetBranchAddress("vLevel", &vLevel);

  // Set branch addresses for gemana
  vector<double> *vXloc     = nullptr;
  vector<double> *vYloc     = nullptr;
  vector<double> *vZloc     = nullptr;
  vector<double> *vXglo     = nullptr;
  vector<double> *vYglo     = nullptr;
  vector<double> *vZglo     = nullptr;
  vector<double> *vTchamber = nullptr;
  vector<int> *vChamber     = nullptr;
  vector<int> *vgPDG        = nullptr;
  vector<int> *vgLevel      = nullptr;

  gemana->SetBranchAddress("vXloc", &vXloc);
  gemana->SetBranchAddress("vYloc", &vYloc);
  gemana->SetBranchAddress("vZloc", &vZloc);
  gemana->SetBranchAddress("vXglo", &vXglo);
  gemana->SetBranchAddress("vYglo", &vYglo);
  gemana->SetBranchAddress("vZglo", &vZglo);
  gemana->SetBranchAddress("vTchamber", &vTchamber);
  gemana->SetBranchAddress("vChamber", &vChamber);
  gemana->SetBranchAddress("vgPDG", &vgPDG);
  gemana->SetBranchAddress("vgLevel", &vgLevel);

  TFile *outfile = TFile::Open(outfile_name.c_str(), "RECREATE");
  if (!outfile || outfile->IsZombie()) {
    std::cerr << "Error creating output file" << std::endl;
    return;
  }
  // Create new tree with only the necessary branches
  int LAD_Hodo_hit_nhits;
  vector<int> *LAD_Hodo_hit_plane      = nullptr;
  vector<int> *LAD_Hodo_hit_paddle     = nullptr;
  vector<double> *LAD_Hodo_hit_xhit    = nullptr;
  vector<double> *LAD_Hodo_hit_yhit    = nullptr;
  vector<double> *LAD_Hodo_hit_zhit    = nullptr;
  vector<double> *LAD_Hodo_hit_tavg    = nullptr;
  vector<double> *LAD_Hodo_hit_beta    = nullptr;
  vector<double> *LAD_Hodo_hit_sumedep = nullptr;

  int LAD_GEM_hit_nhits;
  vector<int> *LAD_GEM_hit_plane    = nullptr;
  vector<double> *LAD_GEM_hit_xin   = nullptr;
  vector<double> *LAD_GEM_hit_yin   = nullptr;
  vector<double> *LAD_GEM_hit_zin   = nullptr;
  vector<double> *LAD_GEM_hit_t     = nullptr;
  vector<double> *LAD_GEM_hit_xout  = nullptr;
  vector<double> *LAD_GEM_hit_yout  = nullptr;
  vector<double> *LAD_GEM_hit_zout  = nullptr;
  vector<double> *LAD_GEM_hit_t_out = nullptr;
  vector<double> *LAD_GEM_hit_edep  = nullptr;

  TTree *outTree = new TTree("T", "Filtered Data");
  outTree->Branch("LAD.Hodo.hit.nhits", &LAD_Hodo_hit_nhits);
  outTree->Branch("LAD.Hodo.hit.plane", &LAD_Hodo_hit_plane);
  outTree->Branch("LAD.Hodo.hit.paddle", &LAD_Hodo_hit_paddle);
  outTree->Branch("LAD.Hodo.hit.xhit", &LAD_Hodo_hit_xhit);
  outTree->Branch("LAD.Hodo.hit.yhit", &LAD_Hodo_hit_yhit);
  outTree->Branch("LAD.Hodo.hit.zhit", &LAD_Hodo_hit_zhit);
  outTree->Branch("LAD.Hodo.hit.tavg", &LAD_Hodo_hit_tavg);
  outTree->Branch("LAD.Hodo.hit.beta", &LAD_Hodo_hit_beta);
  outTree->Branch("LAD.Hodo.hit.sumedep", &LAD_Hodo_hit_sumedep);

  outTree->Branch("LAD.GEM.hit.nhits", &LAD_GEM_hit_nhits);
  outTree->Branch("LAD.GEM.hit.plane", &LAD_GEM_hit_plane);
  outTree->Branch("LAD.GEM.hit.xin", &LAD_GEM_hit_xin);
  outTree->Branch("LAD.GEM.hit.yin", &LAD_GEM_hit_yin);
  outTree->Branch("LAD.GEM.hit.zin", &LAD_GEM_hit_zin);
  outTree->Branch("LAD.GEM.hit.t", &LAD_GEM_hit_t);
  outTree->Branch("LAD.GEM.hit.xout", &LAD_GEM_hit_xout);
  outTree->Branch("LAD.GEM.hit.yout", &LAD_GEM_hit_yout);
  outTree->Branch("LAD.GEM.hit.zout", &LAD_GEM_hit_zout);
  outTree->Branch("LAD.GEM.hit.t_out", &LAD_GEM_hit_t_out);
  outTree->Branch("LAD.GEM.hit.edep", &LAD_GEM_hit_edep);

  // Loop over entries and print the vectors
  Long64_t nEntries = hodoposition->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    // Clear the vectors
    LAD_Hodo_hit_plane->clear();
    LAD_Hodo_hit_paddle->clear();
    LAD_Hodo_hit_xhit->clear();
    LAD_Hodo_hit_yhit->clear();
    LAD_Hodo_hit_zhit->clear();
    LAD_Hodo_hit_tavg->clear();
    LAD_Hodo_hit_beta->clear();
    LAD_Hodo_hit_sumedep->clear();

    LAD_GEM_hit_plane->clear();
    LAD_GEM_hit_xin->clear();
    LAD_GEM_hit_yin->clear();
    LAD_GEM_hit_zin->clear();
    LAD_GEM_hit_t->clear();
    LAD_GEM_hit_xout->clear();
    LAD_GEM_hit_yout->clear();
    LAD_GEM_hit_zout->clear();
    LAD_GEM_hit_t_out->clear();
    LAD_GEM_hit_edep->clear();

    hodoposition->GetEntry(i);
    hodoenergy->GetEntry(i);
    gemana->GetEntry(i);

    // Loop over the hodo hits
    for (int j = 0; j < vXbar->size(); ++j) {
      if (vTrackID->at(j) != 1) {
        continue;
      }
      LAD_Hodo_hit_xhit->push_back(vXbar->at(j) / 1000);
      LAD_Hodo_hit_yhit->push_back(vYbar->at(j) / 1000);
      LAD_Hodo_hit_zhit->push_back(vZbar->at(j) / 1000);
      LAD_Hodo_hit_tavg->push_back(vTbar->at(j));
      int paddle = vPaddle->at(j);
      if (paddle / 100 == 101) {
        cout << "Paddle: " << paddle << endl;
      }
      switch (paddle / 100) {
      case 000:
        LAD_Hodo_hit_plane->push_back(0);
        break;
      case 001:
        LAD_Hodo_hit_plane->push_back(1);
        break;
      case 100:
        LAD_Hodo_hit_plane->push_back(2);
        break;
      case 101:
        LAD_Hodo_hit_plane->push_back(3);
        break;
      case 200:
        LAD_Hodo_hit_plane->push_back(4);
        break;
      default:
        LAD_Hodo_hit_plane->push_back(-1);
        break;
      }
      paddle = paddle % 100;
      LAD_Hodo_hit_paddle->push_back(paddle);

      // hardcode Edep for now, to not have to mess with digitization
      LAD_Hodo_hit_sumedep->push_back((hodo_maxEdep - hodo_minEdep) * ((double)rand() / RAND_MAX) + hodo_minEdep);
    }

    LAD_Hodo_hit_nhits = LAD_Hodo_hit_xhit->size();

    // Loop over the gem hits
    for (int j = 0; j < vXloc->size(); ++j) {

      if (vgLevel->at(j) != 1 || vgPDG->at(j) != 2212) { // 2212 is the PDG ID for protons
        continue;
      }
      if(j != 0 && vChamber->at(j) == vChamber->at(j-1) && vgPDG->at(j) == vgPDG->at(j-1) && vgLevel->at(j) == vgLevel->at(j-1)) {
        continue;
      }
      LAD_GEM_hit_xin->push_back(vXloc->at(j) / 1000);
      LAD_GEM_hit_yin->push_back(vYloc->at(j) / 1000);
      LAD_GEM_hit_zin->push_back(vZloc->at(j) / 1000);
      LAD_GEM_hit_xout->push_back(vXloc->at(j) / 1000);
      LAD_GEM_hit_yout->push_back(vYloc->at(j) / 1000);
      LAD_GEM_hit_zout->push_back((vZloc->at(j) + 5.0) / 1000);
      LAD_GEM_hit_t->push_back(vTchamber->at(j));
      LAD_GEM_hit_t_out->push_back(vTchamber->at(j));
      LAD_GEM_hit_edep->push_back((gem_maxEdep - gem_minEdep) * ((double)rand() / RAND_MAX) + gem_minEdep);
      LAD_GEM_hit_plane->push_back(vChamber->at(j));
    }
    LAD_GEM_hit_nhits = LAD_GEM_hit_xin->size();
    outTree->Fill();
  }

  // Write the new tree to a new file

  outfile->cd();
  outTree->Write();
  // Close the files
  outfile->Close();
  file->Close();
}