#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <string>
#include <vector>
#include <TH1.h>
using namespace std;

int make_test_plots() {
  // Open the ROOT file
  TFile *file = TFile::Open("GMN_SBS8_elastic_job_0_cpy.root");
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file" << std::endl;
    return 1;
  }

  // Get the TTree
  TTree *tree = nullptr;
  file->GetObject("T", tree);
  if (!tree) {
    std::cerr << "Error getting TTree" << std::endl;
    file->Close();
    return 1;
  }

  // Define vectors to hold branch data
  std::vector<double> *x_in  = nullptr;
  std::vector<double> *y_in  = nullptr;
  std::vector<double> *z_in  = nullptr;
  std::vector<double> *x_out = nullptr;
  std::vector<double> *y_out = nullptr;
  std::vector<double> *z_out = nullptr;

  // Set branch addresses
  string prefix = "LAD.GEM.hit.";
  tree->SetBranchAddress((prefix + "xin").c_str(), &x_in);
  tree->SetBranchAddress((prefix + "yin").c_str(), &y_in);
  tree->SetBranchAddress((prefix + "zin").c_str(), &z_in);
  tree->SetBranchAddress((prefix + "xout").c_str(), &x_out);
  tree->SetBranchAddress((prefix + "yout").c_str(), &y_out);
  tree->SetBranchAddress((prefix + "zout").c_str(), &z_out);

  // Create new branches for the differences
  std::vector<double> dx, dy, dz;
  TBranch *b_dx = tree->Branch("dx", &dx);
  TBranch *b_dy = tree->Branch("dy", &dy);
  TBranch *b_dz = tree->Branch("dz", &dz);

  // Loop over the entries in the tree
  Long64_t nentries = tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    tree->GetEntry(i);

    // Calculate differences and fill the new branches
    dx.clear();
    dy.clear();
    dz.clear();
    for (size_t j = 0; j < x_in->size(); ++j) {
      dx.push_back(x_out->at(j) - x_in->at(j));
      dy.push_back(y_out->at(j) - y_in->at(j));
      dz.push_back(z_out->at(j) - z_in->at(j));
    }

    b_dx->Fill();
    b_dy->Fill();
    b_dz->Fill();
  }

  // Write the new branches to the file
  // tree->Write("", TObject::kOverwrite);

  // Create histograms
  double rng = 0.005;
  TH1D *h_dx = new TH1D("h_dx", "dx distribution", 100, -rng, rng);
  TH1D *h_dy = new TH1D("h_dy", "dy distribution", 100, -rng, rng);
  TH1D *h_dz = new TH1D("h_dz", "dz distribution", 100, -rng, rng);

  // Fill histograms
  for (Long64_t i = 0; i < nentries; ++i) {
    tree->GetEntry(i);
    for (size_t j = 0; j < dx.size(); ++j) {
      h_dx->Fill(dx[j]);
      h_dy->Fill(dy[j]);
      h_dz->Fill(dz[j]);
    }
  }

  // Save histograms to a new file
  TFile *outfile = new TFile("histograms.root", "RECREATE");
  h_dx->Write();
  h_dy->Write();
  h_dz->Write();
  outfile->Close();

  // Clean up histograms
  delete h_dx;
  delete h_dy;
  delete h_dz;
  // Clean up
  file->Close();
  delete file;

  return 0;
}