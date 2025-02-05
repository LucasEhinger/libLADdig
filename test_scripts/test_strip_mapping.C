#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <vector>

const int fixed_strip = 98;

void test_strip_mapping() {
  // Open the input file
  TFile *inputFile = TFile::Open("lad_hodo_gem_sim.root", "READ");
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Error opening input file" << std::endl;
    return;
  }

  // Get the tree from the input file
  TTree *inputTree = (TTree*)inputFile->Get("T");
  if (!inputTree) {
    std::cerr << "Error getting tree from input file" << std::endl;
    inputFile->Close();
    return;
  }

  // Create a new file to save the modified tree
  TFile *outputFile = TFile::Open("lad_hodo_gem_sim_strip_test.root", "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "Error creating output file" << std::endl;
    inputFile->Close();
    return;
  }

  // Clone the tree structure
  TTree *outputTree = inputTree->CloneTree(0);

  // Set up the branch to be modified
  int nStrips = 0;
  inputTree->SetBranchAddress("LAD.GEM.dighit.nstrips", &nStrips);
  std::vector<int> *branchValue = nullptr;
  inputTree->SetBranchAddress("LAD.GEM.dighit.strip", &branchValue);

  // Create a new branch in the output tree
  outputTree->Branch("LAD.GEM.dighit.strip", &branchValue);

  // Loop over the entries and modify the branch
  Long64_t nEntries = inputTree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    inputTree->GetEntry(i);

    // Modify the branch value (example modification)
    for (int j = 0; j < nStrips; ++j){
      branchValue->at(j) = fixed_strip;
    }

    // Fill the output tree with the modified entry
    outputTree->Fill();
  }

  cout << "test 1" << endl;
  // Write the output tree to the new file
  outputTree->Write();

  // Close the files
  inputFile->Close();
  outputFile->Close();
}
