// #define g4sbs_tree_cxx
#include "g4sbs_tree.h"
#include <TCanvas.h>
#include <TH2.h>
#include <TStyle.h>
#include <iostream>

// g4sbs_tree constructor: the tree will be the
// the boolean is a flag to consider(true) or ignore(false) the ECal_box and HCal_box data
// g4sbs_tree::g4sbs_tree(TTree *tree, Exp_t expt, bool pythia)
//, bool ecalbox, bool have_hcalbox)
g4sbs_tree::g4sbs_tree(TTree *tree, std::vector<TString> det_list, bool sig_br) : fChain(0) {
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject(
        "/volatile/halla/sbs/efuchey/gmn13.5_elastic_20200228_17/elastic_0.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("/volatile/halla/sbs/efuchey/gmn13.5_elastic_20200228_17/elastic_0.root");
    }
    f->GetObject("T", tree);
  }
  // fExpt = expt;
  // fPythia = pythia;
  // fEcalBox = ecalbox;
  // fHcalBox = have_hcalbox;
  Init(tree, det_list, sig_br);
}

// default destructor
g4sbs_tree::~g4sbs_tree() {
  if (!fChain)
    return;
  delete fChain->GetCurrentFile();
}

// overload of the TTree::GetEntries() function
Int_t g4sbs_tree::GetEntries() {
  // Read contents of entry.
  if (!fChain)
    return 0;
  return fChain->GetEntries();
}
// overload of the TTree::GetEntry(Long64_t) function
Int_t g4sbs_tree::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain)
    return 0;
  return fChain->GetEntry(entry);
}
Long64_t g4sbs_tree::LoadTree(Long64_t entry) {
  // Set the environment to read one entry
  if (!fChain)
    return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0)
    return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void g4sbs_tree::Init(TTree *tree, std::vector<TString> det_list, bool sig_br) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer

  Primaries_PID     = 0;
  Primaries_genflag = 0;
  Primaries_Px      = 0;
  Primaries_Py      = 0;
  Primaries_Pz      = 0;
  Primaries_vx      = 0;
  Primaries_vy      = 0;
  Primaries_vz      = 0;
  Primaries_M       = 0;
  Primaries_E       = 0;
  Primaries_P       = 0;
  Primaries_t       = 0;
  Primaries_theta   = 0;
  Primaries_phi     = 0;

  // Set branch addresses and branch pointers
  if (!tree)
    return;
  fChain   = tree;
  fCurrent = -1;
  // (jc2): Why do we want to make a class again??
  // I disabled this so that we can force the tree to check each
  // SetBranchStatus for matches.
  // fChain->SetMakeClass(1);

  // Setup "Event branch": can be useful
  // fChain->SetBranchAddress("ev", &ev_count, &b_ev);

  for (int k = 0; k < det_list.size(); k++) {
    // GMN/GEN
    if (det_list[k] == "bbhodo") {
      printf(" bbhodo branches set up! \n");
      SetupDetBranch(Earm_BBHodoScint, "LAD.Hodo.hit");
      SetupDetBranch(Earm_BBHodo_Dig, "LAD.Hodo.dighit");
    }
    if (det_list[k] == "bbgem") {
      printf(" bbgem branches set up! \n");
      SetupDetBranch(Earm_BBGEM, "Earm.BBGEM.hit");
      SetupDetBranch(Earm_BBGEM_Dig, "Earm.BBGEM.dighit");
      // if(sig_br)SetupDetBranch(Earm_BBGEM_Dig_sig, "Earm.BBGEM.dighit_sig");
    }
    if (det_list[k] == "sbsgem") {
      printf(" fpp2 branches set up! \n");
      SetupDetBranch(Harm_SBSGEM, "Harm.SBSGEM.hit");
      SetupDetBranch(Harm_SBSGEM_Dig, "Harm.SBSGEM.dighit");
    }
  }

  Notify();
}

Bool_t g4sbs_tree::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void g4sbs_tree::Show(Long64_t entry) {
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain)
    return;
  fChain->Show(entry);
}

Int_t g4sbs_tree::Cut(Long64_t entry) {
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

void g4sbs_tree::Loop() {
  //   In a ROOT session, you can do:
  //      Root > .L gep_tree_with_spin.C
  //      Root > gep_tree_with_spin t
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //

  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  // by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0)
    return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    // if (Cut(ientry) < 0) continue;
  }
}

void g4sbs_tree::SetupDetBranch(TSBSGeant4::VDetData_t &det, const char *prefix) { det.SetupBranches(fChain, prefix); }

void g4sbs_tree::ClearDigBranches() {
  Earm_BBGEM_Dig.ClearBranches();
  Earm_BBHodo_Dig.ClearBranches();
  Harm_SBSGEM_Dig.ClearBranches();
}

void g4sbs_tree::FillDigBranches() {
  Earm_BBGEM_Dig.FillBranches();
  Earm_BBHodo_Dig.FillBranches();
  Harm_SBSGEM_Dig.FillBranches();
}
