#include <TFile.h>
#include <TTree.h>
#include <vector>

using namespace std;

struct HodoscopeHit {
  int plane;
  double x;
  double y;
};

/*
Variables to write
int Earm.BBGEM.hit.nhits
vector<int> Earm.BBGEM.hit.plane
vector<double> Earm.BBGEM.hit.x
vector<double> Earm.BBGEM.hit.y
vector<double> Earm.BBGEM.hit.z
vector<double> Earm.BBGEM.hit.t
vector<double> Earm.BBGEM.hit.beta
vector<double> Earm.BBGEM.hit.edep
*/
void CreateTree() {
  // Create a ROOT file
  TFile *file = new TFile("lad_hodo_sim.root", "RECREATE");

  // Create a TTree
  TTree *tree = new TTree("T", "Hodoscope Hits");

  // Create branches for hodoscope hits
  int nHits;
  vector<int> plane;
  vector<double> x;
  vector<double> y;
  vector<double> z;
  vector<double> t;
  vector<double> beta;
  vector<double> edep;

  tree->Branch("nHits", &nHits, "nHits/I");
  tree->Branch("plane", &plane);
  tree->Branch("x", &x);
  tree->Branch("y", &y);
  tree->Branch("z", &z);
  tree->Branch("t", &t);
  tree->Branch("beta", &beta);
  tree->Branch("edep", &edep);

  // Fill the tree with some example data
  double nEvents = 20;
  for (int i = 0; i < nEvents; i++) {
    nHits = 3;
    plane.resize(nHits);
    x.resize(nHits);
    y.resize(nHits);
    z.resize(nHits);
    t.resize(nHits);
    beta.resize(nHits);
    edep.resize(nHits);
    for (int j = 0; j < nHits; j++) {
      plane[j] = j + 1;
      x[j]     = 0.1 * j;
      y[j]     = 0.2 * j;
      z[j]     = 0.3 * j;
      t[j]     = 0.4 * j;
      beta[j]  = 0.5;
      // 10^-5 to 10^-2
      edep[j]  = 1/1000 * j;
    }
    tree->Fill();
  }

  // Write the tree to the file and close it
  file->Write();
  file->Close();

  delete file;
}

int main() {
  CreateTree();

  return 0;
}