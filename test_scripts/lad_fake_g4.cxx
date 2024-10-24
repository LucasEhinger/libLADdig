#include <TFile.h>
#include <TTree.h>
#include <vector>

using namespace std;

void CreateTree(int nEvents = 20) {
  // Create a ROOT file
  TFile *file = new TFile("lad_hodo_sim_time_walk.root", "RECREATE");

  // Create a TTree
  TTree *tree = new TTree("T", "Hodoscope Hits");

  // Create branches for hodoscope hits
  int nHits;
  vector<int> plane;
  vector<int> paddle;
  vector<double> xhit;
  vector<double> yhit;
  vector<double> zhit;
  vector<double> t;
  vector<double> beta;
  vector<double> edep;

  // Set min and max values for random generation
  int minPlane   = 0;
  int maxPlane   = 4;
  int minPaddle  = 0;
  int maxPaddle  = 10;
  int minNhits   = 1;
  int maxNhits   = 6;
  double minXhit = -1.0;
  double maxXhit = 1.0;
  double minYhit = -1.0;
  double maxYhit = 1.0;
  double minZhit = -1.0;
  double maxZhit = 1.0;
  double minT    = 0.0;
  double maxT    = 0.0;
  double minBeta = 0.0;
  double maxBeta = 1.0;
  double minEdep = 1. / 1000;
  double maxEdep = 1. / 100;

  string prefix = "LAD.Hodo.hit.";

  tree->Branch((prefix + "nhits").c_str(), &nHits, "nhits/I");
  tree->Branch((prefix + "plane").c_str(), &plane);
  tree->Branch((prefix + "paddle").c_str(), &paddle);
  tree->Branch((prefix + "xhit").c_str(), &xhit);
  tree->Branch((prefix + "yhit").c_str(), &yhit);
  tree->Branch((prefix + "zhit").c_str(), &zhit);
  tree->Branch((prefix + "tavg").c_str(), &t);
  tree->Branch((prefix + "beta").c_str(), &beta);
  tree->Branch((prefix + "sumedep").c_str(), &edep);

  // Fill the tree with some example data
  for (int i = 0; i < nEvents; i++) {
    nHits = rand() % (maxNhits - minNhits + 1) + minNhits;
    plane.resize(nHits);
    paddle.resize(nHits);
    xhit.resize(nHits);
    yhit.resize(nHits);
    zhit.resize(nHits);
    t.resize(nHits);
    beta.resize(nHits);
    edep.resize(nHits);
    // Fill the vectors with random values within the specified range
    for (int j = 0; j < nHits; j++) {
      plane[j]  = rand() % (maxPlane - minPlane + 1) + minPlane;
      paddle[j] = rand() % (maxPaddle - minPaddle + 1) + minPaddle;
      xhit[j]   = (maxXhit - minXhit) * ((double)rand() / RAND_MAX) + minXhit;
      yhit[j]   = (maxYhit - minYhit) * ((double)rand() / RAND_MAX) + minYhit;
      zhit[j]   = (maxZhit - minZhit) * ((double)rand() / RAND_MAX) + minZhit;
      // zhit[j]   = 0.0;
      t[j]      = (maxT - minT) * ((double)rand() / RAND_MAX) + minT;
      beta[j]   = (maxBeta - minBeta) * ((double)rand() / RAND_MAX) + minBeta;
      edep[j]   = (maxEdep - minEdep) * ((double)rand() / RAND_MAX) + minEdep;
    }
    tree->Fill();
  }

  // Write the tree to the file and close it
  file->Write();
  file->Close();

  delete file;
}

int lad_fake_g4(double nEvents = 20) {
  CreateTree(nEvents);

  return 0;
}