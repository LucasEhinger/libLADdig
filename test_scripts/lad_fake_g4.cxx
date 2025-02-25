#include <TFile.h>
#include <TTree.h>
#include <vector>

using namespace std;

void CreateTree(int nEvents = 20) {
  // Create a ROOT file
  TFile *file = new TFile("lad_hodo_gem_sim.root", "RECREATE");

  // Create a TTree
  TTree *tree = new TTree("T", "LAD Hits");

  // Create branches for hodoscope hits
  int hodo_nHits;
  vector<int> hodo_plane;
  vector<int> hodo_paddle;
  vector<double> hodo_xhit;
  vector<double> hodo_yhit;
  vector<double> hodo_zhit;
  vector<double> hodo_t;
  vector<double> hodo_beta;
  vector<double> hodo_edep;

  // Create branches for GEM hits
  int gem_nHits;
  vector<int> gem_plane;
  vector<int> gem_strip;
  vector<double> gem_xhit_in;
  vector<double> gem_yhit_in;
  vector<double> gem_zhit_in;
  vector<double> gem_t_in;
  vector<double> gem_xhit_out;
  vector<double> gem_yhit_out;
  vector<double> gem_zhit_out;
  vector<double> gem_t_out;
  vector<double> gem_edep;

  // Set min and max values for random generation
  // Hodo
  int hodo_minPlane   = 0;
  int hodo_maxPlane   = 4;
  int hodo_minPaddle  = 0;
  int hodo_maxPaddle  = 10;
  int hodo_minNhits   = 1;
  int hodo_maxNhits   = 6;
  double hodo_minXhit = -1.0;
  double hodo_maxXhit = 1.0;
  double hodo_minYhit = -1.0;
  double hodo_maxYhit = 1.0;
  double hodo_minZhit = -1.0;
  double hodo_maxZhit = 1.0;
  double hodo_minT    = 0.0;
  double hodo_maxT    = 0.0;
  double hodo_minBeta = 0.0;
  double hodo_maxBeta = 1.0;
  double hodo_minEdep = 1. / 1000;
  double hodo_maxEdep = 1. / 100;

  // GEM
  const double gem_width  = 0.4 * 3072 / 1000;
  const double gem_height = 0.4 * 1536 / 1000; // 0.4mm spacing * n strips in m

  int gem_minPlane   = 0; // Plane naming starts at 0
  int gem_maxPlane   = 1; // includes 1
  int gem_minStrip   = 0;
  int gem_maxStrip   = 100;
  int gem_minNhits   = 2;
  int gem_maxNhits   = 2;
  double gem_minXhit = -gem_width/2; // -gem_width/2;
  double gem_maxXhit = gem_width/2; // gem_width/2;
  double gem_minYhit = -gem_height/2; //-gem_height/2;
  double gem_maxYhit = gem_height/2; // gem_height/2;
  double gem_minZhit = 0.0;
  double gem_maxZhit = 0.0;
  double gem_minT    = 10.0;
  double gem_maxT    = 12.0;
  double gem_minEdep = 9e-6;
  double gem_maxEdep = 9e-6; // origeonally 1e-7 to 1e-6

  // Create branches for Hodo hits
  string hodo_prefix = "LAD.Hodo.hit.";

  tree->Branch((hodo_prefix + "nhits").c_str(), &hodo_nHits); //, "nHits/I");
  tree->Branch((hodo_prefix + "plane").c_str(), &hodo_plane);
  tree->Branch((hodo_prefix + "paddle").c_str(), &hodo_paddle);
  tree->Branch((hodo_prefix + "xhit").c_str(), &hodo_xhit);
  tree->Branch((hodo_prefix + "yhit").c_str(), &hodo_yhit);
  tree->Branch((hodo_prefix + "zhit").c_str(), &hodo_zhit);
  tree->Branch((hodo_prefix + "tavg").c_str(), &hodo_t);
  tree->Branch((hodo_prefix + "beta").c_str(), &hodo_beta);
  tree->Branch((hodo_prefix + "sumedep").c_str(), &hodo_edep);

  // Create branches for GEM hits
  // string gem_prefix = "LAD.GEM.hit.";
  string gem_prefix = "LAD.GEM.hit.";
  tree->Branch((gem_prefix + "nhits").c_str(), &gem_nHits); //, "nHits/I");
  tree->Branch((gem_prefix + "plane").c_str(), &gem_plane);
  tree->Branch((gem_prefix + "strip").c_str(), &gem_strip);
  tree->Branch((gem_prefix + "xin").c_str(), &gem_xhit_in);
  tree->Branch((gem_prefix + "yin").c_str(), &gem_yhit_in);
  tree->Branch((gem_prefix + "zin").c_str(), &gem_zhit_in);
  tree->Branch((gem_prefix + "t").c_str(), &gem_t_in);
  tree->Branch((gem_prefix + "xout").c_str(), &gem_xhit_out);
  tree->Branch((gem_prefix + "yout").c_str(), &gem_yhit_out);
  tree->Branch((gem_prefix + "zout").c_str(), &gem_zhit_out);
  tree->Branch((gem_prefix + "t_out").c_str(), &gem_t_out);
  tree->Branch((gem_prefix + "edep").c_str(), &gem_edep);

  // Fill the tree with some example data
  for (int i = 0; i < nEvents; i++) {
    hodo_nHits = rand() % (hodo_maxNhits - hodo_minNhits + 1) + hodo_minNhits;
    hodo_plane.resize(hodo_nHits);
    hodo_paddle.resize(hodo_nHits);
    hodo_xhit.resize(hodo_nHits);
    hodo_yhit.resize(hodo_nHits);
    hodo_zhit.resize(hodo_nHits);
    hodo_t.resize(hodo_nHits);
    hodo_beta.resize(hodo_nHits);
    hodo_edep.resize(hodo_nHits);
    // Fill the vectors with random values within the specified range
    for (int j = 0; j < hodo_nHits; j++) {
      hodo_plane[j]  = rand() % (hodo_maxPlane - hodo_minPlane + 1) + hodo_minPlane;
      hodo_paddle[j] = rand() % (hodo_maxPaddle - hodo_minPaddle + 1) + hodo_minPaddle;
      hodo_xhit[j]   = (hodo_maxXhit - hodo_minXhit) * ((double)rand() / RAND_MAX) + hodo_minXhit;
      hodo_yhit[j]   = (hodo_maxYhit - hodo_minYhit) * ((double)rand() / RAND_MAX) + hodo_minYhit;
      hodo_zhit[j]   = (hodo_maxZhit - hodo_minZhit) * ((double)rand() / RAND_MAX) + hodo_minZhit;
      // hodo_zhit[j]   = 0.0;
      hodo_t[j]    = (hodo_maxT - hodo_minT) * ((double)rand() / RAND_MAX) + hodo_minT;
      hodo_beta[j] = (hodo_maxBeta - hodo_minBeta) * ((double)rand() / RAND_MAX) + hodo_minBeta;
      hodo_edep[j] = (hodo_maxEdep - hodo_minEdep) * ((double)rand() / RAND_MAX) + hodo_minEdep;
    }

    gem_nHits = rand() % (gem_maxNhits - gem_minNhits + 1) + gem_minNhits;
    gem_plane.resize(gem_nHits);
    gem_strip.resize(gem_nHits);
    gem_xhit_in.resize(gem_nHits);
    gem_yhit_in.resize(gem_nHits);
    gem_zhit_in.resize(gem_nHits);
    gem_t_in.resize(gem_nHits);
    gem_xhit_out.resize(gem_nHits);
    gem_yhit_out.resize(gem_nHits);
    gem_zhit_out.resize(gem_nHits);
    gem_t_out.resize(gem_nHits);
    gem_edep.resize(gem_nHits);
    // Fill the vectors with random values within the specified range
    for (int j = 0; j < gem_nHits; j++) {
      // gem_plane[j]    = rand() % (gem_maxPlane - gem_minPlane + 1) + gem_minPlane;
      gem_plane[j]   = (j == 0) ? 0 : 1;
      gem_strip[j]   = rand() % (gem_maxStrip - gem_minStrip + 1) + gem_minStrip;
      gem_xhit_in[j] = (gem_maxXhit - gem_minXhit) * ((double)rand() / RAND_MAX) + gem_minXhit;
      gem_yhit_in[j] = (gem_maxYhit - gem_minYhit) * ((double)rand() / RAND_MAX) + gem_minYhit;
      // cout << "gem_xhit_in[j] = " << gem_xhit_in[j] << endl;
      // cout << "gem_yhit_in[j] = " << gem_yhit_in[j] << endl;
      gem_zhit_in[j]  = (gem_maxZhit - gem_minZhit) * ((double)rand() / RAND_MAX) + gem_minZhit;
      gem_t_in[j]     = (gem_maxT - gem_minT) * ((double)rand() / RAND_MAX) + gem_minT;
      gem_xhit_out[j] = gem_xhit_in[j]; // (gem_maxXhit - gem_minXhit) * ((double)rand() / RAND_MAX) + gem_minXhit;
      gem_yhit_out[j] = gem_yhit_in[j]; // (gem_maxYhit - gem_minYhit) * ((double)rand() / RAND_MAX) + gem_minYhit;
      gem_zhit_out[j] =
          gem_zhit_in[j] + 0.005; // (gem_maxZhit - gem_minZhit) * ((double)rand() / RAND_MAX) + gem_minZhit;
      gem_t_out[j] = gem_t_in[j];// (gem_maxT - gem_minT) * ((double)rand() / RAND_MAX) + gem_minT;
      gem_edep[j]  = (gem_maxEdep - gem_minEdep) * ((double)rand() / RAND_MAX) + gem_minEdep;
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