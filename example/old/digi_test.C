// Example "replay" script
//#define DEBUG 1
#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "TSystem.h"
#include "TDatime.h"
#include "TSBSGeant4File.h"
#include "TSBSSimHCal.h"
#include "TSBSSimECal.h"
#include "TSBSSimScint.h"
#include "TSBSSimCher.h"
#include "TSBSSimGEM.h"
#include "TSBSDBManager.h"
#include "TSBSSimDigitizer.h"
#include "THaAnalysisObject.h"
#endif

//simple, ready to run macro script with libsbsdig
void digi_test(ULong64_t nentries = 100, int nbkgd = 0, int debuglevel = 1)
{
  printf("\n** This gets called with 'analyzer' and not 'root' **\n");
  printf("** If you're getting missing symbol errors, this is likely the cause **\n\n");

  TDatime run_time = 991231;

  gSystem->Load("../libsbsdig.so");

  ////////////////////////////////////////////////////////////////
  
  TSBSDBManager* manager = TSBSDBManager::GetInstance();
  manager->SetDebug(debuglevel);
  
  if(debuglevel>=1)cout << "About to read database " << endl;

  manager->LoadGenInfo("db_geninfo_gmn.dat");
  
  if(debuglevel>=1)cout << "Setup digitizer " << endl;
  
  // Create the SBS Digitizer (will control the digitization process)
  TSBSSimDigitizer *digitizer = new TSBSSimDigitizer("simdig_test.root");
  digitizer->SetDebug(debuglevel);
  
  if(debuglevel>=1)cout << "Setup input file " << endl;
  
  digitizer->AddInputFile("/volatile/halla/sbs/efuchey/gmn13.5_elastic_20200214_10/elastic_0.root", 0, 1);
  
  ifstream inputfile("BeamBkgd_GMn13.5.txt");
  if(nbkgd && inputfile.good()){
    //TFile *fout = new TFile( outputfilename, "RECREATE" );
    TString currentline;
    while( currentline.ReadLine(inputfile) && !currentline.BeginsWith("endlist") ){
      if( !currentline.BeginsWith("#") ){
	digitizer->AddInputFile(currentline.Data(), 1, -nbkgd);
      }
    }
  }
  
  // It is recommended  to declare the detector with its unique ID (second parameter)
  // See list of unique det IDs defined in src/g4sbs_types.h
  
  if(debuglevel>=1)cout << "Declare detectors and add them to digitizer " << endl;
  
  TSBSSimHCal *hcal = new TSBSSimHCal("hcal", 0);
  hcal->SetDebug(debuglevel);
  digitizer->AddDetector(hcal);
  
  TSBSSimScint *hodo = new TSBSSimScint("hodo", 30);
  hodo->SetDebug(debuglevel);
  digitizer->AddDetector(hodo);
  
  TSBSSimCher *grinch = new TSBSSimCher("grinch", 20);
  grinch->SetDebug(debuglevel);
  digitizer->AddDetector(grinch);
  
  TSBSSimECal *ps = new TSBSSimECal("ps", 10);
  ps->SetDebug(debuglevel);
  digitizer->AddDetector(ps);
  
  TSBSSimECal *sh = new TSBSSimECal("sh", 11);
  sh->SetDebug(debuglevel);
  digitizer->AddDetector(sh);
  
  TSBSSimGEM *bbgem = new TSBSSimGEM("gem", 40);
  bbgem->SetDebug(debuglevel);
  digitizer->AddDetector(bbgem);
  
  if(debuglevel>=1)cout << "About to process digitization for " << nentries << "events " << endl;
  
  digitizer->Process(nentries);
}
