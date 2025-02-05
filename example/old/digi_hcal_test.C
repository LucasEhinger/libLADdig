// Example "replay" script
//#define DEBUG 1
#include "TSystem.h"
#include "TDatime.h"
#include "TSBSGeant4File.h"
#include "TSBSSimHCal.h"
#include "TSBSDBManager.h"
#include "TSBSSimDigitizer.h"


void digi_hcal_test(int runnum = 931, int nentries = 100, int debuglevel = 1)

{
  printf("\n** This gets called with 'analyzer' and not 'root' **\n");
  printf("** If you're getting missing symbol errors, this is likely the cause **\n\n");

  TDatime run_time = 991231;

  gSystem->AddDynamicPath("${SBS_ANALYSIS}");
  gSystem->Load("../libsbsdig.so");

  ////////////////////////////////////////////////////////////////

  TSBSDBManager* manager = TSBSDBManager::GetInstance();
  manager->SetDebug(debuglevel);
  //manager->LoadGeneralInfo(Form("%s/db_generalinfo_grinch.dat",gSystem->Getenv("DB_DIR")));
  manager->LoadGenInfo("db_geninfo_hcal.dat");
  //manager->LoadGeoInfo("g4sbs_grinch");


  // Create the SBS Digitizer (will control the digitization process)
  TSBSSimDigitizer *digitizer = new TSBSSimDigitizer(
      Form("digitized/simdig_%d.root",runnum));

  // First load the input root file
  TSBSGeant4File *f = new TSBSGeant4File(Form("data/sbsin_%d.root",runnum));
  
  TSBSSimHCal *hcal = new TSBSSimHCal("hcal",0);
  digitizer->AddDetector(hcal);

  // Now start the digitization on this file
  digitizer->Process(f,nentries);


  delete hcal;
}
