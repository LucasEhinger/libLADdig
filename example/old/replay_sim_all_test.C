
#include "TSystem.h"
#include "TList.h"
#include "THaRun.h"
#include "THaEvent.h"
#include "THaEvtTypeHandler.h"
#include "THaAnalyzer.h"
#include "THaApparatus.h"
#include "THaVarList.h"

//#include "SBSGEMStand.h"
//#include "SBSBigBite.h"
#include "SBSEArm.h"
#include "SBSCalorimeter.h"
#include "SBSCDet.h"
#include "SBSTimingHodoscope.h"
#include "SBSGRINCH.h"
#include "TSBSDBManager.h"
#include "THaInterface.h"
#include "TSBSSimDecoder.h"
#include "TSBSSimFile.h"

//R__ADD_LIBRARY_PATH($SBS_ANALYSIS)
//R__ADD_LIBRARY_PATH(../)
//R__LOAD_LIBRARY(libsbs.so)
//R__LOAD_LIBRARY(libsbsdig.so)

void replay_sim_all_test(Int_t runnum = 931, Int_t lastEvent = -1){
  // Set up the analyzer - we use the standard one,
  // but this could be an experiment-specific one as well.
  // The Analyzer controls the reading of the data, executes
  // tests/cuts, loops over Apparatus's and PhysicsModules,
  // and executes the output routines.
  THaAnalyzer* analyzer = new THaAnalyzer;

  TSBSDBManager* manager = TSBSDBManager::GetInstance();
  manager->LoadGenInfo("../db/db_geninfo_gmn.dat");
  //manager->LoadROCMap("../db/db_rocmap.dat");
  THaInterface::SetDecoder( TSBSSimDecoder::Class() );



  // A simple event class to be output to the resulting tree.
  // Creating your own descendant of THaEvent is one way of
  // defining and controlling the output.
  THaEvent* event = new THaEvent;

  // Define the run(s) that we want to analyze.
  // We just set up one, but this could be many.
  //  THaRun* run = new THaRun( "prod12_4100V_TrigRate25_4.dat" );
  //THaRun* run = new THaRun(TString::Format("digitized/simin_%d.root",runnum) );
  THaRunBase *run = new TSBSSimFile( "digitized/simdig_test.root" );
  run->SetFirstEvent(0);
  run->SetLastEvent(lastEvent);

  run->SetDataRequired(0);
  run->SetDate(TDatime());

  analyzer->SetVerbosity(100);

  // Define the analysis parameters
  //analyzer->SetEvent( event );
  analyzer->SetOutFile( TString::Format("rootfiles/simout_%d.root",runnum));
  // File to record cuts accounting information
  analyzer->SetSummaryFile("sbs_hcal_test.log"); // optional

  // Change the cratemap to point to the sim one
  analyzer->SetCrateMapFileName("db_sbssim_cratemap");


  SBSCalorimeter *hcal = new SBSCalorimeter("hcal","HCAL");
  hcal->SetWithADCSamples(true);
  hcal->SetWithTDC(true);

  //SBSCDet *cdet = new SBSCDet("cdet","Coordinate Detector");

  SBSEArm *harm = new SBSEArm("sbs","Hadron Arm with HCal");
  harm->AddDetector(hcal);
  //harm->AddDetector(cdet);
  gHaApps->Add(harm);


  SBSCalorimeter *sh = new SBSCalorimeter("sh","BB shower");
  //SBSTimingHodoscope *hodo = new SBSTimingHodoscope("hodo","BB timing hodoscope");
  SBSCalorimeter *ps = new SBSCalorimeter("ps","BB preshower");
  SBSGRINCH *grinch = new SBSGRINCH("grinch","GRINCH");
  SBSEArm *earm = new SBSEArm("bb","Big Bite electron arm");
  earm->AddDetector(sh);
  //earm->AddDetector(hodo);
  earm->AddDetector(ps);
  earm->AddDetector(grinch);
  gHaApps->Add(earm);
  
  analyzer->SetOdefFile("output_hcal_test.def");
  
  TIter next(gHaEvtHandlers);
  THaEvtTypeHandler *obj = 0;
  int ootnum = 0;
  while( (obj = (THaEvtTypeHandler*)next()) ) {
  //for(THaEvtTypeHandler* obj: *gHaEvtHandlers) {
    obj->SetDebugFile(TString::Format("oooot_debug_%d.log",ootnum));
    ootnum++;
    std::cout << "Here!!" << std::endl;
  }


  //gHaEvtHandlers->SetDebugFile("oooot_debug.log");
  std::cout << "About to initialize Analyzer." << std::endl;
  if( analyzer->Init(run) == 0 ) {
    std::cout << "Analyzer Initialized successfully initialized."
     << std::endl << std::endl;
    //analyzer->SetCompressionLevel(0); // turn off compression
  //gHaEvtHandlers->GetListOfPrimitives()->R__FOR_EACH(THaEvtTypeHandler,
  //    SetDebugFile)("oooot_debug.log");


    std::cout << "Analyzer Processing run:" << std::endl;
    int ret = analyzer->Process(run);     // start the actual analysis
    std::cout << "ret = " << ret << std::endl;
    bool fail = (ret < 0 );
    if( fail ) {
      std::cerr << "Terminated with analyzer error = " << ret << std::endl;
    } else {
      std::cout << "Analyzed " << ret << " events" << std::endl;
    }
    analyzer->Close();
  }

  std::cout << std::endl << std::endl;
  std::cout << "Finished processing. All done now." << std::endl;
}
