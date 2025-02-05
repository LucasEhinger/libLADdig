// #ifndef __CINT__

// #include "SBSSpec.h"
// #include "TSBSSimDecoder.h"
// #include "TSolSimFile.h"

// #include "THaInterface.h"
// #include "THaGlobals.h"
// #include "THaTextvars.h"
// #include "THaAnalyzer.h"
// #include "THaDetector.h"

// #include "TSystem.h"
// #include "TList.h"
// #include "TString.h"
// #include "TFile.h"
// #include "TVector3.h"

// #include <iostream>

// #endif

void ReplayMCDigitized(const char* filename = "digitized", 
		       const char* detsuffix = "grinch",//detector suffix: 
		       //"bbgem" for BigBite spectrometer (GMn, GEn, SIDIS);
		       //"FT" for Front Tracker spectrometer (GEp);
		       //"FPP" for Focal Plane Polarimeters (GEp).
		       bool bkgd,// flag to indicate if digitized file includes background or not.
		       Int_t nevent = -1, // number of events to process
		       Int_t nseg = 0, // number of segments
		       bool do_cuts = true )
{
  if( nseg < 0 || nseg > 100 ) {
    cerr << "Invalid number of run segments = " << nseg
	 << ", must be 0-100" << endl;
    return;
  }
  bool do_parts = true;
  if( nseg == 0 ) {
    do_parts = false;
    nseg = 1;
  }
  
  string bg = "nobkgd";
  if(bkgd)bg = "bkgd";
  
  gSystem->Load("../libsbsdig.so");
  gSystem->Load("/u/home/efuchey/g4work/Offline/SBS-offline/libsbs.so");
  //gSystem->Load("libsbs.so");
  
  TSBSDBManager* manager = TSBSDBManager::GetInstance();
  manager->LoadGeneralInfo(Form("db_generalinfo_%s.dat", detsuffix));
  manager->LoadGeoInfo(Form("g4sbs_%s", detsuffix));

  // dde = new TSBSSimDecoder();
  // dde->SetCrateMapName("db_sbssim_cratemap.dat");
  
  THaInterface::SetDecoder( TSBSSimDecoder::Class() );
  
  cout << "Reading " << detsuffix << endl;
  THaApparatus* SBS_BBSpec = new SBSBigBite( "sbs_bb", "SBS / BigBite", 0, true, false, false );
  gHaApps->Add( SBS_BBSpec );
  cout << "Just read " << detsuffix << endl;
  
  SBS_BBSpec->Print("DET");
  
  TString db_prefix = SBS_BBSpec->GetName();
  TString detsuf(detsuffix);
  db_prefix += "."+detsuf;
  cout << db_prefix.Data() << endl;
  gHaTextvars->Add( "DET", db_prefix.Data() );
  gHaTextvars->Add( "APP", SBS_BBSpec->GetName() );

  THaAnalyzer* analyzer = new THaAnalyzer;
  
  TString rootfile(Form("%s_%s_%s", filename, detsuffix, bg.c_str())), infile0(Form("%s_%s_%s", filename, detsuffix, bg.c_str()));
  TString odeffile("sbssim.odef"), cutfile(Form("sbs_%ssim.cuts",detsuffix));
  rootfile.Append("_replayed_new.root");
  analyzer->EnableBenchmarks();
  analyzer->SetOutFile(rootfile);
  analyzer->SetOdefFile(odeffile);
  if( do_cuts ) analyzer->SetCutFile(cutfile);
  analyzer->SetSummaryFile(Form("%s_%s_new.sum", filename, detsuffix));
  analyzer->SetCrateMapFileName("sbssim_cratemap");
  
  //static int Nrun = TMath::Max(nseg,1);
  THaRunBase* run[0];
  TString title0 = "Digitized MC data";
  for( int i=0; i<nseg; ++i ) {
    TString title(title0), infile(infile0);
    if( do_parts ) {
      title.Append(Form(" part %d", i+1));
      infile.Append(Form("_p%d", i+1));
    }
    infile.Append(".root");
    run[i] = new TSBSSimFile(infile,title);
  }
  if( nseg == 1 && nevent > 0 )
    run[0]->SetLastEvent(nevent);

  bool fail = true;
  if( analyzer->Init(run[0]) == 0 ) {
    cout << "initialization successful..." << endl;
    SBSGRINCH* grinch = SBS_BBSpec->GetDetector("grinch");
    cout << "grinch object pointer : " << grinch << endl;
    //grinch->Print("DET");
    
    // Process the runs
    Int_t ret = 0, ntotal = 0;
    for( int i=0; i<nseg && ret >= 0; ++i ) {
      //cout << "processing segment i / nseg : " << i << "/" << nseg << endl;
      if( i>0 )
	run[i]->SetDate(run[0]->GetDate());
      //cout << "Start processing " << endl;
      ret = analyzer->Process(run[i]);
      cout << "ret = " << ret << endl;
      if( ret > 0 )
	ntotal += ret;
    }
    fail = (ret < 0 );
    if( fail )
      cerr << "Terminated with analyzer error = " << ret << endl;
    else
      cout << "Analyzed " << ntotal << " events" << endl;
    analyzer->Close();
  }
  
  for( int i=0; i<nseg; ++i ) {
    delete run[i]; run[i] = 0;
  }
  delete analyzer; analyzer = 0;
  gHaApps->Delete();
  //}
  
  //TFile* f =
  if( !fail )
    new TFile(rootfile,"r");
} 
