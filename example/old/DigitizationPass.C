// Example "replay" script
//#define DEBUG 1

void DigitizationPass(int fspec = 1, // Spectrometer flag: 
		      // 1 for GRINCH
		      // 2 for RICH
		      UInt_t NsigFiles = 1,//Number of signal files to analize
		      UInt_t Nmax = -1, //number of events to digitize
		      UInt_t nbacktoadd = 2, // number of background *files* to add to each event
		      bool testdis = false,// boolean to digitize DIS test files instead of regular signal.
		      bool print = false){
  printf("\n** This gets called with 'analyzer' and not 'root' **\n");
  printf("** If you're getting missing symbol errors, this is likely the cause **\n\n");
  
  TDatime run_time = 991231;
  
  gSystem->Load("../libsbsdig.so");
  
  ////////////////////////////////////////////////////////////////
  
  int Ngood = 0;
      
  TSBSCher *ddy;
  TSBSSpec *dds;
  TSBSSimCherDigitization *ddd;
  
  char* outname;
  string bg = "bkgd";
  if(nbacktoadd==0)bg = "nobkgd";

  string infile_sig_prefix;
  string infile_bkgd_prefix;
  
  TSBSDBManager* manager = TSBSDBManager::GetInstance();
  switch(fspec){
  case(1):
    manager->LoadGeneralInfo("db_generalinfo_grinch.dat");
    manager->LoadGeoInfo("g4sbs_grinch");
    dds = new TSBSSpec ("g4sbs_grinch", "BB spectrometer");
    outname = Form("digitized_grinch_%s.root", bg.c_str());
    infile_sig_prefix = "/volatile/halla/sbs/efuchey/gmn_elastic/gmn13.5_elastic_sig_20171211_08";
    if(testdis)infile_sig_prefix = "/volatile/halla/sbs/efuchey/misc/test_grinch_20170727_09";//"grinch_signal_0.root";
    infile_bkgd_prefix = "/volatile/halla/sbs/efuchey/gmn_beam_bkgd/gmn13.5_beam_bkgd_20170630_14";
    dds->Init(run_time);
    break;
    //case(2):
    //break;
  default:
    cout << "No corresponding geometry; choose: " << endl 
	 << "1 (GRINCH)" << endl << "2 (RICH - not implemented -> quit too) " << endl;
    return;
    break;
  }
  
  for(int i = 0; i<manager->GetNDetectors(); i++){
    ddy = new TSBSCher (Form("cher%d", i), Form("Cherenkov detector %d", i));
    ddy->SetApparatus(dds);
    if( ddy->Init() )
      return;
    dds->AddCher(ddy);
  }
  
  printf("\n");
  
  //cout << outname << " " << &outname << endl;
  
  if(print)dds->Print();
    
  ddd = new TSBSSimCherDigitization (*dds,"ratedig");
  ////////////////////////////////////////////////////////////////
    
  int nevent = 1;
  
  TSBSGeant4File *f;
  
  int  ndata, i;
  TSBSCherData *chd, *chb;
  g4sbsgendata *gen;
  
  cout << "creating file " << outname << endl;
  ddd->InitTree (*dds, outname);
  
  printf("Digitizing events\n");
  
  //Nmax = TMath::Min((Long64_t)Nmax, f->GetEntries());
    
  int hadback = 1;
  
  int N_bg_file_g = 0;
  
  //Add the loop on the signal files
  for(int i_sig = 0; i_sig<NsigFiles; i_sig++){
    if(testdis){
      f = new TSBSGeant4File(Form("%s/gc_signal_%d.root",infile_sig_prefix.c_str(), i_sig));
    }else{
      f = new TSBSGeant4File(Form("%s/elastic_%d.root",infile_sig_prefix.c_str(), i_sig));
    }
    printf("The filename returned is %s\n", f->GetFileName());
    f->SetSource(0);
    
    int res;
    
    res = f->Open();
    
    if( res != 1 ){
      printf("Opening g4sbs file returned %d\n", res);
      return;
    }
    
    int d_flag_readevent = 0;
    while( f->ReadNextEvent(d_flag_readevent) && hadback && nevent<Nmax ){
      
      ndata = 0;
      if(nevent%100==0){
	cout << "Evt " << nevent << endl;
      }
      
      if(f->GetNData()==0){
	if(print)
	  cout << "No hits, skip evt " << nevent << endl;
	nevent++;
	continue;
      }
      
      chd = f->GetCherData();
      if(f->GetNGen()>0){
	gen = f->GetGenData(0);
	Ngood++;
      }else{
	if(print)
	  cout << "No generated data for event " << nevent 
	       << ", skip it (Nhits = " << f->GetNData() << ")" << endl;
	nevent++;
	continue;
      }
      
      if(print)
	cout << "Evt " << nevent << " has hits and generated data " << endl;
      
      ddd->SetTreeEvent((*chd), (*f), nevent);
    
      if(print){
	cout << "number of hits in GRINCH data " << chd->GetNHit() << endl;
	while(ndata<chd->GetNHit()){
	
	  //if(chd->GetParticleID(ndata)>1)continue;
	  chd->Print();
	  cout << "hit number " << ndata << endl;
	  chd->PrintHit(ndata);
	  ndata++;
	}
      }
      ddd->Digitize(*chd, *dds);
      
      // Access to generated vertex and momentum
      // gen->GetV();
      // gen->GetP();
      // gen->GetWeight();
    

      // Add some number of background files...
      int N_bg_file_g_post = N_bg_file_g+nbacktoadd;
 
      if(nbacktoadd){
	for(int Nfile = N_bg_file_g; Nfile < N_bg_file_g_post; Nfile++){
	  //if(print)cout << N_bg_file_g << " <= " << Nfile << " < " << N_bg_file_g+nbacktoadd << endl;
	  TSBSGeant4File *fback = new TSBSGeant4File(Form("%s/beam_bkgd_%d.root",infile_bkgd_prefix.c_str(), Nfile));
	  int open = fback->Open();
	  if(!open){
	    N_bg_file_g_post++;
	    N_bg_file_g++;
	  
	    if(N_bg_file_g>=2000){
	      int n_temp = Nfile;
	      Nfile = N_bg_file_g_post-n_temp;
	      N_bg_file_g_post = nbacktoadd;
	      N_bg_file_g = 0;
	    }
	    //if(print)cout << Form("/group/exjpsi/eric/31722/beam_bkgd_%d.root does not exist", Nfile) << endl;
	    continue;
	  }
	  
	  fback->SetSource(1);
	  //if(print)
	  printf("The background file read is %s\n", fback->GetFileName());
	    
	  int backidx = 0;
	  //while( hadback = fback->ReadNextEvent() && backidx < nbacktoadd ){
	  while( backidx < fback->GetEntries() ){
	    hadback = fback->ReadNextEvent();
	    chb = fback->GetCherData();
	    if(print && chb->GetNHit()>0){
	      cout << "Bkgd evt: " << chb->GetEvent() << ", number of hits " 
		   << chb->GetNHit() << endl;
	      int nback = 0;
	      while(nback<chb->GetNHit()){
		//if(chd->GetParticleID(ndata)>1)continue;
		chb->Print();
		cout << "hit number " << nback << endl;
		chb->PrintHit(nback);
		nback++;
	      }
	    
	    }
	  
	    ddd->AdditiveDigitize(*chb, *dds);
	    
	    // //Randomize times based on gate width
	    // for( int bidx = 0; bidx < chb->GetNHit(); bidx++ ){
	    //   double timeshift = gRandom->Uniform(-ddd->GetGateWidth(), 75.0 );//ns
	    //   chb->SetHitTime(bidx, chb->GetHitTime(bidx) + timeshift );
	    // }	
	    // //chd->AddGEMData(chb);
	    backidx++;
	  }
	  
	  // if( backidx != nbacktoadd ){
	  // printf("Warning:  Not enough background events to be added (%d)\n", backidx);
	  // }
	  
	  fback->Close();
	}
	//if(print)cout << "new number of hits in GEM data " << chd->GetNHit() << endl;
	N_bg_file_g = N_bg_file_g_post;
      }//end if nbacktoadd
    
      if(N_bg_file_g>=2000)N_bg_file_g = 0;
      
      ddd->FillTree();
      
      //if(nevent==7)ddd->GetEvent()->Print("all");
      if(print)ddd->GetEvent()->Print("all");
      //ddd->GetEvent()->Print("clust");
      
      delete chd;
      nevent++;
    }
    
  }
  
  printf("Completed %d events total: %d good events \n", nevent, Ngood);

  ddd->WriteTree();
  ddd->CloseTree();
  
  cout << "Tree closed" << endl;
}
