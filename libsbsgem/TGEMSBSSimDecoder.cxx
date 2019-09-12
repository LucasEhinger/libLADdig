//*-- Author :    Ole Hansen (ole@jlab.org)    9-Dec-2011

/////////////////////////////////////////////////////////////////////
//
//   TGEMSBSSimDecoder
//
//   Decoder for SoLID simulation data
//
//   Interprets event buffer from input as TGEMSBSSimEvent objects
//   (containing digitized simulation data) and unpacks them into
//   crateslot arrays for low-level decoding by detectors.
//
/////////////////////////////////////////////////////////////////////

#include "TGEMSBSSimDecoder.h"
#include "THaCrateMap.h"
#include "THaBenchmark.h"
#include "VarDef.h"
#include "TGEMSBSDBManager.h"
#include "ha_compiledata.h"

#include "TError.h"
#include "TSystem.h"
#include "TMath.h"
#include "TDatabasePDG.h"
#include "TRandom.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <utility>
#include <stdexcept>

using namespace std;
using namespace Podd;

//EFuchey: 2016/12/10: it is necessary to declare the TGEMSBSDBManager as a static instance here 
// (and not make it a member) because it is used by functions whic are defined as "static inline".
static TGEMSBSDBManager* fManager = TGEMSBSDBManager::GetInstance();
static const Int_t kPrimaryType = 1, kPrimarySource = 0;
// Projection types must match the definitions in TreeSearch
enum EProjType { kUPlane = 0, kVPlane =1, kXPlane = 2, kYPlane = 3};

typedef vector<int>::size_type vsiz_t;

// TODO: Have dbconvert write out MAXSLOT (and possibly other parameters)
//  to the database to allow client to understand the generated detector maps.
// FIXME: The number 30 is hardcoded in dbconvert
static const Int_t SIM_MAXSLOT = TMath::Min(Decoder::MAXSLOT,30);

//-----------------------------------------------------------------------------
TGEMSBSSimDecoder::TGEMSBSSimDecoder()
{
  // Constructor

  fMCHits     = new TClonesArray( "TGEMSBSSimGEMHit",    200 );
  fMCTracks   = new TClonesArray( "TGEMSBSSimTrack",       1 );
  fBackTracks = new TClonesArray( "TGEMSBSSimBackTrack",   5 );
  
  DefineVariables();

  gSystem->Load("libEG.so");  // for TDatabasePDG
  
  for (int i=0; i<fManager->GetNSigParticle(); i++){
    fSignalInfo.push_back(SignalInfo(fManager->GetSigPID(i),
                                     fManager->GetSigTID(i)));
  }
}

//-----------------------------------------------------------------------------
TGEMSBSSimDecoder::~TGEMSBSSimDecoder() {

  DefineVariables( THaAnalysisObject::kDelete );

  delete fBackTracks;
  // fMCHits and fMCTracks are deleted by SimDecoder destructor
}

//-----------------------------------------------------------------------------
Int_t TGEMSBSSimDecoder::DefineVariables( THaAnalysisObject::EMode mode )
{
  // Define global variables for the MC quantities. Extends the base
  // class method.

  const char* const here = "TGEMSBSSimDecoder::DefineVariables";
  
  if( mode == THaAnalysisObject::kDefine && fIsSetup )
    return THaAnalysisObject::kOK;
  
  SimDecoder::DefineVariables( mode );

  RVarDef vars[] = {
    // Generated track info
    //{ "tr.n",      "Number of tracks",      "GetNMCTracks()" },  // already defined in Podd::SimDecoder
    { "tr.vx",     "Track origin x (m)",    "fMCTracks.TGEMSBSSimTrack.VX()" },
    { "tr.vy",     "Track origin y (m)",    "fMCTracks.TGEMSBSSimTrack.VY()" },
    { "tr.vz",     "Track origin z (m)",    "fMCTracks.TGEMSBSSimTrack.VZ()" },
    { "tr.p",      "Track momentum (GeV)",  "fMCTracks.TGEMSBSSimTrack.P() "},
    { "tr.theta",  "Track theta_p (rad)",   "fMCTracks.TGEMSBSSimTrack.PTheta()" },
    { "tr.phi",    "Track phi_p (rad)",     "fMCTracks.TGEMSBSSimTrack.PPhi()" },
    { "tr.pid",    "Track PID (PDG)",       "fMCTracks.TGEMSBSSimTrack.fPID" },
    { "tr.num",    "GEANT track number",    "fMCTracks.TGEMSBSSimTrack.fNumber" },
    { "tr.planes", "Bitpattern of planes hit", "fMCTracks.TGEMSBSSimTrack.fHitBits" },
    { "tr.nhits",  "Number of tracker hits","fMCTracks.TGEMSBSSimTrack.fNHits" },
    { "tr.nfound", "Number of hits found",  "fMCTracks.TGEMSBSSimTrack.fNHitsFound" },
    { "tr.flags",  "Reconstruction status", "fMCTracks.TGEMSBSSimTrack.fReconFlags" },

    // Results of fit to MC points - measures multiple scattering
    // Those below are not useful for SBS, which needs X, Y, Xdir, Ydir (unless otherwise demonstrated)
    // refer to comment in TGEMSBSSimEvent.h l. 30-32
    // { "tr.mcfit.r",     "Track x from MC fit [m]", "fMCTracks.TGEMSBSSimTrack.MCFitR()" },
    // { "tr.mcfit.phi",   "Track phi from MC fit [rad]", "fMCTracks.TGEMSBSSimTrack.MCFitPhi()" },
    // { "tr.mcfit.thdir", "Track dir theta from MC fit [rad]", "fMCTracks.TGEMSBSSimTrack.MCFitThetaDir()" },
    // { "tr.mcfit.phdir", "Track x from MC fit [rad]", "fMCTracks.TGEMSBSSimTrack.MCFitPhiDir()" },
    // { "tr.mcfit.x",     "Track x from MC fit [m]",       "fMCTracks.TGEMSBSSimTrack.MCFitX_print()" },
    { "tr.mcfit.x",     "Track x from MC fit [m]",       "fMCTracks.TGEMSBSSimTrack.fMCFitPar[0]" },
    { "tr.mcfit.xdir",  "Track dir x from MC fit [rad]", "fMCTracks.TGEMSBSSimTrack.fMCFitPar[1]" },
    { "tr.mcfit.y",     "Track y from MC fit [rad]",     "fMCTracks.TGEMSBSSimTrack.fMCFitPar[2]" },
    { "tr.mcfit.ydir",  "Track dir y from MC fit [rad]", "fMCTracks.TGEMSBSSimTrack.fMCFitPar[3]" },
    { "tr.mcfit.chi2",  "Chi2 of MC fit",                "fMCTracks.TGEMSBSSimTrack.fMCFitPar[4]" },
    { "tr.mcfit.ndof",  "NDoF of MC fit",                "fMCTracks.TGEMSBSSimTrack.fMCFitPar[5]" },
    { "tr.mcfit.vx",    "Vertex x from MC fit [m]",      "fMCTracks.TGEMSBSSimTrack.fMCFitPar[6]" },
    { "tr.mcfit.vy",    "Vertex y from MC fit [m]",      "fMCTracks.TGEMSBSSimTrack.fMCFitPar[7]" },
    { "tr.mcfit.vz",    "Vertex z from MC fit [m]",      "fMCTracks.TGEMSBSSimTrack.fMCFitPar[8]" },

    // Results of fit to reconstructed MC hits - checks hit resolution effects
    // independent of track finding
    // Those below are not useful for SBS, which needs X, Y, Xdir, Ydir (unless otherwise demonstrated)
    // refer to comment in TGEMSBSSimEvent.h l. 30-32
    // { "tr.fit.r",     "Track x from rec hit fit [m]", "fMCTracks.TGEMSBSSimTrack.RcFitR()" },
    // { "tr.fit.phi",   "Track phi from rec hit fit [rad]", "fMCTracks.TGEMSBSSimTrack.RcFitPhi()" },
    // { "tr.fit.thdir", "Track dir theta from rec hit fit [rad]", "fMCTracks.TGEMSBSSimTrack.RcFitThetaDir()" },
    // { "tr.fit.phdir", "Track x from rec hit fit [rad]", "fMCTracks.TGEMSBSSimTrack.RcFitPhiDir()" },
    { "tr.fit.x",     "Track x from rec hit fit [m]",       "fMCTracks.TGEMSBSSimTrack.fRcFitPar[0]" },
    { "tr.fit.xdir",  "Track dir x from rec hit fit [rad]", "fMCTracks.TGEMSBSSimTrack.fRcFitPar[1]" },
    { "tr.fit.y",     "Track y from rec hit fit [rad]",     "fMCTracks.TGEMSBSSimTrack.fRcFitPar[2]" },
    { "tr.fit.ydir",  "Track dir y from rec hit fit [rad]", "fMCTracks.TGEMSBSSimTrack.fRcFitPar[3]" },
    { "tr.fit.chi2",  "Chi2 of rec hit fit",                "fMCTracks.TGEMSBSSimTrack.fRcFitPar[4]" },
    { "tr.fit.ndof",  "NDoF of rec hit fit",                "fMCTracks.TGEMSBSSimTrack.fRcFitPar[5]" },
    { "tr.fit.vx",    "Vertex x from rec hit fit [m]",      "fMCTracks.TGEMSBSSimTrack.fRcFitPar[6]" },
    { "tr.fit.vy",    "Vertex y from rec hit fit [m]",      "fMCTracks.TGEMSBSSimTrack.fRcFitPar[7]" },
    { "tr.fit.vz",    "Vertex z from rec hit fit [m]",      "fMCTracks.TGEMSBSSimTrack.fRcFitPar[8]" },

    // "Back tracks": hits of the primary particle in the first tracker plane
    { "btr.n",     "Number of back tracks",     "GetNBackTracks()" },
    { "btr.pid",   "Track PID (PDG)",           "fBackTracks.TGEMSBSSimBackTrack.fPID" },
    { "btr.num",   "GEANT particle number",     "fBackTracks.TGEMSBSSimBackTrack.fType" },
    { "btr.planes","Bitpattern of planes hit",  "fBackTracks.TGEMSBSSimBackTrack.fHitBits" },
    { "btr.ufail", "Undigitized u planes",      "fBackTracks.TGEMSBSSimBackTrack.fUfailBits" },
    { "btr.vfail", "Undigitized v planes",      "fBackTracks.TGEMSBSSimBackTrack.fVfailBits" },
    { "btr.sect",  "Sector number",             "fBackTracks.TGEMSBSSimBackTrack.fSector" },
    { "btr.p",     "Track momentum (GeV)",      "fBackTracks.TGEMSBSSimBackTrack.P() "},
    // Track position in Cartesian/TRANSPORT coordinates, optimal for SBS, not for SoLID
    { "btr.x",     "Track pos lab x [m]",       "fBackTracks.TGEMSBSSimBackTrack.X()" },
    { "btr.y",     "Track pos lab y [m]",       "fBackTracks.TGEMSBSSimBackTrack.Y()" },
    { "btr.th",    "Track dir tan(theta)",      "fBackTracks.TGEMSBSSimBackTrack.ThetaT()" },
    { "btr.ph",    "Track dir tan(phi)",        "fBackTracks.TGEMSBSSimBackTrack.PhiT()" },
    // Track position and direction in cylindrical coordinates, good for SoLID
    // { "btr.r",     "Track pos lab r_trans (m)", "fBackTracks.TGEMSBSSimBackTrack.R()" },
    // { "btr.theta", "Track pos lab theta [rad]", "fBackTracks.TGEMSBSSimBackTrack.Theta()" },
    // { "btr.phi",   "Track pos lab phi [rad]",   "fBackTracks.TGEMSBSSimBackTrack.Phi()" },
    // { "btr.thdir", "Track dir theta [rad]",     "fBackTracks.TGEMSBSSimBackTrack.ThetaDir()" },
    // { "btr.phdir", "Track dir phi [rad]",       "fBackTracks.TGEMSBSSimBackTrack.PhiDir()" },
    // Hit coordinates in first tracker plane, relative to plane origin
    { "btr.hx",    "Track pos plane x [m]",     "fBackTracks.TGEMSBSSimBackTrack.HX()" },
    { "btr.hy",    "Track pos plane y [m]",     "fBackTracks.TGEMSBSSimBackTrack.HY()" },

    // Digitized hits registered in the GEMs
    //    { "hit.n",     "Number of MC hits",          "GetNMCHits()" },
    { "hit.id",    "MC hit number",              "fMCHits.TGEMSBSSimGEMHit.fID" },
    { "hit.sect",  "MC hit sector",              "fMCHits.TGEMSBSSimGEMHit.fSector" },
    { "hit.rsect", "MC hit non-mapped sector",   "fMCHits.TGEMSBSSimGEMHit.fRealSector" },
    { "hit.plane", "MC hit plane",               "fMCHits.TGEMSBSSimGEMHit.fPlane" },
    { "hit.src",   "MC data set source",         "fMCHits.TGEMSBSSimGEMHit.fSource" },
    { "hit.type",  "MC hit GEANT counter",       "fMCHits.TGEMSBSSimGEMHit.fType" },
    { "hit.pid",   "MC hit PID (PDG)",           "fMCHits.TGEMSBSSimGEMHit.fPID" },
    { "hit.p",     "MC hit particle mom [GeV]",  "fMCHits.TGEMSBSSimGEMHit.P()" },
    { "hit.x",     "MC hit lab x position [m]",  "fMCHits.TGEMSBSSimGEMHit.X()" },
    { "hit.y",     "MC hit lab y position [m]",  "fMCHits.TGEMSBSSimGEMHit.Y()" },
    { "hit.z",     "MC hit lab z position [m]",  "fMCHits.TGEMSBSSimGEMHit.Z()" },
    // Hit position in cylindrical/spherical coordinates, good for SoLID
    // { "hit.r",     "MC hit lab r [m]",           "fMCHits.TGEMSBSSimGEMHit.R()" },
    // { "hit.theta", "MC hit lab theta [rad]",     "fMCHits.TGEMSBSSimGEMHit.Theta()" },
    // { "hit.phi",   "MC hit lab phi [rad]",       "fMCHits.TGEMSBSSimGEMHit.Phi()" },
    { "hit.charge","MC hit cluster charge",      "fMCHits.TGEMSBSSimGEMHit.fCharge" },
    { "hit.time",  "MC hit time offset [s]",     "fMCHits.TGEMSBSSimGEMHit.fTime" },
    { "hit.usz",   "MC hit u cluster size",      "fMCHits.TGEMSBSSimGEMHit.fUSize" },
    { "hit.ustart","MC hit u cluster 1st strip", "fMCHits.TGEMSBSSimGEMHit.fUStart" },
    { "hit.upos",  "MC hit u cluster center [m]","fMCHits.TGEMSBSSimGEMHit.fUPos" },
    { "hit.vsz",   "MC hit v cluster size",      "fMCHits.TGEMSBSSimGEMHit.fVSize" },
    { "hit.vstart","MC hit v cluster 1st strip", "fMCHits.TGEMSBSSimGEMHit.fVStart" },
    { "hit.vpos",  "MC hit v cluster center [m]","fMCHits.TGEMSBSSimGEMHit.fVPos" },
    
    { "pt.fmctrk", "MC point track number",      "fMCPoints.Podd::MCTrackPoint.fMCTrack" },
    
    { 0 }
  };

  return THaAnalysisObject::
    DefineVarsFromList( vars, THaAnalysisObject::kRVarDef,
			mode, "", this, MC_PREFIX, here );
}

//-----------------------------------------------------------------------------
void TGEMSBSSimDecoder::Clear( Option_t* opt )
{
  // Clear track and plane data

  SimDecoder::Clear(opt);   // clears fMCHits, fMCTracks and fMCPoints

  fBackTracks->Clear(opt);
  fStripMap.clear();
}

//-----------------------------------------------------------------------------
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
int TGEMSBSSimDecoder::LoadEvent(const UInt_t* evbuffer )
#else
int TGEMSBSSimDecoder::LoadEvent(const Int_t* evbuffer )
#endif
{
  // Wrapper around DoLoadEvent so we can conveniently stop the benchmark
  // counter in case of errors

  int ret = DoLoadEvent( evbuffer );

  if( fDoBench ) fBench->Stop("physics_decode");

  return ret;
}

//-----------------------------------------------------------------------------
static inline
void StripToROC( Int_t s_plane, Int_t s_sector, Int_t s_proj,
		 Int_t s_chan,
		 Int_t& crate, Int_t& slot, Int_t& chan )
{
  // Convert location parameters (plane,sector,proj,chan) of the given strip
  // to hardware channel (crate,slot,chan)
  // The (crate,slot,chan) assignment must match the detmap definition in
  // the database!  See TreeSearch/dbconvert.cxx
  
  // cout << "Chan per slot ? " << fManager->GetChanPerSlot() << endl;
  // cout << "Module per readout ? " << fManager->GetModulesPerReadOut() << endl;
  // cout << "N readout ? " << fManager->GetNReadOut() << ", N Chambers ? " << fManager->GetNChamber() << endl;
  // cout << "Chambers per crate ? " << fManager->GetChambersPerCrate() << endl;
  // cout << "Module per readout ? " << fManager->GetModulesPerChamber() << endl;
  
  div_t d = div( s_chan, fManager->GetChanPerSlot() );
  Int_t module = d.quot;
  chan = d.rem;
  Int_t ix = module +
    fManager->GetModulesPerReadOut()*( s_proj + fManager->GetNReadOut()*( s_plane + fManager->GetNChamber()*s_sector ));
  
  //cout << "StripToROC: module " << module << ", ix " << ix << endl;
  
  d = div( ix, fManager->GetChambersPerCrate()*fManager->GetModulesPerChamber() );
  crate = d.quot;
  slot  = d.rem;
}

//-----------------------------------------------------------------------------
static inline
void StripToROCD( Int_t s_plane, Int_t s_module, Int_t s_proj,
		 Int_t s_chan,
		 Int_t& crate, Int_t& slot, Int_t& chan )
{
  div_t d = div( s_chan, fManager->GetChanPerSlot() );
  //  Int_t module = d.quot;
  chan = d.rem;
  //total slot id
  Int_t ix = s_proj + 2*( s_module + fManager->GetNModule(s_plane-1)*s_plane );
  
  //  cout << "StripToROC: module " << module << ", ix " << ix << Decoder::MAXSLOT<<endl;
  
  d = div( ix, SIM_MAXSLOT);//fManager->GetChambersPerCrate()*fManager->GetModulesPerChamber() );
  crate = d.quot;
  slot  = d.rem;
}

//-----------------------------------------------------------------------------
static inline
Int_t MakeROCKey( Int_t crate, Int_t slot, Int_t chan )
{
  return chan +
    fManager->GetChanPerSlot()*( slot + SIM_MAXSLOT*crate );
}

//-----------------------------------------------------------------------------
Int_t TGEMSBSSimDecoder::StripFromROC( Int_t crate, Int_t slot, Int_t chan ) const
{
  // Return index of digitized strip correspomding to hardware channel
  // (crate,slot,chan)

  if( fStripMap.empty() )
    return -1;

  StripMap_t::const_iterator found = fStripMap.find( MakeROCKey(crate,slot,chan) );
  if( found == fStripMap.end() )
    return -1;

  return found->second;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<Double_t>> TGEMSBSSimDecoder::GetAllMCHits() const
{
  std::vector<std::vector<Double_t>> hits;
  std::vector<Double_t> vtemp = {0,0,0,0,0,0};//v[0]--posx, v[1]--posy, v[2]--charge, v[3]--planeID v[4]--moduleID v[5]time_zero
  assert( buffer );       // Must still have the event buffer
  const TGEMSBSSimEvent* simEvent = reinterpret_cast<const TGEMSBSSimEvent*>(buffer);
  for(size_t i=0;i<simEvent->fGEMClust.size();i++){
    const TGEMSBSSimEvent::GEMCluster& clust = simEvent->fGEMClust[i];
    if(clust.fSource!=0){continue;}
    vtemp[0] = clust.fMCpos.X();
    vtemp[1] = clust.fMCpos.Y();
    vtemp[2] = clust.fCharge;
    vtemp[3] = clust.fPlane;
    vtemp[4] = clust.fModule;
    vtemp[5] = clust.fTime;
    // cout<<"######## "<<vtemp[1]<<endl;
    hits.push_back(vtemp);
  }

  return hits;
}

//-----------------------------------------------------------------------------
TGEMSBSMCHitInfo TGEMSBSSimDecoder::GetSBSMCHitInfo( Int_t crate, Int_t slot, Int_t chan ) const
{
  // Get MC truth info for the given hardware channel
//  const char* const here = "TGEMSBSSimDecoder::GetSBSMCHitInfo";

  Int_t istrip = StripFromROC( crate, slot, chan );
  assert( istrip >= 0 );  // else logic error in caller or bad fStripMap
  
  
  assert( buffer );       // Must still have the event buffer
  const TGEMSBSSimEvent* simEvent = reinterpret_cast<const TGEMSBSSimEvent*>(buffer);
  
  assert( static_cast<vsiz_t>(istrip) < simEvent->fGEMStrips.size() );
  const TGEMSBSSimEvent::DigiGEMStrip& strip = simEvent->fGEMStrips[istrip];
  assert( strip.fProj >= 0 && strip.fProj < fManager->GetNReadOut() );
  
  TGEMSBSMCHitInfo mc;
  mc.fSigType = strip.fSigType;
  // cout<<kPrimaryStrip<<" "<<kSecondaryStrip<<" "<<kInducedStrip<<endl;getchar();0 1 2
    // if(strip.fProj==0 && strip.fPlane==4 && strip.fTime1>50.0)
  //   printf("%f \n", strip.fTime1);
  
  //for cross talk
  if (TESTBIT(strip.fSigType, kInducedStrip) && !TESTBIT(strip.fSigType, kPrimaryStrip) &&
      !TESTBIT(strip.fSigType, kSecondaryStrip) ){
    mc.fMCTrack = 0;
    mc.fMCPos = fManager->GetPosFromModuleStrip(strip.fProj, strip.fPlane, strip.fModule, strip.fChan);
    mc.fMCTime = strip.fTime1;
    
    //cout << "strip = " << strip.fChan << ", time = " << mc.fMCTime << ", pos = " <<  mc.fMCPos << endl;
    return mc;
  }

  mc.fMCCharge = strip.fCharge;
  // cout<<"strip: "<<istrip<<" cc: "<<mc.fMCCharge<<endl;
  Double_t nOverlapSignal = 0.;
  // cout<<strip.fClusters.GetSize()<<" # "<<strip.fClusterRatio[0].GetSize()<<" # "<<strip.fClusterRatio[1].GetSize()<<" # "<<strip.fClusterRatio[2].GetSize()<<" # "<<strip.fClusterRatio[3].GetSize()<<endl;getchar();
  for( Int_t i = 0; i<strip.fClusters.GetSize(); ++i ) {
    Int_t iclust = strip.fClusters[i] - 1;  // yeah, array index = clusterID - 1

   
    //   cout<<mc.vClusterID.size()<<" : "<<mc.vClusterADC[1].size()<<endl;getchar(); 

    assert( iclust >= 0 && static_cast<vsiz_t>(iclust) < simEvent->fGEMClust.size() );
    const TGEMSBSSimEvent::GEMCluster& c = simEvent->fGEMClust[iclust];
    mc.vClusterID.push_back(iclust);

    //work here! remove after finish!
    // add cluster type info here to tell whether a cluster is a primary or background using "fSource" in simevent.h !!
    // then go to get b_overtotal in gemplane.cxx and so u can tell whether a cluster has background in it. and how much! great!
    mc.vClusterType.push_back(c.fSource);




    //
    mc.vClusterPeakTime.push_back(c.fTime);
    mc.vClusterPos.push_back(c.fHitpos);
    mc.vClusterCharge.push_back(c.fCharge);

    mc.vClusterStripWeight.push_back(strip.fStripWeightInCluster[i]);
    // cout<<strip.fStripWeightInCluster[i]<<endl;getchar();
    for(Int_t its=0;its<6;its++)
      {
	mc.vClusterADC[its].push_back(strip.fClusterRatio[its][i]);
      }

    assert( c.fID == iclust+1 );
    assert( strip.fPlane == c.fPlane && strip.fSector == c.fSector );
    Int_t signalID = -1;
    for (unsigned int ii = 0; ii<fSignalInfo.size(); ii++){
      if (c.fType == fSignalInfo.at(ii).tid && c.fPID == fSignalInfo.at(ii).pid) // cluster_type(primary or secondary) == type_requested(primary) && particle == partical_requested
	signalID = ii;
    }
    // cout << "Plane " << strip.fPlane << ", proj (x: 0, y: 1) " << strip.fProj 
    //  	 << ": pos[proj] = "  << c.fXProj[strip.fProj] << endl;
    if( signalID >= 0 && c.fSource == kPrimarySource ) {
      if( mc.fMCTrack > 0 ) {
        //this means that there two signal hits overlapping
        //for now I keep the fMCTrack to the first one, by average the fMCPos nad fMCTime
        //Weizhi Xiong
        //assert(manager->GetNSigParticle() > 1); //otherwise should not happen
        
        mc.fMCPos += c.fXProj[strip.fProj]+(1-strip.fProj)*fManager->GetXOffset(c.fPlane,c.fModule);
        mc.fMCTime += c.fTime; 
      }else{
        // Strip contains a contribution from a primary particle hit :)
        mc.fMCTrack = fSignalInfo.at(signalID).tid; 
        mc.fMCPos   = c.fXProj[strip.fProj]+(1-strip.fProj)*fManager->GetXOffset(c.fPlane,c.fModule);
        mc.fMCTime  = c.fTime;
      }
      nOverlapSignal++;
    } else {
      ++mc.fContam;
      if( mc.fMCTrack == 0 ) {
	mc.fMCPos += c.fXProj[strip.fProj]+(1-strip.fProj)*fManager->GetXOffset(c.fPlane,c.fModule);
      }
    }
  }
  assert( strip.fClusters.GetSize() == 0 || mc.fMCTrack > 0 || mc.fContam > 0 );
  
  if( mc.fMCTrack == 0 ) {
    if( mc.fContam > 1 ) {
      // If only background hits, report the mean position of all those hits
      mc.fMCPos /= static_cast<Double_t>(mc.fContam);
    }
    mc.fMCTime = strip.fTime1;
  }else{
    mc.fMCPos /= nOverlapSignal;
    mc.fMCTime /= nOverlapSignal;
  }
  
  return mc;
}

//-----------------------------------------------------------------------------
static inline Int_t NumberOfSetBits( UInt_t v )
{
  // Count number of bits set in 32-bit integer. From
  // http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel

  v = v - ((v >> 1) & 0x55555555);
  v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
  return (((v + (v >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

//-----------------------------------------------------------------------------
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
Int_t TGEMSBSSimDecoder::DoLoadEvent(const UInt_t* evbuffer )
#else
Int_t TGEMSBSSimDecoder::DoLoadEvent(const Int_t* evbuffer )
#endif
{
  // cout<<"Tsbs sim decoder DoLoadEvent"<<endl;
  // Fill crateslot structures with Monte Carlo event data in 'evbuffer'
  static const char* const here = "TGEMSBSSimDecoder::LoadEvent";

#if ANALYZER_VERSION_CODE < ANALYZER_VERSION(1,6,0)
  Bool_t fNeedInit = fgNeedInit;
#endif
  assert( fMap || fNeedInit );

  // Local copy of evbuffer pointer, used in GetMCHitInfo
  buffer = evbuffer;

  // Cast the evbuffer pointer back to exactly the event type that is present
  // in the input file (in TGEMSBSSimFile). The pointer-to-unsigned integer is
  // needed compatibility with the standard decoder.
  const TGEMSBSSimEvent* simEvent = reinterpret_cast<const TGEMSBSSimEvent*>(buffer);
 
  /* 
  cout<<simEvent->fGEMClust.size()<<endl;getchar();
  for(Int_t i=0;i<simEvent->fGEMClust.size();i++){
    cout<<"Plane: "<<simEvent->fGEMClust[i].fPlane<<" type: "<<simEvent->fGEMClust[i].fType<<"   x: "<<simEvent->fGEMClust[i].fMCpos.X()
	<<"   y: "<<simEvent->fGEMClust[i].fMCpos.Y()<<endl;
  }
  */

  Int_t ret = HED_OK;
  if (first_decode || fNeedInit) {
    if( (ret = init_cmap()) != HED_OK )
      return ret;
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
    if( (ret = init_slotdata()) != HED_OK)
#else
    if( (ret = init_slotdata(fMap)) != HED_OK)
#endif
      return ret;
    first_decode = false;
  }

  if( fDoBench ) fBench->Begin("clearEvent");
  Clear();
  for( int i=0; i<fNSlotClear; i++ )
    crateslot[fSlotClear[i]]->clearEvent();
  if( fDoBench ) fBench->Stop("clearEvent");

  // FIXME: needed?
  evscaler = 0;
  event_length = 0;
  
  event_type = 1;
  event_num = simEvent->fEvtID;
  recent_event = event_num;

  // Event weight
  fWeight = simEvent->fWeight;

  //
  if( fDoBench ) fBench->Begin("physics_decode");

  // Decode the digitized strip data.  Populate crateslot array.
  for( vector<TGEMSBSSimEvent::DigiGEMStrip>::size_type i = 0;
       i < simEvent->fGEMStrips.size(); i++) {
    //cout << "i " << i << endl;
    const TGEMSBSSimEvent::DigiGEMStrip& s = simEvent->fGEMStrips[i];
    Int_t crate, slot, chan;
    //cout << "striptoroc: " << endl;
    //StripToROC( s.fPlane, s.fSector, s.fProj, s.fChan, crate, slot, chan );
    StripToROCD( s.fPlane, s.fModule, s.fProj, s.fChan, crate, slot, chan );
    //cout<<"Plane: "<<s.fPlane<<" Module:  "<<s.fModule<<" projection:  "<<s.fProj<<" channel:  "<<s.fChan<<" sigType: "<<s.fSigType<<endl;
    //cout << "crate = " << crate << ", slot = " << slot << ", chan " << chan << endl;getchar();
    //cout << "samples: " << endl;
    for( Int_t k = 0; k < s.fNsamp; k++ ) { 
      Int_t raw = s.fADC[k];
      //cout << raw << " ### ";
      
      if( crateslot[idx(crate,slot)]->loadData("adc",chan,raw,raw) == SD_ERR )
	return HED_ERR;
    }
    //    cout << endl<<endl;
    //cout << "stripmap : " << endl;
    // Build map from ROC address to strip index. This is needed to extract
    // the MC truth info later in the tracking detector decoder via GetMCChanInfo.
#ifndef NDEBUG
    pair<StripMap_t::const_iterator,bool> ins =
#endif
      fStripMap.insert( make_pair( MakeROCKey(crate,slot,chan), i ) );
    //cout<<crate<<" "<<slot<<" "<<chan<<endl;
    // cout<<MakeROCKey(crate,slot,chan)<<endl;
    // getchar();
    
    // cout << "ROC key inserted in strip map " << endl;
    //  cout << "ins.second ? " << ins.second << endl;
    assert( ins.second );
  }
  
  // Create lists of two types of tracks:
  // 1) Physics tracks, as generated at the target
  // 2) "Back tracks": hits in any GEM plane from the primary particle

  // Physics tracks. We need to copy them here so we can export them as global
  // variables.
  TClonesArray* tracks = simEvent->fMCTracks;
  assert( tracks );
  for( Int_t i = 0; i < tracks->GetLast()+1; i++ ) {
    TGEMSBSSimTrack* trk = static_cast<TGEMSBSSimTrack*>(tracks->UncheckedAt(i));
   new( (*fMCTracks)[i] ) TGEMSBSSimTrack(*trk);
  }
  assert( GetNMCTracks() > 0 );

  // MC hit data ("clusters") and "back tracks"
  Int_t best_primary = -1, best_primary_plane = fManager->GetNChamber(), primary_sector = -1;
  UInt_t primary_hitbits = 0, ufail = 0, vfail = 0;
  for( vector<TGEMSBSSimEvent::GEMCluster>::size_type i = 0;
       i < simEvent->fGEMClust.size(); ++i ) {
    const TGEMSBSSimEvent::GEMCluster& c = simEvent->fGEMClust[i];

    if( c.fPlane < 0 || c.fPlane >= fManager->GetNChamber() ) {
      Error( here, "Illegal plane number = %d in cluster. "
	     "Should never happen. Call expert.", c.fPlane );
      simEvent->Print("clust");
      return HED_FATAL;
    }

    // Save hits in the GEMs
    new( (*fMCHits)[GetNMCHits()] ) TGEMSBSSimGEMHit(c);

    // Extra bookkeeping for primary tracks, used for making back tracks below
    if( c.fType == kPrimaryType && c.fSource == kPrimarySource ) {
      // Record the primary track's points for access via the SimDecoder interface.
      // Record one point per projection so that we can study residuals.
      Int_t itrack = 1;
      primary_sector = c.fSector;
      MCTrackPoint* upt = // kUPlane changed to kXPlane: necessary to match TreeSearch EProjType
	new( (*fMCPoints)[GetNMCPoints()] ) MCTrackPoint( itrack,
							  c.fPlane, kXPlane,
							  c.fMCpos, c.fP );
      upt->fMCTime = c.fTime;
      MCTrackPoint* vpt =// kVPlane changed to kYPlane: necessary to match TreeSearch EProjType
	new( (*fMCPoints)[GetNMCPoints()] ) MCTrackPoint( itrack,
							  c.fPlane, kYPlane,
							  c.fMCpos, c.fP );
      vpt->fMCTime = c.fTime;
      
      // //debug...
      // cout << "TGEMSBSSimDecoder.cxx: Print MC points " << endl;
      // cout << "kXplane ? " << kXPlane << endl;
      // upt->Print("");
      // cout << "kVYlane ? " << kYPlane << endl;
      // vpt->Print("");
      
      // Keep bitpattern of planes crossed by this primary
      SETBIT(primary_hitbits,c.fPlane);

      //cout << "Plane number " << c.fPlane << ", primary hitbits " << primary_hitbits << endl; 
      
      // Save index of the primary particle hit closest to plane 0
      if( c.fPlane < best_primary_plane ) {
	best_primary = i;
	best_primary_plane = c.fPlane;
      }
      // Determine digitization hit inefficiency: Check if this MC hit
      // activated GEM strips in both readout planes
      if( c.fSize[0] == 0 ) {
	SETBIT(ufail, c.fPlane);
	CLRBIT(upt->fStatus, MCTrackPoint::kDigitized);
      } else {
	SETBIT(upt->fStatus, MCTrackPoint::kDigitized);
      }
      if( c.fSize[1] == 0 ) {
	SETBIT(vfail, c.fPlane);
	CLRBIT(vpt->fStatus, MCTrackPoint::kDigitized);
      } else {
	SETBIT(vpt->fStatus, MCTrackPoint::kDigitized);
      }
    }
  }

  // Sort fMCPoints by type (u,v) and plane number, then calculate plane-to-plane
  // differences. The following assumes that all points are from the same track
  // (ensured above). If that is no longer so one day, fMCPoints will need to
  // be sorted by track number as well, and the algo below needs to be changed.
  fMCPoints->Sort();
  Double_t mass = 0;
  TGEMSBSSimTrack* trk = static_cast<TGEMSBSSimTrack*>(fMCTracks->UncheckedAt(0));
  assert(trk);
  if( TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(trk->fPID) )
    mass = particle->Mass();
  else
    Warning( "LoadEvent", "No enrty in PDG database for PID = %d", trk->fPID );

  MCTrackPoint* prev_pt = 0;
  for( Int_t i = 0; i < GetNMCPoints(); ++i ) {
    MCTrackPoint* pt = static_cast<MCTrackPoint*>( fMCPoints->UncheckedAt(i) );
    assert(pt);
    if( prev_pt && prev_pt->fType == pt->fType ) {
      assert( pt->fMCTrack == prev_pt->fMCTrack );
      if( prev_pt->fPlane+1 == pt->fPlane ) {
	pt->fDeltaE = TMath::Sqrt(prev_pt->fMCP.Mag2() + mass*mass) -
	  TMath::Sqrt(pt->fMCP.Mag2() + mass*mass);
	pt->fDeflect = prev_pt->fMCP.Angle(pt->fMCP);
	pt->fToF = pt->fMCTime - prev_pt->fMCTime;
      }
    }
    prev_pt = pt;
  }

  // Keep statistics in the MC track
  trk->fNHits = 2*NumberOfSetBits(primary_hitbits);
  trk->fHitBits = primary_hitbits;

  // "Back tracks"
  // Record the apparent track from the primary particle
  // of the signal data here, i.e. type == 1 and source == 0.
  // There is only ever one primary particle per event.
  if( best_primary >= 0 ) {

    Int_t nback = GetNBackTracks();
    assert( nback == 0 );

    TGEMSBSSimBackTrack* btr = new( (*fBackTracks)[nback] )
      TGEMSBSSimBackTrack(simEvent->fGEMClust[best_primary]);

    //cout << "Backtrack primary hitbits " << primary_hitbits << endl;
    
    btr->SetHitBits(primary_hitbits);
    btr->SetUfailBits(ufail);
    btr->SetVfailBits(vfail);

    // Use the back track to emulate calorimeter hits.
    // Assumptions:
    // - Only tracks crossing all fManager->GetNChamber() GEMs (points in all planes)
    //   make a calorimeter hit. This is a crude model for the trigger.
    // - The track propagates without deflection from the last GEM plane
    //   to the front of the emulated calorimeter.
    // - The measured calorimeter position is independent of the incident
    //   track angle.
    if( fManager->DoCalo() && trk->fNHits == 2*fManager->GetNChamber() ) {
      // Retrieve last MC track point
      assert( GetNMCPoints() == 2*fManager->GetNChamber() );
      MCTrackPoint* pt =
	static_cast<MCTrackPoint*>( fMCPoints->UncheckedAt(2*fManager->GetNChamber()-1) );
      assert( pt );
      const TVector3& pos = pt->fMCPoint;
      TVector3 dir = pt->fMCP.Unit();
      if( fManager->GetCaloZ() <= pos.Z() ) {
	Error( here, "Calorimeter z = %lf less than z of last GEM plane = %lf. "
	       "Set correct value with SetCaloZ() or turn off calo emulation.",
	       fManager->GetCaloZ(), pos.Z() );
	return HED_FATAL;
      }
      if( TMath::Abs(dir.Z()) < 1e-6 ) {
	Error( here, "Illegal primary track direction (%lf,%lf,%lf). "
	       "Should never happen. Call expert.", dir.X(), dir.Y(), dir.Z() );
	return HED_ERR;
      }
      dir *= 1.0/dir.Z();  // Make dir a transport vector
      TVector3 hitpos = pos + (fManager->GetCaloZ()-pos.Z()) * dir;

      // Smear the position with the given resolution
      // Assumes z-axis normal to calorimeter plane. Otherwise we would have to
      // get the plane's fXax and fYax
      TVector3 res( gRandom->Gaus(0.0, fManager->GetCaloRes()),
		    gRandom->Gaus(0.0, fManager->GetCaloRes()), 0.0 );
      hitpos += res;

      // Encode the raw hit data for the dummy GEM planes.
      // The actual coordinate transformation to u or v takes place in each
      // plane's Decode() where all the required geometry information is at hand.
      // This bypasses any type of digitization. That should be fine for dummy
      // planes where we want to inject known lab hit positions into the tracking.
      //
      // Because of the way the detector map is layed out at the moment,
      // we place the calorimeter in fake sector 31, so the data are in two slots
      // (for u and v coordinates, respectively) in the ROC immediately
      // following the GEM trackers for sector 30. In each slot, channels
      // 0-29 correspond to the sector of the MC track sector (should always be 0
      // if mapping sectors. Each "hit" corresponds to one measured position.
      // Currently, there is only ever one hit per channel since there is only
      // one MC track. The hit's raw data are hitpos.X(), the data, hitpos.Y(),
      // each a Float_t value interpreted as Int_t.
      assert( primary_sector == 0 );

      union FloatIntUnion {
	Float_t f;
	Int_t   i;
      } datx, daty;
      datx.f = static_cast<Float_t>(hitpos.X());
      daty.f = static_cast<Float_t>(hitpos.Y());

      Int_t crate, slot, chan;
      //StripToROC( 0, fManager->GetNSector(), kUPlane, primary_sector, crate, slot, chan );
      StripToROC( 0, fManager->GetNSector(), kXPlane, primary_sector, crate, slot, chan );
       if( crateslot[idx(crate,slot)]->loadData("adc",chan,datx.i,daty.i)
	  == SD_ERR )
	 {
	   return HED_ERR;
	 }
	//StripToROC( 0, fManager->GetNSector(), kVPlane, primary_sector, crate, slot, chan );
      StripToROC( 0, fManager->GetNSector(), kYPlane, primary_sector, crate, slot, chan );
      if( crateslot[idx(crate,slot)]->loadData("adc",chan,datx.i,daty.i)
	  == SD_ERR )
	{
	  return HED_ERR;
	}
    }
  }

  // DEBUG:
  //cout << "SimDecoder: nTracks = " << GetNMCTracks() << endl;
  //fMCTracks.Print();
  // cout<<" \n\n###\n\n"<<endl;
  return HED_OK;
}

//-----------------------------------------------------------------------------
TGEMSBSSimGEMHit::TGEMSBSSimGEMHit( const TGEMSBSSimEvent::GEMCluster& c )
  : fID(c.fID), fSector(c.fSector), fPlane(c.fPlane),
    fRealSector(c.fRealSector), fSource(c.fSource), fType(c.fType),
    fPID(c.fPID), fP(c.fP), fXEntry(c.fXEntry), fMCpos(c.fMCpos),
    fHitpos(c.fHitpos), fCharge(c.fCharge), fTime(c.fTime),
    fUSize(c.fSize[0]), fUStart(c.fStart[0]), fUPos(c.fXProj[0]),
    fVSize(c.fSize[1]), fVStart(c.fStart[1]), fVPos(c.fXProj[1])
{
  // Construct hit from cluster
}

//-----------------------------------------------------------------------------
void TGEMSBSSimGEMHit::Print( const Option_t* ) const
{
  // Print TGEMSBSSimGEMHit info

}

//-----------------------------------------------------------------------------
TGEMSBSSimBackTrack::TGEMSBSSimBackTrack( const TGEMSBSSimEvent::GEMCluster& c )
  : fType(c.fType), fPID(c.fPID), fSector(c.fSector), fSource(c.fSource),
    fHitBits(0), fUfailBits(0), fVfailBits(0)
{
  // Construct track from cluster info

  Update( c );
  SetHitBit( c.fPlane );
}

//-----------------------------------------------------------------------------
Int_t TGEMSBSSimBackTrack::Update( const TGEMSBSSimEvent::GEMCluster& c )
{
  // Project track coordinates to first tracker plane

  static const char* const here = "TGEMSBSSimBackTrack::Update";

  // Currently not needed since Update only called from constructor
  // if( fType != c.fType || fPID != c.fPID || fSector != c.fSector ) {
  //   Error( here, "Updating with inconsistent GEMCluster data: "
  // 	   "type = %d/%d, pid = %d/%d, sector = %d/%d.\n"
  // 	   "Should never happen. Call expert.",
  // 	   fType, c.fType, fPID, c.fPID, fSector, c.fSector );
  //   return -1;
  // }

  if( c.fPlane > 0 ) {
    Double_t dz = c.fMCpos.Z() - fManager->GetZ0();
    if( dz <= 0 ) {
      Error( here, "Illegal fMCpos z-coordinate in plane = %d. "
	     "Should never happen. Call expert.", c.fPlane );
      c.fMCpos.Print();
      return -2;
    }
    fOrigin = c.fMCpos - dz * c.fP.Unit();
  } else {
    fOrigin = c.fMCpos;
  }
  fHitpos = c.fHitpos; // FIXME: project this, too?
  fMomentum = c.fP;

  return 0;
}

//-----------------------------------------------------------------------------
void TGEMSBSSimBackTrack::Print( const Option_t* ) const
{
  // Print TGEMSBSSimBackTrack info

  cout << "track: type = " << fType
       << ", PID = "       << fPID
       << ", sector = "    << fSector
       << endl;
  cout << "  Origin    = ";  fOrigin.Print();
  cout << "  Momentum  = ";  fMomentum.Print();
  cout << "  Hitpos    = ";  fHitpos.Print();
  cout << "  HitBits   = " << fHitBits << endl;
}

//-----------------------------------------------------------------------------
ClassImp(TGEMSBSSimDecoder)
ClassImp(TGEMSBSSimGEMHit)
ClassImp(TGEMSBSSimBackTrack)