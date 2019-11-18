#include "TSBSSimCher.h"
#include <iostream>
#include <TSBSSimData.h>
#include <TSBSSimEvent.h>
#include "TSBSDBManager.h"
//#include "sbs_types.h"

TSBSSimCher::TSBSSimCher(const char* name, short id)
{
  fName = name;
  SetUniqueDetID(id);
  Init();
}

TSBSSimCher::~TSBSSimCher()
{
}

void TSBSSimCher::Init()
{
  TSBSSimDetector::Init();
  if(fDebug>=1)
    cout << "Cherenkov detector with UniqueDetID = " << UniqueDetID() << ": TSBSSimCher::Init() " << endl;
  
  // Get the Detector info
  fDetInfo = fDBmanager->GetDetInfo(fName.Data());
  
  // Get all necessary info to parameterize the PMT pulse shape.
  double tau = fDetInfo.DigInfo().SPE_Tau();
  double sigma = fDetInfo.DigInfo().SPE_Sigma();
  double tmin = -fDetInfo.DigInfo().GateWidth()/2.0;
  double tmax = +fDetInfo.DigInfo().GateWidth()/2.0;
  double t0 = 0.0;//+fDetInfo.DigInfo().SPE_TransitTime()-fDetInfo.DigInfo().TriggerOffset();
  
  fSPE = new TSPEModel(fName.Data(), tau, sigma, t0, tmin, tmax);
  
  //Configure the PMT signals array
  fSignals.resize(fDetInfo.NChan());
  for(size_t i_ch = 0; i_ch<fDetInfo.NChan(); i_ch++){
    fSignals[i_ch].SetNpeChargeConv(fDetInfo.DigInfo().NpeChargeConv(i_ch));
  }
}


void TSBSSimCher::LoadEventData(const std::vector<g4sbshitdata*> &evbuffer)
{
  Clear();
  LoadAccumulateData(evbuffer);
  /*
  bool signal = false;
  int chan = 0;
  int type = 0;
  double time = 0;
  double data = 0;
  
  // for( const g4sbshitdata *ev: evbuffer) {
  for(std::vector<g4sbshitdata*>::const_iterator it = evbuffer.begin(); it!= evbuffer.end(); ++it ) {
    g4sbshitdata* ev = (*it);
    // Only get detector data for Cherenkov
    // new detector ID convetion proposal: UniqueDetID = 10*DetType+DetID
    if(ev->GetDetUniqueID() == UniqueDetID()) {
      signal = (ev->GetData(0)==0);
      chan = ev->GetData(1);
      type = ev->GetData(2);
      time = ev->GetData(3)+fDetInfo.DigInfo().SPE_TransitTime()-fDetInfo.DigInfo().TriggerOffset()+fDetInfo.DigInfo().TriggerJitter();//add 
      data = ev->GetData(4);
      
      if(fDebug>=3)
	cout << "Detector " << UniqueDetID() << " chan = " << chan << " Npe " << data << " Evt Time: "<< ev->GetData(3) << " " << time << endl;
      
      if(type == 0) {
        //std::cout << "Filling data for chan: " << ev->GetData(0) << ", t=" << 
        // ev->GetData(1) - 60. << std::endl;
        //if(ev->GetData(1)<mint)
        //  mint = ev->GetData(1);
	//fSPE->SetNpe(data);
        //fSignals[chan].Fill(chan, fNPE, data);//
	fSignals[chan].Fill(fSPE, data, fDetInfo.DigInfo().Threshold(chan), time, signal);//
      }
    }
  }
  */
}

void TSBSSimCher::LoadAccumulateData(const std::vector<g4sbshitdata*> &evbuffer)
{
  bool signal = false;
  int chan = 0;
  int type = 0;
  double time = 0;
  double data = 0;
  
  // for( const g4sbshitdata *ev: evbuffer) {
  for(std::vector<g4sbshitdata*>::const_iterator it = evbuffer.begin(); it!= evbuffer.end(); ++it ) {
    g4sbshitdata* ev = (*it);
    // Only get detector data for Cherenkov
    // new detector ID convetion proposal: UniqueDetID = 10*DetType+DetID
    if(ev->GetDetUniqueID() == UniqueDetID()) {
      signal = (ev->GetData(0)==0);
      chan = ev->GetData(1);
      type = ev->GetData(2);
      time = ev->GetData(3)+fDetInfo.DigInfo().SPE_TransitTime()-fDetInfo.DigInfo().TriggerOffset() + fTimeZero;//+fDetInfo.DigInfo().TriggerJitter()
      data = ev->GetData(4);
      
      //if(fabs(time)>fDetInfo.DigInfo().GateWidth()/2.0)continue;
      
      if(fDebug>=3)
	cout << "Detector " << UniqueDetID() << " chan = " << chan << " Npe " << data << " Evt Time: "<< ev->GetData(3) << " " << time << endl;
      
      if(type == 0) {
        //std::cout << "Filling data for chan: " << ev->GetData(0) << ", t=" << 
        // ev->GetData(1) - 60. << std::endl;
        //if(ev->GetData(1)<mint)
        //  mint = ev->GetData(1);
	//fSPE->SetNpe(data);
        //fSignals[chan].Fill(chan, fNPE, data);//
	fSignals[chan].Fill(fSPE, data, fDetInfo.DigInfo().Threshold(chan), time, signal);//
      }
    }
  }
}

void TSBSSimCher::Digitize(TSBSSimEvent &event)
{
  bool any_events = false;

  if(fDebug>=3)cout << "TSBSSimCher::Digitize() : Unique Det ID " << UniqueDetID() << " signal size = " << fSignals.size() << endl;
  
  //TSBSSimEvent::DetectorData data;
  //TSBSSimEvent::SimDetectorData simdata;
  
  //UInt_t VETROCword;
  
  //bool header[8] = {0, 0, 0, 0, 0, 0, 1, 1};
  //bool channel[8];
  //bool tdc[16];
  //bool trail;
  //short edgebitpos = 26;
  
  SimEncoder::adc_data adc_data;
  std::vector<uint32_t> data;
  UInt_t tdcval;
  //std::vector<double> simdata;
  int mult = 0;
  for(size_t m = 0; m < fSignals.size(); m++) {
    //data.fData.clear();
    //simdata.fData.clear();
    data.clear();
    //simdata.clear();
    fSignals[m].Digitize(fDetInfo.DigInfo(), m);
    if(fSignals[m].Npe() > 0) {
      any_events = true;
      //data.fDetID = UniqueDetID();
      //data.fChannel = m;
      
      if(fDebug>=4)cout << " = > fSignals[" << m << "].TDCSize() " << fSignals[m].TDCSize() << endl;
      
      /*
      //event.SimDetID.push_back(Short_t(UniqueDetID()));
      event.SimDetChannel[fDetInfo.DetFullName()].push_back(Short_t(m));
      event.SimDetDataType[fDetInfo.DetFullName()].push_back(1);
      event.SimDetNData[fDetInfo.DetFullName()].push_back(1);
      simdata.push_back(fSignals[m].Npe());
      event.SimDetData[fDetInfo.DetFullName()].push_back(simdata);
      event.NSimDetData[fDetInfo.DetFullName()]++;
      simdata.clear();
      
      //event.SimDetID.push_back(Short_t(UniqueDetID()));
      event.SimDetChannel[fDetInfo.DetFullName()].push_back(Short_t(m));
      event.SimDetDataType[fDetInfo.DetFullName()].push_back(2);
      event.SimDetNData[fDetInfo.DetFullName()].push_back(Short_t(fSignals[m].LeadTimesSize()+fSignals[m].TrailTimesSize()));
      for(size_t i = 0; i<fSignals[m].LeadTimesSize(); i++){
	simdata.push_back(fSignals[m].LeadTime(i));
	
	if(fDebug>=3)cout << " leadtime " << i << " = " << fSignals[m].LeadTime(i) << endl;
      }
      for(size_t i = 0; i<fSignals[m].TrailTimesSize(); i++){
	simdata.push_back(fSignals[m].TrailTime(i));
	if(fDebug>=3)cout << " trail time " << i << " = " << fSignals[m].TrailTime(i) << endl;;
      }
      event.SimDetData[fDetInfo.DetFullName()].push_back(simdata);
      event.NSimDetData[fDetInfo.DetFullName()]++;
      simdata.clear();
      */
      
      if(fDebug>=4){
	cout << fDetInfo.DetFullName() << ": check MC vec size" << endl;
	fSignals[m].check_vec_size();
      }
      for(uint i_mc = 0; i_mc<fSignals[m].MCHitSize(); i_mc++){
	event.NSimDetHits[fDetInfo.DetFullName()]++;
	event.SimDetChannel[fDetInfo.DetFullName()].push_back(Short_t(m));
	event.SimDetNpe[fDetInfo.DetFullName()].push_back(fSignals[m].MCHitNpe(i_mc));
	event.SimDetTime[fDetInfo.DetFullName()].push_back(fSignals[m].MCHitTime(i_mc));
	if(fEncoderTDC){
	  if(fDebug>=4)cout << fSignals[m].MCHitLeadTime(i_mc) << " " << fSignals[m].MCHitTrailTime(i_mc) << endl;
	  event.SimDetLeadTime[fDetInfo.DetFullName()].push_back(fSignals[m].MCHitLeadTime(i_mc));
	  event.SimDetTrailTime[fDetInfo.DetFullName()].push_back(fSignals[m].MCHitTrailTime(i_mc));
	}
      }
       
      
      //define convention for type:
      // 0: ADC
      // 1: TDC
      mult = 0;

      // Fill ADC
      if(fEncoderADC) {
        adc_data.integral=fSignals[m].ADC();
        fEncoderADC->EncodeADC(adc_data,fEncBuffer,fNEncBufferWords);
        CopyEncodedData(fEncoderADC,mult++,data);//.fData);
	
	for(uint i = 0; i<data.size(); i++){
	  event.NDetHits[fDetInfo.DetFullName()]++;
	  event.DetChannel[fDetInfo.DetFullName()].push_back(Short_t(m));
	  event.DetDataWord[fDetInfo.DetFullName()].push_back(data.at(i));
	  if(i==0){//header
	    event.DetADC[fDetInfo.DetFullName()].push_back(-1000000);
	    if(fEncoderTDC){
	      event.DetTDC_L[fDetInfo.DetFullName()].push_back(-1000000);
	      event.DetTDC_T[fDetInfo.DetFullName()].push_back(-1000000);
	    }
	  }else{
	    event.DetADC[fDetInfo.DetFullName()].push_back(fSignals[m].ADC()-fDetInfo.DigInfo().Pedestal(m));
	    if(fEncoderTDC){
	      event.DetTDC_L[fDetInfo.DetFullName()].push_back(-1000000);
	      event.DetTDC_T[fDetInfo.DetFullName()].push_back(-1000000);
	    }
	  }
	}
	data.clear();
      }
      
      // Fill TDC 
      if(fEncoderTDC) {
 	if(fDebug>=4)cout << fSignals[m].TDC(0) << " " << fSignals[m].TDC(1) << endl;
	fEncoderTDC->EncodeTDC(fSignals[m].TDCData(),fEncBuffer,
            fNEncBufferWords);
        CopyEncodedData(fEncoderTDC,mult++,data);//.fData);
	
	for(uint i = 0; i<data.size(); i++){
	  event.NDetHits[fDetInfo.DetFullName()]++;
	  event.DetChannel[fDetInfo.DetFullName()].push_back(Short_t(m));
	  event.DetDataWord[fDetInfo.DetFullName()].push_back(data.at(i));
	  if(i==0){//header
	    if(fEncoderADC)event.DetADC[fDetInfo.DetFullName()].push_back(-1000000);
	    event.DetTDC_L[fDetInfo.DetFullName()].push_back(-1000000);
	    event.DetTDC_T[fDetInfo.DetFullName()].push_back(-1000000);
	  }else{
	    if(fEncoderADC)event.DetADC[fDetInfo.DetFullName()].push_back(-1000000);
	    if( fSignals[m].TDC(i-1) & ( 1 << (31) ) ){
	      tdcval = fSignals[m].TDC(i-1);
	      if(fDebug>=4)cout << " T: " << tdcval << " => ";
	      tdcval ^= ( -0 ^ tdcval) & ( 1 << (31) );
	      if(fDebug>=4)cout << tdcval << endl;
	      event.DetTDC_L[fDetInfo.DetFullName()].push_back(-1000000);
	      event.DetTDC_T[fDetInfo.DetFullName()].push_back(tdcval);
	    }else{
	      event.DetTDC_L[fDetInfo.DetFullName()].push_back(fSignals[m].TDC(i-1));
	      if(fDebug>=4)cout << " L: " << fSignals[m].TDC(i-1) << endl;
	      event.DetTDC_T[fDetInfo.DetFullName()].push_back(-1000000);
	    }
	  }
	}
	data.clear();
      }
      
            /*
      //event.DetID.push_back(Short_t(UniqueDetID()));
      event.DetChannel[fDetInfo.DetFullName()].push_back(Short_t(m));
      event.DetNData[fDetInfo.DetFullName()].push_back(Short_t(data.size()));
      event.DetData[fDetInfo.DetFullName()].push_back(data);
      event.NDetData[fDetInfo.DetFullName()]++;

      //event.fDetectorData.push_back(data);
      //data.fData.clear();
      data.clear();
      */
    }
  }
  SetHasDataFlag(any_events);
}

// Clear signals in array
void TSBSSimCher::Clear(Option_t *)
{
  for(size_t i = 0; i < fSignals.size(); i++ ) {
    fSignals[i].Clear();
  }
}
ClassImp(TSBSSimCher) // Implements TSBSSimCher
