#include "g4sbs_data.h"
#include <TTree.h>
#include <TBranch.h>
#include <iostream>


namespace TSBSGeant4 {
  template<typename T>
  int VDetData_t::SetupBranch(TTree *tree, const char* prefix,
      const char* varname, T &var)
  {
    TString branchname = TString::Format("%s.%s",prefix,varname);
    if(!tree)
      return 1;
    var = 0;
    int ret = tree->SetBranchAddress(branchname.Data(),&var);
    if( ret != 0 ) {
      std::cerr << "Unable to set branch '" << branchname
        << "' failed with error code: " << ret << std::endl;
      return 1;
    }

    return 0;
  }

  //Hodo
  bool CalData_t::SetupBranches(TTree *tree, const char* prefix)
  {
    int ret = 0;
    ret += SetupBranch(tree,prefix,"nhits", nhits);
    // ret += SetupBranch(tree,prefix,"row", row);
    // ret += SetupBranch(tree,prefix,"col", col);
    ret += SetupBranch(tree,prefix,"plane", plane);
    ret += SetupBranch(tree,prefix,"paddle", paddle);
    // ret += SetupBranch(tree,prefix,"xcell", xcell);
    // ret += SetupBranch(tree,prefix,"ycell", ycell);
    // ret += SetupBranch(tree,prefix,"zcell", zcell);
    // ret += SetupBranch(tree,prefix,"xcellg", xcellg);
    // ret += SetupBranch(tree,prefix,"ycellg", ycellg);
    // ret += SetupBranch(tree,prefix,"zcellg", zcellg);
    ret += SetupBranch(tree,prefix,"xhit", xhit);
    ret += SetupBranch(tree,prefix,"yhit", yhit);
    ret += SetupBranch(tree,prefix,"zhit", zhit);
    // ret += SetupBranch(tree,prefix,"xhitg", xhitg);
    // ret += SetupBranch(tree,prefix,"yhitg", yhitg);
    // ret += SetupBranch(tree,prefix,"zhitg", zhitg);
    ret += SetupBranch(tree,prefix,"sumedep", sumedep);
    ret += SetupBranch(tree,prefix,"tavg", tavg);
    // ret += SetupBranch(tree,prefix,"trms", trms);
    // ret += SetupBranch(tree,prefix,"tmin", tmin);
    // ret += SetupBranch(tree,prefix,"tmax", tmax);
    return (ret ==0);
  }
  

  bool GEMData_t::SetupBranches(TTree* tree, const char *prefix)
  {
    int ret = 0;
    ret += SetupBranch(tree,prefix,"nhits",nhits);
    ret += SetupBranch(tree,prefix,"plane",plane);
    ret += SetupBranch(tree,prefix,"strip",strip);
    ret += SetupBranch(tree,prefix,"x",x);
    ret += SetupBranch(tree,prefix,"y",y);
    ret += SetupBranch(tree,prefix,"z",z);
    ret += SetupBranch(tree,prefix,"polx",polx);
    ret += SetupBranch(tree,prefix,"poly",poly);
    ret += SetupBranch(tree,prefix,"polz",polz);
    ret += SetupBranch(tree,prefix,"t",t);
    ret += SetupBranch(tree,prefix,"trms",trms);
    ret += SetupBranch(tree,prefix,"tmin",tmin);
    ret += SetupBranch(tree,prefix,"tmax",tmax);
    ret += SetupBranch(tree,prefix,"tx",tx);
    ret += SetupBranch(tree,prefix,"ty",ty);
    ret += SetupBranch(tree,prefix,"txp",txp);
    ret += SetupBranch(tree,prefix,"typ",typ);
    ret += SetupBranch(tree,prefix,"xg",xg);
    ret += SetupBranch(tree,prefix,"yg",yg);
    ret += SetupBranch(tree,prefix,"zg",zg);
    ret += SetupBranch(tree,prefix,"trid",trid);
    ret += SetupBranch(tree,prefix,"mid",mid);
    ret += SetupBranch(tree,prefix,"pid",pid);
    ret += SetupBranch(tree,prefix,"vx",vx);
    ret += SetupBranch(tree,prefix,"vy",vy);
    ret += SetupBranch(tree,prefix,"vz",vz);
    ret += SetupBranch(tree,prefix,"p",p);
    ret += SetupBranch(tree,prefix,"edep",edep);
    ret += SetupBranch(tree,prefix,"beta",beta);
    ret += SetupBranch(tree,prefix,"xin",xin);
    ret += SetupBranch(tree,prefix,"yin",yin);
    ret += SetupBranch(tree,prefix,"zin",zin);
    ret += SetupBranch(tree,prefix,"xout",xout);
    ret += SetupBranch(tree,prefix,"yout",yout);
    ret += SetupBranch(tree,prefix,"zout",zout);
    return (ret==0);
  }
  
  
  bool DigCalData_t::SetupBranches(TTree* tree, const char *prefix)
  {
    if(!tree)return(false);
    chan = new std::vector<int>;
    adc = new std::vector<int>;
    b_nchan = tree->Branch(Form("%s.nchan", prefix), &nchan);
    b_chan = tree->Branch(Form("%s.chan", prefix), &chan);
    b_adc = tree->Branch(Form("%s.adc", prefix), &adc);
    // b_tdc = tree->Branch(Form("%s.tdc", prefix), &tdc);
    // b_amp = tree->Branch(Form("%s.amp", prefix), &amp);
    // b_ped = tree->Branch(Form("%s.ped", prefix), &ped);
    return true;
  }
  
  void DigCalData_t::ClearBranches()
  {
    if(chan){//if one var is defined they all are
      nchan = 0;
      chan->clear();
      adc->clear();
      // tdc->clear();
      // amp->clear();
      // ped->clear();
    }
  }
  
  void DigCalData_t::FillBranches()
  {
    if(b_nchan){//if one branch is defined they all are
      b_nchan->Fill();
      b_chan->Fill();
      b_adc->Fill();
      // b_tdc->Fill();
      // b_amp->Fill();
      // b_ped->Fill();
    }
  }
  
  bool DigCalFADC7Data_t::SetupBranches(TTree* tree, const char *prefix)
  {
    if(!tree)return(false);
    chan = new std::vector<int>;
    adc = new std::vector<int>;
    tdc = new std::vector<int>;
    amp = new std::vector<int>;
    ped = new std::vector<int>;
    b_nchan = tree->Branch(Form("%s.nchan", prefix), &nchan);
    b_chan = tree->Branch(Form("%s.chan", prefix), &chan);
    b_adc = tree->Branch(Form("%s.adc", prefix), &adc);
    b_tdc = tree->Branch(Form("%s.tdc", prefix), &tdc);
    b_amp = tree->Branch(Form("%s.amp", prefix), &amp);
    b_ped = tree->Branch(Form("%s.ped", prefix), &ped);
    return true;
  }
  
  void DigCalFADC7Data_t::ClearBranches()
  {
    if(chan){//if one var is defined they all are
      nchan = 0;
      chan->clear();
      adc->clear();
      tdc->clear();
      amp->clear();
      ped->clear();
    }
  }
  
  void DigCalFADC7Data_t::FillBranches()
  {
    if(b_nchan){//if one branch is defined they all are
      b_nchan->Fill();
      b_chan->Fill();
      b_adc->Fill();
      b_tdc->Fill();
      b_amp->Fill();
      b_ped->Fill();
    }
  }
  
  bool DigTimingData_t::SetupBranches(TTree* tree, const char *prefix)
  {
    if(!tree)return(false);
    chan = new std::vector<int>;
    adc = new std::vector<int>;
    tdc_l = new std::vector<int>;
    tdc_t = new std::vector<int>;
    
    b_nchan = tree->Branch(Form("%s.nchan", prefix), &nchan);
    b_chan = tree->Branch(Form("%s.chan", prefix), &chan);
    b_adc = tree->Branch(Form("%s.adc", prefix), &adc);
    b_tdc_l = tree->Branch(Form("%s.tdc_l", prefix), &tdc_l);
    b_tdc_t = tree->Branch(Form("%s.tdc_t", prefix), &tdc_t);
    return true;
  }
  
  void DigTimingData_t::ClearBranches()
  {
    if(chan){
      nchan = 0;
      chan->clear();
      adc->clear();
      tdc_l->clear();
      tdc_t->clear();
    }
  }
  
  void DigTimingData_t::FillBranches()
  {
    if(b_nchan){
      b_nchan->Fill();
      b_chan->Fill();
      b_adc->Fill();
      b_tdc_l->Fill();
      b_tdc_t->Fill();
    }
  }
  
  bool DigSampCalData_t::SetupBranches(TTree* tree, const char *prefix)
  {
    if(!tree)return(false);
    chan = new std::vector<int>;
    /*
    adc_0 = new std::vector<int>;
    adc_1 = new std::vector<int>;
    adc_2 = new std::vector<int>;
    adc_3 = new std::vector<int>;
    adc_4 = new std::vector<int>;
    adc_5 = new std::vector<int>;
    adc_6 = new std::vector<int>;
    adc_7 = new std::vector<int>;
    adc_8 = new std::vector<int>;
    adc_9 = new std::vector<int>;
    adc_10 = new std::vector<int>;
    adc_11 = new std::vector<int>;
    adc_12 = new std::vector<int>;
    adc_13 = new std::vector<int>;
    adc_14 = new std::vector<int>;
    adc_15 = new std::vector<int>;
    adc_16 = new std::vector<int>;
    adc_17 = new std::vector<int>;
    adc_18 = new std::vector<int>;
    adc_19 = new std::vector<int>;
    */
    adc = new std::vector<int>;
    samp = new std::vector<int>;
    tdc = new std::vector<int>;
    
    b_nchan = tree->Branch(Form("%s.nchan", prefix), &nchan);
    b_chan = tree->Branch(Form("%s.chan", prefix), &chan);
    /*
    b_adc_0 = tree->Branch(Form("%s.adc_0", prefix), &adc_0);
    b_adc_1 = tree->Branch(Form("%s.adc_1", prefix), &adc_1);
    b_adc_2 = tree->Branch(Form("%s.adc_2", prefix), &adc_2);
    b_adc_3 = tree->Branch(Form("%s.adc_3", prefix), &adc_3);
    b_adc_4 = tree->Branch(Form("%s.adc_4", prefix), &adc_4);
    b_adc_5 = tree->Branch(Form("%s.adc_5", prefix), &adc_5);
    b_adc_6 = tree->Branch(Form("%s.adc_6", prefix), &adc_6);
    b_adc_7 = tree->Branch(Form("%s.adc_7", prefix), &adc_7);
    b_adc_8 = tree->Branch(Form("%s.adc_8", prefix), &adc_8);
    b_adc_9 = tree->Branch(Form("%s.adc_9", prefix), &adc_9);
    b_adc_10 = tree->Branch(Form("%s.adc_10", prefix), &adc_10);
    b_adc_11 = tree->Branch(Form("%s.adc_11", prefix), &adc_11);
    b_adc_12 = tree->Branch(Form("%s.adc_12", prefix), &adc_12);
    b_adc_13 = tree->Branch(Form("%s.adc_13", prefix), &adc_13);
    b_adc_14 = tree->Branch(Form("%s.adc_14", prefix), &adc_14);
    b_adc_15 = tree->Branch(Form("%s.adc_15", prefix), &adc_15);
    b_adc_16 = tree->Branch(Form("%s.adc_16", prefix), &adc_16);
    b_adc_17 = tree->Branch(Form("%s.adc_17", prefix), &adc_17);
    b_adc_18 = tree->Branch(Form("%s.adc_18", prefix), &adc_18);
    b_adc_19 = tree->Branch(Form("%s.adc_19", prefix), &adc_19);
    */
    b_adc = tree->Branch(Form("%s.adc", prefix), &adc);
    b_samp = tree->Branch(Form("%s.samp", prefix), &samp);
    b_tdc = tree->Branch(Form("%s.tdc", prefix), &tdc);
    return true;
  }

  void DigSampCalData_t::ClearBranches()
  {
    if(chan){
      nchan = 0;
      chan->clear();
      /*
      adc_0->clear();
      adc_1->clear();
      adc_2->clear();
      adc_3->clear();
      adc_4->clear();
      adc_5->clear();
      adc_6->clear();
      adc_7->clear();
      adc_8->clear();
      adc_9->clear();
      adc_10->clear();
      adc_11->clear();
      adc_12->clear();
      adc_13->clear();
      adc_14->clear();
      adc_15->clear();
      adc_16->clear();
      adc_17->clear();
      adc_18->clear();
      adc_19->clear();
      */
      adc->clear();
      samp->clear();
      tdc->clear();
    }
  }
  
  void DigSampCalData_t::FillBranches()
  {
    if(b_nchan){
      b_nchan->Fill();
      b_chan->Fill();
      /*
      b_adc_0->Fill();
      b_adc_1->Fill();
      b_adc_2->Fill();
      b_adc_3->Fill();
      b_adc_4->Fill();
      b_adc_5->Fill();
      b_adc_6->Fill();
      b_adc_7->Fill();
      b_adc_8->Fill();
      b_adc_9->Fill();
      b_adc_10->Fill();
      b_adc_11->Fill();
      b_adc_12->Fill();
      b_adc_13->Fill();
      b_adc_14->Fill();
      b_adc_15->Fill();
      b_adc_16->Fill();
      b_adc_17->Fill();
      b_adc_18->Fill();
      b_adc_19->Fill();
      */
      b_adc->Fill();
      b_samp->Fill();
      b_tdc->Fill();
    }
  }
  
  bool DigGEMData_t::SetupBranches(TTree* tree, const char *prefix)
  {
    if(!tree)return(false);
    
    module = new std::vector<int>;
    strip = new std::vector<int>;
    adc = new std::vector<int>;
    samp = new std::vector<int>;
    
    b_nstrips = tree->Branch(Form("%s.nstrips", prefix), &nstrips);
    b_module = tree->Branch(Form("%s.module", prefix), &module);
    b_strip = tree->Branch(Form("%s.strip", prefix), &strip);
    /*
    b_adc_0 = tree->Branch(Form("%s.adc_0", prefix), &adc_0);
    b_adc_1 = tree->Branch(Form("%s.adc_1", prefix), &adc_1);
    b_adc_2 = tree->Branch(Form("%s.adc_2", prefix), &adc_2);
    b_adc_3 = tree->Branch(Form("%s.adc_3", prefix), &adc_3);
    b_adc_4 = tree->Branch(Form("%s.adc_4", prefix), &adc_4);
    b_adc_5 = tree->Branch(Form("%s.adc_5", prefix), &adc_5);
    */
    b_adc = tree->Branch(Form("%s.adc", prefix), &adc);
    b_samp = tree->Branch(Form("%s.samp", prefix), &samp);
    return true;
  }
  
  void DigGEMData_t::ClearBranches()
  {
    if(strip){
      nstrips = 0;
      strip->clear();
      module->clear();
      /*
      adc_0->clear();
      adc_1->clear();
      adc_2->clear();
      adc_3->clear();
      adc_4->clear();
      adc_5->clear();
      */
      adc->clear();
      samp->clear();
    }
  }
  
  void DigGEMData_t::FillBranches()
  {
    if(b_nstrips){
      b_nstrips->Fill();
      b_module->Fill();
      b_strip->Fill();
      /*
      b_adc_0->Fill();
      b_adc_1->Fill();
      b_adc_2->Fill();
      b_adc_3->Fill();
      b_adc_4->Fill();
      b_adc_5->Fill();
      */
      b_adc->Fill();
      b_samp->Fill();
    }
  } 

}

