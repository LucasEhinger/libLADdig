#include "SBSDigGEMDet.h"
#include "TMath.h"

using namespace std;

SBSDigGEMDet::SBSDigGEMDet()
{
}

SBSDigGEMDet::SBSDigGEMDet(UShort_t uniqueid, UInt_t nplanes, int* layer, int* nstrips, double* offset, double* roangle, int nsamp, double zsup_thr):
  fUniqueID(uniqueid), fNPlanes(nplanes)
{
  //for(uint i = 0; i<fNPlanes; i++)GEMPlanes[i] = SBSDigGEMPlane(nstrips[i], nsamp, zsup_thr);
  for(uint i = 0; i<fNPlanes; i++){
    cout << i << " " << layer[i] << " " << nstrips[i] << " " << offset[i] << " " << roangle[i] << endl; 
    GEMPlanes.push_back(SBSDigGEMPlane(layer[i], i/2, nstrips[i], nsamp, zsup_thr, offset[i], roangle[i]));
  }
}

SBSDigGEMDet::~SBSDigGEMDet()
{
  
}

void SBSDigGEMDet::Clear()
{
  for(uint i = 0; i<fNPlanes; i++)GEMPlanes[i].Clear();
  fGEMhits.clear();
}
