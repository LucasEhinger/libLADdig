#ifndef _G4SBS_DATA_H
#define _G4SBS_DATA_H

#include <Rtypes.h>
class TTree;
class TBranch;

namespace TSBSGeant4 {

// Purely virtual data structure for detector data in a G4SBS ROOT tree.
struct VDetData_t {
  VDetData_t() {};
  virtual ~VDetData_t() {};
  // All sub-classes *must* implement a concrete instance of this
  virtual bool SetupBranches(TTree *t, const char *prefix) = 0;
  template <typename T> int SetupBranch(TTree *tree, const char *prefix, const char *varname, T &var);
};

struct CalData_t : public VDetData_t {
  int nhits;
  std::vector<int> *row;
  std::vector<int> *col;
  std::vector<int> *plane;
  std::vector<int> *paddle;

  std::vector<double> *xhit;
  std::vector<double> *yhit;
  std::vector<double> *zhit;

  std::vector<double> *sumedep;
  std::vector<double> *tavg;
  std::vector<double> *trms;
  std::vector<double> *tmin;
  std::vector<double> *tmax;
  CalData_t()
      : nhits(0), row(0), col(0), plane(0), paddle(0), xhit(0), yhit(0), zhit(0), sumedep(0), tavg(0), trms(0), tmin(0),
        tmax(0) {}
  virtual ~CalData_t() {};
  virtual bool SetupBranches(TTree *t, const char *prefix);
};



struct GEMData_t : public VDetData_t {
  Int_t nhits;
  std::vector<int> *plane;
  std::vector<int> *strip;
  std::vector<double> *x;
  std::vector<double> *y;
  std::vector<double> *z;
  std::vector<double> *polx;
  std::vector<double> *poly;
  std::vector<double> *polz;
  std::vector<double> *t;
  std::vector<double> *trms;
  std::vector<double> *tmin;
  std::vector<double> *tmax;
  std::vector<double> *tx;
  std::vector<double> *ty;
  std::vector<double> *txp;
  std::vector<double> *typ;
  std::vector<double> *xg;
  std::vector<double> *yg;
  std::vector<double> *zg;
  std::vector<int> *trid;
  std::vector<int> *mid;
  std::vector<int> *pid;
  std::vector<double> *vx;
  std::vector<double> *vy;
  std::vector<double> *vz;
  std::vector<double> *p;
  std::vector<double> *edep;
  std::vector<double> *beta;
  std::vector<double> *xin;
  std::vector<double> *yin;
  std::vector<double> *zin;
  std::vector<double> *xout;
  std::vector<double> *yout;
  std::vector<double> *zout;

  GEMData_t()
      : nhits(0), plane(0), strip(0), x(0), y(0), z(0), polx(0), poly(0), polz(0), t(0), trms(0), tmin(0), tmax(0),
        tx(0), ty(0), txp(0), typ(0), xg(0), yg(0), zg(0), trid(0), mid(0), pid(0), vx(0), vy(0), vz(0), p(0), edep(0),
        beta(0), xin(0), yin(0), zin(0), xout(0), yout(0), zout(0) {};
  virtual ~GEMData_t() {};
  virtual bool SetupBranches(TTree *t, const char *prefix);
};

struct DigCalData_t : public VDetData_t {
  Int_t nchan;
  std::vector<Int_t> *chan;
  std::vector<Int_t> *adc;

  TBranch *b_nchan; //!
  TBranch *b_chan;  //!
  TBranch *b_adc;   //!

  DigCalData_t() : nchan(0), chan(0), adc(0), b_nchan(0), b_chan(0), b_adc(0) {};
  virtual ~DigCalData_t() {};
  virtual bool SetupBranches(TTree *t, const char *prefix);
  void ClearBranches();
  void FillBranches();
};

struct DigCalFADC7Data_t : public VDetData_t {
  Int_t nchan;
  std::vector<Int_t> *chan;
  std::vector<Int_t> *adc;
  std::vector<Int_t> *tdc;
  std::vector<Int_t> *amp;
  std::vector<Int_t> *ped;

  TBranch *b_nchan; //!
  TBranch *b_chan;  //!
  TBranch *b_adc;   //!
  TBranch *b_tdc;   //!
  TBranch *b_amp;   //!
  TBranch *b_ped;   //!

  DigCalFADC7Data_t()
      : nchan(0), chan(0), adc(0), tdc(0), amp(0), ped(0), b_nchan(0), b_chan(0), b_adc(0), b_tdc(0), b_amp(0),
        b_ped(0) {};
  virtual ~DigCalFADC7Data_t() {};
  virtual bool SetupBranches(TTree *t, const char *prefix);
  void ClearBranches();
  void FillBranches();
};

struct DigTimingData_t : public VDetData_t {
  Int_t nchan;
  std::vector<Int_t> *chan;
  std::vector<Int_t> *adc;
  std::vector<Int_t> *tdc_l;
  std::vector<Int_t> *tdc_t;
  std::vector<Int_t> *amp;

  TBranch *b_nchan; //!
  TBranch *b_chan;  //!
  TBranch *b_adc;   //!
  TBranch *b_tdc_l; //!
  TBranch *b_tdc_t; //!
  TBranch *b_amp;   //!

  DigTimingData_t()
      : nchan(0), chan(0), adc(0), tdc_l(0), tdc_t(0), amp(0), b_nchan(0), b_chan(0), b_adc(0), b_tdc_l(0), b_tdc_t(0), b_amp(0) {};
  virtual ~DigTimingData_t() {};
  virtual bool SetupBranches(TTree *t, const char *prefix);
  void ClearBranches();
  void FillBranches();
};

struct DigSampCalData_t : public VDetData_t {
  Int_t nchan;
  std::vector<Int_t> *chan;
  /*
  std::vector<Int_t> *adc_0;
  std::vector<Int_t> *adc_1;
  std::vector<Int_t> *adc_2;
  std::vector<Int_t> *adc_3;
  std::vector<Int_t> *adc_4;
  std::vector<Int_t> *adc_5;
  std::vector<Int_t> *adc_6;
  std::vector<Int_t> *adc_7;
  std::vector<Int_t> *adc_8;
  std::vector<Int_t> *adc_9;
  std::vector<Int_t> *adc_10;
  std::vector<Int_t> *adc_11;
  std::vector<Int_t> *adc_12;
  std::vector<Int_t> *adc_13;
  std::vector<Int_t> *adc_14;
  std::vector<Int_t> *adc_15;
  std::vector<Int_t> *adc_16;
  std::vector<Int_t> *adc_17;
  std::vector<Int_t> *adc_18;
  std::vector<Int_t> *adc_19;
  */
  std::vector<Int_t> *adc;
  std::vector<Int_t> *samp;
  std::vector<Int_t> *tdc;

  TBranch *b_nchan; //!
  TBranch *b_chan;  //!
  /*
  TBranch *b_adc_0;   //!
  TBranch *b_adc_1;   //!
  TBranch *b_adc_2;   //!
  TBranch *b_adc_3;   //!
  TBranch *b_adc_4;   //!
  TBranch *b_adc_5;   //!
  TBranch *b_adc_6;   //!
  TBranch *b_adc_7;   //!
  TBranch *b_adc_8;   //!
  TBranch *b_adc_9;   //!
  TBranch *b_adc_10;   //!
  TBranch *b_adc_11;   //!
  TBranch *b_adc_12;   //!
  TBranch *b_adc_13;   //!
  TBranch *b_adc_14;   //!
  TBranch *b_adc_15;   //!
  TBranch *b_adc_16;   //!
  TBranch *b_adc_17;   //!
  TBranch *b_adc_18;   //!
  TBranch *b_adc_19;   //!
  */
  TBranch *b_adc;  //!
  TBranch *b_samp; //!
  TBranch *b_tdc;  //!

  DigSampCalData_t()
      : nchan(0), chan(0), adc(0), samp(0),
        tdc(0), // adc_0(0), adc_1(0), adc_2(0), adc_3(0), adc_4(0), adc_5(0), adc_6(0), adc_7(0), adc_8(0), adc_9(0),
                // adc_10(0), adc_11(0), adc_12(0), adc_13(0), adc_14(0), adc_15(0), adc_16(0), adc_17(0), adc_18(0),
                // adc_19(0), tdc(0), b_nchan(0), b_chan(0), b_adc_0(0), b_adc_1(0), b_adc_2(0), b_adc_3(0), b_adc_4(0),
                // b_adc_5(0), b_adc_6(0), b_adc_7(0), b_adc_8(0), b_adc_9(0), b_adc_10(0), b_adc_11(0), b_adc_12(0),
                // b_adc_13(0), b_adc_14(0), b_adc_15(0), b_adc_16(0), b_adc_17(0), b_adc_18(0), b_adc_19(0),
        b_nchan(0), b_chan(0), b_adc(0), b_samp(0), b_tdc(0) {};
  virtual ~DigSampCalData_t() {};
  virtual bool SetupBranches(TTree *t, const char *prefix);
  void ClearBranches();
  void FillBranches();
};

struct DigGEMData_t : public VDetData_t {
  Int_t nstrips;
  // std::vector<Int_t> *planeID;
  std::vector<Int_t> *module;
  std::vector<Int_t> *strip;
  /*
  std::vector<Int_t> *adc_0;
  std::vector<Int_t> *adc_1;
  std::vector<Int_t> *adc_2;
  std::vector<Int_t> *adc_3;
  std::vector<Int_t> *adc_4;
  std::vector<Int_t> *adc_5;
  */
  std::vector<Int_t> *adc;
  std::vector<Int_t> *samp;

  TBranch *b_nstrips; //!
  TBranch *b_module;  //!
  TBranch *b_strip;   //!
  /*
  TBranch *b_adc_0;   //!
  TBranch *b_adc_1;   //!
  TBranch *b_adc_2;   //!
  TBranch *b_adc_3;   //!
  TBranch *b_adc_4;   //!
  TBranch *b_adc_5;   //!
  */
  TBranch *b_adc;  //!
  TBranch *b_samp; //!

  DigGEMData_t()
      : nstrips(0), module(0),
        strip(0), // adc_0(0), adc_1(0), adc_2(0), adc_3(0), adc_4(0), adc_5(0), b_nstrips(0), b_module(0), b_strip(0),
                  // b_adc_0(0), b_adc_1(0), b_adc_2(0), b_adc_3(0), b_adc_4(0), b_adc_5(0)
        adc(0), samp(0), b_nstrips(0), b_module(0), b_strip(0), b_adc(0), b_samp(0) {};

  virtual ~DigGEMData_t() {};
  virtual bool SetupBranches(TTree *t, const char *prefix);
  void ClearBranches();
  void FillBranches();
};

} // namespace TSBSGeant4

#endif // _G4SBS_DATA_H
