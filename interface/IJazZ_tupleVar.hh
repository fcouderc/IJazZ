//---------------------------------------------------------------------------//
// Description: This class is a wrapper so any Z tree can be used as input
//              The public variables should exist and be properly filled
//                          via the loadEvent() function
//                          any selection that are ztree input dependent should
//                          be implemented in the passPrivateSelection();
//
//
//
// Author: Fabrice Couderc, fabrice.couderc@cea.fr 
// History: 2012/11/27 - v1
//
//---------------------------------------------------------------------------//
#ifndef ijazz_tupleVar_hh__
#define ijazz_tupleVar_hh__

#include "interface/EcalUtils.hh"
#include "interface/IJazZAxisND.hh"

#include <math.h>
#include <vector>

#define nEleInTuple 3

class TTree;
class TChain;
class TFile;
class TH3F;
class TH1F;

class IJazZdiEle {
public:
  IJazZdiEle(TChain *c);
  IJazZdiEle(TTree  *c);
  ~IJazZdiEle();

  void isMC   (bool ismc = true   ) { _isMC = ismc; }
  void set7TeV(bool is7TeV = true ) { _7TeV = is7TeV; }

  //--- absolutely needed by IJazZ, they have to be filled properly
  unsigned int  runId;
  unsigned char nPV;
  unsigned long long evtId;
  //  unsigned long evtId;
  float esOverEcal[nEleInTuple];
  float eSC[nEleInTuple], pt[nEleInTuple], R9[nEleInTuple];
  float std_eSC[nEleInTuple];
  float scEta[nEleInTuple], scPhi[nEleInTuple];
  float eta[nEleInTuple], phi[nEleInTuple];
  float seedEnF[nEleInTuple];
  int evtQuality;

  //--- will be filled by IJazZ
  float weight;
  
  //--- potentially needed depending on the config options
  float eSC_err[nEleInTuple];
  float pho_eSC[nEleInTuple], pho_eSC_err[nEleInTuple];
  
  //--- required for MC oversmearing checks (config options)
  float energyMCEle[nEleInTuple],etaMCEle[nEleInTuple],phiMCEle[nEleInTuple];

  //--- not called by IJazZ but used either for selection or 
  //    to fill the required variables
  unsigned int ID[nEleInTuple];
  unsigned char tracker[nEleInTuple];
  short q[nEleInTuple];
  float avgLC[nEleInTuple];
  float raw_eSC[nEleInTuple], raw_eES[nEleInTuple], raw_eSeed[nEleInTuple];
  float raw_eES1[nEleInTuple], raw_eES2[nEleInTuple]; 
  float iEta[nEleInTuple], iPhi[nEleInTuple], iDee[nEleInTuple], iHarn[nEleInTuple];
  float iX[nEleInTuple], iY[nEleInTuple];
  short seedIx[nEleInTuple], seedIy[nEleInTuple];
  float deltaR_MC[nEleInTuple];
  float corr[nEleInTuple];
  float _mee_MC, mZ_MC;
  float rho;




  //------ needed for E/p ntuples format
  int isZ;
  float raw_eSC_ele1,raw_eES_ele1,raw_e3x3_ele1,std_eSC_ele1,raw_eSeed_ele1,eta_ele1,phi_ele1,scEta_ele1,scPhi_ele1,avgLC_ele1,q_ele1;
  int seedIx_ele1,seedIy_ele1,seedIEta_ele1,seedIPhi_ele1;
  float raw_eSC_ele2,raw_eES_ele2,raw_e3x3_ele2,std_eSC_ele2,raw_eSeed_ele2,eta_ele2,phi_ele2,scEta_ele2,scPhi_ele2,avgLC_ele2,q_ele2;
  int seedIx_ele2,seedIy_ele2,seedIEta_ele2,seedIPhi_ele2;




  double mee() { return sqrt(2*pt[0]*pt[1]*(cosh(eta[0]-eta[1])-cos(phi[0]-phi[1]))); }
  void corrResp( double resp[nEleInTuple] ) { 
    for( int iele = 0; iele < 2; iele++ ) { 
      corr[iele] = 1./resp[iele];
      eSC[iele] /= resp[iele]; pt[iele] /= resp[iele]; 
      std_eSC[iele] /= resp[iele];
    } 
  }

  //-- IJazZ input tree
  void setBranchAddress(bool addMCtruth);
  int  getEntries(void); 
  int  loadEvent(int ievt);
  void setEnergyCorrectionType(int type) {_scEnergyCorrType=type;}

  //-- IJazZ does not make much selections they should be done here
  bool passPrivateSelection();
  
  //-- IJazZ output tree (user defined)
  void createOutputTree( TFile * fileout );
  void fillOutputTree(void);
  void writeTree(void);

  //-- axis binder: find axis bin for IJazZ
  void fill_Xvar_forBinFinder( IJazZAxisND<double> *eCorr, std::vector<double> *x_eCorr );

  int _ncat;
private: 
  bool _7TeV;
  bool _isMC;
  TTree *_ztree;
  TTree *_ztreeOut;
  TH3F  *_hMCTruthStudy;
  TH1F  *_hMeeTruthStudyMul[20];
  TH1F  *_hMeeTruthStudyAdd[20];
  TH1F  *_hMeeRecoStudy[20];
  TH1F  *_hMeeTruth[20];
  

  float _mee_raw;
  float _mee_corr;
  float _fZ;
  int  _scEnergyCorrType;
  harnessEcal _harnessDef;
  
};



#endif
