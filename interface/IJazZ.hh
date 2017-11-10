//---------------------------------------------------------------------------//
// Description: IJazZ Main class: handles the option and the event looper
//
//
// Author: Fabrice Couderc, fabrice.couderc@cea.fr 
// History: 2012/11/28 - v1
//
//---------------------------------------------------------------------------//

#ifndef IJazZ_hh_
#define IJazZ_hh_

#include "interface/IJazZAxis.hh"
#include "interface/IJazZAxisND.hh"

#include <string>
#include <vector>

namespace ROOT { namespace Math { class Minimizer;}; };

class TFile;
class TChain;
class TTree;
class TH1F;
class ZFitterMinuit2;
class TGraphErrors;

TChain* addInputFiles( std::string filelist, std::string treename );

class IJazZ {
public:
  IJazZ(void);
  ~IJazZ(void);

  /// Fitter ran over MC / Data
  void isMC( bool ismc = true ) { _isMC = int(ismc); }

  /// i/o
  void doFit( bool fit ) { _doFIT = fit;}
  void createOutputTree(bool create) { _outputTree = create;}

  /// create output file name
  void setVersion( std::string vers ) {_version = vers; }
  void setUserDir( std::string uDir ) {_userDir = uDir; }  
  std::string createOutputFileName( std::string ecalRegion );

  /// setup the files containing the tree
  void addInputFilesMC(   std::string filelist, std::string treeName = "selected" ) { _ztreeMC   = addInputFiles( filelist, treeName ); }
  void addInputFilesData( std::string filelist, std::string treeName = "selected" ) { _ztreeData = addInputFiles( filelist, treeName ); }

  /// define run ranges
  void defineRunRange(int rmin, int rmax) { _runMin = rmin; _runMax = rmax; }

  /// correct for already fitter corrections
  void addResponseCorrection( std::vector<std::string> respCorrFiles ) { _file_eRespCorr = respCorrFiles; }
  void addResoDataFile(       std::vector<std::string> resoDataFile  ) { _resoDataFile   = resoDataFile;  }
  void addResoMCFile(         std::vector<std::string> resoMCFile    ) { _resoMCFile     = resoMCFile;    }
  
  /// EB or EE fit
  void isEE(bool isee = true ) { _isEEfit = isee; }

  /// energy correction method
  void scEnergyCorrection( int type) { _scEnergyCorrType = type; }

  /// define eta bining (return true is eta axis is positive only, ie EE fit instead of EE- and EE+)
  bool setBiningND( std::string fileAxisResp = "undefined",  std::string fileAxisReso = "undefined");

  void useEtaScaleBining(void) { _analysis =  1; }
  void useExternalBining(void) { _analysis =  0; }
  void useTestBining(void)     { _analysis = -1; }

  /// IJazZ classic mode (removing low mass tail)
  void setIIazZClassicMode( bool ijazzClassic = true ) { _ijazzClassic = ijazzClassic; }
  
  /// run on a subsample only
  void setAllEvenOddEvents(int allEvenOddEvents ) { _allEvenOddEvents = allEvenOddEvents; }

  /// cuts
  void setPtCut( float ptcut ) { _ptCut = ptcut; }

  /// prepare vtx reweighting histograms
  void setupNvtxReweighting(void);
  float getNvtxWeight( int nVtx );

  /// setup fitter & minimizer
  void mcOversmearingStudy(int mcStudy ) { _testOversmearingMC = mcStudy; }
  void setupZFitter( int ncpu);
  void eventSelectionZFitter(std::string ecalRegion);
  void minimize(std::string ecalp);
  void saveFitResults(ROOT::Math::Minimizer*, std::string* parName);
  
  /// debugging printout (default level = 0)
  void setPrintLevel( int level = 1 ) { _debugLevel = level; }


  /// analysis output result
  static  void etaHistoStyle( TGraphErrors *g, const IJazZAxis<double> *axis, int color = 1, std::string yTitle = "",
			      int yMode = -1, float yMin = -1, float yMax = -1 );
  void etaScaleFromZeeAnaFit(std::string mcVersion="notdefined");

  TGraphErrors* graphMerger( int yBin, bool ER, bool isIeta = false ); /// FIX me: can certainly remove that, obsolete

private:
  /// axis
  IJazZAxisND<double>  _respAxisND, _resoAxisND;
  std::vector<IJazZAxisND<double> >  _tailAxisND;

  TH1F *_hNvtxMC;
  TH1F *_hNvtxData;
  TH1F *_hNvtxWeight;
  
  TChain *_ztreeData;
  TChain *_ztreeMC  ;
  std::string _ztreeName;
  
  std::string _version;
  std::string _userDir;
  std::string _ecalp;
  std::string _EEringsDefinitionFile[2];
  
  bool _isEEfit;
  int  _isMC;
  unsigned int  _runMin,_runMax;
  int  _analysis;
  bool _bin3D;
  bool _ijazzClassic;
  int  _allEvenOddEvents;
  //--- mc oversmearing test
  int _testOversmearingMC;

  //--- addtional cuts
  float _ptCut;

  //--- actual ZFitter
  int _nParameters;
  ZFitterMinuit2 *_zfitterFab;
  float _meeFitMin,_meeFitMax;

  //--- energy correction type
  int  _scEnergyCorrType;
  bool _fitRelResoForRegrEnCorr;

  //--- outputs
  bool  _outputTree;
  bool  _doFIT;
  TFile *_fileout;

  //--- debugging
  int _debugLevel;

  //--- already fitted corrections
  std::vector<std::string>  _file_eRespCorr;
  std::vector<std::string> _resoDataFile, _resoMCFile;


  std::vector<double> _tailPar0Rotate_EB,_tailPar0Rotate_EE;
};


//--- FIX Me should be somewhere else, this is just for test
void checkTestMCSmearing(void);

#endif
