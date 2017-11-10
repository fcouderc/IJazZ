#ifndef __energyCorrectionAndSmearingHgg__hh__
#define __energyCorrectionAndSmearingHgg__hh__
#include <string>
#include <vector>

class EnergyScaleReader {
 public:
  EnergyScaleReader( void );
  ~EnergyScaleReader(void);

  int EcalPart( float r9, float scEta );
  int EcalPart( std::string ecal );

 
  std::string EcalPartString( float r9, float scEta );

  int runRange( int run, int iEcal );

  bool setup( std::string file );


  //-------------- new, legacy
  int EcalPart_EtDep( float r9, float scEta, float Et  );
  int EcalPart_EtDep( std::string ecal ); 
  std::string EcalPartString_EtDep( float r9, float scEta, float Et );
  bool setup_EtDep( std::string file );
  //--------------


  //----- common energy scale for EtDep and not EtDep correction
  float energyScale(   float r9, float scEta, float Et, int run, int systShift = 0 );

  private:
  bool _isEtDep;
  bool _isSetup;
  std::vector<int>    _runStart[24];
  std::vector<int>    _runStop[24];
  std::vector<std::string> _ecalPart;
  std::vector<float>  _eScale[24];
  std::vector<float>  _err_eScale[24];

  int _nWarnings;

  /// private so nobody can access it... use common function 
  float energyScale_noEtDep( float r9, float scEta, int run, int systShift = 0 );
  float energyScale_EtDep(   float r9, float scEta, float Et, int run, int systShift = 0 );


};


#include <TRandom3.h> 
class photonOverSmearing {
 public:
  enum category { EBlowEtaGold_CM = 0, EBlowEtaGold_GAP = 1, 
		  EBlowEtaBad_CM  = 2, EBlowEtaBad_GAP  = 3,
		  EBhighEtaGold   = 4, EBhighEtaBad     = 5,
		  EElowEtaGold    = 6, EElowEtaBad      = 7,
		  EEhighEtaGold   = 8, EEhighEtaBad     = 9,
		  UNKNOWN = 10 };

  photonOverSmearing( ) {};
  photonOverSmearing::category photonOverSmearingCategory( float scEta, float r9, bool isEBGap );
  void initialize( std::string setupType );
  void setSmearSeed(int evtNum, int runNum, int lumisec, float SCPhi, float SCEta);
  float meanOverSmearing( float scEta, float r9, bool isEBGap, int syst = 0 );
  float meanOverSmearing( float scEta, float r9, float et,bool isEBGap, int syst = 0 );
  float randOverSmearing( float scEta, float r9, float et, bool isEBGap, int syst = 0 );

 private:
  std::vector<float> _oversmearing;
  std::vector<float> _oversmearing_err;
  std::vector<float> _oversmearingstoch_rho;
  std::vector<float> _oversmearingstoch_rhoerr;
  std::vector<float> _oversmearingstoch_EMean;
  std::vector<float> _oversmearingstoch_EMeanerr;
  std::vector<float> _oversmearingstoch_phi;
  std::vector<float> _oversmearingstoch_phierr;
  TRandom3 overSmearingGen;
  std::string setupType_;
};


#endif
