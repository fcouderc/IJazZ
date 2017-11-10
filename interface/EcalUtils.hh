#ifndef ecalUtils__hh__
#define ecalUtils__hh__

#include <vector>
#include <string>

#include "interface/IJazZAxisND.hh"
#include "interface/IJazZAxisViewer.hh"

#define NETARING_EB 85
#define NETARING_EE 39

class TLine;
class TH1F;
class TH2F;
class TGraphErrors;

//----- run range checker: from an IJazZ fit output file, extract the run range from the directory path
//                         this assumes std output path building from IJazZ: a.k.a ../IJazZ/IJazZ_output/FitResults/<UserDir>/run1-run2/xxxx.Resp
class runRangeChecker {
public:
  runRangeChecker(std::string iJazZfile);
  runRangeChecker(const runRangeChecker &r) { _runMin = r._runMin; _runMax = r._runMax; }
  bool inRange( int run ) { return (run >= _runMin && run <= _runMax); }
  
  int rMin() {return _runMin;}
  int rMax() {return _runMax;}
private:
  int _runMin;
  int _runMax;
};


class harnessEcal {
public:
  harnessEcal();

  TH2F* energyScaleEB( std::string calibFile );
  TH2F* energyScaleEE( std::string calibFile, int plusOrMinus, bool ixiy = false );

  int deeEE( int ix  , int iy , int iside  );
  int harnessEE( int ix  , int iy, int iside   );
  int ietaEE( float ix, float iy, int side );
  
  /// ad'hoc function in EB, just for IJazZ
  int deeEB(     float ieta );
  int harnessEB( float ieta, float iphi );

  int ieta(    float scEta, float ix, float iy );
  int dee(     float scEta, float ix, float iy );
  int harness( float scEta, float ix, float iy );
  int ix( float scEta, float ix );
  int iy( float scEta, float iy );

  /// vice et versa
  static float ietaFromDeeEB(float dee );
  static float iphiDeeAndHarnessEB( float dee, float harness );
private:
  int _eering_harness_EEm[150][150];
  int _eering_dee_EEm[150][150];
  int _eering_harness_EEp[150][150];
  int _eering_dee_EEp[150][150];
  int _eering_eta[150][150];

  std::string _defHarnessEE;
  std::string _defEtaRingEE;
};

IJazZAxisND<double> combineEtaScale(  std::string  ecalFilePattern );

double ijazzEtaScaleFun(double *x, double *p);

std::vector<TLine*> ecalModuleEdges(float ymin, float ymax, bool IETA = true );
void addPlotEmbellishement( TH1F *frame, std::string plotTitle, bool addLineOne = false ); 
void printEtaScaleToFile( std::string calibOutput, IJazZAxisND<double> &response, float r9 );
void printEtaScaleToFile_EE5x5( std::string calibOutput, IJazZAxisND<double> &response );

void ijazzEtaScaleFitter( TGraphErrors *gr, TF1 *f1 );

void etaScaleFromZee_AnaFit( std::string dataScaleFile, std::string  mcScaleFile, bool saveScale = false );
void plotIOVFile( std::string file);
//boost::tuple<IJazZAxisViewer<double>, std::vector<std::vector<string> >, std::vector<std::vector<float> > > historyPlot( std::string userDir, std::string version, std::string ESorER = "Resp" );

class IJazZhistory {
public:
  std::vector<IJazZAxisViewer<double>  > chViewer;
  std::vector<std::vector<std::string> > chAxisName;
  std::vector<std::vector<float> >       chAxisVal;
  std::string channelTitle( int ich );
};

IJazZhistory historyPlot( std::string userDir, std::string version, std::string ESorER = "Resp" );
IJazZAxisND<double> fwhmOver2( std::string rootFileName );
IJazZAxisND<double> meanResoPdf( std::string rootFileName );
IJazZAxisND<double> maxResoPdf( std::string rootFileName );
IJazZAxisND<double> rmsResoPdf( std::string rootFileName );
IJazZAxisND<double> widthResoPdf( std::string rootFileName );

float convert_ieta_eta( float ieta );
float convert_eta_ieta( float eta  );
float convert_phi_iphi( float phi  );
#endif
