#ifdef __CLING__

//#pragma link off all globals;
//#pragma link off all classes;
//#pragma link off all functions;

#pragma link C++ class IJazZAxisViewer<double>+;
#pragma link C++ class IJazZAxisND<double>+;
#pragma link C++ class IJazZAxis<double>+;
#pragma link C++ class ijazzAxisND;
#pragma link C++ class ijazzViewer;

#pragma link C++ class harnessEcal+;
#pragma link C++ class IJazZhistory+;

//#pragma link C++ class vector<TH1F*>+;
//#pragma link C++ class vector<TH2F*>+;
//#pragma link C++ class vector<double>+;

#pragma link C++ function combineEtaScale( std::string );
#pragma link C++ function etaScaleFromZee_AnaFit( std::string, std::string, bool );
#pragma link C++ function itostr( int ); 
#pragma link C++ function EtaHistoStyle( TGraphErrors*,float,float,float,float,int, std::string, std::string);
#pragma link C++ function AverageGraph(TGraphErrors*,TGraphErrors*, std::vector<double> *) ;
#pragma link C++ function DivideGraph( TGraphErrors *, TGraphErrors *);
#pragma link C++ function DivideGraphPerAverage( TGraphErrors *, TGraphErrors *);
#pragma link C++ function SubtractQuadGraph(  TGraphErrors *, TGraphErrors *, float);
#pragma link C++ function SubtractQuadGraphs( TGraphErrors *, std::vector<TGraphErrors *>);
#pragma link C++ function selectGraphOp( TGraphErrors *, TGraphErrors *, int, float);
#pragma link C++ function addAxisTitles( TH1F *, std::string, std::string, double, double);
#pragma link C++ function addPlotTitle( TH1F *, std::string);
#pragma link C++ function addPlotTitle( TH2F *, std::string);
#pragma link C++ function randRootName();
#pragma link C++ function ijazzEtaScaleFun(double *x, double *p);
#pragma link C++ function ijazzEtaScaleFitter( TGraphErrors *gr, TF1 *f1 );
#pragma link C++ function plotIOVFile(std::string);
#pragma link C++ function historyPlot(std::string,std::string,std::string);
#pragma link C++ function printEtaScaleToFile_EE5x5( std::string, IJazZAxisND<double>&);
#pragma link C++ function printEtaScaleToFile(       std::string, IJazZAxisND<double>&);
#pragma link C++ function ecalModuleEdges(float,float, bool);
#pragma link C++ function fwhmOver2(    std::string );
#pragma link C++ function meanResoPdf(  std::string );
#pragma link C++ function maxResoPdf(   std::string );
#pragma link C++ function rmsResoPdf(   std::string );
#pragma link C++ function widthResoPdf( std::string );
#pragma link C++ function convert_ieta_eta(float);
#pragma link C++ function convert_eta_ieta(float);

#endif

