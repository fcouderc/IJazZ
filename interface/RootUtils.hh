#ifndef rootUtils__hh__
#define rootUtils__hh__

#include <string>
#include <vector>

class TChain;
class TGraphErrors;
class TH1F;
class TH2F;

//----- Utils function
TChain* addInputFiles( std::string filelist, std::string treename );
std::string itostr(int i ); 
std::string ftostr(float i ); 

int  FindXindexInGraph( float x, TGraphErrors* g ); 
void EtaHistoStyle( TGraphErrors *g, 
		    float xmin, float xmax, 
		    float ymin, float ymax, 
		    int color, std::string ytitle, std::string xtitle);
TGraphErrors *AverageGraph( TGraphErrors *g1, TGraphErrors *g2, std::vector<double> *correl = 0);
TGraphErrors *DivideGraph( TGraphErrors *g1, TGraphErrors *g2);
TGraphErrors *DivideGraphPerAverage( TGraphErrors *g1, TGraphErrors *g2);
TGraphErrors *SubtractQuadGraph( TGraphErrors *g1, TGraphErrors *g2, float correlation = 0);
TGraphErrors *SubtractQuadGraphs( TGraphErrors *g1, std::vector<TGraphErrors*> g2);

TGraphErrors *selectGraphOp( TGraphErrors *g1, TGraphErrors *g2, int sel, float correl =0 );

void addAxisTitles( TH1F *frame, std::string xTitle, std::string yTitle, double yMin = 0, double yMax = -1 );
void addPlotTitle( TH1F *frame, std::string title );
void addPlotTitle( TH2F *frame, std::string title );

TH1F * yAxisSpread(TH1F * h, float &errExp);
std::string randRootName();

#endif
