#include "interface/RootUtils.hh"

#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLine.h>

#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
using namespace std;


TChain* addInputFiles( string filelist, string treename ) {
  TChain *ztree = new TChain( treename.c_str() );
  
  ifstream input(filelist.c_str());  
  while ( input.good() && !input.eof() ) {
    string file;
    getline(input,file,'\n');    
    /// comment
    if( file.find("#") != string::npos ) continue;
    if( file.size() < 2 ) continue;
    cout << " adding file: " << file << endl;
    ztree->Add( file.c_str() );
  }
  return ztree;
}

string itostr(int i ) {
  std::stringstream out;
  out << i;
  return out.str();
}

string ftostr(float i ) {
  std::stringstream out;
  out << i;
  return out.str();
}

void addAxisTitles( TH1F *frame, string xTitle, string yTitle, double yMin, double yMax  ) {
  frame->GetXaxis()->SetTitle(xTitle.c_str() ); frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->SetTitle(yTitle.c_str() ); frame->GetYaxis()->CenterTitle();
  if( yMax > yMin ) {
    frame->SetMinimum(yMin);
    frame->SetMaximum(yMax);
  }
}


void addPlotTitle( TH2F *frame, string plotTitle) {
  float xmin = frame->GetXaxis()->GetBinLowEdge(frame->GetXaxis()->GetFirst());
  float ymin = frame->GetYaxis()->GetBinLowEdge( frame->GetYaxis()->GetFirst());; 
  float ymax = frame->GetYaxis()->GetBinUpEdge( frame->GetYaxis()->GetLast());; 
  TLatex tex;
  tex.SetTextColor(kBlack); tex.SetTextSize(0.05); 
  tex.DrawLatex(xmin*0.998,ymax+(ymax-ymin)*0.02, plotTitle.c_str() );
}

void addPlotTitle( TH1F *frame, string plotTitle) {
  float xmin = frame->GetXaxis()->GetBinLowEdge(frame->GetXaxis()->GetFirst());
  float ymin = frame->GetMinimum();
  float ymax = frame->GetMaximum();
  TLatex tex;
  tex.SetTextColor(kBlack); tex.SetTextSize(0.05); 
  tex.DrawLatex(xmin*0.998,ymax+(ymax-ymin)*0.02, plotTitle.c_str() );
}


void EtaHistoStyle(TGraphErrors *g, float xmin, float xmax, float ymin, float ymax, int color, string ytitle, string xtitle) {
  g->GetHistogram()->GetXaxis()->SetRangeUser(xmin,xmax);
  g->GetHistogram()->GetXaxis()->SetTitle(xtitle.c_str());
  g->GetHistogram()->GetXaxis()->CenterTitle();
  g->GetHistogram()->GetYaxis()->SetTitle(ytitle.c_str());
  g->GetHistogram()->GetYaxis()->CenterTitle();
  g->GetHistogram()->SetMinimum(ymin);
  g->GetHistogram()->SetMaximum(ymax);
  g->GetHistogram()->SetTitle("");
  g->SetMarkerStyle(8);
  g->SetMarkerSize(0.8);
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetLineWidth(2);
}



int FindXindexInGraph( float x, TGraphErrors* g ) {
  float dxmin = 999999;
  int ipmin = -1;
  for( int ip = 0 ; ip < g->GetN(); ip++ ) {
    float dx = fabs(g->GetX()[ip]-x);
    if( dx < dxmin ) {
      dxmin = dx;
      ipmin = ip;
    }
  }
  
  return ipmin;
}


TGraphErrors *AverageGraph( TGraphErrors *g1, TGraphErrors *g2, vector<double> *correl) {
  TGraphErrors *gres = 0;
  if( g1->GetN() != g2->GetN() ) {
    cout << " g1 and g2 have different numbers of point... I will crash" << endl;
    return (TGraphErrors*) gres;
  }
  gres = (TGraphErrors*) g1->Clone(); 

  for( int ip = 0 ; ip < g1->GetN(); ++ip ) {

    double x1,y1,ex1,ey1;
    double x2,y2,ey2;
    
    g1->GetPoint(ip,x1,y1);
    g2->GetPoint(ip,x2,y2);
    ex1 = g1->GetErrorX(ip);
    ey1 = g1->GetErrorY(ip);
    ey2 = g2->GetErrorY(ip);
    
    float cor = 0;
    if( correl != 0 && correl->size() > 0 ) cor = correl->at(ip);
    double d = 0;
    float ed = 0;
    if( ey1 > 0.00001 && ey2 > 0.00001 ) {
      d =
	( y1*(ey2*ey2) + y2*(ey1*ey1) - cor*ey1*ey2*(y1+y2) ) /
	(  1*(ey2*ey2) +  1*(ey1*ey1) - cor*ey1*ey2*( 1+ 1) ) ;
      ed = sqrt( (ey1*ey1)*(ey2*ey2)*(1-cor*cor) /
		 (  1*(ey2*ey2) +  1*(ey1*ey1) - cor*ey1*ey2*( 1+ 1) ) ) ;
      
      /* d  = ( y1/(ey1*ey1) + y2/(ey2*ey2) ) / (1./(ey1*ey1) + 1./(ey2*ey2)); */
      /* ed = sqrt( 1. / (1./(ey1*ey1) + 1./(ey2*ey2)) ); */
    }
    gres->SetPoint(     ip, x1, d);
    gres->SetPointError(ip,ex1,ed);
  }

  return gres;
}


TGraphErrors *DivideGraph( TGraphErrors *g1, TGraphErrors *g2) {
  TGraphErrors *gres = 0;
  if( g1->GetN() != g2->GetN() ) {
    cout << " g1 and g2 have different numbers of point... I will crash" << endl;
    cout << " N[g1] = " << g1->GetN()  << " vs N[g2] = " << g2->GetN() << endl;
    for( int in = 0; in < g2->GetN(); in++ ) {
      cout << " X[i1] = " << g1->GetX()[in] << " ; X[i2] = " <<  g2->GetX()[in] << endl;
    }
    
    return (TGraphErrors*) gres;
  }

  gres = (TGraphErrors*) g1->Clone();   
  for( int ip = 0 ; ip < g1->GetN(); ++ip ) {

    double x1,y1,ex1,ey1;
    double x2,y2,ey2;
    
    g1->GetPoint(ip,x1,y1);
    g2->GetPoint(ip,x2,y2);
    ex1 = g1->GetErrorX(ip);
    ey1 = g1->GetErrorY(ip);
    ey2 = g2->GetErrorY(ip);
    
    double d = y1/y2;
    double ed = d*sqrt(ey1*ey1/(y1*y1) + ey2*ey2/(y2*y2));
    if( fabs( y1 ) < 0.00001 || fabs(y2) < 0.00001 ) { d = 0; ed = 0; }
    gres->SetPoint(     ip, x1, d);
    gres->SetPointError(ip,ex1,ed);
  }

  return gres;
}

TGraphErrors* DivideGraphPerAverage(TGraphErrors* gr1,TGraphErrors* gr2 ) {

  TGraphErrors *grout = new TGraphErrors( gr1->GetN() );
  for( int ip = 0 ; ip < gr1->GetN() ; ip++ ) {
    float e1 = gr1->GetEY()[ip];
    float e2 = gr2->GetEY()[ip];
    float c1 = gr1->GetY()[ip];
    float c2 = gr2->GetY()[ip];
    

    float w1 = 0.5;
    float w2 = 0.5;
    if( e1 > 0.0001 && e2 > 0.0001 )  {
      w1 = 1./( 1.+(e1/e2*e1/e2) );
      w2 = 1./( 1.+(e2/e1*e2/e1) );
    }

    float av =  (w1*c1+w2*c2);
    float r = c1 / av;
    /// Error... assuming poissonian stat.
    float er = e1*c2 / (av*av);
    if( fabs( av ) < 0.000001 ) {r=1;er=0;}

    grout->SetPoint(     ip,gr1->GetX()[ip] , r);
    grout->SetPointError(ip,gr1->GetEX()[ip],er);
  }

  return grout;

}

TGraphErrors *SubtractQuadGraph( TGraphErrors *g1, TGraphErrors *g2, float correlation ) {
  TGraphErrors *gres = 0;
  int igclone = 1;
  if( g1->GetN() != g2->GetN() ) {
    cout << " g1 and g2 have different numbers of point... This might crash" << endl;
    //    return (TGraphErrors*) gres;
    if( g1->GetN() < g2->GetN() ) igclone = 2;
  }
  //  gres = new TGraphErrors( g1->GetN() );
  if( igclone == 1 ) gres = (TGraphErrors*) g1->Clone(); 
  if( igclone == 2 ) gres = (TGraphErrors*) g2->Clone(); 
  
  for( int ip = 0 ; ip < gres->GetN(); ++ip ) {
    
    double x = gres->GetX()[ip];
    int ip1 = FindXindexInGraph(x,g1);
    int ip2 = FindXindexInGraph(x,g2);

    double x1,y1,ex1,ey1;
    double x2,y2,ey2;
    
    g1->GetPoint(ip1,x1,y1);
    g2->GetPoint(ip2,x2,y2);
    ex1 = g1->GetErrorX(ip1);
    ey1 = g1->GetErrorY(ip1);
    ey2 = g2->GetErrorY(ip2);
    double d = y1*y1 - y2*y2;
    if( d > 0 )  d = +sqrt(fabs(d)); 
    else         d = -sqrt(fabs(d));
    double ed = sqrt(ey1*ey1*y1*y1 + ey2*ey2*y2*y2-2*y1*y2*ey1*ey2*correlation)/fabs(d);
    if( fabs( y1 ) < 0.00001 || fabs(y2) < 0.00001 ) { d = 0; ed = 0; }
    //    if( d <= 0.001 ) ed = 0; 
    gres->SetPoint(     ip, x1, d);
    gres->SetPointError(ip,ex1,ed);
  }

  return gres;
}



TGraphErrors *SubtractQuadGraphs( TGraphErrors *g1, vector<TGraphErrors*> g2 ) {
  TGraphErrors *gres = 0;
  gres = (TGraphErrors*) g1->Clone(); 
  for( int ip = 0 ; ip < gres->GetN(); ++ip ) {
    double x1,y1,ex1,ey1;
    double x2,y2,ey2(1);
    
    double x = gres->GetX()[ip];
    int ip1 = FindXindexInGraph(x,g1);
    g1->GetPoint(ip1,x1,y1);
    ex1 = g1->GetErrorX(ip1);
    ey1 = g1->GetErrorY(ip1);
    
    float os2 = 0;
    for( int imc = int(g2.size()-1) ; imc >= 0; imc-- ) {
      int ip2 = FindXindexInGraph(x,g2[imc]);
      g2[imc]->GetPoint(ip2,x2,y2);
      if( imc == int(g2.size()-1) ) ey2 = g2[imc]->GetErrorY(ip2);      
      os2 += y1*y1 - y2*y2;
    }
    
    double d = os2;
    if( d > 0 )  d = +sqrt(fabs(d)); 
    else         d = -sqrt(fabs(d));
    double ed = sqrt(ey1*ey1*y1*y1 + ey2*ey2*y2*y2)/fabs(d);
    if( fabs( y1 ) < 0.00001 || fabs(y2) < 0.00001 ) { d = 0; ed = 0; }
    if( d <= 0.001 ) ed = 0; 

    gres->SetPoint(     ip, x1, d);
    gres->SetPointError(ip,ex1,ed);
  }

  return gres;
}

TGraphErrors *selectGraphOp( TGraphErrors *g1, TGraphErrors *g2, int sel, float correl ) {
  if     ( sel == 0 ) return g1;
  else if( sel == 1 ) return g2;
  else if( sel == 2 ) return AverageGraph( g1, g2 );
  else if( sel == 3 ) return DivideGraph(  g1, g2 );
  else if( sel == 4 ) return DivideGraphPerAverage(  g1, g2 );
  else if( sel == 5 ) return SubtractQuadGraph( g1, g2, correl );
  else {
    cerr << "TGraphErrors *chooseGraph ERROR, sel " << sel << " not permitted: " <<  endl
	 << "  - sel = 0: return g1 " << endl
	 << "  - sel = 1: return g2 " << endl
	 << "  - sel = 2: return average(g1,g2) " << endl
	 << "  - sel = 3: return g1 / g2 " << endl
	 << "  - sel = 4: return g1 / (g1+g2) " << endl
	 << "  - sel = 5: return \"sqrt(g1*g1 - g2*g2)\" " << endl
	 << "   .... I will bail out soon ... " << endl;
  }
  
  return (TGraphErrors*) 0;
}



///----------------- Histo func -----------------------------///
int irandName = 0;
string randRootName() {
  return "RootNameRand_" + itostr(irandName++);
}

TH1F * yAxisSpread(TH1F * h, float &errExp ) {
  vector<float>  bC;
  vector<float> ebC;
  errExp = 0;
  for( int ix = 1; ix <= h->GetXaxis()->GetNbins(); ix++ ) {
    bC.push_back( h->GetBinContent(ix) );
    ebC.push_back( h->GetBinError(ix) );
    errExp += h->GetBinError(ix);
  }
  errExp /= ebC.size();
  float minimum = *min_element(bC.begin(),bC.end());
  float maximum = *max_element(bC.begin(),bC.end());
  float l = maximum-minimum;
  
  TH1F *hout = new TH1F(randRootName().c_str(),"unknown",20,minimum-l/10.,maximum+l/10.);
  for( unsigned i = 0 ; i < bC.size(); i++ ) hout->Fill(bC[i]);
  
  return hout;
}
