#include "interface/IJazZ.hh"
#include "interface/IJazZAxisViewer.hh"
#include "interface/RootUtils.hh"
#include "interface/EcalUtils.hh"

#include <TGraphErrors.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFile.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>

#include <iostream>
#include <vector>
using namespace std;


int ijazzEtaScale_nPar = 7;
double xmean = 0;

double ijazzEtaScaleFun( double *x, double *p ) {
  double xx = 1;
  double xmean = p[ijazzEtaScale_nPar-1];

  // xx*=xx;
  double res = 0;
  //  cout << " npar = " << npar << endl;
  
  //  for( int i = ijazzEtaScale_nPar-1; i >=0 ; i-- ) res = p[i] + res*xx;
  for( int i = 0; i < ijazzEtaScale_nPar-1; i++ ) {
    res += p[i] * xx;
    xx *= (x[0]-xmean)/10.;    
  }

  
  return res + 1;
}

void ijazzEtaScaleFitter( TGraphErrors *gr, TF1 *f1 ) {

  ijazzEtaScale_nPar = f1->GetNpar();
 
  for( int ip = 0; ip < ijazzEtaScale_nPar-1; ip++ ) {
    f1->SetParameter(ip,0.001);
    f1->SetParLimits(ip,-0.5,+0.5);
    f1->SetLineColor( gr->GetLineColor()-1 );
    f1->SetLineWidth(2);
  }
  double x0(-1),x1(1);
  f1->GetRange(x0,x1);
  xmean = (x0+x1)/2.;
  //f1->AddParameter( "xmean", xmean );
  f1->SetParName(   ijazzEtaScale_nPar-1,"xmean");
  f1->FixParameter( ijazzEtaScale_nPar-1, xmean);
  cout << "  npar F1: " << f1->GetNpar() << endl;
  gr->Fit(f1, "mhe r");
}

  
//  string dataVersion = _version;
//  if( mcVersion == "notdefined" ) mcVersion = _version;

void etaScaleFromZee_AnaFit( string dataScaleFile, string  mcScaleFile, bool saveScale ) {
  string plotTitle = "CMS 2012 prelim., #sqrt{s}= 8 TeV, L = 19.5 fb^{-1}";
  vector<TGraphErrors*> grES[2];
  int colors[2][5] = { {kBlue+1, kAzure+1 , kBlue , kAzure -2, kAzure -2 },  /// Data colors
		       {kRed +1, kOrange+1, kRed  , kOrange-1, kRed }   /// MC colors
		       };  

  IJazZAxisViewer<double> viewer[2];
  viewer[0] = IJazZAxisViewer<double>(combineEtaScale(dataScaleFile));
  viewer[1] = IJazZAxisViewer<double>(combineEtaScale(  mcScaleFile));
  IJazZAxis<double> *xeta = viewer[0].getAxisND()->getAxis(0);
  for( int ismc = 0 ; ismc <= 1; ismc++ ) {
    IJazZAxis<double> *xR9 = viewer[ismc].getAxisND()->getAxis(1);
    for( xR9->itBegin(); !xR9->itEnd(); xR9->itForward() ) {
      grES[ismc].push_back( viewer[ismc].proj1D(0,xR9->itBinCenter() ) );
      IJazZ::etaHistoStyle(  grES[ismc][xR9->itBinNumber()], xeta, colors[ismc][xR9->itBinNumber()], "r( i#eta )" );
    }
  }

  /// incl scale
  int iIncl = grES[0].size();
  for( int ismc = 0 ; ismc <= 1; ismc++ ) {
    grES[ismc].push_back( AverageGraph( grES[ismc][0], grES[ismc][1] ) );
    IJazZ::etaHistoStyle( grES[ismc][grES[ismc].size()-1], xeta,colors[ismc][grES[ismc].size()-1], "r( i#eta )_{incl}" );
  }

  /// gold ratio
  int iGold = grES[0].size();
  for( int ismc = 0 ; ismc <= 1; ismc++ ) {
    grES[ismc].push_back( DivideGraphPerAverage( grES[ismc][1], grES[ismc][0] ) );
    IJazZ::etaHistoStyle( grES[ismc][grES[ismc].size()-1], xeta,colors[ismc][grES[ismc].size()-1], "r( i#eta )_{gold} / r( i#eta )_{incl}" );
  }

  /// bad ratio
  int iBad = grES[0].size();
  for( int ismc = 0 ; ismc <= 1; ismc++ ) {
    grES[ismc].push_back( DivideGraphPerAverage( grES[ismc][0], grES[ismc][1] ) );
    IJazZ::etaHistoStyle( grES[ismc][grES[ismc].size()-1], xeta, colors[ismc][grES[ismc].size()-1], "r( i#eta )_{bad} / r( i#eta )_{incl}" );
  }


  TH1F *frame = grES[0][0]->GetHistogram();  
  string pdfout = dataScaleFile + ".etaScalePlots.pdf";

  unsigned nCanvas = grES[0].size()+3;
  vector<TCanvas*> can;
  for( unsigned ic = 0 ; ic < nCanvas; ic++ ) {
    can.push_back(  new TCanvas( ("canvas" + itostr(ic)).c_str(),"canvas",800,600) );
    can[ic]->SetTopMargin(0.08);
  }

  can[0]->Print( string(pdfout + "[").c_str());
  /// 
  unsigned ican = 0;
  for( ican = 0 ; ican <  grES[0].size(); ican++ ) {
    can[ican]->cd();
    grES[0][ican]->Draw("AP"); grES[1][ican]->Draw("P");
    addPlotEmbellishement( frame, plotTitle);
  } 
  

  TGraphErrors *grScaleIncl = DivideGraph(grES[0][iIncl], grES[1][iIncl] );
  TGraphErrors *grScaleGold = DivideGraph(grES[0][iIncl], grES[1][iIncl] );
  TGraphErrors *grScaleBad  = DivideGraph(grES[0][iIncl], grES[1][iIncl] );
  TGraphErrors *grScaleGoldWithErrors = DivideGraph(grES[0][iIncl], grES[1][iIncl] );
 
  /// incl -> gold fit

  float xminEB = 0.0;
  float xmaxEB = 1.49;
  float xminEE = 1.49;
  float xmaxEE = 2.5;
  if(  viewer[0].getAxisND()->getAxis(0)->getName() == "IEta" ) {
    xminEB =   0.5;
    xmaxEB =  84.5;
    xminEE =  87.5;
    xmaxEE = 120.5;
  }

  TF1 *fGold[] = { new TF1("fGold1",ijazzEtaScaleFun, -xmaxEB,-xminEB, 7),
		   new TF1("fGold2",ijazzEtaScaleFun, +xminEB,+xmaxEB, 7),
		   new TF1("fGold3",ijazzEtaScaleFun, -xmaxEE,-xminEE, 7),
		   new TF1("fGold4",ijazzEtaScaleFun, +xminEE,+xmaxEE, 7)};

  TF1 *fBad[] = { new TF1("fBad1",ijazzEtaScaleFun, -xmaxEB,-xminEB, 7),
		  new TF1("fBad2",ijazzEtaScaleFun, +xminEB,+xmaxEB, 7),
		  new TF1("fBad3",ijazzEtaScaleFun, -xmaxEE,-xminEE, 7),
		  new TF1("fBad4",ijazzEtaScaleFun, +xminEE,+xmaxEE, 7)};

  TGraphErrors *ratioGold = DivideGraph( grES[0][iGold], grES[1][iGold]);
  IJazZ::etaHistoStyle(ratioGold, xeta, kOrange-1, "r_{gold} / r_{incl} : ratio Data/MC" );
  can[ican++]->cd(); ratioGold->Draw("AP");
  for( unsigned iec = 0 ; iec < 4; iec++ ) {
    fGold[iec]->SetLineColor( kYellow+1);
    fGold[iec]->SetParLimits(0,-0.10,0.10);
    ratioGold->Fit( fGold[iec], "mhe r" );
    fGold[iec]->Draw("same");
  }
  addPlotEmbellishement( frame, plotTitle, 1);

  TGraphErrors *ratioBad = DivideGraph( grES[0][iBad], grES[1][iBad]);
  IJazZ::etaHistoStyle(ratioBad, xeta, kOrange, "r_{bad} / r_{incl} : ratio Data/MC" );
  can[ican++]->cd(); ratioBad->Draw("AP");
  for( unsigned iec = 0 ; iec < 4; iec++ ) {
    fBad[iec]->SetLineColor( kYellow+1);
    fBad[iec]->SetParLimits(0,0.90,1.10);
    ratioBad->Fit( fBad[iec], "mhe r" );
    fBad[iec]->Draw("same");
  }
  addPlotEmbellishement( frame, plotTitle, 1);


  if(  viewer[0].getAxisND()->getAxis(0)->getName() == "IEta" ||
       viewer[0].getAxisND()->getAxis(0)->getName() == "SCEta"
       ) {
    for( int ip = 0; ip < grScaleGold->GetN(); ip++ ) {
      float x = grScaleIncl->GetX()[ip];
      float y1 = grScaleIncl->GetY()[ip];
      float y2 = grScaleIncl->GetY()[ip];
      float ey = grScaleIncl->GetEY()[ip];
      if     ( x <= -xmaxEB ) y1 *= fGold[2]->Eval(x);
      else if( x <= -xminEB ) y1 *= fGold[0]->Eval(x);
      else if( x <= +xmaxEB ) y1 *= fGold[1]->Eval(x);
      else if( x <= +120.5  ) y1 *= fGold[3]->Eval(x);
      
      if     ( x <= -xmaxEB ) y2 *= fBad[2]->Eval(x);
      else if( x <= -xminEB ) y2 *= fBad[0]->Eval(x);
      else if( x <= +xmaxEB ) y2 *= fBad[1]->Eval(x);
      else if( x <= +120.5  ) y2 *= fBad[3]->Eval(x);
      grScaleGold->SetPoint( ip, x, y1 );
      grScaleBad ->SetPoint( ip, x, y2 );
      grScaleGoldWithErrors->SetPoint( ip, x, y1 );
      if( ey > 0 )  grScaleGold->SetPointError( ip, 0., 0. );
      if( ey > 0 )  grScaleBad ->SetPointError( ip, 0., 0. );
      if( ey > 0 )  grScaleGoldWithErrors->SetPointError( ip, 0., ey*1.20 );

      //      if( (fabs( x ) < 88.5 && fabs(x) > 83.5) ||  fabs( x ) <   xminEB+0.000001 || fabs( x ) > 118.5 )  {      
      if( false ) {
	grScaleGold->GetY()[ip]  = 1;
	grScaleGold->GetEY()[ip] = 0.05;
	grScaleBad->GetY()[ip]  = 1;
	grScaleBad->GetEY()[ip] = 0.05;
	grScaleGoldWithErrors->GetY()[ip]  = 1;
	grScaleGoldWithErrors->GetEY()[ip] = 0.05;
	grScaleIncl->GetY()[ip]  = 1;
	grScaleIncl->GetEY()[ip] = 0.05;
      }
    }
  } else {
    /// here we do not want the ratio to the inclusieve but the absolute ones to get MC/Data
    iGold = 1;
    iBad  = 0;
    grScaleIncl = DivideGraph(grES[0][iIncl], grES[1][iIncl] );
    grScaleGold = DivideGraph(grES[0][iGold], grES[1][iGold] );
    grScaleBad  = DivideGraph(grES[0][iBad ], grES[1][iBad ] );    
    grScaleGoldWithErrors = DivideGraph(grES[0][iGold], grES[1][iGold] );
  }


  IJazZ::etaHistoStyle(grScaleIncl, xeta, kGreen+1, "scale incl." );
  IJazZ::etaHistoStyle(grScaleGold, xeta, kGray+1 , "scale gold " );
  IJazZ::etaHistoStyle(grScaleBad , xeta, kYellow+1 , "scale bad " );

  can[nCanvas-1]->cd();
  grScaleIncl->Draw("AP");
  grScaleGold->Draw("P");
  grScaleBad->Draw("P");
  addPlotEmbellishement( frame, plotTitle, 1);



  /// Add legends
  vector<TLegend*> leg;
  for( unsigned ic = 0 ; ic < nCanvas; ic++ ) {
    leg.push_back( new TLegend(0.32,0.75,0.75,0.89) );
    leg[ic]->SetBorderSize(0);
    leg[ic]->SetFillColor(0);
    string legTitle;
    if( ic == 0 ) legTitle = " R_{9} < 0.94";
    if( ic == 1 ) legTitle = " R_{9} #geq 0.94";
    if( ic == 2 ) legTitle = " Incl. ";
    if( ic == 3 ) legTitle = " Gold / Incl ";
    if( ic == 4 ) legTitle = " Bad / Incl ";
    if( ic == 5 ) legTitle = " Gold / Incl: Data/MC ";
    if( ic == 6 ) legTitle = " Bad / Incl: Data/MC ";
    if( ic == 7 ) legTitle = " #eta scale ";


    if( ic < 5 ) {
      for( int ismc = 0 ; ismc <= 1; ismc++ ) {
	string dataORmc = "Data: ";
	if( ismc == 1 ) dataORmc = "MC  : " ;      
	leg[ic]->AddEntry( grES[ismc][ic], (dataORmc+legTitle).c_str(), "lp" );
      }
    }
    if( ic == 5 ) leg[ic]->AddEntry( ratioGold, legTitle.c_str(), "lp" );
    if( ic == 6 ) leg[ic]->AddEntry( ratioBad , legTitle.c_str(), "lp" );
    if( ic == 7 ) leg[ic]->AddEntry( grScaleIncl, (legTitle+" incl.").c_str(),"lp");
    if( ic == 7 ) leg[ic]->AddEntry( grScaleGold, (legTitle+" gold ").c_str(),"lp");
    if( ic == 7 ) leg[ic]->AddEntry( grScaleBad , (legTitle+" bad ").c_str(),"lp");
    
    can[ic]->cd(); leg[ic]->Draw();    
    can[ic]->Print( pdfout.c_str() );
  }

  

  can[0]->Print( string(pdfout + "]").c_str() );
  cout << " All Plots saved in: " << pdfout << endl;

  IJazZAxisND<double> dataresp = *(viewer[0].getAxisND());
  for( int ip = 0 ; ip < grScaleGold->GetN(); ip++ ) {
    vector<double> xx; xx.resize(2,0);
    xx[0] = grScaleGold->GetX()[ip];
    xx[1] = 0.95;
    dataresp.value(xx) = grScaleGold->GetY()[ip] ;
    dataresp.error(xx) = grScaleGold->GetEY()[ip];
    xx[0] = grScaleBad->GetX()[ip];
    xx[1] = 0.85;
    dataresp.value(xx) = grScaleBad->GetY()[ip]; 
    dataresp.error(xx) = grScaleBad->GetEY()[ip];
  }
  dataresp.saveToFile( dataScaleFile + ".etaScaleDataOverMC" );
  cout << " Saved EtaScale to: " <<  dataScaleFile + ".etaScaleDataOverMC" << endl;

  if( ! saveScale ) return;
  string calibIncl = "scaleFactorsZee_allR9.txt" ;
  string calibGold = "scaleFactorsZee_highR9.txt" ;
  printEtaScaleToFile( calibIncl, dataresp, 0.80);
  printEtaScaleToFile( calibGold, dataresp, 0.95 );

}



///------------ eta scale print out ------------///
#include <sstream>
#include <fstream>
void printEtaScaleToFile( string calibOutput, IJazZAxisND<double> &response, float r9 ) {
  harnessEcal h;

#define netaRING_EB 85
#define netaRING_EE 35
  string calibSkeleton = "etc/geom/skeleton.txt";
  ifstream skeleton(calibSkeleton.c_str());
  ofstream outputfile( calibOutput.c_str() );
  while ( skeleton.good() && !skeleton.eof() ) {
    string line;
    getline(skeleton,line,'\n');
    istringstream isstream(line);
    int ix,iy,iz;
    float ic, eic;
    isstream >> ix >> iy >> iz >> ic >> eic;
    ic = 1; eic = 0;
    int ieta = -1;
    if( iz ==  0 ) ieta = ix ;
    else           ieta = h.ietaEE(ix,iy,iz);

    if( ieta < -120 || ieta > 120 ) { ic = 1;  eic = 0.; }
    else {
      vector<double> x; x.resize(2);
      x[0] = ieta; x[1] = r9;
      ic      = 1./response.value(x);
      eic     = response.error(x);
    }
    outputfile  << ix << " " << iy << " " << iz << " " << ic     << " " << eic     << endl;
  }
  outputfile.close();
  skeleton.close();
}



void printEtaScaleToFile_EE5x5( string calibOutput, IJazZAxisND<double> &response ) {
  harnessEcal h;

#define netaRING_EB 85
#define netaRING_EE 35
  string calibSkeleton = "etc/skeleton.txt";
  ifstream skeleton(calibSkeleton.c_str());
  ofstream outputfile( calibOutput.c_str() );
  while ( skeleton.good() && !skeleton.eof() ) {
    string line;
    getline(skeleton,line,'\n');
    istringstream isstream(line);
    int ix,iy,iz;
    float ic, eic;
    isstream >> ix >> iy >> iz >> ic >> eic;
    ic = 1; eic = 0;
    int ieta = -1;
    if( iz ==  0 ) ieta = ix ;
    else           ieta = h.ietaEE(ix,iy,iz);

    if     ( ieta < -120 || ieta > 120 ) { ic = 1;  eic = 1.; }
    else if( iz == 0 ) { ic = 1; eic = 1;} 
    else {
      vector<double> x; x.resize(3);
      x[0] = ieta; x[1] = ix; x[2] = iy;
      ic      = 1./response.value(x);
      eic     = response.error(x);
    }
    outputfile  << ix << " " << iy << " " << iz << " " << ic     << " " << eic     << endl;
  }
  outputfile.close();
  skeleton.close();
}
