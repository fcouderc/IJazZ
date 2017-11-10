#include "interface/EcalUtils.hh"
#include "interface/RootUtils.hh"
#include "interface/IJazZAxisViewer.hh"

#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>

#include <vector>
#include <string>
using namespace std;

vector<TGraphErrors*> ShervinOS(void);

void EcalPaperPlot(string dirafs = "./", bool doResp = false) {
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadLeftMargin(0.12);
  string outputpdf = "EcalPlotsMoriond2013.pdf";
  outputpdf = "EcalPlotsCompare2016.pdf";

  string version1 = "prompt";
  string version2 = "IC_piZero_v1";
  string version3 = "IC_piZero_v2";
  string version4 = "ref";
  string userDir = "user";

  userDir  = "2016_July";


  string dir = "IJazZ_output/FitResultsClassic/" + userDir + "/";
  vector<string> resofiles, title;
  vector<int>    colors;
  dir = dirafs + dir;

  string xtitle = "SuperCluster | #eta |";
  string ytitle = "#sigma_{E} / E";
  string ESorER = ".fittedReso";
  float y0 = 0;
  float y1 = 0.06;
  if( doResp ) {
    ESorER = ".fittedResp";
    ytitle = "r(#eta)";
    y0 = 0.96;
    y1 = 1.04;
  }

  float yMaxOS = 0.03;
  string dir2 = "IJazZ_output/FitResults/user/";
  dir2 = dir;
  resofiles.push_back( dir + "0-999999/IJazZ_EB_MC_"   + version1 + ".root" ); colors.push_back( kRed+1  ); title.push_back("MC");
  resofiles.push_back( dir + "0-999999/IJazZ_EB_Data_" + version1 + ".root" ); colors.push_back( kBlack); title.push_back("Data prompt");

  if( version2 != "notdefined" ){
    resofiles.push_back( dir2 + "0-999999/IJazZ_EB_Data_" + version2 + ".root" ); colors.push_back( kGray+1  ); title.push_back("Data ref ");
  }
  if( version3 != "notdefined" ){
    resofiles.push_back( dir2 + "0-999999/IJazZ_EB_Data_" + version3 + ".root" ); colors.push_back( kGreen+2  ); title.push_back("Data ReReco + Zee");
  }
  if( version4 != "notdefined" ){
    resofiles.push_back( dir + "0-999999/IJazZ_EB_Data_" + version4 + ".root" ); colors.push_back( kAzure    ); title.push_back("new MC runD");
  }

  // resofiles.push_back( dir + "0-999999/IJazZ_EE_Data_" + version1 + ".root" ); colors.push_back( kBlack  ); title.push_back("Data Final 2012");
  // resofiles.push_back( dir + "0-999999/IJazZ_EE_MC_"   + version1 + ".root" ); colors.push_back( kRed+1  ); title.push_back("MC std madgraph");


  int imc    =  resofiles.size()-1; imc   = 0;
  int idata  =  resofiles.size()-2; idata = 1;// idata = 0;
  int imc2   = 0; 
  int idata2 = 1;


  TLegend *leg[2] = { new TLegend(0.12,0.63,0.50,0.9),
		      new TLegend(0.12,0.63,0.50,0.9) };
  for( int ir9 = 0 ; ir9 < 2; ir9++ ) {
    leg[ir9]->SetBorderSize(0);
    leg[ir9]->SetFillColor(0);
  }

  vector<TGraphErrors*> gr[2];

  double r9[2] = {0.85, 0.95};
  string titlePlot = "CMS 2012 preliminary: L = 19.5 fb^{-1}, #sqrt{s} = 8TeV";
  //string titlePlot = "CMS 2011: L = 4.98 fb^{-1}, #sqrt{s} = 7TeV";
  //  string titleR9[2]   = { ", inclusive", ", R_{9} #geq 0.94"};
  string titleR9[2]   = { ", R_{9} < 0.94", ", R_{9} #geq 0.94"};

  TCanvas *c[2] = { new TCanvas("can_lowR9","low  R9",800,600),
		    new TCanvas("can_higR9","high R9",800,600) };
  for( int ir9 = 0 ; ir9 < 2; ir9++ ) c[ir9]->SetTopMargin(0.085);


  std::vector<TLine*> modEdges = ecalModuleEdges(y0, y1, false);

  TBox *box = new TBox(1.45,y0,1.560,y1);
  box->SetLineWidth(0);
  box->SetFillColor(kGray);


  for( unsigned i = 0 ; i < resofiles.size(); i++ ) {
    for( int ir9 = 0 ; ir9 < 2 ;ir9++ ) {
    c[ir9]->SetGridy();

      c[ir9]->cd();
      if( i > 10   ) {
	if( ESorER.find("Resp") != string::npos )  gr[ir9].push_back( IJazZAxisViewer<double>( maxResoPdf(   resofiles[i] ) ).proj1D(0,r9[ir9]) );
	if( ESorER.find("Reso") != string::npos )  gr[ir9].push_back( IJazZAxisViewer<double>( widthResoPdf( resofiles[i] ) ).proj1D(0,r9[ir9]) );

	//	if( ir9 == 0 ) gr[ir9].push_back( IJazZAxisViewer<double>( maxResoPdf( resofiles[i] ) ).proj1D(0,r9[ir9]) );
	//	else           gr[ir9].push_back( IJazZAxisViewer<double>( maxResoPdf( resofiles[i] ) ).proj1D(0,r9[ir9]) );
      } else {
	if( ir9 == 0 ) gr[ir9].push_back( IJazZAxisViewer<double>( combineEtaScale( resofiles[i] + ESorER) ).proj1D(0,r9[ir9]) );
	else           gr[ir9].push_back( IJazZAxisViewer<double>( combineEtaScale( resofiles[i] + ESorER) ).proj1D(0,r9[ir9]) );
      }
      EtaHistoStyle( gr[ir9][i], 0., 2.5, y0, y1, colors[i], ytitle, xtitle );
      for( int ip = 0 ; ip < gr[ir9][i]->GetN(); ip++ ) if( fabs(gr[ir9][i]->GetX()[ip]-1.5) < 0.01 ) gr[ir9][i]->RemovePoint(ip);
      

      string optDraw = "P"; if( i == 0 )  optDraw = "AP";   
      gr[ir9][i]->Draw(optDraw.c_str());
      leg[ir9]->AddEntry( gr[ir9][i], (title[i]+titleR9[ir9]).c_str(),"lp");
      addPlotTitle( gr[ir9][i]->GetHistogram(), titlePlot.c_str() );
      
    }
  }
  for( int ir9 = 0 ; ir9 < 2; ir9++ ) {
    c[ir9]->cd();
    //    box->Draw();
    for( unsigned m = 0 ; m < modEdges.size(); m++ ) modEdges[m]->Draw("same");
  }



  ///------------------------ plot oversmearings ---------------------------------------///
  TCanvas *cOS[] = { new TCanvas("canOS_lowR9","low  R9",800,600) ,
		     new TCanvas("canOS_highR9","high R9",800,600) };


  for( int ir9 = 0 ; ir9 < 2; ir9++ ) cOS[ir9]->SetTopMargin(0.085);
  std::vector<TLine*> modEdgesOS = ecalModuleEdges(0, yMaxOS, false);
  int colorsOld[] = { kAzure -1, kRed - 3 };
  int colorsNew[] = { kBlue + 2, kRed + 2 };
  TLegend *legOS =  new TLegend(0.12,0.63,0.50,0.9);
  legOS->SetFillColor(0); legOS->SetBorderSize(0);
  std::vector<TGraphErrors*> shOS = ShervinOS();
  for( int ir9 = 0 ; ir9 < 2; ir9++ ) {
    cOS[ir9]->SetGridy();
    cOS[ir9]->cd();
    //    box->Draw();
    
    if( !doResp && idata >= 0  ) {
      TGraphErrors *grSub = SubtractQuadGraph( gr[ir9][idata], gr[ir9][imc] );
      EtaHistoStyle( grSub, 0., 2.5, 0.0, yMaxOS, kGreen + 1, ytitle, xtitle );
      grSub->SetFillColor( colorsOld[ir9] );
      grSub->SetLineColor( colorsOld[ir9] );
      grSub->SetFillStyle( 3004 );
      if( ir9 <= 1 ) grSub->Draw("AE3");
      else           grSub->Draw("E3");
      grSub->Draw("LX");
      shOS[ir9]->Draw("P");
      if( ir9 == 0 ) legOS->AddEntry( shOS[ir9], "Shervin  R_{9} < 0.94","LP");
      if( ir9 == 0 ) legOS->AddEntry( grSub    , "OldIJazZ R_{9} < 0.94","LF");
      if( ir9 == 1 ) legOS->AddEntry( shOS[ir9], "Shervin  R_{9} #geq 0.94","LP");
      if( ir9 == 1 ) legOS->AddEntry( grSub    , "OldIJazZ R_{9} #geq 0.94","LF");
    }


    if( !doResp && idata2 >= 0  ) {
      TGraphErrors *grSub = SubtractQuadGraph( gr[ir9][idata2], gr[ir9][imc2] );
      EtaHistoStyle( grSub, 0., 2.5, 0.0, 0.04, kGreen + 1, ytitle, xtitle );
      grSub->SetFillColor( colorsNew[ir9] );
      grSub->SetLineColor( colorsNew[ir9] );
      grSub->SetMarkerColor( colorsNew[ir9] );
      grSub->SetLineStyle( kDashed );

      for( int ip = 0 ; ip < grSub->GetN(); ip++ ) {
	grSub->GetEX()[ip] = 0.1;
	//	grSub->GetX()[ip]  = grSub->GetX()[ip]-0.05;
      }
      grSub->SetFillStyle( 3005 );
      //      grSub->Draw("E3");
      grSub->Draw("LX");
      //      grSub->Draw("P>");
      if( ir9 == 0 ) legOS->AddEntry( grSub    , "NewIJazZ R_{9} < 0.94","LP");
      if( ir9 == 1 ) legOS->AddEntry( grSub    , "NewIJazZ R_{9} #geq 0.94","LP");
    }
    for( unsigned m = 0 ; m < modEdges.size(); m++ ) modEdgesOS[m]->Draw("same");
  }

  TCanvas *cOSImprovment = 0; 
  if( idata2 > 0 && ESorER == ".fittedReso" ) {
    int cr9[2] = { kAzure, kRed+2 };
    //// ---- improvement in oversmearing ----- /////  
    cOSImprovment = new TCanvas("cOverSMC5","Improvement",700,500);
    TLine *l = new TLine(0,0.0,2.5,0.0); l->SetLineWidth(2);
    for( int ir9 = 0 ; ir9 < 2 ;ir9++ ) {
      cOSImprovment->cd();
      cOSImprovment->SetGridy();
      
      TGraphErrors *gtmp = SubtractQuadGraph( gr[ir9][idata2], gr[ir9][idata], 1 ) ;
      EtaHistoStyle( gtmp, 0., 2.5, -0.03, 0.03, cr9[ir9], "OverSmearing gain", xtitle );
      gtmp->SetFillColor(cr9[ir9]);
      if( ir9 == 0 ) gtmp->SetFillStyle(3002);
      else  gtmp->SetFillStyle(3006);
      string optDraw = "AE3 LP";
      if( ir9 != 0 ) optDraw = "E3 LP";      
      gtmp->Draw(optDraw.c_str());
      if( ir9 == 0 )    addPlotTitle( gtmp->GetHistogram(), "CMS 2012 RunABCD, L = 19.5 fb^{-1}" );
    }
    { l->Draw(); cOSImprovment->RedrawAxis();  }
   
  }

  ///------------------- polishing plots ---------------------------///

  c[0]->Print( (outputpdf + "[").c_str() );

  for( int ir9 = 0; ir9 < 2; ir9++ ) { 
    c[ir9]->cd(); leg[ir9]->Draw(); 
    c[ir9]->RedrawAxis(); 
    c[ir9]->Print( outputpdf.c_str() ); 
  }
  for( int ir9 = 0; ir9 < 2; ir9++ ) { 
    cOS[ir9]->cd(); legOS->Draw();
    cOS[ir9]->RedrawAxis();
    cOS[ir9]->Print( outputpdf.c_str() ); 
  }
  if( cOSImprovment != 0 ) cOSImprovment->Print( outputpdf.c_str() );

    c[0]->Print( (outputpdf +"]").c_str() );

    // c[0]->Print( "EcalScaleEta_Incl.eps");
  // c[1]->Print( "EcalScaleEta_LowBrem.eps");
  // c[0]->Print( "EcalScaleEta_Incl.C");
  // c[1]->Print( "EcalScaleEta_LowBrem.C");


}



vector<TGraphErrors*> ShervinOS(void) {
  vector<TGraphErrors*> grExtraS_shervin; grExtraS_shervin.resize(2 );
  int dr9[] = {kAzure, kOrange+1 };

  for( int ir9 = 0 ; ir9 < 2 ;ir9++ ) grExtraS_shervin[ir9] = new TGraphErrors(4);
  int  ip = 0; float perCent = 0.01;
  grExtraS_shervin[0]->SetPoint(ip,0.55,0.88*perCent); grExtraS_shervin[0]->SetPointError(ip++,0.1,0.02*perCent);
  grExtraS_shervin[0]->SetPoint(ip,1.30,1.87*perCent); grExtraS_shervin[0]->SetPointError(ip++,0.1,0.02*perCent);
  grExtraS_shervin[0]->SetPoint(ip,1.80,2.32*perCent); grExtraS_shervin[0]->SetPointError(ip++,0.1,0.03*perCent);
  grExtraS_shervin[0]->SetPoint(ip,2.30,2.69*perCent); grExtraS_shervin[0]->SetPointError(ip++,0.1,0.04*perCent);
  grExtraS_shervin[0]->SetMarkerColor( dr9[0] );
  grExtraS_shervin[0]->SetLineColor(   dr9[0] );
  ip = 0;
  grExtraS_shervin[1]->SetPoint(ip,0.55,0.97*perCent); grExtraS_shervin[1]->SetPointError(ip++,0.1,0.02*perCent);
  grExtraS_shervin[1]->SetPoint(ip,1.30,1.31*perCent); grExtraS_shervin[1]->SetPointError(ip++,0.1,0.02*perCent);
  grExtraS_shervin[1]->SetPoint(ip,1.80,2.25*perCent); grExtraS_shervin[1]->SetPointError(ip++,0.1,0.07*perCent);
  grExtraS_shervin[1]->SetPoint(ip,2.30,2.54*perCent); grExtraS_shervin[1]->SetPointError(ip++,0.1,0.03*perCent);
  grExtraS_shervin[1]->SetMarkerColor( dr9[1] );
  grExtraS_shervin[1]->SetLineColor(   dr9[1] );
   return grExtraS_shervin;
}
