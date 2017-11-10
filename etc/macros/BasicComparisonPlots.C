#include<TLine.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TMath.h>
#include <TStyle.h>
#include <TString.h>
#include <TChain.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TSpline.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TLatex.h>

#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooDataSet.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooDataHist.h>

#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooPolynomial.h>
#include <RooPolyVar.h>
#include <RooCBShape.h>
#include <RooArgusBG.h>
#include <RooBreitWigner.h>
#include <RooNumConvPdf.h>
#include <RooFFTConvPdf.h>
#include <RooChebychev.h>

#include <iostream>
#include <vector>
#include <map>

//#include "Utils.h"

using namespace std;
using namespace RooFit;

float _scaleMC = 1;

TChain * SetData(TString file) {
  TChain *ch = new TChain("selected");
  ch->Add(file);
  return ch;
}

int ipoint_MC = 0;
int ipoint_Data = 0;
int ipoint = -1;
TGraphErrors *grAvsR = 0;
TGraphErrors *grAvsR_Data = new TGraphErrors(6);
TGraphErrors *grAvsR_MC   = new TGraphErrors(6);

float previousAlpha = 1;
bool fixAlpha       = 0;

float Get_Escale_Esigma_FromZ_rootfit( TH1F *hmass, 
				       float &escale, float &e_escale,
				       float &Sigma , float &e_Sigma ,
				       TPad *pad = 0,float meanZ = 91.1876 );


TH1F* CompareVariables( TChain *ch, int color1,  TCut cut, TString v1, TString v2,
		       int nbin = 180, float xmin = 60, float xmax = 150);

TH1F* MakePlot( TChain *chdata, TChain *chmc, int datacolor, TCut cut, TString variable,
		int nbin = 180, float xmin = 60, float xmax = 150,
		TString v1 = "unknown", TString v2 = "unknown", TString vmc = "unknown",
		bool setScale = false); 


void BasicComparisonPlots( TString filedata = "roottuples/FabTreeData_reg_ReCalib0.root",
			   TString filemc   = "roottuples/FabTreeMC_reg_ReCalib0.root"
			   ) {
  TChain *chData = SetData(filedata);
  TChain *chMC   = SetData(filemc);
  chData->SetLineColor( kBlue - 3 );
  chMC->SetLineColor( kRed  - 2 );
  chData->SetMarkerColor( kBlue - 3 );
  chMC->SetMarkerColor( kRed  - 2 );

  bool save = true;
  TCut cut_e1 = "abs(etaSCEle[0])<2.5";
  TCut cut_e2 = "abs(etaSCEle[1])<2.5";
  TCut cut_R9 = "R9Ele[0]>0.94&&R9Ele[1]>0.94";
  TCut cut_notR9 = "R9Ele[0]<0.94&&R9Ele[1]<0.94";

  
  //// basic distri ee
 
  TString corPV = "*1";//91.9083/(91.4413+0.084201*nPV-0.00293408*nPV*nPV)";

  TCut cut_mee = "mee>60 && mee < 130";
  TCut cut_EBEB  = "abs(etaSCEle[0])<1.5&&abs(etaSCEle[1])<1.5";
    TCut cut_M4M4  = "(abs(etaSCEle[0])>1.0&&abs(etaSCEle[0])<1.5&&abs(etaSCEle[1])<1.5)||(abs(etaSCEle[1])>1.0&&abs(etaSCEle[1])<1.5&&abs(etaSCEle[0])<1.5)";
  //  TCut cut_M4M4  = "(abs(etaSCEle[0])>1.0&&abs(etaSCEle[0])<1.5)&&(abs(etaSCEle[1])>1.0&&abs(etaSCEle[1])<1.5)";
  TCut cut_EEEE1 = "(abs(etaSCEle[0])>1.5&&abs(etaSCEle[0])<2.0&&abs(etaSCEle[1])<1.5)||(abs(etaSCEle[1])>1.5&&abs(etaSCEle[1])<2.0&&abs(etaSCEle[0])<1.5)";
  TCut cut_EEEE2 = "(abs(etaSCEle[0])>2.0&&abs(etaSCEle[0])<2.5&&abs(etaSCEle[1])<1.5)||(abs(etaSCEle[1])>2.0&&abs(etaSCEle[1])<2.5&&abs(etaSCEle[0])<1.5)";

  TCut cut_alleta ="";

  int datacolor_noR9 = kBlack;
  int datacolor_R9 = kBlue - 3;
  TCut laser_cut = "";
  vector<TCut>    cut_list;
  vector<TString> cutname_list;
  cut_list.push_back(  cut_alleta  && cut_e1 && cut_e2 ); cutname_list.push_back( "AllEta" );
  cut_list.push_back(  cut_EBEB    && cut_e1 && cut_e2 ); cutname_list.push_back( "EBEB" );
  //  cut_list.push_back(  cut_EBEE   && cut_e1 && cut_e2 ); cutname_list.push_back( "EBEE" );
  cut_list.push_back(  cut_EEEE1   && cut_e1 && cut_e2 ); cutname_list.push_back( "EEEE1" );
  cut_list.push_back(  cut_EEEE2   && cut_e1 && cut_e2 ); cutname_list.push_back( "EEEE2" );
  cut_list.push_back(  cut_M4M4    && cut_e1 && cut_e2 ); cutname_list.push_back( "M4M4" );

  gSystem->Exec(" mkdir -p plots_zee");
  TString basename = "plots_zee/BasicPlots_";
  
  TCut cutKin = cut_list[0];// || cut_list[4];

  if( filemc.Contains("MC") ) {
    TCanvas ctmp; ctmp.cd();
    MakePlot( chData, chMC,  datacolor_noR9, cutKin , "nPV" ,1,0,100,"unknown","unknown","unknown",1);  
  }
  TCanvas *cVars1 = new TCanvas("Zee_Var1","variables",1200,800);

  //// do MC scaling for EBEB pair.
  cVars1->Divide(3,2); 
  cVars1->cd(1); MakePlot( chData, chMC,  datacolor_noR9, cutKin , "ptEle[0]" ,100,   0,100);
  cVars1->cd(2); MakePlot( chData, chMC,  datacolor_noR9, cutKin , "etaEle[0]",100,-2.5,2.5);
  cVars1->cd(4); MakePlot( chData, chMC,  datacolor_noR9, cutKin , "ptEle[1]" ,100,   0,100);
  cVars1->cd(5); MakePlot( chData, chMC,  datacolor_noR9, cutKin , "etaEle[1]",100,-2.5,2.5);
  cVars1->cd(3); MakePlot( chData, chMC,  datacolor_noR9, cutKin , "nPV"      , 36, -0.5,35.5);
  cVars1->cd(6); MakePlot( chData, chMC,  datacolor_noR9, cutKin , "mee_corr" ,100, 70,120);
  if( save )  cVars1->SaveAs( basename + "kin.gif");
  
  
  TLatex text;
  text.SetNDC(); text.SetTextSize(0.07);
  
 for( unsigned icut = 0 ; icut < cut_list.size(); icut++ ) {
   TString canName = "MeeAll_" + cutname_list[icut];

    canName = "Mee2_" + cutname_list[icut];    
    TCanvas *cmee3 = new TCanvas(canName,canName,900,700);
    cmee3->Divide(2,2);

    cmee3->cd(2); CompareVariables( chMC,  kBlue, cut_list[icut]&&cut_notR9   , "mee_corr","unknown", 320,60, 140 );
    fixAlpha = false;
    cmee3->cd(1); CompareVariables( chData,   kRed, cut_list[icut]&&cut_notR9 , "mee_corr","unknown", 320,60, 140 );
    fixAlpha = false;
    cmee3->cd(4); CompareVariables( chMC,  kBlue, cut_list[icut]&&cut_R9      , "mee_corr","unknown", 320,60, 140 );
    fixAlpha = false; 
    cmee3->cd(3); CompareVariables( chData ,   kRed, cut_list[icut]&&cut_R9   , "mee_corr","unknown", 320,60, 140 );
    fixAlpha = false;
    cmee3->Update();

    if( save ) cmee3->SaveAs( basename + canName + ".gif" );

  }
}

float Get_Escale_Esigma_FromZ_rootfit( TH1F *hmass, 
				       float &escale, float &e_escale,
				       float &sigma , float &e_sigma ,
				       TPad *pad, float meanZ) {
  
  escale = 1;
  e_escale = sigma = e_sigma = -1;
  if( hmass->Integral() < 10 ) return -1;
  if( meanZ < 70 ) return -1;
  float mz = 91.1876;
  RooRealVar     M( "M_{ee}", "M_{ee}" ,65. , 115.);
  //  RooRealVar     MeanBW("MeanZ","MeanZ",meanZ,85,94);
  RooRealVar     MeanBW("MeanZ","MeanZ",meanZ,91.188);
  RooRealVar     WidthBW("WidthBW","WidthBW",2.4952);
  RooBreitWigner BW("BW","BW", M, MeanBW, WidthBW);

  RooRealVar     Mean("Mean","Mean",0.0,-1,1);
  RooRealVar     Sigma("Sigma","Sigma",1.,0.1,3.);
  RooRealVar     N("N","N",5);//,0,10);
  // RooRealVar     N("N","N",5,0.1,10);
  RooRealVar     Alpha("Alpha","Alpha",3.0,1.0,10.);
  if( fixAlpha ) {
    Alpha.setVal( previousAlpha );
    Alpha.setConstant();
  }
  RooCBShape     sig("CB","CB",M, Mean ,Sigma, Alpha, N);
  RooFFTConvPdf  Conv1("Conv1","Conv",M,BW,sig);

 // Build Chebychev polynomial p.d.f.  
  RooRealVar a0("a0","a0", -1./M.getMax()) ;
  RooPolynomial bkg("bkg1","Background",M,RooArgSet(a0)) ;

  // full distri (fix Nbkg to 0)
  RooRealVar Nbkg("Nbkg","Nbkg",0) ;
  RooRealVar Nsig("Nsig","Nsig",hmass->Integral()*0.95,0.,1.0*hmass->Integral()*3) ;

  RooAddPdf   model1("model1","bkg+sig",RooArgList(bkg,Conv1),RooArgList(Nbkg,Nsig)) ;
  RooDataHist Mdist("Mdist","Mdist", M, hmass);

  if( pad != 0 ) pad->cd();
  RooPlot* plot=M.frame();
  Mdist.plotOn(plot); 

  model1.fitTo(Mdist,Verbose(0),PrintLevel(1),Warnings(0),PrintEvalErrors(-1),Extended(),Minos(false) );
  model1.plotOn(plot,LineColor(kBlue-3));
  model1.plotOn(plot,LineStyle(kDashed),Components(bkg));
  model1.plotOn(plot,LineColor(kBlue-3),LineStyle(kDashed),Components(bkg));

  previousAlpha = Alpha.getVal();
  //  escale   = MeanBW.getVal()   / mz;
  //  e_escale = MeanBW.getError() / mz;
  //  escale   = (MeanBW.getVal() + Mean.getVal() ) / mz;
  //  e_escale = Mean.getError() / mz;
  //  sigma    = Sigma.getVal()    / mz;
  //  e_sigma  = Sigma.getError()  / mz;
  escale   =  Mean.getVal() ;
  e_escale = Mean.getError() ;
  sigma    = Sigma.getVal()  ;
  e_sigma  = Sigma.getError();
  double fbkg   = Nbkg.getVal() / (Nbkg.getVal()+Nsig.getVal());
  double e_fbkg = Nbkg.getError() / (Nbkg.getVal()+Nsig.getVal());

  plot->Draw();
  if( pad != 0 ) pad->Update();

  return Nsig.getVal();
}


TH1F* CompareVariables( TChain *ch, int color1,  TCut cut, TString v1, TString v2,
		       int nbin, float xmin, float xmax) {
 
  TH1F *hmassD = new TH1F("hmassD","Z->ee mass", nbin,xmin,xmax);
  TH1F *hmassM = new TH1F("hmassM","Z->ee mass", nbin,xmin,xmax);
  hmassD->Sumw2();
  hmassM->Sumw2();


  ch->Draw(v1 + ">>hmassD",cut ,"goff");
  if( v2 != "unknown" ) ch->Draw(v2 + ">>hmassM",cut, "goff");
  hmassM->SetLineColor(kGray + 2);
  hmassM->SetFillColor(kGray + 1);
  hmassD->SetLineColor(color1);
  hmassD->GetXaxis()->SetTitle(v1);
  hmassM->GetXaxis()->SetTitle(v2);
  hmassM->GetXaxis()->CenterTitle();
  hmassD->GetXaxis()->CenterTitle();

  hmassD->SetMarkerStyle(8);
  hmassD->SetMarkerColor(color1);
  hmassD->SetMarkerSize(0.8);
  float max = hmassM->GetMaximum();
  if( hmassD->GetMaximum() > max ) max = hmassD->GetMaximum();
  hmassM->SetMaximum( max *1.10);
  hmassM->SetMinimum( 0);

  if( v2 != "unknown" ) {
    hmassM->DrawCopy("HIST");
    hmassD->DrawCopy("e same");
  } else {
    hmassD->DrawCopy("e");
  }
//   /// fill mass histogram
  TPad * pad = (TPad*) gPad->GetPad(0);

  float mzee_cor(0),s_cor(0);
  float e_cor(0),e_s_cor(0);
  float nsig =  Get_Escale_Esigma_FromZ_rootfit(hmassD, mzee_cor,e_cor,s_cor,e_s_cor,pad);
  
  pad->cd();
  TLatex text; text.SetNDC();
  char str[1000];text.SetTextSize(0.06);
  //  sprintf( str,"M = %2.2f +/- %1.2f GeV", mzee_cor876, e_cor*91.1876 );
  sprintf( str,"M = %2.2f +/- %1.2f GeV", mzee_cor, e_cor*91.1876 );
  text.DrawLatex(0.18,0.8,str);
  float esigma = s_cor/mzee_cor * sqrt( e_s_cor*e_s_cor/s_cor/s_cor+e_cor*e_cor/mzee_cor/mzee_cor);
  //  sprintf( str,"#sigma_{M} / M = %2.2f +/- %2.2f %%", s_cor/mzee_cor*100,  esigma*100);
  sprintf( str,"#sigma_{M} / M = %2.2f +/- %2.2f GeV", s_cor ,  esigma);
  text.DrawLatex(0.18,0.7,str);
  sprintf( str,"N_{Z} = %2.1f", nsig );
  text.DrawLatex(0.18,0.6,str);
  
  return hmassD;
}

TH1F* MakePlot( TChain *chdata, TChain *chmc, int datacolor, TCut cut, TString variable,	       
		int nbin, float xmin, float xmax,
		TString v1, TString v2, TString vmc , bool setScale ) {


  TPad *p0 = (TPad*) gPad->GetPad(0);
  string pUpName = p0->GetName(); pUpName += "_Up"; 
  string pDoName = p0->GetName(); pDoName += "_Do"; 
  double ysep = 0.30;
  TPad *pUp = new TPad( pUpName.c_str(), pUpName.c_str(), 0,ysep,1,1,0,0,0);
  TPad *pDo = new TPad( pDoName.c_str(), pDoName.c_str(), 0,0,1,ysep,0,0,0);
  pUp->SetTopMargin(0.005);  pUp->SetBottomMargin(0.00);  
  pDo->SetTopMargin(0.000);  pDo->SetBottomMargin(0.21);  

  TH1F *hmassD  = new TH1F("hmassD" ,"Z->ee mass", nbin,xmin,xmax);
  TH1F *hmassD1 = new TH1F("hmassD1","Z->ee mass", nbin,xmin,xmax);
  TH1F *hmassD2 = new TH1F("hmassD2","Z->ee mass", nbin,xmin,xmax);
  TH1F *hmassM  = new TH1F("hmassM" ,"Z->ee mass", nbin,xmin,xmax);
  hmassD->Sumw2();
 

  TString mass_str = variable;
  chdata->Draw(mass_str + ">>hmassD",cut ,"goff");
  if( v1 != "unknown" ) chdata->Draw(v1 + ">>hmassD1",cut ,"goff");
  if( v2 != "unknown" ) chdata->Draw(v2 + ">>hmassD2",cut ,"goff");

  // if( mass_str.Contains("R9") ) mass_str = "1.0035*"+mass_str;
  if( vmc == "unknown" ) vmc = variable;
  chmc  ->Draw( vmc + ">>hmassM","weight"*cut, "goff");
  hmassM->SetLineColor(kGray + 2);
  hmassM->SetFillColor(kGray + 1);
  hmassD->SetLineColor(datacolor);
  hmassD1->SetLineColor(kGreen+3);
  hmassD2->SetLineColor(kRed -3);
  hmassD1->SetMarkerColor(kGreen+3);
  hmassD2->SetMarkerColor(kRed -3);
  hmassD1->SetMarkerSize(0.8);
  hmassD2->SetMarkerSize(0.8);

  hmassM->GetXaxis()->SetTitle(variable);
  hmassD->GetXaxis()->SetTitle(variable);
  hmassM->GetXaxis()->CenterTitle();
  hmassD->GetXaxis()->CenterTitle();

  hmassD->SetMarkerStyle(8);
  hmassD->SetMarkerColor(datacolor);
  hmassD->SetMarkerSize(0.8);
  if( setScale && hmassM->Integral()  > 0 ) {
    // int ib1 = hmassD->GetXaxis()->FindBin(75);
    // int ib2 = hmassD->GetXaxis()->FindBin(110);
    int ib1 = hmassD->GetXaxis()->FindBin(85);
    int ib2 = hmassD->GetXaxis()->FindBin(95);
    _scaleMC = hmassD->Integral(ib1,ib2) / hmassM->Integral(ib1,ib2) ;
    //    _scaleMC = 1;
  }
  
  hmassM->Scale(_scaleMC );
  float max = hmassM->GetMaximum();
  if( hmassD->GetMaximum() > max ) max = hmassD->GetMaximum();
  hmassM->SetMaximum( max *1.10);
  hmassM->SetMinimum( 0.0001);
  hmassD->SetMaximum( max *1.10);
  hmassD->SetMinimum( 0.0001);

  pUp->cd();
  hmassM->DrawCopy("HIST");
  if( v2 != "unknown" )  hmassD2->DrawCopy("e same");
  if( v1 != "unknown" )  hmassD1->DrawCopy("e same");
  hmassD->DrawCopy("e same");

  if( ( variable == "PVz" || variable == "nPV" ) && false ) {
    TFile *f = new TFile("calib_files/nPV_calib.root","update");
    hmassM->Write( variable + "_mc"  ,TObject::kOverwrite );
    hmassD->Write( variable + "_data",TObject::kOverwrite );
    cout << " variable: " << variable << " saved in file: " << f->GetName() << endl;;
    f->Close();
    //    return (TH1F*) 0;
  }
  
  hmassD->Add( hmassM,-1);
  // if( v2 != "unknown" )  { hmassD2->Add(hmassM,-1); SetMaximum( hmassD, hmassD2 ); }
  // if( v1 != "unknown" )  { hmassD1->Add(hmassM,-1); SetMaximum( hmassD, hmassD1 ); }
  hmassD->SetLabelSize(0.10,"xy");
  hmassD->SetTitleSize(0.11,"x");
  pDo->cd();
  hmassD->DrawCopy("e");
  if( v1 != "unknown" )  hmassD1->DrawCopy("e same");
  if( v2 != "unknown" )  hmassD2->DrawCopy("e same");


  p0->cd(0);
  pUp->Draw();
  pDo->Draw();
  p0->Update();

  return hmassD;
}



//-------------------------------- Make FZ measurement --------------------------------//
TGraphErrors* FzEnergyDependence( TString file ) {
  TTree *t = (TTree*) TFile::Open( file, "read" )->Get("selected");
  TCut cut0 = "ptEle[0] > 35 && ptEle[1]> 30";
  TCanvas *c = new TCanvas;
  int np = 40;
  TGraphErrors *gr = new TGraphErrors(np);
  for( int it = 0 ; it < np; it++ ) {
    TH1F hmass("hm","mass [GeV]",80,70,110);
    float fzmin = 0.5+it*1./np;
    float fzmax = 0.5+(it+1)*1./np;
    TString cutstr = TString::Format("fZ>=%f && fZ < %f",fzmin,fzmax);
    TCut    cut(cutstr);
    t->Draw("mee_corr>>hm",cut,"goff");
    float escale(1), e_escale(0.01);
    float sigma(0.01), e_sigma(0.001);
    Get_Escale_Esigma_FromZ_rootfit( &hmass, escale, e_escale,
				     sigma , e_sigma , c );
    gr->SetPoint(it,(fzmin+fzmax)/2.,escale);
    gr->SetPointError(it,(fzmax-fzmin)/2.,e_escale);
  }

  gr->Draw("AP");
  return gr;
}

void FzDataMC( TString data, TString mc )  {

  TGraphErrors *fzData = FzEnergyDependence(data);
  TGraphErrors *fzMC   = FzEnergyDependence( mc );

  fzMC  ->SetLineColor(kRed+1);
  fzMC  ->SetMarkerColor(kRed+1);
  fzData->SetLineColor(kAzure+1);
  fzData->SetMarkerColor(kAzure+1);

  fzData->Draw("AP");
  fzMC  ->Draw("P");

}
