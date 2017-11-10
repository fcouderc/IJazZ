#include "interface/EcalUtils.hh"
#include "interface/RootUtils.hh"

#include <TProfile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TMath.h>
#include <TLatex.h>
#include <TString.h>
#include <TFile.h>
#include <TProfile.h>
#include <TF1.h>

#include <string>
#include <iostream>
#include <map>
using namespace std;

harnessEcal harnGeom;
map<int,map<int, int> > nXtalEE[2];


int getXstalPerHarness( int x, int y, int z ) {
  int idee     = harnGeom.deeEE(x,y,z);
  int iharness = harnGeom.harnessEE(x,y,z);
  if( idee < 0 || iharness < 0 ) return 0;
  
  if( nXtalEE[(z+1)/2].find(idee) != nXtalEE[(z+1)/2].end() &&
      nXtalEE[(z+1)/2][idee].find(iharness) != nXtalEE[(z+1)/2][idee].end() ) 
    return nXtalEE[(z+1)/2][idee][iharness];
    
  
  int nXtal = 0;
  for( int ix = 1; ix <= 100; ix++ )
  for( int iy = 1; iy <= 100; iy++ ) {
    int jdee     = harnGeom.deeEE(ix,iy,z);
    int jharness = harnGeom.harnessEE(ix,iy,z);
    if( idee == jdee && iharness == jharness ) nXtal++;
  }

  if( nXtalEE[(z+1)/2].find(idee) != nXtalEE[(z+1)/2].end() ) nXtalEE[(z+1)/2][idee][iharness] = nXtal;
  else {
    map<int,int> tmp; tmp[iharness] = nXtal;
    nXtalEE[(z+1)/2][idee] = tmp;
  }

  return nXtal;
}

void HarnessLevelEnergyScale( string calibFile ) {
  TString outputpdf = "HarnessStudy.pdf";
  gStyle->SetPadRightMargin(0.15);
  harnessEcal harnessViewer;

  TH2F *hEB  = harnessViewer.energyScaleEB( calibFile );
  TH2F *hEEp = harnessViewer.energyScaleEE( calibFile, +1 );
  TH2F *hEEm = harnessViewer.energyScaleEE( calibFile, -1 );

  TCanvas *canEB = new TCanvas("canEB","EB scale",700,600); canEB->cd();
  canEB->Print( outputpdf + "[" );
  hEB->DrawCopy("colz");
  addPlotTitle( hEB, "EB IC per laser harness");
  
  TCanvas *canEEm = new TCanvas("canEEm","EEm scale",700,600); canEEm->cd();
  hEEm->DrawCopy("colz");
  addPlotTitle( hEEm, "EE- IC per laser harness");
  
  TCanvas *canEEp = new TCanvas("canEEp","EEp scale",700,600); canEEp->cd();
  hEEp->DrawCopy("colz");
  addPlotTitle( hEEp, "EE+ IC per laser harness");
  
  canEB ->Print(outputpdf);
  canEEm->Print(outputpdf);
  canEEp->Print(outputpdf);

  vector<TCanvas*> cEB_phi;
  //  cEB_phi->Divide(5,2);
  TLatex tex;
  TH1F *hEB_ICspread = new TH1F("hEB_ICspread","EB spread",100,0.95,1.05);
  float eb_ErrExp = 0;
  for( int ibx = 1; ibx <= hEB->GetXaxis()->GetNbins(); ibx++ ) {
    cEB_phi.push_back( new TCanvas(TString::Format("cEB_iphi%d",ibx),"",600,800) );
    cEB_phi[ibx-1]->Divide(1,2);
    float ieta = hEB->GetXaxis()->GetBinCenter(ibx);
    TH1F *h = (TH1F*)hEB->ProjectionY("toto",ibx,ibx);
    h->SetMaximum(1.02);
    h->SetMinimum(0.98);
    cEB_phi[ibx-1]->cd(1); h->DrawCopy("E");
    addPlotTitle(h, TString::Format("IC #varphi proj. for <i#eta> = %2.1f ",ieta).Data());
    
    float errExp = -1;
    TH1F *hspread = yAxisSpread( h, errExp ); 
    for( int ip = 1; ip <= h->GetNbinsX(); ip++ ) hEB_ICspread->Fill( h->GetBinContent(ip) );
    eb_ErrExp += errExp;
    hspread->SetMaximum( hspread->GetMaximum()*1.20);
    cEB_phi[ibx-1]->cd(2); hspread->DrawCopy("hist"); 
    tex.SetTextSize(0.06);
    tex.DrawLatex(hspread->GetXaxis()->GetXmin(),hspread->GetMaximum(0.95),
		  TString::Format(" <mean> = %1.4f , RMS = %1.4f, pull = %1.1f",
				  hspread->GetMean(),hspread->GetRMS(), hspread->GetRMS()/errExp).Data() );
    cEB_phi[ibx-1]->Print(outputpdf);
  }
  eb_ErrExp /= hEB->GetXaxis()->GetNbins();
  TCanvas* cEB_ICspread = new TCanvas("canEb_ICspread","EB ICspread",600,400); cEB_ICspread->SetLogy();
  hEB_ICspread->SetMaximum(hEB_ICspread->GetMaximum()*1.20);
  hEB_ICspread->GetXaxis()->SetTitle("IC/harness Z#rightarrow ee");
  hEB_ICspread->GetXaxis()->CenterTitle();
  hEB_ICspread->DrawCopy("hist");
  tex.DrawLatex(hEB_ICspread->GetXaxis()->GetXmin(),hEB_ICspread->GetMaximum(0.95),
		  TString::Format(" <mean> = %1.4f , RMS = %1.4f, pull = %1.1f",
				  hEB_ICspread->GetMean(),hEB_ICspread->GetRMS(), hEB_ICspread->GetRMS()/eb_ErrExp).Data() );
  cEB_ICspread->Print(outputpdf);
  
  double bins[] = {-85, -65, -45 ,-25, -5, 0, 5, 25, 45 ,65, 85 };
  TH1F *hICPrecision = new TH1F("hICPrecision","precision vs #eta",10,bins);
  for( int ibx = 1; ibx <= hEB->GetXaxis()->GetNbins(); ibx++ ) {
    float ieta = hEB->GetXaxis()->GetBinCenter(ibx);
    TH1F *h = (TH1F*)hEB->ProjectionY("titi",ibx,ibx);
    float err = 0;   
    for( int ib = 1; ib <= h->GetNbinsX(); ib++ ) err += h->GetBinError(ib);
    err /= h->GetNbinsX();
    if( fabs(ieta) < 5 ) err *= sqrt(100);
    else                 err *= sqrt(200);
    int ibpreci = hICPrecision->FindBin(ieta);
    hICPrecision->SetBinContent(ibpreci,err);
    hICPrecision->SetBinError(ibpreci,0.00001);    
  }
  hICPrecision->GetXaxis()->SetTitle("i#eta_{SC}");
  hICPrecision->GetXaxis()->CenterTitle();
  hICPrecision->GetYaxis()->SetTitle("IC/harness error #times #sqrt{N_{crys/harness}}");
  hICPrecision->SetMarkerStyle(8);
  hICPrecision->SetMarkerSize(1.5);
  hICPrecision->SetMaximum(0.02);
  hICPrecision->SetMinimum(0.00);
  TCanvas *cPreciEB = new TCanvas("cPreciEB","hPrecision",800,500); cPreciEB->cd();
  hICPrecision->DrawCopy("E");
  cPreciEB->Print(outputpdf);

  TProfile *hICPrecisionEEp = new TProfile("hICPrecisionEEp","precision vs #eta",60,10,60);
  TProfile *hICPrecisionEEm = new TProfile("hICPrecisionEEm","precision vs #eta",60,10,60);
  for( int ibx = 1 ; ibx <= hEEp->GetXaxis()->GetNbins(); ibx++ )
  for( int iby = 1 ; iby <= hEEp->GetYaxis()->GetNbins(); iby++ ) {
    int ix = hEEp->GetXaxis()->GetBinCenter(ibx);
    int iy = hEEp->GetYaxis()->GetBinCenter(iby);
    int nXtalp = getXstalPerHarness(ix,iy,+1);
    int nXtalm = getXstalPerHarness(ix,iy,-1);
    //    hEEp->SetBinError(ibx,iby,hEEp->GetBinError(ibx,iby)*sqrt(nXtalp) );
    //    hEEm->SetBinError(ibx,iby,hEEm->GetBinError(ibx,iby)*sqrt(nXtalm) );
    float r = sqrt((ix-50)*(ix-50) + (iy-50)*(iy-50));
    hICPrecisionEEm->Fill(r, hEEp->GetBinError(ibx,iby)*sqrt(nXtalp));
    hICPrecisionEEm->Fill(r, hEEm->GetBinError(ibx,iby)*sqrt(nXtalm));
  }
  // hICPrecisionEEm->SetLineColor(kRed+1);
  // hICPrecisionEEm->SetMarkerColor(kRed+1);
  hICPrecisionEEm->SetMarkerStyle(8);
  hICPrecisionEEm->SetMarkerSize(1.5);
  hICPrecisionEEm->SetMaximum(0.03);
  hICPrecisionEEm->SetMinimum(0.00);
  hICPrecisionEEm->GetXaxis()->SetTitle("R");
  hICPrecisionEEm->GetXaxis()->CenterTitle();
  hICPrecisionEEm->GetYaxis()->SetTitle("IC/harness error #times #sqrt{N_{crys/harness}}");
  hICPrecisionEEp->SetMarkerStyle(8);
  hICPrecisionEEp->SetMarkerSize(1.5);
  TCanvas *cPreciEE = new TCanvas("cPreciEE","hPrecision",800,500); cPreciEE->cd();
  hICPrecisionEEm->DrawCopy("E");
  //  hICPrecisionEEp->DrawCopy("E same");
  cPreciEE->Print(outputpdf);

  canEB->Print(outputpdf +"]");
}


void EE_IC( string calibFile ) {
  TString outputpdf = "EE5x5_IC.pdf";
  
  gStyle->SetPadRightMargin(0.15);
  harnessEcal harnessViewer;

  TH2F *hEEp = harnessViewer.energyScaleEE( calibFile, +1,1 );
  TH2F *hEEm = harnessViewer.energyScaleEE( calibFile, -1,1 );

  TCanvas *canEEm = new TCanvas("canEEm","EEm scale",700,600); canEEm->cd();
  canEEm->Print( outputpdf + "[" );
  hEEm->DrawCopy("colz");
  addPlotTitle( hEEm, "EE- IC per 5x5 crystal tower");
  
  TCanvas *canEEp = new TCanvas("canEEp","EEp scale",700,600); canEEp->cd();
  hEEp->DrawCopy("colz");
  addPlotTitle( hEEp, "EE+ IC per 5x5 crystal tower");
  
  TH1F *hPull[2] = {
    new TH1F("hPullEEp","pullEEm",25,-10,10),
    new TH1F("hPullEEm","pullEEp",25,-10,10) };

  TH1F *hError[2] = {
    new TH1F("hErrorEEp","pullEEm",25,0,5),
    new TH1F("hErrorEEm","pullEEp",25,0,5) };

  IJazZAxisND<double> fitEE = combineEtaScale( calibFile );
  for( int iv = 0 ; iv < fitEE.nBinsND(); iv++ ) {
    int ih = -1;
    if( fitEE.getBinCenterDimN(iv,0) < -85.1 ) ih = 0;
    if( fitEE.getBinCenterDimN(iv,0) > +85.1 ) ih = 1;
    
    if( ih < 0 ) continue;
    if( fabs(fitEE.value(iv) -1 ) < 0.0001 &&  fabs(fitEE.error(iv) -0.001 ) < 0.0001 ) continue;
    hPull[ih]->Fill( (fitEE.value(iv)-1)/ fitEE.error(iv) ); 
    hError[ih]->Fill( fitEE.error(iv) * 100 *5 ); 
  }
  
  TCanvas *cPull = new TCanvas("canPull","pull",900,800);
  cPull->Divide(2,2);
  string pulltitle[2] = {"EE-","EE+"};
  TF1 *f = new TF1("fgaus","gaus",-10,10);
  TLatex tex; tex.SetNDC(); tex.SetTextSize(0.06);
  for( int ih = 0 ; ih < 2; ih++ ) { 
    f->SetParameters(50,0,1);
    cPull->cd(ih+1); 
    hPull[ih]->GetXaxis()->SetTitle("Pull");
    hPull[ih]->GetXaxis()->CenterTitle();
    hPull[ih]->SetMaximum(55);
    hPull[ih]->Fit(f,"MLHE R");
    hPull[ih]->DrawCopy(); 
    tex.DrawLatex(0.20,0.82,TString::Format("pull = %1.2f #pm %1.2f",hPull[ih]->GetRMS(),f->GetParError(2)) );
    addPlotTitle( hPull[ih], "CMS Pull Z IC, " + pulltitle[ih] );

    cPull->cd(ih+3); 
    hError[ih]->SetMaximum(60);
    hError[ih]->GetXaxis()->SetTitle("Error/Crystal [%]");
    hError[ih]->GetXaxis()->CenterTitle();
    hError[ih]->DrawCopy(); 
    addPlotTitle( hError[ih], "CMS IC5x5 Fit uncertainty x 5, " + pulltitle[ih] );
  }

  canEEm->Print(outputpdf);
  canEEp->Print(outputpdf);
  cPull->Print(outputpdf);
  canEEm->Print( outputpdf + "]" );



}
