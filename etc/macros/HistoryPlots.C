#include "interface/EcalUtils.hh"

#include <TCanvas.h>
#include <TH1F.h>
#include <TLine.h>

#include <iostream>
using namespace std;

void ResponseHistoryPlots( string userDir, string version ) {
  string outputpdf[3] = {
    "ResponseHistoryPlotIJazZ_d" + userDir + "_v" + version + "_EEm.pdf",
    "ResponseHistoryPlotIJazZ_d" + userDir + "_v" + version + "_EB.pdf",
    "ResponseHistoryPlotIJazZ_d" + userDir + "_v" + version + "_EEp.pdf"
  };
  IJazZhistory history = historyPlot(userDir, version, "Resp" );
  cout << " # of channels = " << history.chViewer.size() << endl;

  TCanvas *c = new TCanvas("historyCanvas","canvas history",1000,500);
  c->SetRightMargin(0.02);
  c->SetBottomMargin(0.22);

  TLine *l1 = 0;

  c->Print( (outputpdf[0] +"[").c_str() );
  c->Print( (outputpdf[1] +"[").c_str() );
  c->Print( (outputpdf[2] +"[").c_str() );

  for( unsigned ich = 0 ; ich < history.chViewer.size(); ich++ ) {
    TH1F *h = history.chViewer[ich].hist1D(0);
    for( int ib = 1 ; ib <= h->GetXaxis()->GetNbins(); ib++ ) {
      if( h->GetBinContent(ib) <= 0.0000001 ) continue;
      int runmin = int(h->GetXaxis()->GetBinLowEdge(ib));
      int runmax = int(h->GetXaxis()->GetBinUpEdge( ib));
      h->GetXaxis()->SetBinLabel( ib, (itostr(runmin) +"-"+itostr(runmax)).c_str() );
    }
    if( l1 == 0 ) {
      l1 = new TLine( h->GetXaxis()->GetXmin(), 1, h->GetXaxis()->GetXmax(),1);
      l1->SetLineStyle(kDashed);
    }
    string title = "       " + history.channelTitle(ich);
    int ical = 1;
    if( title.find("EE-") != string::npos ) ical = 0;
    if( title.find("EE+") != string::npos ) ical = 2;
    h->SetLabelSize(0.035);
    h->SetMinimum(0.97);
    h->SetMaximum(1.03);
    h->Draw("e"); l1->Draw();
    addPlotTitle(h,title);
    c->Update();
    c->Print( outputpdf[ical].c_str() );
    h->SetMinimum(0.90);
    h->SetMaximum(1.10);
    h->Draw("e"); l1->Draw();
    addPlotTitle(h,title);
    c->Update();
    c->Print( outputpdf[ical].c_str() );


  }
  c->Print( (outputpdf[0] +"]").c_str() );
  c->Print( (outputpdf[1] +"]").c_str() );
  c->Print( (outputpdf[2] +"]").c_str() );
}
