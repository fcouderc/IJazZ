#include "interface/IJazZAxis.hh"
#include "interface/IJazZAxisND.hh"
#include "interface/IJazZAxisViewer.hh"

#include "interface/EcalUtils.hh"
#include "interface/RootUtils.hh"

#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TLine.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TSystem.h>

#include <math.h>
#include <iostream>

#include "boost/algorithm/string.hpp"
#include "boost/algorithm/string/split.hpp"
using namespace std;

IJazZhistory historyPlot( string userDir, string version, string ESorER ){
  assert( ESorER == "Resp" || ESorER == "Reso" || ESorER == "Resp.etaScaleDataOverMC" );

  /// I will assume here std IJazZ path building to create the input filelist
  string filelist = "/tmp/ijazz.filelist";
  string dir = "IJazZ_output/FitResults/" + userDir + "/*/*EB*Data*" + version + ".root.fitted" + ESorER + " > " + filelist;
  gSystem->Exec( ("ls -1 "+ dir).c_str() );

  IJazZAxis<double> rRange; rRange.setName("run");
  vector<IJazZAxisND<double> > fitResults;
  vector<double> runMean;
  ifstream listOfIJazZfit(filelist.c_str());
  int rMin = -1; 
  int rMax = -1; 
  while (  listOfIJazZfit.good() && !listOfIJazZfit.eof() ) {
    string line;
    getline( listOfIJazZfit,line,'\n');
    if( line.size() < 1 ) continue;
    if( line.find("IJazZ_output") == string::npos ) continue;
    runRangeChecker runchecker(line);
    rMin = runchecker.rMin();
    rMax = runchecker.rMax();
    if( rMin >= 0 && rMax <= 999999 ) {
      rRange.addBin( runchecker.rMin() );
      rRange.addBin( runchecker.rMax() );
      fitResults.push_back( combineEtaScale( line ) );
      runMean.push_back( double(rMin+rMax)/2.) ;
    }
  }
  
  rRange.print();
  /// now loop over bins
  vector<IJazZAxisND<double> >     channels;
  vector<IJazZAxisViewer<double> > channelsViewer;
  vector<vector<string> > channelsNames;
  vector<vector<float> >  channelsBinCenter;
  for( int ich = 0 ; ich < fitResults[0].nBinsND() ; ich++ ) {
    IJazZAxisND<double> channel; channel.addAxis( rRange );
    int nTimes1 = 0;
    for( unsigned it = 0 ; it < fitResults.size(); it++ ) {
      vector<double> run; run.resize(1,runMean[it]);
      channel.value( channel.findBin( run ) ) = fitResults[it].value(ich);
      channel.error( channel.findBin( run ) ) = fitResults[it].error(ich);
      if( fabs( fitResults[it].value(ich) -1    ) < 0.00001 ||
	  fabs( fitResults[it].value(ich) -0.02 ) < 0.00001 ) nTimes1++;
    }
    channels.push_back(channel);
    if( nTimes1 > 3 ) continue;
    channelsViewer.push_back( IJazZAxisViewer<double>(channels[ich]) );
    vector<string> xnames;
    vector<float>  xval;
    for( unsigned id = 0; id < fitResults[0].getND(); id++ ) {
      xnames.push_back( fitResults[0].getAxis(id)->getName() );
      xval  .push_back( (float) fitResults[0].getBinCenterDimN( ich, id ));
    }
    channelsNames.push_back(xnames);
    channelsBinCenter.push_back(xval);
  }
  IJazZhistory out;
  out.chViewer   = channelsViewer;
  out.chAxisName = channelsNames ;
  out.chAxisVal  = channelsBinCenter ;
  return out;
}

string IJazZhistory::channelTitle( int ich ) {
  string title;
  if( chAxisName[ich].size() == 3 && 
      chAxisName[ich][0] == "IEta" && chAxisName[ich][1] == "IDee" && chAxisName[ich][2] == "IHarness" ) {
    if     ( chAxisVal[ich][0] < -85 ) title += "EE- ,";
    else if( chAxisVal[ich][0] > +85 ) title += "EE+ ,";
    else                                {
      title += "EB  , i#eta: "; title += ftostr( harnessEcal::ietaFromDeeEB(chAxisVal[ich][1]) );
      title += " , i#varphi: "; title += ftostr( harnessEcal::iphiDeeAndHarnessEB(chAxisVal[ich][1],chAxisVal[ich][2]) );
      return title;
    }
    title += " iDee = "    ;  title += ftostr( chAxisVal[ich][1] );
    title += " iHarness = ";  title += ftostr( chAxisVal[ich][2] );
    return title;
  }
  if( chAxisName[ich].size() == 2 && 
      chAxisName[ich][0] == "AbsSCEta" && chAxisName[ich][1] == "R9" ) {
    if     ( chAxisVal[ich][0] < 1.  ) title += "EB-absEta_0_1";
    else if( chAxisVal[ich][0] < 1.5 ) title += "EB-absEta_1_1.4442";
    else if( chAxisVal[ich][0] < 2.0 ) title += "EE-absEta_1.566_2";
    else if( chAxisVal[ich][0] < 2.5 ) title += "EE-absEta_2_2.5";
    title += "-";
    if( chAxisVal[ich][1] < 0.94 ) title += "bad";
    else                           title += "gold";
    return title;
  }


  /// generic name
  for( unsigned id = 0 ; id < chAxisName[ich].size(); id++ ) {
    title += chAxisName[ich][id] + ": ";
    title += ftostr(chAxisVal[ich][id]);
    if( id != chAxisName[ich].size()-1 ) title += " , ";
  }
  return title;
}

///------------------ run range checker ------------------------------///
runRangeChecker::runRangeChecker( string fileCorr ) {
  vector<string> breakFilePath;
  boost::split( breakFilePath, fileCorr, boost::is_any_of("/"), boost::token_compress_on );
  vector<string>::iterator p = find( breakFilePath.begin(), breakFilePath.end(), "IJazZ_output" );
  if( p == breakFilePath.end() ) {
    cerr << "---- IJazZ WARNING: correction file (" << fileCorr << ") path has wrong format, " << endl
	 << "                    can not find run range for this correction..." << endl;
    cerr << "                    I will apply it for every event" << endl;
    _runMin = -1;
    _runMax = 999999;
    return;
  }
  p += 3;
  vector<string> runs;
  boost::split( runs, *p, boost::is_any_of("-"), boost::token_compress_on );
  if( runs.size() != 2 ) {
    cerr << "IJazZ::runRangeChecker(): correction file: " << fileCorr << endl
	 <<" has not the expected path format, i.e.: /.../calib_files/xxx/runMin-runMax/...." << endl;
    cerr << "    ... bail out ... " << endl;
  }
  _runMin = atoi(runs[0].c_str());
  _runMax = atoi(runs[1].c_str());
}



IJazZAxisND<double> combineEtaScale(  string  ecalFilePattern ) {
  int n_cFiles = 3;
  
  string cFiles[3];
  int i1 = ecalFilePattern.find("_E");
  string ecalp = ecalFilePattern.substr(i1+1,2);
  if     ( ecalp == "EB"  ) cFiles[0] = ecalFilePattern.replace(i1+1,2,"EB" );
  else if( ecalp == "EE"  ) cFiles[0] = ecalFilePattern.replace(i1+1,2,"EB" );
  else if( ecalp == "EEm" ) cFiles[0] = ecalFilePattern.replace(i1+1,3,"EB" );
  else if( ecalp == "EEp" ) cFiles[0] = ecalFilePattern.replace(i1+1,3,"EB" );
  else {
    //--- non std format, do not attempt to change anything
    IJazZAxisND<double> out; out.readFromFile(ecalFilePattern);
    return out;
  }
  cFiles[0] = ecalFilePattern;

  IJazZAxisND<double> calibEB; calibEB.readFromFile(cFiles[0]);
  if( calibEB.getAxis(0)->getName().find("Abs") != string::npos ) n_cFiles = 2;
  if( n_cFiles == 3 ) {
    cFiles[1] = ecalFilePattern.replace(i1+1,2,"EEp");
    cFiles[2] = ecalFilePattern.replace(i1+1,3,"EEm");    
  } else {
    cFiles[1] = ecalFilePattern.replace(i1+1,2,"EE");
    cFiles[2] = "unknown";    
  }
  
  IJazZAxisND<double> calib[3]; 
  bool open[3]; for( int i = 0 ; i < n_cFiles; i++ ) open[i] = calib[i].readFromFile(cFiles[i]);
  IJazZAxisND<double> out = calibEB;
  double CCstop = 1.502;
  if( calibEB.getAxis(0)->getName().find("IEta") != string::npos ) CCstop = 85.2;

  for( int ib = 0 ; ib < out.nBinsND(); ib++ ) {
    int icalib = -1;
    if     ( out.getBinCenterDimN(ib,0) < -CCstop ) icalib = 2;
    else if( out.getBinCenterDimN(ib,0) > +CCstop ) icalib = 1;
    else icalib = 0;
    if( !open[icalib] ) continue;
    if( icalib == 2 && n_cFiles == 2 ) {
      cerr << " This can not happen: eta = " << out.getBinCenterDimN(ib,0) << " is < 0 while axis name contains Abs" << endl;
    }
    out.value(ib) = calib[icalib].value(ib);
    out.error(ib) = calib[icalib].error(ib);
  }
  return out;
}

///// helper for final plots
vector<TLine*> ecalModuleEdges(float ymin, float ymax, bool IETA ) {
  vector<float> xModEdge;
  xModEdge.push_back(-85.5);
  xModEdge.push_back(-65.5);
  xModEdge.push_back(-45.5);
  xModEdge.push_back(-25.5);
  xModEdge.push_back(  0.0);
  xModEdge.push_back(+25.5);
  xModEdge.push_back(+45.5);
  xModEdge.push_back(+65.5);
  xModEdge.push_back(+85.5);

  vector<float> xModEdgeEta;
  xModEdgeEta.push_back(-1.50);
  xModEdgeEta.push_back(-1.137);
  xModEdgeEta.push_back(-0.789);
  xModEdgeEta.push_back(-0.440);
  xModEdgeEta.push_back(  0.0);
  xModEdgeEta.push_back(+0.440);
  xModEdgeEta.push_back(+0.789);
  xModEdgeEta.push_back(+1.137);
  xModEdgeEta.push_back(+1.50);


  vector<TLine*> edges; edges.resize(xModEdge.size());
  for( unsigned iedge = 0 ; iedge < xModEdge.size(); iedge++ ) {
    if( IETA ) edges[iedge] = new TLine( xModEdge[iedge]   , ymin, xModEdge[iedge]   , ymax );
    else       edges[iedge] = new TLine( xModEdgeEta[iedge], ymin, xModEdgeEta[iedge], ymax );
    edges[iedge]->SetLineWidth(2);
    edges[iedge]->SetLineStyle(kDashed);
  }
  return edges;
}


void addPlotEmbellishement( TH1F *frame, string plotTitle, bool addLineOne  ) {
  float xmin = frame->GetXaxis()->GetBinLowEdge(frame->GetXaxis()->GetFirst());
  float xmax = frame->GetXaxis()->GetBinUpEdge( frame->GetXaxis()->GetLast() );

  float ymax = frame->GetMaximum(); 
  vector<TLine*> modEdges = ecalModuleEdges(frame->GetMinimum(),frame->GetMaximum() );
  TLatex tex;
  tex.SetTextColor(kBlack); tex.SetTextSize(0.05); 
  tex.DrawLatex(xmin*0.998,(ymax+0.001)*1.002, plotTitle.c_str() );

  for( unsigned ie = 0 ;ie < modEdges.size(); ++ie ) modEdges[ie]->Draw();

  if( addLineOne ) {
    TLine *l1 = new TLine( xmin,1, xmax, 1 );
    l1->SetLineStyle( kDashed );
    l1->Draw();
    TLine *l2 = new TLine( xmin,0, xmax, 0 );
    l2->Draw();
  }
}



#include <TFile.h>
void plotIOVFile( string iofile ) {
  TH2F* hEcal[3] = { new TH2F( "EEm","EEp Energy scale", 100, 0.5, 100.5, 100, 0.5, 100.5 ),
		     new TH2F( "EB" ,"EB  Energy scale", 171, -85.5, 85.5, 360, 0.5, 360.5 ),
		     new TH2F( "EEp","EEm Energy scale", 100, 0.5, 100.5, 100, 0.5, 100.5  ) };

  ifstream skeleton(iofile.c_str());
  while ( skeleton.good() && !skeleton.eof() ) {
    string line;
    getline(skeleton,line,'\n');
    istringstream isstream(line);
    int ix,iy,iz;
    float ic, eic;
    isstream >> ix >> iy >> iz >> ic >> eic;
    int ibx =  hEcal[iz+1]->GetXaxis()->FindBin(ix);
    int iby =  hEcal[iz+1]->GetYaxis()->FindBin(iy);
    hEcal[iz+1]->SetBinContent( ibx,iby,ic);
    hEcal[iz+1]->SetBinError  ( ibx,iby,eic);
  }
  hEcal[0]->GetXaxis()->SetTitle("ix");  hEcal[0]->GetYaxis()->SetTitle("iy");
  hEcal[2]->GetXaxis()->SetTitle("ix");  hEcal[2]->GetYaxis()->SetTitle("iy");
  hEcal[1]->GetXaxis()->SetTitle("i#eta");  hEcal[1]->GetYaxis()->SetTitle("i#varphi");
  for( int i = 0 ; i < 3; i++ ) {
    hEcal[i]->SetMinimum(0.97);
    hEcal[i]->SetMaximum(1.03);
    hEcal[i]->GetXaxis()->CenterTitle();
    hEcal[i]->GetYaxis()->CenterTitle();
  }

  TCanvas *cEB = new TCanvas("cEB","EB",800,600); cEB->cd();
  cEB->SetRightMargin(0.15);
  hEcal[1]->DrawCopy("colz"); addPlotTitle(hEcal[1],"EB, ICs");
  TCanvas *cEE = new TCanvas("cEE","EE",1100,500); 
  cEE->Divide(2,1);
  for( int ipad = 0 ; ipad < 2; ipad++ ) cEE->GetPad(ipad+1)->SetRightMargin(0.15);
  cEE->cd(1); hEcal[0]->DrawCopy("colz"); addPlotTitle(hEcal[0],"EE-, ICs");
  cEE->cd(2); hEcal[2]->DrawCopy("colz"); addPlotTitle(hEcal[2],"EE+, ICs");
  
  TFile f( (iofile + ".dump.root").c_str(),"recreate");
  f.cd();
  for( int iecal = 0 ; iecal < 3; iecal++ ) hEcal[iecal]->Write();
  f.Close();
}

//---------------------------------------- FWHM --------------------------------------//
#include <TFile.h>
#include <TVectorD.h>
IJazZAxisND<double> meanResoPdf( string rootFileName ) {  
  IJazZAxisND<double> m; m.readFromFile( rootFileName + ".fittedResp" );
  IJazZAxisND<double> s; s.readFromFile( rootFileName + ".fittedReso" );
  IJazZAxisND<double> d; d.readFromFile( rootFileName + ".fittedTailP0" );
  IJazZAxisND<double> r; r.readFromFile( rootFileName + ".fittedTailP1" );
  IJazZAxisND<double> f; f.readFromFile( rootFileName + ".fittedTailP2" );
  
  IJazZAxisND<double> out = m;
  for( int ib = 0 ; ib < out.nBinsND(); ib++ ) {
    vector<double> x( out.getND() );
    vector<double> x0( 1 );
    for( unsigned id = 0; id < out.getND(); id++ ) x[id] = out.getBinCenterDimN( ib, id );
    x0[0] = x[0];
    float param[5];
    if( m.getND() == 1 ) param[0] = m.value(x0);
    else                 param[0] = m.value(x );
    if( s.getND() == 1 ) param[1] = s.value(x0);
    else                 param[1] = s.value(x );
    if( d.getND() == 1 ) param[2] = d.value(x0);
    else                 param[2] = d.value(x );
    if( r.getND() == 1 ) param[3] = r.value(x0);
    else                 param[3] = r.value(x );
    if( f.getND() == 1 ) param[4] = f.value(x0);
    else                 param[4] = f.value(x );

    //    double mean = param[4] * param[0] + ( 1 - param[4] ) * (  param[0] - param[2]*param[4]*param[1] );
    double mean = param[0] - ( 1 - param[4] ) * param[2]*param[1]/((1-param[4])*5) ;
	
    out.value(x) = mean;
    out.error(x) = m.error(x);
  }
  //  fin->Close();
  return out;
}

IJazZAxisND<double> rmsResoPdf( string rootFileName ) {  
  IJazZAxisND<double> m; m.readFromFile( rootFileName + ".fittedResp" );
  IJazZAxisND<double> s; s.readFromFile( rootFileName + ".fittedReso" );
  IJazZAxisND<double> d; d.readFromFile( rootFileName + ".fittedTailP0" );
  IJazZAxisND<double> r; r.readFromFile( rootFileName + ".fittedTailP1" );
  IJazZAxisND<double> f; f.readFromFile( rootFileName + ".fittedTailP2" );
  
  IJazZAxisND<double> out = m;
  for( int ib = 0 ; ib < out.nBinsND(); ib++ ) {
    vector<double> x( out.getND() );
    vector<double> x0( 1 );
    for( unsigned id = 0; id < out.getND(); id++ ) x[id] = out.getBinCenterDimN( ib, id );
    x0[0] = x[0];
    float param[5];
    if( m.getND() == 1 ) param[0] = m.value(x0);
    else                 param[0] = m.value(x );
    if( s.getND() == 1 ) param[1] = s.value(x0);
    else                 param[1] = s.value(x );
    if( d.getND() == 1 ) param[2] = d.value(x0);
    else                 param[2] = d.value(x );
    if( r.getND() == 1 ) param[3] = r.value(x0);
    else                 param[3] = r.value(x );
    if( f.getND() == 1 ) param[4] = f.value(x0);
    else                 param[4] = f.value(x );

    //    double mean = param[4] * param[0] + ( 1 - param[4] ) * (  param[0] - param[2]*param[4]*param[1] );
    double f1   = param[4];
    double m1   = param[0];
    double sig1 = param[1];
    double m2   = m1 - (1-f1)*sig1*f1*param[2];
    double sig2 = param[3]*sig1/(2.*sqrt(2.*log(2.)));
    double rms = sqrt(f1*sig1*sig1+(1-f1)*sig2*sig2+f1*(1-f1)*(m1-m2)*(m1-m2));
	
    out.value(x) = rms;
    out.error(x) = s.error(x)/f1;
  }
  //  fin->Close();
  return out;
}

IJazZAxisND<double> maxResoPdf( string rootFileName ) {  
  IJazZAxisND<double> m; m.readFromFile( rootFileName + ".fittedResp" );
  IJazZAxisND<double> s; s.readFromFile( rootFileName + ".fittedReso" );
  IJazZAxisND<double> d; d.readFromFile( rootFileName + ".fittedTailP0" );
  IJazZAxisND<double> r; r.readFromFile( rootFileName + ".fittedTailP1" );
  IJazZAxisND<double> f; f.readFromFile( rootFileName + ".fittedTailP2" );
  

  string g1x  = "TMath::Gaus(x,[0],[1]*[4],1)*[4]";
  string g2x  = "TMath::BreitWigner(x,[0]-[2]*[1]*[4],[1]*[3])*(1-[4])";
  string formulax = "(" + g1x + "+" + g2x + ")";
  TF1  ftmp("gtmp",formulax.c_str(),0.90,1.10);
  TH1F htmp("htmp","tmp",1,-1,1);
  htmp.SetBinContent(1,0);
  htmp.SetBinError(1,0.0001);

  IJazZAxisND<double> out = s;
  // TCanvas *c = new TCanvas;
  for( int ib = 0 ; ib < out.nBinsND(); ib++ ) {
    vector<double> x( out.getND() );
    vector<double> x0( 1 );
    for( unsigned id = 0; id < out.getND(); id++ ) x[id] = out.getBinCenterDimN( ib, id );
    x0[0] = x[0];
    float param[5];
    if( m.getND() == 1 ) param[0] = m.value(x0);
    else                 param[0] = m.value(x );
    if( s.getND() == 1 ) param[1] = s.value(x0);
    else                 param[1] = s.value(x );
    if( d.getND() == 1 ) param[2] = d.value(x0);
    else                 param[2] = d.value(x );
    if( r.getND() == 1 ) param[3] = r.value(x0);
    else                 param[3] = r.value(x );
    if( f.getND() == 1 ) param[4] = f.value(x0);
    else                 param[4] = f.value(x );

    for( int ip = 0 ; ip <= 4; ip++ ) ftmp.FixParameter( ip, param[ip] );
    double xMax = ftmp.GetMaximumX(0.90,1.10); 
    out.value(x) = xMax;
    out.error(x) = m.error(x);
    //    ftmp.Draw();
    //   c->Update();
  }
  return out;
}
IJazZAxisND<double> widthResoPdf( string rootFileName ) {  
  IJazZAxisND<double> m; m.readFromFile( rootFileName + ".fittedResp" );
  IJazZAxisND<double> s; s.readFromFile( rootFileName + ".fittedReso" );
  IJazZAxisND<double> d; d.readFromFile( rootFileName + ".fittedTailP0" );
  IJazZAxisND<double> r; r.readFromFile( rootFileName + ".fittedTailP1" );
  IJazZAxisND<double> f; f.readFromFile( rootFileName + ".fittedTailP2" );
  

  string g1x  = "TMath::Gaus(x,[0],[1]*[4],1)*[4]";
  string g2x  = "TMath::BreitWigner(x,[0]-[2]*[1]*[4],[1]*[3])*(1-[4])";
  string formulax = "(" + g1x + "+" + g2x + ")";
  TF1  ftmp("gtmp",formulax.c_str(),0.80,1.20);
  TF1  bwtmp("bwtmp",g2x.c_str(),0.80,1.20);

  IJazZAxisND<double> out = s;
  //  TCanvas *c = 0;
  for( int ib = 0 ; ib < out.nBinsND(); ib++ ) {
    vector<double> x( out.getND() );
    vector<double> x0( 1 );
    for( unsigned id = 0; id < out.getND(); id++ ) x[id] = out.getBinCenterDimN( ib, id );
    x0[0] = x[0];
    float param[5];
    if( m.getND() == 1 ) param[0] = m.value(x0);
    else                 param[0] = m.value(x );
    if( s.getND() == 1 ) param[1] = s.value(x0);
    else                 param[1] = s.value(x );
    if( d.getND() == 1 ) param[2] = d.value(x0);
    else                 param[2] = d.value(x );
    if( r.getND() == 1 ) param[3] = r.value(x0);
    else                 param[3] = r.value(x );
    if( f.getND() == 1 ) param[4] = f.value(x0);
    else                 param[4] = f.value(x );

    for( int ip = 0 ; ip <= 4; ip++ ) ftmp.FixParameter( ip, param[ip] );
    for( int ip = 0 ; ip <= 4; ip++ ) bwtmp.FixParameter( ip, param[ip] );
    double yMax = ftmp.GetMaximum(0.90,1.10); 
    if( yMax > 0.0001 ) out.value(x) = 1./(sqrt(TMath::TwoPi())*yMax);
    else                out.value(x) = 0;
    double err = 0;
    if( s.value(x) > 0.001 && f.value(x) > 0.001 ) {
      err = s.error(x)*s.error(x)/(s.value(x)*s.value(x)) + f.error(x)*f.error(x)/(f.value(x)*f.value(x));
      err = s.value(x)/f.value(x) * sqrt(err);
    }
    out.error(x) = err;
    //    if( c == 0 ) c = new TCanvas;
    // ftmp.Draw();
    // bwtmp.SetLineStyle(kDashed);
    // bwtmp.Draw("same");
    // c->Update();
    // getchar();
  }
  return out;
}


IJazZAxisND<double> fwhmOver2( string rootFileName ) {  
  // TFile *fin = TFile::Open("etc/geom/IJazZFit.RotateTail.root","read");
  // assert( fin->IsOpen() );
  // TVectorD *rotEB = (TVectorD*) fin->Get("RotateTail0_ijazz_EB"); 
  // TVectorD *rotEE = (TVectorD*) fin->Get("RotateTail0_ijazz_EE"); 
  // TVectorD *rot;

  IJazZAxisND<double> m; m.readFromFile( rootFileName + ".fittedResp" );
  IJazZAxisND<double> s; s.readFromFile( rootFileName + ".fittedReso" );
  IJazZAxisND<double> d; d.readFromFile( rootFileName + ".fittedTailP0" );
  IJazZAxisND<double> r; r.readFromFile( rootFileName + ".fittedTailP1" );
  IJazZAxisND<double> f; f.readFromFile( rootFileName + ".fittedTailP2" );
  
  string g1x  = "TMath::Gaus(x,[0],[1],1)*[4]";
  string g2x  = "TMath::BreitWigner(x,[0]-[2]*[1],[1]*[3])*(1-[4])";
  string g1   = "TMath::Gaus([5],[0],[1],1)*[4]";
  string g2   = "TMath::BreitWigner([5],[0]-[2]*[1],[1]*[3])*(1-[4])";
  string max  = "(TMath::Gaus([6],[0],[1],1)*[4]";
  max += " + TMath::BreitWigner([6],[0]-[2]*[1],[1]*[3])*(1-[4]))*0.5";
  string formula = "abs(" + g1 + "+" + g2 + "-" + max + ")";
  string formula2 = "abs(" + g1x + "+" + g2x + "-" + max + ")";
  cout << formula << endl;
  TF1  ftmp("ftmp",formula.c_str(),-1,1);
  TF1  gtmp("gtmp",(g1x+"+"+g2x).c_str(),0.90,1.10);
  TF1  gttp("gttp",formula2.c_str(),0.90,1.10);
  TH1F htmp("htmp","tmp",1,-1,1);
  htmp.SetBinContent(1,0);
  htmp.SetBinError(1,0.0001);

  IJazZAxisND<double> out = s;
  // TCanvas *c = new TCanvas;
  for( int ib = 0 ; ib < out.nBinsND(); ib++ ) {
    vector<double> x( out.getND() );
    vector<double> x0( 1 );
    for( unsigned id = 0; id < out.getND(); id++ ) x[id] = out.getBinCenterDimN( ib, id );
    x0[0] = x[0];
    float param[5];
    if( m.getND() == 1 ) param[0] = m.value(x0);
    else                 param[0] = m.value(x );
    if( s.getND() == 1 ) param[1] = s.value(x0);
    else                 param[1] = s.value(x );
    if( d.getND() == 1 ) param[2] = d.value(x0);
    else                 param[2] = d.value(x );
    if( r.getND() == 1 ) param[3] = r.value(x0);
    else                 param[3] = r.value(x );
    if( f.getND() == 1 ) param[4] = f.value(x0);
    else                 param[4] = f.value(x );
    // ---- need to rotate parameter 2
    param[2] = param[2]*param[4];
    for( int ip = 0 ; ip <= 4; ip++ ) ftmp.FixParameter( ip, param[ip] );
    for( int ip = 0 ; ip <= 4; ip++ ) gtmp.SetParameter( ip, param[ip] );
    for( int ip = 0 ; ip <= 4; ip++ ) gttp.SetParameter( ip, param[ip] );

    //    param[0] = (param[0]-1)/0.01;
    //    param[1] = param[1]/0.01;
    // double tail0 = 0;
    // if( fabs(x[0]) < 1.5 ) rot = rotEB;
    // else                   rot = rotEE;
    // rot->Print();
    // for( int i = 0 ; i < rot->GetNrows(); i++ ) tail0 += (*rot)(i) *  param[i];
    // param[2] = tail0*0.01;

    //    cout << " Delta MU = " << tail0 << endl;
    // int ip = 2;
    // ftmp.FixParameter( ip, param[ip] );
    // gtmp.FixParameter( ip, param[ip] );
    // gttp.FixParameter( ip, param[ip] );
    float xMax = gtmp.GetMaximumX(0.9,1.1);
    ftmp.FixParameter(6,xMax);
    double x0_p, x0_m; 
    ftmp.SetParameter(5,1. + s.value(x));  htmp.Fit( &ftmp, "MHE q 0n" ); x0_p = ftmp.GetParameter(5);
    ftmp.SetParameter(5,1. - s.value(x));  htmp.Fit( &ftmp, "MHE q 0n" ); x0_m = ftmp.GetParameter(5);
    //    cout << " x0+ = " << x0_p << " x0- = " << x0_m << " -- s = " << (x0_p-x0_m)/2. << endl;
    double fwhmOver2 = (x0_p-x0_m)/(2.*sqrt(2.*log(2.)));
    out.value(x) = fwhmOver2;
    out.error(x) = s.error(x)/param[4];
    // gtmp.Draw();
    // gttp.Draw("same");
    // c->Update();
    // getchar();
  }
  //  fin->Close();
  return out;
}




//------------------------------------- harness plotted -------------------------------//
harnessEcal::harnessEcal() {
  _defEtaRingEE = "etc/geom/eering.dat";
  _defHarnessEE = "etc/geom/Numbering_EE_All.dat";

  /// EE+ and EE-
  for( int i1 = 0 ; i1 < 150 ; i1++ ) 
  for( int i2 = 0 ; i2 < 150 ; i2++ ) {
    _eering_harness_EEm[i1][i2] = _eering_dee_EEm[i1][i2] = -1;
    _eering_harness_EEp[i1][i2] = _eering_dee_EEp[i1][i2] = -1;
    _eering_eta[i1][i2] = -1;

  }

  ifstream eeDefinition(_defHarnessEE.c_str() );
  int iline = 0 ;
  while ( eeDefinition.good() && !eeDefinition.eof() ) {
    string line;
    getline(eeDefinition,line,'\n');
    istringstream isstream(line);
    int i1,i2,i3,i4,i5,i6,i7,i8,i9;
    isstream >> i1 >> i2 >> i3 >> i4 >> i5 >> i6 >> i7 >> i8 >> i9;
    if( iline != 0 ) {
      if( i3 < 3 ) {
	_eering_harness_EEp[i1][i2] = i7;
	_eering_dee_EEp[i1][i2]     = i3;
      } else {
	_eering_harness_EEm[i1][i2] = i7;
	_eering_dee_EEm[i1][i2]     = i3;
      }
    } else iline++;
  }


  ifstream eenumdata (_defEtaRingEE.c_str());
  while ( eenumdata.good() && !eenumdata.eof() ) {
    string line;
    getline(eenumdata,line,'\n');
    istringstream isstream(line);
    int i1,i2,i3;
    isstream >> i1 >> i2 >> i3;
    _eering_eta[i1][i2] = i3+1;
  }

}

TH2F*  harnessEcal::energyScaleEB( std::string calibFile ) {
  IJazZAxisND<double> histo = combineEtaScale( calibFile) ;
  float yMinEB = 0.98;
  float yMaxEB = 1.02;

  double ebEtaBin[11] = { -85, -65, -45, -25, -5, 0, 5, 25, 45, 65, 85 };
  double ebPhiBin[37]; for( int iphi = 0; iphi < 37; iphi++ ) ebPhiBin[iphi] = iphi*10;
  TH2F *hEB = new TH2F( ("harnessEB"+calibFile).c_str(),"Harness Energy Scale", 10, ebEtaBin, 36, ebPhiBin );
  
  for( int ieta = 1; ieta <= hEB->GetXaxis()->GetNbins(); ieta++ ) 
  for( int iphi = 1; iphi <= hEB->GetYaxis()->GetNbins(); iphi++ ) {
    vector<double> x; x.resize(3);
    x[0] = 0;
    float ietaf = hEB->GetXaxis()->GetBinCenter(ieta);
    float iphif = hEB->GetYaxis()->GetBinCenter(iphi);
    x[1] = deeEB(ietaf);
    x[2] = harnessEB(ietaf,iphif);
    if( x[1] <= 0.5 ) continue; 
    hEB->SetBinContent(ieta,iphi,histo.value(x));
    hEB->SetBinError(  ieta,iphi,histo.error(x));
  }
  addAxisTitles( (TH1F*)hEB, "i#eta", "i#varphi", yMinEB,yMaxEB );
  return hEB;
}

TH2F* harnessEcal::energyScaleEE( std::string calibFile, int plusOrMinus, bool ixiy ) {
  IJazZAxisND<double> histo = combineEtaScale( calibFile);
  float yMinEE = 0.96;
  float yMaxEE = 1.04;
  string hname = "hEEm"; if( plusOrMinus > 0 ) hname = "hEEp";

  TH2F *hEE = new TH2F( (hname+calibFile).c_str(),"Harness Energy scale", 100, 0.5, 100.5, 100, 0.5, 100.5 );
  for( int ix = 1; ix <= hEE->GetXaxis()->GetNbins(); ix++ ) 
  for( int iy = 1; iy <= hEE->GetYaxis()->GetNbins(); iy++ ) {
    vector<double> x; x.resize(3);
    x[0] = plusOrMinus*90;
    int ixf = int( hEE->GetXaxis()->GetBinCenter(ix));
    int iyf = int( hEE->GetYaxis()->GetBinCenter(iy));
    if( plusOrMinus < 0 ) {
      x[1] = _eering_dee_EEm[    ixf][iyf];
      x[2] = _eering_harness_EEm[ixf][iyf];
    } else {
      x[1] = _eering_dee_EEp[    ixf][iyf];
      x[2] = _eering_harness_EEp[ixf][iyf];
    }
    if( ixiy ) { x[1] = ixf; x[2] = iyf; }

    if( x[1] <= 0 || x[2] <= 0 )  continue;
    if( _eering_dee_EEp[ixf][iyf] < 0 || _eering_harness_EEp[ixf][iyf] < 0 ) continue;
    /// remove not fitted point;
    if( fabs(histo.value(x) -1 ) < 0.0001 &&  fabs(histo.error(x) -0.001 ) < 0.0001 ) continue;
    hEE->SetBinContent(ix,iy,histo.value(x));
    hEE->SetBinError(  ix,iy,histo.error(x));    
  }
  addAxisTitles( (TH1F*)hEE, "ix", "iy",  yMinEE,yMaxEE );
  return hEE;
}

int harnessEcal::deeEE(     int ix, int iy, int iside  ) { return iside > 0 ? _eering_dee_EEp[ix][iy]:_eering_dee_EEm[ix][iy]; }
int harnessEcal::harnessEE( int ix, int iy, int iside  ) { return iside > 0 ? _eering_harness_EEp[ix][iy]:_eering_harness_EEm[ix][iy]; }
int harnessEcal::deeEB(     float ieta ) {   
  if     ( ieta > -85.5 && ieta < -65.5 ) return  1;
  else if( ieta > -65.5 && ieta < -45.5 ) return  2;
  else if( ieta > -45.5 && ieta < -25.5 ) return  3;
  else if( ieta > -25.5 && ieta <  -5.5 ) return  4;
  else if( ieta > - 5.5 && ieta < - 0.5 ) return  5;
  else if( ieta > + 0.5 && ieta < + 5.5 ) return  6;
  else if( ieta > + 5.5 && ieta < +25.5 ) return  7;
  else if( ieta > +25.5 && ieta < +45.5 ) return  8;
  else if( ieta > +45.5 && ieta < +65.5 ) return  9;
  else if( ieta > +65.5 && ieta < +85.5 ) return 10;
  
  return -1;
}

int harnessEcal::harnessEB( float ieta, float iphi ) { return fabs(ieta)>5.001? TMath::Ceil(iphi/10. ):TMath::Ceil(iphi/20. ) ; }

int harnessEcal::ietaEE( float ix, float iy, int side ) { 
  return side*(_eering_eta[int(ix)][int(iy)] + NETARING_EB); 

}


int  harnessEcal::ieta( float scEta, float ix, float iy ) {
  if( fabs( scEta ) > 1.5001 ) {
    if( !( int(ix) >= 0 && int(ix) <= 149 && int(iy) >= 0 && int(iy) <= 149 ) ) {
      cout << "Something wrong in etaEE calculation: ix = " << ix << " - iy = " << iy << " while scEta = " << scEta <<  endl;
      return 120;
    }
  }
  
  if     ( scEta < -1.50001 ) return ietaEE(ix,iy,-1);
  else if( scEta > +1.50001 ) return ietaEE(ix,iy,+1);
  return ix;
}

int  harnessEcal::dee( float scEta, float ix, float iy ) {
  if     ( scEta < -1.50001 ) return deeEE(ix,iy,-1);
  else if( scEta > +1.50001 ) return deeEE(ix,iy,+1);
  return deeEB(ix);
}

int  harnessEcal::harness( float scEta, float ix, float iy ) {
  if     ( scEta < -1.50001 ) return harnessEE(ix,iy,-1);
  else if( scEta > +1.50001 ) return harnessEE(ix,iy,+1);
  return harnessEB(ix,iy);
}

int  harnessEcal::ix( float scEta, float ix ) {
  if     ( scEta < -1.50001 || scEta > +1.50001 ) return ix;
  return 1;
}

int  harnessEcal::iy( float scEta, float iy ) {
  if     ( scEta < -1.50001 || scEta > +1.50001 ) return iy;
  return 1;
}


float harnessEcal::ietaFromDeeEB( float dee ) {
  if( dee == 1 ) return -75;
  if( dee == 2 ) return -55;
  if( dee == 3 ) return -35;
  if( dee == 4 ) return -15;
  if( dee == 5 ) return -2.5;
  if( dee == 6 ) return +2.5;
  if( dee == 7 ) return +15;
  if( dee == 8 ) return +35;
  if( dee == 9 ) return +55;
  if( dee ==10 ) return +75;

  return -1;
}


float harnessEcal::iphiDeeAndHarnessEB( float dee, float harness ) {
  if( dee == 5 || dee == 6 ) return harness*20-10;
  return harness*10 - 5;
}




float convert_ieta_eta( float ieta ) {
  double pEB[] = {-4.75731e-03,5.73445e+01};
  /// --> [0] + (x/[1])^2
  //double pEE[] = {4.11086e-01,8.33415e+01};
  /// next one give better result: --> [0] + (x/[1])^4
  double pEE[] = {1.13288e+00,1.09372e+02};

  float sign = +1;
  if( ieta < 0 ) sign = -1;
  ieta = fabs(ieta);
  if( ieta < 85 ) return sign*(pEB[0]+ieta/pEB[1]);
  float r =  ieta/pEE[1];
  //return sign*(pEE[0] + r*r);
  return sign*(pEE[0] + r*r*r*r);
}


float convert_eta_ieta( float eta ) {
  double pEB[] = {-4.75731e-03,5.73445e+01};
  double pEE[] = {7.36055e-01,1.89612e+02,5.37159e+01};

  float sign = +1;
  if( eta < 0 ) sign = -1;
  eta = fabs(eta);
  if( eta < 1.46 ) return sign*((eta-pEB[0])*pEB[1]);
  return sign*( sqrt(eta-pEE[0])*pEE[1] - pEE[2]*eta);
}
