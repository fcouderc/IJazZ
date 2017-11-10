//---------------------------------------------------------------------------//
// Description: a few utils functions for IJazZ:
//              - constructor: default value 
//              - vertex reweighting handler
//              - graphical interface
//
//
// Author: Fabrice Couderc, fabrice.couderc@cea.fr 
// History: 2012/11/28 - v1
//
//---------------------------------------------------------------------------//

#include "interface/IJazZ.hh"
#include "interface/RootUtils.hh"
#include "interface/EcalUtils.hh"

#include <TChain.h>
//#include <TFitterMinuit.h>
#include <TH1.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <TLine.h>
#include <TFile.h>
#include <TLatex.h>
#include <TVectorD.h>

#include <map>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;


//--- output file name handling
string IJazZ::createOutputFileName( string ecalp ) {
  _ecalp = ecalp;
  string dir = "IJazZ_output/";
  if( _doFIT )  {
    if( _ijazzClassic )    dir += "FitResultsClassic/";
    else                   dir += "FitResults/";
  } else if( _outputTree ) dir += "CorrectedTrees/";
  if( _doFIT ) {
    if     ( _analysis == -1 ) dir += "testIJazZ/";
    else if( _analysis == +0 ) dir += _userDir + "/";
    else if( _analysis == +1 ) dir += "ecalEtaScale/";
  }
  dir += itostr(_runMin) + "-" + itostr(_runMax) + "/";

  string nameOutFilePref  = "IJazZ_";
  nameOutFilePref += ecalp + "_";
  if( _isMC ) nameOutFilePref += "MC_";
  else        nameOutFilePref += "Data_";
  nameOutFilePref += _version;
  
  string mkdir = "mkdir -p " + dir;
  gSystem->Exec( mkdir.c_str() );
  return dir + nameOutFilePref + ".root";
}

//--- constructor
IJazZ::IJazZ(void) {
  _debugLevel  = 0;
  _hNvtxMC     = 0;
  _hNvtxData   = 0;
  _hNvtxWeight = 0;
  _ztreeName  = "selected";
  // _ztreeName  = "simpleNtupleEoverP/SimpleNtupleEoverP";
  // _ztreeName  = "SimpleNtupleEoverP";
  _isMC = false;
  _runMin = -1;
  _runMax = 9999999;
  _bin3D = false;
  _analysis = -1;
  _isEEfit = false;
  _meeFitMin =  89;
  _meeFitMax = 100;

  _EEringsDefinitionFile[0] = "etc/geom/eering.dat";
  _EEringsDefinitionFile[1] = "etc/geom/Numbering_EE_All.dat";
  _fitRelResoForRegrEnCorr = false;  /// can fit only the ratio reso / resoExpRegr

  _ptCut = 25;
  _scEnergyCorrType = 0;

  _allEvenOddEvents = 0;

  _zfitterFab = 0;
  _testOversmearingMC = 0;

  _doFIT        = false;
  _outputTree   = false;
  _ijazzClassic = false; // default mode for IJazZ is now fit with tails

  _version = "Test";
  _userDir = "user";
  _fileout = (TFile*)0;
  TFile *f = TFile::Open("etc/geom/IJazZFit.RotateTail.root","read");
  assert( f->IsOpen() );
  TVectorD *pointer = (TVectorD*) f->Get("RotateTail0_ijazz_EB"); 
  _tailPar0Rotate_EB.resize(pointer->GetNrows());
  for( int i = 0 ; i < pointer->GetNrows(); i++ )
    _tailPar0Rotate_EB[i] = (*pointer)(i);

  pointer = (TVectorD*) f->Get("RotateTail0_ijazz_EE"); 
  _tailPar0Rotate_EE.resize(pointer->GetNrows());
  for( int i = 0 ; i < pointer->GetNrows(); i++ )
    _tailPar0Rotate_EE[i] = (*pointer)(i);
  f->Close();
}



IJazZ::~IJazZ(void) {
}


void IJazZ::setupNvtxReweighting(void) {
  string varNvtx = "PV_n";
  _hNvtxMC     = new TH1F("hNvtxMC"    ,"# of vertex MC"  ,99,0.5,99.5);
  _hNvtxData   = new TH1F("hNvtxData"  ,"# of vertex Data",99,0.5,99.5);
  _hNvtxWeight = new TH1F("hNvtxWeight","weight for # of vertex",99,0.5,99.5);

  string strRunRange = "runId >= " + itostr(_runMin ) + " && runId <= " + itostr(_runMax);
  if( _ztreeName == "selected" ) {
    strRunRange = "runNumber >= " + itostr(_runMin ) + " && runNumber <= " + itostr(_runMax);
    varNvtx     = "nPV";
  } 
  TCut cutRunRange( strRunRange.c_str() );
  cout << " -- cut = " << cutRunRange << endl;
  string todraw;
  todraw = varNvtx + ">>hNvtxMC";
  _ztreeMC  ->Draw(todraw.c_str(),"","goff");
  todraw = varNvtx + ">>hNvtxData";
  _ztreeData->Draw(todraw.c_str(),cutRunRange,"goff");

  _hNvtxWeight->Divide(_hNvtxData,_hNvtxMC,
		       1. / _hNvtxData->Integral(),
		       1. / _hNvtxMC  ->Integral() );
}


float IJazZ::getNvtxWeight( int nVtx ) {
  if( !_isMC ) return 1;
  if( !_hNvtxWeight ) {
    cerr << "WARNING: IJazZ::getNvtxWeight(...): For MC one should set the reweighting histograms"
	 <<" before calling this function" << endl
	 << "     - return 1" << endl;
    return 1;
  }
 
  int ibin = _hNvtxWeight->FindBin( nVtx );
  if( ibin <= 1 ) ibin = 1;
  if( ibin >= _hNvtxWeight->GetNbinsX() ) ibin = _hNvtxWeight->GetNbinsX();
  
  float weight = _hNvtxWeight->GetBinContent(ibin);
  if( weight >= 5 ) weight = 5;

  return weight;
}


void IJazZ::etaScaleFromZeeAnaFit(string mcVersion ) {
  int ismc = _isMC;
  string version = _version;
  _isMC = 0;  string datafile = createOutputFileName("EB");
  if( mcVersion != "notdefined" ) _version = mcVersion;
  _isMC = 1;  string   mcfile = createOutputFileName("EB");
  
  etaScaleFromZee_AnaFit(datafile, mcfile );

  _isMC = ismc;
  _version = version;
}




void IJazZ::etaHistoStyle( TGraphErrors *gr, const IJazZAxis<double> *axis,int color, string yTitle, 		
					    int yMode, float yMin, float yMax ) {

  string grName = gr->GetName();
  bool ER = false;
  if( grName.find("ER") != string::npos ) ER = true;
  if( yMode == 0 ) ER = false;
  if( yMode == 1 ) ER = true;

  float xmin = -2.5;
  float xmax = +2.5;
  string xTitle = "#eta";
  if( axis->getName() == "IEta" ) {
    xmin = -120.5;
    xmax = +120.5;
    xTitle = "i#eta ring";
  } 
  
  if( axis->getName().find("Abs") != string::npos ) {
    xmin = 0;
    xTitle = "| #eta |";
  }

  float ymin = 0.97;
  float ymax = 1.04;
  if( ER ) {
    ymin = 0.00;
    ymax = 0.06;
  }
  if( yMin != -1 ) ymin = yMin;
  if( yMax != -1 ) ymax = yMax;

  EtaHistoStyle( gr, xmin, xmax, ymin, ymax, color, yTitle, xTitle );
}




TGraphErrors* IJazZ::graphMerger( int yBin, bool ER, bool isIeta  ) {
  string grName   = "grES";
  if( ER ) grName = "grER";
  grName += "_ir9:" + itostr(yBin);

  vector<TGraphErrors*> gr;
  vector<string> ecalp;
  ecalp.push_back("EB"); 
  ecalp.push_back("EEm");
  ecalp.push_back("EEp");
  
  for( unsigned iecal = 0; iecal < ecalp.size(); ++iecal ) {
    string fin = createOutputFileName(ecalp[iecal]);
    cout << "IJazZ::GraphMerger Open: " << fin << endl;

    gr.push_back( (TGraphErrors*) 0 );
    TFile *f = 0; f = TFile::Open( fin.c_str() , "read" );

    if( f && f->IsOpen() ) gr[iecal] = (TGraphErrors*) f->Get( grName.c_str() );
 
   // if ( doieta == 1)
   //  for( int ip = 0; ip < gr[iEcal]->GetN(); ip++ )
   //    gr[iecal]->GetX()[ip] = convert_eta_ieta( gr[iecal]->GetX()[ip] );

    if( gr[iecal] )
    for( int ip = 0; ip < gr[iecal]->GetN(); ip++ )
     if( fabs(gr[iecal]->GetX()[ip]) > 120.01 ) {
       cout << " Removing point: " << gr[iecal]->GetX()[ip] << " for iEcal = " << ecalp[iecal] << endl;
       gr[iecal]->RemovePoint(ip);    
     }
  }
  
  
  TGraphErrors *grout = (TGraphErrors*) gr[0]->Clone();
  float CCstop = 85.5;
  if( ! isIeta ) CCstop = 1.506;

  for( int ip = 0; ip < grout->GetN(); ip++ ) {
    float eta = grout->GetX()[ip];
    int iEcal = 0;
    if( eta <= -CCstop && gr[1] != 0 ) { 
      iEcal = 1;
      int iq = 0;
      if( grName.find("ER") != string::npos ) {
	for( iq = 0; iq < grout->GetN(); iq++ ) 
	  if( grout->GetX()[ip] == -grout->GetX()[iq] ) break;
      } else iq = ip;
      grout->GetY()[ip]  = gr[iEcal]->GetY()[iq];
      grout->GetEY()[ip] = gr[iEcal]->GetEY()[iq];
    }
    if( eta >= +CCstop && gr[2] != 0 ) {
      iEcal = 2;    
      grout->GetY()[ip]  = gr[iEcal]->GetY()[ip];
      grout->GetEY()[ip] = gr[iEcal]->GetEY()[ip];
    } 
  }
 
 
  //// remove points  
  if( grName.find( "grER") == string::npos )
  for( int ip = 0; ip < grout->GetN(); ip++ ) {
    float eta  = grout->GetX()[ip];
    float deta = grout->GetEX()[ip];

    if(   ( (fabs( eta ) < 0.5 && deta < 1.5 ) || fabs( eta ) > 118.5  )
	  || ( fabs( eta ) >= 83.999  && fabs( eta ) < 88 ) 
	  ) {
      grout->GetEX()[ip] = 0;
      grout->GetY()[ip]  = 1;
      grout->GetEY()[ip] = 0;
    }
  }
  /// miror resolution in EB
  if( grName.find( "grER") != string::npos ){
    for( int ip = 0; ip < grout->GetN(); ip++ ) {
      float eta = grout->GetX()[ip];
      if(   fabs( eta ) < CCstop &&  eta < 0 ) {
	int iq = 0;
	for( iq = 0; iq < grout->GetN(); iq++ ) {
	  if( fabs( grout->GetX()[ip] + grout->GetX()[iq] ) < 0.00001) break;
	}
	
	/// mirror resolution in EB
	grout->GetEX()[ip] = grout->GetEX()[iq];
	grout->GetY()[ip]  = grout->GetY()[iq];
	grout->GetEY()[ip] = grout->GetEY()[iq];
      }
    }
  }
  return grout;
}




