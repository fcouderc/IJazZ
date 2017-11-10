#include "interface/IJazZ.hh"
#include "interface/ZFitterMinuit2_ND_MThread.hh"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "interface/RootUtils.hh"

//#include <TFitterMinuit.h>
#include <TFitter.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TTree.h>

using namespace std;
bool IJazZ::setBiningND( string axisRespDef, string axisResoDef ) {

  //--- define here the X and Y bining as well as the X and Y variables
  _respAxisND.reset();
  _resoAxisND.reset();


  IJazZAxis<double> xBinResp, xBinReso;
  IJazZAxis<double> yBinResp, yBinReso;

  if        ( _analysis == 0 ) {
    _respAxisND.readFromFile(axisRespDef);
    _resoAxisND.readFromFile(axisResoDef);
  } else if( _analysis == 1 ) {
    xBinResp.setName( "IEta" );
    xBinReso.setName( "AbsIEta" );
    yBinResp.setName( "R9" );
    yBinReso.setName( "R9" );

    /// eta bining
    float eta_loc  = -120.5;
    while( eta_loc < 120.9 ) { 
      xBinResp.addBin(eta_loc);
      eta_loc += 1.0;
    }
    
    float deta_loc = 10;
    eta_loc  = -125.5;
    while( eta_loc < 125.999 ) {
      xBinReso.addBin(eta_loc);
      eta_loc += deta_loc;
      if( eta_loc > 0 && eta_loc < deta_loc/2. - 0.0001 ) eta_loc = deta_loc/2. + 0.5;
    }
    
    /// r9 bining
    double r9TMP[] = {-100,0.94,100};
    yBinResp.setBining(2,r9TMP);
    yBinReso.setBining(2,r9TMP);    
  } else {
    xBinResp.setName( "AbsSCEta" );
    xBinReso.setName( "AbsSCEta" );
    yBinResp.setName( "R9" );
    yBinReso.setName( "R9" );
   
    double etaTMP[] = {0,0.75,1.5,2.00,2.5};
    double r9TMP[]  = {-100,0.94,+100};
    xBinResp.setBining(  4, etaTMP );
    yBinResp .setBining( 2, r9TMP  );
    xBinReso = xBinResp;
    yBinReso  = yBinResp ;
  }

  if( _analysis != 0 ) {
    _respAxisND.addAxis( xBinResp );
    _resoAxisND.addAxis( xBinReso );
    _respAxisND.addAxis( yBinResp );
    _resoAxisND.addAxis( yBinReso );
  }

  IJazZAxis<double> xSingleBin, ySingleBin, single;
  xSingleBin.setName("AbsSCEta");
  xSingleBin.addBin(0.0);
  xSingleBin.addBin(0.5);
  xSingleBin.addBin(1.0);
  xSingleBin.addBin(1.5);
  xSingleBin.addBin(2.0);
  xSingleBin.addBin(2.5);
  ySingleBin.setName("R9");
  ySingleBin.addBin(0.0);
  ySingleBin.addBin(0.94);
  ySingleBin.addBin(5);
  single.setName("AbsSCEta");
  single.addBin(0.0);
  single.addBin(1.5);
  single.addBin(2.5);

  _tailAxisND.resize(7);
  /// for the tail axis, only use eta
  // _tailAxisND[0].addAxis( xSingleBin );
  // //  _tailAxisND[0].addAxis( ySingleBin );
  // //  _tailAxisND[1].addAxis( single );
  // //  _tailAxisND[2].addAxis( single );
  // //  _tailAxisND[1].addAxis( xSingleBin );
  // //  _tailAxisND[1].addAxis( ySingleBin );
  // _tailAxisND[2].addAxis( xSingleBin );
  // _tailAxisND[2].addAxis( ySingleBin );
  // if( _resoAxisND.nBinsND() < _tailAxisND[2].nBinsND() ) _tailAxisND[2] = _resoAxisND;
  // if( _resoAxisND.nBinsND() < _tailAxisND[0].nBinsND() ) _tailAxisND[0] = _resoAxisND;


  _tailAxisND[0].readFromFile("etc/geom/MCTruthParameters.bin.fittedResp");
  _tailAxisND[1].readFromFile("etc/geom/MCTruthParameters.bin.fittedReso");
  _tailAxisND[2].readFromFile("etc/geom/MCTruthParameters.bin.fittedTailP0");
  _tailAxisND[3].readFromFile("etc/geom/MCTruthParameters.bin.fittedTailP1");
  _tailAxisND[4].readFromFile("etc/geom/MCTruthParameters.bin.fittedTailP2");
  _tailAxisND[5].addAxis( xSingleBin );
  _tailAxisND[5].addAxis( ySingleBin );
  _tailAxisND[6].addAxis( xSingleBin );
  _tailAxisND[6].addAxis( ySingleBin );

  for( int ib = 0 ; ib < _tailAxisND[5].nBinsND(); ib++ ) _tailAxisND[5].value(ib)  = 1;
  for( int ib = 0 ; ib < _tailAxisND[6].nBinsND(); ib++ ) _tailAxisND[6].value(ib)  = 1;

  cout << " =================== Axis for Response   parameter =====================" << endl;  _respAxisND.print();
  cout << " =================== Axis for Resolution parameter =====================" << endl;  _resoAxisND.print();
  //  cout << " =================== Axis for Tail0      parameter =====================" << endl;  _tailAxisND[0].print();
  //  cout << " =================== Axis for Tail2      parameter =====================" << endl;  _tailAxisND[2].print();

  if( _respAxisND.getAxis(0)->getName().find("Abs") != string::npos ||
      _respAxisND.getAxis(0)->getName().find("abs") != string::npos ) return true;
  return false;
}



//---------------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------ Saving fit result --------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------//
#include <TH2F.h>
void IJazZ::saveFitResults( ROOT::Math::Minimizer * fitter, string * parName) {    
  _fileout->cd();
  _fileout->Print();  
  double * par0, *parE;
  par0 = new double[_nParameters];
  parE = new double[_nParameters];

  for( int  ip = 0 ; ip < _nParameters; ip++ ) {
    par0[ip] = fitter->X()[ip];  
    parE[ip] = fitter->Errors()[ip];
    //    parE[ip] = sqrt(fitter->CovMatrix(ip, ip ));
  }

  for( int ibin = 0 ; ibin < _respAxisND.nBinsND(); ibin++ ) { 
    _respAxisND.value(ibin) = 1+par0[ibin]*0.01;
    _respAxisND.error(ibin) = parE[ibin]*0.01;
  }

  for( int ibin = 0 ; ibin < _resoAxisND.nBinsND(); ibin++ ) { 
    int iibin = ibin + _respAxisND.nBinsND();
    _resoAxisND.value(ibin) = par0[iibin]*0.01;
    _resoAxisND.error(ibin) = parE[iibin]*0.01;
  }


  if( ! _ijazzClassic ) {
  int shift =  _respAxisND.nBinsND() + _resoAxisND.nBinsND();
  for( unsigned itail = 0 ; itail < _tailAxisND.size(); itail++ ) {
    for( int ibin = 0 ; ibin < _tailAxisND[itail].nBinsND(); ibin++ ) { 
      int iibin =   ibin + shift;
      _tailAxisND[itail].value(ibin) = par0[iibin];
      _tailAxisND[itail].error(ibin) = parE[iibin];
    }
    shift += _tailAxisND[itail].nBinsND();   
  }
  }


  IJazZAxis<double> esResp;
  esResp.setName("ESreponse");
  esResp.addBin(0.5); esResp.addBin(1.5); esResp.addBin(2.5);
  IJazZAxisND<double> esRespAxis; esRespAxis.addAxis( esResp );
  vector<double> x(1); 
  x[0] = 1;  esRespAxis.value(x) = par0[ _nParameters-2]; esRespAxis.error(x) = parE[_nParameters-2];
  x[0] = 2;  esRespAxis.value(x) = par0[ _nParameters-1]; esRespAxis.error(x) = parE[_nParameters-1];
    
  TCanvas *cCrossCheck = _zfitterFab->fitCrossCheck(par0);
  cCrossCheck->Write( cCrossCheck->GetName(),  TObject::kOverwrite );
  cout << "---- IJazZ: crosscheck can saved to: " << _fileout->GetName() << endl;
  cout << "     IJazZ: reponse    fit saved to: " << string(_fileout->GetName())+".fittedResp" << endl;
  cout << "     IJazZ: resolution fit saved to: " << string(_fileout->GetName())+".fittedReso" << endl;
  cout << "     IJazZ: responseES fit saved to: " << string(_fileout->GetName())+".fittedRespES" << endl;

  _respAxisND.saveToFile(  string(_fileout->GetName())+".fittedResp");
  _resoAxisND.saveToFile(  string(_fileout->GetName())+".fittedReso");
  esRespAxis .saveToFile(  string(_fileout->GetName())+".fittedRespES");
  /*
    for( unsigned itail = 0 ; itail < _tailAxisND.size(); itail++ )
    _tailAxisND[itail].saveToFile(  string(_fileout->GetName())+".fittedTailP" + itostr(itail) );
  */

  
  /// save correlation matrix
  TH2F *hCorrelation = new TH2F( "hCorrelation","correlation",
				 fitter->NFree(),0.5,fitter->NFree()+0.5,
				 fitter->NFree(),0.5,fitter->NFree()+0.5 );
  TH2F *hCovariance  = new TH2F( "hCovariance","covariance",
				 fitter->NFree(),0.5,fitter->NFree()+0.5,
				 fitter->NFree(),0.5,fitter->NFree()+0.5 );	
  hCorrelation->SetDirectory(0);
  hCovariance ->SetDirectory(0);
  
  /// _nParameters can be much larger than # of floated parameters
  /// make a correspondence 
  vector<int> parAllToParFree(_nParameters); 
  int iFree = 0;
  for( int ip = 0 ; ip < _nParameters; ip++ )
   if( fitter->CovMatrix(ip,ip) > 0.00000001 ) {
     parAllToParFree[ip] = iFree; iFree++; 
   } else parAllToParFree[ip] = -1;
  
  for( int i = 0 ; i < _nParameters; i++ )
  for( int j = i ; j < _nParameters; j++ )
   if( parAllToParFree[i] >= 0 &&  parAllToParFree[j] >= 0 ) {
     hCorrelation->SetBinContent( parAllToParFree[i]+1, 
				  parAllToParFree[j]+1,  
				  fitter->Correlation(i,j) );
     hCovariance ->SetBinContent( parAllToParFree[i]+1, 
				  parAllToParFree[j]+1,  
				  fitter->CovMatrix(i,j) );
     hCorrelation->SetBinContent( parAllToParFree[j]+1, 
				  parAllToParFree[i]+1,  
				  fitter->Correlation(i,j) );
     hCovariance ->SetBinContent( parAllToParFree[j]+1, 
				  parAllToParFree[i]+1,  
				  fitter->CovMatrix(i,j) );
     
  }
  for( int i = 0 ; i < _nParameters; i++ )
   if( parAllToParFree[i] >= 0 ) {
    hCorrelation->GetXaxis()->SetBinLabel(  parAllToParFree[i]+1, parName[i].c_str() );
    hCorrelation->GetYaxis()->SetBinLabel(  parAllToParFree[i]+1, parName[i].c_str() );
    hCovariance->GetXaxis()->SetBinLabel(  parAllToParFree[i]+1, parName[i].c_str() );
    hCovariance->GetYaxis()->SetBinLabel(  parAllToParFree[i]+1, parName[i].c_str() );

  }

  delete[] par0;
  delete[] parE;
  hCorrelation->Write();
  hCovariance ->Write();
  _fileout->Close();
}


#include <TRandom3.h>
//---------------------------------------------------------------------------------------------------------------------------//
//------------------------------------------- setup and minize with ZFitter -------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------//
void IJazZ::setupZFitter( int ncpu ) {
  _meeFitMin =  80.0;
  _meeFitMax = 100.0;
  if( _isEEfit ) { 
    _meeFitMin =  80;
    _meeFitMax = 100;
  }
  if( _ijazzClassic ) {
    _meeFitMin =  89.0;
    _meeFitMin =  82.0; /// first 2015 data
    _meeFitMax = 100.0;
    if( _isEEfit ) { 
      _meeFitMin =  88;
      _meeFitMin =  75;  /// first 2015 data
      _meeFitMax = 100;
      if( _allEvenOddEvents == 1 ) _meeFitMin =  84; /// add more evts for dedicated scale fit 
    }
    if( _scEnergyCorrType == 2 ) { _meeFitMin -= 1.0; _meeFitMax -= 1.0; }
    //    _meeFitMin =  82.5;
  }

  _nParameters = _respAxisND.nBinsND() + _resoAxisND.nBinsND() + 2;
  if( !_ijazzClassic ) for( unsigned ip = 0 ; ip <  _tailAxisND.size(); ip++ ) _nParameters += _tailAxisND[ip].nBinsND();

  _zfitterFab = new ZFitterMinuit2();
  _zfitterFab->setNParams(_nParameters);
  _zfitterFab->setNumCPU(ncpu);
  _zfitterFab->setFitRange(  _meeFitMin, _meeFitMax);
  _zfitterFab->setTailParRotation(  _tailPar0Rotate_EB, _tailPar0Rotate_EE );
  if( _ijazzClassic ) _zfitterFab->setFitterModeClassic();
}



void IJazZ::minimize( string ecalp ) {
  TCanvas *cAdjust = _zfitterFab->fitAutoAdjustRange(_meeFitMin,_meeFitMax);
  _zfitterFab->setFitRange(  _meeFitMin, _meeFitMax);

  _fileout->cd();
  _fileout->Print();
  cAdjust->Write( cAdjust->GetName(),  TObject::kOverwrite );
  //  _fileout->Close();
  //  return;


  string tailParName[] = { "Mu", "Sig", "dMu", "rSig",  "fMain", "KdMu", "KrSig" };
  


  double *par0, *parE, *parMin, *parMax;
  string *parName;
  bool   *parFixed;
  par0     = new double[_nParameters];
  parE     = new double[_nParameters];
  parMin   = new double[_nParameters];
  parMax   = new double[_nParameters];
  parName  = new string[_nParameters];
  parFixed = new bool  [_nParameters];

  for( int ip = 0; ip < _nParameters; ip++ ) parFixed[ip] = false;  
  cout << "IJazZ: # of parameters in fit (including fixed param): " << _nParameters << endl;
  int iip = _nParameters-2;
  parName[iip] = "rES1"; par0[iip] = 1.0; parE[iip] = 0.001; parMin[iip] = 0.2; parMax[iip] = 5.0; parFixed[iip] = 1;
  iip = _nParameters-1;
  parName[iip] = "rES2"; par0[iip] = 1.0; parE[iip] = 0.001; parMin[iip] = 0.2; parMax[iip] = 5.0; parFixed[iip] = 1;
  
  
  TRandom3 gen;
  for( int ib = 0 ; ib < _respAxisND.nBinsND(); ib++ ) {
    parName[ib] = "response_bin:" + itostr(ib);
    //    par0[ib] = gen.Uniform(-2,2); parE[ib] = 0.02; parMin[ib] = -0.05/0.01; parMax[ib] = 0.05/0.01;
    par0[ib] = 0; parE[ib] = 0.002; parMin[ib] = -0.5/0.01; parMax[ib] = 0.5/0.01;
  }
  
  for( int ib = 0 ; ib < _resoAxisND.nBinsND(); ib++ ) {
    int iib = ib +  _respAxisND.nBinsND();
    parName[iib] = "resolution_bin:" + itostr(ib);
    par0[iib] = 0.03/0.01; parE[iib] = 0.05; parMin[iib] = +0.003/0.01; parMax[iib] = 1.0/0.01; 
    if( _fitRelResoForRegrEnCorr ) { par0[iib] = 1; parE[iib] = 0.0001; parMin[iib] = 0.5; parMax[iib] = 2.0; }
  }
  
  int shift =  _respAxisND.nBinsND() + _resoAxisND.nBinsND();

  if( ! _ijazzClassic )
  for( unsigned itail = 0 ; itail < _tailAxisND.size(); itail++ ) {
    for( int ib    = 0 ; ib < _tailAxisND[itail].nBinsND(); ib++ ) {
      int iib = ib +  shift;
      parName[iib] = tailParName[itail] + "_bin:"   + itostr(ib);
      par0[iib]    = _tailAxisND[itail].getValue(ib); parE[iib] = 0.01;
      parMin[iib] = -1000.00; parMax[iib] = +1000.0; 
      if( itail > 5 ) parMin[iib] = 0.5; parMax[iib] = 5;
      // if( itail == 0 ) { par0[iib] =   2.49; parE[iib] = 0.0001; parMin[iib] = +1.00; parMax[iib] = +3.0; }
      // if( itail == 1 ) { par0[iib] =   3.7; parE[iib] = 0.0001; parMin[iib] = +0.10; parMax[iib]  =  5.0; }
      // if( itail == 2 ) { par0[iib] =  0.80; parE[iib] = 0.0001; parMin[iib] = +0.40; parMax[iib] =  0.95; }       
      // if( itail == 1 ) par0[iib] = _tailAxisND[itail].getValue(ib);
    }
    shift += _tailAxisND[itail].nBinsND();
  }
  
  //--- for EE fit, fix the EB values from EB fit
  if( _isEEfit ) {
    string fileFit = createOutputFileName("EB");
    IJazZAxisND<double> EBfitResp; EBfitResp.readFromFile( fileFit + ".fittedResp" );
    IJazZAxisND<double> EBfitReso; EBfitReso.readFromFile( fileFit + ".fittedReso" );
    for( int ib = 0 ; ib < EBfitResp.nBinsND(); ib++ )  par0[ib                    ] = (EBfitResp[ib]-1)/0.01;
    for( int ib = 0 ; ib < EBfitReso.nBinsND(); ib++ )  par0[ib+EBfitResp.nBinsND()] = (EBfitReso[ib])/0.01;
    // for( int ib = 0 ; ib < EBfitResp.nBinsND(); ib++ )  parE[ib                    ] = EBfitResp.error(ib)/0.01;
    // for( int ib = 0 ; ib < EBfitReso.nBinsND(); ib++ )  parE[ib+EBfitResp.nBinsND()] = EBfitReso.error(ib)/0.01;
  }
  

  float CCstop = 1.479;
  if( _respAxisND.getAxis(0)->getName().find( "IEta" ) != string::npos ) CCstop = 85.2;
  
  if( _debugLevel >= 1 ) cout << "IJazZ:  ====== fixing parameters ======= " << endl;
  if( _debugLevel >= 1 ) cout << "CC stop = " << CCstop << endl;
  for( int ib = 0 ; ib < _respAxisND.nBinsND(); ib++ ) {
    float etaBinCenter = _respAxisND.getBinCenterDimN(ib,0);
    bool CC =  etaBinCenter < CCstop+0.00001 && etaBinCenter > -CCstop-0.00001;
    bool fixParam = false;
    if(   CC &&  _isEEfit ) fixParam = true;
    if(  !CC && !_isEEfit ) fixParam = true;
    if( etaBinCenter > +CCstop+0.00001 && ecalp == "EEm" ) fixParam = true;
    if( etaBinCenter < -CCstop-0.00001 && ecalp == "EEp" ) fixParam = true;
    
    if( etaBinCenter > +CCstop+0.00001 && ecalp == "EEm" && 
  	( _respAxisND.getAxis(0)->getName().find("Abs") != string::npos ) ) fixParam = false;
    
    if( fixParam ) parFixed[ib] =  1;
    
    if( _debugLevel >= 3 ) {
      cout << "    - response  ,  par{ " << ib << " -> Bin[1d] = " << etaBinCenter;
      for( unsigned id = 1; id < _respAxisND.getND();id++ ) cout << ", Bin[" << id+1 <<"d] = " << _respAxisND.getBinCenterDimN(ib,id);
      cout << "} ";
      if( fixParam ) cout << " fixed to: " << par0[ib] << endl; 
      else           cout << " free" << endl;
    }
  }
  
  CCstop = 1.479;
  if( _resoAxisND.getAxis(0)->getName().find( "IEta" ) != string::npos ) CCstop = 85.2;
  for( int ib = 0 ; ib < _resoAxisND.nBinsND(); ib++ ) {
    float etaBinCenter = _resoAxisND.getBinCenterDimN(ib,0);
    bool CC =  etaBinCenter < CCstop+0.00001 && etaBinCenter > -CCstop-0.00001;
    bool fixParam = false;
    if(  CC &&  _isEEfit )  fixParam = true;
    if( !CC && !_isEEfit )  fixParam = true; 
    //    if( etaBinCenter < -0.00001 ) fixParam = true;
    
    if( fixParam ) parFixed[ib+_respAxisND.nBinsND()] = 1;
    
    if( _debugLevel >= 3 ) {
      if( fixParam && etaBinCenter < -0.00001 ) { cout << " res[eta<0] are fixed to 1 (only |eta|): " ;}
      cout << "    - resolution,  par{ " << ib+ _respAxisND.nBinsND()  << " -> Bin[1d] = " << etaBinCenter;
      for( unsigned id = 1; id < _resoAxisND.getND();id++ ) cout << ", Bin[" << id+1 <<"d] = " << _resoAxisND.getBinCenterDimN(ib,id);
      cout << "} ";
      if( fixParam ) cout << " fixed to: " << par0[ib+ _respAxisND.nBinsND()] << endl; 
      else   	     cout << " free" << endl;
    }
  }
  
  if( !_ijazzClassic ) {
    shift =  _respAxisND.nBinsND() + _resoAxisND.nBinsND();
    for( unsigned itail = 0 ; itail < _tailAxisND.size(); itail++ ) {
      CCstop = 1.479;
      if( _tailAxisND[itail].getAxis(0)->getName().find( "IEta" ) != string::npos ) CCstop = 85.2;
      for( int ib    = 0 ; ib < _tailAxisND[itail].nBinsND(); ib++ ) {
	int iParNumber =  ib + shift;
	float etaBinCenter = _tailAxisND[itail].getBinCenterDimN(ib,0);
	bool CC =  etaBinCenter < CCstop+0.00001 && etaBinCenter > -CCstop-0.00001;
	bool fixParam = false;
	if(  CC &&  _isEEfit )  fixParam = true;
	if( !CC && !_isEEfit )  fixParam = true; 
	if( itail < 5 ) fixParam = true;
	//      fixParam = true;
	// if( etaBinCenter < -0.00001 ) fixParam = true;
	
	// if( itail == 0 ) fixParam = true;
	// if( itail == 1 ) fixParam = true;
	// if( itail == 0 && !CC )
	// 	minimizer.SetParameter( iParNumber, parName[iParNumber].c_str(),
	// 				2.31 , 0.0001, 1.50000, 2.500001 );
	// if( itail == 0 &&  CC && !_isEEfit )
	// 	/// do not change EB value for EE fits since this is fixed from previous fit
	// 	minimizer.SetParameter( iParNumber, parName[iParNumber].c_str(),
	// 				2.49 , 0.0001, 0.5000, 2.900001 );
	if( _ijazzClassic ) fixParam = true;
	if( fixParam ) parFixed[iParNumber] = 1;
	
	if( _debugLevel >= 3 ) {      
	  cout << "    - " << tailParName[itail] << ",  par{ " << iParNumber 
	       << " -> Bin[1d] = " << etaBinCenter;
	  for( unsigned id = 1; id < _tailAxisND[itail].getND();id++ ) cout << ", Bin[" << id+1 <<"d] = " << _tailAxisND[itail].getBinCenterDimN(ib,id);
	  cout << "} ";
	  if( fixParam ) cout << " fixed to: " << par0[iParNumber] << endl; 
	  else   	     cout << " free" << endl;
	}
      }
      shift += _tailAxisND[itail].nBinsND();
    }
  }
  
  vector<int> nEvtPerParam = _zfitterFab->nEvtPerParam();
  for( unsigned ip = 0; ip < nEvtPerParam.size()-2; ip++ )
    if( !parFixed[ip] && nEvtPerParam[ip] < 10 ) {
      if( _debugLevel >= 2 ) cout << " - Fix parameter(" << ip << ") because it has a low statistic, only nEvt: " << nEvtPerParam[ip] <<endl;
      parFixed[ip] = 1;
    }
  

 
  /// define the function to minimize
  //string minimizerName = "Minuit2"; string minimizerAlgo  = "Fumili2";
  string minimizerName = "Minuit2"; string minimizerAlgo  = "Migrad";
  //  string minimizerName = "Minuit";  string minimizerAlgo  = "Migrad";
  // string minimizerName = "Minuit";  string minimizerAlgo  = "Simplex";
  //  string minimizerName = "GSLMultiMin";  string minimizerAlgo  = "BFGS2";
  //  string minimizerName = "GSLSimAn"; string minimizerAlgo  = "";
  ROOT::Math::Functor FCN(_zfitterFab,&ZFitterMinuit2::eval,_nParameters);
  ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer(minimizerName.c_str(),minimizerAlgo.c_str());
  minimizer->SetFunction( FCN );

  cout << "======================= setting up fitter with free and fixed variables ========================" << endl;
  for( int ip = 0 ; ip < _nParameters; ip++ ) {
    if( parFixed[ip] ) minimizer->SetFixedVariable(   ip, parName[ip].c_str(), par0[ip] );
    else               {
      minimizer->SetLimitedVariable( ip, parName[ip].c_str(), par0[ip], parE[ip], parMin[ip], parMax[ip] );
      if( _debugLevel >= 1 ) cout << " Par{ " << parName[ip] << " } free, init to: " << par0[ip] << " step: " << parE[ip] << endl;
    }
  }


  /// --- count number of free parameters in this fit:
  unsigned nParFitted = 0;
  for( unsigned ip = 0; ip < nEvtPerParam.size(); ip++ ) 
    if( !parFixed[ip] ) nParFitted++;
  
  cout << "---- IJazZ: # of param free to float in this fit: " << nParFitted << " / " << _nParameters << endl;

  double tolerance = 0.1;
  minimizer->SetPrintLevel(_debugLevel);
  //minimizer->SetPrintLevel(2);
  minimizer->SetMaxFunctionCalls(500000); // for Minuit/Minuit2
  minimizer->SetMaxIterations(500);  // for GSL
  minimizer->SetStrategy(1);   
  // minimizer->SetStrategy(0);
  if( nParFitted > 100 )  minimizer->SetStrategy(0);   
  //  if( nParFitted > 50 )  tolerance = 10;
  //  if( nParFitted > 50 )  minimizer->SetValidError(0);
  minimizer->SetTolerance(tolerance);    
  cout << "Minimizing with:  " << minimizerName << " (algo: " 
       << minimizerAlgo << "), tolerance: " << minimizer->Tolerance() << " , Up: " << minimizer->ErrorDef() << endl
       << " ==> ValidError?: " << minimizer->IsValidError() << endl
       << " ==> precision: " << minimizer->Precision() << endl;
  minimizer->Minimize();
  cout << "entering save fit" << endl;
  // run with Hesse for errors if not many parameters.
  // if( nParFitted < 100 )  minimizer->Hesse();
  // minimizer->PrintResults();
  saveFitResults(minimizer, parName);

  delete[] par0;
  delete[] parE;
  delete[] parMin;
  delete[] parMax;
  delete[] parName;
  delete[] parFixed;
  delete minimizer;

}

