//---------------------------------------------------------------------------//
// Description: IJazZ event looper, very little selections are actually perform
//              here, just what is necessary for a potential fit
//              it loops over IJazZ diele class that can be defined by the user
//              Actually all the i/o are defined by the user class IJazZ diele
//                       user event selection criteria are also defined there.
//
//
// Author: Fabrice Couderc, fabrice.couderc@cea.fr 
// History: 2012/11/28 - v1
//
//---------------------------------------------------------------------------//

#include "interface/IJazZ.hh"
#include "interface/IJazZ_tupleVar.hh"
#include "interface/ZFitterMinuit2_ND_MThread.hh"
#include "interface/EcalUtils.hh"
#include "interface/EnergyCorrectionAndSmearingHgg.hh"

#include <TFile.h>
#include <TChain.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
//#include <TFitterMinuit.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TH1F.h>
using namespace std;

float mcSmearingTest( float sceta, float r9 );


double generateLorentZ(TRandom3 & gen) {

  double mZ     = 91.1876;
  double gammaZ = 2.4952;
  double gZO2 = 0.5*gammaZ; gZO2 *= 0.5*gammaZ;

  double mllOK = -1;
  bool goON = true;
  while( goON ) {

    double mll = gen.Uniform(mZ-10,mZ+10);
    
    double FllOverFmax = 1 / ( 1 + (mll-mZ)*(mll-mZ)/gZO2 );
    double y = gen.Uniform(0,1);
    if( y < FllOverFmax ) {
      mllOK = mll;
      goON = false;
    }
  }

  if( mllOK < 0 ) cout << "mllOK < 0, this can not be true! " << endl;

  return mllOK;

}

void IJazZ::eventSelectionZFitter( string ecalp ) {
  TRandom3 gen;
  //--- response correction (data only)
  vector<IJazZAxisND<double> > eRespCorr; 
  eRespCorr.resize(_file_eRespCorr.size());
  vector<runRangeChecker> runCheck;
  for( unsigned icorr = 0 ; icorr < _file_eRespCorr.size(); icorr++ ) {
    cout << "---- IJazZ apply electron correction number: " << icorr << endl;
    eRespCorr[icorr] = combineEtaScale( _file_eRespCorr[icorr] );
    runCheck.push_back( runRangeChecker(_file_eRespCorr[icorr]) );
  }
  
  //--- ES response
  double rES1(1), rES2(1);
  if(  _file_eRespCorr.size() > 0 ) {
    IJazZAxisND<double> eRespCorrES;
    string fileForES =  _file_eRespCorr[0] + "ES"; 
    eRespCorrES.readFromFile( fileForES);
    vector<double> x(1);
    x[0] = 1; rES1 = eRespCorrES.value(x);
    x[0] = 2; rES2 = eRespCorrES.value(x);
    cout << " *** ES plane 1 weight = " << rES1 << endl;
    cout << " *** ES plane 2 weight = " << rES2 << endl;
  }

  //--- oversmearing (MC only)
  bool overSmearMC = false;
  IJazZAxisND<double> eResoData;
  vector<IJazZAxisND<double> > eResoMC;
  if( _resoDataFile.size() > 0 && _resoMCFile.size() > 0 ) {
    for( unsigned i = 0 ; i < _resoMCFile.size(); i++ )
      eResoMC.push_back( combineEtaScale(_resoMCFile[i]  ) );
    eResoData   = combineEtaScale(_resoDataFile[_resoDataFile.size()-1]);
    overSmearMC = true;
  }

  //--- Fix ME should be a class variable tunable in the config file
  bool _7TeV = false;

  TChain *ztree = 0;
  if( _isMC )  ztree = _ztreeMC  ;
  else         ztree = _ztreeData;
  
  //--- input diele variables
  IJazZdiEle diele(ztree);
  diele.set7TeV(_7TeV);
  diele.isMC(_isMC);
  diele.setEnergyCorrectionType(_scEnergyCorrType);
  diele.setBranchAddress( _testOversmearingMC );
  diele.setBranchAddress( 1 );

  //--- output (potential output tree)
  string outName = createOutputFileName(ecalp);
  if( _outputTree ) outName = createOutputFileName("Ecal");
  cout << "---- IJazZ create (overwrite if necessary) file: " << outName << endl;
  _fileout = TFile::Open( outName.c_str(), "recreate");
  if( _outputTree ) diele.createOutputTree( _fileout );
  
  
  //--- Hgg corrections
  EnergyScaleReader scaleHgg;
  scaleHgg.setup("etc/geom/EtDependentCorr.txt");
  photonOverSmearing overSmearHgg;
  overSmearHgg.initialize("Legacy2013_8TeV");

  int printNevents = 500000;
  if( _debugLevel >= 1 )  printNevents = 10000;
  int fittedEvents = 0;
  int nTotEvents   = 0;
  float CCstop = 1.479;
  int nentries = diele.getEntries();
  //  nentries /= 5;
  // nentries = 1000000;
  cout << "ZFitterFabrice: Total number of events: " << nentries << endl;
  
  vector<int>    nEvtCut;
  vector<string> nameCut;
  for( int ievt = 0; ievt < nentries; ievt++ ) {
    unsigned icut = 0;

    diele.loadEvent(ievt);
    if( nEvtCut.size() < icut+1 ) {
      nEvtCut.push_back(0);
      nameCut.push_back("- nEvt[Total    ]: " );
    }
    nEvtCut[icut++]++;

    if( ievt%printNevents == 0) cout << " Event: " << ievt << "\r" << flush ;
    if( !_isMC && ( diele.runId < _runMin || diele.runId > _runMax ) ) continue;
    nTotEvents++;   
    if( nEvtCut.size() < icut+1 ) {
      nEvtCut.push_back(0);
      nameCut.push_back("- nEvt[run      ]: " );
    }
    nEvtCut[icut++]++;

    
    if( ! diele.passPrivateSelection() ) continue;
    if( nEvtCut.size() < icut+1 ) {
      nEvtCut.push_back(0);
      nameCut.push_back("- nEvt[cuts base]: " );
    }
    nEvtCut[icut++]++;


    if( nEvtCut.size() < icut+1 ) {
      nEvtCut.push_back(0);
      nameCut.push_back("- nEvt[evtQuality]: " );
    }
    if( _allEvenOddEvents != 1 && diele.evtQuality < 5 ) continue; /// not fit on subsample (tight selection)
    nEvtCut[icut++]++;


    if( nEvtCut.size() < icut+1 ) {
      nEvtCut.push_back(0);
      nameCut.push_back("- nEvt[odd/even]: " );
    }
    if(  _allEvenOddEvents == 1 && 
     	 (  diele.evtQuality >= 5 && diele.evtId%2 == 0 ) ) continue; /// fit
    if(  _allEvenOddEvents == 2   && diele.evtId%2 != 0   ) continue; /// test    
    nEvtCut[icut++]++;


    if( nEvtCut.size() < icut+1 ) {
      nEvtCut.push_back(0);
      nameCut.push_back("- nEvt[nPV<50  ]: " );
    }
    if( diele.nPV > 50 ) continue;
    nEvtCut[icut++]++;

    diele.weight = 1;
    if( _isMC ) diele.weight *= getNvtxWeight( diele.nPV );
    
    //----  define i1 and i2 for a potential asymmetric cut (but never used for now)
    int i1 = 0, i2 = 1;
    if( diele.eSC[1]/cosh(diele.eta[1]) > diele.eSC[0]/cosh(diele.eta[0]) ) { i1 = 1; i2 = 0; } 
    
    /// when not in eta scale mode, do not correct for ES
    //    if( _analysis != 1 ) for( int iele = 0; iele < 2; iele++ ) diele.esOverEcal[iele] = 0;
    
    
    ///  potentially correct electron energy (only for data)
    if( !_isMC || true )  {
      for( unsigned icorr = 0 ; icorr < eRespCorr.size(); icorr++ ) 
       if( runCheck[icorr].inRange( diele.runId ) ) {
	vector<double> xERespCorr[2]; 
	double respEle[2] = {1,1};
	diele.fill_Xvar_forBinFinder( & eRespCorr[icorr], xERespCorr );
	for( int iele = 0 ; iele < 2; iele++ ) {
	  //	  respEle[iele] = 
	  //	    ( eRespCorr[icorr].value(xERespCorr[iele]) * diele.raw_eSC[iele]
	  //	      + rES1 * diele.raw_eES1[iele] + rES2 * diele.raw_eES2[iele] ) / 
	  //	    ( diele.raw_eSC[iele] + diele.raw_eES1[iele] + diele.raw_eES2[iele] )	    ;
	  respEle[iele] = (diele.esOverEcal[iele] + 1 ) / 
	    (diele.esOverEcal[iele] + 1./eRespCorr[icorr].value(xERespCorr[iele]) );
	}
	diele.corrResp( respEle );
       }
    }
    
    if( _isMC && overSmearMC ) {
      double resoData[2] = {1,1};
      double resoMC[2]   = {1,1};
      for( int iele = 0 ; iele < 2; iele++ ) {
	vector<double> xEResoData[2]; 
	diele.fill_Xvar_forBinFinder( & eResoData, xEResoData );
	resoData[iele] = eResoData.value( xEResoData[iele] );
	vector<double> xEResoMC[2]; 
	double overSmearing2 = 0;
	diele.fill_Xvar_forBinFinder( & eResoMC[0]  , xEResoMC );
	for( int ios =  eResoMC.size()-1; ios >= 0; ios-- ) {
	  resoMC[iele]   = eResoMC[ios]  .value( xEResoMC[iele] );
	  overSmearing2 += resoData[iele]*resoData[iele] - resoMC[iele]*resoMC[iele];
	}
	if( overSmearing2 >= 0   ) {
	  diele.eSC[iele] *= gen.Gaus(1.,sqrt(overSmearing2));
	  diele.pt[iele]   = diele.eSC[iele]/cosh(diele.eta[iele]);
	}
      }
    }


    // ---- Hgg corrections
    if( false ) {
    if( !_isMC ) {
      /// correct data scale with Hgg official corrections
      double respEle[2] = {1,1};
      for( int iele = 0 ; iele < 2; iele++ ) 
	respEle[iele] = 1./scaleHgg.energyScale( diele.R9[iele], diele.scEta[iele],
						 diele.eSC[iele]/cosh(diele.eta[iele]), diele.runId, 0);
      diele.corrResp( respEle );
    } else {
      /// oversmear MC with Hgg official smearing
      for( int iele = 0 ; iele < 2; iele++ ) {
	diele.eSC[iele] *= overSmearHgg.randOverSmearing( diele.scEta[iele],diele.R9[iele],diele.eSC[iele]/cosh(diele.eta[iele]),0);
	diele.pt[iele]   = diele.eSC[iele]/cosh(diele.eta[iele]);
      }
    }
    }

    double mee_true = diele.mZ_MC;
    /*
    if(  _testOversmearingMC > 0 ) {
      for( int iele = 0 ; iele < 2; iele++ ) {
	diele.eSC[iele] = diele.energyMCEle[iele];
	float smear = 0.025;
	if( _testOversmearingMC == 4 ) smear =  mcSmearingTest(diele.scEta[iele],diele.R9[iele]);
	if( _testOversmearingMC >= 3 ) diele.eSC[iele] *= gen.Gaus(1+0.015*cos(diele.scEta[iele]*TMath::TwoPi()/2.5),smear);
	diele.pt[iele]  = diele.eSC[iele]/cosh(diele.etaMCEle[iele]);
	diele.eta[iele] = diele.etaMCEle[iele];
	diele.phi[iele] = diele.phiMCEle[iele];
      }      
    }
    */
    if(  _testOversmearingMC > 0 ) {
      double mee_gen = generateLorentZ(gen);
      double mee_rec = diele.mee();
      double corr_resp[2] = { 1./ ( mee_gen / mee_rec ),  1./ ( mee_gen / mee_rec ) };
      double smear = 0;
      double scale = 0;
      for( int iele = 0 ; iele < 2; iele++ ) {
	if(                                  fabs(diele.scEta[iele]) < 1.0 ) { smear = 0.015; scale = +0.01; }
	if( fabs(diele.scEta[iele]) > 0.5 && fabs(diele.scEta[iele]) < 1.0 ) { smear = 0.020; scale = -0.02; }
	if( fabs(diele.scEta[iele]) > 1.0 && fabs(diele.scEta[iele]) < 1.5 ) { smear = 0.030; scale = -0.045; }
	if( fabs(diele.scEta[iele]) > 1.5 && fabs(diele.scEta[iele]) < 2.5 ) { smear = 0.045; scale = -0.05; }
	
	double c = gen.Gaus(1+scale,smear);
	corr_resp[iele] *= 1./c;
      }
      diele.corrResp(corr_resp);
    }


    ///---- start selection cut
    float    R9min = 0.10;
    //    cout << " eSC 1 = " << diele.eSC[i1] << " eta = " << diele.eta[i1] << endl;
    if( nEvtCut.size() < icut+1 ) {
      nEvtCut.push_back(0);
      nameCut.push_back("- nEvt[ptCuts && |etaSC|< 2.5]: " );
    }
    //    cout << " scEta1: " << diele.scEta[i1] << " - scEta2 : " << diele.scEta[i2] << endl;
    //    cout << " eta1: " << diele.eta[i1] << " - eta2 : " << diele.eta[i2] << endl;
    //    cout << " eng1: " << diele.eSC[i1] << " - eng2 : " << diele.eSC[i2] << endl;
    //    cout << " pt1 : " << diele.eSC[i1]/cosh(diele.eta[i1]) << " - pt2  : " << diele.eSC[i2]/cosh(diele.eta[i2]) << endl;

    if( diele.eSC[i1]/cosh(diele.eta[i1]) < _ptCut || fabs( diele.scEta[i1] ) > 2.5 ) continue;
    if( diele.eSC[i2]/cosh(diele.eta[i2]) < _ptCut || fabs( diele.scEta[i2] ) > 2.5 ) continue;
    nEvtCut[icut++]++;

    //---- fiducial cuts except when computing the scale
    if( nEvtCut.size() < icut+1 ) {
      nEvtCut.push_back(0);
      nameCut.push_back("- nEvt[fiducial cuts]: " );
    }
    if( _analysis != 1 )
      if( ( fabs(diele.scEta[i1]) < 1.566 && fabs(diele.scEta[i1]) > 1.4446 ) ||
	  ( fabs(diele.scEta[i2]) < 1.566 && fabs(diele.scEta[i2]) > 1.4446 ) ) continue;
    nEvtCut[icut++]++;

    
    if( nEvtCut.size() < icut+1 ) {
      nEvtCut.push_back(0);
      nameCut.push_back("- nEvt[r9 > r9min    ]: " );
    }

    if( diele.R9[i1] < R9min || diele.R9[i2] < R9min ) continue;
    nEvtCut[icut++]++;
    if( _doFIT) {
      /// when doing FIT only select EB-EB or EB-EE+EE-EE pairs depending on the fit
      if( !_isEEfit && ( fabs(diele.scEta[i1]) > CCstop || fabs(diele.scEta[i2]) > CCstop )  ) continue;

      /// EE calibration: remove central gap and mod4
      if( _isEEfit ) {
	if(   fabs(diele.scEta[i1]) < 0.05 || fabs(diele.scEta[i2]) < 0.05 )   continue;
	if( ( fabs(diele.scEta[i1]) < CCstop && fabs(diele.scEta[i1]) > 1.38 ) ||
	    ( fabs(diele.scEta[i2]) < CCstop && fabs(diele.scEta[i2]) > 1.38 ) ) continue;
	if( fabs(diele.scEta[i1]) < +CCstop && fabs(diele.scEta[i2]) < +CCstop ) continue;
	if( ecalp == "EEm" && (diele.scEta[i1] > +1.444 || diele.scEta[i2] > +1.444 ) ) continue;
	if( ecalp == "EEp" && (diele.scEta[i1] < -1.444 || diele.scEta[i2] < -1.444 ) ) continue;
      }

      if( nEvtCut.size() < icut+1 ) {
	nEvtCut.push_back(0);
	nameCut.push_back("- nEvt[EB || EE sel ]: " );
      }
      nEvtCut[icut++]++;
    }
     
    if( std::isnan(diele.eSC[i1]) || std::isnan(diele.eSC[i2]) ) continue;
    if( nEvtCut.size() < icut+1 ) {
      nEvtCut.push_back(0);
      nameCut.push_back("- nEvt[NAN eSC check  ]: " );
    }
    nEvtCut[icut++]++;
    
    double mee_corr = diele.mee();
    if(  _testOversmearingMC == 1 )  mee_corr = mee_true;

    if( std::isnan(mee_corr) ) continue;
    if( nEvtCut.size() < icut+1 ) {
      nEvtCut.push_back(0);
      nameCut.push_back("- nEvt[NAN mee check  ]: " );
    }
    nEvtCut[icut++]++;

    if( _doFIT ) {
      //--- add event to fitter 
      //      if( mee_corr < _meeFitMin-0.5 || mee_corr > _meeFitMax+0.5 ) continue;
      if( mee_corr < 50 || mee_corr > 150 ) continue;
      if( nEvtCut.size() < icut+1 ) {
	nEvtCut.push_back(0);
	nameCut.push_back("- nEvt[fit window  ]: " );
      }
      nEvtCut[icut++]++;

      float expReso[2] = {-1,-1};
      for( int iele = 0 ; iele < 2; iele++ ) 
	expReso[iele] = diele.eSC_err[iele]/diele.eSC[iele];
      bool passTestRelReso = true;
      if( _fitRelResoForRegrEnCorr ) 
	for( int iele = 0 ; iele < 2; iele++ ) {
	  if( expReso[iele] <= 0.0001 ) 
	    cout << " expReso[" << iele << "] < 0:  resolution: " << diele.eSC_err[iele]
		 << " / E = " << diele.eSC[iele] << " eta: " << diele.scEta[iele] << " r9: " << diele.R9[iele] << endl;
	  passTestRelReso = false;
	  break;
	}
      if( !passTestRelReso ) continue;
      
      int iresp[2], ireso[2];
      vector<int> itailPar[2]; for( int iele = 0 ; iele < 2; iele++ ) itailPar[iele].resize( _tailAxisND.size() );
      vector<double> xResp[2], xReso[2];
      vector<vector<double>* > xTail(_tailAxisND.size() ) ; 
      for( unsigned itail = 0; itail < _tailAxisND.size(); itail++ ) xTail[itail] = new vector<double>[2];


      diele.fill_Xvar_forBinFinder( & _respAxisND, xResp );
      diele.fill_Xvar_forBinFinder( & _resoAxisND, xReso );
      for( unsigned itail = 0; itail < _tailAxisND.size(); itail++ ) 
	diele.fill_Xvar_forBinFinder( & _tailAxisND[itail], xTail[itail] );

      for( int iele = 0 ; iele < 2; iele++ ) {
	iresp[iele] = _respAxisND.findBin( xResp[iele] );
	ireso[iele] = _resoAxisND.findBin( xReso[iele] ) + _respAxisND.nBinsND();
	
	int iParTail0 = _respAxisND.nBinsND() + _resoAxisND.nBinsND();
	for( unsigned itailp = 0 ; itailp < _tailAxisND.size(); itailp++ ) {
	  itailPar[iele][itailp] = _tailAxisND[itailp].findBin(xTail[itailp][iele]) + iParTail0; 
	  iParTail0 += _tailAxisND[itailp].nBinsND();
	}
      }
      for( unsigned itail = 0; itail < _tailAxisND.size(); itail++ ) delete[] xTail[itail];
	    
      
      bool fixResp[2] = {false,false};
      bool isEB[2] = { true, true };
      for( int iele = 0 ; iele < 2; iele++ )
	if( fabs(diele.scEta[iele]) > CCstop ) isEB[iele] = false;
      if( nEvtCut.size() < icut+1 ) {
	nEvtCut.push_back(0);
	nameCut.push_back("- nEvt[final number ]: " );
      }
      nEvtCut[icut++]++;
      fittedEvents++;
      _zfitterFab->addEvent( mee_corr, iresp, ireso,  itailPar, isEB,
			     diele.raw_eSC , diele.raw_eES ,
			     diele.raw_eES1, diele.raw_eES2,
			     expReso, diele.weight,mee_true,fixResp );

      
    } else if( _outputTree ) {
      diele.fillOutputTree();
      if( nEvtCut.size() < icut+1 ) {
	nEvtCut.push_back(0);
	nameCut.push_back("- nEvt[final number ]: " );
      }
      nEvtCut[icut++]++;
    }

  }

  for( unsigned icut = 0 ; icut < nEvtCut.size(); icut++ ) 
    cout << nameCut[icut] << nEvtCut[icut] << endl;

  cout << "---- IJazZ: # of events in run range    :  " << nTotEvents   << endl;
  cout << "---- IJazZ: # of events selected for fit:  " << fittedEvents << endl;

  if( _isMC ) _hNvtxWeight->Write();
  if( _outputTree ) {
    diele.writeTree();
    cout << "---- IJazZ: output tree saved to: " << _fileout->GetName() << endl;
    _fileout->Close();
  }

}






///--- Fix ME: put me somewhere else
///-------------- smearing mc for testing purpose ----------------///
float mcSmearingTest( float sceta, float r9 ) {
  float r9_bad_eb[]  = { 8.91085e-03, 6.60522e-03, -2.86797e-02, 4.73330e-02, -1.95607e-02 };
  float r9_gold_eb[] = { 7.62684e-03, 1.13788e-02, -4.14171e-02, 5.57636e-02, -1.93949e-02 };
  
  float r9_gold_ee[] = { -4.64302e-01, 9.20859e-01, -5.54852e-01, 1.07274e-01, 0 };
  float r9_bad_ee[]  = { -1.47432e-01, 2.22487e-01, -1.26847e-02,-6.83499e-02, 1.99454e-02 };


  float *par = 0;
  if( fabs(sceta) <  1.5 && r9 >  0.94 ) par = r9_gold_eb;
  if( fabs(sceta) <  1.5 && r9 <= 0.94 ) par = r9_bad_eb;
  if( fabs(sceta) >= 1.5 && r9 >  0.94 ) par = r9_gold_ee;
  if( fabs(sceta) >= 1.5 && r9 <= 0.94 ) par = r9_bad_ee;

  float res = 0;
  for( int ip = 4 ; ip >= 0; ip-- ) res = par[ip] + fabs(sceta)*res; 

  double offset = 0.;
  offset = 0.005*(1.+1./1.5*fabs(sceta));
  if( r9 < 0.94 ) offset += 0.005;

  return res + offset;
}
