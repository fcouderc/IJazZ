//---------------------------------------------------------------------------//
// Description: loop other events, compute the LLH, this is the function
//               called by the minimizer
//               supports MultiThread processing
//
// Warning: by defaults all the CPU on the machin are used (control via setNumCPU)
//
//
// Author: Fabrice Couderc, fabrice.couderc@cea.fr 
// History: 2012/11/19 - v1
//
//---------------------------------------------------------------------------//

#include "interface/ZFitterMinuit2_ND_MThread.hh"
#include "interface/IJazZ_voigt.hh"

#include <TMath.h>
#include <TStopwatch.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TThread.h>

#include <iostream>
#include <iomanip>
#include <vector>

#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <pthread.h>

#include <Minuit2/FCNBase.h>
using namespace std;

TStopwatch zfitterTimer;

int _niter;


ZFitterMinuit2::ZFitterMinuit2() {
  _mZ     = 91.1876;
  _gammaZ = 2.4952;
  _niter  = 0;

  _mmin = 85;
  _mmax = 100;
  _fUp = 1.0;
  
  _numCPU = sysconf( _SC_NPROCESSORS_ONLN );
  if( _numCPU <= 0 ) _numCPU = 1;
  _gaussOnly = false;
  _fitterClassic = false;
}




ZFitterMinuit2::ZFitterMinuit2(const ZFitterMinuit2& es) {
  _mZ     = 91.1876;
  _gammaZ = 2.4952;
  
  _fUp = es._fUp;
  
  _mee_evt    = es. _mee_evt;
  _mref_evt   = es. _mref_evt;
  _iResp1_evt = es._iResp1_evt;
  _iResp2_evt = es._iResp2_evt;
  _iReso1_evt = es._iReso1_evt;
  _iReso2_evt = es._iReso2_evt;
  _iTail1_evt = es._iTail1_evt;
  _iTail2_evt = es._iTail2_evt;
  _numCPU     = es._numCPU;
  _weight_evt = es. _weight_evt;
  _rawSC1_evt = es._rawSC1_evt;
  _rawSC2_evt = es._rawSC2_evt;
  _rawES1_evt = es._rawES1_evt;
  _rawES2_evt = es._rawES2_evt;
  _rawES11_evt = es._rawES11_evt;
  _rawES12_evt = es._rawES12_evt;
  _rawES21_evt = es._rawES21_evt;
  _rawES22_evt = es._rawES22_evt;


  _fitterClassic = es._fitterClassic;
  
 }

void ZFitterMinuit2::setNParams( int npars ) { 
  _nevtParam.resize(npars,0);
}

void  ZFitterMinuit2::clearEvents(void) {
  _mee_evt.clear();
  _mref_evt.clear();
  _weight_evt.clear();

  _esOverEcal1_evt.clear();
  _esOverEcal2_evt.clear();
  _iResp1_evt.clear();
  _iResp2_evt.clear();
  _iReso1_evt.clear();
  _iReso2_evt.clear();
  _iTail1_evt.clear();
  _iTail2_evt.clear();
  _isEB1_evt.clear();
  _isEB2_evt.clear();

  _rawES1_evt.clear();
  _rawES2_evt.clear();
  _rawES11_evt.clear();
  _rawES12_evt.clear();
  _rawES21_evt.clear();
  _rawES22_evt.clear();
  _rawSC1_evt.clear();
  _rawSC2_evt.clear();

  _fixResp1_evt.clear();
  _fixResp2_evt.clear();
  _res1_evt.clear();
  _res2_evt.clear();
  _nevtParam.clear();
}

ZFitterMinuit2::~ZFitterMinuit2() {
  clearEvents();
}

void ZFitterMinuit2::addEvent(float mee, int iresp[2], int ireso[2],
			      vector<int>  tailP[2],
			      bool isEB[2],
			      float raw_eSC[2] , float raw_eES[2] ,
			      float raw_eES1[2], float raw_eES2[2],
			      float exp_ResE[2],
			      float weight, float mref,
			      bool fixResp[2]) {

  if( _nevtParam.size() == 0 ) {
    cerr << "ZFitterMinuit2::addEvent: bins number not set, I will crash..." << endl;
    return;
  }


  _mee_evt.push_back(mee/_mZ);
  _mref_evt.push_back(mref);
  _weight_evt.push_back(weight);
  //  _esOverEcal1_evt.push_back(raw_eES[0]/raw_eSC[0]);
  //  _esOverEcal2_evt.push_back(raw_eES[1]/raw_eSC[1]);
  
  _esOverEcal1_evt.push_back( 0 );
  _esOverEcal2_evt.push_back( 0 );

  _rawSC1_evt.push_back(raw_eSC[0]);
  _rawSC2_evt.push_back(raw_eSC[1]);
  _rawES1_evt.push_back(raw_eES[0]);
  _rawES2_evt.push_back(raw_eES[1]);
  _rawES11_evt.push_back(raw_eES1[0]);
  _rawES12_evt.push_back(raw_eES1[1]);
  _rawES21_evt.push_back(raw_eES2[0]);
  _rawES22_evt.push_back(raw_eES2[1]);
  
  _iResp1_evt.push_back(iresp[0]);
  _iResp2_evt.push_back(iresp[1]);  
  _iReso1_evt.push_back(ireso[0]);
  _iReso2_evt.push_back(ireso[1]);

  
  _iTail1_evt.push_back( tailP[0] );
  _iTail2_evt.push_back( tailP[1] );

  _isEB1_evt.push_back(isEB[0]);
  _isEB2_evt.push_back(isEB[1]);

  for( int iele = 0 ; iele < 2; iele++ ) {
    _nevtParam[iresp[iele]]++;
    _nevtParam[ireso[iele]]++;
    if( ! _fitterClassic ) {
    for( unsigned itail = 0 ;itail < tailP[0].size(); itail++ )
      _nevtParam[ tailP[iele][itail] ]++;
    }
  }
 
  _res1_evt.push_back(exp_ResE[0]);
  _res2_evt.push_back(exp_ResE[1]);

  _fixResp1_evt.push_back( fixResp[0] );
  _fixResp2_evt.push_back( fixResp[1] );
}


void ZFitterMinuit2::setFitRange(double mmin,double mmax ) {
  _mmin = mmin;
  _mmax = mmax;
}


void ZFitterMinuit2::llhCompute( thread_data_t *data ) const {
  double LLH    = 0;

  float nAccept = 0;
  voigtian ijazzVoigt;

  int nContrib = 4;
  if( _fitterClassic ) nContrib = 1;

  double r1(1), r2(1), s1(0.015), s2(0.015), y1(1), y2(1);
  //  double rES(1);
  //  double rES1(1), rES2(1);
  double f1(1), f2(1), d1(0), d2(0), k1(1), k2(1),m1(0),m2(0);
  //  double exp1(1), exp2(1);
  //  double mref(_mZ);
  double integral(1);
  double relGammaZ = _gammaZ / _mZ;
  double mmin = _mmin / _mZ;
  double mmax = _mmax / _mZ;
  double integrals[4], scales[4], frac[4], sigmas[4];
  double gammas[4];
  for( unsigned ievt = data->ifirst ; ievt < data->ilast; ++ievt ) {
    if(  _mee_evt[ievt] < mmin || _mee_evt[ievt] > mmax )  continue;  
    r1 = data->param[_iResp1_evt[ievt]];
    r2 = data->param[_iResp2_evt[ievt]];
    y1 = data->param[_iReso1_evt[ievt]];
    y2 = data->param[_iReso2_evt[ievt]];
    if( !_fitterClassic ) {
      m1 = data->param[_iTail1_evt[ievt][0]];
      m2 = data->param[_iTail2_evt[ievt][0]];
      s1 = data->param[_iTail1_evt[ievt][1]];
      s2 = data->param[_iTail2_evt[ievt][1]];
      d1 = data->param[_iTail1_evt[ievt][2]];
      d2 = data->param[_iTail2_evt[ievt][2]];
      k1 = data->param[_iTail1_evt[ievt][3]];
      k2 = data->param[_iTail2_evt[ievt][3]];
      f1 = data->param[_iTail1_evt[ievt][4]];
      f2 = data->param[_iTail2_evt[ievt][4]];
      double k1dmu =  data->param[_iTail1_evt[ievt][5]]; 
      double k2dmu =  data->param[_iTail2_evt[ievt][5]];
      d1 *= k1dmu;
      d2 *= k2dmu;
      double k1rsig =  data->param[_iTail1_evt[ievt][6]]; 
      double k2rsig =  data->param[_iTail2_evt[ievt][6]];
      f1 *= k1rsig; if( f1 >= 1 ) f1 = 1;
      f2 *= k2rsig; if( f2 >= 1 ) f2 = 1;
    } else {
      s1 = y1;
      s2 = y2;
    }
    
    if( !_fitterClassic ) {
      r1 = 1 + (r1+m1)*0.01; 
      r2 = 1 + (r2+m2)*0.01; 
      d1 = -d1*s1*0.01;
      d2 = -d2*s2*0.01;
    } else {
      r1 = r1*0.01+1;
      r2 = r2*0.01+1;
      //      rES1 = data->param[_nevtParam.size()-2];
      //      rES2 = data->param[_nevtParam.size()-1];
      //r1 = (r1*_rawSC1_evt[ievt]+rES1*_rawES11_evt[ievt]+rES2*_rawES21_evt[ievt])/(_rawSC1_evt[ievt]+_rawES11_evt[ievt]+_rawES21_evt[ievt]);
      //r2 = (r2*_rawSC2_evt[ievt]+rES1*_rawES12_evt[ievt]+rES2*_rawES22_evt[ievt])/(_rawSC2_evt[ievt]+_rawES12_evt[ievt]+_rawES22_evt[ievt]);

      /// potentially for eta scale correct for ES
      r1 = (_esOverEcal1_evt[ievt] + 1 ) / (_esOverEcal1_evt[ievt] + 1./r1 );
      r2 = (_esOverEcal2_evt[ievt] + 1 ) / (_esOverEcal2_evt[ievt] + 1./r2 );
    }

    double s1M_2 =  s1*s1;
    double s2M_2 =  s2*s2; 
    scales[0] = sqrt( r1    * r2) ; frac[0] = 1; 
    sigmas[0] = 0.5 * sqrt( s1M_2 + s2M_2 )*0.01; gammas[0] = scales[0]*relGammaZ;    
    if( !_fitterClassic ) {
      frac[0] =    f1 *   f2 ;
      scales[1] = sqrt( r1    * (r2+d2)); frac[1] =    f1 *(1-f2);
      scales[2] = sqrt((r1+d1)* r2)     ; frac[2] = (1-f1)*   f2 ;
      scales[3] = sqrt((r1+d1)* (r2+d2)); frac[3] = (1-f1)*(1-f2);
      if( _gaussOnly ) {
	double s1W_2 = k1*k1*s1*s1;
	double s2W_2 = k2*k2*s2*s2;
	double addOS = y1+y2 ;
	sigmas[0] = addOS > -s1M_2-s2M_2 ? 0.5 * sqrt( s1M_2 + s2M_2 + addOS )*0.01 : 0;
	sigmas[1] = addOS > -s1W_2-s2M_2 ? 0.5 * sqrt( s1M_2 + s2W_2 + addOS )*0.01 : 0; 
	sigmas[2] = addOS > -s1M_2-s2W_2 ? 0.5 * sqrt( s1W_2 + s2M_2 + addOS )*0.01 : 0; 
	sigmas[3] = addOS > -s1W_2-s2W_2 ? 0.5 * sqrt( s1W_2 + s2W_2 + addOS )*0.01 : 0; 
	gammas[0] = scales[0]*relGammaZ;
	gammas[1] = scales[1]*relGammaZ;
	gammas[2] = scales[2]*relGammaZ;
	gammas[3] = scales[3]*relGammaZ;
      } else {
	double addOS = y1+y2 ;
	sigmas[0] = s1M_2 + s2M_2 + addOS > 0 ? 0.5 * sqrt( s1M_2 + s2M_2 + addOS )*0.01 : 0; 
	sigmas[1] = s1M_2 + addOS > 0 ? 0.5 * sqrt( s1M_2 + addOS )*0.01 : 0; 
	sigmas[2] = s2M_2 + addOS > 0 ? 0.5 * sqrt( s2M_2 + addOS )*0.01 : 0; 
	sigmas[3] = 0     + addOS > 0 ? 0.5 * sqrt( addOS )*0.01 : 0;
	gammas[0] = scales[0]*relGammaZ;    
	gammas[1] = scales[1]*relGammaZ + k2*s2*0.5*0.01;
	gammas[2] = scales[2]*relGammaZ + k1*s1*0.5*0.01;
	gammas[3] = scales[3]*relGammaZ + (k1*s1*+k2*s2)*0.5*0.01;
      }
    }

    double mu_mee = 0;
    integral = 0;
    for( int i = 0 ; i < nContrib; i++ ) {
      ijazzVoigt.setParams( scales[i], sigmas[i], gammas[i]);
      integrals[i] = ijazzVoigt.integral( mmin, mmax );
      integral += integrals[i] * frac[i];
    }

    for( int i = 0 ; i < nContrib; i++ ) {
      ijazzVoigt.setParams( scales[i], sigmas[i], gammas[i]);
      mu_mee   += ijazzVoigt.val( _mee_evt[ievt] )* frac[i];
    }
    mu_mee /= integral;

    for( int i = 0 ; i < nContrib; i++ )
    if( integrals[i] <= 0.00000001 || integrals[i] > 1 ) 
      cout << " integral [ i = " << i << "] is wrong: " << integrals[i] << endl
    	   << " mass scale: " << sqrt(r1*r2) << " mass resolution: " << sigmas[0] << endl
    	   << "   ==> EB1: " << _isEB1_evt[ievt] << "; EB2: " << _isEB2_evt[ievt] << endl
    	   << "   ==> r1: " << r1 << "; r2: " << r2 << endl
	//    	   << "   ==> d1: " << d1 << "; d2: " << d2 << endl
	//    	   << "   ==> d1(before): " << d1 << "; d2: " << d2 << endl
	//    	   << "   ==> f1: " << f1 << "; f2: " << f2 << endl
	   << "   ==> s1: " << s1 << "; s2: " << s2 << endl
	//    	   << "   ==> k1: " << k1 << "; k2: " << k2 << endl
    	   << "   ==> mee: " << _mee_evt[ievt] << " for ievt = " << ievt << endl
    	   << "   ==> iresp1: " << _iResp1_evt[ievt] << endl
    	   << "   ==> iresp2: " << _iResp2_evt[ievt] << endl;

    LLH     += log(mu_mee) *_weight_evt[ievt];
    nAccept += _weight_evt[ievt];

  }

  data->llh = LLH;
  data->n   = nAccept;
}


void* ZFitterMinuit2::llhComputeStatic( void *d )  {
  thread_data_t *data = (thread_data_t*) d;
  data->that->llhCompute(data);
  return 0;
}
  

double ZFitterMinuit2::eval(const double *  param)   {
  //  if( _vIntegrals[0] == 0 || _vIntegrals[1] == 0 ) cout << " --- I will crash... setup the mass range befre." << endl;  
  if( _niter == 0      ) zfitterTimer.Reset();

  zfitterTimer.Start(0);
  //------------------------------------------------------------------------------------------------//
  //----------  sending the LLH computation over several CPUs
  double LLH = 0;
  vector<thread_data_t> threadData; threadData.resize( _numCPU );
  vector<pthread_t>       mtp; mtp.resize(_numCPU);
  
  unsigned nEvt = _mee_evt.size()/_numCPU;
  for( int icpu = 0 ; icpu < _numCPU; icpu++ ) {
    threadData[icpu].ifirst = 0    + nEvt*icpu;
    threadData[icpu].ilast  = nEvt + nEvt*icpu;
    threadData[icpu].param  = param;
    threadData[icpu].that   = this;
  }
  threadData[_numCPU-1].ilast  = _mee_evt.size();

  for( int icpu = 0 ; icpu < _numCPU; icpu++ ) {
    int rc = pthread_create( &mtp[icpu], 0, llhComputeStatic, (void*) &threadData[icpu] );
    assert( rc == 0 );
  }

  for( int icpu = 0 ; icpu < _numCPU; icpu++ ) pthread_join(mtp[icpu],NULL);
  for( int icpu = 0 ; icpu < _numCPU; icpu++ ) LLH += threadData[icpu].llh;  
  //------------------------------------------------------------------------------------------------//
  zfitterTimer.Stop();
  
  if( _niter == 10 || _niter%1000 == 0 ) {    
    cout << "**** nevts = " << _mee_evt.size() 
	 << "  ( FCN = " << std::setprecision(10) <<  -2*LLH << " ) - iter: " << _niter
	 << "  ; cpu time calc / iter: " << zfitterTimer.CpuTime() / _niter
      //       << "  ; Total CPU  time: " << zfitterTimer.CpuTime() 
	 << "  ; Total Real time: " << zfitterTimer.RealTime() 
      //       << "   - nCPU = " << _numCPU
	 << "\r" << flush;
    //    getchar();
  }  
  _niter++;
  
  return -2*(LLH);
}





#include <TCanvas.h>
#include <TProfile.h>
#include <TF1.h>
#include <TH1.h>
TCanvas * ZFitterMinuit2::fitCrossCheck( const double* param ) {
  double mmin = _mmin / _mZ; //0.85
  double mmax = _mmax / _mZ; //1.15
  int nbin = 150;
  TProfile *hmass_fitT    =   new TProfile("hmass_fitT" ,"hmass",nbin,mmin,mmax);
  TProfile *hmass_fitW[]  = { new TProfile("hmass_fitW1","hmass",nbin,mmin,mmax),
			      new TProfile("hmass_fitW2","hmass",nbin,mmin,mmax),
			      new TProfile("hmass_fitW3","hmass",nbin,mmin,mmax) };

  TProfile *hmass_bkg   = new TProfile("hmass_bkg","hmass",nbin,mmin,mmax);
  TH1F     *hmass_data  = new TH1F(   "hmass_data","hmass",nbin,mmin,mmax);
  //  TH1F     *hmass_fitH  = new TH1F(   "hmass_fitH","hmass",nbin,mmin,mmax);

  hmass_data->SetLineColor(kBlack);
  hmass_data->SetLineWidth(2);
  hmass_data->GetXaxis()->SetTitle( "M_{ee} [GeV]" );
  hmass_data->GetXaxis()->CenterTitle();

  hmass_fitT->SetLineColor(kAzure);
  //  hmass_fitWa->SetLineColor(kRed+3);
  //  hmass_fitWa->SetLineStyle(kDashed);
  hmass_fitW[0]->SetLineColor(kGreen+1);
  hmass_fitW[0]->SetLineStyle(kDashed);
  hmass_fitW[1]->SetLineColor(kYellow+2);
  hmass_fitW[1]->SetLineStyle(kDashed);
  hmass_fitW[2]->SetLineColor(kCyan+3);
  hmass_fitW[2]->SetLineStyle(kDashed);
  hmass_fitT->SetLineWidth(2);
  //  hmass_fitWa->SetLineWidth(2);
  for( int i = 0 ; i < 3; i++ ) hmass_fitW[i]->SetLineWidth(2);
  hmass_bkg->SetLineColor(kRed);

  voigtian ijazzVoigt;
  int nContrib = 4;
  if( _fitterClassic ) nContrib = 1;

  double integral(1);
  double relGammaZ = _gammaZ / _mZ;
  double integrals[4], scales[4], frac[4], sigmas[4],gammas[4];
  double r1(1), r2(1), f1(1), f2(0), d1(0), d2(0), k1(1), k2(1), s1(1), s2(1), y1(1), y2(1);
  //  double rES(1);
  //  double rES1(1), rES2(1);
  double m1(0), m2(0);
  /*
  _tailAxisND[0] = _respAxisND;
  _tailAxisND[1] = _resoAxisND;
  _tailAxisND[2] = _resoAxisND;
  _tailAxisND[3] = _resoAxisND;
  _tailAxisND[4] = _resoAxisND;

  for( int ipt = 0 ; ipt <  5; ipt++ )
  for( int ib  = 0 ; ib  < _tailAxisND[ipt].nBinsND(); ib++ )  _tailAxisND[ipt].value(ib) = 0;
  */

  double llh = 0;
  cout << " NEVT in xCheck: " <<  _mee_evt.size() << endl;
  for( unsigned ievt = 0 ; ievt < _mee_evt.size(); ++ievt ) {
    if(  _mee_evt[ievt] < mmin || _mee_evt[ievt] > mmax )  continue;

    r1 = param[_iResp1_evt[ievt]];
    r2 = param[_iResp2_evt[ievt]];
    y1 = param[_iReso1_evt[ievt]];
    y2 = param[_iReso2_evt[ievt]];
    if( !_fitterClassic ) {
      m1 = param[_iTail1_evt[ievt][0]];
      m2 = param[_iTail1_evt[ievt][0]];
      s1 = param[_iTail1_evt[ievt][1]];
      s2 = param[_iTail2_evt[ievt][1]];
      d1 = param[_iTail1_evt[ievt][2]];
      d2 = param[_iTail2_evt[ievt][2]];
      k1 = param[_iTail1_evt[ievt][3]];
      k2 = param[_iTail2_evt[ievt][3]];
      f1 = param[_iTail1_evt[ievt][4]];
      f2 = param[_iTail2_evt[ievt][4]];
      double k1dmu = param[_iTail1_evt[ievt][5]]; 
      double k2dmu = param[_iTail2_evt[ievt][5]];
      d1 *= k1dmu;
      d2 *= k2dmu;

    } else {
      s1 = y1;
      s2 = y2;
    }

    if( !_fitterClassic ) {
      r1 = 1+(r1+m1)*0.01; 
      r2 = 1+(r2+m2)*0.01; 
      d1 = -d1*s1*0.01;
      d2 = -d2*s2*0.01;
    } else {
      r1 = r1*0.01+1;
      r2 = r2*0.01+1;
      //      rES = param[_nevtParam.size()-2];	    
      //      r1 = (r1*_rawSC1_evt[ievt]+rES*_rawES1_evt[ievt])/(_rawSC1_evt[ievt]+_rawES1_evt[ievt]);
      //      r2 = (r1*_rawSC2_evt[ievt]+rES*_rawES2_evt[ievt])/(_rawSC2_evt[ievt]+_rawES2_evt[ievt]);
      
      //      rES1 = param[_nevtParam.size()-2];	    
      //      rES2 = param[_nevtParam.size()-1];
      //r1 = (r1*_rawSC1_evt[ievt]+rES1*_rawES11_evt[ievt]+rES2*_rawES21_evt[ievt])/(_rawSC1_evt[ievt]+_rawES11_evt[ievt]+_rawES21_evt[ievt]);
      //r2 = (r2*_rawSC2_evt[ievt]+rES1*_rawES12_evt[ievt]+rES2*_rawES22_evt[ievt])/(_rawSC2_evt[ievt]+_rawES12_evt[ievt]+_rawES22_evt[ievt]);

      // /// potentially for eta scale correct for ES
      r1 = (_esOverEcal1_evt[ievt] + 1 ) / (_esOverEcal1_evt[ievt] + 1./r1 );
      r2 = (_esOverEcal2_evt[ievt] + 1 ) / (_esOverEcal2_evt[ievt] + 1./r2 );
    }
    
    double s1M_2 =  s1*s1;
    double s2M_2 =  s2*s2; 
    scales[0] = sqrt( r1    * r2) ; frac[0] = 1; 
    sigmas[0] = 0.5 * sqrt( s1M_2 + s2M_2 )*0.01; 
    gammas[0] = scales[0]*relGammaZ;    
    if( !_fitterClassic ) {
      scales[0] = sqrt( r1    * r2)     ; frac[0] =    f1 *   f2 ; 
      scales[1] = sqrt( r1    * (r2+d2)); frac[1] =    f1 *(1-f2);
      scales[2] = sqrt((r1+d1)* r2)     ; frac[2] = (1-f1)*   f2 ;
      scales[3] = sqrt((r1+d1)* (r2+d2)); frac[3] = (1-f1)*(1-f2);

      if( _gaussOnly ) {
	double s1W_2 = k1*k1*s1*s1;
	double s2W_2 = k2*k2*s2*s2;
	double addOS = y1+y2 ;
	sigmas[0] = addOS > -s1M_2-s2M_2 ? 0.5 * sqrt( s1M_2 + s2M_2 + addOS )*0.01 : 0;
	sigmas[1] = addOS > -s1W_2-s2M_2 ? 0.5 * sqrt( s1M_2 + s2W_2 + addOS )*0.01 : 0; 
	sigmas[2] = addOS > -s1M_2-s2W_2 ? 0.5 * sqrt( s1W_2 + s2M_2 + addOS )*0.01 : 0; 
	sigmas[3] = addOS > -s1W_2-s2W_2 ? 0.5 * sqrt( s1W_2 + s2W_2 + addOS )*0.01 : 0; 
	gammas[1] = scales[1]*relGammaZ;
	gammas[2] = scales[2]*relGammaZ;
	gammas[3] = scales[3]*relGammaZ;
	
      } else {
	double addOS = y1+y2 ;
	sigmas[0] = s1M_2 + s2M_2 + addOS > 0 ? 0.5 * sqrt( s1M_2 + s2M_2 + addOS )*0.01 : 0; 
	sigmas[1] = s1M_2 + addOS > 0 ? 0.5 * sqrt( s1M_2 + addOS )*0.01 : 0; 
	sigmas[2] = s2M_2 + addOS > 0 ? 0.5 * sqrt( s2M_2 + addOS )*0.01 : 0; 
	sigmas[3] = 0     + addOS > 0 ? 0.5 * sqrt( addOS )*0.01 : 0;
	gammas[0] = scales[0]*relGammaZ;    
	gammas[1] = scales[1]*relGammaZ + k2*s2*0.5*0.01;
	gammas[2] = scales[2]*relGammaZ + k1*s1*0.5*0.01;
	gammas[3] = scales[3]*relGammaZ + (k1*s1*+k2*s2)*0.5*0.01;
      }
    }

    integral = 0;
    for( int i = 0 ; i < nContrib; i++ ) {
      ijazzVoigt.setParams( scales[i], sigmas[i], gammas[i]);
      integrals[i] = ijazzVoigt.integral( mmin, mmax );
      integral += integrals[i] * frac[i];
    }                    
    
    hmass_data ->Fill( _mee_evt[ievt],         _weight_evt[ievt]);
    for( int ib = 1; ib <= hmass_data->GetXaxis()->GetNbins(); ib++ ) {      
      double mee_b = hmass_data->GetXaxis()->GetBinCenter(ib);
      for( int i = 0 ; i < nContrib; i++ ) {
	ijazzVoigt.setParams( scales[i], sigmas[i], gammas[i]);
	double c = ijazzVoigt.val( mee_b ) / integrals[i]*frac[i];  //// check formula when there is more than 1  contribution
	hmass_fitT ->Fill( mee_b ,  c*_weight_evt[ievt]);
	if( i > 0 ) hmass_fitW[i-1]->Fill(mee_b ,  c*_weight_evt[ievt]);
      }
    }
    //    hmass_fitWa->Fill( _mee_evt[ievt], mu_wide,_weight_evt[ievt]);

    //    llh += log(mu_mee) *_weight_evt[ievt];
  }


  TF1 *fTest = new TF1("ff","[0]*TMath::Voigt(x-[1],[2]/91.1876,2.4952/91.1876)",mmin,mmax);
  fTest->SetParameters(10000,1,0.02);
  cout << " integral fit = " << hmass_fitT->Integral() << endl;
  float scale = hmass_data->Integral() / hmass_fitT->Integral();
  hmass_fitT->Scale( scale );
  //  hmass_fitWa->Scale( scale );

  for( int i = 1 ; i < 4; i++ )   hmass_fitW[i-1]->Scale( scale );//hmass_bkg->Scale( scale );



  TCanvas *c = new TCanvas("canCrossCheckFit","fit cross check",800,600);
  c->cd();
  hmass_data->Fit(fTest,"mlhe r 0 q");
  //  hmass_data->Fit(f2,"mlhe r 0");
  hmass_data->DrawCopy("e");
  fTest->Draw("same");
  // f2->Draw("same");
  for( int i = 1 ; i < nContrib; i++ )  hmass_fitW[i-1]->DrawCopy("hist same");
  //  hmass_fitWa->DrawCopy("hist same");
  hmass_fitT->DrawCopy("hist same");
  hmass_data->DrawCopy("e same");
  cout << "   ZFitterMinuit2: cross check with simple 1D binned fit: " << endl;
  cout << "        - ROOT Fit ener scale: " << fTest->GetParameter(1) << endl;
  cout << "        - ROOT Fit resolution: " << fTest->GetParameter(2) / _mZ *sqrt(2)<< endl;

  double chi2FitND = 0;
  int ndof = 0;
  for( int ib  = 1 ; ib <= hmass_data->GetXaxis()->GetNbins(); ib++ ) {
    double vD = hmass_data->GetBinContent(ib);
    double eD = hmass_data->GetBinError(ib);
    double vFitND =  hmass_fitT->GetBinContent(ib);
    double chi = (vD - vFitND);
    if( eD > 0.000001 ) { 
      chi2FitND += chi*chi / (eD*eD);
      ndof += 1;
    }
  }
  cout << " chi2[ND] = " << chi2FitND << " / ndof = " << ndof << " --  [-2llh] = " << -2*llh <<  endl;
  
  c->Update();
  delete hmass_fitT;
  //  delete hmass_fitH;
  //  delete hmass_fitWa;
  for( int i = 0; i < 3; i++ ) delete hmass_fitW[i];
  delete hmass_bkg;
  delete hmass_data;
  delete fTest;
  return c;
}




TCanvas * ZFitterMinuit2::fitAutoAdjustRange( float &xMin, float &xMax ) {
  double mmin =  70; //0.85
  double mmax = 120; //1.15
  int nbin = 150;
  TH1F     *hmass_data  = new TH1F(   "hmass_data_autoAdjust","hmass adjusting fit range",nbin,mmin,mmax);
  hmass_data->SetLineColor(kBlack);
  hmass_data->SetLineWidth(2);
  hmass_data->GetXaxis()->SetTitle( "M_{ee} [GeV]" );
  hmass_data->GetXaxis()->CenterTitle();

  
  int nEvtAccept = 0;
  for( unsigned ievt = 0 ; ievt < _mee_evt.size(); ++ievt ) {
    if(  _mee_evt[ievt]*_mZ < mmin || _mee_evt[ievt]*_mZ > mmax )  continue;
    hmass_data ->Fill( _mee_evt[ievt]*_mZ,         _weight_evt[ievt]);
    nEvtAccept++;
  }
  cout << " Total nEvts: " << _mee_evt.size() << " Accept: " << nEvtAccept << endl;
  double max = hmass_data->GetXaxis()->GetBinCenter(hmass_data->GetMaximumBin());
  double rms = hmass_data->GetRMS();


  TF1 *fTest = new TF1("ff","[0]*TMath::Voigt(x-[1],[2],2.4952)",max-0.5*rms,max+1.5*rms);
  fTest->SetParameters(10000,max,0.8*rms);
  fTest->SetParLimits(1,max-2*rms,max+2*rms);
  fTest->SetParLimits(2,0.05*rms,1.5*rms);

  fTest->SetLineColor(kRed);
  TCanvas *c = new TCanvas("canvasFitAutoAdjust","fit autoadjust range",800,600);
  c->cd();
  hmass_data->Fit( fTest, "mlhe r 0");
  hmass_data->DrawCopy("e");
  fTest->Draw("same");
  hmass_data->DrawCopy("e same");
  c->Update();


  xMin = max-1.0*rms;
  xMax = max+2.0*rms;

  double sigG   =fTest->GetParameter(2);
  double fwhmG = sigG*sqrt(log(2)*2);
  double fwhmL = 2.4952; 
  double phiLG = fwhmL/fwhmG;
  double c0 = 2.0056; double c1 = 1.0593;
  double fwhmV = fwhmG*(1-c0*c1+sqrt(phiLG*phiLG+2*c1*phiLG+c0*c0*c1*c1));
  double rmsVapprox = fwhmV / sqrt(log(2)*2);
  double mean       = fTest->GetParameter(1);
  xMin = mean - 1*rmsVapprox;
  xMax = mean + 3*rmsVapprox;
  cout << "   ZFitterMinuit2: cross check with simple 1D binned fit: " << endl;
  cout << "        - ROOT Fit ener scale: " << fTest->GetParameter(1) << endl;
  cout << "        - ROOT Fit resolution: " << fTest->GetParameter(2) << endl;
  cout << "        - x @ maximum = " << max  << " rms = " << rms << endl;
  cout << "        - x @ maximum = " << mean << " rms = " << rmsVapprox << endl;
  cout << "        - Deduce fit range:  " << xMin << " --> " << xMax << "" << endl;


  delete hmass_data;
  ///  delete fTest;
  return c;
}
