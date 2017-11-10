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

#ifndef globalfit_z_hh_
#define globalfit_z_hh_

#include <vector>

#include <interface/IJazZAxisND.hh>

class TH2F;
class TFile;
class TCanvas;
class TThread;

class ZFitterMinuit2 {
  
 public:
  ZFitterMinuit2();
  ZFitterMinuit2(const ZFitterMinuit2&);

  ~ZFitterMinuit2();
  void setNParams( int nparams );
  

  void setTailParRotation( const std::vector<double> &rotEB, const std::vector<double> & rotEE ) { 
    _tailPar0Rotate_EB = rotEB;  _tailPar0Rotate_EE = rotEE;  }
  
  void setFitterModeClassic(bool classic = true ) { _fitterClassic = classic; }
  void addEvent(float mee, 
		int iresp[2], int ireso[2], 
		std::vector<int> tailP[2],
		bool isEB[2], 
		//float esOverEcal[2],
		float raw_eSC[2], float raw_eES[2],
		float raw_eES1[2], float raw_eES2[2],
		
		float exp_ResE[2],
		float weight,
		float mref,
		bool fixResp[2] );
  std::vector<int> nEvtPerParam( void ) {return _nevtParam;}
  void setFitRange(double mmin,double mmax );
  void setNumCPU( int nCPU = 1 ) { if( nCPU >= 1 ) _numCPU = nCPU; }

  /// evaluate objective function 
  //  double operator()(const std::vector<double>& param) const; 
  double eval(const double *param) ; 
  
  double Up() const { return _fUp; }

  void clearEvents(void);
  
 
  TCanvas* fitCrossCheck( const double * param );
  TCanvas* fitAutoAdjustRange( float &xMin, float &xMax );
 private:
  
  bool _gaussOnly;
  bool _fitterClassic;
 /// DeltaChi2 == N sigma
  double   _fUp;
  double _mmin, _mmax;
  std::vector<float> _mee_evt;
  std::vector<float> _mref_evt;
  std::vector<float> _weight_evt;
  std::vector<float> _esOverEcal1_evt;
  std::vector<float> _esOverEcal2_evt;
  std::vector<float> _rawES1_evt;
  std::vector<float> _rawES2_evt;
  std::vector<float> _rawES11_evt;
  std::vector<float> _rawES12_evt;
  std::vector<float> _rawES21_evt;
  std::vector<float> _rawES22_evt;
  std::vector<float> _rawSC1_evt;
  std::vector<float> _rawSC2_evt;

  std::vector<float> _res1_evt;
  std::vector<float> _res2_evt;
  std::vector<bool> _isEB1_evt;
  std::vector<bool> _isEB2_evt;

  std::vector<int> _iResp1_evt;
  std::vector<int> _iResp2_evt;
  std::vector<int> _iReso1_evt;
  std::vector<int> _iReso2_evt;
  std::vector<std::vector<int> > _iTail1_evt;
  std::vector<std::vector<int> > _iTail2_evt;

  std::vector<int> _fixResp1_evt;
  std::vector<int> _fixResp2_evt;
  std::vector<int> _nevtParam;

  float _mZ;
  float _gammaZ;


  IJazZAxisND<float> _integral;
  std::vector<double> _tailPar0Rotate_EB, _tailPar0Rotate_EE;
  

  struct thread_data_t {
    unsigned ifirst;
    unsigned ilast;
    double llh;
    double n;
    const ZFitterMinuit2 *that;
    const double *param;
  };

  int _numCPU;
  
  void llhCompute( thread_data_t *data ) const;
  static void* llhComputeStatic( void *d );
 };

#endif
