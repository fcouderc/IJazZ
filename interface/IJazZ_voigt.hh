#ifndef __IJazZ_voigt
#define __IJazZ_voigt

#include <TMath.h>
#include <vector>
#include <interface/IJazZ_voigt.hh>
#include <map>
#include <iostream>

class paramDef {
public: 
  bool operator==( const paramDef &d ) const { 
    if( fabs(d._gamma - _gamma ) < _equal && 
	fabs(d._sigma - _sigma ) < _equal && 
	fabs(d._mean  - _mean  ) < _equal )
      return true; 
    return false;
  }

  bool operator<=( const paramDef &d ) const { 
    if( fabs( _mean - d._mean ) < _equal  ) {
      if( fabs( _sigma - d._sigma) < _equal ) {
	if( fabs( _gamma - d._gamma) > _equal ) {
	  if( _gamma < d._gamma ) return true;
	  else                    return false;
	} else return true;
      } else {
	if( _sigma < d._sigma ) return true;
	else                    return false;
      } 
    } else {
      if( _mean < d._mean  ) return true;
      else                   return false;
    }
    
    //    std::cout << " ! a <= b " << std:: endl;
    return false;

  }

  bool operator<( const paramDef &d  ) const { 
    //    std::cout << " operator <: "; print();  std:: cout << "           vs "; d.print();

    if( fabs( _mean - d._mean ) < _equal  ) {
      if( fabs( _sigma - d._sigma) < _equal ) {
	if( fabs( _gamma - d._gamma) > _equal ) {
	  if( _gamma < d._gamma ) return true;
	  else                    return false;
	} else return false;
      } else {
	if( _sigma < d._sigma ) return true;
	else                    return false;
      } 
    } else {
      if( _mean < d._mean ) return true;
      else                   return false;
    }
    
    //    std::cout << " ! a < b " << std:: endl;
    return false;
  }

  bool operator>( const paramDef &d  ) const { 
    return !( d <= *const_cast<paramDef*>(this) );
  }

  bool operator>=( const paramDef &d ) const { 
    return !( d < *const_cast<paramDef*>(this) );
  }
  
  paramDef & operator=( const paramDef &d ) { 
    _sigma = d._sigma;
    _gamma = d._gamma;
    _mean  = d._mean ;
    _equal = d._equal;
    return *this;
  }
  
  void print(void) const { std::cout << " mean: " << _mean << " sigma: " << _sigma << " ; gamma: " << _gamma << std::endl; }
  
  paramDef(void);
  paramDef( double mean, double sigma, double gamma );
  paramDef( const paramDef & d ) { *this = d; }
  

  double _gamma;
  double _sigma;
  double _mean;
  
  double _integral;
  double _equal;
};

paramDef::paramDef(void) : _equal(0.000001) {}
paramDef::paramDef( double mean, double sigma, double gamma ): _equal(0.000001) {
  _mean = mean; _sigma = sigma; _gamma = gamma; 
}


class voigtian {
  
public: 
  voigtian(void);
  voigtian(double mean, double sigma, double gamma );


  void   setParams( double mean, double sigma, double gamma );
  double val( double x );
  double primitive( double x );
  double integral(  double x1, double x2 );

  double valLorentz( double x );
  double integralLorentz( double x1, double x2 );

  void disableIntegralStorage( bool disable = true ) { _disableStorage = disable; }
  unsigned nIntegralStored(void) { return _integralList.size(); }
  
private:
  void init(void);
  std::vector<double> _Lv;
  std::vector<double> _Kv;
  std::vector<double> _Mv;
  std::vector<double> _Nv;
  std::vector<double> _Gv;

  double _mean;
  double _sigma;
  double _gamma;
  double _gOs;
  double _integral;

  std::map<paramDef,double> _integralList;
  bool _disableStorage;
};

voigtian::voigtian(void) {init();}
voigtian::voigtian(double mean, double sigma, double gamma ) { setParams(mean,sigma,gamma); init();}



void voigtian::init(void) {
  _Lv.resize(3);
  _Kv.resize(3);
  _Mv.resize(3);
  _Nv.resize(3);
  _Gv.resize(3);

  _Lv[0] = 0        ; _Kv[0] = 0      ; _Mv[0] = 1.32272; _Nv[0] = 0.081905;
  _Lv[1] = 0.090227 ; _Kv[1] = 1.09148; _Mv[1] = 1.29081; _Nv[1] = 0.0093116;
  _Lv[2] = 0.0035776; _Kv[2] = 2.30556; _Mv[2] = 1.17417; _Nv[2] = -0.0116099;

  _integral = 0;
  for( unsigned i = 0 ; i < _Gv.size(); i++ ) _integral += _Nv[i];
  _integral *= TMath::TwoPi();
  _disableStorage = false;
}

void voigtian::setParams(double mean, double sigma, double gamma ) {
  _mean  = mean; 
  _sigma = sigma;
  _gamma = gamma;
  _gOs   = _gamma/(2*_sigma);

  for( unsigned i = 0 ; i < _Gv.size(); i++ )
    _Gv[i] = _gOs + _Mv[i];
}


double voigtian::val( double x ) {
  if( fabs(_sigma) < 0.0001 ) return valLorentz(x);

  double v = (x-_mean)/_sigma;

  double Vx = 0;
  for( unsigned i = 0; i < _Kv.size(); i++ ) {
    double Gi2 = _Gv[i]*_Gv[i];
    double Vpi = (v + _Kv[i]);
    double Vmi = (v - _Kv[i]);
    double Lpi = 1./(Vpi*Vpi+Gi2);    
    double Lmi = 1./(Vmi*Vmi+Gi2);

    Vx += _Nv[i]*_Gv[i]*( Lmi + Lpi);
    if( i > 0 ) Vx += _Lv[i]*(Vpi*Lpi - Vmi*Lmi);
  }
  return Vx / (_integral*_sigma) ;
}


double voigtian::primitive(double x) {
  double v = (x-_mean)/_sigma;

  double Px = 0;
  for( unsigned i = 0; i < _Kv.size(); i++ ) {
    double Gi2 = _Gv[i]*_Gv[i];
    double Vpi = (v + _Kv[i]);
    double Vmi = (v - _Kv[i]);
    double Lpi = 1./(Vpi*Vpi+Gi2);    
    double Lmi = 1./(Vmi*Vmi+Gi2);

    Px += _Nv[i]*(atan(Vpi/_Gv[i])+atan(Vmi/_Gv[i]));
    if( i > 0 ) Px += _Lv[i]/2 * log(Lmi/Lpi);
  }

  return Px / _integral ;
}


double voigtian::integral( double x1, double x2 ) {
  if( _disableStorage ) {
    if( fabs(_sigma) < 0.0001 ) return integralLorentz(x1,x2) ;
    else                        return primitive(x2) - primitive(x1);
  }

  /// auto disable the integral storage if the size is becoming comparable to 
  /// the anticipated number of events in loop
  /// this has to determined in a better way
  if( _integralList.size() > 1000 ) disableIntegralStorage();

  paramDef p(_mean,_sigma,_gamma);
  std::map<paramDef,double>::const_iterator integral = _integralList.find( p );
  if( integral != _integralList.end() ) return integral->second;
  else {
    //   std::cout << " adding integral n = : " << _integralList.size()  << " ... for p: "; p.print();
    double inttmp = 0;
    if( fabs(_sigma) < 0.0001 ) inttmp = integralLorentz(x1,x2) ;
    else                        inttmp = primitive(x2) - primitive(x1);    
    _integralList.insert( std::make_pair( p, inttmp ) );
    
    return inttmp;
  }
  return 1;				   
}


double voigtian::valLorentz( double x ) {
  double gO2 = _gamma/2.;
  double v = (x-_mean)/gO2;
  
  /// this is the correct formula:
  /// gO2 really at the denominator because v is normalised
  return 1./(gO2*TMath::Pi()) * 1./(v*v+1);
}

double voigtian::integralLorentz( double x1, double x2 ) {
  double gO2 = _gamma/2.;
  double v1 = (x1-_mean)/gO2;
  double v2 = (x2-_mean)/gO2;
  
  return (atan(v2)-atan(v1))/TMath::Pi();
}

#endif
