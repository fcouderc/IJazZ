//---------------------------------------------------------------------------//
// Description: Transform IJazZAxisND histogram in a TGraphErrors ROOT format
//              to visualize the result or into a TH2F for 2D visualisation
//
//
// Author: Fabrice Couderc, fabrice.couderc@cea.fr 
// History:
//
//---------------------------------------------------------------------------//

#ifndef fabAxisViewer__hh__
#define fabAxisViewer__hh__

#include "interface/IJazZAxis.hh"
#include "interface/IJazZAxisND.hh"
#include "interface/RootUtils.hh"

#include <TGraphErrors.h>
#include <TH2F.h>
#include <assert.h>


template <typename T>
class IJazZAxisViewer {
public:
  IJazZAxisViewer(void) {};
  IJazZAxisViewer(const IJazZAxisND<T> &axis)       { _axis = axis;}
  IJazZAxisViewer(const IJazZAxisViewer<T> &viewer) { _axis = viewer._axis;}
  IJazZAxisViewer(std::string axisFile)           { _axis.readFromFile(axisFile); }
  IJazZAxisViewer& operator=(const IJazZAxisViewer<T> &viewer)   { _axis = viewer._axis; return *this;}

  ~IJazZAxisViewer(void) {}
  IJazZAxisND<T>* getAxisND() { return &_axis; }
  TGraphErrors *proj1D( void );
  TGraphErrors *proj1D( unsigned dim, T x );
  TGraphErrors *proj1D( unsigned dim, T x, T y );
  TGraphErrors *proj1D( unsigned dim, const std::vector<T> &xMinus1 = std::vector<T>() );
  TH1F *hist1D( unsigned dim, T x );
  TH1F *hist1D( unsigned dim, T x, T y );
  TH1F *hist1D( unsigned dim, const std::vector<T> &xMinus1 = std::vector<T>() );
  TH2F *proj2D( unsigned d1 = 0, unsigned d2 = 1, const std::vector<T> &xMinus2 = std::vector<T>() );
  TH2F *proj2D( unsigned d1, unsigned d2, T z );
  TH2F *proj2D( unsigned d1, unsigned d2, T z, T t );
  
private:
  IJazZAxisND<T> _axis;
  
};

template <typename T>
TGraphErrors *IJazZAxisViewer<T>::proj1D( void ) {
  assert( _axis.getND() == 1 );
  typename std::vector<T> xMinus1; 
  return proj1D(0,xMinus1);
}

template <typename T>
TGraphErrors *IJazZAxisViewer<T>::proj1D( unsigned dim, T x ) {
  assert( _axis.getND() == 2 );
  typename std::vector<T> xMinus1; xMinus1.push_back(x);
  return proj1D(dim,xMinus1);
}


template <typename T>
TGraphErrors *IJazZAxisViewer<T>::proj1D( unsigned dim, T x, T y ) {
  assert( _axis.getND() == 3 );
  typename std::vector<T> xMinus1; 
  xMinus1.push_back(x);
  xMinus1.push_back(y);
  return proj1D(dim,xMinus1);
}

template <typename T>
TGraphErrors *IJazZAxisViewer<T>::proj1D( unsigned dim, const std::vector<T> &xMinus1 ) {
  IJazZAxis<T> *bX = _axis.getAxis(dim);
  TGraphErrors *out = new TGraphErrors( bX->nBins() );
  int ip = 0;
  for( bX->itBegin(); ! bX->itEnd(); bX->itForward() ) {
    double v  = 0;
    double ev2 = 0;
    for( int ib = 0 ; ib < _axis.nBinsND(); ib++ ) {
      bool addBin = false;
      if( fabs( _axis.getBinCenterDimN(ib,dim) - bX->itBinCenter() ) < 0.0001 ) addBin = true;
      if( addBin && xMinus1.size() == _axis.getND() - 1 ) {
	for( unsigned d = 0 ; d < _axis.getND(); d++ ) {
	  if(  d != dim && fabs( _axis.getBinCenterDimN(ib,d) - _axis.getAxis(d)->binCenter(xMinus1[d<dim?d:d-1]) ) >  0.0001 ) addBin = false;	  
	  if( !addBin ) break;
	}
      }
      if( addBin ) {
	double err2 = 999999;
	if( _axis.error(ib) > 0.00001 ) err2 = _axis.error(ib)*_axis.error(ib);
	v   += _axis[ib]/err2;
	ev2 += 1./err2;
      }
    }
    if( ev2 < 0.00001 )  ev2 = 1;
    out->SetPoint(      ip  , bX->itBinCenter(), v/ev2 );
    out->SetPointError( ip++, bX->itBinWidth() , 1./sqrt(ev2));
  }
  return out;
}


template <typename T>
TH2F* IJazZAxisViewer<T>::proj2D( unsigned d1, unsigned d2, const std::vector<T> &xMinus2 ) {
  assert( _axis.getND() - xMinus2.size()  == 2 );
  
  IJazZAxis<T> *bX = _axis.getAxis(d1);
  IJazZAxis<T> *bY = _axis.getAxis(d2);
  float *binX = new float[bX->nBins()+1];
  float *binY = new float[bY->nBins()+1];
  for( bX->itBegin(); ! bX->itEnd();  bX->itForward() ) {
    binX[bX->itBinNumber()+0] = bX->itBinCenter() - bX->itBinWidth();
    binX[bX->itBinNumber()+1] = bX->itBinCenter() + bX->itBinWidth();
  }
  for( bY->itBegin(); ! bY->itEnd();  bY->itForward() ) {
    binY[bY->itBinNumber()+0] = bY->itBinCenter() - bY->itBinWidth();
    binY[bY->itBinNumber()+1] = bY->itBinCenter() + bY->itBinWidth();
  }

  TH2F *h = new TH2F("hProj2D","proj 2D",bX->nBins(),binX,bY->nBins(),binY);
  h->Sumw2();
  for( bX->itBegin(); ! bX->itEnd();  bX->itForward() )
  for( bY->itBegin(); ! bY->itEnd();  bY->itForward() ) {
    std::vector<T> xy; 
    int dMinus2 = 0;
    for( unsigned d = 0 ; d < _axis.getND(); d++ ) {
      if( d == d1 || d == d2 ) xy.push_back( _axis.getAxis(d)->itBinCenter() ); 
      else                     xy.push_back( xMinus2[dMinus2++] ); 
    }
    int ibin = _axis.findBin( xy );
    h->SetBinContent(bX->itBinNumber()+1,bY->itBinNumber()+1,_axis.value(ibin));
    h->SetBinError(  bX->itBinNumber()+1,bY->itBinNumber()+1,_axis.error(ibin));
  }
  h->GetXaxis()->SetTitle(bX->getName().c_str());
  h->GetYaxis()->SetTitle(bY->getName().c_str());
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  return h;
}


template <typename T>
TH2F* IJazZAxisViewer<T>::proj2D( unsigned d1, unsigned d2, T z ) {
  typename std::vector<T> xMinus2; xMinus2.push_back(z);
  return proj2D(d1,d2,xMinus2);
}

template <typename T>
TH2F* IJazZAxisViewer<T>::proj2D( unsigned d1, unsigned d2, T z, T t ) {
  typename std::vector<T> xMinus2; 
  xMinus2.push_back(z);
  xMinus2.push_back(t);
  return proj2D(d1,d2,xMinus2);
}



template <typename T>
TH1F *IJazZAxisViewer<T>::hist1D( unsigned dim, T x ) {
  assert( _axis.getND() == 2 );
  typename std::vector<T> xMinus1; xMinus1.push_back(x);
  return hist1D(dim,xMinus1);
}


template <typename T>
TH1F *IJazZAxisViewer<T>::hist1D( unsigned dim, T x, T y ) {
  assert( _axis.getND() == 3 );
  typename std::vector<T> xMinus1; 
  xMinus1.push_back(x);
  xMinus1.push_back(y);
  return hist1D(dim,xMinus1);
}

template <typename T>
TH1F *IJazZAxisViewer<T>::hist1D( unsigned dim, const std::vector<T> &xMinus1 ) {
  
  IJazZAxis<T> *bX = _axis.getAxis(dim);
  float *binX = new float[bX->nBins()+1];
  for( bX->itBegin(); ! bX->itEnd();  bX->itForward() ) {
    binX[bX->itBinNumber()+0] = bX->itBinCenter() - bX->itBinWidth();
    binX[bX->itBinNumber()+1] = bX->itBinCenter() + bX->itBinWidth();
  }

  TH1F *h = new TH1F(randRootName().c_str(),"proj 1D",bX->nBins(),binX );
  h->Sumw2();
  for( bX->itBegin(); ! bX->itEnd(); bX->itForward() ) {
    double v  = 0;
    double ev2 = 0;
    for( int ib = 0 ; ib < _axis.nBinsND(); ib++ ) {
      bool addBin = false;
      if( fabs( _axis.getBinCenterDimN(ib,dim) - bX->itBinCenter() ) < 0.0001 ) addBin = true;
      if( addBin && xMinus1.size() == _axis.getND() - 1 ) {
	for( unsigned d = 0 ; d < _axis.getND(); d++ ) {
	  if(  d != dim && fabs( _axis.getBinCenterDimN(ib,d) - _axis.getAxis(d)->binCenter(xMinus1[d<dim?d:d-1]) ) >  0.0001 ) addBin = false;	  
	  if( !addBin ) break;
	}
      }
      if( addBin ) {
	double err2 = 999999;
	if( _axis.error(ib) > 0.00001 ) err2 = _axis.error(ib)*_axis.error(ib);
	v   += _axis[ib]/err2;
	ev2 += 1./err2;
      }
    }
    if( ev2 < 0.00001 )  ev2 = 1;
    h->SetBinContent( bX->itBinNumber()+1, v/ev2 );
    h->SetBinError( bX->itBinNumber()+1, 1./sqrt(ev2));
  }
  return h;
}

typedef IJazZAxisViewer<double> ijazzViewer;

#endif
