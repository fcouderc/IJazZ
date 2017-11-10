//---------------------------------------------------------------------------//
// Description: Simple axis classes optimized to search 
//              for a bin number knowing the x coordinate,
//              not to access directly the bin knowing its bin number
//              This is a multi dimension axis, potential content per bin supported
//
// Warning:    C++ , template function have to be defined inline only
//
// Author: Fabrice Couderc, fabrice.couderc@cea.fr 
// History: 2012/11/19 - v1
//
//---------------------------------------------------------------------------//

#ifndef fabAxisND__hh__
#define fabAxisND__hh__

#include "interface/IJazZAxis.hh"

#include <math.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>



//-----------------------------------------------------------------------------//
//------------------------------ Multi D axis ---------------------------------//
//-----------------------------------------------------------------------------//
template <typename T> 
class IJazZAxisND {
public:
  IJazZAxisND();
  ~IJazZAxisND();
  IJazZAxisND(const IJazZAxisND &) ;
  IJazZAxisND& operator=(const IJazZAxisND &);
  
  bool addAxis( const IJazZAxis<T> & f);
  void print(int iprint = 0);

  /// bin number in 1D
  int findBin( unsigned dim, T x ) const;
  
  /// bin number in ND
  int findBin( const std::vector<T> & x ) const;

  /// total number of bin
  int nBinsND(void) const;

  void reset(void);
  
  IJazZAxis<T>* getAxis( unsigned d ) { return d >= _axis_nd.size() ? (IJazZAxis<T>*) 0: &_axis_nd[d]; }
  unsigned getND(void) const { return _axis_nd.size(); }
  
  double& operator[]( int ib) { return ib < int(_content.size()) ? _content[ib]: _tmp; }
  double& value(int ib) { return ib < int(_content.size()) ? _content[ib]: _tmp; }
  double& error(int ib) { return ib < int(_error.size())   ? _error[ib]: _tmp; }
  double& value(const std::vector<T> &x) { return value(findBin(x)); }
  double& error(const std::vector<T> &x) { return error(findBin(x)); }

  void setValue(const std::vector<T> &x, double val) { int ib = findBin(x);  ib < int(_content.size()) ? _content[ib] = val : _tmp;  }
  void setError(const std::vector<T> &x, double val) { int ib = findBin(x);  ib < int(_error  .size()) ? _error[ib]   = val : _tmp;  }
  

  double getValue(int ib) const { return ib < int(_content.size()) ? _content[ib]: _tmp; }
  double getError(int ib) const { return ib < int(_error.size())   ? _error[ib]: _tmp; }
  double getValue(const std::vector<T> &x) const { return getValue(findBin(x)); }
  double getError(const std::vector<T> &x) const { return getError(findBin(x)); }


  void saveToFile( std::string file );
  bool readFromFile( std::string file );
  double getBinCenterDimN(  int ibin, unsigned dim );
  double getBinCenterDimN(  int ibin, int dim ) { return dim < 0 ? -99999 : getBinCenterDimN(ibin,unsigned(dim)); }
  unsigned getBinDimN( unsigned ibin, unsigned dim );

private:
  std::vector<IJazZAxis<T> >     _axis_nd;
  typename std::vector<double> _content;
  typename std::vector<double> _error;
  double _tmp;
};




/////////////////////////////////// MultiD axis functions ///////////////////////////////////////
template <typename T> 
IJazZAxisND<T>::IJazZAxisND(void) {}

template <typename T> 
IJazZAxisND<T>::~IJazZAxisND(void) {}

template <typename T> 
IJazZAxisND<T>::IJazZAxisND(const IJazZAxisND &ndAxis) { 
  _axis_nd = ndAxis._axis_nd; 
  _content = ndAxis._content; 
  _error   = ndAxis._error;
}

template <typename T> 
IJazZAxisND<T>& IJazZAxisND<T>::operator=(const IJazZAxisND &ndAxis) {
  _axis_nd = ndAxis._axis_nd;
  _content = ndAxis._content; 
  _error   = ndAxis._error;
  return *this;
}

template <typename T> 
bool IJazZAxisND<T>::addAxis( const IJazZAxis<T> &axis ) { 
    _axis_nd.push_back(axis);
    _content.resize( nBinsND() );
    _error  .resize( nBinsND() );
    return true; 
}

template <typename T> 
void IJazZAxisND<T>::print(int iprint)  {
  std::cout << " IJazZAxisND has " << _axis_nd.size() << " Dimensions. Print each axis" << std::endl;
  for( unsigned id = 0 ; id < _axis_nd.size(); id++ ) {
    std::cout << " D" << id+1 << ": "; _axis_nd[id].print();
  }
  if( iprint >= 1 )
  for( unsigned ib = 0 ; ib < _content.size(); ib++ ) {
    std::cout << "bin[" << ib << "] = " << value(ib) << " +/- " << error(ib) << std::endl;
  } 
}

template <typename T> 
int IJazZAxisND<T>::findBin( unsigned d, T x ) const {
  if( d >= _axis_nd.size() ) {
    std::cerr << "IJazZAxisND::findBin: dimension " << d << " does not exist, D max = " <<  _axis_nd.size()-1 << std::endl;
    return -1;
  }
  return _axis_nd[d].findBin(x);
}

template <typename T> 
int IJazZAxisND<T>::findBin( const std::vector<T> &x ) const {
  if( x.size() != _axis_nd.size() ) {
    std::cerr << "IJazZAxisND::findBin: vector dimension " << x.size() << " should be equal to number of axis " <<  _axis_nd.size() << std::endl;
    return -1;
  }
  int out = 0;
  unsigned ProdD = 1;
  for( unsigned d = 0 ; d < _axis_nd.size(); ++d ) { 
    int locBin =  _axis_nd[d].findBin( x[d] );
    if( locBin <= 0 ) locBin = 0;
    if( locBin ==  int(_axis_nd[d].nBins()) ) locBin = _axis_nd[d].nBins()-1;
    
    out   += locBin*ProdD;
    ProdD *= _axis_nd[d].nBins();
  }

  return out;
}

template <typename T> 
unsigned IJazZAxisND<T>::getBinDimN( unsigned ibin, unsigned dim ) {
  unsigned ProdD = nBinsND();
  unsigned loc = ibin;
  for( unsigned d = _axis_nd.size()-1 ; d > dim; --d ) {
    ProdD /= _axis_nd[d].nBins();
    loc=loc%ProdD;
  }
  ProdD /= _axis_nd[dim].nBins();
  return loc/ProdD;
}


template <typename T> 
double IJazZAxisND<T>::getBinCenterDimN( int ibin, unsigned dim ) {
  return _axis_nd[dim].getBinCenter( getBinDimN(ibin,dim) );
}



/// total number of bin
template <typename T> 
int IJazZAxisND<T>::nBinsND(void) const {
    unsigned ProdD = 1;
  for( unsigned d = 0 ; d < _axis_nd.size(); ++d ) { 
    ProdD *= _axis_nd[d].nBins();
  }
  return ProdD;
}

template <typename T> 
void IJazZAxisND<T>::reset(void)  {
  for( unsigned d = 0 ; d < _axis_nd.size(); ++d )
    _axis_nd[d].reset();
  _axis_nd.clear();
}

template <typename T> 
void IJazZAxisND<T>::saveToFile(std::string file)  {
  std::ofstream output1(file.c_str(),std::ios_base::out);
  output1 << "## FabHisto ND class plain ascii output" << std::endl;
  output1.close();
  for( unsigned d = 0 ; d < _axis_nd.size(); ++d )
    _axis_nd[d].saveToFile(file,d);

  std::ofstream output(file.c_str(),std::ios_base::app);
  output.precision( 7);
  for( unsigned ib = 0 ; ib < _content.size(); ib++ )
    output << _content[ib] << " " << _error[ib] << std::endl;
  output.close();
}

template <typename T> 
bool IJazZAxisND<T>::readFromFile(std::string file)  {
  std::cout << " IJazZAxisND::reading: " << file << std::endl;
  std::ifstream input(file.c_str());
  if( !input.good() ) {
    std::cout << "    * IJazZ can not open file: " << file << std::endl;
    return false;
  }
  while ( input.good() && !input.eof() ) {
    std::string line;
    getline(input,line,'\n');
    if( line.find("#") != std::string::npos ) continue;
    if( line.size() > 5 && line.substr(0,4) == "axis" ) {
      IJazZAxis<T> axis; axis.readAxis(line);
      addAxis( axis );
      _content.clear();
      _error.clear();
      continue;
    }
    if( line.size() <= 2 ) continue;
    std::istringstream parse(line);
    double c, e;
    parse >> c >> e;
    _content.push_back(c);
    _error  .push_back(e);
  }
  if( _content.size() != unsigned(nBinsND()) ) _content.resize( nBinsND() );
  if( _error  .size() != unsigned(nBinsND()) ) _error  .resize( nBinsND() );

  input.close();
  return true;
}


typedef IJazZAxisND<double> ijazzAxisND;
#endif
