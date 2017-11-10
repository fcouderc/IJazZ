//---------------------------------------------------------------------------//
// Description: Simple axis classes optimized to search 
//              for a bin number knowing the x coordinate,
//              not to access directly the bin knowing its bin number
//
// Warning:    C++ , template function have to be defined inline only
//
//
// Author: Fabrice Couderc, fabrice.couderc@cea.fr 
// History: 2012/11/19 - v1
//
//---------------------------------------------------------------------------//

#ifndef fabAxis__hh__
#define fabAxis__hh__

#include <math.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>


//-----------------------------------------------------------------------------//
//-------------------------------- 1D axis ------------------------------------//
//-----------------------------------------------------------------------------//

template <typename T> 
class IJazZAxis {
public:
  IJazZAxis();
  ~IJazZAxis();
  IJazZAxis(const IJazZAxis &);
  IJazZAxis& operator=(const IJazZAxis &);
  
  bool addBin(T binEdge);
  bool setBining(unsigned nbin, T* binEdges);
  bool setBining(std::vector<T> binEdges );
  
  inline void setName(std::string name) { _name = name; }
  inline std::string getName(void) const { return _name; }


  inline unsigned nBins(void) const { return _axis.size()>0?_axis.size()-1:0; }
  int findBin(T x) const;
  T binLowEdge(T x) const;
  T binUpEdge( T x) const;
  double binCenter( T x) const;

  void print(void) const;

  inline void reset(void) { _axis.clear(); }
  
  inline void itBegin() { _bin = _axis.begin(); }
  inline unsigned itBinNumber() { return _bin->second; }
  inline void itForward() { ++_bin; }
  bool   itEnd();
  double itBinCenter() ;
  double itBinWidth()  ;

  double getBinCenter( unsigned ib );
  
  void saveToFile(std::string file, unsigned dim = 0);
  void readAxis( std::string axisdef );
private:
  std::string _name;
  std::map<T,unsigned> _axis;
  typename std::map<T,unsigned>::const_iterator _bin;
  void reshuffle();
};



/////////////////////////////////// 1D axis functions ///////////////////////////////////////
template <typename T>
IJazZAxis<T>::IJazZAxis(): _name("notDefined") {}

template <typename T>
IJazZAxis<T>::~IJazZAxis() {}

template <typename T>
IJazZAxis<T>::IJazZAxis(const IJazZAxis &axis) { _axis = axis._axis; _name = axis._name;  _bin  = axis._bin; }

template <typename T>
IJazZAxis<T>& IJazZAxis<T>::operator=(const IJazZAxis<T> &axis) {
 _axis = axis._axis;
 _name = axis._name;
 _bin  = axis._bin;
 return *this;
}

template <typename T>
void IJazZAxis<T>::reshuffle() {
  /// order bins
  typename std::map<T,unsigned>::iterator bin = _axis.begin();
  int ibin = 0;
  for( ; bin != _axis.end(); ++bin ) bin->second = ibin++;
}

template <typename T>
bool IJazZAxis<T>::addBin(T binEdge) {
  _axis.insert( std::make_pair(binEdge,_axis.size()-1) );
  reshuffle();
  return true;
}

template <typename T>
bool IJazZAxis<T>::setBining(unsigned nbins,  T* bins) {
  for( unsigned ib = 0 ; ib <= nbins; ib++ )  
    _axis.insert( std::make_pair(bins[ib],ib) );
  reshuffle();
  return true;
}

template <typename T>
bool IJazZAxis<T>::setBining(std::vector<T> binEdges) {
  for( unsigned ib = 0 ; ib < binEdges.size(); ib++ )  
    _axis.insert( std::make_pair(binEdges[ib],ib) );
  reshuffle();
  return true;
}

template <typename T>
int IJazZAxis<T>::findBin(T x) const {
  typename std::map<T,unsigned>::const_iterator b = _axis.upper_bound( x );
  if( b == _axis.end() ) return nBins();
  return b->second - 1;  
}

template <typename T>
T IJazZAxis<T>::binLowEdge(T x) const {
 typename std::map<T,unsigned>::const_iterator b = _axis.upper_bound( x );
 if( b != _axis.end() && b->second == 0 ) return -999999;
  b--; 
  return b->first;
}

template <typename T>
T IJazZAxis<T>::binUpEdge( T x) const {
  typename std::map<T,unsigned>::const_iterator b = _axis.upper_bound( x );
   if( b == _axis.end() ) return +999999;

  return b->first;
}

template <typename T>
double IJazZAxis<T>::binCenter( T x) const {
  return double( binUpEdge(x) + binLowEdge(x) )/2.;
}


template <typename T>
void IJazZAxis<T>::print(void) const {
 typename std::map<T,unsigned>::const_iterator bin;
 std::cout << "Axis[nbins= " << nBins() << "," << _name <<  "]: ";
 for( bin = _axis.begin(); bin != _axis.end(); ++bin ) 
   std::cout << bin->first << " , ";
 std::cout << std::endl;
}

template <typename T> 
void IJazZAxis<T>::saveToFile(std::string file, unsigned dim)  {
  std::ofstream output(file.c_str(),std::ios_base::app);
  output << "axis[ " << _name << " " << dim << " ]: NBIN= " << _axis.size()<< " : ";
  typename std::map<T,unsigned>::const_iterator bin;
  for( bin = _axis.begin(); bin != _axis.end(); ++bin ) 
    output << bin->first << " ";
  output << std::endl;
  output.close();
}

template <typename T> 
void IJazZAxis<T>::readAxis(std::string axisdef)  {
  std::istringstream parse(axisdef);
  std::string dummy;
  unsigned  nbin(-1);
  parse >> dummy >> _name >> dummy >> dummy >> dummy >> nbin >> dummy;
  for( unsigned ib = 0 ; ib < nbin; ib++ ) {
    T x; parse >> x;
    _axis.insert( std::make_pair(x,ib ) );
  }
}



template <typename T>
bool IJazZAxis<T>::itEnd() { 
  typename std::map<T,unsigned>::const_iterator binNext = _bin;
  binNext++;
  return binNext == _axis.end(); 
}

template <typename T>
double IJazZAxis<T>::itBinCenter() {
  typename std::map<T,unsigned>::const_iterator binNext = _bin; binNext++;
  if( binNext == _axis.end() ) return -999999.;

  return ( double(binNext->first) + double(_bin->first) )/2.;
}
template <typename T>
double IJazZAxis<T>::itBinWidth() {
  return fabs( itBinCenter() - double( _bin->first ) );
}

template <typename T>
double IJazZAxis<T>::getBinCenter( unsigned ib ) {
  typename std::map<T,unsigned>::const_iterator bin = _bin;
  itBegin();
  for( unsigned i = 0; i < ib; i++ ) itForward();
  double c = itBinCenter();
  _bin = bin;
  return c;
}

#endif
