#include "interface/IJazZAxis.hh"
#include "interface/IJazZAxisND.hh"
#include <iostream>

using namespace std;
int main(void) {
  double x = 4.8;
  double y = 122;
  IJazZAxis<double> axis;
  double xx[6] = {-3, -6, 7, 9.4, 5, 12 };
  axis.setBining( 5, xx );
  axis.print();
  axis.setName("axisX");
  cout << " axis name: " << axis.getName() << endl;
  cout << " bin(3.5) = " << axis.findBin(3.5) << endl;
  cout << " -low edge(3.5) = " <<  axis.binLowEdge(3.5) << endl;
  cout << " -up  edge(3.5) = " <<  axis.binUpEdge(3.5)  << endl;

  cout << " bin(-12) = " << axis.findBin(-12) << endl;
  cout << " -low edge(-12) = " <<  axis.binLowEdge(-12) << endl;
  cout << " -up  edge(-12) = " <<  axis.binUpEdge(-12)  << endl;


  cout << " bin(123) = " << axis.findBin(123) << endl;
  cout << " -low edge(123) = " <<  axis.binLowEdge(123) << endl;
  cout << " -up  edge(123) = " <<  axis.binUpEdge(123)  << endl;

  axis.addBin( 125 );
  axis.print();

  cout << " bin(123) = " << axis.findBin(123) << endl;
  cout << " -low edge(123) = " <<  axis.binLowEdge(123) << endl;
  cout << " -up  edge(123) = " <<  axis.binUpEdge(123)  << endl;

  cout << endl;
  cout << endl;
  cout << " =========== ND axis checks ========= " << endl;
  IJazZAxis<double> axis_y;
  double yy[8] = {-55, -12.65, 4444, 9.4, 5, 12, 904, 6.543 };
  axis_y.setBining( 7, yy );
  IJazZAxisND<double> AxisND;
  AxisND.addAxis( axis );
  AxisND.addAxis( axis_y );
  AxisND[5] = 145.;
  AxisND[34] = -8.;
  AxisND[128] = -120.;
  //  AxisND.print();
  AxisND.saveToFile("toto2D.out");
  vector<double> xy; xy.resize(2);
  xy[0] = x;
  xy[1] = y;
  cout << " binX = " << axis.findBin(x) << endl;
  cout << " binY = " << axis.findBin(y) << endl;
  cout << " ND bins(" << x << ", " << y <<" ) = " << AxisND.findBin( xy ) << endl; 

  int mybin = AxisND.findBin( xy ); 
  if( mybin >= 0 ) {
    cout << "    --> binX = " << AxisND.getBinDimN( mybin, 0 ) << endl;
    cout << "    --> binY = " << AxisND.getBinDimN( mybin, 1 ) << endl;
  }

  cout << endl;
  cout << endl;
  cout << " ======== check bin iteration ====== " << endl;
  for( axis.itBegin(); ! axis.itEnd(); axis.itForward() ) {
    cout << " Bin[" << axis.itBinNumber() << "]: " << axis.itBinCenter() << " width: " << axis.itBinWidth() << endl;
  }
  for( unsigned ib = 0 ; ib < axis.nBins(); ib++ )
    cout << " Bin[" << ib << "]: " << axis.getBinCenter(ib) << endl;

  cout << " one more time... : " << endl;
  for( axis.itBegin(); ! axis.itEnd(); axis.itForward() ) {
    cout << " BinCenter: " << axis.itBinCenter() << " width: " << axis.itBinWidth() << endl;
  }
  
  IJazZAxis<double> xAxis;
  xAxis.readAxis("axis[ notDefined 1 ]: NBIN= 8 : -55 -12.65 5 6.543 9.4 12 904 4444 ");

  IJazZAxisND<double> histo;
  histo.readFromFile("toto2D.out");

  return 1;
}

