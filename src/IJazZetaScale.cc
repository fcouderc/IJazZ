#include "interface/IJazZ.hh"
#include "interface/RootUtils.hh"
#include "interface/EcalUtils.hh"

#include <TSystem.h>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <boost/program_options.hpp>
using namespace std;

void getRunRange_ECALELF(string runrange, int &runMin,int &runMax);

int main( int nargc, char **argv ) {

  //-------------------- Parse IJazZ options -------------------------//
  string version;
  string versMC;
  string userDir;
  string runRange;
  bool save(false), ijazzClassic(false);
  namespace po = boost::program_options;
  po::options_description generic("Generic options");
  generic.add_options()
     ("help,h", "produce help message")
     ("classic"     , "use IJazZ classic mode (removes low mass tail)")
     ("version,v"   ,po::value<string>(& version )->default_value("undefined") , "Versioning the output")
     ("versMC"      ,po::value<string>(& versMC  )->default_value("undefined") , "Version for MC if different from Data")
     ("userDir,u"   ,po::value<string>(&userDir  )->default_value("undefined") , "user output sub-dir")
     ("runRange,r"  ,po::value<string>(&runRange )->default_value("0-999999" ) , "run range")
     ("saveIOV,s"   ,po::value<bool>(  &save )->default_value(false) , "save in IOV format")
     ;

   po::options_description cmdline_options;
   cmdline_options.add(generic);

   po::options_description visible("Allowed options");
   visible.add(generic);
   
   po::variables_map vm;
   po::store(po::command_line_parser(nargc, argv).
	     options(cmdline_options).run(), vm);
   po::notify(vm);

   if( vm.count("help") ) {
     cout << visible << endl;
     cout << "Usage: " << endl
	  << "./bin/IJazZetaScale [options] configFile" << endl;
     return 0;
   }

   if( vm.count("classic") ) ijazzClassic = true;
   else                      ijazzClassic = false;
   
   if( version == "undefined" ) {
     cout << "version not defined... " << endl;
     cout << "./bin/IJazZetaScale --help" << endl;
     return 1;
   }

   if( userDir == "undefined" ) {
     cout << "userDir not defined... " << endl;
     cout << "./bin/IJazZetaScale --help" << endl;
     return 1;
   }

  if( versMC == "undefined" ) versMC = version;
  
    
   
   //-------------------- Read runRange and input file -------------------------//   
   int runmin =     -1;
   int runmax = 999999;
   getRunRange_ECALELF(runRange,runmin,runmax);
   
   
   //-------------------- Setup IJazZ options to get filenames ------------------//
   IJazZ fitter;
   fitter.defineRunRange(runmin,runmax);
   fitter.setUserDir(userDir);
   fitter.useExternalBining();
   fitter.setIIazZClassicMode( ijazzClassic );
   fitter.doFit(1);

   fitter.setVersion(version);
   fitter.isMC(false);
   string dataRespFile = fitter.createOutputFileName("EB") + ".fittedResp";
   fitter.setVersion(versMC );
   fitter.isMC(true);
   string mcRespFile   = fitter.createOutputFileName("EB") + ".fittedResp";
   cout << "etaScale from data: " << dataRespFile << endl;
   cout << "etaScale from MC  : " <<   mcRespFile << endl;
   etaScaleFromZee_AnaFit( dataRespFile, mcRespFile, save);

   return 0;
}

#include <stdlib.h>
void getRunRange_ECALELF(string runrange, int &runMin,int &runMax) {
  string sep = "-";
  int posSep = runrange.find(sep);
  string runMinStr = runrange.substr(0,posSep);
  string runMaxStr = runrange.substr(posSep+1,runrange.size());
  runMin = atoi( runMinStr.c_str() );
  runMax = atoi( runMaxStr.c_str() );
}
