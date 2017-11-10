#include "interface/IJazZ.hh"
#include "interface/RootUtils.hh"

#include <TSystem.h>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <boost/program_options.hpp>
using namespace std;

void createDataAndMCinput_ECALELF( string ecalelf, string dataIn, string mcIn ); 
void getRunRange_ECALELF(string runrange, int &runMin,int &runMax);


int main( int nargc, char **argv ) {


  //-------------------- Parse IJazZ options -------------------------//
  bool isMC, ijazzClassic(true);
  int  analysis, mcTest, numCPU;
  bool doFit, doAnaFit, outputTree;
  int  debugLevel;
  int  enCorrType;
  string axisDefResp, axisDefReso, runRange, ecalElfCfg, version, fileListData, fileListMC;
  string userDir;
  vector<string> resoDataFile,resoMCFile;
  vector<string> dataERespCorr, mcERespCorr;
  int allEvenOddEvents = 0;
  int noEBfit = 0;
  namespace po = boost::program_options;
  string config_file;
  po::options_description config("Configuration");
  config.add_options()
    ("eleEnCorr"   ,po::value<int>(&enCorrType) ->default_value(0),"sc cluster energy correction (0=EleRegr, 1=EleStd, 2=PhoRegr,...)")
    ("axisDefReso" ,po::value<string>(&axisDefReso) ->default_value("undefined"),"Resolution axis definition")
    ("axisDefResp" ,po::value<string>(&axisDefResp) ->default_value("undefined"),"Response   axis definition")
    ("numCPU"      ,po::value<int>(&numCPU)         ->default_value(-1)         ,"Number of CPU (-1 = all CPU available)")
    ("fileListData",po::value<string>(&fileListData), "List of input files for data")
    ("fileListMC"  ,po::value<string>(&fileListMC  ), "List of input files for mc  ")
    ("mcTest"      ,po::value<int>(&mcTest)         ->default_value(0), "Oversmear truth MC to perform test")
    ("ecalRespCorr",po::value<vector<string> >(&dataERespCorr), "correct  ecal response (data), several corrections are possible (use IJazZ axis format, like a fit result)")
    ("ecalRespCorrMC",po::value<vector<string> >(&mcERespCorr), "correct  ecal response (MC), several corrections are possible (use IJazZ axis format, like a fit result)")
    ("resoData"    ,po::value<vector<string> >(&resoDataFile),"Resolution Fitted in Data (For MC OverSmearing), IJazZ axis format")
    ("resoMC"      ,po::value<vector<string> >(&resoMCFile  ),"Resolution Fitted in MC   (For MC OverSmearing), IJazZ axis format")
    ("doAnaFit"    ,po::value<bool>(  &doAnaFit    )->default_value(false), "produce Eta Scale for analysis +1 [NOT SUPPORTED]")
    ;
   po::options_description generic("Generic options");
   generic.add_options()
     ("help,h", "produce help message")
     ("isMC"         , "Data or MC")
     ("classic"      ,po::value<bool>(  &ijazzClassic     )->default_value(true) ,  "1: classic IJazZ, 0: not yet available ")
     ("noEB"         ,po::value<int>(   &noEBfit          )->default_value(0)    ,  "0: all; 1: noEB; 2: EB+EEm; 3: EB+EEp")
     ("analysis,a"   ,po::value<int>(   &analysis         )->default_value(-1)   , "Analysis type -1: test with 2 params only, 0: user defined axis, 1: eta scale axis")
     ("doFit"        ,po::value<bool>(  &doFit            )->default_value(true) , "perform Z fit")
     ("outputTree"   ,po::value<bool>(  &outputTree       )->default_value(false), "output includes a tree with corrected energy scale ")
     ("debugLevel,d" ,po::value<int>(   &debugLevel       )->default_value(0)    , "Debug level")
     ("version,v"    ,po::value<string>(& version         )->default_value("vTest") , "Versioning the output")
     ("ecalElf"      ,po::value<string>(&ecalElfCfg       )->default_value("undefined"), "EcalElf input trees")
     ("runRange,r"   ,po::value<string>(&runRange         )->default_value("0-999999") , "run range")
     ("userDir,u"    ,po::value<string>(&userDir          )->default_value("user")     , "user output sub-dir")
     ("subSampleOnly",po::value<int>(   &allEvenOddEvents )->default_value(0)     , "all/fit/test evts (0/1/2) default: 0")
     ;

   po::options_description hidden("Hidden options");
   hidden.add_options()
     ("config-file", po::value<string>(&config_file), "config file")
     ;

   po::positional_options_description p;
   p.add("config-file", -1);


   po::options_description cmdline_options;
   cmdline_options.add(generic).add(config).add(hidden);

   po::options_description config_file_options;
   config_file_options.add(config).add(generic);
   
   po::options_description visible("Allowed options");
   visible.add(generic).add(config);
   
   po::variables_map vm;
   po::store(po::command_line_parser(nargc, argv).
	     options(cmdline_options).positional(p).run(), vm);
   po::notify(vm);

   if( vm.count("help") ) {
     cout << visible << endl;
     cout << "Usage: " << endl
	  << "./bin/IJazZexe [options] configFile" << endl;
     return 0;
   }

   ifstream ifs(config_file.c_str());
   bool setupOK = true;
   if (!ifs) {
     cout << "can not open config file: " << config_file << endl;
     setupOK = false;
   } else {
     po::store(parse_config_file(ifs, config_file_options), vm);
     notify(vm);
   }
    
   if( !doFit && !outputTree ) {
     cout << "---- IJazZ: I ve nothing to do, either ask for a fit or just for an output tree." << endl;
     setupOK = false;
   }

   if( !setupOK ) {
     cout << "---- IJazZ bail out: " << endl
	  << "* For more info, try:     " << endl
	  << "./bin/IJazZexe -h     " << endl << endl;
     cout << "* For even more info try: " << endl
	  << "        fabrice.couderc@cea.fr" << endl;
     return 0;
   }

   if( vm.count("isMC")    ) isMC = true;
   else                      isMC = false;
   // if( vm.count("noEB")    ) noEBfit = true;
   // else                      noEBfit = false;
   //-------------------- Read runRange and input file -------------------------//   
   int runmin =     -1;
   int runmax = 999999;
   getRunRange_ECALELF(runRange,runmin,runmax);
   
   if( ecalElfCfg != "undefined" ) {
     int pid = gSystem->GetPid();
     fileListData = "/tmp/input.zfitterFab." + itostr(pid) + ".data" ;
     fileListMC   = "/tmp/input.zfitterFab." + itostr(pid) + ".mc";
     createDataAndMCinput_ECALELF(ecalElfCfg,fileListData,fileListMC);
   }
   
   //------------------------- Setup IJazZ options -----------------------------//
   IJazZ fitter;
   if(  isMC ) fitter.isMC(true);
   
   fitter.defineRunRange(runmin,runmax);
   fitter.setVersion(version);
   fitter.setPrintLevel(debugLevel);
   fitter.setUserDir(userDir);
   if( isMC ) fitter.addResponseCorrection( mcERespCorr   );   
   else       fitter.addResponseCorrection( dataERespCorr );   
   fitter.addResoDataFile( resoDataFile );
   fitter.addResoMCFile(   resoMCFile   );
   fitter.setIIazZClassicMode( ijazzClassic );
   fitter.doFit(doFit);
   fitter.createOutputTree(outputTree);
   fitter.scEnergyCorrection(enCorrType);
   fitter.setAllEvenOddEvents(allEvenOddEvents);

   if     ( analysis == -1 ) fitter.useTestBining();
   else if( analysis ==  0 ) fitter.useExternalBining();
   else if( analysis ==  1 ) fitter.useEtaScaleBining();
   bool absEtaAxis = fitter.setBiningND(axisDefResp,axisDefReso);

   fitter.mcOversmearingStudy(mcTest);
   //   string treename = "simpleNtupleEoverP/SimpleNtupleEoverP";
   //   string treename = "SimpleNtupleEoverP";
   string treename = "selected";
   fitter.addInputFilesMC(   fileListMC   , treename );
   fitter.addInputFilesData( fileListData , treename );
   if( isMC ) fitter.setupNvtxReweighting();	

   
   //--------------------------- fit or simple run ------------------------------//
   if( doFit ) {
     vector<string> ecalp;
     if( noEBfit != 1 ) ecalp.push_back("EB");
     if( absEtaAxis   ) ecalp.push_back("EE");
     else {
       if     ( noEBfit <= 1 && noEBfit >= 0 ) {
	 ecalp.push_back("EEm");
	 ecalp.push_back("EEp");  
       } 
       else if( noEBfit == 2 ) ecalp.push_back("EEm");
       else if( noEBfit == 3 ) ecalp.push_back("EEp");       
     }
     for( unsigned ip = 0; ip < ecalp.size(); ++ip ) {
       if( ecalp[ip] != "EB" ) fitter.isEE();
       else                    fitter.isEE(false);
       if( ecalp[ip] != "EB" && false ) {
	 /// add EB correction and remove potential EE corrections
	 dataERespCorr.push_back( fitter.createOutputFileName("EB") + ".fittedResp" );
	 gSystem->Exec( ("rm -f " + fitter.createOutputFileName(ecalp[ip]) + ".fittedResp").c_str() );
	 fitter.addResponseCorrection( dataERespCorr ); 
       }
       fitter.setupZFitter(numCPU);
       fitter.eventSelectionZFitter(ecalp[ip]);
       fitter.minimize(ecalp[ip]);
    }
     if( doAnaFit    ) fitter.etaScaleFromZeeAnaFit();
   } else {
     cout << "---- IJazZ: no Fit required, I just produce an output tree" << endl;
     fitter.eventSelectionZFitter("Ecal");
   }
   
   if( ecalElfCfg != "undefined" ) {
     string rmdata = "rm -f " + fileListData;
     string rmmc   = "rm -f " + fileListMC  ;
     gSystem->Exec( rmdata.c_str() );
     gSystem->Exec( rmmc.c_str()   );
   }
     
   return 0;
}


void createDataAndMCinput_ECALELF( string ecalelf, string dataIn, string mcIn ) {
  
  cout << " Opening ECALF config: " << ecalelf << endl;
  ifstream input(ecalelf.c_str());
  ofstream outputMC(mcIn.c_str());
  ofstream outputData(dataIn.c_str());

  while ( input.good() && !input.eof() ) {
    string line;
    getline(input,line,'\n');
    if( line.find("#") != string::npos ) continue;
    istringstream parse(line);
    string type, treename,file;
    parse >> type >> treename >> file;
    if( treename == "selected" )  {
      if( type == "d" ) outputData << file << endl;
      if( type == "s" ) outputMC   << file << endl;
    }
  }
  input.close();
  outputData.close();
  outputMC.close();
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
