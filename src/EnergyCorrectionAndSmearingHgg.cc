#include "interface/EnergyCorrectionAndSmearingHgg.hh"

using namespace std;


#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>


int EnergyScaleReader::EcalPart( string ecal ) {
  
  /// internal code convention
  if( ecal == "EBlowEtaBad"   ) return 0;
  if( ecal == "EBlowEtaGold"  ) return 1;
  if( ecal == "EBhighEtaBad"  ) return 2;
  if( ecal == "EBhighEtaGold" ) return 3;
  if( ecal == "EElowEtaBad"   ) return 4;
  if( ecal == "EElowEtaGold"  ) return 5;
  if( ecal == "EEhighEtaBad"  ) return 6;
  if( ecal == "EEhighEtaGold" ) return 7;

  /// read other shervin's table format
  if( ecal == "EB-absEta_0_1-bad"      ) return 0;
  if( ecal == "EB-absEta_0_1-gold"     ) return 1;
  if( ecal == "EB-absEta_1_1.4442-bad" ) return 2;
  if( ecal == "EB-absEta_1_1.4442-gold") return 3;
  if( ecal == "EE-absEta_1.566_2-bad"  ) return 4;
  if( ecal == "EE-absEta_1.566_2-gold" ) return 5;
  if( ecal == "EE-absEta_2_2.5-bad"    ) return 6;
  if( ecal == "EE-absEta_2_2.5-gold"   ) return 7;
  
  return -1;
}

string EnergyScaleReader::EcalPartString( float r9, float scEta ) {
  string ecal = "Unknown";
  if     ( fabs(scEta) <= 1.00 ) ecal = "EBlowEta";
  else if( fabs(scEta) <= 1.45 ) ecal = "EBhighEta";
  else if( fabs(scEta) <= 2.00 ) ecal = "EElowEta";
  else if( fabs(scEta) <= 2.50 ) ecal = "EEhighEta";
  
  if( r9 > 0.94 ) ecal += "Gold";
  else            ecal += "Bad";
  
  return ecal;
}

int EnergyScaleReader::EcalPart(float r9, float scEta ) {
  return EcalPart( EcalPartString( r9, scEta ) );
}


bool EnergyScaleReader::setup( string energyscaleName ) {

   if( energyscaleName.find( "EtDep" ) != string::npos ) {
     _isEtDep = true;
     return setup_EtDep(energyscaleName);
   }
   _isEtDep = false;

   ifstream input(energyscaleName.c_str());

   //// default energyscale values


   cout << " Opening energy scale file: " << energyscaleName<< endl;
   while ( input.good() && !input.eof() ) {
     string line;
    getline(input,line,'\n');
    
    /// comment
    if( line.find("#") != string::npos ) continue;

    istringstream isstream1(line);
    istringstream isstream2(line);
    string ecal;
    int runMin,runMax;
    float eScale, err_eScale;
    isstream1 >> ecal >> runMin >> runMax >> eScale >> err_eScale ;

    if( runMin <= 0 ) {
      string dummy;
      /// most likely we are reading the new Format
      isstream2 >> ecal >> dummy >> runMin >> runMax >> eScale >> err_eScale ;
      
    }
    int iEcal = EcalPart(ecal);
    //cout <<" iEcal "<<iEcal <<endl;

    if( iEcal < 0 ) continue;
    _runStart[  iEcal].push_back( runMin );
    _runStop[   iEcal].push_back( runMax );
    _eScale[    iEcal].push_back( eScale );
    _err_eScale[iEcal].push_back( err_eScale );
   }
   _isSetup = true;
   return true;
}

EnergyScaleReader::EnergyScaleReader(void) {_isSetup = false; _nWarnings = 0;}
EnergyScaleReader::~EnergyScaleReader(void) {}




// *************************************  combine ******************************************

float EnergyScaleReader::energyScale( float r9, float scEta, float Et, int run, int systShift ) {
  
  if( _isEtDep ) 
    return energyScale_EtDep( r9, scEta, Et, run, systShift );
  return energyScale_noEtDep( r9, scEta, run, systShift );
}


// *********************************energy dependent correction******************************************

float EnergyScaleReader::energyScale_EtDep( float r9, float scEta, float Et, int run, int systShift ) {
  
  if( !_isSetup ) {
    cout << " EnergyScaleReader is not setup! use function setup with the ad hoc setup file" << endl;
    return 1;
  }

  int iEcal = EcalPart_EtDep( r9, scEta, Et );
  if( iEcal < 0 ) {
    cout << " EnergyScaleReader unknown Ecal part for: r9 = " << r9 << " eta = " << scEta << " et = " << Et << endl;
    return 1;
  }
  
  int irun = runRange( run, iEcal );

  if( irun < 0 ) {
    cout << " EnergyScaleReader: can not find run: " << run << " for ECAL: " << EcalPartString( r9, scEta ) << endl;
    return 1;
  }

  return _eScale[iEcal][irun] + systShift * _err_eScale[iEcal][irun];
}

// *********************************energy independent correction************************************

float EnergyScaleReader::energyScale_noEtDep( float r9, float scEta, int run, int systShift ) {
  
  if( !_isSetup ) {
    cout << " EnergyScaleReader is not setup! use function setup with the ad hoc setup file" << endl;
    return 1;
  }

  int iEcal = EcalPart( r9, scEta );

  if( iEcal < 0 ) {
    cout << " EnergyScaleReader unknown Ecal part for: r9 = " << r9 << " eta = " << scEta << endl;
    return 1;
  }
    

  int irun = runRange( run, iEcal );
  if( irun < 0 ) {
    cout << " EnergyScaleReader: can not find run: " << run << " for ECAL: " << EcalPartString( r9, scEta ) << endl;
    return 1;
  }
  
  //  cout <<"inside energy scale --->   r9, scEta "<<  r9 <<", "<<  scEta  <<"  iEcal  "<< iEcal << "   _eScale[iEcal][irun] + systShift * _err_eScale[iEcal][irun] "<<  _eScale[iEcal][irun] + systShift * _err_eScale[iEcal][irun] <<endl;

  return _eScale[iEcal][irun] + systShift * _err_eScale[iEcal][irun];
}

// *****************************************************************************************

int EnergyScaleReader::runRange( int run, int iEcal ) {
  for( unsigned irun = 0; irun < _runStart[iEcal].size(); irun++ )
    if( run >= _runStart[iEcal][irun] && run <= _runStop[iEcal][irun] ) return int(irun);


  if( run > _runStop[iEcal][_runStart[iEcal].size()-1] ) {
    if( _nWarnings < 20 ) {
      cout << "run " << run << " > last run known: " << _runStop[iEcal][_runStart[iEcal].size()-1] << " for iEcal " << iEcal << endl;
      _nWarnings++;
    }
    return _runStart[iEcal].size()-1;
  }

  return -1;
}
  
//=====================================================================
//================  Et dependent corrections ================
///=====================================================================
int EnergyScaleReader::EcalPart_EtDep( string ecal ) {


  if( ecal == "absEta_0_1-gold-Et_20_35"      ) return 0;
  if( ecal == "absEta_0_1-gold-Et_35_43"      ) return 1;
  if( ecal == "absEta_0_1-gold-Et_43_50"      ) return 2;
  if( ecal == "absEta_0_1-gold-Et_50_55"      ) return 3;
  if( ecal == "absEta_0_1-gold-Et_55_100"    ) return 4;

  if( ecal == "absEta_0_1-bad-Et_20_33"       ) return 5;
  if( ecal == "absEta_0_1-bad-Et_33_39"       ) return 6;
  if( ecal == "absEta_0_1-bad-Et_39_45"       ) return 7;
  if( ecal == "absEta_0_1-bad-Et_45_50"       ) return 8;
  if( ecal == "absEta_0_1-bad-Et_50_58"       ) return 9;
  if( ecal == "absEta_0_1-bad-Et_58_100"     ) return 10;


  if( ecal == "absEta_1_1.4442-gold-Et_20_40"      ) return 11;
  if( ecal == "absEta_1_1.4442-gold-Et_40_50"      ) return 12;
  if( ecal == "absEta_1_1.4442-gold-Et_50_100"    ) return 13;


  if( ecal == "absEta_1_1.4442-bad-Et_20_33"       ) return 14;
  if( ecal == "absEta_1_1.4442-bad-Et_33_39"       ) return 15;
  if( ecal == "absEta_1_1.4442-bad-Et_39_45"       ) return 16;
  if( ecal == "absEta_1_1.4442-bad-Et_45_50"       ) return 17;
  if( ecal == "absEta_1_1.4442-bad-Et_50_58"       ) return 18;
  if( ecal == "absEta_1_1.4442-bad-Et_58_100"     ) return 19;

  if( ecal == "absEta_1.566_2-gold" ) return 20;
  if( ecal == "absEta_1.566_2-bad"  ) return 21;
  if( ecal == "absEta_2_2.5-gold"   ) return 22;
  if( ecal == "absEta_2_2.5-bad"    ) return 23;
  
  return -1;
}


string EnergyScaleReader::EcalPartString_EtDep( float r9, float scEta, float Et ) {
  string ecal = "Unknown";
  
  if     ( fabs(scEta) <= 1.00 ) ecal = "absEta_0_1";
  else if( fabs(scEta) <= 1.45 ) ecal = "absEta_1_1.4442";
  else if( fabs(scEta) <= 2.00 ) ecal = "absEta_1.566_2";
  else if( fabs(scEta) <= 3.00 ) ecal = "absEta_2_2.5";
  
  if( r9 > 0.94 ) ecal += "-gold";
  else            ecal += "-bad";
  
  if(ecal == "absEta_0_1-gold"){
    if(                Et<35)ecal +="-Et_20_35";
    else if (Et>=35 && Et<43)ecal +="-Et_35_43";
    else if (Et>=43 && Et<50)ecal +="-Et_43_50";
    else if (Et>=50 && Et<55)ecal +="-Et_50_55";
    else if (Et>=55)ecal +="-Et_55_100";
  }else if(ecal == "absEta_0_1-bad" || ecal == "absEta_1_1.4442-bad" ){
    if(                Et<33)ecal +="-Et_20_33";
    else if (Et>=33 && Et<39)ecal +="-Et_33_39";
    else if (Et>=39 && Et<45)ecal +="-Et_39_45";
    else if (Et>=45 && Et<50)ecal +="-Et_45_50";
    else if (Et>=50 && Et<58)ecal +="-Et_50_58";
    else if (Et>=58)ecal +="-Et_58_100";
  }else if(ecal == "absEta_1_1.4442-gold" ){
    if(                Et<40)ecal +="-Et_20_40";
    else if (Et>=40 && Et<50)ecal +="-Et_40_50";
    else if (Et>=50)ecal +="-Et_50_100";
  }

  //  cout <<"ecal "<<   ecal <<endl;
  return ecal;
}

int EnergyScaleReader::EcalPart_EtDep(float r9, float scEta, float Et ) {
  return EcalPart_EtDep( EcalPartString_EtDep( r9, scEta, Et ) );
}

bool EnergyScaleReader::setup_EtDep( string energyscaleName ) {
   ifstream input(energyscaleName.c_str());

   //// default energyscale values  
   if(!input.good()){
     std::cerr << "[ERROR] file " << energyscaleName.c_str() << " not readable" << std::endl;
   }


   cout << " Opening energy scale file: " << energyscaleName << endl;
   while ( input.good() && !input.eof() ) {
     string line;
     getline(input,line,'\n');
    
     
     /// comment
     if( line.find("#") != string::npos ) continue;
     
     //  istringstream isstream1(line);
     istringstream isstream2(line);
     string ecal;
     int runMin,runMax;
     float eScale, err_eScale;
     
     string dummy;
     /// most likely we are reading the new Format
     isstream2 >> ecal >> dummy >> runMin >> runMax >> eScale >> err_eScale ;
     
     
     int iEcal = EcalPart_EtDep(ecal);
     
     
     if( iEcal < 0 ) continue;
     _runStart[  iEcal].push_back( runMin );
     _runStop[   iEcal].push_back( runMax );
     _eScale[    iEcal].push_back( eScale );
     _err_eScale[iEcal].push_back( err_eScale );

   }
   _isSetup = true;
   return true;
}





/// ------------------------------------------------------------------------------ ///

void photonOverSmearing::initialize( string setupType ) {
  _oversmearing.resize(10,-1);
  _oversmearing_err.resize(10,-1);
  _oversmearingstoch_rho.resize(6,-1);
  _oversmearingstoch_rhoerr.resize(6,-1);
  _oversmearingstoch_phi.resize(6,-1);
  _oversmearingstoch_phierr.resize(6,-1);
  _oversmearingstoch_EMean.resize(6,-1);
  _oversmearingstoch_EMeanerr.resize(6,-1);
  setupType_= setupType;
  /// given in percents

  /// Jan16 smearing numbers
  _oversmearing_err[photonOverSmearing::EBlowEtaGold_CM ] = 0.36;
  _oversmearing_err[photonOverSmearing::EBlowEtaGold_GAP] = 0.26;
  _oversmearing_err[photonOverSmearing::EBlowEtaBad_CM  ] = 0.25;
  _oversmearing_err[photonOverSmearing::EBlowEtaBad_GAP ] = 0.25;
  _oversmearing_err[photonOverSmearing::EBhighEtaGold   ] = 0.73;
  _oversmearing_err[photonOverSmearing::EBhighEtaBad    ] = 0.61;
  _oversmearing_err[photonOverSmearing::EElowEtaGold    ] = 0.95;
  _oversmearing_err[photonOverSmearing::EElowEtaBad     ] = 0.34;
  _oversmearing_err[photonOverSmearing::EEhighEtaGold   ] = 0.37;
  _oversmearing_err[photonOverSmearing::EEhighEtaBad    ] = 0.55;

  if( setupType == "Prompt2012_ichep" ) {
    //FOR MORIOND RUN ABCD SMEARING 
    _oversmearing[photonOverSmearing::EBlowEtaGold_GAP] = 1.11;
    _oversmearing[photonOverSmearing::EBlowEtaGold_CM ] = 1.11;
    _oversmearing[photonOverSmearing::EBlowEtaBad_GAP ] = 1.07;
    _oversmearing[photonOverSmearing::EBlowEtaBad_CM  ] = 1.07;
    _oversmearing[photonOverSmearing::EBhighEtaGold   ] = 1.55;
    _oversmearing[photonOverSmearing::EBhighEtaBad    ] = 1.94;
    
    _oversmearing[photonOverSmearing::EElowEtaGold    ] = 2.95;
    _oversmearing[photonOverSmearing::EElowEtaBad     ] = 2.76;
    _oversmearing[photonOverSmearing::EEhighEtaGold   ] = 3.70;    
    _oversmearing[photonOverSmearing::EEhighEtaBad    ] = 3.71;
    _oversmearing_err[photonOverSmearing::EBlowEtaGold_GAP] = sqrt(0.07*0.07+0.22*0.22);
    _oversmearing_err[photonOverSmearing::EBlowEtaGold_CM ] = sqrt(0.07*0.07+0.22*0.22);
    _oversmearing_err[photonOverSmearing::EBlowEtaBad_GAP ] = sqrt(0.06*0.06+0.24*0.24);
    _oversmearing_err[photonOverSmearing::EBlowEtaBad_CM  ] = sqrt(0.06*0.06+0.24*0.24);
    _oversmearing_err[photonOverSmearing::EBhighEtaGold   ] = sqrt(0.40*0.40+0.60*0.60);
    _oversmearing_err[photonOverSmearing::EBhighEtaBad    ] = sqrt(0.11*0.11+0.59*0.59);
    _oversmearing_err[photonOverSmearing::EElowEtaGold    ] = sqrt(0.25*0.25+0.90*0.90);
    _oversmearing_err[photonOverSmearing::EElowEtaBad     ] = sqrt(0.13*0.13+0.30*0.30);
    _oversmearing_err[photonOverSmearing::EEhighEtaGold   ] = sqrt(0.11*0.11+0.34*0.34);
    _oversmearing_err[photonOverSmearing::EEhighEtaBad    ] = sqrt(0.16*0.16+0.52*0.52);

  }

  if( setupType == "oversmear_hcp2012" ) {
    /// HCP2012 numbers
    _oversmearing[photonOverSmearing::EBlowEtaGold_GAP] = 1.10;
    _oversmearing[photonOverSmearing::EBlowEtaGold_CM ] = 1.00;
    _oversmearing[photonOverSmearing::EBlowEtaBad_GAP ] = 1.06;
    _oversmearing[photonOverSmearing::EBlowEtaBad_CM  ] = 1.06;
    _oversmearing[photonOverSmearing::EBhighEtaGold   ] = 1.86;
    _oversmearing[photonOverSmearing::EBhighEtaBad    ] = 1.96;

    _oversmearing[photonOverSmearing::EElowEtaGold    ] = 2.83;
    _oversmearing[photonOverSmearing::EElowEtaBad     ] = 2.67;
    _oversmearing[photonOverSmearing::EEhighEtaGold   ] = 3.43;    
    _oversmearing[photonOverSmearing::EEhighEtaBad    ] = 3.45;
    _oversmearing_err[photonOverSmearing::EBlowEtaGold_GAP] = sqrt(0.14*0.14+0.22*0.22);
    _oversmearing_err[photonOverSmearing::EBlowEtaGold_CM ] = sqrt(0.28*0.28+0.22*0.22);
    _oversmearing_err[photonOverSmearing::EBlowEtaBad_GAP ] = sqrt(0.08*0.08+0.24*0.24);
    _oversmearing_err[photonOverSmearing::EBlowEtaBad_CM  ] = sqrt(0.08*0.08+0.24*0.24);
    _oversmearing_err[photonOverSmearing::EBhighEtaGold   ] = sqrt(0.42*0.42+0.60*0.60);
    _oversmearing_err[photonOverSmearing::EBhighEtaBad    ] = sqrt(0.14*0.14+0.59*0.59);
    _oversmearing_err[photonOverSmearing::EElowEtaGold    ] = sqrt(0.31*0.31+0.90*0.90);
    _oversmearing_err[photonOverSmearing::EElowEtaBad     ] = sqrt(0.17*0.17+0.30*0.30);
    _oversmearing_err[photonOverSmearing::EEhighEtaGold   ] = sqrt(0.15*0.15+0.34*0.34);
    _oversmearing_err[photonOverSmearing::EEhighEtaBad    ] = sqrt(0.18*0.18+0.52*0.52);

  }



  if( setupType == "oversmear_moriond2013" ) {
    _oversmearing[photonOverSmearing::EBlowEtaGold_GAP] = 1.11;
    _oversmearing[photonOverSmearing::EBlowEtaGold_CM ] = 1.11;
    _oversmearing[photonOverSmearing::EBlowEtaBad_GAP ] = 1.07;
    _oversmearing[photonOverSmearing::EBlowEtaBad_CM  ] = 1.07;
    _oversmearing[photonOverSmearing::EBhighEtaGold   ] = 1.55;
    _oversmearing[photonOverSmearing::EBhighEtaBad    ] = 1.94;
    
    _oversmearing[photonOverSmearing::EElowEtaGold    ] = 2.95;
    _oversmearing[photonOverSmearing::EElowEtaBad     ] = 2.76;
    _oversmearing[photonOverSmearing::EEhighEtaGold   ] = 3.70;    
    _oversmearing[photonOverSmearing::EEhighEtaBad    ] = 3.71;
    _oversmearing_err[photonOverSmearing::EBlowEtaGold_GAP] = sqrt(0.07*0.07+0.22*0.22);
    _oversmearing_err[photonOverSmearing::EBlowEtaGold_CM ] = sqrt(0.07*0.07+0.22*0.22);
    _oversmearing_err[photonOverSmearing::EBlowEtaBad_GAP ] = sqrt(0.06*0.06+0.24*0.24);
    _oversmearing_err[photonOverSmearing::EBlowEtaBad_CM  ] = sqrt(0.06*0.06+0.24*0.24);
    _oversmearing_err[photonOverSmearing::EBhighEtaGold   ] = sqrt(0.40*0.40+0.60*0.60);
    _oversmearing_err[photonOverSmearing::EBhighEtaBad    ] = sqrt(0.11*0.11+0.59*0.59);
    _oversmearing_err[photonOverSmearing::EElowEtaGold    ] = sqrt(0.25*0.25+0.90*0.90);
    _oversmearing_err[photonOverSmearing::EElowEtaBad     ] = sqrt(0.13*0.13+0.30*0.30);
    _oversmearing_err[photonOverSmearing::EEhighEtaGold   ] = sqrt(0.11*0.11+0.34*0.34);
    _oversmearing_err[photonOverSmearing::EEhighEtaBad    ] = sqrt(0.16*0.16+0.52*0.52);

  }

  if( setupType == "Legacy2013_8TeV" ) {
    
    _oversmearing[photonOverSmearing::EBlowEtaGold_GAP] = 0.75;
    _oversmearing[photonOverSmearing::EBlowEtaGold_CM ] = 0.75;
    _oversmearing[photonOverSmearing::EBlowEtaBad_GAP ] = 0.86;
    _oversmearing[photonOverSmearing::EBlowEtaBad_CM  ] = 0.86;
    _oversmearing[photonOverSmearing::EBhighEtaGold   ] = 1.22;
    _oversmearing[photonOverSmearing::EBhighEtaBad    ] = 1.88;
    _oversmearing[photonOverSmearing::EElowEtaGold    ] = 1.63;
    _oversmearing[photonOverSmearing::EElowEtaBad     ] = 1.98;
    _oversmearing[photonOverSmearing::EEhighEtaGold   ] = 1.86;
    _oversmearing[photonOverSmearing::EEhighEtaBad    ] = 1.92;
    _oversmearing_err[photonOverSmearing::EBlowEtaGold_GAP] = sqrt(0.03*0.03+0.23*0.23);
    _oversmearing_err[photonOverSmearing::EBlowEtaGold_CM ] = sqrt(0.03*0.03+0.23*0.23);
    _oversmearing_err[photonOverSmearing::EBlowEtaBad_GAP ] = sqrt(0.02*0.02+0.25*0.25);
    _oversmearing_err[photonOverSmearing::EBlowEtaBad_CM  ] = sqrt(0.02*0.02+0.25*0.25);
    _oversmearing_err[photonOverSmearing::EBhighEtaGold   ] = sqrt(0.09*0.09+0.72*0.72);
    _oversmearing_err[photonOverSmearing::EBhighEtaBad    ] = sqrt(0.02*0.02+0.60*0.60);
    _oversmearing_err[photonOverSmearing::EElowEtaGold    ] = 0.903;
    _oversmearing_err[photonOverSmearing::EElowEtaBad     ] = 0.301;
    _oversmearing_err[photonOverSmearing::EEhighEtaGold   ] = 0.341;
    _oversmearing_err[photonOverSmearing::EEhighEtaBad    ] = 0.522;

   //above are previously derived numbers for old smearing method 


   //here compute stochastic smearing term  rho*meansqrtEt*cosphi
   _oversmearingstoch_rho[photonOverSmearing::EBlowEtaGold_GAP] =0.74;
   _oversmearingstoch_rhoerr[photonOverSmearing::EBlowEtaGold_GAP]=0.019;
   _oversmearingstoch_rho[photonOverSmearing::EBlowEtaGold_CM]=0.74;
   _oversmearingstoch_rhoerr[photonOverSmearing::EBlowEtaGold_CM]=0.019;
   _oversmearingstoch_rho[photonOverSmearing::EBlowEtaBad_GAP]=0.77;
   _oversmearingstoch_rhoerr[photonOverSmearing::EBlowEtaBad_GAP]=0.0145;
   _oversmearingstoch_rho[photonOverSmearing::EBlowEtaBad_CM]=0.77;
   _oversmearingstoch_rhoerr[photonOverSmearing::EBlowEtaBad_CM]=0.0145;
   
   _oversmearingstoch_rho[photonOverSmearing::EBhighEtaGold   ]=1.12;
   _oversmearingstoch_rhoerr[photonOverSmearing::EBhighEtaGold   ]=0.0704;
   _oversmearingstoch_rho[photonOverSmearing::EBhighEtaBad   ]=1.26;
   _oversmearingstoch_rhoerr[photonOverSmearing::EBhighEtaBad ]=0.0204;
      
   _oversmearingstoch_phi[photonOverSmearing::EBlowEtaGold_GAP]=0.0;
   _oversmearingstoch_phierr[photonOverSmearing::EBlowEtaGold_GAP]=15.6;
   _oversmearingstoch_phi[photonOverSmearing::EBlowEtaGold_CM]=0.0;
   _oversmearingstoch_phierr[photonOverSmearing::EBlowEtaGold_CM]=15.6;
   _oversmearingstoch_phi[photonOverSmearing::EBlowEtaBad_GAP]=0.0;
   _oversmearingstoch_phierr[photonOverSmearing::EBlowEtaBad_GAP]=16.2;
   _oversmearingstoch_phi[photonOverSmearing::EBlowEtaBad_CM]=0.0;
   _oversmearingstoch_phierr[photonOverSmearing::EBlowEtaBad_CM]=16.2;   

   _oversmearingstoch_phi[photonOverSmearing::EBhighEtaGold   ]=0.0;
   _oversmearingstoch_phierr[photonOverSmearing::EBhighEtaGold   ]=22.3;
   _oversmearingstoch_phi[photonOverSmearing::EBhighEtaBad   ]=0.0;
   _oversmearingstoch_phierr[photonOverSmearing::EBhighEtaBad ]=7.17;

//these are not in percent absolute energy
   _oversmearingstoch_EMean[photonOverSmearing::EBlowEtaGold_GAP]=6.60119;
   _oversmearingstoch_EMeanerr[photonOverSmearing::EBlowEtaGold_GAP]=0.31923;
   _oversmearingstoch_EMean[photonOverSmearing::EBlowEtaGold_CM]=6.60119;
   _oversmearingstoch_EMeanerr[photonOverSmearing::EBlowEtaGold_CM]=0.31923;
   _oversmearingstoch_EMean[photonOverSmearing::EBlowEtaBad_GAP]=6.7276;
   _oversmearingstoch_EMeanerr[photonOverSmearing::EBlowEtaBad_GAP]=0.347156;
   _oversmearingstoch_EMean[photonOverSmearing::EBlowEtaBad_CM]=6.7276;
   _oversmearingstoch_EMeanerr[photonOverSmearing::EBlowEtaBad_CM]=0.347156;

   _oversmearingstoch_EMean[photonOverSmearing::EBhighEtaGold   ]=6.52397;
   _oversmearingstoch_EMeanerr[photonOverSmearing::EBhighEtaGold   ]=0.584685;
   _oversmearingstoch_EMean[photonOverSmearing::EBhighEtaBad   ]=6.73261;
   _oversmearingstoch_EMeanerr[photonOverSmearing::EBhighEtaBad ]=0.29007;  
   
  }

  if( setupType == "Legacy2013_7TeV" ) {
    _oversmearing[photonOverSmearing::EBlowEtaGold_GAP] = 0.68;
    _oversmearing[photonOverSmearing::EBlowEtaGold_CM ] = 0.68;
    _oversmearing[photonOverSmearing::EBlowEtaBad_GAP ] = 0.96;
    _oversmearing[photonOverSmearing::EBlowEtaBad_CM  ] = 0.96;
    _oversmearing[photonOverSmearing::EBhighEtaGold   ] = 1.01;
    _oversmearing[photonOverSmearing::EBhighEtaBad    ] = 1.85;

    _oversmearing[photonOverSmearing::EElowEtaGold    ] = 1.58;
    _oversmearing[photonOverSmearing::EElowEtaBad     ] = 1.85;
    _oversmearing[photonOverSmearing::EEhighEtaGold   ] = 1.83;
    _oversmearing[photonOverSmearing::EEhighEtaBad    ] = 2.01;
    //store error (absolute proofed against globe file)
    _oversmearing_err[photonOverSmearing::EBlowEtaGold_GAP] = 0.221;
    _oversmearing_err[photonOverSmearing::EBlowEtaGold_CM ] = 0.221;
    _oversmearing_err[photonOverSmearing::EBlowEtaBad_GAP ] = 0.241;
    _oversmearing_err[photonOverSmearing::EBlowEtaBad_CM  ] = 0.241;
    _oversmearing_err[photonOverSmearing::EBhighEtaGold   ] = 0.616;
    _oversmearing_err[photonOverSmearing::EBhighEtaBad    ] = 0.591;
    _oversmearing_err[photonOverSmearing::EElowEtaGold    ] = 0.903;
    _oversmearing_err[photonOverSmearing::EElowEtaBad     ] = 0.301;
    _oversmearing_err[photonOverSmearing::EEhighEtaGold   ] = 0.341;
    _oversmearing_err[photonOverSmearing::EEhighEtaBad    ] = 0.522;
  }


  

  /// absolute:
  for( unsigned ios = 0 ; ios < _oversmearing.size()    ; ios++ ) _oversmearing[ios]     /= 100.;
  for( unsigned ios = 0 ; ios < _oversmearing_err.size(); ios++ ) _oversmearing_err[ios] /= 100.;
  for( unsigned ios = 0 ; ios <_oversmearingstoch_rho.size(); ios++ ) {
  _oversmearingstoch_rho[ios]/=100.;
  _oversmearingstoch_rhoerr[ios]/=100.;
  _oversmearingstoch_phi[ios]/=100.;
  _oversmearingstoch_phierr[ios]/=100.;
}
}

void  photonOverSmearing::setSmearSeed(int evtNum, int runNum, int lumisec, float SCPhi, float SCEta){
        UInt_t seedBase = (UInt_t) evtNum + (UInt_t) runNum + (UInt_t) lumisec;
        UInt_t seed1    = seedBase + 100000*(UInt_t) (TMath::Abs(100.*SCPhi)) +1000*(UInt_t) (TMath::Abs(100.*SCEta));
       overSmearingGen.SetSeed(seed1);       
}

photonOverSmearing::category photonOverSmearing::photonOverSmearingCategory( float scEta, float r9, bool isEBGap  ) {
  scEta = fabs(scEta);
  if     ( scEta < 1.0 ) {
    if( r9 > 0.94 ) {
      if( isEBGap ) return photonOverSmearing::EBlowEtaGold_GAP;
      else          return photonOverSmearing::EBlowEtaGold_CM;
    } else {
      if( isEBGap ) return photonOverSmearing::EBlowEtaBad_GAP;
      else          return photonOverSmearing::EBlowEtaBad_CM;
    }
  } else if( scEta < 1.45 )  {
    if( r9 > 0.94 ) return photonOverSmearing::EBhighEtaGold;
    else            return photonOverSmearing::EBhighEtaBad;
  } else if( scEta < 2.0 ) {
    if( r9 > 0.94 ) return photonOverSmearing::EElowEtaGold;
    else            return photonOverSmearing::EElowEtaBad;
  } else if( scEta < 2.5 ) {
    if( r9 > 0.94 ) return photonOverSmearing::EEhighEtaGold;
    else            return photonOverSmearing::EEhighEtaBad;
  }

  return  photonOverSmearing::UNKNOWN;
}


float photonOverSmearing::meanOverSmearing( float scEta, float r9, bool isEBGap, int syst) {
  photonOverSmearing::category cat = photonOverSmearingCategory(scEta, r9, isEBGap);
  if( cat == photonOverSmearing::UNKNOWN ) return 0;  
  
  return _oversmearing[cat] + syst * _oversmearing_err[cat];
}

float photonOverSmearing::meanOverSmearing( float scEta, float r9, float et, bool isEBGap, int syst) {
  photonOverSmearing::category cat = photonOverSmearingCategory(scEta, r9, isEBGap);
  if( cat == photonOverSmearing::UNKNOWN ) return 0;
  float phi   = _oversmearingstoch_phi[cat]+ syst*_oversmearingstoch_phierr[cat];
  float rho   = _oversmearingstoch_rho[cat]+ syst*_oversmearingstoch_rhoerr[cat];
  float EMean = _oversmearingstoch_EMean[cat]+ syst*_oversmearingstoch_EMeanerr[cat];

  float constSmear = rho*sin(phi);
  float stochSmear = rho*cos(phi)*EMean;   
  if(cat == photonOverSmearing::EElowEtaGold  || cat == photonOverSmearing::EElowEtaBad ||
     cat == photonOverSmearing::EEhighEtaGold ||cat==photonOverSmearing::EEhighEtaBad) return _oversmearing[cat] + syst * _oversmearing_err[cat];
  return sqrt( constSmear*constSmear + (stochSmear*stochSmear)/et );
}


float photonOverSmearing::randOverSmearing( float scEta, float r9,  float et, bool isEBGap, int syst) {
  float oversmearing = meanOverSmearing( scEta, r9, isEBGap, syst );
  if(setupType_== "Legacy2013_8TeV"){
  oversmearing=meanOverSmearing( scEta, r9, et,isEBGap, syst );
}

  double scale = -1.0;
  while(scale<0.0){
    scale=overSmearingGen.Gaus( 1.0, oversmearing);
  }

  return scale;
}



