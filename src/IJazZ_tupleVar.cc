 //---------------------------------------------------------------------------//
 // Description: This class is a wrapper so any Z tree can be used as input
 //              The public variables should exist and be properly filled
 //                          via the loadEvent() function
 //                          any selection that are ztree input dependent should
 //                          be implemented in the passPrivateSelection();
 //
 //
 //
 // Author: Fabrice Couderc, fabrice.couderc@cea.fr 
 // History: 2012/11/27 - v1
 //
 //---------------------------------------------------------------------------//

 #include "interface/IJazZ_tupleVar.hh"
 #include "interface/EcalUtils.hh"
 #include "interface/IJazZAxisND.hh"
 #include "interface/eleIDMap.h"

 #include <TChain.h>
 #include <TTree.h>
 #include <TFile.h>
 #include <TH3.h>
 #include <TMath.h>

 #include <iostream>
 #include <vector>
 #include <string>
 #include <cmath>

 //-------------------- function used to find the correction bin, any axis variable can be added
 //-------------------- provided it has handled by IJazZdiEle
 void IJazZdiEle::fill_Xvar_forBinFinder( IJazZAxisND<double> *eCorr, std::vector<double> *x_eCorr ) {
   for( int iele = 0 ; iele < 2; iele++ )
   for( unsigned id = 0 ; id < eCorr->getND(); id++ ) {
     std::string aName = eCorr->getAxis(id)->getName();
     //--- just check the axis name and fill the coordinate accordingly
     if     ( aName == "AbsSCEta"   ) x_eCorr[iele].push_back( fabs(scEta[iele]) );
     else if( aName == "AbsEta"     ) x_eCorr[iele].push_back( fabs(eta[iele])   );
     else if( aName == "AbsIEta"    ) x_eCorr[iele].push_back( fabs(iEta[iele])  );
     else if( aName == "SCEta"      ) x_eCorr[iele].push_back( scEta[iele]       );
     else if( aName == "SCPhi"      ) x_eCorr[iele].push_back( scPhi[iele]       );
     else if( aName == "EtaCrys"    ) x_eCorr[iele].push_back( fmod(convert_eta_ieta(scEta[iele]), 1) );
     else if( aName == "PhiMod"     ) x_eCorr[iele].push_back( fmod(iPhi[iele],20)     );
     else if( aName == "Eta"        ) x_eCorr[iele].push_back( eta[iele]         );
     else if( aName == "IX"         ) x_eCorr[iele].push_back( iX[iele]          );
     else if( aName == "IY"         ) x_eCorr[iele].push_back( iY[iele]          );
     else if( aName == "IEta"       ) x_eCorr[iele].push_back( iEta[iele]        );
     else if( aName == "IPhi"       ) x_eCorr[iele].push_back( iPhi[iele]        );
     else if( aName == "IDee"       ) x_eCorr[iele].push_back( iDee[iele]        );
     else if( aName == "IHarness"   ) x_eCorr[iele].push_back( iHarn[iele]       );
     else if( aName == "LogLC"      ) x_eCorr[iele].push_back( log(avgLC[iele])  );
     else if( aName == "R9"         ) x_eCorr[iele].push_back( R9[iele]          );
     else if( aName == "Rho"        ) x_eCorr[iele].push_back( rho               );
     else if( aName == "nPV"        ) x_eCorr[iele].push_back( nPV               );
     else if( aName == "ExpResEle"  ) x_eCorr[iele].push_back( eSC_err[iele]/eSC[iele] );
     else if( aName == "ExpResPho"  ) x_eCorr[iele].push_back( pho_eSC_err[iele]/pho_eSC[iele] );
     else if( aName == "seedEnFrac" ) x_eCorr[iele].push_back( seedEnF[iele] );
     else {
       std::cerr << " IJazZ::eventSelectionND is not supposed to handle axis: " <<  aName << std::endl;
     }
   }
 }


 //-------------------- IJazZ diEle actual code
 IJazZdiEle::IJazZdiEle(TChain *c): 
   _7TeV(false),_isMC(false),
   _scEnergyCorrType(0)
 { _ztree = dynamic_cast<TTree*>(c); }

 IJazZdiEle::IJazZdiEle(TTree  *c): 
   _7TeV(false),_isMC(false),
   _scEnergyCorrType(0)
 { _ztree = c;}

 IJazZdiEle::~IJazZdiEle(){
   //-- do nothing, trees are extern
 } 

 int  IJazZdiEle::getEntries(void) { return _ztree->GetEntries(); }

 void  IJazZdiEle::setBranchAddress(bool addMCtruth) {
   /// input branch address
   _ztree->SetBranchAddress("runNumber", &runId  );
   _ztree->SetBranchAddress("eventNumber", &evtId  );
   _ztree->SetBranchAddress("nPV"      , &nPV    );
   //   _ztree->SetBranchAddress("weight"   , &weight );
   _ztree->SetBranchAddress("recoFlagsEle"    , tracker);
   _ztree->SetBranchAddress("rawEnergySCEle"  , raw_eSC );
   _ztree->SetBranchAddress("esEnergySCEle"   , raw_eES );
   //   _ztree->SetBranchAddress("energySCEle_corr", std_eSC );
   //   _ztree->SetBranchAddress("seedEnergySCEle" , raw_eSeed );
   if     ( _scEnergyCorrType == 0 ) _ztree->SetBranchAddress("rawEnergySCEle"  , eSC );
   else if( _scEnergyCorrType == 1 ) _ztree->SetBranchAddress("energy_ECAL_ele" , eSC );

   _ztree->SetBranchAddress("rho"     , &rho);
   //   _ztree->SetBranchAddress("es1Ele"  , raw_eES1);
   //   _ztree->SetBranchAddress("es2Ele"  , raw_eES2);
   _ztree->SetBranchAddress("R9Ele"   , R9);
   _ztree->SetBranchAddress("etaEle"  , eta);
   _ztree->SetBranchAddress("phiEle"  , phi);
   _ztree->SetBranchAddress("etaSCEle", scEta);
   _ztree->SetBranchAddress("phiSCEle", scPhi);
   //   _ztree->SetBranchAddress("avgLCSCEle", avgLC  );
   _ztree->SetBranchAddress("eleID", ID);

   _ztree->SetBranchAddress("xSeedSC", seedIx);
   _ztree->SetBranchAddress("ySeedSC", seedIy);
   _ztree->SetBranchAddress("chargeEle" , q);
   if( addMCtruth ) {
     //     _ztree->SetBranchAddress("invMass_MC", &mZ_MC);
     //     _ztree->SetBranchAddress("energyMCEle", energyMCEle );
     //     _ztree->SetBranchAddress("phiMCEle"   , phiMCEle    );
     //     _ztree->SetBranchAddress("etaMCEle"   , etaMCEle    );
   }
 }

 int IJazZdiEle::loadEvent(int ievt) {
   int nb = _ztree->GetEntry(ievt);
   evtQuality = -1;
   if( tracker[0] <= 1 || tracker[1] <= 1 ) return nb;
   if( abs(q[0])  != 1 || abs(q[1]) != 1  ) return nb;
   evtQuality =  0;

   for( int iele = 0 ; iele < 2; iele++ ) {    
     if( _isMC && _7TeV ) {
       scEta[iele] = eta[iele];
       //      scPhi[iele] = phi[iele];
     }
     if( _scEnergyCorrType == 2 ) {
       eSC[iele] = pho_eSC[iele];
       eSC_err[iele] = pho_eSC_err[iele];
     }
     //     raw_eES1[iele] = 0;
     //     raw_eES2[iele] = 0;
     
     seedEnF[iele] = raw_eSeed[iele] / raw_eSC[iele];

     //     if( fabs(scEta[iele]) < 1.49 ) eSC[iele] *= 1.000;
     //     else                           eSC[iele] *= 1.050;

     pt[iele]    = eSC[iele] / cosh(eta[iele]);      
     //     esOverEcal[iele] =  raw_eES[iele] / (raw_eSC[iele]-raw_eES[iele]);
     //     esOverEcal[iele] =  raw_eES[iele] / raw_eSC[iele];
     //     esOverEcal[iele] = 0;

     //// to be changed when switching to actual regression
     //     raw_eES[iele] = 0;
     if     ( _scEnergyCorrType == 0 ) raw_eSC[iele] = eSC[iele];

     if( R9[iele] > 1.0 ) R9[iele] = 0.9999999;

     /// Variables I want to use for corrections so I compute them here.
     iEta[iele]  = _harnessDef.ieta(   scEta[iele],seedIx[iele], seedIy[iele]);
     iPhi[iele]  = seedIy[iele];
     iDee[iele]  = _harnessDef.dee(    scEta[iele],seedIx[iele], seedIy[iele]);
     iHarn[iele] = _harnessDef.harness(scEta[iele],seedIx[iele], seedIy[iele]);
     iX[iele]    = _harnessDef.ix( scEta[iele], seedIx[iele]);
     iY[iele]    = _harnessDef.iy( scEta[iele], seedIy[iele]);

     /// potentially fill here the variables needed by IJazZ
   }


   /*
   unsigned int bitLoose  = 1 << 1;
   unsigned int bitMedium = 1 << 2;
   unsigned int bitTight  = 1 << 3;
   unsigned int bitWP90   = 1 << 4;
   unsigned int bitWP80   = 1 << 5;
   unsigned int bitFid    = 1;

   /// apply WP80
   if(   ( ID[0] & bitLoose  ) != 0 && ( ID[1] & bitLoose  ) != 0 )    evtQuality = 0;
   if( ( ( ID[0] & bitLoose  ) != 0 && ( ID[1] & bitMedium ) != 0 ) ||
       ( ( ID[1] & bitLoose  ) != 0 && ( ID[0] & bitMedium ) != 0 ) )  evtQuality = 1;
   if(   ( ID[0] & bitMedium ) != 0 && ( ID[1] & bitMedium ) != 0 )    evtQuality = 2;
   if( ( ( ID[0] & bitMedium ) != 0 && ( ID[1] & bitTight  ) != 0 ) ||
       ( ( ID[1] & bitMedium ) != 0 && ( ID[0] & bitTight  ) != 0 ) )  evtQuality = 3;
   if(   ( ID[0] & bitTight  ) != 0 && ( ID[1] & bitTight  ) != 0 )    evtQuality = 4;
   if(   ( ID[0] & bitWP80   ) != 0 && ( ID[1] & bitWP80   ) != 0 ) {
     if( evtQuality >= 2 ) evtQuality =  5;
     else                  evtQuality = 10; /// high quality evts not medium-medium
   }
   */
   eleIDMap  idMap;
   //   if( (ID[0] & idMap.eleIDmap["tight25nsRun2"] ) != 0  && (ID[1] & idMap.eleIDmap["tight25nsRun2"] ) != 0  ) evtQuality = 5;
   if( (ID[0] & idMap.eleIDmap["eleID80X-loose"] ) != 0  && 
       (ID[1] & idMap.eleIDmap["eleID80X-loose"] ) != 0  ) evtQuality = 5;

   //   evtQuality = 5;
   /// --- mee before IJazZ potential corrections
   _mee_raw = mee();
   return nb;
 }


bool  IJazZdiEle::passPrivateSelection(void) {
  // std::cout << "evt quality = " <<  evtQuality 
  // 	    << " * q0   = " << q[0]  << " * q1   = " << q[1] 
  // 	    << " * trk0 = " << (int) tracker[0]  << " * trk1 = " << (int) tracker[1] 
  // 	    << " * ID0  = " << ID[0] << " * ID1  = " << ID[1] << std::endl; 

  if( evtQuality < 1 ) return false;
  //  if( R9[0] < 0.94 || R9[1] < 0.94 ) return false;
  // if( fabs(scEta[0]) > 0.8 ) return false;
  // if( fabs(scEta[1]) > 0.8 ) return false;
    for( int iele = 0 ; iele < 2; iele ++ )
     if( fabs(scEta[iele]) >= 1.5 && (seedIx[iele] < 0 || seedIx[iele]> 100 || seedIy[iele] < 0 || seedIy[iele]> 100) ) {
       std::cout << " * iele: " << iele << " - scEta = " << scEta[iele] << "  tracker: " << tracker[iele] << " q: " << q[iele] 
		 << " - Ix  = " << seedIx[iele] << " Iy = " << seedIy[iele] <<  std::endl;
     }
  

  ///  check charge
    if( abs(q[0])  != 1 || abs(q[1]) != 1   ) return false;

  /// require tracker electron
  if( q[0]*q[1] > 0 ) return false;  

  /// for MC purpose
  if( (fabs(eta[0]) < 0.0001 && fabs(phi[0]) < 0.0001) ||
      (fabs(eta[1]) < 0.0001 && fabs(phi[1]) < 0.0001) ) return false;
  return true;
}




void  IJazZdiEle::createOutputTree( TFile *fileout ) {
  evtQuality =  - evtQuality * q[0]*q[1] ;
  /// put the tree where you want in the file
  /// here at the root of directory
  _ztreeOut = new TTree("selected","output IJazZ tree");
  _ztreeOut->SetDirectory(fileout->GetDirectory(""));

  /// output branch address
  _ztreeOut->Branch("runNumber", &runId  ,"runNumber/I"   );
  _ztreeOut->Branch("nPV"      , &nPV    ,"nPV/I"         );
  _ztreeOut->Branch("wei"   , &weight ,"weight/F"      );
  _ztreeOut->Branch("evtQuality"   , &evtQuality ,"evtQuality/I"      );
  _ztreeOut->Branch("energyEle", eSC     ,"energyEle[2]/F");
  _ztreeOut->Branch("ptEle"    , pt      ,"ptEle[2]/F"    );
  _ztreeOut->Branch("energySigmaEle", eSC_err, "energySigmaEle[2]/F");
  _ztreeOut->Branch("R9Ele"    , R9      , "R9Ele[2]/F"    );
  _ztreeOut->Branch("etaEle"   , eta     , "etaEle[2]/F"   );
  _ztreeOut->Branch("corrEle"  , corr    , "corrEle[2]/F"  );
  _ztreeOut->Branch("seedEnFracEle"  , seedEnF    , "seedEnFracEle[2]/F"  );
  _ztreeOut->Branch("avgLCSCEle", avgLC  , "avgLCSCEle[2]/F"   );
  _ztreeOut->Branch("phiEle"   , phi     , "phiEle[2]/F"   );
  _ztreeOut->Branch("etaSCEle" , scEta   , "etaSCEle[2]/F" );
  _ztreeOut->Branch("seedIEtaSCEle", iEta, "seedIEtaSCEle[2]/F");
  _ztreeOut->Branch("mee_raw"  , &_mee_raw , "mee_raw/F");
  _ztreeOut->Branch("mee_corr" , &_mee_corr, "mee_corr/F");
  _ztreeOut->Branch("mee_MC"   , &_mee_MC  , "mee_MC/F");
  _ztreeOut->Branch("mZ_MC"    , &mZ_MC  , "mZ_MC/F");
  _ztreeOut->Branch("fZ"       , &_fZ    , "fZ/F");

  _ztreeOut->Branch("eSC", eSC         ,"eSC[2]/F");
  _ztreeOut->Branch("eMC", energyMCEle ,"eMC[2]/F");
  _ztreeOut->Branch("dR" , deltaR_MC   ,"dR[2]/F");

  _ztreeOut->Branch("Ix" , seedIx , "Ix[2]/F");
  _ztreeOut->Branch("Iy" , seedIy , "Iy[2]/F");
  _ztreeOut->Branch("std_eSC", std_eSC, "std_eSC[2]/F");
  _ztreeOut->Branch("raw_eES" , raw_eES, "raw_eES[2]/F");
  _ztreeOut->Branch("raw_eES1", raw_eES1, "raw_eES1[2]/F");
  _ztreeOut->Branch("raw_eES2", raw_eES2, "raw_eES2[2]/F");
  _ztreeOut->Branch("raw_eSC", raw_eSC, "raw_eSC[2]/F");



  double mean = 0;
  _hMCTruthStudy = new TH3F("hMCTruthStudy","MCTruthStudy",
			    300,0.7-mean,1.3-mean, 
			    250,0.0,2.5,
			    25 ,0.5,1.0);

  _ncat = 12;
  std::string hName;
  for( int i = 0 ; i < _ncat; i++ ) {
    hName = "hMeeTruthStudyMul_cat" + itostr(i);
    _hMeeTruthStudyMul[i] = new TH1F(hName.c_str(),hName.c_str(),300,0.7-mean,1.3-mean );

    hName = "hMeeTruthStudyAdd_cat" + itostr(i);
    _hMeeTruthStudyAdd[i] = new TH1F(hName.c_str(),hName.c_str(),300,0.7-mean,1.3-mean );

    hName = "hMeeRecoStudy_cat" + itostr(i);
    _hMeeRecoStudy[i] = new TH1F(hName.c_str(),hName.c_str(),300,0.7-mean,1.3-mean );

    hName = "hMeeTruth_cat" + itostr(i);
    _hMeeTruth[i] = new TH1F(hName.c_str(),hName.c_str(),300,0.7-mean,1.3-mean );
  }

}

void IJazZdiEle::fillOutputTree() {

  float angle = (cosh(eta[0]-eta[1])-cos(phi[0]-phi[1]))/(cosh(eta[0])*cosh(eta[1]));
  _fZ = (eSC[0]+eSC[1])/(2*mee())*angle;
  /// --- mee after IJazZ potential corrections
  _mee_corr = mee();
  _mee_MC   = sqrt(2*energyMCEle[0]/cosh(etaMCEle[0])*energyMCEle[1]/cosh(etaMCEle[1])*(cosh(etaMCEle[0]-etaMCEle[1])-cos(phiMCEle[0]-phiMCEle[1])));
  _ztreeOut->Fill();  

  double dphi = -1;
  double deta = -1;
  dphi = acos(cos(phi[0]-phiMCEle[0]));
  deta = eta[0]-etaMCEle[0];
  deltaR_MC[0] = sqrt( dphi*dphi + deta*deta);
  dphi = acos(cos(phi[1]-phiMCEle[1]));
  deta = eta[1]-etaMCEle[1];
  deltaR_MC[1] = sqrt( dphi*dphi + deta*deta);

  if( deltaR_MC[0] > 0.05 || deltaR_MC[1] > 0.05 ) return;
  // if( ( fabs( scEta[0]) < 1.5 && int(iPhi[0])%20 < 2.5 ) || 
  //     ( fabs( scEta[1]) < 1.5 && int(iPhi[1])%20 < 2.5 ) ||
  //     fabs( fabs(iEta[0]) - 25 ) < 2.5 || fabs( fabs(iEta[1]) - 25 ) < 2.5 ||
  //     fabs( fabs(iEta[0]) - 45 ) < 2.5 || fabs( fabs(iEta[1]) - 45 ) < 2.5 ||
  //     fabs( fabs(iEta[0]) - 65 ) < 2.5 || fabs( fabs(iEta[1]) - 65 ) < 2.5 ) return;
  
  //  double angReco  = cosh(eta[0]-eta[1])-cos(phi[0]-phi[1]);
  //  double meeMatchMC = sqrt(2*energyMCEle[0]/cosh(etaMCEle[0])*energyMCEle[1]/cosh(etaMCEle[1])*(cosh(etaMCEle[0]-etaMCEle[1])-cos(phiMCEle[0]-phiMCEle[1])));
  //  double angTruth = mZ_MC*mZ_MC/(2*energyMCEle[0]/cosh(etaMCEle[0])*energyMCEle[1]/cosh(etaMCEle[1]));

  double scaleToTruthMult =  sqrt( eSC[0]/energyMCEle[0] * eSC[1]/energyMCEle[1] ) ;
  double scaleToTruthAdd  =  0.5*( eSC[0]/energyMCEle[0] + eSC[1]/energyMCEle[1] ) ;

  int etaCatEle[2];
  int r9CatEle[ 2];

  for( int iele = 0 ; iele < 2; iele++ ) {
    etaCatEle[iele] = 0; 
    r9CatEle[iele]  = 0; 
    if     ( fabs(scEta[iele] ) >= 0.00 && fabs(scEta[iele] ) < 0.40 ) etaCatEle[iele] = 0; 
    else if( fabs(scEta[iele] ) >= 0.40 && fabs(scEta[iele] ) < 0.80 ) etaCatEle[iele] = 1; 
    else if( fabs(scEta[iele] ) >= 0.80 && fabs(scEta[iele] ) < 1.20 ) etaCatEle[iele] = 2; 
    else if( fabs(scEta[iele] ) >= 1.20 && fabs(scEta[iele] ) < 1.55 ) etaCatEle[iele] = 3; 
    else if( fabs(scEta[iele] ) >= 1.55 && fabs(scEta[iele] ) < 2.00 ) etaCatEle[iele] = 4; 
    else if( fabs(scEta[iele] ) >= 2.00 && fabs(scEta[iele] ) < 2.50 ) etaCatEle[iele] = 5;
    if( R9[iele] > 0.94 ) r9CatEle[iele]  = 1;
  }

  int sumR9 = r9CatEle[0] + r9CatEle[1];  
  bool cat[12];
  cat[ 0] = etaCatEle[0] == 0 && etaCatEle[1] == 0 && sumR9 == 2;
  cat[ 1] = etaCatEle[0] == 0 && etaCatEle[1] == 0 && sumR9 == 1;
  cat[ 2] = etaCatEle[0] == 1 && etaCatEle[1] == 1 && sumR9 == 2;
  cat[ 3] = etaCatEle[0] == 1 && etaCatEle[1] == 1 && sumR9 == 1;
  cat[ 4] = etaCatEle[0] == 2 && etaCatEle[1] == 2 && sumR9 == 2;
  cat[ 5] = etaCatEle[0] == 2 && etaCatEle[1] == 2 && sumR9 == 1;
  cat[ 6] = etaCatEle[0] == 0 && etaCatEle[1] == 2 && sumR9 == 0;
  cat[ 7] = etaCatEle[0] == 3 && etaCatEle[1] == 3 && sumR9 == 1;
  cat[ 8] = etaCatEle[0] == 4 && etaCatEle[1] == 4 && sumR9 == 2;
  cat[ 9] = etaCatEle[0] == 4 && etaCatEle[1] == 4 && sumR9 == 1;
  cat[10] = etaCatEle[0] == 5 && etaCatEle[1] == 5 && sumR9 == 2;
  cat[11] = etaCatEle[0] == 5 && etaCatEle[1] == 5 && sumR9 == 1;

  for( int ic = 0 ; ic < _ncat; ic++ ) {
    if( cat[ic] ) _hMeeTruthStudyMul[ic]->Fill( scaleToTruthMult, weight );
    if( cat[ic] ) _hMeeTruthStudyAdd[ic]->Fill( scaleToTruthAdd , weight );
  }

  for( int ic = 0 ; ic < _ncat; ic++ )
    if( cat[ic] ) _hMeeRecoStudy[ic]->Fill( mee() / 91.188 , weight );

  for( int ic = 0 ; ic < _ncat; ic++ )
    if( cat[ic] ) _hMeeTruth[ic]->Fill( mee()/scaleToTruthMult / 91.188, weight );

  int pair = 0; // 0 = EBEB; 1 = EBEE; 2 = EEEE
  for( int iele = 0 ; iele < 2; iele++ ) 
    if( fabs( scEta[iele] ) > 1.5 ) pair++;
    
  for( int iele = 0 ; iele < 2; iele++ ) 
    if( fabs(scEta[iele]) > 1.566 || fabs(scEta[iele]) < 1.4446 ) 
    if( pair == 0 || ( pair == 1 && fabs( scEta[iele] ) > 1.5 ) )
	_hMCTruthStudy->Fill( eSC[iele]/energyMCEle[iele] , fabs(scEta[iele]), R9[iele], weight );

}


void IJazZdiEle::writeTree() {
  _ztreeOut->Write();
  _hMCTruthStudy->Write();
  for( int ic = 0 ; ic < _ncat; ic++ ) _hMeeTruthStudyMul[ic]->Write();
  for( int ic = 0 ; ic < _ncat; ic++ ) _hMeeTruthStudyAdd[ic]->Write();
  for( int ic = 0 ; ic < _ncat; ic++ ) _hMeeRecoStudy[ic]->Write();
  for( int ic = 0 ; ic < _ncat; ic++ ) _hMeeTruth[ic]->Write();
}
