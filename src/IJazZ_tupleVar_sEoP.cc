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
   _ztree->SetBranchAddress("runId"    , &runId  );
   //   _ztree->SetBranchAddress("eventId"  , &evtId  );
   _ztree->SetBranchAddress("PV_n"     , &nPV    );
   //   _ztree->SetBranchAddress("weight"   , &weight );   
   _ztree->SetBranchAddress("isZ"      , &isZ    );   

   /// ele1  variables
   _ztree->SetBranchAddress("ele1_scERaw"  , &raw_eSC_ele1 );   
   _ztree->SetBranchAddress("ele1_es"      , &raw_eES_ele1 );   
   _ztree->SetBranchAddress("ele1_e3x3"    , &raw_e3x3_ele1 );   
   _ztree->SetBranchAddress("ele1_scE"     , &std_eSC_ele1 );
   if( _scEnergyCorrType == 1 )   _ztree->SetBranchAddress("ele1_scE_regression"     , &std_eSC_ele1 );
   _ztree->SetBranchAddress("ele1_seedE" , &raw_eSeed_ele1 );
   _ztree->SetBranchAddress("ele1_eta"   , &eta_ele1);
   _ztree->SetBranchAddress("ele1_phi"   , &phi_ele1);
   _ztree->SetBranchAddress("ele1_scEta" , &scEta_ele1);
   _ztree->SetBranchAddress("ele1_scPhi" , &scPhi_ele1);
   _ztree->SetBranchAddress("ele1_scLaserCorr", &avgLC_ele1  );
   _ztree->SetBranchAddress("ele1_seedIx"  , &seedIx_ele1);
   _ztree->SetBranchAddress("ele1_seedIy"  , &seedIy_ele1);
   _ztree->SetBranchAddress("ele1_seedIeta", &seedIEta_ele1);
   _ztree->SetBranchAddress("ele1_seedIphi", &seedIPhi_ele1);
   _ztree->SetBranchAddress("ele1_charge"  , &q_ele1);


   /// ele2  variables
   _ztree->SetBranchAddress("ele2_scERaw"  , &raw_eSC_ele2 );   
   _ztree->SetBranchAddress("ele2_es"      , &raw_eES_ele2 );   
   _ztree->SetBranchAddress("ele2_e3x3"    , &raw_e3x3_ele2 );   
   _ztree->SetBranchAddress("ele2_scE"     , &std_eSC_ele2 );
   if( _scEnergyCorrType == 1 )   _ztree->SetBranchAddress("ele2_scE_regression"     , &std_eSC_ele2 );
   _ztree->SetBranchAddress("ele2_seedE" , &raw_eSeed_ele2 );
   _ztree->SetBranchAddress("ele2_eta"   , &eta_ele2);
   _ztree->SetBranchAddress("ele2_phi"   , &phi_ele2);
   _ztree->SetBranchAddress("ele2_scEta" , &scEta_ele2);
   _ztree->SetBranchAddress("ele2_scPhi" , &scPhi_ele2);
   _ztree->SetBranchAddress("ele2_scLaserCorr", &avgLC_ele2  );
   _ztree->SetBranchAddress("ele2_seedIx"  , &seedIx_ele2);
   _ztree->SetBranchAddress("ele2_seedIy"  , &seedIy_ele2);
   _ztree->SetBranchAddress("ele2_seedIeta", &seedIEta_ele2);
   _ztree->SetBranchAddress("ele2_seedIphi", &seedIPhi_ele2);
   _ztree->SetBranchAddress("ele2_charge"  , &q_ele2);


   
 }

 int IJazZdiEle::loadEvent(int ievt) {
   int nb = _ztree->GetEntry(ievt);

    /// ele1 variables
   raw_eSC[0] = raw_eSC_ele1;
   raw_eES[0] = raw_eES_ele1;
   eSC[0]   = std_eSC_ele1;
   raw_eSeed[0] = raw_eSeed_ele1;
   R9[0]    = raw_e3x3_ele1/raw_eSC_ele1;
   eta[0]   = eta_ele1;
   phi[0]   = phi_ele1;
   scEta[0] = scEta_ele1;
   scPhi[0] = scPhi_ele1;
   avgLC[0] = avgLC_ele1;
   seedIx[0] = seedIx_ele1; if( seedIx_ele1 < 0 ) seedIx[0] = seedIEta_ele1;
   seedIy[0] = seedIy_ele1; if( seedIy_ele1 < 0 ) seedIy[0] = seedIPhi_ele1;
   q[0]      = q_ele1;

   /// ele2 variables
   raw_eSC[1] = raw_eSC_ele2;
   raw_eES[1] = raw_eES_ele2;
   eSC[1]   = std_eSC_ele2;
   raw_eSeed[1] = raw_eSeed_ele2;
   R9[1]    = raw_e3x3_ele2/raw_eSC_ele2;
   eta[1]   = eta_ele2;
   phi[1]   = phi_ele2;
   scEta[1] = scEta_ele2;
   scPhi[1] = scPhi_ele2;
   avgLC[1] = avgLC_ele2;
   seedIx[1] = seedIx_ele2; if( seedIx_ele2 < 0 ) seedIx[1] = seedIEta_ele2;
   seedIy[1] = seedIy_ele2; if( seedIy_ele2 < 0 ) seedIy[1] = seedIPhi_ele2;
   q[1]      = q_ele2;

   /// apply WP80
   evtQuality = -1;
   if( isZ == 1 ) evtQuality = 5;

   if( evtQuality < 0 ) return -1;

   for( int iele = 0 ; iele < 2; iele++ ) {    
     if( _isMC && _7TeV ) {
       scEta[iele] = eta[iele];
       //      scPhi[iele] = phi[iele];
     }
     if( _scEnergyCorrType == 2 ) {
       eSC[iele] = pho_eSC[iele];
       eSC_err[iele] = pho_eSC_err[iele];
     }
     raw_eES1[iele] = 0;
     raw_eES2[iele] = 0;
     
     seedEnF[iele] = raw_eSeed[iele] / raw_eSC[iele];

     pt[iele]    = eSC[iele] / cosh(eta[iele]);      
     esOverEcal[iele] =  raw_eES[iele] / raw_eSC[iele];

     if( R9[iele] > 1.0 ) R9[iele] = 0.9999999;
     //     std::cout << " e: " << eSC[iele] <<  " eta: " << eta[iele] << " phi: " << phi[iele] << " r9: " << R9[iele] << std::endl;
     /// Variables I want to use for corrections so I compute them here.
     iEta[iele]  = _harnessDef.ieta(   scEta[iele],seedIx[iele], seedIy[iele]);
     iPhi[iele]  = seedIy[iele];
     iDee[iele]  = _harnessDef.dee(    scEta[iele],seedIx[iele], seedIy[iele]);
     iHarn[iele] = _harnessDef.harness(scEta[iele],seedIx[iele], seedIy[iele]);
     iX[iele]    = _harnessDef.ix( scEta[iele], seedIx[iele]);
     iY[iele]    = _harnessDef.iy( scEta[iele], seedIy[iele]);
     /// potentially fill here the variables needed by IJazZ
   }

   /// --- mee before IJazZ potential corrections
    _mee_raw = mee();
    //   std::cout << " Mee = " << _mee_raw << std::endl;
   return nb;
 }


bool  IJazZdiEle::passPrivateSelection(void) {

  if( evtQuality < 1 ) return false;
  //  if( R9[0] < 0.94 || R9[1] < 0.94 ) return false;
  // if( fabs(scEta[0]) > 0.8 ) return false;
  // if( fabs(scEta[1]) > 0.8 ) return false;
  
  /// require tracker electron
  //  if( tracker[0] == 1 || tracker[1] == 1 ) return false;

  //  if( q[0]*q[1] > 0 ) return false;  

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
  
  double angReco  = cosh(eta[0]-eta[1])-cos(phi[0]-phi[1]);
  double meeMatchMC = sqrt(2*energyMCEle[0]/cosh(etaMCEle[0])*energyMCEle[1]/cosh(etaMCEle[1])*(cosh(etaMCEle[0]-etaMCEle[1])-cos(phiMCEle[0]-phiMCEle[1])));
  double angTruth = mZ_MC*mZ_MC/(2*energyMCEle[0]/cosh(etaMCEle[0])*energyMCEle[1]/cosh(etaMCEle[1]));

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

  bool pair = 0; // 0 = EBEB; 1 = EBEE; 2 = EEEE
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
