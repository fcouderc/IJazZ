void CloneEoPFormat(void) {
  TString inputDir  = "~/eos/cms/store/group/alca_ecalcalib/ecalMIBI/amartell/MaterialBudget/";
  TString outputDir = "EOP_data/";
 
  //  TString subDir  = "DYEE_8TeV_53X_MatBudget_Summer12ExtendedGeo14DR-PU2012ABCD_START53_V29A-v3_GEN-SIM-RECODEBUG";
  //  TString subDir  = "MC_GEN-SIM-RECODEBUG";
  //  TString subDir  = "DoubleElectron_Run2012ABCD-22Jan2013-v1_AOD" ;
  TString subDir = "DYToEE_M20_powheg-Summer12-START53-ZSkim-runDependent";  
TString outFile = outputDir + subDir + ".skimmed.root"; 
  TString treeName  = "simpleNtupleEoverP/SimpleNtupleEoverP";
  //  TFile *fOld = new TFile(treeInput,"read");
  TChain *tOld  = new TChain(treeName);
  cout << " adding: " <<  inputDir + subDir  + "/*root" << endl;
  tOld->Add( inputDir + subDir  + "/*root");
  //  tOld->Add( inputDir + subDir + "/DYToEE_M_20_TuneZ2star_8TeV_pythia6_Summer12_DR53X-NominalGeo_PU_S10_START53_V15-v1_GEN-SIM-RECODEBUG.root");

  tOld->SetBranchStatus("*",0);
  tOld->SetBranchStatus("runId",1);
  tOld->SetBranchStatus("eventId",1);
  tOld->SetBranchStatus("PV_n",1);
  tOld->SetBranchStatus("isZ",1);
  tOld->SetBranchStatus("ele1_scERaw",1);
  tOld->SetBranchStatus("ele1_es",1);
  tOld->SetBranchStatus("ele1_e3x3",1);
  tOld->SetBranchStatus("ele1_scE",1);
  tOld->SetBranchStatus("ele1_scE_regression",1);
  tOld->SetBranchStatus("ele1_seedE",1);
  tOld->SetBranchStatus("ele1_eta",1);
  tOld->SetBranchStatus("ele1_phi",1);
  tOld->SetBranchStatus("ele1_scEta",1);
  tOld->SetBranchStatus("ele1_scPhi",1);
  tOld->SetBranchStatus("ele1_scLaserCorr",1);
  tOld->SetBranchStatus("ele1_seedIx",1);
  tOld->SetBranchStatus("ele1_seedIy",1);
  tOld->SetBranchStatus("ele1_seedIeta",1);
  tOld->SetBranchStatus("ele1_seedIphi",1);
  tOld->SetBranchStatus("ele1_charge",1);
  tOld->SetBranchStatus("ele2_scERaw",1);
  tOld->SetBranchStatus("ele2_es",1);
  tOld->SetBranchStatus("ele2_e3x3",1);
  tOld->SetBranchStatus("ele2_scE",1);
  tOld->SetBranchStatus("ele2_scE_regression",1);
  tOld->SetBranchStatus("ele2_seedE",1);
  tOld->SetBranchStatus("ele2_eta",1);
  tOld->SetBranchStatus("ele2_phi",1);
  tOld->SetBranchStatus("ele2_scEta",1);
  tOld->SetBranchStatus("ele2_scPhi",1);
  tOld->SetBranchStatus("ele2_scLaserCorr",1);
  tOld->SetBranchStatus("ele2_seedIx",1);
  tOld->SetBranchStatus("ele2_seedIy",1);
  tOld->SetBranchStatus("ele2_seedIeta",1);
  tOld->SetBranchStatus("ele2_seedIphi",1);
  tOld->SetBranchStatus("ele2_charge",1);

  cout << " Cloning in NEW FILE: " << outFile << endl;
  TFile *fNew = new TFile(outFile,"recreate");
  TTree *tNew = tOld->CloneTree(); 
  fNew->Write();
  delete fNew;




}
