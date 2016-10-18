//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 12 14:17:59 2016 by ROOT version 6.06/01
// from TTree selected/selected
// found on file: root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/ecalelf/ntuples/13TeV/ALCARERECO/Cal_Sep2016_final_v3/SingleElectron-Run2016G-WSkim-Prompt-v1-weightsReco/278817-279588/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016G-WSkim-Prompt-v1-weightsReco-278817-279588.root
//////////////////////////////////////////////////////////

#ifndef EoverP_shervin_h
#define EoverP_shervin_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class EoverP_shervin {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           runNumber;
   ULong64_t       eventNumber;
   Int_t           lumiBlock;
   UInt_t          runTime;
   Int_t           nBX;
   Float_t         mcGenWeight;
   Char_t          HLTfire;
   Int_t           nPU;
   Float_t         rho;
   Int_t           nPV;
   UInt_t          eleID[3];
   Int_t           chargeEle[3];
   Float_t         etaSCEle[3];
   Float_t         phiSCEle[3];
   Float_t         etaEle[3];
   Float_t         phiEle[3];
   Int_t           classificationEle[3];
   Int_t           recoFlagsEle[3];
   Float_t         PtEle[3];
   Float_t         fbremEle[3];
   Float_t         seedXSCEle[3];
   Float_t         seedYSCEle[3];
   Float_t         seedEnergySCEle[3];
   Float_t         seedLCSCEle[3];
   Float_t         avgLCSCEle[3];
   Float_t         alphaSeedSCEle[3];
   UChar_t         gainEle[3];
   Float_t         energyMCEle[3];
   Float_t         energyEle[3];
   Float_t         energySCEle[3];
   Float_t         energySCEle_must[3];
   Float_t         rawEnergySCEle[3];
   Float_t         rawEnergySCEle_must[3];
   Float_t         energySCEle_must_regrCorr_ele[3];
   Float_t         energySigmaSCEle_must_regrCorr_ele[3];
   Float_t         energySCEle_pho_regrCorr[3];
   Float_t         energySigmaSCEle_pho_regrCorr[3];
   Float_t         esEnergySCEle[3];
   Float_t         esEnergyPlane2SCEle[3];
   Float_t         esEnergyPlane1SCEle[3];
   Float_t         rawESEnergyPlane2SCEle[3];
   Float_t         rawESEnergyPlane1SCEle[3];
   Float_t         energySCEle_corr[3];
   Float_t         R9Ele[3];
   Float_t         e5x5SCEle[3];
   Float_t         efull5x5SCEle[3];
   Float_t         pModeGsfEle[3];
   Float_t         pAtVtxGsfEle[3];
   Float_t         pNormalizedChi2Ele[3];
   Float_t         trackMomentumErrorEle[3];
   Float_t         invMass;
   Float_t         invMass_SC;
   Float_t         invMass_SC_must;
   Float_t         invMass_SC_must_regrCorr_ele;
   Float_t         invMass_SC_pho_regrCorr;
   Float_t         invMass_e5x5;
   Float_t         invMass_efull5x5;
   Float_t         invMass_rawSC;
   Float_t         invMass_rawSC_must;
   Float_t         invMass_rawSC_esSC;
   Float_t         invMass_SC_corr;
   Float_t         invMass_MC;
   Float_t         invMass_mumu;
   Float_t         etaMCEle[3];
   Float_t         phiMCEle[3];
   Int_t           nHitsSCEle[3];

   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_runTime;   //!
   TBranch        *b_nBX;   //!
   TBranch        *b_mcGenWeight;   //!
   TBranch        *b_HLTfire;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_eleID;   //!
   TBranch        *b_chargeEle;   //!
   TBranch        *b_etaSCEle;   //!
   TBranch        *b_phiSCEle;   //!
   TBranch        *b_etaEle;   //!
   TBranch        *b_phiEle;   //!
   TBranch        *b_classificationEle;   //!
   TBranch        *b_recoFlagsEle;   //!
   TBranch        *b_PtEle;   //!
   TBranch        *b_fbremEle;   //!
   TBranch        *b_seedXSCEle;   //!
   TBranch        *b_seedYSCEle;   //!
   TBranch        *b_seedEnergySCEle;   //!
   TBranch        *b_seedLCSCEle;   //!
   TBranch        *b_avgLCSCEle;   //!
   TBranch        *b_alphaSeedSCEle;   //!
   TBranch        *b_gainEle;   //!
   TBranch        *b_energyMCEle;   //!
   TBranch        *b_energyEle;   //!
   TBranch        *b_energySCEle;   //!
   TBranch        *b_energySCEle_must;   //!
   TBranch        *b_rawEnergySCEle;   //!
   TBranch        *b_rawEnergySCEle_must;   //!
   TBranch        *b_energySCEle_must_regrCorr_ele;   //!
   TBranch        *b_energySigmaSCEle_must_regrCorr_ele;   //!
   TBranch        *b_energySCEle_pho_regrCorr;   //!
   TBranch        *b_energySigmaSCEle_pho_regrCorr;   //!
   TBranch        *b_esEnergySCEle;   //!
   TBranch        *b_esEnergyPlane2SCEle;   //!
   TBranch        *b_esEnergyPlane1SCEle;   //!
   TBranch        *b_rawESEnergyPlane2SCEle;   //!
   TBranch        *b_rawESEnergyPlane1SCEle;   //!
   TBranch        *b_energySCEle_corr;   //!
   TBranch        *b_R9Ele;   //!
   TBranch        *b_e5x5SCEle;   //!
   TBranch        *b_efull5x5SCEle;   //!
   TBranch        *b_pModeGsfEle;   //!
   TBranch        *b_pAtVtxGsfEle;   //!
   TBranch        *b_pNormalizedChi2Ele;   //!
   TBranch        *b_trackMomentumErrorEle;   //!
   TBranch        *b_invMass;   //!
   TBranch        *b_invMass_SC;   //!
   TBranch        *b_invMass_SC_must;   //!
   TBranch        *b_invMass_SC_must_regrCorr_ele;   //!
   TBranch        *b_invMass_SC_pho_regrCorr;   //!
   TBranch        *b_invMass_e5x5;   //!
   TBranch        *b_invMass_efull5x5;   //!
   TBranch        *b_invMass_rawSC;   //!
   TBranch        *b_invMass_rawSC_must;   //!
   TBranch        *b_invMass_rawSC_esSC;   //!
   TBranch        *b_invMass_SC_corr;   //!
   TBranch        *b_invMass_MC;   //!
   TBranch        *b_invMass_mumu;   //!
   TBranch        *b_etaMCEle;   //!
   TBranch        *b_phiMCEle;   //!
   TBranch        *b_nHitsSCEle;   //!

   EoverP_shervin(TTree *tree=0);
   virtual ~EoverP_shervin();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(const std::string, const std::vector<Float_t> &, const std::string &);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef EoverP_shervin_cxx
EoverP_shervin::EoverP_shervin(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/ecalelf/ntuples/13TeV/ALCARERECO/Cal_Sep2016_final_v3/SingleElectron-Run2016G-WSkim-Prompt-v1-weightsReco/278817-279588/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016G-WSkim-Prompt-v1-weightsReco-278817-279588.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/ecalelf/ntuples/13TeV/ALCARERECO/Cal_Sep2016_final_v3/SingleElectron-Run2016G-WSkim-Prompt-v1-weightsReco/278817-279588/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016G-WSkim-Prompt-v1-weightsReco-278817-279588.root");
      }
      f->GetObject("selected",tree);

   }
   Init(tree);
}

EoverP_shervin::~EoverP_shervin()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EoverP_shervin::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EoverP_shervin::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EoverP_shervin::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("runTime", &runTime, &b_runTime);
   fChain->SetBranchAddress("nBX", &nBX, &b_nBX);
   fChain->SetBranchAddress("mcGenWeight", &mcGenWeight, &b_mcGenWeight);
   fChain->SetBranchAddress("HLTfire", &HLTfire, &b_HLTfire);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("eleID", eleID, &b_eleID);
   fChain->SetBranchAddress("chargeEle", chargeEle, &b_chargeEle);
   fChain->SetBranchAddress("etaSCEle", etaSCEle, &b_etaSCEle);
   fChain->SetBranchAddress("phiSCEle", phiSCEle, &b_phiSCEle);
   fChain->SetBranchAddress("etaEle", etaEle, &b_etaEle);
   fChain->SetBranchAddress("phiEle", phiEle, &b_phiEle);
   fChain->SetBranchAddress("classificationEle", classificationEle, &b_classificationEle);
   fChain->SetBranchAddress("recoFlagsEle", recoFlagsEle, &b_recoFlagsEle);
   fChain->SetBranchAddress("PtEle", PtEle, &b_PtEle);
   fChain->SetBranchAddress("fbremEle", fbremEle, &b_fbremEle);
   fChain->SetBranchAddress("seedXSCEle", seedXSCEle, &b_seedXSCEle);
   fChain->SetBranchAddress("seedYSCEle", seedYSCEle, &b_seedYSCEle);
   fChain->SetBranchAddress("seedEnergySCEle", seedEnergySCEle, &b_seedEnergySCEle);
   fChain->SetBranchAddress("seedLCSCEle", seedLCSCEle, &b_seedLCSCEle);
   fChain->SetBranchAddress("avgLCSCEle", avgLCSCEle, &b_avgLCSCEle);
   fChain->SetBranchAddress("alphaSeedSCEle", alphaSeedSCEle, &b_alphaSeedSCEle);
   fChain->SetBranchAddress("gainEle", gainEle, &b_gainEle);
   fChain->SetBranchAddress("energyMCEle", energyMCEle, &b_energyMCEle);
   fChain->SetBranchAddress("energyEle", energyEle, &b_energyEle);
   fChain->SetBranchAddress("energySCEle", energySCEle, &b_energySCEle);
   fChain->SetBranchAddress("energySCEle_must", energySCEle_must, &b_energySCEle_must);
   fChain->SetBranchAddress("rawEnergySCEle", rawEnergySCEle, &b_rawEnergySCEle);
   fChain->SetBranchAddress("rawEnergySCEle_must", rawEnergySCEle_must, &b_rawEnergySCEle_must);
   fChain->SetBranchAddress("energySCEle_must_regrCorr_ele", energySCEle_must_regrCorr_ele, &b_energySCEle_must_regrCorr_ele);
   fChain->SetBranchAddress("energySigmaSCEle_must_regrCorr_ele", energySigmaSCEle_must_regrCorr_ele, &b_energySigmaSCEle_must_regrCorr_ele);
   fChain->SetBranchAddress("energySCEle_pho_regrCorr", energySCEle_pho_regrCorr, &b_energySCEle_pho_regrCorr);
   fChain->SetBranchAddress("energySigmaSCEle_pho_regrCorr", energySigmaSCEle_pho_regrCorr, &b_energySigmaSCEle_pho_regrCorr);
   fChain->SetBranchAddress("esEnergySCEle", esEnergySCEle, &b_esEnergySCEle);
   fChain->SetBranchAddress("esEnergyPlane2SCEle", esEnergyPlane2SCEle, &b_esEnergyPlane2SCEle);
   fChain->SetBranchAddress("esEnergyPlane1SCEle", esEnergyPlane1SCEle, &b_esEnergyPlane1SCEle);
   fChain->SetBranchAddress("rawESEnergyPlane2SCEle", rawESEnergyPlane2SCEle, &b_rawESEnergyPlane2SCEle);
   fChain->SetBranchAddress("rawESEnergyPlane1SCEle", rawESEnergyPlane1SCEle, &b_rawESEnergyPlane1SCEle);
   fChain->SetBranchAddress("energySCEle_corr", energySCEle_corr, &b_energySCEle_corr);
   fChain->SetBranchAddress("R9Ele", R9Ele, &b_R9Ele);
   fChain->SetBranchAddress("e5x5SCEle", e5x5SCEle, &b_e5x5SCEle);
   fChain->SetBranchAddress("efull5x5SCEle", efull5x5SCEle, &b_efull5x5SCEle);
   fChain->SetBranchAddress("pModeGsfEle", pModeGsfEle, &b_pModeGsfEle);
   fChain->SetBranchAddress("pAtVtxGsfEle", pAtVtxGsfEle, &b_pAtVtxGsfEle);
   fChain->SetBranchAddress("pNormalizedChi2Ele", pNormalizedChi2Ele, &b_pNormalizedChi2Ele);
   fChain->SetBranchAddress("trackMomentumErrorEle", trackMomentumErrorEle, &b_trackMomentumErrorEle);
   fChain->SetBranchAddress("invMass", &invMass, &b_invMass);
   fChain->SetBranchAddress("invMass_SC", &invMass_SC, &b_invMass_SC);
   fChain->SetBranchAddress("invMass_SC_must", &invMass_SC_must, &b_invMass_SC_must);
   fChain->SetBranchAddress("invMass_SC_must_regrCorr_ele", &invMass_SC_must_regrCorr_ele, &b_invMass_SC_must_regrCorr_ele);
   fChain->SetBranchAddress("invMass_SC_pho_regrCorr", &invMass_SC_pho_regrCorr, &b_invMass_SC_pho_regrCorr);
   fChain->SetBranchAddress("invMass_e5x5", &invMass_e5x5, &b_invMass_e5x5);
   fChain->SetBranchAddress("invMass_efull5x5", &invMass_efull5x5, &b_invMass_efull5x5);
   fChain->SetBranchAddress("invMass_rawSC", &invMass_rawSC, &b_invMass_rawSC);
   fChain->SetBranchAddress("invMass_rawSC_must", &invMass_rawSC_must, &b_invMass_rawSC_must);
   fChain->SetBranchAddress("invMass_rawSC_esSC", &invMass_rawSC_esSC, &b_invMass_rawSC_esSC);
   fChain->SetBranchAddress("invMass_SC_corr", &invMass_SC_corr, &b_invMass_SC_corr);
   fChain->SetBranchAddress("invMass_MC", &invMass_MC, &b_invMass_MC);
   fChain->SetBranchAddress("invMass_mumu", &invMass_mumu, &b_invMass_mumu);
   fChain->SetBranchAddress("etaMCEle", etaMCEle, &b_etaMCEle);
   fChain->SetBranchAddress("phiMCEle", phiMCEle, &b_phiMCEle);
   fChain->SetBranchAddress("nHitsSCEle", nHitsSCEle, &b_nHitsSCEle);
   Notify();
}

Bool_t EoverP_shervin::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EoverP_shervin::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EoverP_shervin::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef EoverP_shervin_cxx
