#define EoverP_cxx
#include "EoverP.h"

//#include "makeEcorrWeight.C"

#include <TROOT.h>
#include <TAttFill.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TVirtualFitter.h>

#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h                                                           
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>      // std::istringstream ; to read array of numbers from a line in a file                          
#include <string>
#include <vector>
#include <iomanip> //for input/output manipulators 

#include <Rtypes.h> // to use kColor

#include "histoFunc.h"
#include "Option.h"

#define DATA2016 1
#define SKIM_1LEP1JET_80X 0
#define SKIM_1TLEP30 1 // overrides the previous, but better to have only one with 1 
#define FIT_2SIDE_CB 1      // 0 for single tail Crystal Ball for the fit, 1 for double tail
#define SET_SCALE_ON_Y 0    // select default or user defined ranges for y axis
#define USE_E 1 // 0 for ET and 1 for E in the binning
#define USE_RAWE 0 // when 1, use raw SC energy instead of regression corrected ECAL energy
#define NO_MC_WEIGHT 0   // if you don't want to weight MC with scale factors, PU ecc...
#define ELE_ETA_MAX 1.0  // cut on eta
#define MC_NLO_NOHTBIN 1  // to use NLO MC for W which is not binned in HT, it works with SKIM_1TLEP30 set to 1


using namespace std;


//=====================================================================

// class Option {

// public:
//   Option();
//   ~Option();
//   Int_t Get_data2016() const { return data2016;}
//   Int_t Get_skim1lep1jet80X() const { return skim1lep1jet80X;}
//   Int_t Get_fit2sideCB() const { return fit2sideCB;}
//   Int_t Get_setScaleOnY() const { return setScaleOnY;}
//   Int_t Get_useE() const { return useE;}
//   string Get_dirName() const { return dirName;}

//   void Set_data2016(const Int_t &value) { data2016 = value; }
//   void Set_skim1lep1jet80X(const Int_t &value) { skim1lep1jet80X = value; }
//   void Set_fit2sideCB(const Int_t &value) { fit2sideCB = value; }
//   void Set_setScaleOnY(const Int_t &value) { setScaleOnY = value; }
//   void Set_useE(const Int_t &value) { useE = value; }
//   void Set_dirName(const string &value) { dirName = value; }


// private:
//   Int_t data2016;
//   Int_t skim1lep1jet80X;
//   Int_t fit2sideCB;
//   Int_t setScaleOnY;
//   Int_t useE;
//   string dirName;

// };  

// Option::Option() {
//   data2016 = 1;
//   skim1lep1jet80X = 1;
//   fit2sideCB = 0;
//   setScaleOnY = 0;
//   useE = 0;
//   dirName = "default";
// }

// Option::~Option() {
//   cout<<"Option::~Option() called"<<endl;
// }

//=====================================================================

void getEoP_templateMC(TH1F *hinput = NULL, const string& dirName = "") {

  string fileName = dirName + "EoverP_WJetsToLNu.root";  // template is created when running on MC, and is stored in this file

  TH1::AddDirectory(kFALSE);  // with this line the histogram taken with TFile::Get() is no more associated to the file (no need to Clone anymore)
  // note that even though I clone the object, if the file is closed the clone is lost, as it happens for the object retrieved with TFile::Get()

  TH1F *htmp = NULL;

  cout << "Opening file " << fileName << endl;

  TFile* f = TFile::Open(fileName.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<fileName<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  // get E/P template from file
  htmp = (TH1F*)f->Get("hEoP_template");  
  if (!htmp) {
    cout << "Error: histogram " << htmp->GetName() << " not found in file ' " << fileName << "'. End of programme." << endl;
    exit(EXIT_FAILURE);
  }
  htmp->SetDirectory(0);
  //  hinput = (TH1F*) htmp->Clone();
  hinput = new TH1F(*htmp);  // it looks like that after exiting the function hinput is gone
  if (!hinput) {
    cout << "Error in function getEoP_templateMC(TH1F * hinput, const string& dirName): hinput is NULL. End of programme." << endl;
    exit(EXIT_FAILURE);
  }

  // do NOT close the file, or hinput is lost 
  //  f->Close();  

}


//=====================================================================


void getLowAndUpBinGivenRange(const TH1F* histo, const Double_t &lowerX, const Double_t &upperX, Int_t &lowBin, Int_t &upBin) {
  lowBin = histo->GetXaxis()->FindFixBin(lowerX);
  upBin = histo->GetXaxis()->FindFixBin(upperX);
}

Int_t getLowBinGivenRange(const TH1F* histo, const Double_t &lowerX) {
  Int_t lowBin = histo->GetXaxis()->FindFixBin(lowerX);
  return lowBin;
}

Int_t getUpBinGivenRange(const TH1F* histo, const Double_t &upperX) {
  Int_t upBin = histo->GetXaxis()->FindFixBin(upperX);
  return upBin;
}

//=====================================================================

Double_t myCrystalBallRightTail(double* x, double* par) {

  Double_t xcur = x[0];
  Double_t alpha = par[0];
  Double_t n = par[1];
  Double_t mu = par[2];
  Double_t sigma = par[3];
  Double_t N = par[4];
  Double_t t = (xcur-mu)/sigma;
  Double_t absAlpha = fabs((Double_t)alpha);
  Double_t invAbsAlpha = 1./absAlpha;

  if ( t <= absAlpha)  {   // would be t >= -absAlpha for left tail
    return N*exp(-0.5*t*t);
  } else {
    Double_t A = TMath::Power(n*invAbsAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    Double_t B = n*invAbsAlpha - absAlpha;
    return N*A * TMath::Power(B+t,-n);  // would be B-t for left tail
  }

}

//=====================================================================

Double_t myCrystalBallLeftTail(double* x, double* par) {

  Double_t xcur = x[0];
  Double_t alpha = par[0];
  Double_t n = par[1];
  Double_t mu = par[2];
  Double_t sigma = par[3];
  Double_t N = par[4];
  Double_t t = (xcur-mu)/sigma;
  Double_t absAlpha = fabs((Double_t)alpha);
  Double_t invAbsAlpha = 1./absAlpha;

  if ( t >= -absAlpha)  {   // would be t <= absAlpha for right tail
    return N*exp(-0.5*t*t);
  } else {
    Double_t A = TMath::Power(n*invAbsAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    Double_t B = n*invAbsAlpha - absAlpha;
    return N*A * TMath::Power(B-t,-n);  // would be B+t for right tail
  }

}

//=====================================================================

Double_t my2sideCrystalBall(double* x, double* par) {

  //a priori we allow for different shape of right and left tail, thus two values of alpha and n 

  Double_t xcur = x[0];
  Double_t alphaL = par[0];
  Double_t nL = par[1];
  Double_t mu = par[2];
  Double_t sigma = par[3];
  Double_t N = par[4];
  Double_t alphaR = par[5];
  Double_t nR = par[6];
  Double_t t = (xcur-mu)/sigma;
  Double_t absAlphaL = fabs((Double_t)alphaL);
  Double_t invAbsAlphaL = 1./absAlphaL;
  Double_t absAlphaR = fabs((Double_t)alphaR);
  Double_t invAbsAlphaR = 1./absAlphaR;

  
  if ( t<-absAlphaL ) {
    //cout<<"checkpoint dscb left"<<endl;
    Double_t AL = TMath::Power(nL*invAbsAlphaL,nL)*exp(-0.5*absAlphaL*absAlphaL);
    Double_t BL = nL*invAbsAlphaL - absAlphaL;
    return N*AL*TMath::Power(BL-t,-nL);
  } else if ( t <= absAlphaR )  {
    //cout<<"checkpoint dscb gaussian"<<endl;
    return N*exp(-0.5*t*t);
  } else {
    //cout<<"checkpoint dscb right"<<endl;
    Double_t AR = TMath::Power(nR*invAbsAlphaR,nR)*exp(-0.5*absAlphaR*absAlphaR);
    Double_t BR = nR*invAbsAlphaR - absAlphaR;
    return N*AR*TMath::Power(BR+t,-nR);
  }

}

//=====================================================================

void setHistColor(vector<Int_t> &histColor, const Int_t nObject) {
  
  Int_t colorList[] = {kBlue, kRed, kGreen, kOrange+1, kCyan, kViolet};  
  // the first color is for the main object. This array may contain more values than nSamples  
  Int_t vectorIsEmpty = ((Int_t)histColor.size()) == 0 ? 1 : 0;
                                                          
  for (Int_t i = 0; i < nObject; i++) {   // now color are assigned in reverse order (the main contribution is the last object in the sample array)         

    if (vectorIsEmpty) histColor.push_back(colorList[i]);
    else histColor[i] = colorList[i];    // if histColor is just uninitialized but not really empty

  }
  
}


void buildChainWithFriend(TChain* chain, TChain* chFriend, string sampleName, TChain* chsfFriend = NULL) {
  
  cout << "Creating chain ..." << endl;
  
  vector<string> subSampleNameVector;

  if (sampleName == "DATA") {

    // 2016 2.6fb^-1
    if (DATA2016) {
      if (SKIM_1TLEP30) {
	subSampleNameVector.push_back("SingleElectron_Run2016B_PromptReco_v1_runs_272023_273146");
	subSampleNameVector.push_back("SingleElectron_Run2016B_PromptReco_v2_runs_273150_275376");
	subSampleNameVector.push_back("SingleElectron_Run2016C_PromptReco_v2_runs_275420_276283");
	subSampleNameVector.push_back("SingleElectron_Run2016D_PromptReco_v2_runs_276315_276811");
	subSampleNameVector.push_back("SingleElectron_Run2016E_PromptReco_v2_runs_276830_277420");
	subSampleNameVector.push_back("SingleElectron_Run2016F_PromptReco_v1_runs_277820_278808");
	subSampleNameVector.push_back("SingleElectron_Run2016G_PromptReco_v1_runs_278817_279931");	
	subSampleNameVector.push_back("SingleElectron_Run2016H_PromptReco_v1_runs_281085_281201");	
      } else if (SKIM_1LEP1JET_80X) {
	subSampleNameVector.push_back("SingleElectron_Run2016B_PromptReco_v1_runs_272023_273146");
	subSampleNameVector.push_back("SingleElectron_Run2016B_PromptReco_v2_runs_273150_275376");
	subSampleNameVector.push_back("SingleElectron_Run2016C_PromptReco_v2_runs_275420_276283");
	subSampleNameVector.push_back("SingleElectron_Run2016D_PromptReco_v2_runs_276315_276811");
	subSampleNameVector.push_back("SingleElectron_Run2016E_PromptReco_v2_runs_276830_277420");
	subSampleNameVector.push_back("SingleElectron_Run2016F_PromptReco_v1_runs_277820_278808");
	subSampleNameVector.push_back("SingleElectron_Run2016G_PromptReco_v1_runs_278817_279931");	
      } else {
	subSampleNameVector.push_back("SingleElectron_PromptReco_v1_runs_272021_273149");
	subSampleNameVector.push_back("SingleElectron_PromptReco_v2_runs_273150_274443");
      }
    } else {
    // 2015 2.32 fb^-1
      subSampleNameVector.push_back("SingleElectron_Run2015C_16Dec_runs_254227_254914");
      subSampleNameVector.push_back("SingleElectron_Run2015D_16Dec_runs_256630_260627");
    }

  } else if (sampleName == "WJetsToLNu") {

    if (SKIM_1TLEP30) {
      if (MC_NLO_NOHTBIN) {
	subSampleNameVector.push_back("WJetsToLNu");	
      } else {
	subSampleNameVector.push_back("WJetsToLNu_HT100to200");
	//subSampleNameVector.push_back("WJetsToLNu_HT100to200_ext");
	subSampleNameVector.push_back("WJetsToLNu_HT200to400");
	//subSampleNameVector.push_back("WJetsToLNu_HT200to400_ext");
	subSampleNameVector.push_back("WJetsToLNu_HT400to600");
	//subSampleNameVector.push_back("WJetsToLNu_HT400to600_ext");
	subSampleNameVector.push_back("WJetsToLNu_HT600to800");
	subSampleNameVector.push_back("WJetsToLNu_HT800to1200");
	//subSampleNameVector.push_back("WJetsToLNu_HT800to1200_ext");
	subSampleNameVector.push_back("WJetsToLNu_HT1200to2500");
	subSampleNameVector.push_back("WJetsToLNu_HT2500toInf");
      }
    } else if (SKIM_1LEP1JET_80X) {
      subSampleNameVector.push_back("WJetsToLNu_HT100to200");
      //subSampleNameVector.push_back("WJetsToLNu_HT100to200_ext");
      subSampleNameVector.push_back("WJetsToLNu_HT200to400");
      //subSampleNameVector.push_back("WJetsToLNu_HT200to400_ext");
      subSampleNameVector.push_back("WJetsToLNu_HT400to600");
      //subSampleNameVector.push_back("WJetsToLNu_HT400to600_ext");
      subSampleNameVector.push_back("WJetsToLNu_HT600to800");
      subSampleNameVector.push_back("WJetsToLNu_HT800to1200");
      //subSampleNameVector.push_back("WJetsToLNu_HT800to1200_ext");
      subSampleNameVector.push_back("WJetsToLNu_HT1200to2500");
      subSampleNameVector.push_back("WJetsToLNu_HT2500toInf");
    } else {    
      subSampleNameVector.push_back("WJetsToLNu_HT100to200");
      subSampleNameVector.push_back("WJetsToLNu_HT200to400");
      subSampleNameVector.push_back("WJetsToLNu_HT400to600");
      subSampleNameVector.push_back("WJetsToLNu_HT600to800");
      subSampleNameVector.push_back("WJetsToLNu_HT800to1200");
      if (!DATA2016) subSampleNameVector.push_back("WJetsToLNu_HT1200to2500"); // note, 1200to2500 bin missing for 2016 MC
      subSampleNameVector.push_back("WJetsToLNu_HT2500toInf"); 
    }

  } else {

    cout << "Error: unknown sampleName " << sampleName <<". End of programme" << endl;
    exit(EXIT_FAILURE);

  }

  //2016 trees
  string treePath = "";
  if (SKIM_1TLEP30){
    if (MC_NLO_NOHTBIN) treePath = "/afs/cern.ch/work/m/mciprian/EoverP_study_new/CMSSW_8_0_10/src/eleEoverP/TREES_1tightEle30_skim/";
    else treePath = "root://eoscms//eos/cms/store/group/phys_exotica/monojet/mciprian/trees_80X/TREES_1TIGHTELE30SKIM_4EoP/";
  else if (SKIM_1LEP1JET_80X) treePath = "root://eoscms//eos/cms/store/group/phys_exotica/monojet/mciprian/trees_80X/TREES_1TLEP1JET_80X_4EoP/"; // 2016 trees with skim 1 lep1jet
  else {
    if (DATA2016) treePath = "root://eoscms//eos/cms/store/cmst3/group/susy/emanuele/monox/trees/TREES_1LEPSKIM_80X/"; // 2016 trees
    else treePath = "root://eoscms//eos/cms/store/cmst3/group/susy/emanuele/monox/trees/TREES_25ns_1LEPSKIM_76X/";   //2015 trees
  }

  for(UInt_t i = 0; i < subSampleNameVector.size(); i++) {
  
    string treeRootFile = treePath + subSampleNameVector[i] + "_treeProducerDarkMatterMonoJet_tree.root";
    string friend_treeRootFile = "";
    string sffriend_treeRootFile = "";
    if (SKIM_1LEP1JET_80X || SKIM_1TLEP30) {
      friend_treeRootFile = treePath + "evVarFriend_" + subSampleNameVector[i]+ ".root";
      if (SKIM_1TLEP30) sffriend_treeRootFile = treePath + "sfFriend_" + subSampleNameVector[i]+ ".root";
    } else {
      if (DATA2016) friend_treeRootFile = treePath + "friends_evVarFriend_" + subSampleNameVector[i]+ ".root";
      else friend_treeRootFile = treePath + "evVarFriend_" + subSampleNameVector[i]+ ".root";
    }

    chain->Add(TString(treeRootFile.c_str()));
    chFriend->Add(TString(friend_treeRootFile.c_str()));
    if (SKIM_1TLEP30) chsfFriend->Add(TString(sffriend_treeRootFile.c_str()));

  }

  cout << "Adding friend to chain ..." << endl;
  chain->AddFriend(chFriend);  //adding whole friend chain as friend                                                           
  if (SKIM_1TLEP30 && chsfFriend != NULL) chain->AddFriend(chsfFriend);  //adding whole friend chain as friend                

  if(!chain || !chFriend) {
    cout << "Error: chain not created. End of programme" << endl;
    exit(EXIT_FAILURE);
  }

  cout << chain->GetEntries() << endl;

}

//============================================================   

Int_t getBinNumber(const Float_t value, const vector<Float_t> &binEdgesVector) {

  // return invalid value if bin not found
  // return -1 if value > binEdgesVector.last(),
  // return -2 if value < binEdgesVector.first(),

  if (value > binEdgesVector[0]) {

    Int_t bin = 0;
    Int_t howManyBins = binEdgesVector.size() - 1;

    while (bin < howManyBins) {
      if (value < binEdgesVector[bin+1]) {
	return bin;
      }
      bin++;
    }      

    return -1;

  } else return -2;

}

//============================================================

void EoverP::Loop(const string sampleName, const vector<Float_t> &energybinEdges, const string& dirName, const Int_t reweightMC = 1)
{


   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);

   fChain->SetBranchStatus("nEle10V",1);  // # of electrons passing loose selection for electron veto
   fChain->SetBranchStatus("nEle40T",1);
   fChain->SetBranchStatus("nLepGood",1);
   fChain->SetBranchStatus("LepGood_pdgId",1);  // must be 13 for muons ( -13 for mu+), 11 for electrons and 15 for taus                                               
   fChain->SetBranchStatus("LepGood_pt",1);
   fChain->SetBranchStatus("LepGood_eta",1);
   fChain->SetBranchStatus("LepGood_correctedEcalEnergy",1);
   fChain->SetBranchStatus("LepGood_superCluster_rawEnergy",1);
   fChain->SetBranchStatus("LepGood_eSuperClusterOverP",1);
   fChain->SetBranchStatus("LepGood_full5x5_r9",1);
   fChain->SetBranchStatus("met_pt",1);
   fChain->SetBranchStatus("LepGood_gsfTrackP",1);

   // for additional studies
   fChain->SetBranchStatus("LepGood_full5x5_sigmaIetaIeta",1);

   // branches for MC study
   if (sampleName != "DATA"){

     fChain->SetBranchStatus("ngenLep",1);
     fChain->SetBranchStatus("genLep_motherId",1);
     fChain->SetBranchStatus("genLep_pdgId",1);
     fChain->SetBranchStatus("genLep_pt",1);
     fChain->SetBranchStatus("genLep_eta",1);
     fChain->SetBranchStatus("genLep_phi",1);
     fChain->SetBranchStatus("genLep_mass",1);

     fChain->SetBranchStatus("weight",1);
     fChain->SetBranchStatus("puw",1);
     fChain->SetBranchStatus("SF_BTag",1);

     fChain->SetBranchStatus("SF_trig1lep",1);
     fChain->SetBranchStatus("SF_LepTight",1);
     fChain->SetBranchStatus("SF_NLO_QCD",1);
     fChain->SetBranchStatus("SF_NLO_EWK",1);

   }

   gStyle->SetStatStyle(0);

   // passing from main
   // vector<Float_t> energybinEdges;
   // energybinEdges.push_back(25.0);
   // energybinEdges.push_back(50.0);
   // energybinEdges.push_back(75.0);
   // energybinEdges.push_back(100.0);
   // energybinEdges.push_back(150.0);
   // energybinEdges.push_back(200.0);
   // energybinEdges.push_back(250.0);
   // energybinEdges.push_back(300.0);
   // energybinEdges.push_back(350.0);
   // energybinEdges.push_back(450.0);
   // energybinEdges.push_back(600.0);
   // energybinEdges.push_back(750.0);
   // energybinEdges.push_back(900.0);

   TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 

   Int_t nEnergyBins = energybinEdges.size() -1;

   string rootfileName = dirName + "EoverP_" + sampleName + ".root";

   TFile *rootFile = new TFile((rootfileName).c_str(),"RECREATE");
   if (!rootFile || !rootFile->IsOpen()) {
     cout << "Error: file \"" << rootfileName << "\" was not opened." << endl;
     exit(EXIT_FAILURE);
   }


   string fwgtname = dirName;
   if (USE_E) fwgtname += "hEcorrWeight_dataMCratio.root";
   else fwgtname += "hETcorrWeight_dataMCratio.root";

   TFile* fdataMCratioWeight = NULL;
   TH1F *hwgt = NULL;
   if (reweightMC && sampleName != "DATA") {
     fdataMCratioWeight =  TFile::Open(fwgtname.c_str(),"READ");
     if (!fdataMCratioWeight || !fdataMCratioWeight->IsOpen()) {
       cout<<"*******************************"<<endl;
       cout<<"Error: could not find file to reweight MC based on data/MC ratio for Ecorr distribution."<<endl;
       cout<<"*******************************"<<endl;
       exit(EXIT_FAILURE);
     }
     hwgt = (TH1F*) fdataMCratioWeight->Get("hdataMCratio");
     //hwgt->SetDirectory(0);
     if (!hwgt) {
       cout << "Error: histogram not found in file '" << fwgtname << "'. End of programme." << endl;
       exit(EXIT_FAILURE);
     }
   } else if (sampleName != "DATA") {
     cout << "This time MC will not be reweighted" << endl;
   }

   rootFile->cd(); // get back to this file

   vector<TH1F*> hEoverP_energyBin(nEnergyBins,NULL);
   // histogram to store the mean energy in a given bin, useful to plot E/P asf of the mean energy in the bin with a TGraph
   TH1F* hMeanEnergyInEnergyBin = new TH1F("hMeanEnergyInEnergyBin","",nEnergyBins,energybinEdges.data());

   vector<TH1F*> hEcorrOverEtrue_energyBin(nEnergyBins,NULL);
   vector<TH1F*> hErawOverEtrue_energyBin(nEnergyBins,NULL);
   vector<TH1F*> hPtrackOverEtrue_energyBin(nEnergyBins,NULL);

   // for additional studies
   TH1F* hSigmaIetaIeta_lowE = new TH1F("hSigmaIetaIeta_lowE","",96,0.0,0.012);
   TH1F* hSigmaIetaIeta_hotE = new TH1F("hSigmaIetaIeta_hotE","",96,0.0,0.012);
   TH1F* hR9_lowE = new TH1F("hR9_lowE","",90,0.2,1.1);
   TH1F* hR9_hotE = new TH1F("hR9_hotE","",90,0.2,1.1);
   TH1F* hPt_lowE = new TH1F("hPt_lowE","",140,0.0,700.0);
   TH1F* hPt_hotE = new TH1F("hPt_hotE","",140,0.0,700.0);   

   TH1F* hElePt = new TH1F("hElePt","",200,0,1000);
   TH1F* hEleEcorr = new TH1F("hEleEcorr","",200,0,1000);
   TH1F* hEleEraw = new TH1F("hEleEraw","",200,0,1000);

   for (Int_t i = 0; i < nEnergyBins; i++) {
     hEoverP_energyBin[i] = new TH1F(Form("hEoverP_energyBin%1.0fTo%1.0f",energybinEdges[i],energybinEdges[i+1]),"",200,0.05,2.05);
   }

   if (sampleName != "DATA") {
     for (Int_t i = 0; i < nEnergyBins; i++) {
       hEcorrOverEtrue_energyBin[i] = new TH1F(Form("hEcorrOverEtrue_energyBin%1.0fTo%1.0f",energybinEdges[i],energybinEdges[i+1]),"",200,0.55,1.55);
       hErawOverEtrue_energyBin[i] = new TH1F(Form("hErawOverEtrue_energyBin%1.0fTo%1.0f",energybinEdges[i],energybinEdges[i+1]),"",200,0.55,1.55);
       hPtrackOverEtrue_energyBin[i] = new TH1F(Form("hPtrackOverEtrue_energyBin%1.0fTo%1.0f",energybinEdges[i],energybinEdges[i+1]),"",200,0.05,2.05);
     }
   }

   TH1F* hEoP_template = new TH1F("hEoP_template","template",200,0.05,2.05);  // a template for E/P in the whole energy spectrum (not used at the moment)

   // Float_t lowMeanCut = 0.85;
   // Float_t upMeanCut = 1.15;

   // used only for MC study
   /////////////////////////
   Int_t MCtruthMatchFound = 0;
   Float_t Etrue = 0.0; 
   /////////////////////////

   Double_t wgt = 1.0;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
       
     Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%500000 == 0) cout << jentry << endl;
  
      // if (met_pt < 50) continue;
      if (!( nEle10V == 1 && nEle40T == 1) ) continue;
      if (!( fabs(LepGood_pdgId[0]) == 11 && LepGood_pt[0] > 30 && fabs(LepGood_eta[0]) < ELE_ETA_MAX) ) continue;

      Double_t energyToUse = -1.0;
      if (USE_RAWE) energyToUse = LepGood_superCluster_rawEnergy[0];
      else energyToUse = LepGood_correctedEcalEnergy[0];

      if (!USE_E) {
	// could use ET = E /cosh(eta), maybe it's faster
	Double_t theta = 2. * atan(exp(-LepGood_eta[0])); 
	energyToUse *= sin(theta);
      } 

      if (sampleName != "DATA") {
	if (reweightMC) wgt = (Double_t) hwgt->GetBinContent(hwgt->FindBin(energyToUse));
	else if (NO_MC_WEIGHT) wgt = 1.0;
	else {
	  wgt = 24.4 * weight * puw * SF_NLO_QCD * SF_NLO_EWK * SF_trig1lep * SF_LepTight;
	  if (!MC_NLO_NOHTBIN) wgt /= 1.21;
	}
      }
      else wgt = 1.0;


      // fill these before cutting on r9
      if (energyToUse > 99.9 && energyToUse < 225.0) {
	hR9_lowE->Fill(LepGood_full5x5_r9[0],wgt);
      } else if (energyToUse > 349.9 && energyToUse < 450.0) {
	hR9_hotE->Fill(LepGood_full5x5_r9[0],wgt);
      }

      if (!(LepGood_full5x5_r9[0] > 0.94)) continue;

      if (energyToUse > 99.9 && energyToUse < 225.0) {
	hSigmaIetaIeta_lowE->Fill(LepGood_full5x5_sigmaIetaIeta[0],wgt);
	hPt_lowE->Fill(LepGood_pt[0],wgt);
      } else if (energyToUse > 349.9 && energyToUse < 450.0) {
	hSigmaIetaIeta_hotE->Fill(LepGood_full5x5_sigmaIetaIeta[0],wgt);
	hPt_hotE->Fill(LepGood_pt[0],wgt);
      }

      hElePt->Fill(LepGood_pt[0],wgt);
      hEleEcorr->Fill(LepGood_correctedEcalEnergy[0],wgt);  
      hEleEraw->Fill(LepGood_superCluster_rawEnergy[0],wgt); 

      Double_t EoverP_toUse = -1.0;
      // LepGood_eSuperClusterOverP[0] is E_regrCorr/P_track. If we need E_raw/P_track then divide by E_egrCorr and multiply by E_raw
      if (USE_RAWE) EoverP_toUse = LepGood_superCluster_rawEnergy[0] / LepGood_gsfTrackP[0];
      else EoverP_toUse = LepGood_correctedEcalEnergy[0] / LepGood_gsfTrackP[0];


      //if (!(fabs(LepGood_pdgId[0]) == 11 && LepGood_pt[0] > 40 && fabs(LepGood_eta[0]) > 1.0 && fabs(LepGood_eta[0]) < 1.479)) continue;

      // here goes with the algorithm to match gen to reco electrons and in case fill histograms

      if (sampleName != "DATA") {

	MCtruthMatchFound = 0;  // reset to 0 
	Etrue = 0.0;            // reset this too
	Int_t i = 0;

	while (!MCtruthMatchFound && i < ngenLep) {

	  // start asking for e+/- (|pdgID| = 11) coming from a W+/- (|pdgID| = 24)
	    
	  if (fabs(genLep_pdgId[i]) == 11 && fabs(genLep_motherId[i]) == 24) {

	    TVector3 eleReco;
	    eleReco.SetPtEtaPhi(LepGood_pt[0],LepGood_eta[0],LepGood_phi[0]);
	    TVector3 eleGen;
	    eleGen.SetPtEtaPhi(genLep_pt[i],genLep_eta[i],genLep_phi[i]);

	    // now match is ok if deltaR < 0.3
	    // in this case, compute Etrue and set MCtruthMatchFound flag as 1

	    if (eleGen.DeltaR(eleReco) < 0.3) {

	      Etrue = eleGen.Mag();  // neglect mass since we have electrons
	      MCtruthMatchFound = 1;  

	    }
	    
	  }

	  i++;

	}

      } 

      // end of loop to match gen and reco electrons

      hEoP_template->Fill(EoverP_toUse,wgt);

      // look for the bin in the LepGood_correctedEcalEnergy variable      
      Int_t bin = getBinNumber(energyToUse,energybinEdges);  // this function returns negative value if bin not found
      
      if (bin >= 0) {

	hEoverP_energyBin[bin]->Fill(EoverP_toUse,wgt);

	//sum the energy to the bin content in the bin it belongs to (at the end we will divide by the number of entries in each bin)
	// using bin+1 because the histogram bin number goes from 1 to number of bins, while "bin" variable starts from 0
	hMeanEnergyInEnergyBin->SetBinContent(bin+1, energyToUse + hMeanEnergyInEnergyBin->GetBinContent(bin+1));  

      } else if (bin == -1) {
        // fill last bin with overflows to gain in statistics                                           
        hEoverP_energyBin[nEnergyBins-1]->Fill(EoverP_toUse,wgt);
      }



      Int_t etrueBin = getBinNumber(Etrue,energybinEdges);

      if (etrueBin >= 0) {

	if (MCtruthMatchFound) {
	  
	  Double_t invEtrue = 1./Etrue;
	  hEcorrOverEtrue_energyBin[etrueBin]->Fill(LepGood_correctedEcalEnergy[0] * invEtrue);
	  hErawOverEtrue_energyBin[etrueBin]->Fill(LepGood_superCluster_rawEnergy[0] * invEtrue);
	  hPtrackOverEtrue_energyBin[etrueBin]->Fill(LepGood_gsfTrackP[0] * invEtrue);    
	  // Ptrack/Etrue = corrE * (EoverP)^-1 / Etrue
	}
      }


   }  // end of loop on entries

   for (Int_t i = 0; i < nEnergyBins; i++) {

     //hEoverP_energyBin[i]->SetStats(0);  // keep stat box in file, and remove in directly the function that plots it on canvas
     // keep the following axis settings
     hEoverP_energyBin[i]->GetXaxis()->SetTitle("E / P");
     hEoverP_energyBin[i]->GetXaxis()->SetTitleSize(0.06);
     hEoverP_energyBin[i]->GetXaxis()->SetTitleOffset(0.8);
     hEoverP_energyBin[i]->GetYaxis()->SetTitle("events");
     hEoverP_energyBin[i]->GetYaxis()->SetTitleSize(0.055);
     hEoverP_energyBin[i]->GetYaxis()->SetTitleOffset(0.8);

     hMeanEnergyInEnergyBin->SetBinContent(i+1, hMeanEnergyInEnergyBin->GetBinContent(i+1)/hEoverP_energyBin[i]->GetEntries());

   }

   string whichEnergy = "";
   if (USE_RAWE) whichEnergy = "raw";
   else whichEnergy = "corrected";

   if (USE_E) {
     hMeanEnergyInEnergyBin->GetXaxis()->SetTitle((whichEnergy + " E [GeV]").c_str());
     hMeanEnergyInEnergyBin->GetYaxis()->SetTitle("mean E [GeV]");
   } else {
     hMeanEnergyInEnergyBin->GetXaxis()->SetTitle((whichEnergy + " E_{T} [GeV]").c_str());
     hMeanEnergyInEnergyBin->GetYaxis()->SetTitle("mean E_{T} [GeV]");
   }
   hMeanEnergyInEnergyBin->GetXaxis()->SetTitleSize(0.06);
   hMeanEnergyInEnergyBin->GetXaxis()->SetTitleOffset(0.8);
   hMeanEnergyInEnergyBin->GetYaxis()->SetTitleSize(0.055);
   hMeanEnergyInEnergyBin->GetYaxis()->SetTitleOffset(0.8);

   if (sampleName != "DATA") {

     for (Int_t i = 0; i < nEnergyBins; i++) {

       hEcorrOverEtrue_energyBin[i]->GetXaxis()->SetTitle("E_{corr} / E_{true}");
       hEcorrOverEtrue_energyBin[i]->GetXaxis()->SetTitleSize(0.06);
       hEcorrOverEtrue_energyBin[i]->GetXaxis()->SetTitleOffset(0.75);
       hEcorrOverEtrue_energyBin[i]->GetYaxis()->SetTitle("events");
       hEcorrOverEtrue_energyBin[i]->GetYaxis()->SetTitleSize(0.055);
       hEcorrOverEtrue_energyBin[i]->GetYaxis()->SetTitleOffset(0.8);

       hErawOverEtrue_energyBin[i]->GetXaxis()->SetTitle("E_{raw} / E_{true}");
       hErawOverEtrue_energyBin[i]->GetXaxis()->SetTitleSize(0.06);
       hErawOverEtrue_energyBin[i]->GetXaxis()->SetTitleOffset(0.75);
       hErawOverEtrue_energyBin[i]->GetYaxis()->SetTitle("events");
       hErawOverEtrue_energyBin[i]->GetYaxis()->SetTitleSize(0.055);
       hErawOverEtrue_energyBin[i]->GetYaxis()->SetTitleOffset(0.8);

       hPtrackOverEtrue_energyBin[i]->GetXaxis()->SetTitle("P_{track} / E_{true}");
       hPtrackOverEtrue_energyBin[i]->GetXaxis()->SetTitleSize(0.06);
       hPtrackOverEtrue_energyBin[i]->GetXaxis()->SetTitleOffset(0.75);
       hPtrackOverEtrue_energyBin[i]->GetYaxis()->SetTitle("events");
       hPtrackOverEtrue_energyBin[i]->GetYaxis()->SetTitleSize(0.055);
       hPtrackOverEtrue_energyBin[i]->GetYaxis()->SetTitleOffset(0.8);
    
     }

   }

   rootFile->Write();
   rootFile->Close();
   delete rootFile;


}

//=========================================================================

void getTexMCSampleName(const string &MCSampleName, string &texMCSampleName) {
  
  if (MCSampleName == "WJetsToLNu") texMCSampleName = "W(l#nu)+jets";
  else if (MCSampleName == "DYJetsToLL") texMCSampleName = "Z(ll)+jets";

}

//=========================================================================

void plotDistribution(const string &sampleName, const vector<Float_t> &energybinEdges, TH1F* hPeak, TH1F* hSigma, const string &hNameID, const string &dirName, TH1F* hPeakMeanInRange = NULL, const Double_t lowRangeValue = 0.0, const Double_t upRangeValue = 2.0, TH1F* hPeakShift_MCdata = NULL, TH1F* hScaleFactor_MCdata = NULL) {

  cout << endl;

  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2() 
  TVirtualFitter::SetDefaultFitter("Minuit");

  string fileName = dirName + "EoverP_" + sampleName + ".root";

  TCanvas *c = new TCanvas("c",""); 

  TH1F* htmp = NULL;
  TH1F* hist = NULL;    

  TFile* f = TFile::Open(fileName.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<fileName<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  TH1F* hEoP_template = NULL;
  // get E/P template from MC file

  // -----------------------------------------------------------
  // get MC template. If running on data, open MC file, else the file is already the right one and I just need to get histograms

  TFile* fMC = NULL;  // to get MC file when running on data
  string fileNameMC = dirName + "EoverP_WJetsToLNu.root";  // template is created when running on MC, and is stored in this file

  if (sampleName == "DATA") {

    //    getEoP_templateMC(hEoP_template, dirName);
    cout << "Running on data. Opening file " << fileNameMC << " to get template." << endl;

    fMC = TFile::Open(fileNameMC.c_str(),"READ");
    if (!fMC || !fMC->IsOpen()) {
      cout<<"*******************************"<<endl;
      cout<<"Error opening file \""<<fileNameMC<<"\".\nApplication will be terminated."<<endl;
      cout<<"*******************************"<<endl;
      exit(EXIT_FAILURE);
    }
    
  //   // get E/P template from file
  //   htmp = (TH1F*)fMC->Get("hEoP_template");  
  //   if (!htmp) {
  //     cout << "Error: histogram " << htmp->GetName() << " not found in file ' " << fileNameMC << "'. End of programme." << endl;
  //     exit(EXIT_FAILURE);
  //   }
  //   hEoP_template = (TH1F*) htmp->Clone();
    
  } else {

  //   htmp = (TH1F*)f->Get("hEoP_template");  
  //   if (!htmp) {
  //     cout << "Error: histogram " << htmp->GetName() << " not found in file ' " << fileName << "'. End of programme." << endl;
  //     exit(EXIT_FAILURE);
  //   }
  //   hEoP_template = (TH1F*) htmp->Clone();

  // }

  // if (!hEoP_template) {
  //   cout << "Error: hEoP_template is NULL! End of programme." << endl;
  //   exit(EXIT_FAILURE);
  }

  // end of access of MC file to get templates
  // -----------------------------------------------------------

  UInt_t nBins = energybinEdges.size() -1;

  for (UInt_t i = 0; i < nBins; i++) {

    cout << " Energy bin: [" << energybinEdges[i] << ", " << energybinEdges[i+1] << "]" << endl;

    htmp = (TH1F*)f->Get(Form("h%s_energyBin%1.0fTo%1.0f",hNameID.c_str(),energybinEdges[i],energybinEdges[i+1]));  
    if (!htmp) {
      cout << "Error: histogram " << Form("h%s_energyBin%1.0fTo%1.0f",hNameID.c_str(),energybinEdges[i],energybinEdges[i+1]) 
	   << " not found in file ' " << fileName << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
    // now, if using data I will fit with template obtained from the same distribution in MC file. Since the histogram's name is the same in both file, I 
    // got a strange behaviour, that is, when using fMC->Get(<name>) I actually get the same histogram in data. A solution is deleting the previous histogram
    // got from data file before getting that for MC. In this case, better to give the clone a different name to make explicit it is for data 
    if (sampleName == "DATA") hist = (TH1F*)htmp->Clone(Form("%s_data",htmp->GetName()));
    else hist = (TH1F*)htmp->Clone();

    //    hist->SetStats(0);
    gStyle->SetOptStat(10);
    gStyle->SetOptFit(112);
    hist->Draw("HE");
    if (sampleName == "DATA") {
      if (energybinEdges[i] > 449.9) hist->Rebin(5); 
      else if (energybinEdges[i] > 349.9) hist->Rebin(4);
      else if (energybinEdges[i] > 274.9) hist->Rebin(3); 
      else if (energybinEdges[i] > 174.9) hist->Rebin(2); 
      //else if (energybinEdges[i] > 199.9) hist->Rebin(2); // when using E instead of Et 249.9 is ok
    } else if (hNameID != "EoverP") {
      if (hNameID == "PtrackOverEtrue") {
	if (energybinEdges[i] > 449.9) hist->Rebin(5); 
	else if (energybinEdges[i] > 349.9) hist->Rebin(4); 
	else if (energybinEdges[i] > 274.9) hist->Rebin(3); 
	else if (energybinEdges[i] > 174.9) hist->Rebin(2); 
      } else {
	hist->GetXaxis()->SetRangeUser(0.85,1.15);
      }
    } else {
      if (energybinEdges[i] > 549.9) hist->Rebin(4); 
      else if (energybinEdges[i] > 349.9) hist->Rebin(3); 
      else if (energybinEdges[i] > 224.9) hist->Rebin(2);      
    }
    c->Update();

    // fitting:
    // do a first fit with a simple gaussian in the core
    Double_t gaussEdgeL = 0.9;  //left side of the gaussian to be used in the fit (I use a variable so that I change this value only once)
    Double_t gaussEdgeR = 1.1;  //right side ...
    // some general settings
    if (energybinEdges[i] > 349.9) {
      gaussEdgeL = 0.75;
      gaussEdgeR = 1.25;
    }
    // if MC, use more specific settings (decided after looking at plots, there is no a priori reason for them)
    if (sampleName != "DATA") {
      if (hNameID != "EoverP") {
	if (hNameID == "PtrackOverEtrue") {
	  gaussEdgeL = 0.9;
	  gaussEdgeR = 1.1;
	} else {
	  gaussEdgeL = 0.95;
	  gaussEdgeR = 1.05;	
	}
      } else {
	if (energybinEdges[i] > 349.9) {
	  gaussEdgeL = 0.75;
	  gaussEdgeR = 1.25;
	}
      }
    }
    hist->Fit("gaus","E L I Q 0","",gaussEdgeL,gaussEdgeR);  // L: loglikelihood method, 0: do not plot this fit, Q: quiet mode (minimum printing)
    Double_t gaussNorm = hist->GetFunction("gaus")->GetParameter(0);
    Double_t gaussMean = hist->GetFunction("gaus")->GetParameter(1);
    //Double_t gaussMeanError = hist->GetFunction("gaus")->GetParError(1);
    Double_t gaussSigma = hist->GetFunction("gaus")->GetParameter(2);
    // now use crystal ball with right tail
    Double_t funcRangeLeft = gaussEdgeL;
    Double_t funcRangeRight = 2.0;
    if (FIT_2SIDE_CB) {
      if (hNameID != "EoverP") {
	if (hNameID != "PtrackOverEtrue") {
	  funcRangeLeft = 0.8;
	  funcRangeRight = gaussEdgeR;
	} else {
	  funcRangeLeft = 0.2;
	  funcRangeRight = 1.3;
	}
      } else {
	funcRangeLeft = 0.5;
      }
    }
    else {
      if (hNameID != "EoverP") {
	if (hNameID != "PtrackOverEtrue") funcRangeLeft = 0.8;
	else funcRangeLeft = 0.2;
	funcRangeRight = gaussEdgeR;
      }
    }

    TF1 *cb1;

    if (FIT_2SIDE_CB) {

      cb1 = new TF1("cb1",&my2sideCrystalBall,funcRangeLeft,funcRangeRight,7);
      cb1->SetParNames("alphaL","nL","mu","sigma","N","alphaR","nR");  
      cb1->SetParLimits(cb1->GetParNumber("nL"),0.1,15);
      cb1->SetParLimits(cb1->GetParNumber("nR"),0.1,15);
      cb1->SetParLimits(cb1->GetParNumber("alphaL"),-10.,-0.01); 
      cb1->SetParLimits(cb1->GetParNumber("alphaR"),0.01,10);
      cb1->SetParameters((gaussEdgeL-gaussMean)/gaussSigma,5,gaussMean,gaussSigma,gaussNorm,(gaussEdgeR-gaussMean)/gaussSigma,5);
      cb1->SetLineColor(kRed);

    } else {

      if (hNameID != "EoverP") cb1 = new TF1("cb1",&myCrystalBallLeftTail,funcRangeLeft,funcRangeRight,5);  // last parameter is the number of free parameters
      else cb1 = new TF1("cb1",&myCrystalBallRightTail,funcRangeLeft,funcRangeRight,5);
      cb1->SetParNames("alpha","n","mu","sigma","N");  
      cb1->SetParLimits(cb1->GetParNumber("n"),0.1,15); 
      if (hNameID != "EoverP") {
	cb1->SetParLimits(cb1->GetParNumber("alpha"),-10.0,-0.01);
	cb1->SetParameters((gaussEdgeL-gaussMean)/gaussSigma,5,gaussMean,gaussSigma,gaussNorm);
      } else {
	cb1->SetParLimits(cb1->GetParNumber("alpha"),0.01,10);
	cb1->SetParameters((gaussEdgeR-gaussMean)/gaussSigma,5,gaussMean,gaussSigma,gaussNorm);
      }
      // with the following (or some of them) it looks like the fit doesn't work well, tipically the value from fit is out of the range
      // cb1->SetParLimits(cb1->GetParNumber("mu"),gaussMean-3.0*gaussSigma,gaussMean+3.0*gaussSigma);
      // cb1->SetParLimits(cb1->GetParNumber("sigma"),0.1*gaussSigma,10.0*gaussSigma);
      // cb1->SetParLimits(cb1->GetParNumber("N"),0.1*gaussNorm,10.0*gaussNorm);
    }

    // build TLegend
    TLegend *leg = new TLegend(0.11,0.8,0.4,0.89);
    TLegend *legFitFunction = new TLegend(0.11,0.5,0.35,0.8);

    // do the fit
    // before fitting, clone hist . Will use the clone to do the fit with the template. This is needed to plot 2 statistic boxes for crystal Ball and template fit
    TH1F* histClone = (TH1F*) hist->Clone("histClone");
    TFitResultPtr frp1 = hist->Fit(cb1,"E L I S Q B R","HE",funcRangeLeft,funcRangeRight);
    leg->AddEntry(hist,"distribution","l");
    legFitFunction->SetHeader("fit functions:");
    legFitFunction->AddEntry(cb1,"Crystal Ball","l");
    TFitResultPtr frp_template;
    TF1* f_EoP_template = NULL;

    if (sampleName == "DATA" && hNameID == "EoverP") {

      // ----------------------------------------------------------------
      // fitting with template (only for data and for E/P
      
      // ------------------------------------------------------------
      // get binned MC histograms to use them as template
	
      //    getEoP_templateMC(hEoP_template, dirName);
    
      // get E/P template from file

      // before getting the template, delete the same named histogram from data file, otherwise I am getting the same data histogram despite the fact that
      // I'm using the pointer to MC file
      delete htmp;
      fMC->cd();
      TH1F* htmpMC = NULL;
      htmpMC = (TH1F*)fMC->Get(Form("hEoverP_energyBin%1.0fTo%1.0f",energybinEdges[i],energybinEdges[i+1]));  
      if (!htmpMC) {
	cout << "Error: histogram " << Form("hEoverP_energyBin%1.0fTo%1.0f",energybinEdges[i],energybinEdges[i+1]) 
	     << " not found in file ' " << fileNameMC << "'. End of programme." << endl;
	exit(EXIT_FAILURE);
      } else {
	//cout << "Energy ["<< energybinEdges[i] <<"," <<energybinEdges[i+1]<< "] \t max htmpMC: " << htmpMC->GetBinContent(htmpMC->GetMaximumBin()) << endl;
      }

      hEoP_template = (TH1F*) htmpMC->Clone(Form("hEoverP_energyBin%1.0fTo%1.0f_MC",energybinEdges[i],energybinEdges[i+1]));
      if (!hEoP_template) {
	cout << "Error: hEoP_template is NULL! End of programme." << endl;
	exit(EXIT_FAILURE);
      }	
      // rebinning if needed
      if (energybinEdges[i] > 449.9) hEoP_template->Rebin(4); 
      else if (energybinEdges[i] > 274.9) hEoP_template->Rebin(3);      
      else if (energybinEdges[i] > 174.9) hEoP_template->Rebin(2);      

      // end of access of MC file to get templates
      // -----------------------------------------------------------

      // class defining the template function
      histoFunc* templateHistoFunc = new histoFunc(hEoP_template); 
      //cout << "templateHistoFunc->GetIntegral() : " << templateHistoFunc->GetIntegral() << endl;
      Double_t templateFitRangeLow = 0.7;  // 0.7
      Double_t templateFitRangeUp = 1.8;   // 2.0
      if (energybinEdges[i] > 274.9) {
	templateFitRangeLow = 0.75;
	templateFitRangeUp = 1.8; 
      }
      
      f_EoP_template = new TF1("f_EoP_template", templateHistoFunc, templateFitRangeLow, templateFitRangeUp, 3, "histoFunc");
      f_EoP_template -> SetParName(0,"Norm"); 
      f_EoP_template -> SetParName(1,"Scale factor");
      f_EoP_template -> SetParName(2,"Shift of x"); 
      //      f_EoP_template -> SetLineWidth(1); 
      f_EoP_template -> SetNpx(5000);  // granularity of template fit (template shuld be as smooth as possible, so either change this or rebin the MC histogram)
      //Double_t xNorm = histClone->GetEntries()/hEoP_template->GetEntries() * histClone->GetBinWidth(1)/hEoP_template->GetBinWidth(1); // GetEntries includes under/overflow

      //// PAR 0 /////
      Double_t xNorm = 0.0;
      // xnorm = (histClone->Integral(getLowBinGivenRange(histClone,templateFitRangeLow),getUpBinGivenRange(histClone,templateFitRangeUp)) /
      // 	       hEoP_template->Integral(getLowBinGivenRange(hEoP_template,templateFitRangeLow),getUpBinGivenRange(hEoP_template,templateFitRangeUp))) * 
      // 	( histClone->GetBinWidth(1) / hEoP_template->GetBinWidth(1) ); 
      xNorm = (histClone->Integral(getLowBinGivenRange(histClone,templateFitRangeLow),getUpBinGivenRange(histClone,templateFitRangeUp),"width") /
	       hEoP_template->Integral(getLowBinGivenRange(hEoP_template,templateFitRangeLow),getUpBinGivenRange(hEoP_template,templateFitRangeUp),"width"));
      // f_EoP_template -> SetParameter(0, xNorm);
      // f_EoP_template -> SetParLimits(0, 0.98 * xNorm, 1.02 * xNorm);
      f_EoP_template -> FixParameter(0, xNorm);

      //// PAR 1 /////      
      f_EoP_template -> SetParameter(1, 0.9);  // maybe this should be the inverse of the width ratio
      //f_EoP_template -> SetParLimits(1,0.1,10.0); 
      f_EoP_template -> SetParLimits(1,0.8,1.1); 
      
      //// PAR 2 /////
      // I think the shift should be negative: x --> x' = a (x - dx)
      // This means the data are shifter toward left with respect to MC. Also, a is slightly less than 1, so that distribution is broader in data than in MC
      // Suppose f(x) for MC is a gaussian centered in 0 and sigma = 1. Then, g(x') = g( a(x-dx) ) is a gaussian peaked in dx and with sigma = 1/a
      
      // Double_t MC_data_peakDeltaX = hEoP_template->GetBinCenter(hEoP_template->GetMaximumBin()) - frp1->Parameter(2);
      // cout << "MC_data_peakDeltaX = " << MC_data_peakDeltaX << endl;
      //f_EoP_template -> SetParameter(2, MC_data_peakDeltaX);  // will set this to Dx between peak in data and MC ( MC - data )
      f_EoP_template -> SetParameter(2, - 0.05);  // will set this to Dx between peak in data and MC ( MC - data )
      //f_EoP_template -> SetParLimits(2, 0.1 * MC_data_peakDeltaX, 10.0 * MC_data_peakDeltaX);  // will set this to Dx between peak in data and MC ( MC - data )
      f_EoP_template -> SetParLimits(2, -0.5, 0.1);  
      //Double_t MC_data_peakDeltaX = 0.2;
      //f_EoP_template -> FixParameter(2, MC_data_peakDeltaX);  // will set this to Dx between peak in data and MC ( MC - data )
      f_EoP_template -> SetLineColor(kGreen+2);
      // -----------------------------------------------------------------

      frp_template = histClone->Fit("f_EoP_template","E L I S Q B R","HE SAMES", templateFitRangeLow,templateFitRangeUp);
      legFitFunction->AddEntry(f_EoP_template,"MC template","l");

    }

    leg->Draw();
    leg->SetMargin(0.3); 
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    legFitFunction->Draw();
    legFitFunction->SetMargin(0.3); 
    legFitFunction->SetBorderSize(0);
    legFitFunction->SetFillStyle(0);

    c->Update();
    // box for fit with Crystal Ball
    TPaveStats *stat = (TPaveStats*)(hist->FindObject("stats"));
    if(stat) {
      // stat->SetTextColor(kBlue);
      // stat1->SetTextColor(kGreen);
      float width = stat->GetX2NDC() - stat->GetX1NDC();
      // make stat box bigger
      stat->SetX1NDC(stat->GetX1NDC() - 0.5 * width);     
      if (FIT_2SIDE_CB) stat->SetY1NDC(stat->GetY2NDC() - 1.5 * (stat->GetY2NDC() - stat->GetY1NDC()));
      stat->Draw();
    }
    // Box for fit with template (data only)
    TPaveStats *stat1 = NULL;
    if (sampleName == "DATA") {
      histClone->Draw("HE SAME");
      c->Update();
      stat1 = (TPaveStats*)(histClone->FindObject("stats"));
      if (stat1) {
	float height = stat1->GetY2NDC() - stat1->GetY1NDC();
	stat1->SetY1NDC(stat->GetY1NDC() - height);
	stat1->SetY2NDC(stat->GetY1NDC() );
	stat1->SetX1NDC(stat->GetX1NDC());
	stat1->Draw();
      }
    }

    if (hPeak != NULL) {
      hPeak->SetBinContent(i+1, frp1->Parameter(2)); // 2 is mu (starts with alpha, which is parameter number 0)  
      hPeak->SetBinError(i+1, frp1->ParError(2)); // 2 is mu (starts with alpha, which is parameter number 0)  
    }
    if (hPeakMeanInRange != NULL) {
      // Clone the histogram and set its range so that TH1F::GetMean() computes the mean in that range
      TH1F* hCloneForMeanInRange = (TH1F*) hist->Clone();
      TAxis* axisForMeanInRange = (TAxis*) hCloneForMeanInRange->GetXaxis();
      axisForMeanInRange->SetRange(axisForMeanInRange->FindFixBin(lowRangeValue),axisForMeanInRange->FindFixBin(upRangeValue));
      hPeakMeanInRange->SetBinContent(i+1, hCloneForMeanInRange->GetMean());
      hPeakMeanInRange->SetBinError(i+1, hCloneForMeanInRange->GetMeanError());
    }
    if (hScaleFactor_MCdata  != NULL) {
      hScaleFactor_MCdata ->SetBinContent(i+1, frp_template->Parameter(1)); // 1 is scale factor 
      hScaleFactor_MCdata ->SetBinError(i+1, frp_template->ParError(1)); // 1 is scale factor
    }
    if (hPeakShift_MCdata  != NULL) {
      hPeakShift_MCdata ->SetBinContent(i+1, frp_template->Parameter(2)); // 2 is deltaX 
      hPeakShift_MCdata ->SetBinError(i+1, frp_template->ParError(2)); // 2 is deltaX 
    }
    if (hSigma != NULL) {
      hSigma->SetBinContent(i+1, frp1->Parameter(3)); // 2 is mu (starts with alpha, which is parameter number 0)  
      hSigma->SetBinError(i+1, frp1->ParError(3)); // 2 is mu (starts with alpha, which is parameter number 0)  
    }

    string whichEnergy = "";
    if (USE_RAWE) whichEnergy = "raw E";
    else whichEnergy = "corrected E";
    if ((sampleName != "DATA") && (hNameID != "EoverP")) whichEnergy = "true E";

    if (USE_E) {      
      hist->SetTitle(Form("%1.0f < %s [GeV] < %1.0f",energybinEdges[i],whichEnergy.c_str(),energybinEdges[i+1]));
      c->SaveAs(Form("%s%sdistribution_E%1.0fTo%1.0f_%s.pdf",dirName.c_str(),hNameID.c_str(),energybinEdges[i],energybinEdges[i+1],sampleName.c_str()));
      c->SaveAs(Form("%s%sdistribution_E%1.0fTo%1.0f_%s.png",dirName.c_str(),hNameID.c_str(),energybinEdges[i],energybinEdges[i+1],sampleName.c_str()));
    } else {
      hist->SetTitle(Form("%1.0f < %s_{T}[GeV] < %1.0f",energybinEdges[i],whichEnergy.c_str(),energybinEdges[i+1]));
      c->SaveAs(Form("%s%sdistribution_ET%1.0fTo%1.0f_%s.pdf",dirName.c_str(),hNameID.c_str(),energybinEdges[i],energybinEdges[i+1],sampleName.c_str()));
      c->SaveAs(Form("%s%sdistribution_ET%1.0fTo%1.0f_%s.png",dirName.c_str(),hNameID.c_str(),energybinEdges[i],energybinEdges[i+1],sampleName.c_str()));
    }

  }

  //cout << "CHECK 1" << endl;
  delete c;
  //cout << "CHECK 2" << endl;
  //delete htmp;
  //cout << "CHECK 3" << endl;
  //delete hist;  // if I delete hist I get segmentation fault (it didn't happen before trying to fit with template)
  //cout << "CHECK 4" << endl;
  
  cout << endl;

}

//==========================================================================================


void drawPlotDataMC(TH1F* hdata, TH1F* hmc, const string& MCSampleName, const string& xAxisName, const string& yAxisName, const string& canvasName, const string &dirName, const string &canvasTitle = "", const string& hNameToTriggerOptions = "", const Int_t setLogy = 0) {

  TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 

  string hOriginalTitle = "";

  TCanvas *cfit = new TCanvas("cfit","",700,700);

  // now here we go with the canvas                                                                                                                                    
  TPad *subpad_1 = NULL;  // will use it to access specific subpad in canvas                                                                                           
  TPad *subpad_2 = NULL;
  TLegend *leg = new TLegend(0.5,0.7,0.79,0.90);

  subpad_1 = new TPad("pad_1","",0.0,0.28,1.0,1.0);
  subpad_2 = new TPad("pad_2","",0.0,0.0,1.0,0.32);
  subpad_2->SetGridy();
  subpad_2->SetBottomMargin(0.3);
  subpad_1->Draw();
  subpad_2->Draw();
  subpad_1->cd();

  Double_t maximumYaxisValue = -100000.0;
  Double_t minimumYaxisValue = 100000.0;

  // Double_t value = hdata->GetBinContent(hdata->GetMaximumBin()) + hdata->GetBinError(hdata->GetMaximumBin()); 
  // if ( value > maximumYaxisValue) maximumYaxisValue = value; 
  // value =  hmc->GetBinContent(hmc->GetMaximumBin()) + hmc->GetBinError(hmc->GetMaximumBin());
  // if ( value > maximumYaxisValue) maximumYaxisValue = value;
  // value = hdata->GetBinContent(hdata->GetMinimumBin()) - hdata->GetBinError(hdata->GetMinimumBin());
  // if ( value < minimumYaxisValue) minimumYaxisValue = value;
  // value = hmc->GetBinContent(hmc->GetMinimumBin()) - hmc->GetBinError(hmc->GetMinimumBin());
  // if ( value < minimumYaxisValue) minimumYaxisValue = value;

  maximumYaxisValue = TMath::Max(hdata->GetBinContent(hdata->GetMaximumBin()),hmc->GetBinContent(hmc->GetMaximumBin()));
  minimumYaxisValue = TMath::Min(hdata->GetBinContent(hdata->GetMinimumBin()),hmc->GetBinContent(hmc->GetMinimumBin()));
  
  Double_t diff = maximumYaxisValue - minimumYaxisValue;
  minimumYaxisValue -= diff * 0.1;
  maximumYaxisValue += diff * 0.1;

  if (setLogy != 0) {
    if (minimumYaxisValue < 0.0000001) minimumYaxisValue = 0.00001; // if compatible with 0 or negative, set it to value close to but different from zero
    subpad_1->SetLogy(); 
  } 
  hdata->SetStats(0);
  hdata->SetLineColor(kRed);
  hdata->Draw("HE");
  hOriginalTitle = string(hdata->GetTitle());
  hdata->SetTitle(canvasTitle.c_str());
  hdata->GetXaxis()->SetLabelSize(0.45);
  hdata->GetYaxis()->SetTitle(yAxisName.c_str());
  hdata->GetYaxis()->SetTitleSize(0.06);
  hdata->GetYaxis()->SetTitleOffset(0.8);
  hdata->SetMinimum(minimumYaxisValue);
  hdata->SetMaximum(maximumYaxisValue);
  // string hdataName(hdata->GetName());
  // if (hdataName == hNameToTriggerOptions) {
  if (hNameToTriggerOptions == "hPeakEoverPdata_templateFit") {
    hdata->SetMinimum(0.93);
    hdata->SetMaximum(1.0);
  } else if (hNameToTriggerOptions == "distributions") {
    hdata->SetMinimum(0.0001);
    hdata->SetMaximum(1.0);
  }

  if (SET_SCALE_ON_Y) {
    if (yAxisName == "peak(E/P)") {
      hdata->SetMinimum(0.94);
      hdata->SetMinimum(0.965);
      hdata->SetMaximum(1.005);
    } else if (yAxisName == "#sigma(E/P)") {
      hdata->SetMinimum(0.02);
      hdata->SetMaximum(0.15);
    }
  }
  hmc->SetLineColor(kBlue);
  hmc->Draw("HE SAME");

  string texMCSampleName = "";
  getTexMCSampleName(MCSampleName, texMCSampleName);

  leg->AddEntry(hdata,"data","lf");
  leg->AddEntry(hmc,Form("%s",texMCSampleName.c_str()),"lf");
  leg->Draw();
  leg->SetMargin(0.3);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);  // transparent legend

  //ratio plot

  subpad_2->cd();
  TH1F * ratioplot = NULL; // will use it for the ratio plots                                                                                                   
  ratioplot = new TH1F(*hdata);  // to have ratioplot with the same x axis range as hdata (or hmc), I copy it from hdata created above, then I substitute bin content with hdata/hmc                                                                                                                                                   
  ratioplot->Divide(hdata,hmc);
  ratioplot->SetStats(0);
  ratioplot->SetTitle("");
  ratioplot->GetXaxis()->SetLabelSize(0.10);
  ratioplot->GetXaxis()->SetTitle(xAxisName.c_str());
  ratioplot->GetXaxis()->SetTitleSize(0.14);
  ratioplot->GetXaxis()->SetTitleOffset(0.8);
  ratioplot->GetYaxis()->SetLabelSize(0.10);
  ratioplot->GetYaxis()->SetTitle("data / MC");
  ratioplot->GetYaxis()->SetTitleSize(0.12);
  ratioplot->GetYaxis()->SetTitleOffset(0.45);
  ratioplot->GetYaxis()->CenterTitle();
  ratioplot->GetYaxis()->SetNdivisions(011);
  ratioplot->SetMarkerStyle(8);  //medium dot  
  TH1F *ratioplotCopy = (TH1F*) ratioplot->DrawCopy("E");
  //ratioplot->Draw("E");   
  if (hNameToTriggerOptions == "hPeakEoverPdata_templateFit") {
    ratioplotCopy->SetMinimum(0.93);
    ratioplotCopy->SetMaximum(1.05);
  } else if (hNameToTriggerOptions == "distributions") {
    ratioplotCopy->SetMinimum(0.5);
    ratioplotCopy->SetMaximum(1.5);
  }
  if (SET_SCALE_ON_Y) {
    if (yAxisName == "peak(E/P)") {
      ratioplotCopy->SetMinimum(0.97);
      ratioplotCopy->SetMinimum(0.985);
      ratioplotCopy->SetMaximum(1.015);
    } else if (yAxisName == "#sigma(E/P)") {
      // ratioplotCopy->SetMinimum(0.9);
      // ratioplotCopy->SetMaximum(0.14);
    }
  }

  maximumYaxisValue = ratioplotCopy->GetBinContent(ratioplotCopy->GetMaximumBin());
  minimumYaxisValue = ratioplotCopy->GetBinContent(ratioplotCopy->GetMinimumBin());
  diff = maximumYaxisValue - minimumYaxisValue;
  minimumYaxisValue -= diff * 0.1;
  maximumYaxisValue += diff * 0.1;
  if (hNameToTriggerOptions == "distributions") {
    minimumYaxisValue = 0.0;
    maximumYaxisValue = 2.0;
  }
  ratioplotCopy->SetMaximum(maximumYaxisValue);
  ratioplotCopy->SetMinimum(minimumYaxisValue);
  
  if( canvasName== "modeEoverPfromFit" ) {

    // save ratioplot in root file and also produce a TGraphError associating each point to the mean energy in the bin. This information is stored in the files with all histgrams
    string filename = dirName + "EoverP_" + MCSampleName + ".root";	
    TFile *histfile = new TFile((filename).c_str(),"READ");
    if (!histfile || !histfile->IsOpen()) {
      cout << "Error: file \"" << filename << "\" was not opened." << endl;
      exit(EXIT_FAILURE);
    }

    TH1F* htmp = NULL;
    htmp = (TH1F*) histfile->Get("hMeanEnergyInEnergyBin");  
    if (!htmp) {
      cout << "Error: histogram not found in file ' " << filename << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
    TH1F* hMeanEnergy = (TH1F*) htmp->Clone(); 

    string tmpString = "";
    if (USE_RAWE) tmpString = "Raw";
    else tmpString = "RegrCorr";
    TH1F * histo = (TH1F*) ratioplot->Clone( ("EoverP_vs" + tmpString + "Energy_dataMCRatio").c_str());

    TGraphErrors *graph = new TGraphErrors(); 
    for (Int_t i = 0; i < hMeanEnergy->GetNbinsX(); i++) {
      graph->SetPoint(i,hMeanEnergy->GetBinContent(i+1),histo->GetBinContent(i+1));
      graph->SetPointError(i,hMeanEnergy->GetBinError(i+1),histo->GetBinError(i+1));
    }

    TFile *f = new TFile( (dirName + "EoverP_vs" + tmpString + "Energy_dataMCRatio.root").c_str(),"RECREATE");
    f->cd();
    histo->Write();
    string graphName(histo->GetName());
    graphName = "graph" + graphName; 
    graph->Write(graphName.c_str());
    f->Close();
    delete f;

  }


  cfit->SaveAs( (dirName + canvasName + ".pdf").c_str() );
  cfit->SaveAs( (dirName + canvasName + ".png").c_str() );

  delete cfit;

  hdata->SetTitle(hOriginalTitle.c_str());

}

//==========================================================================================


void drawPlotOnlyMC(vector<TH1F*> &hmcVector, const vector<string> &legEntryName, const string& MCSampleName, const string& xAxisName, const string& yAxisName, const string& canvasName, const string &dirName, const string &canvasTitle = "", const Double_t miny = -1.0, const Double_t maxy = -1.0) {

  TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2()

  // to draw the canvas title in the frame, we actually use TH1::GetTitle, setting it to canvasTitle argument passed to this function.
  // after saving the canvas, we set the histogram title to its original value
  string hOriginalTitle = "";

  TCanvas *cfitMC = new TCanvas("cfitMC","");

  // now here we go with the canvas                                                                                                                                    
  // TPad *subpad_1 = NULL;  // will use it to access specific subpad in canvas                                                                                        
  // TPad *subpad_2 = NULL;
  TLegend *leg = new TLegend(0.5,0.6,0.79,0.90);

  // subpad_1 = new TPad("pad_1","",0.0,0.28,1.0,1.0);
  // subpad_2 = new TPad("pad_2","",0.0,0.0,1.0,0.32);
  // subpad_2->SetGridy();
  // subpad_2->SetBottomMargin(0.3);
  // subpad_1->Draw();
  // subpad_2->Draw();
  // subpad_1->cd();

  vector<Int_t> histColor;
  setHistColor(histColor,(Int_t)hmcVector.size());
  
  Double_t maximumYaxisValue = -100000.0;
  Double_t minimumYaxisValue = 100000.0;
  
  for (UInt_t i = 0; i < hmcVector.size(); i++) {

    Double_t value = hmcVector[i]->GetBinContent(hmcVector[i]->GetMaximumBin()); 
    if ( value > maximumYaxisValue) maximumYaxisValue = value; 
    value = hmcVector[i]->GetBinContent(hmcVector[i]->GetMinimumBin());
    if ( value < minimumYaxisValue) minimumYaxisValue = value;

  } 

  // maximumYaxisValue *= 1.02; // slightly increase the y scale for maximum
  // minimumYaxisValue *= 0.98; // slightly decrease the y scale for minimum

  Double_t diff = maximumYaxisValue - minimumYaxisValue;
  minimumYaxisValue -= diff * 0.1;
  maximumYaxisValue += diff * 0.2;

  for (UInt_t i = 0; i < hmcVector.size(); i++) {

    hmcVector[i]->SetStats(0);
    hmcVector[i]->SetLineColor(histColor[i]);
    if ( i == 0 ) {
      hmcVector[i]->SetMaximum(maximumYaxisValue);
      hmcVector[i]->SetMinimum(minimumYaxisValue);
      if (hmcVector[i]->GetMinimum() < 0.02) hmcVector[i]->SetMinimum(0.0);
      if (miny > 0.0) hmcVector[i]->SetMinimum(miny);
      if (maxy > 0.0) hmcVector[i]->SetMaximum(maxy);
      if (FIT_2SIDE_CB) {
	if (canvasName == "modeMCstudy") {
	  hmcVector[i]->SetMinimum(0.96);
	  hmcVector[i]->SetMaximum(maximumYaxisValue);
	} else if (canvasName == "sigmaMCstudy") {
	  hmcVector[i]->SetMinimum(0.0);
	  hmcVector[i]->SetMaximum(0.15);
	}
      }
      // if (canvasName == "modeMCstudy") hmcVector[i]->SetMinimum(0.9);
      // else if (canvasName == "sigmaMCstudy") hmcVector[i]->SetMaximum(0.2);
      hmcVector[i]->Draw("HE");
      hOriginalTitle = string(hmcVector[i]->GetTitle());
      hmcVector[i]->SetTitle(canvasTitle.c_str());
      hmcVector[i]->GetXaxis()->SetTitle(xAxisName.c_str());
      hmcVector[i]->GetXaxis()->SetLabelSize(0.05);
      hmcVector[i]->GetXaxis()->SetTitleSize(0.06);
      hmcVector[i]->GetXaxis()->SetTitleOffset(0.8);
      hmcVector[i]->GetYaxis()->SetTitle(yAxisName.c_str());
      hmcVector[i]->GetYaxis()->SetTitleSize(0.06);
      hmcVector[i]->GetYaxis()->SetTitleOffset(0.8);
      if (SET_SCALE_ON_Y) {
	if (yAxisName == "peak position") hmcVector[i]->GetYaxis()->SetRangeUser(0.98,1.04);
	else if (yAxisName == "#sigma of distribution") hmcVector[i]->GetYaxis()->SetRangeUser(0.0,0.11);
      }
    } else {
      hmcVector[i]->Draw("HE SAME");
    }

  }

  string texMCSampleName = "";
  getTexMCSampleName(MCSampleName, texMCSampleName);

  leg->AddEntry((TObject*)0,texMCSampleName.c_str(),"");
  for (UInt_t i = 0; i < hmcVector.size(); i++) {
    leg->AddEntry(hmcVector[i],Form("%s",legEntryName[i].c_str()),"lf");
  }
  leg->Draw();
  leg->SetMargin(0.3);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);  // transparent legend

  cfitMC->SaveAs( (dirName + canvasName + ".pdf").c_str() );
  cfitMC->SaveAs( (dirName + canvasName + ".png").c_str() );

  delete cfitMC;

  hmcVector[0]->SetTitle(hOriginalTitle.c_str());  // reset the title of the histogram used to draw the canvas title (it was set in this function)


}

//===========================================================================

void  plotHistogramsFromFile(const string &dataSampleName, const string &MCSampleName, const string& dirName) { 

  /*
  void drawPlotDataMC(TH1F* hdata, 
		      TH1F* hmc, 
		      const string& MCSampleName, 
		      const string& xAxisName, 
		      const string& yAxisName, 
		      const string& canvasName, 
		      const string &dirName, 
		      const string &canvasTitle = "", 
		      const string& hNameToTriggerOptions = "") {
  */

  TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2()  

  string fileNameData = dirName + "EoverP_" + dataSampleName + ".root";
  string fileNameMC = dirName + "EoverP_" + MCSampleName + ".root";

  TH1F* hvar = NULL;  // to get histogram from file
  TH1F* hvar2 = NULL;  // to get histogram from file

  vector<string> hNameList;
  hNameList.push_back("hSigmaIetaIeta_lowE");
  hNameList.push_back("hSigmaIetaIeta_hotE");
  hNameList.push_back("hR9_lowE");
  hNameList.push_back("hR9_hotE");
  hNameList.push_back("hPt_lowE");
  hNameList.push_back("hPt_hotE");
  hNameList.push_back("hElePt");
  hNameList.push_back("hEleEcorr");
  hNameList.push_back("hEleEraw");
  vector<string> xAxisName;
  xAxisName.push_back("electron #sigma_{i#eta i#eta}");
  xAxisName.push_back("electron #sigma_{i#eta i#eta}");
  xAxisName.push_back("electron R9");
  xAxisName.push_back("electron R9");
  xAxisName.push_back("electron p_{T} [GeV]");
  xAxisName.push_back("electron p_{T} [GeV]");
  xAxisName.push_back("electron p_{T} [GeV]");
  xAxisName.push_back("electron corrected energy [GeV]");
  xAxisName.push_back("electron raw energy [GeV]");
  vector<string> canvasName;
  canvasName.push_back("sigmaIetaIeta_lowE");
  canvasName.push_back("sigmaIetaIeta_hotE");
  canvasName.push_back("R9_lowE");
  canvasName.push_back("R9_hotE");
  canvasName.push_back("pT_lowE");
  canvasName.push_back("pT_hotE");
  canvasName.push_back("pT");
  canvasName.push_back("corrE");
  canvasName.push_back("rawE");
  vector<string> canvasTitle;

  string whichEnergy = "";
  if (USE_RAWE) whichEnergy = "raw";
  else whichEnergy = "corrected";

  if (USE_E) {
    canvasTitle.push_back("100 < " +  whichEnergy + " E [GeV] < 225");
    canvasTitle.push_back("350 < " +  whichEnergy + " E [GeV] < 450");
    canvasTitle.push_back("100 < " +  whichEnergy + " E [GeV] < 225");
    canvasTitle.push_back("350 < " +  whichEnergy + " E [GeV] < 450");
    canvasTitle.push_back("100 < " +  whichEnergy + " E [GeV] < 225");
    canvasTitle.push_back("350 < " +  whichEnergy + " E [GeV] < 450");
    canvasTitle.push_back("");
    canvasTitle.push_back("");
    canvasTitle.push_back("");
  } else {
    canvasTitle.push_back("100 < " +  whichEnergy + " E_{T} [GeV] < 225");
    canvasTitle.push_back("325 < " +  whichEnergy + " E_{T} [GeV] < 425");
    canvasTitle.push_back("100 < " +  whichEnergy + " E_{T} [GeV] < 225");
    canvasTitle.push_back("325 < " +  whichEnergy + " E_{T} [GeV] < 425");
    canvasTitle.push_back("100 < " +  whichEnergy + " E_{T} [GeV] < 225");
    canvasTitle.push_back("325 < " +  whichEnergy + " E_{T} [GeV] < 425");
    canvasTitle.push_back("");
    canvasTitle.push_back("");
    canvasTitle.push_back("");
  }    

  TH1F* hData = NULL;
  TH1F* hMC = NULL;

  cout << "check here" << endl;

  TFile* fData = TFile::Open(fileNameData.c_str(),"READ");
  if (!fData || !fData->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<fileNameData<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }
  
  TFile* fMC = TFile::Open(fileNameMC.c_str(),"READ");
  if (!fMC || !fMC->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<fileNameMC<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }
  
  for (UInt_t i = 0; i < hNameList.size(); i++) {
    
    fData->cd();
    // reading data file
    hvar = (TH1F*)fData->Get(hNameList[i].c_str());
    if (!hvar) {
      cout << "Error: histogram not found in file ' " << fileNameData << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
    hData = (TH1F*) hvar->Clone();

    fMC->cd();
    // reading MC file
    hvar2 = (TH1F*)fMC->Get(hNameList[i].c_str());
    if (!hvar2) {
      cout << "Error: histogram not found in file ' " << fileNameMC << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
    hMC = (TH1F*) hvar2->Clone();

    /////////////////////////
    // Now plotting
    ////////////////////////
    // to compare, normalize MC to same area
    // since no weights are applied, I can use GetEntries() to get integral (this includes underflow and overflow events)
    hData->Sumw2();
    hData->Scale(1./hData->Integral());
    hMC->Sumw2();
    hMC->Scale(1./hMC->Integral());
    drawPlotDataMC(hData, hMC, MCSampleName, xAxisName[i], "a.u.",canvasName[i],dirName,canvasTitle[i],"distributions",1);
    //////////////////////////

  }

}


//=================================================================================

void plotFromFit(const string &dataSampleName, const string &MCSampleName, const vector<Float_t> &energybinEdges, const string & dirName) {

  TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 
  
  string xAxisName = "";

  string whichEnergy = "";
  if (USE_RAWE) whichEnergy = "raw";
  else whichEnergy = "corrected";

  if (USE_E) xAxisName = whichEnergy + " E [GeV]";
  else xAxisName = whichEnergy + " E_{T} [GeV]";

  Int_t nEnergyBins = energybinEdges.size() -1;

  TH1F* hPeakEoverPdata = new TH1F("hPeakEoverPdata","",nEnergyBins,energybinEdges.data());
  TH1F* hSigmaEoverPdata = new TH1F("hSigmaEoverPdata","",nEnergyBins,energybinEdges.data());
  TH1F* hPeakEoverPmc = new TH1F("hPeakEoverPmc","",nEnergyBins,energybinEdges.data());
  TH1F* hSigmaEoverPmc = new TH1F("hSigmaEoverPmc","",nEnergyBins,energybinEdges.data());

  TH1F* hPeakShift_MCdata = new TH1F("hPeakShift_MCdata","",nEnergyBins,energybinEdges.data()); 
  TH1F* hScaleFactor_MCdata = new TH1F("hScaleFactor_MCdata","",nEnergyBins,energybinEdges.data()); 

  // valuse used to define the range when getting peak position from TH1::GetMean()
  double_t lowRangeValueForMean = 0.8; 
  Double_t upRangeValueForMean = 1.2; 
  string titleForPlotsWithMeanInRange = string(Form("peak from TH1F::GetMean() in [%1.1f, %1.1f]",lowRangeValueForMean,upRangeValueForMean));

  TH1F* hPeakEoverPdata_meanInRange = new TH1F("hPeakEoverPdata_meanInRange","",nEnergyBins,energybinEdges.data());
  TH1F* hPeakEoverPmc_meanInRange = new TH1F("hPeakEoverPmc_meanInRange","",nEnergyBins,energybinEdges.data());

  plotDistribution(dataSampleName, energybinEdges, hPeakEoverPdata, hSigmaEoverPdata, "EoverP",dirName, 
		   hPeakEoverPdata_meanInRange,lowRangeValueForMean, upRangeValueForMean, hPeakShift_MCdata, hScaleFactor_MCdata);
  plotDistribution(MCSampleName, energybinEdges, hPeakEoverPmc, hSigmaEoverPmc, "EoverP",dirName,
		   hPeakEoverPmc_meanInRange,lowRangeValueForMean, upRangeValueForMean);

  TH1F* hPeakEoverPdata_templateFit = NULL;
  if (hPeakShift_MCdata && hScaleFactor_MCdata) hPeakEoverPdata_templateFit = new TH1F("hPeakEoverPdata_templateFit","",nEnergyBins,energybinEdges.data());
  if (hPeakEoverPdata_templateFit) {
    // the following is because we get peak of data distribution from template by taking MC, dividing by scale factor in each bin and adding the shift
    // this is because the template acts by changing the x in MC, x_MC, to x in data, x_D in the following way:
    // x_MC = (x_D - mu_D)/sigma_D, where mu_D is the shift and 1/sigma_D is the scale factor
    hPeakEoverPdata_templateFit->Sumw2();
    hPeakEoverPdata_templateFit->Divide(hPeakEoverPmc,hScaleFactor_MCdata,1.0,1.0);
    hPeakEoverPdata_templateFit->Add(hPeakShift_MCdata,1.0);  // data = MC/scale_factor + shift (with shift < 0 since data are shifted to the left of MC)
  }

  drawPlotDataMC(hPeakEoverPdata, hPeakEoverPmc, MCSampleName, xAxisName, "peak(E/P) from fit", "modeEoverPfromFit",dirName,"fit with Crystal Ball");
  if (hPeakEoverPdata_templateFit) drawPlotDataMC(hPeakEoverPdata_templateFit, hPeakEoverPmc, MCSampleName, xAxisName, "peak(E/P) from fit", 
						  "modeEoverPfromTemplateFit",dirName,"fit with template (data) and Crystal Ball (MC)","hPeakEoverPdata_templateFit");
  drawPlotDataMC(hSigmaEoverPdata, hSigmaEoverPmc, MCSampleName, xAxisName, "#sigma(E/P) from fit", "sigmaEoverPfromFit",dirName,"fit with Crystal Ball");
  if ( (hPeakEoverPdata_meanInRange != NULL) && (hPeakEoverPmc_meanInRange != NULL) ) {
    drawPlotDataMC(hPeakEoverPdata_meanInRange, hPeakEoverPmc_meanInRange, MCSampleName, xAxisName, "peak(E/P)", "modeEoverP_usingMean",dirName,
		   titleForPlotsWithMeanInRange);
  }

  TH1F* hPeakEcorrOverEtrue = new TH1F("hPeakEcorrOverEtrue","",nEnergyBins,energybinEdges.data());
  TH1F* hSigmaEcorrOverEtrue = new TH1F("hSigmaEcorrOverEtrue","",nEnergyBins,energybinEdges.data());
  TH1F* hPeakErawOverEtrue = new TH1F("hPeakErawOverEtrue","",nEnergyBins,energybinEdges.data());
  TH1F* hSigmaErawOverEtrue = new TH1F("hSigmaErawOverEtrue","",nEnergyBins,energybinEdges.data());
  TH1F* hPeakPtrackOverEtrue = new TH1F("hPeakPtrackOverEtrue","",nEnergyBins,energybinEdges.data());
  TH1F* hSigmaPtrackOverEtrue = new TH1F("hSigmaPtrackOverEtrue","",nEnergyBins,energybinEdges.data());

  TH1F* hPeakEcorrOverEtrue_meanInRange = new TH1F("hPeakEcorrOverEtrue_meanInRange","",nEnergyBins,energybinEdges.data());
  TH1F* hPeakErawOverEtrue_meanInRange = new TH1F("hPeakErawOverEtrue_meanInRange","",nEnergyBins,energybinEdges.data());
  TH1F* hPeakPtrackOverEtrue_meanInRange = new TH1F("hPeakPtrackOverEtrue_meanInRange","",nEnergyBins,energybinEdges.data());

  plotDistribution(MCSampleName, energybinEdges, hPeakEcorrOverEtrue, hSigmaEcorrOverEtrue, "EcorrOverEtrue",dirName, 
		   hPeakEcorrOverEtrue_meanInRange,lowRangeValueForMean, upRangeValueForMean);
  plotDistribution(MCSampleName, energybinEdges, hPeakErawOverEtrue, hSigmaErawOverEtrue, "ErawOverEtrue",dirName, 
		   hPeakErawOverEtrue_meanInRange,lowRangeValueForMean, upRangeValueForMean);
  plotDistribution(MCSampleName, energybinEdges, hPeakPtrackOverEtrue, hSigmaPtrackOverEtrue, "PtrackOverEtrue",dirName, 
		   hPeakPtrackOverEtrue_meanInRange,lowRangeValueForMean, upRangeValueForMean);

  vector<TH1F*> hPeakVectorMC;
  hPeakVectorMC.push_back(hPeakEcorrOverEtrue);
  hPeakVectorMC.push_back(hPeakErawOverEtrue);
  hPeakVectorMC.push_back(hPeakPtrackOverEtrue);

  vector<TH1F*> hSigmaVectorMC;
  hSigmaVectorMC.push_back(hSigmaEcorrOverEtrue);
  hSigmaVectorMC.push_back(hSigmaErawOverEtrue);
  hSigmaVectorMC.push_back(hSigmaPtrackOverEtrue);

  vector<TH1F*> hPeakVectorMC_meanInRange;
  hPeakVectorMC_meanInRange.push_back(hPeakEcorrOverEtrue_meanInRange);
  hPeakVectorMC_meanInRange.push_back(hPeakErawOverEtrue_meanInRange);
  hPeakVectorMC_meanInRange.push_back(hPeakPtrackOverEtrue_meanInRange);

  vector<string> legEntryName;
  legEntryName.push_back("E_{corr}/E_{true}");
  legEntryName.push_back("E_{raw}/E_{true}");
  legEntryName.push_back("P_{track}/E_{true}");

  xAxisName = "true E [GeV]";
  drawPlotOnlyMC(hPeakVectorMC, legEntryName, MCSampleName, xAxisName, "peak position from fit", "modeMCstudy",dirName,"fit with Crystall Ball"); 
  drawPlotOnlyMC(hSigmaVectorMC, legEntryName, MCSampleName, xAxisName, "#sigma of distribution from fit", "sigmaMCstudy",dirName, "fit with Crystall Ball"); 

  drawPlotOnlyMC(hPeakVectorMC_meanInRange, legEntryName, MCSampleName, xAxisName, "peak position", "peakMCstudy_usingMean",dirName,
		 titleForPlotsWithMeanInRange); 

}

//======================================================================


Int_t main(Int_t argc, char* argv[]) {

  Int_t doAll_flag = 1;
  Int_t doLoop_flag = 1;
  Int_t doLoopMC_flag = 1;
  Int_t doLoopData_flag = 1; 
  Int_t doPlot_flag = 1;
  Int_t doLoop4weight_flag = 0;

  string dirName = "";   // will be a path like "plot/<name>/" . Note the ending "/" 

 
  if (argc > 1) {

    for (Int_t i = 1; i < argc; i++) {

      string thisArgument(argv[i]);
      if (thisArgument == "-nl") {
	cout << "Passing option -nl: skip Loop on ntuples" << endl;
	doLoop_flag = 0;  // -nl --> no Loop
	doLoopMC_flag = 0;  // -nl --> no Loop, so it implies no loop on MC
	doLoopData_flag = 0;  // -nl --> no Loop, so it implies no loop on data
	doAll_flag = 0;
      }	else if (thisArgument == "-nlMC") {
	cout << "Passing option -nlMC: skip Loop on MC ntuples (but do loop for data)." << endl; // Useful when new data are available but MC didn't change.
	doLoopMC_flag = 0; 
      } else if (thisArgument == "-nldata") {
	cout << "Passing option -nldata: skip Loop on data ntuples (but do loop for MC)." << endl; // Useful when new data are available but MC didn't change.
	doLoopData_flag = 0; 
      } else if (thisArgument == "-np") {
	cout << "Passing option -np: skip creation of plots" << endl;
	doPlot_flag = 0;  // -np --> no plots
	doAll_flag = 0;
      } else if (thisArgument == "-lw") {
	cout << "Passing option -lw: running with data/MC weight." << endl;
	doLoop4weight_flag = 1;  // 
      } else if (thisArgument == "-dn") {   // -dn --> directory name
	cout << "Passing option -dn: passing name for directories to be created" << endl;
	dirName = string(argv[i+1]);
	//option->Set_dirName(string(argv[i+1]));
	cout << "Saving output in '" << dirName << "'" << endl;
	i++;
      }

    }

  }

  Option* option = new Option();
  option->Set_data2016(1);
  option->Set_skim1lep1jet80X(1);
  option->Set_fit2sideCB(0);
  option->Set_setScaleOnY(1);
  option->Set_useE(1);
  option->Set_dirName(dirName);

  vector<string> sampleName;
  sampleName.push_back("DATA");
  sampleName.push_back("WJetsToLNu");

  vector<Float_t> energybinEdges;
  // this binning looks ok
  energybinEdges.push_back(25.0);
  energybinEdges.push_back(50.0);
  energybinEdges.push_back(75.0);
  energybinEdges.push_back(100.0);
  energybinEdges.push_back(125.0);
  energybinEdges.push_back(175.0);
  energybinEdges.push_back(225.0);
  if (USE_E) {
    energybinEdges.push_back(275.0);
    energybinEdges.push_back(350.0);
    energybinEdges.push_back(450.0);
    energybinEdges.push_back(650.0);
    energybinEdges.push_back(900.0);
  } else {
    energybinEdges.push_back(275.0);
    //energybinEdges.push_back(350.0);
    energybinEdges.push_back(325.0);
    energybinEdges.push_back(425.0);
    energybinEdges.push_back(600.0);
    energybinEdges.push_back(900.0);
  }    

  // energybinEdges.push_back(25.0);
  // energybinEdges.push_back(50.0);
  // energybinEdges.push_back(75.0);
  // energybinEdges.push_back(100.0);
  // energybinEdges.push_back(125.0);
  // energybinEdges.push_back(150.0);
  // energybinEdges.push_back(200.0);
  // energybinEdges.push_back(250.0);
  // energybinEdges.push_back(325.0);
  // energybinEdges.push_back(450.0);
  // energybinEdges.push_back(900.0);

  if(doAll_flag || doLoop_flag) {

    for (UInt_t i = 0; i < sampleName.size(); i++) {
    
      if (!doLoopMC_flag && (sampleName[i] != "DATA")) continue;
      if (!doLoopData_flag && (sampleName[i] == "DATA")) continue;

      // create chain                                                                                                         
    
      TChain* chain = new TChain("tree");
      TChain* chFriend = new TChain("mjvars/t");
      TChain* chsfFriend = new TChain("sf/t");

      buildChainWithFriend(chain, chFriend, sampleName[i], chsfFriend);
    
      if(!chain || !chFriend) {
	cout << "Error: chain not created. End of programme" << endl;
	exit(EXIT_FAILURE);
      }

      EoverP tree(chain);
      tree.Loop(sampleName[i], energybinEdges, dirName, 0);
      if (sampleName[i] != "DATA" && doLoop4weight_flag) {
	cout <<"Creating data/MC weight histogram to reweight MC" << endl;
	system("root -l -q 'makeEcorrWeight.C++(\"DATA\",\"WJetsToLNu\",dirName.c_str())");
	cout <<"Looping again on MC with data/MC weight for MC" << endl;
	tree.Loop(sampleName[i], energybinEdges, dirName);
      }

      delete chain;
      delete chFriend;
      delete chsfFriend;

    }

  }

  //  cout << "check" << endl;

  if(doAll_flag || doPlot_flag) {

    plotFromFit(sampleName[0],sampleName[1], energybinEdges, dirName);
    plotHistogramsFromFile(sampleName[0],sampleName[1], dirName);

  }

  // plotEoverP(sampleName[0],sampleName[1],"hMeanEoverP","< E/P >");
  // plotEoverP(sampleName[0],sampleName[1],"hModeEoverP"," mode of E/P");
  // plotEoverPdistribution(sampleName[0], energybinEdges);
  // plotEoverPdistribution(sampleName[1], energybinEdges);

  return 0;

}



////////////////////////////////////////
//TO BE IMPLEMENTED
//
//Plots of E/P asaf mean energy
//
