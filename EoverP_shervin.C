#define EoverP_shervin_cxx
#include "EoverP_shervin.h"

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
#include <TLegendEntry.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TProfile.h>
#include <TTreeIndex.h>
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

#define FIT_2SIDE_CB 1      // 0 for single tail Crystal Ball for the fit, 1 for double tail
#define SET_SCALE_ON_Y 0    // select default or user defined ranges for y axis
#define USE_E 1 // 0 for ET and 1 for E in the binning
#define USE_RAWE 0 // when 1, use raw SC energy instead of regression corrected ECAL energy
#define DO_TEMPLATE_FIT 0  // when 0 will not do template fit
#define READ_FROM_LOCAL 1 // I skimmed and copied Shervin's ntuples. If the selection is the same as or tighter than the skim, better to read the skimmed version
#define USE_P_AT_VTX 1  // choose p at vtx or mode of p (two variables in Sehrvin's ntuples
#define BIN_P 0 // decide to bin in P or E (or Eraw depending on other defines)
#define ELE_ETA_MAX 1.0

using namespace std;


// must use this variables
// recoFlagsEle[0] > 1   (says if electron is tracker or ECAL driven)
// eleID[0], to be used as -->  (eleID[0] & 0x0008) == 0x0008   (see https://github.com/GiuseppeFasanella/ECALELF/blob/miniAOD/ZNtupleDumper/interface/eleIDMap.h)
// in this way we use a tight ID on the electron

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

//===========================================================================

void buildChain(TChain* chain, string sampleName) {
  
  cout << "Creating chain ..." << endl;
  
  vector<string> subSampleNameVector;

  if (sampleName == "DATA") {

    if (READ_FROM_LOCAL) {
      
      subSampleNameVector.push_back("/afs/cern.ch/work/m/mciprian/EoverP_study_new/CMSSW_8_0_10/src/eleEoverP/TREES_shervin_skim/DATA_shervin.root");
 
    } else {

      subSampleNameVector.push_back("SingleElectron-Run2016B-WSkim-Prompt_v2/273150-275376/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016B-WSkim-Prompt_v2-273150-275376.root");
      subSampleNameVector.push_back("SingleElectron-Run2016C-WSkim-Prompt_v2/275420-276283/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016C-WSkim-Prompt_v2-275420-276283.root");
      subSampleNameVector.push_back("SingleElectron-Run2016D-WSkim-Prompt_v2/276315-276811/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016D-WSkim-Prompt_v2-276315-276811.root");
      subSampleNameVector.push_back("SingleElectron-Run2016E-WSkim-Prompt-v2/276830-277420/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016E-WSkim-Prompt-v2-276830-277420.root");
      subSampleNameVector.push_back("SingleElectron-Run2016F-WSkim-Prompt-v1/277820-278808/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016F-WSkim-Prompt-v1-277820-278808.root");
      subSampleNameVector.push_back("SingleElectron-Run2016G-WSkim-Prompt-v1/278817-279588/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016G-WSkim-Prompt-v1-278817-279588.root");
      subSampleNameVector.push_back("SingleElectron-Run2016G-WSkim-Prompt-v1/279589-279931/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016G-WSkim-Prompt-v1-279589-279931.root");

    }
    
  } else if (sampleName == "DATA_w") {


    if (READ_FROM_LOCAL) {
      
      subSampleNameVector.push_back("/afs/cern.ch/work/m/mciprian/EoverP_study_new/CMSSW_8_0_10/src/eleEoverP/TREES_shervin_skim/DATA_w_shervin.root");
 
    } else {

      subSampleNameVector.push_back("SingleElectron-Run2016B-WSkim-Prompt_v2-weightsReco/273150-275376/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016B-WSkim-Prompt_v2-weightsReco-273150-275376.root");
      subSampleNameVector.push_back("SingleElectron-Run2016C-WSkim-Prompt_v2-weightsReco/275420-276283/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016C-WSkim-Prompt_v2-weightsReco-275420-276283.root");
      subSampleNameVector.push_back("SingleElectron-Run2016D-WSkim-Prompt_v2-weightsReco/276315-276811/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016D-WSkim-Prompt_v2-weightsReco-276315-276811.root");
      subSampleNameVector.push_back("SingleElectron-Run2016E-WSkim-Prompt-v2-weightsReco/276830-277420/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016E-WSkim-Prompt-v2-weightsReco-276830-277420.root");
      subSampleNameVector.push_back("SingleElectron-Run2016F-WSkim-Prompt-v1-weightsReco/277820-278808/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016F-WSkim-Prompt-v1-weightsReco-277820-278808.root");
      subSampleNameVector.push_back("SingleElectron-Run2016G-WSkim-Prompt-v1-weightsReco/278817-279588/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016G-WSkim-Prompt-v1-weightsReco-278817-279588.root");
      subSampleNameVector.push_back("SingleElectron-Run2016G-WSkim-Prompt-v1-weightsReco/279589-279931/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016G-WSkim-Prompt-v1-weightsReco-279589-279931.root");

    }
  
  } else {

    cout << "Error: unknown sampleName " << sampleName <<". End of programme" << endl;
    exit(EXIT_FAILURE);

  }

  string treePath = "";
  treePath = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/ecalelf/ntuples/13TeV/ALCARERECO/Cal_Sep2016_final_v3/";

  for(UInt_t i = 0; i < subSampleNameVector.size(); i++) {
  
    string treeRootFile = treePath + subSampleNameVector[i];
    if (READ_FROM_LOCAL) treeRootFile = subSampleNameVector[i];
    chain->Add(TString(treeRootFile.c_str()));

  }

  if(!chain) {
    cout << "Error: chain not created. End of programme" << endl;
    exit(EXIT_FAILURE);
  }
  cout << chain->GetEntries() << endl;

}

//===============================================================================================

void buildOtherDataChain(TChain* chain, string sampleName) {
  
  // this functions is used in loop to open DATA_w when looping on DATA or viceversa.
  // the reason is that I need to run on the same events (after the selection the number won't be the same anymore)

  cout << "Creating other data chain ..." << endl;
  
  vector<string> subSampleNameVector;

  if (sampleName == "DATA_w") {
    // here reading DATA, because in the loop we are reading DATA_w

    if (READ_FROM_LOCAL) {

      subSampleNameVector.push_back("/afs/cern.ch/work/m/mciprian/EoverP_study_new/CMSSW_8_0_10/src/eleEoverP/TREES_shervin_skim/DATA_shervin.root"); 

    } else {

      subSampleNameVector.push_back("SingleElectron-Run2016B-WSkim-Prompt_v2/273150-275376/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016B-WSkim-Prompt_v2-273150-275376.root");
      subSampleNameVector.push_back("SingleElectron-Run2016C-WSkim-Prompt_v2/275420-276283/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016C-WSkim-Prompt_v2-275420-276283.root");
      subSampleNameVector.push_back("SingleElectron-Run2016D-WSkim-Prompt_v2/276315-276811/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016D-WSkim-Prompt_v2-276315-276811.root");
      subSampleNameVector.push_back("SingleElectron-Run2016E-WSkim-Prompt-v2/276830-277420/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016E-WSkim-Prompt-v2-276830-277420.root");
      subSampleNameVector.push_back("SingleElectron-Run2016F-WSkim-Prompt-v1/277820-278808/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016F-WSkim-Prompt-v1-277820-278808.root");
      subSampleNameVector.push_back("SingleElectron-Run2016G-WSkim-Prompt-v1/278817-279588/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016G-WSkim-Prompt-v1-278817-279588.root");
      subSampleNameVector.push_back("SingleElectron-Run2016G-WSkim-Prompt-v1/279589-279931/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016G-WSkim-Prompt-v1-279589-279931.root");

    }
    
  } else if (sampleName == "DATA") {

    if (READ_FROM_LOCAL) {
      
      subSampleNameVector.push_back("/afs/cern.ch/work/m/mciprian/EoverP_study_new/CMSSW_8_0_10/src/eleEoverP/TREES_shervin_skim/DATA_w_shervin.root");
 
    } else {

      subSampleNameVector.push_back("SingleElectron-Run2016B-WSkim-Prompt_v2-weightsReco/273150-275376/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016B-WSkim-Prompt_v2-weightsReco-273150-275376.root");
      subSampleNameVector.push_back("SingleElectron-Run2016C-WSkim-Prompt_v2-weightsReco/275420-276283/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016C-WSkim-Prompt_v2-weightsReco-275420-276283.root");
      subSampleNameVector.push_back("SingleElectron-Run2016D-WSkim-Prompt_v2-weightsReco/276315-276811/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016D-WSkim-Prompt_v2-weightsReco-276315-276811.root");
      subSampleNameVector.push_back("SingleElectron-Run2016E-WSkim-Prompt-v2-weightsReco/276830-277420/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016E-WSkim-Prompt-v2-weightsReco-276830-277420.root");
      subSampleNameVector.push_back("SingleElectron-Run2016F-WSkim-Prompt-v1-weightsReco/277820-278808/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016F-WSkim-Prompt-v1-weightsReco-277820-278808.root");
      subSampleNameVector.push_back("SingleElectron-Run2016G-WSkim-Prompt-v1-weightsReco/278817-279588/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016G-WSkim-Prompt-v1-weightsReco-278817-279588.root");
      subSampleNameVector.push_back("SingleElectron-Run2016G-WSkim-Prompt-v1-weightsReco/279589-279931/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016G-WSkim-Prompt-v1-weightsReco-279589-279931.root");

    }
  
  } else {

    cout << "Error: unknown sampleName " << sampleName <<". End of programme" << endl;
    exit(EXIT_FAILURE);

  }

  string treePath = "";
  treePath = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/ecalelf/ntuples/13TeV/ALCARERECO/Cal_Sep2016_final_v3/";

  for(UInt_t i = 0; i < subSampleNameVector.size(); i++) {
  
    string treeRootFile = treePath + subSampleNameVector[i];
    if (READ_FROM_LOCAL) treeRootFile = subSampleNameVector[i];
    chain->Add(TString(treeRootFile.c_str()));

  }

  if(!chain) {
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

void EoverP_shervin::Loop(const string sampleName, const vector<Float_t> &energybinEdges, const string& dirName)
{


   if (fChain == 0) return;
   fChain->SetBranchStatus("*",1);

   gStyle->SetStatStyle(0);

   TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 

   // variables to use branch in other chain (when reading data, we have to different samples and I want those events common in both)
   Float_t PtEle_other[3], etaEle_other[3], R9Ele_other[3]; // for selection
   UInt_t eleID_other[3];
   Float_t energySCEle_must_regrCorr_ele_other[3], rawEnergySCEle_must_other[3], pGsfEle_other[3]; // for analysis
   Int_t runNumber_other, lumiBlock_other;
   ULong64_t eventNumber_other;
   TChain* otherDataChain = NULL;
   if (sampleName == "DATA" || sampleName == "DATA_w") {
     otherDataChain = new TChain("selected");
     buildOtherDataChain(otherDataChain, sampleName);
     // set branch to be used for comparison
     otherDataChain->SetBranchAddress("runNumber",&runNumber_other);
     otherDataChain->SetBranchAddress("eventNumber",&eventNumber_other);
     otherDataChain->SetBranchAddress("lumiBlock",&lumiBlock_other);
     otherDataChain->SetBranchAddress("PtEle",PtEle_other);
     otherDataChain->SetBranchAddress("etaEle",etaEle_other);
     otherDataChain->SetBranchAddress("R9Ele",R9Ele_other);
     otherDataChain->SetBranchAddress("eleID",eleID_other);
     otherDataChain->SetBranchAddress("energySCEle_must_regrCorr_ele",energySCEle_must_regrCorr_ele_other);
     otherDataChain->SetBranchAddress("rawEnergySCEle_must",rawEnergySCEle_must_other);
     if (USE_P_AT_VTX) otherDataChain->SetBranchAddress("pAtVtxGsfEle",pGsfEle_other);
     else otherDataChain->SetBranchAddress("pModeGsfEle",pGsfEle_other);
     fChain->AddFriend(otherDataChain);
   }

   Int_t nEnergyBins = energybinEdges.size() -1;

   string rootfileName = dirName + "EoverP_" + sampleName + ".root";

   TFile *rootFile = new TFile((rootfileName).c_str(),"RECREATE");
   if (!rootFile || !rootFile->IsOpen()) {
     cout << "Error: file \"" << rootfileName << "\" was not opened." << endl;
     exit(EXIT_FAILURE);
   }

   vector<TH1F*> hEoverP_energyBin(nEnergyBins,NULL);
   vector<TH1F*> hEoverP_gain6and12_energyBin(nEnergyBins,NULL);
   vector<TH1F*> hEoverP_gain12_energyBin(nEnergyBins,NULL);
   // histogram to store the mean energy in a given bin, useful to plot E/P asf of the mean energy in the bin with a TGraph
   TH1F* hMeanEnergyInEnergyBin = new TH1F("hMeanEnergyInEnergyBin","",nEnergyBins,energybinEdges.data());

   vector<TH1F*> hEcorrOverEtrue_energyBin(nEnergyBins,NULL);
   vector<TH1F*> hErawOverEtrue_energyBin(nEnergyBins,NULL);
   vector<TH1F*> hPtrackOverEtrue_energyBin(nEnergyBins,NULL);

   // for additional studies
   TH1F* hSigmaIetaIeta_lowE = new TH1F("hSigmaIetaIeta_lowE","",96,0.0,0.012);
   TH1F* hSigmaIetaIeta_bumpE = new TH1F("hSigmaIetaIeta_bumpE","",96,0.0,0.012);
   TH1F* hR9_lowE = new TH1F("hR9_lowE","",90,0.2,1.1);
   TH1F* hR9_bumpE = new TH1F("hR9_bumpE","",90,0.2,1.1);
   TH1F* hPt_lowE = new TH1F("hPt_lowE","",140,0.0,700.0);
   TH1F* hPt_bumpE = new TH1F("hPt_bumpE","",140,0.0,700.0);   

   TH1F* hElePt = new TH1F("hElePt","",200,0,1000);
   TH1F* hEleEcorr = new TH1F("hEleEcorr","",200,0,1000);
   TH1F* hEleEraw = new TH1F("hEleEraw","",200,0,1000);
   TH1F* hpAtVtxGsfEle = new TH1F("hpAtVtxGsfEle","",45,0,900);
   TH1F* hpModeGsfEle = new TH1F("hpModeGsfEle","",45,0,900);

   for (Int_t i = 0; i < nEnergyBins; i++) {
     hEoverP_energyBin[i] = new TH1F(Form("hEoverP_energyBin%1.0fTo%1.0f",energybinEdges[i],energybinEdges[i+1]),"",200,0.05,2.05);
     hEoverP_gain6and12_energyBin[i] = new TH1F(Form("hEoverP_gain6and12_energyBin%1.0fTo%1.0f",energybinEdges[i],energybinEdges[i+1]),"",200,0.05,2.05);
     hEoverP_gain12_energyBin[i] = new TH1F(Form("hEoverP_gain12_energyBin%1.0fTo%1.0f",energybinEdges[i],energybinEdges[i+1]),"",200,0.05,2.05);
   }

   if (sampleName != "DATA" && sampleName != "DATA_w") {
     for (Int_t i = 0; i < nEnergyBins; i++) {
       hEcorrOverEtrue_energyBin[i] = new TH1F(Form("hEcorrOverEtrue_energyBin%1.0fTo%1.0f",energybinEdges[i],energybinEdges[i+1]),"",200,0.55,1.55);
       hErawOverEtrue_energyBin[i] = new TH1F(Form("hErawOverEtrue_energyBin%1.0fTo%1.0f",energybinEdges[i],energybinEdges[i+1]),"",200,0.55,1.55);
       hPtrackOverEtrue_energyBin[i] = new TH1F(Form("hPtrackOverEtrue_energyBin%1.0fTo%1.0f",energybinEdges[i],energybinEdges[i+1]),"",200,0.05,2.05);
     }
   }

   TH1F* hEoP_template = new TH1F("hEoP_template","template",200,0.05,2.05);  // a template for E/P in the whole energy spectrum (not used at the moment)

   vector<Double_t> energyBinTH2 = {25.,50.,75.,100.,125.,175.,225., 275.,350., 450.,550.,650.,800.,1000.};
   Int_t NbinsTH2 = energyBinTH2.size() - 1;
   TH2F* h2EnergyPtrack = new TH2F("h2EnergyPtrack","",NbinsTH2,energyBinTH2.data(),NbinsTH2,energyBinTH2.data()); 

   // Float_t lowMeanCut = 0.85;
   // Float_t upMeanCut = 1.15;

   // used only for MC study
   /////////////////////////
   // Int_t MCtruthMatchFound = 0;
   // Float_t Etrue = 0.0; 
   /////////////////////////

   // Build index to associates events in both data samples
   //fChain->BuildIndex("runNumber","eventNumber");  // I think this is not needed
   //cout << "Check 0 \n" << flush;
   if (otherDataChain && otherDataChain != NULL ) otherDataChain->BuildIndex("runNumber","eventNumber");
   //cout << "check \n" << flush;

   // TTreeIndex *TreeIndex_fChain = (TTreeIndex*)fChain->GetTreeIndex();
   // // Long64_t* indexValuesStream = TreeIndex_fChain->GetIndexValues();
   // Long64_t* indexStream = TreeIndex_fChain->GetIndex();
   // cout << "TreeIndex_fChain->GetN() = " << TreeIndex_fChain->GetN() << endl;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   Long64_t finalEntries = 0;
   Long64_t entriesWithEnergyLowerThanFirstBin = 0;
   Long64_t entriesWithEnergyHigherThanLastBin = 0;
   Long64_t entriesWithPLowerThanFirstBin = 0;
   Long64_t entriesWithPHigherThanLastBin = 0;
   Long64_t entriesWithGain1 = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   // for (Long64_t jentry=0; jentry<TreeIndex_fChain->GetN();jentry++) {
 
     // Long64_t localStream = LoadTree(indexStream[jentry]);
     // fChain->GetEntry(localStream);

     Long64_t ientry = LoadTree(jentry);   // could load the entry corresponding to index, but I cut on entries for the other tree, so it should be the same
     if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%500000 == 0) cout << jentry << endl;
  
      if ( otherDataChain != NULL && !(runNumber == runNumber_other && eventNumber == eventNumber_other)) continue;

      // if (met_pt < 50) continue;
      // if (!( nEle10V == 1 && nEle40T == 1) ) continue;
      if ((eleID[0] & 0x0008) != 0x0008) continue;
      if ( !(PtEle[0] > 40 && fabs(etaEle[0]) < ELE_ETA_MAX) ) continue;
      if ( otherDataChain != NULL && !(PtEle_other[0] > 40 && fabs(etaEle_other[0]) < ELE_ETA_MAX && ((eleID_other[0] & 0x0008) == 0x0008) ))  continue;

      Double_t energyToUse = -1.0;
      if (USE_RAWE) energyToUse = rawEnergySCEle_must[0];
      else energyToUse = energySCEle_must_regrCorr_ele[0];

      if (!USE_E) {
	Double_t theta = 2. * atan(exp(-etaEle[0])); 
	energyToUse *= sin(theta);
      } 

      // fill these before cutting on r9 // but the skim had R9 > 0.94
      if (energyToUse > 99.9 && energyToUse < 225.0) {
	hR9_lowE->Fill(R9Ele[0]);
      } else if (energyToUse > 349.9 && energyToUse < 450.0) {
	hR9_bumpE->Fill(R9Ele[0]);
      }

      if (!(R9Ele[0] > 0.94)) continue;
      if (otherDataChain != NULL && !(R9Ele_other[0] > 0.94)) continue;

      if (energyToUse > 99.9 && energyToUse < 225.0) {
	hSigmaIetaIeta_lowE->Fill(0);
	hPt_lowE->Fill(PtEle[0]);
      } else if (energyToUse > 349.9 && energyToUse < 450.0) {
	hSigmaIetaIeta_bumpE->Fill(0);
	hPt_bumpE->Fill(PtEle[0]);
      }

      hElePt->Fill(PtEle[0]);
      hEleEcorr->Fill(energySCEle_must_regrCorr_ele[0]);  
      hEleEraw->Fill(rawEnergySCEle_must[0]); 
      hpAtVtxGsfEle->Fill(pAtVtxGsfEle[0]);
      hpModeGsfEle->Fill(pModeGsfEle[0]);

      Double_t EoverP_toUse = -1.0;
      Double_t pToUse = 0.0;
      if (USE_P_AT_VTX) pToUse = pAtVtxGsfEle[0];
      else pToUse = pModeGsfEle[0];

      Double_t varToUseForBin = -1.0;
      if (BIN_P) varToUseForBin = pToUse;
      else varToUseForBin = energyToUse;      

      EoverP_toUse = energyToUse / pToUse;

      hEoP_template->Fill(EoverP_toUse);

      // look for the bin in the proper variable      
      Int_t bin = getBinNumber(varToUseForBin,energybinEdges);  // this function returns negative value if bin not found
      //Int_t bin = (BIN_P ? getBinNumber(pToUse,energybinEdges) : getBinNumber(energyToUse,energybinEdges));  // this function returns negative value if bin not found
      
      if (bin >= 0) {

	hEoverP_energyBin[bin]->Fill(EoverP_toUse);
	if (gainEle[0] != 2) {
	  hEoverP_gain6and12_energyBin[bin]->Fill(EoverP_toUse);
	  if (gainEle[0] != 1) hEoverP_gain12_energyBin[bin]->Fill(EoverP_toUse);
	}
	//sum the energy to the bin content in the bin it belongs to (at the end we will divide by the number of entries in each bin)
	// using bin+1 because the histogram bin number goes from 1 to number of bins, while "bin" variable starts from 0
	hMeanEnergyInEnergyBin->SetBinContent(bin+1, varToUseForBin + hMeanEnergyInEnergyBin->GetBinContent(bin+1));  

      } else if (bin == -1) {
	// fill last bin with overflows to gain in statistics
	hEoverP_energyBin[nEnergyBins-1]->Fill(EoverP_toUse);
	if (gainEle[0] != 2) {
	  hEoverP_gain6and12_energyBin[nEnergyBins-1]->Fill(EoverP_toUse);
	  if (gainEle[0] != 1) hEoverP_gain12_energyBin[nEnergyBins-1]->Fill(EoverP_toUse);
	}
	// however, the mean is done without them, because overflows can be far outlier and bias the mean
	//hMeanEnergyInEnergyBin->SetBinContent(nEnergyBins, varToUseForBin + hMeanEnergyInEnergyBin->GetBinContent(nEnergyBins));
      }

      h2EnergyPtrack->Fill(pToUse,energyToUse);

      finalEntries++;
      if (energyToUse < energybinEdges[0]) entriesWithEnergyLowerThanFirstBin++;
      if (energyToUse > energybinEdges.back()) entriesWithEnergyHigherThanLastBin++;
      if (energyToUse < energybinEdges[0]) entriesWithPLowerThanFirstBin++;
      if (energyToUse > energybinEdges.back()) entriesWithPHigherThanLastBin++;
      if (gainEle[0] == 2) entriesWithGain1++;

   }  // end of loop on entries

   cout << "Loop ended: final entries = " << finalEntries << endl;
   cout << "Entries with E < " << energybinEdges[0] << " = " << entriesWithEnergyLowerThanFirstBin << endl;
   cout << "Entries with P < " << energybinEdges[0] << " = " << entriesWithPLowerThanFirstBin << endl;
   cout << "Entries with E > " << energybinEdges.back() << " = " << entriesWithEnergyHigherThanLastBin << endl;
   cout << "Entries with P > " << energybinEdges.back() << " = " << entriesWithPHigherThanLastBin << endl;
   cout << "Entries with gain 1 = " << entriesWithGain1 << endl;
   cout << "hEoP_template->GetEntries() = " << hEoP_template->GetEntries() << endl;
   cout << "hEoP_template->Integral() = " << hEoP_template->Integral() << endl;

   for (Int_t i = 0; i < nEnergyBins; i++) {

     //hEoverP_energyBin[i]->SetStats(0);  // keep stat box in file, and remove in directly the function that plots it on canvas
     // keep the following axis settings
     hEoverP_energyBin[i]->GetXaxis()->SetTitle("E / P");
     hEoverP_energyBin[i]->GetXaxis()->SetTitleSize(0.06);
     hEoverP_energyBin[i]->GetXaxis()->SetTitleOffset(0.8);
     hEoverP_energyBin[i]->GetYaxis()->SetTitle("events");
     hEoverP_energyBin[i]->GetYaxis()->SetTitleSize(0.055);
     hEoverP_energyBin[i]->GetYaxis()->SetTitleOffset(0.8);

     hEoverP_gain6and12_energyBin[i]->GetXaxis()->SetTitle("E / P");
     hEoverP_gain6and12_energyBin[i]->GetXaxis()->SetTitleSize(0.06);
     hEoverP_gain6and12_energyBin[i]->GetXaxis()->SetTitleOffset(0.8);
     hEoverP_gain6and12_energyBin[i]->GetYaxis()->SetTitle("events");
     hEoverP_gain6and12_energyBin[i]->GetYaxis()->SetTitleSize(0.055);
     hEoverP_gain6and12_energyBin[i]->GetYaxis()->SetTitleOffset(0.8);


     // remember that the last bin is somewhat special, because it holds the overflows
     // however, the mean is done without them, because overflows can be far outlier and bias the mean
     hMeanEnergyInEnergyBin->SetBinContent(i+1, hMeanEnergyInEnergyBin->GetBinContent(i+1)/hEoverP_energyBin[i]->GetEntries());
     
   }

   string whichEnergy = "";
   if (USE_RAWE) whichEnergy = "raw";
   else whichEnergy = "corrected";

   if (USE_E) {
     hMeanEnergyInEnergyBin->GetXaxis()->SetTitle((whichEnergy + " E [GeV]").c_str());
     hMeanEnergyInEnergyBin->GetYaxis()->SetTitle("mean E [GeV]");
     h2EnergyPtrack->GetYaxis()->SetTitle((whichEnergy + " E [GeV]").c_str());  //the lin below can use P, but this is E (because it is E vs P)
   } else {
     hMeanEnergyInEnergyBin->GetXaxis()->SetTitle((whichEnergy + " E_{T} [GeV]").c_str());
     hMeanEnergyInEnergyBin->GetYaxis()->SetTitle("mean E_{T} [GeV]");
     h2EnergyPtrack->GetYaxis()->SetTitle((whichEnergy + " E_{T} [GeV]").c_str());
   }
   if (BIN_P) {
     hMeanEnergyInEnergyBin->GetXaxis()->SetTitle("P [GeV]");
     hMeanEnergyInEnergyBin->GetYaxis()->SetTitle("mean P [GeV]");
   }


   hMeanEnergyInEnergyBin->GetXaxis()->SetTitleSize(0.06);
   hMeanEnergyInEnergyBin->GetXaxis()->SetTitleOffset(0.8);
   hMeanEnergyInEnergyBin->GetYaxis()->SetTitleSize(0.055);
   hMeanEnergyInEnergyBin->GetYaxis()->SetTitleOffset(0.8);

   h2EnergyPtrack->GetXaxis()->SetTitle("p(track) [GeV]");
   h2EnergyPtrack->GetXaxis()->SetTitleSize(0.06);
   h2EnergyPtrack->GetXaxis()->SetTitleOffset(0.8);
   h2EnergyPtrack->GetYaxis()->SetTitleSize(0.06);
   h2EnergyPtrack->GetYaxis()->SetTitleOffset(0.8);

   rootFile->Write();
   rootFile->Close();
   delete rootFile;


}

//=========================================================================

void getTexSampleName(const string &SampleName, string &texSampleName) {
  
  if (SampleName == "WJetsToLNu") texSampleName = "W(l#nu)+jets";
  else if (SampleName == "DYJetsToLL") texSampleName = "Z(ll)+jets";
  else if (SampleName == "DATA_w") texSampleName = "data(weight)";

}

//=========================================================================

void plotDistribution(const string &sampleName, const vector<Float_t> &energybinEdges, TH1F* hPeak, TH1F* hSigma, const string &hNameID, const string &dirName, TH1F* hPeakMeanInRange = NULL, const Double_t lowRangeValue = 0.0, const Double_t upRangeValue = 2.0, TH1F* hPeakShift_MCdata = NULL, TH1F* hScaleFactor_MCdata = NULL, const string plotNameTag = "") {

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

  if (DO_TEMPLATE_FIT) {

    if (sampleName == "DATA" || sampleName == "DATA_w") {

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

  }

  UInt_t nBins = energybinEdges.size() -1;

  for (UInt_t i = 0; i < nBins; i++) {

    cout << " Energy bin: [" << energybinEdges[i] << ", " << energybinEdges[i+1] << "]" << endl;

    htmp = (TH1F*)f->Get(Form("h%s_energyBin%1.0fTo%1.0f",(hNameID+plotNameTag).c_str(),energybinEdges[i],energybinEdges[i+1]));  
    if (!htmp) {
      cout << "Error: histogram " << Form("h%s_energyBin%1.0fTo%1.0f",(hNameID+plotNameTag).c_str(),energybinEdges[i],energybinEdges[i+1]) 
	   << " not found in file ' " << fileName << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
    // now, if using data I will fit with template obtained from the same distribution in MC file. Since the histogram's name is the same in both file, I 
    // got a strange behaviour, that is, when using fMC->Get(<name>) I actually get the same histogram in data. A solution is deleting the previous histogram
    // got from data file before getting that for MC. In this case, better to give the clone a different name to make explicit it is for data 
    if (sampleName == "DATA" || sampleName == "DATA_w") hist = (TH1F*)htmp->Clone(Form("%s_data",htmp->GetName()));
    else hist = (TH1F*)htmp->Clone();

    //    hist->SetStats(0);
    gStyle->SetOptStat(10);
    gStyle->SetOptFit(112);
    hist->Draw("HE");
    if (sampleName == "DATA" || sampleName == "DATA_w") {
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
      gaussEdgeL = 0.8;
      gaussEdgeR = 1.2;
    } else if (energybinEdges[i] > 449.9) {
      gaussEdgeL = 0.85;
      gaussEdgeR = 1.15;
    }
    // if MC, use more specific settings (decided after looking at plots, there is no a priori reason for them)
    if (sampleName != "DATA" && sampleName != "DATA_w") {
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

    if (DO_TEMPLATE_FIT) {

      if ( (sampleName == "DATA" || sampleName == "DATA_w") && hNameID == "EoverP") {

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
      if (FIT_2SIDE_CB) stat->SetY1NDC(stat->GetY2NDC() - 1.8 * (stat->GetY2NDC() - stat->GetY1NDC()));     
      stat->Draw();
    }
    // Box for fit with template (data only)
    TPaveStats *stat1 = NULL;
    if (DO_TEMPLATE_FIT) {
      if (sampleName == "DATA" || sampleName == "DATA_w") {
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

    if (BIN_P) {
      whichEnergy = "track P";
      hist->SetTitle(Form("%1.0f < %s [GeV] < %1.0f",energybinEdges[i],whichEnergy.c_str(),energybinEdges[i+1]));
      c->SaveAs(Form("%s%sdistribution_P%1.0fTo%1.0f_%s.pdf",dirName.c_str(),(hNameID+plotNameTag).c_str(),energybinEdges[i],energybinEdges[i+1],sampleName.c_str()));
      c->SaveAs(Form("%s%sdistribution_P%1.0fTo%1.0f_%s.png",dirName.c_str(),(hNameID+plotNameTag).c_str(),energybinEdges[i],energybinEdges[i+1],sampleName.c_str()));
    } else if (USE_E) {      
      hist->SetTitle(Form("%1.0f < %s [GeV] < %1.0f",energybinEdges[i],whichEnergy.c_str(),energybinEdges[i+1]));
      c->SaveAs(Form("%s%sdistribution_E%1.0fTo%1.0f_%s.pdf",dirName.c_str(),(hNameID+plotNameTag).c_str(),energybinEdges[i],energybinEdges[i+1],sampleName.c_str()));
      c->SaveAs(Form("%s%sdistribution_E%1.0fTo%1.0f_%s.png",dirName.c_str(),(hNameID+plotNameTag).c_str(),energybinEdges[i],energybinEdges[i+1],sampleName.c_str()));
    } else {
      hist->SetTitle(Form("%1.0f < %s_{T}[GeV] < %1.0f",energybinEdges[i],whichEnergy.c_str(),energybinEdges[i+1]));
      c->SaveAs(Form("%s%sdistribution_ET%1.0fTo%1.0f_%s.pdf",dirName.c_str(),(hNameID+plotNameTag).c_str(),energybinEdges[i],energybinEdges[i+1],sampleName.c_str()));
      c->SaveAs(Form("%s%sdistribution_ET%1.0fTo%1.0f_%s.png",dirName.c_str(),(hNameID+plotNameTag).c_str(),energybinEdges[i],energybinEdges[i+1],sampleName.c_str()));
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
  subpad_1->SetGridy();
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
  } else if (hNameToTriggerOptions == "gain12") {
    hdata->SetMinimum(0.95);
    hdata->SetMaximum(1.02);
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
  getTexSampleName(MCSampleName, texMCSampleName);

  leg->AddEntry(hdata,"data(multifit)","lf");
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
  if (MCSampleName == "DATA_w") ratioplot->GetYaxis()->SetTitle("multifit / weight");
  else ratioplot->GetYaxis()->SetTitle("data / MC");
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
  ratioplotCopy->SetMaximum(maximumYaxisValue);
  ratioplotCopy->SetMinimum(minimumYaxisValue);
  if (hNameToTriggerOptions == "gain12") {
    ratioplotCopy->SetMinimum(0.95);
    ratioplotCopy->SetMaximum(1.05);
  }

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

    string outfname = "";
    if (MCSampleName == "DATA_w") outfname = dirName + "EoverP_vs" + tmpString + "Energy_dataDataWgtRatio.root";
    else outfname = dirName + "EoverP_vs" + tmpString + "Energy_dataMCRatio.root";
    TFile *f = new TFile( outfname.c_str(),"RECREATE");
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

//========================================================
 
 void drawMap(TH2F* h2, const string dirName, const string canvasName) {

  TCanvas *c = new TCanvas("c","",700,700);
  c->SetLogz();
  h2->SetStats(0);
  h2->Draw("COLZ");
  h2->GetXaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetTitleOffset(0.8);
  h2->GetYaxis()->SetTitleSize(0.035);
  h2->GetYaxis()->CenterTitle();
  h2->GetYaxis()->SetTitleOffset(1.3);
  gPad->Update();

  TProfile* mapProfile = h2->ProfileX(Form("%s_pfx",h2->GetName()));
  mapProfile->SetMarkerColor(kBlack);
  mapProfile->SetMarkerStyle(20);
  mapProfile->SetMarkerSize(1);
  mapProfile->Draw("EPsame");

  TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  //  TLegendEntry *l1 = leg->AddEntry((TObject*)0,Form("Correlation = %.2f",h2->GetCorrelationFactor()),"");
  leg->SetHeader(Form("Correlation = %.2f",h2->GetCorrelationFactor()));
  leg->SetTextColor(0);
  leg->Draw("same");

  TPaletteAxis *palette = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");
  if (!palette || palette == NULL) {
    cout << "Error in function drawMap(): palette not found. ABORT" << endl;
    exit(EXIT_FAILURE);
  }
  // the following lines move the palette. Choose the values you need for the position.                                                                              
  palette->SetX1NDC(0.91);
  palette->SetX2NDC(0.94);
  //palette->SetY1NDC(0.2);                                                                                                                                          
  //palette->SetY2NDC(0.8);                                                                                                                                          
  gPad->Modified();
  gPad->Update();

  c->SaveAs( (dirName + canvasName + ".pdf").c_str() );
  c->SaveAs( (dirName + canvasName + ".png").c_str() );

  delete c;
  delete leg;



}

//==========================================================================================

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

  TH2F* h2var = NULL;  // to get histogram from file
  TH2F* h2var2 = NULL;  // to get histogram from file


  vector<string> hNameList;
  hNameList.push_back("hSigmaIetaIeta_lowE");
  hNameList.push_back("hSigmaIetaIeta_bumpE");
  hNameList.push_back("hR9_lowE");
  hNameList.push_back("hR9_bumpE");
  hNameList.push_back("hPt_lowE");
  hNameList.push_back("hPt_bumpE");
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
  canvasName.push_back("sigmaIetaIeta_bumpE");
  canvasName.push_back("R9_lowE");
  canvasName.push_back("R9_bumpE");
  canvasName.push_back("pT_lowE");
  canvasName.push_back("pT_bumpE");
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

  TH2F* h2Data = NULL;
  TH2F* h2MC = NULL;

  //cout << "check here" << endl;

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
    hData->Scale(1./hData->GetEntries());
    hMC->Sumw2();
    hMC->Scale(1./hMC->GetEntries());
    drawPlotDataMC(hData, hMC, MCSampleName, xAxisName[i], "a.u.",canvasName[i],dirName,canvasTitle[i],"distributions",1);
    //////////////////////////

  }


  fData->cd();
  // reading data file
  h2var = (TH2F*)fData->Get("h2EnergyPtrack");
  if (!h2var) {
    cout << "Error: histogram not found in file ' " << fileNameData << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
  }
  h2Data = (TH2F*) h2var->Clone();

  fMC->cd();
  // reading MC file
  h2var2 = (TH2F*)fMC->Get("h2EnergyPtrack");
  if (!h2var2) {
    cout << "Error: histogram not found in file ' " << fileNameMC << "'. End of programme." << endl;
    exit(EXIT_FAILURE);
  }
  h2MC = (TH2F*) h2var2->Clone();

  string canvasTH2Name = "";
  if (USE_RAWE) canvasTH2Name = "ErawPtrack";
  else canvasTH2Name = "EcorrPtrack";
  drawMap(h2Data,dirName,(canvasTH2Name+"_multifit").c_str());
  drawMap(h2MC,dirName,(canvasTH2Name+"_weight").c_str());

}


//=================================================================================

void plotFromFit(const string &dataSampleName, const string &dataWgtSampleName, const vector<Float_t> &energybinEdges, const string & dirName) {

  TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 
  
  string xAxisName = "";

  string whichEnergy = "";
  if (USE_RAWE) whichEnergy = "raw";
  else whichEnergy = "corrected";

  if (USE_E) xAxisName = whichEnergy + " E [GeV]";
  else xAxisName = whichEnergy + " E_{T} [GeV]";

  xAxisName = "P(track) [GeV]";
  Int_t nEnergyBins = energybinEdges.size() -1;

  TH1F* hPeakEoverPdata = new TH1F("hPeakEoverPdata","",nEnergyBins,energybinEdges.data());
  TH1F* hSigmaEoverPdata = new TH1F("hSigmaEoverPdata","",nEnergyBins,energybinEdges.data());
  TH1F* hPeakEoverPdataWgt = new TH1F("hPeakEoverPdataWgt","",nEnergyBins,energybinEdges.data());
  TH1F* hSigmaEoverPdataWgt = new TH1F("hSigmaEoverPdataWgt","",nEnergyBins,energybinEdges.data());

  // valuse used to define the range when getting peak position from TH1::GetMean()
  double_t lowRangeValueForMean = 0.8; 
  Double_t upRangeValueForMean = 1.2; 
  string titleForPlotsWithMeanInRange = string(Form("peak from TH1F::GetMean() in [%1.1f, %1.1f]",lowRangeValueForMean,upRangeValueForMean));

  TH1F* hPeakEoverPdata_meanInRange = new TH1F("hPeakEoverPdata_meanInRange","",nEnergyBins,energybinEdges.data());
  TH1F* hPeakEoverPdataWgt_meanInRange = new TH1F("hPeakEoverPdataWgt_meanInRange","",nEnergyBins,energybinEdges.data());

  plotDistribution(dataSampleName, energybinEdges, hPeakEoverPdata, hSigmaEoverPdata, "EoverP",dirName, 
		   hPeakEoverPdata_meanInRange,lowRangeValueForMean, upRangeValueForMean);
  plotDistribution(dataWgtSampleName, energybinEdges, hPeakEoverPdataWgt, hSigmaEoverPdataWgt, "EoverP",dirName,
		   hPeakEoverPdataWgt_meanInRange,lowRangeValueForMean, upRangeValueForMean);

  string canvasTitle = "";
  if (FIT_2SIDE_CB) canvasTitle = "fit with double Crystal Ball";
  else canvasTitle = "fit with Crystal Ball";

  drawPlotDataMC(hPeakEoverPdata, hPeakEoverPdataWgt, dataWgtSampleName, xAxisName, "peak(E/P) from fit", "modeEoverPfromFit",dirName,canvasTitle);
  drawPlotDataMC(hSigmaEoverPdata, hSigmaEoverPdataWgt, dataWgtSampleName, xAxisName, "#sigma(E/P) from fit", "sigmaEoverPfromFit",dirName,canvasTitle);
  if ( (hPeakEoverPdata_meanInRange != NULL) && (hPeakEoverPdataWgt_meanInRange != NULL) ) {
    drawPlotDataMC(hPeakEoverPdata_meanInRange, hPeakEoverPdataWgt_meanInRange, dataWgtSampleName, xAxisName, "peak(E/P)", "modeEoverP_usingMean",dirName,
		   titleForPlotsWithMeanInRange);
  }

  // study only gain 6 and 12 together
  TH1F* hPeakEoverP_gain6and12_data = new TH1F("hPeakEoverP_gain6and12_data","",nEnergyBins,energybinEdges.data());
  TH1F* hPeakEoverP_gain6and12_dataWgt = new TH1F("hPeakEoverP_gain6and12_dataWgt","",nEnergyBins,energybinEdges.data());
  TH1F* hPeakEoverP_gain6and12_data_meanInRange = new TH1F("hPeakEoverPdata_meanInRange","",nEnergyBins,energybinEdges.data());
  TH1F* hPeakEoverP_gain6and12_dataWgt_meanInRange = new TH1F("hPeakEoverPdataWgt_meanInRange","",nEnergyBins,energybinEdges.data());

  plotDistribution(dataSampleName, energybinEdges, hPeakEoverP_gain6and12_data, NULL, "EoverP",dirName, 
		   hPeakEoverP_gain6and12_data_meanInRange,lowRangeValueForMean, upRangeValueForMean, NULL, NULL, "_gain6and12");
  plotDistribution(dataWgtSampleName, energybinEdges, hPeakEoverP_gain6and12_dataWgt, NULL, "EoverP",dirName,
		   hPeakEoverP_gain6and12_dataWgt_meanInRange,lowRangeValueForMean, upRangeValueForMean, NULL, NULL, "_gain6and12");
  drawPlotDataMC(hPeakEoverP_gain6and12_data, hPeakEoverP_gain6and12_dataWgt, dataWgtSampleName, 
		 xAxisName, "peak(E/P) from fit", "modeEoverP_gain6and12_fromFit",dirName,canvasTitle);
  if ( (hPeakEoverPdata_meanInRange != NULL) && (hPeakEoverPdataWgt_meanInRange != NULL) ) {
    drawPlotDataMC(hPeakEoverP_gain6and12_data_meanInRange, hPeakEoverP_gain6and12_dataWgt_meanInRange, 
		   dataWgtSampleName, xAxisName, "peak(E/P)", "modeEoverP_gain6and12_usingMean",dirName,
		   titleForPlotsWithMeanInRange);
  }

  // study only 12
  TH1F* hPeakEoverP_gain12_data = new TH1F("hPeakEoverP_gain12_data","",nEnergyBins,energybinEdges.data());
  TH1F* hPeakEoverP_gain12_dataWgt = new TH1F("hPeakEoverP_gain12_dataWgt","",nEnergyBins,energybinEdges.data());
  TH1F* hPeakEoverP_gain12_data_meanInRange = new TH1F("hPeakEoverPdata_meanInRange","",nEnergyBins,energybinEdges.data());
  TH1F* hPeakEoverP_gain12_dataWgt_meanInRange = new TH1F("hPeakEoverPdataWgt_meanInRange","",nEnergyBins,energybinEdges.data());

  plotDistribution(dataSampleName, energybinEdges, hPeakEoverP_gain12_data, NULL, "EoverP",dirName, 
		   hPeakEoverP_gain12_data_meanInRange,lowRangeValueForMean, upRangeValueForMean, NULL, NULL, "_gain12");
  plotDistribution(dataWgtSampleName, energybinEdges, hPeakEoverP_gain12_dataWgt, NULL, "EoverP",dirName,
		   hPeakEoverP_gain12_dataWgt_meanInRange,lowRangeValueForMean, upRangeValueForMean, NULL, NULL, "_gain12");
  drawPlotDataMC(hPeakEoverP_gain12_data, hPeakEoverP_gain12_dataWgt, dataWgtSampleName, 
		 xAxisName, "peak(E/P) from fit", "modeEoverP_gain12_fromFit",dirName,canvasTitle,"gain12");
  if ( (hPeakEoverPdata_meanInRange != NULL) && (hPeakEoverPdataWgt_meanInRange != NULL) ) {
    drawPlotDataMC(hPeakEoverP_gain12_data_meanInRange, hPeakEoverP_gain12_dataWgt_meanInRange, 
		   dataWgtSampleName, xAxisName, "peak(E/P)", "modeEoverP_gain12_usingMean",dirName,
		   titleForPlotsWithMeanInRange,"gain12");
  }

}


//======================================================================

Int_t main(Int_t argc, char* argv[]) {

  Int_t doAll_flag = 1;
  Int_t doLoop_flag = 1;
  Int_t doLoopMC_flag = 1;
  Int_t doLoopData_flag = 1; 
  Int_t doPlot_flag = 1;

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
  sampleName.push_back("DATA_w");

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
    //energybinEdges.push_back(650.0);
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
    
      if (!doLoopMC_flag &&  ( (sampleName[i] == "DATA") || (sampleName[i] == "DATA_w") )  ) continue;
      if (  !doLoopData_flag && ( (sampleName[i] == "DATA") || (sampleName[i] == "DATA_w") )  ) continue;

      // create chain                                                                                                         
    
      TChain* chain = new TChain("selected");

      buildChain(chain, sampleName[i]);
    
      if(!chain) {
	cout << "Error: chain not created. End of programme" << endl;
	exit(EXIT_FAILURE);
      }

      EoverP_shervin tree(chain);
      tree.Loop(sampleName[i], energybinEdges, dirName);
      delete chain;

    }

  }

  //  cout << "check" << endl;

  if(doAll_flag || doPlot_flag) {

    plotFromFit(sampleName[0],sampleName[1], energybinEdges, dirName);
    plotHistogramsFromFile(sampleName[0],sampleName[1], dirName);

  }

  return 0;

}



////////////////////////////////////////
//TO BE IMPLEMENTED
//
//Plots of E/P asaf mean energy
//
