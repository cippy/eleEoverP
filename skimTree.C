#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TList.h>

using namespace std;

void buildChain(TChain* chain, string sampleName) {
  
  cout << "Creating chain ..." << endl;
  
  vector<string> subSampleNameVector;

  if (sampleName == "DATA") {

    subSampleNameVector.push_back("SingleElectron-Run2016B-WSkim-Prompt_v2/273150-275376/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016B-WSkim-Prompt_v2-273150-275376.root");
    subSampleNameVector.push_back("SingleElectron-Run2016C-WSkim-Prompt_v2/275420-276283/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016C-WSkim-Prompt_v2-275420-276283.root");
    subSampleNameVector.push_back("SingleElectron-Run2016D-WSkim-Prompt_v2/276315-276811/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016D-WSkim-Prompt_v2-276315-276811.root");
    subSampleNameVector.push_back("SingleElectron-Run2016E-WSkim-Prompt-v2/276830-277420/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016E-WSkim-Prompt-v2-276830-277420.root");
    subSampleNameVector.push_back("SingleElectron-Run2016F-WSkim-Prompt-v1/277820-278808/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016F-WSkim-Prompt-v1-277820-278808.root");
    subSampleNameVector.push_back("SingleElectron-Run2016G-WSkim-Prompt-v1/278817-279588/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016G-WSkim-Prompt-v1-278817-279588.root");
    subSampleNameVector.push_back("SingleElectron-Run2016G-WSkim-Prompt-v1/279589-279931/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016G-WSkim-Prompt-v1-279589-279931.root");
    
  } else if (sampleName == "DATA_w") {

    subSampleNameVector.push_back("SingleElectron-Run2016B-WSkim-Prompt_v2-weightsReco/273150-275376/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016B-WSkim-Prompt_v2-weightsReco-273150-275376.root");
    subSampleNameVector.push_back("SingleElectron-Run2016C-WSkim-Prompt_v2-weightsReco/275420-276283/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016C-WSkim-Prompt_v2-weightsReco-275420-276283.root");
    subSampleNameVector.push_back("SingleElectron-Run2016D-WSkim-Prompt_v2-weightsReco/276315-276811/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016D-WSkim-Prompt_v2-weightsReco-276315-276811.root");
    subSampleNameVector.push_back("SingleElectron-Run2016E-WSkim-Prompt-v2-weightsReco/276830-277420/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016E-WSkim-Prompt-v2-weightsReco-276830-277420.root");
    subSampleNameVector.push_back("SingleElectron-Run2016F-WSkim-Prompt-v1-weightsReco/277820-278808/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016F-WSkim-Prompt-v1-weightsReco-277820-278808.root");
    subSampleNameVector.push_back("SingleElectron-Run2016G-WSkim-Prompt-v1-weightsReco/278817-279588/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016G-WSkim-Prompt-v1-weightsReco-278817-279588.root");
    subSampleNameVector.push_back("SingleElectron-Run2016G-WSkim-Prompt-v1-weightsReco/279589-279931/271036_279931-Cal_Sep2016/withExtras/SingleElectron-Run2016G-WSkim-Prompt-v1-weightsReco-279589-279931.root");
  
  } else {

    cout << "Error: unknown sampleName " << sampleName <<". End of programme" << endl;
    exit(EXIT_FAILURE);

  }

  string treePath = "";
  treePath = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/ecalelf/ntuples/13TeV/ALCARERECO/Cal_Sep2016_final_v3/";

  for(UInt_t i = 0; i < subSampleNameVector.size(); i++) {
  
    string treeRootFile = treePath + subSampleNameVector[i];
    chain->Add(TString(treeRootFile.c_str()));

  }

  if(!chain) {
    cout << "Error: chain not created. End of programme" << endl;
    exit(EXIT_FAILURE);
  }
  cout << chain->GetEntries() << endl;

}

//=================================================0

void doSkimTree(const string sampleName = "", const std::string selection = "") {

  TChain* chain = new TChain("selected");

  buildChain(chain, sampleName);
    
  if(!chain) {
    cout << "Error: chain not built. End of programme" << endl;
    exit(EXIT_FAILURE);
  }

  std::string outfileName = "/u2/mciprian/TREES_shervin_skim/" + sampleName + "_shervin.root";

  TFile* outfile = new TFile(outfileName.c_str(), "RECREATE");
  if (!outfile || outfile->IsZombie()) {
    std::cout << "Cannot open file " << outfileName << std::endl;
    exit(EXIT_FAILURE);
  }

  TTree* outtree = chain->CopyTree(selection.c_str());

  cout << "Entries before skim = " << chain->GetEntries() << endl;
  cout << "Entries after skim = " << outtree->GetEntries() << "\t(" << (Double_t) outtree->GetEntries()/chain->GetEntries() << " % efficiency)" << endl;

  outfile->Write();
  outfile->Close();

  delete chain;

}

//======================================================0

void getChainEntries(const std::string selection = "") {

  vector<string> sampleName;
  sampleName.push_back("DATA");
  sampleName.push_back("DATA_w");

  for (UInt_t i = 0; i < sampleName.size(); i++) {

    TChain* chain = new TChain("selected");
    buildChain(chain, sampleName[i]);
    
    if(!chain) {
      cout << "Error: chain not built. End of programme" << endl;
      exit(EXIT_FAILURE);
    }
    
    cout << sampleName[i] << " :\t N(events) = " << chain->GetEntries(selection.c_str()) << endl;
    
    delete chain;

  }

}


//=============================

void skimTree(const std::string selection = "") {

  doSkimTree("DATA", selection);
  doSkimTree("DATA_w", selection);


}
