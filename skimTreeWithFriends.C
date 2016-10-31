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

void buildChainWithFriend(TChain* chain, TChain* chainevfr, TChain* chainsffr, string sampleName) {
  
  cout << "Creating chain ..." << endl;
  
  vector<string> subSampleNameVector;

  if (sampleName == "DATA") {

    // 2016 2.6fb^-1
    subSampleNameVector.push_back("SingleElectron_Run2016B_PromptReco_v1_runs_272023_273146");
    subSampleNameVector.push_back("SingleElectron_Run2016B_PromptReco_v2_runs_273150_275376");
    subSampleNameVector.push_back("SingleElectron_Run2016C_PromptReco_v2_runs_275420_276283");
    subSampleNameVector.push_back("SingleElectron_Run2016D_PromptReco_v2_runs_276315_276811");
    subSampleNameVector.push_back("SingleElectron_Run2016E_PromptReco_v2_runs_276830_277420");
    subSampleNameVector.push_back("SingleElectron_Run2016F_PromptReco_v1_runs_277820_278808");
    subSampleNameVector.push_back("SingleElectron_Run2016G_PromptReco_v1_runs_278817_279931");	
    subSampleNameVector.push_back("SingleElectron_Run2016H_PromptReco_v1_runs_281085_281201");	

  } else if (sampleName == "WJetsToLNu") {

    subSampleNameVector.push_back("WJetsToLNu_HT100to200");
    subSampleNameVector.push_back("WJetsToLNu_HT100to200_ext");
    subSampleNameVector.push_back("WJetsToLNu_HT200to400");
    subSampleNameVector.push_back("WJetsToLNu_HT200to400_ext");
    subSampleNameVector.push_back("WJetsToLNu_HT400to600");
    subSampleNameVector.push_back("WJetsToLNu_HT400to600_ext");
    subSampleNameVector.push_back("WJetsToLNu_HT600to800");
    subSampleNameVector.push_back("WJetsToLNu_HT800to1200");
    subSampleNameVector.push_back("WJetsToLNu_HT800to1200_ext");
    subSampleNameVector.push_back("WJetsToLNu_HT1200to2500");
    subSampleNameVector.push_back("WJetsToLNu_HT1200to2500_ext");
    subSampleNameVector.push_back("WJetsToLNu_HT2500toInf");
    subSampleNameVector.push_back("WJetsToLNu_HT2500toInf_ext");

  } else {

    cout << "Error: unknown sampleName " << sampleName <<". End of programme" << endl;
    exit(EXIT_FAILURE);

  }

  //2016 trees
  string treePath = "root://eoscms//eos/cms/store/group/phys_exotica/monojet/mciprian/trees_80X/TREES_1TIGHTELE30SKIM_4EoP/";

  for(UInt_t i = 0; i < subSampleNameVector.size(); i++) {
  
    string treeRootFile = treePath + subSampleNameVector[i] + "_treeProducerDarkMatterMonoJet_tree.root";
    string friend_treeRootFile = treePath + "evVarFriend_" + subSampleNameVector[i]+ ".root";
    string sffriend_treeRootFile = treePath + "sfFriend_" + subSampleNameVector[i]+ ".root";

    chain->Add(TString(treeRootFile.c_str()));
    chainevfr->Add(TString(friend_treeRootFile.c_str()));
    if (sampleName != "DATA") chainsffr->Add(TString(sffriend_treeRootFile.c_str()));

  }

  // attaching friends to main chain because will skim main tree using also variables in friends
  cout << "Adding friend to chain ..." << endl;
  chain->AddFriend(chainevfr);  //adding whole friend chain as friend                                                    
  if (sampleName != "DATA") chain->AddFriend(chainsffr);  //adding whole friend chain as friend                                    

  if(!chain || !chainevfr) {
    cout << "Error: chain not created. End of programme" << endl;
    exit(EXIT_FAILURE);
  }
  if(sampleName != "DATA" && !chainsffr) {
    cout << "Error: chain not created. End of programme" << endl;
    exit(EXIT_FAILURE);
  }

  cout << chain->GetEntries() << endl;
  cout << chainevfr->GetEntries() << endl;
  if (sampleName != "DATA") cout << chainsffr->GetEntries() << endl;

}

//=================================================0

void doSkimTree(const string sampleName = "") {

  TChain* chain = new TChain("tree");
  TChain* chainevfr = new TChain("mjvars/t");
  TChain* chainsffr = new TChain("sf/t");

  cout << "sampleName: " << sampleName << endl;

  buildChainWithFriend(chain, chainevfr, chainsffr, sampleName);
    
  if (!chain) {
    cout << "Error: chain not built. End of programme" << endl;
    exit(EXIT_FAILURE);
  }
  if (!chainevfr) {
    cout << "Error: chainevfr not built. End of programme" << endl;
    exit(EXIT_FAILURE);
  }
  if (sampleName != "DATA" && !chainsffr) {
    cout << "Error: chainsffr not built. End of programme" << endl;
    exit(EXIT_FAILURE);
  }

  cout << "main tree" << endl;

  std::string outfileName = "/u2/mciprian/TREES_1tightEle30_skimAgain_exactly1looseAndTightEle_met50/" + sampleName + ".root";

  TFile* outfile = new TFile(outfileName.c_str(), "RECREATE");
  if (!outfile || outfile->IsZombie()) {
    std::cout << "Cannot open file " << outfileName << std::endl;
    exit(EXIT_FAILURE);
  }

  TTree* outtree = chain->CopyTree("met_pt > 50 && nEle10V == 1 && nEle40T == 1 && fabs(LepGood_pdgId[0]) == 11");

  // used to skim friends
  Long64_t nentries = chain->GetEntries();

  cout << "Entries before skim = " << nentries << endl;
  cout << "Entries after skim = " << outtree->GetEntries() << "\t(" << (Double_t) outtree->GetEntries()/nentries << " % efficiency)" << endl;

  outfile->Write();
  outfile->Close();

  //////////////////////////////
  // now skimming friends
  //////////////////////////////

  Float_t met_pt;
  Float_t nEle10V;
  Float_t nEle40T;
  Int_t LepGood_pdgId[4];

  chain->SetBranchAddress("met_pt",&met_pt);
  chain->SetBranchAddress("nEle10V",&nEle10V);
  chain->SetBranchAddress("nEle40T",&nEle40T);
  chain->SetBranchAddress("LepGood_pdgId",LepGood_pdgId);

  ///////////////////////
  // evVarFriend

  cout << "evVarFriend" << endl;

  std::string outfileFriendName = "/u2/mciprian/TREES_1tightEle30_skimAgain_exactly1looseAndTightEle_met50/evVarFriend_" + sampleName + ".root";

  TFile* outfileFriend = new TFile(outfileFriendName.c_str(), "RECREATE");
  if (!outfileFriend || outfileFriend->IsZombie()) {
    std::cout << "Cannot open file " << outfileFriendName << std::endl;
    exit(EXIT_FAILURE);
  }
  TDirectory *dir = outfileFriend->mkdir("mjvars");
  dir->cd();
  TTree* outfriend = chainevfr->CloneTree(0);

  for (Long64_t i=0; i<nentries; i++) {
    chain->GetEntry(i);
    if (met_pt > 50 && nEle10V == 1 && nEle40T == 1 && fabs(LepGood_pdgId[0]) == 11) {
      chainevfr->GetEntry(i);
      outfriend->Fill();
    }
  }

  outfileFriend->Write();
  outfileFriend->Close();

  ///////////////////////
  // sfFriend

  cout << "sfFriend" << endl;

  if (sampleName != "DATA") {

    std::string outfileSfFriendName = "/u2/mciprian/TREES_1tightEle30_skimAgain_exactly1looseAndTightEle_met50/sfFriend_" + sampleName + ".root";

    TFile* outfileSfFriend = new TFile(outfileSfFriendName.c_str(), "RECREATE");
    if (!outfileSfFriend || outfileSfFriend->IsZombie()) {
      std::cout << "Cannot open file " << outfileSfFriendName << std::endl;
      exit(EXIT_FAILURE);
    }
    TDirectory *dirsf = outfileFriend->mkdir("sf");
    dirsf->cd();
    TTree* outsffriend = chainsffr->CloneTree(0);

    for (Long64_t i=0; i<nentries; i++) {
      chain->GetEntry(i);
      if (met_pt > 50 && nEle10V == 1 && nEle40T == 1 && fabs(LepGood_pdgId[0]) == 11) {
	chainsffr->GetEntry(i);
	outsffriend->Fill();
      }
    }

    outfileSfFriend->Write();
    outfileSfFriend->Close();

  }

  delete chain;
  delete chainevfr;
  delete chainsffr;


}

//======================================================0

void getChainEntries(const std::string selection = "") {

  vector<string> sampleName;
  sampleName.push_back("DATA");
  sampleName.push_back("WJetsToLNu");

  for (UInt_t i = 0; i < sampleName.size(); i++) {

    TChain* chain = new TChain("tree");
    TChain* chainevfr = new TChain("mjvars/t");
    TChain* chainsffr = new TChain("sf/t");
    buildChainWithFriend(chain, chainevfr, chainsffr, sampleName[i]) ;
    
    if(!chain || !chainevfr || (sampleName[i] != "DATA" && !chainsffr)) {
      cout << "Error: chain not built. End of programme" << endl;
      exit(EXIT_FAILURE);
    }
    
    cout << sampleName[i] << " :\t N(events) = " << chain->GetEntries(selection.c_str()) << endl;
    
    delete chain;

  }

}


//=============================

void skimTreeWithFriends() {

  //  doSkimTree("DATA");
  doSkimTree("WJetsToLNu");


}
