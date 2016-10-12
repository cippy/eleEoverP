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


void makeEcorrWeight(const string &dataSampleName = "DATA", const string &MCSampleName = "WJetsToLNu", const string& dirName = "plot/2016_singleEleRunBtoG_skim1lep1jet_usingE_fitTemplate/") {

  TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2()                         

  string fileNameData = dirName + "EoverP_" + dataSampleName + ".root";
  string fileNameMC = dirName + "EoverP_" + MCSampleName + ".root";

  TH1F* hvar = NULL;  // to get histogram from file  
  TH1F* hvar2 = NULL;  // to get histogram from file

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

  fData->cd();
  // reading data file                                                                                                                                                
  hvar = (TH1F*)fData->Get("hEleEcorr");
  if (!hvar) {
    cout << "Error: histogram not found in file ' " << fileNameData << "'. End of programme." << endl;
    exit(EXIT_FAILURE);
  }
  hData = (TH1F*) hvar->Clone();

  fMC->cd();
  // reading MC file                                                                                                                                                  
  hvar2 = (TH1F*)fMC->Get("hEleEcorr");
  if (!hvar2) {
    cout << "Error: histogram not found in file ' " << fileNameMC << "'. End of programme." << endl;
    exit(EXIT_FAILURE);
  }
  hMC = (TH1F*) hvar2->Clone();

  /////////////////////////                                                                                                                                           
  // Now creating ratio      
  ////////////////////////                                                                                                                                            
  // to compare, normalize MC to same area                                                                                                                            
  // since no weights are applied, I can use GetEntries() to get integral (this includes underflow and overflow events)                                               
  hData->Sumw2();
  hData->Scale(1./hData->Integral());
  hMC->Sumw2();
  hMC->Scale(1./hMC->Integral());

  string outname = dirName + "hEcorrWeight_dataMCratio.root";
  TFile* fout = TFile::Open(outname.c_str(),"RECREATE");
  if (!fout || !fout->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file hEcorrWeight_dataMCratio.\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  TH1F* hdataMCratio = (TH1F*) hData->Clone("hdataMCratio");
  hdataMCratio->Divide(hMC);

  hdataMCratio->Rebin(2);
  fout->Write();
  fout->Close();
  delete fout;

  //////////////////////////                                                                                                                                          


}
