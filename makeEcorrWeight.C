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
  // hData->Scale(1./hData->Integral());
  // hMC->Scale(1./hMC->Integral());
  // hData->Sumw2();
  // hMC->Sumw2();
  Double_t error = 0.0;  //dummy variable to use YH1::IntegralAndError()
  hData->Scale(1./hData->IntegralAndError(0,hData->GetNbinsX()+1,error));
  hMC->Scale(1./hMC->IntegralAndError(0,hMC->GetNbinsX()+1,error));

  //////////////////////////////////////////////
  // draw distribution

  TCanvas *c = new TCanvas("c","",700,700);

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

  maximumYaxisValue = TMath::Max(hData->GetBinContent(hData->GetMaximumBin()),hMC->GetBinContent(hMC->GetMaximumBin()));
  minimumYaxisValue = TMath::Min(hData->GetBinContent(hData->GetMinimumBin()),hMC->GetBinContent(hMC->GetMinimumBin()));
  
  Double_t diff = maximumYaxisValue - minimumYaxisValue;
  minimumYaxisValue -= diff * 0.1;
  maximumYaxisValue += diff * 0.1;

  if (minimumYaxisValue < 0.0000001) minimumYaxisValue = 0.00001; // if compatible with 0 or negative, set it to value close to but different from zero
  subpad_1->SetLogy(); 

  hData->SetStats(0);
  hData->SetLineColor(kRed);
  hData->Draw("HE");
  hData->SetTitle("distributions before reweighting");
  hData->GetXaxis()->SetLabelSize(0.45);
  hData->GetYaxis()->SetTitle("a.u.");
  hData->GetYaxis()->SetTitleSize(0.06);
  hData->GetYaxis()->SetTitleOffset(0.8);
  hData->SetMinimum(minimumYaxisValue);
  hData->SetMaximum(maximumYaxisValue);
  hData->SetMinimum(0.00001);
  hData->SetMaximum(1.0);
  hMC->SetLineColor(kBlue);
  hMC->Draw("HE SAME");

  string texMCSampleName = "";
  if (MCSampleName == "WJetsToLNu") texMCSampleName = "W(l#nu)+jets";
  else if (MCSampleName == "DYJetsToLL") texMCSampleName = "Z(ll)+jets";

  leg->AddEntry(hData,"data","lf");
  leg->AddEntry(hMC,Form("%s",texMCSampleName.c_str()),"lf");
  leg->Draw();
  leg->SetMargin(0.3);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);  // transparent legend

  //ratio plot

  subpad_2->cd();
  TH1F * ratioplot = NULL; // will use it for the ratio plots                                                                                                   
  TH1F *tmp = (TH1F*) hData->Clone("cloneOfHdata");
  ratioplot = new TH1F(*tmp);  // to have ratioplot with the same x axis range as hData (or hMC), I copy it from hData created above, then I substitute bin content with hData/hMC                                                                                                                                                   
  ratioplot->Divide(hData,hMC);
  ratioplot->SetStats(0);
  ratioplot->SetTitle("");
  ratioplot->GetXaxis()->SetLabelSize(0.10);
  ratioplot->GetXaxis()->SetTitle("electron corrected energy [GeV]");
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
  ratioplotCopy->SetMinimum(0.0);
  ratioplotCopy->SetMaximum(2.0);
  
  c->SaveAs( (dirName + "EcorrBeforeReweighting.pdf").c_str() );
  c->SaveAs( (dirName + "EcorrBeforeReweighting.png").c_str() );

  delete c;

  ///////////////////////////////////////////////

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
  hdataMCratio->SetTitle("data/MC ratio to reweight MC");

  hdataMCratio->Rebin(2);
  fout->Write();
  fout->Close();
  delete fout;

  //////////////////////////                                                                                                                                          


}
