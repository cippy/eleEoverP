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

#define FIT_2SIDE_CB 1      // 0 for single tail Crystal Ball for the fit, 1 for double tail 
#define USE_E 1 // 0 for ET and 1 for E in the binning                                       
#define USE_RAWE 0 // when 1, use raw SC energy instead of regression corrected ECAL energy  
#define READ_FROM_LOCAL 1
#define USE_P_AT_VTX 1
#define ELE_ETA_MAX 1.47

string finalDirOutName = "/afs/cern.ch/user/m/mciprian/www/EoverP/plot/2016_singleEleRunBtoG_Eregr_pAtVtx_fit2sideCB_weightReco_eta1p47/";

//=======================================================

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

//===============================================================================================

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

//=========================================================================

void plotDistribution(vector<TH1F*> hToFit, TH1F * hToFill, const string &dirName, const string &hNameID, const vector<Float_t> &energybinEdges, const string whatBin = "Ecorr", const string whatRatio = "wgtOverMf") {

  cout << endl;

  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2() 
  TVirtualFitter::SetDefaultFitter("Minuit");

  TCanvas *cdistr = NULL;

  UInt_t nBins = energybinEdges.size() -1;

  for (UInt_t i = 0; i < nBins; i++) {

    cdistr = new TCanvas("c",""); 

    // cout << "check -2" << endl;
    cout << " Energy bin: [" << energybinEdges[i] << ", " << energybinEdges[i+1] << "]" << endl;

    //    hToFit[i]->SetStats(0);
    // cout << "check -1" << endl;
    gStyle->SetOptStat(10);
    gStyle->SetOptFit(0112);
    // cout << "check :(" << endl;
    hToFit[i]->Draw("HE");
    if (whatRatio == "wgtOverMf") hToFit[i]->GetXaxis()->SetRangeUser(0.9,1.15);
    else hToFit[i]->GetXaxis()->SetRangeUser(0.85,1.1);
    //cout << "check o_O" << endl;
    // if (energybinEdges[i] > 449.9) hToFit[i]->Rebin(5);
    // else if (energybinEdges[i] > 349.9) hToFit[i]->Rebin(4);
    // else if (energybinEdges[i] > 274.9) hToFit[i]->Rebin(3);
    // else if (energybinEdges[i] > 174.9) hToFit[i]->Rebin(2);
    cdistr->Update();

    // cout << "check 0" << endl;

    // fitting:
    // do a first fit with a simple gaussian in the core
    Double_t gaussEdgeL = 0.99;  //left side of the gaussian to be used in the fit (I use a variable so that I change this value only once)
    Double_t gaussEdgeR = 1.03;  //right side ...
    if (whatRatio == "wgtOverMf") {
      gaussEdgeL = 0.99;  //left side of the gaussian to be used in the fit (I use a variable so that I change this value only once)
      gaussEdgeR = 1.03;  //right side ...
    } else {
      gaussEdgeL = 0.97;  //left side of the gaussian to be used in the fit (I use a variable so that I change this value only once)
      gaussEdgeR = 1.01;  //right side ...
    }
    // some general settings
    // if (energybinEdges[i] > 349.9) {
    //   gaussEdgeL = 0.75;
    //   gaussEdgeR = 1.25;
    // }
    // if MC, use more specific settings (decided after looking at plots, there is no a priori reason for them)

    hToFit[i]->Fit("gaus","E L I Q 0","",gaussEdgeL,gaussEdgeR);  // L: loglikelihood method, 0: do not plot this fit, Q: quiet mode (minimum printing)
    Double_t gaussNorm = hToFit[i]->GetFunction("gaus")->GetParameter(0);
    Double_t gaussMean = hToFit[i]->GetFunction("gaus")->GetParameter(1);
    //Double_t gaussMeanError = hToFit[i]->GetFunction("gaus")->GetParError(1);
    Double_t gaussSigma = hToFit[i]->GetFunction("gaus")->GetParameter(2);
    // now use crystal ball with right tail
    Double_t funcRangeLeft = gaussEdgeL;
    Double_t funcRangeRight = 1.13;
    if (whatRatio == "wgtOverMf") {
      funcRangeLeft = gaussEdgeL;
      funcRangeRight = 1.13;
    } else {
      funcRangeLeft = 0.87;
      funcRangeRight = gaussEdgeR;
    }

    if (FIT_2SIDE_CB) {
      if (whatRatio == "wgtOverMf") funcRangeLeft = 0.98;
      else funcRangeRight = 1.02;
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

      cb1 = new TF1("cb1",&myCrystalBallRightTail,funcRangeLeft,funcRangeRight,5);
      cb1->SetParNames("alpha","n","mu","sigma","N");  
      cb1->SetParLimits(cb1->GetParNumber("n"),0.1,15); 
      cb1->SetParLimits(cb1->GetParNumber("alpha"),0.01,10);
      cb1->SetParameters((gaussEdgeR-gaussMean)/gaussSigma,5,gaussMean,gaussSigma,gaussNorm);
      // with the following (or some of them) it looks like the fit doesn't work well, tipically the value from fit is out of the range
      // cb1->SetParLimits(cb1->GetParNumber("mu"),gaussMean-3.0*gaussSigma,gaussMean+3.0*gaussSigma);
      // cb1->SetParLimits(cb1->GetParNumber("sigma"),0.1*gaussSigma,10.0*gaussSigma);
      // cb1->SetParLimits(cb1->GetParNumber("N"),0.1*gaussNorm,10.0*gaussNorm);
    }

    // build TLegend
    TLegend *leg = new TLegend(0.11,0.8,0.4,0.89);
    TLegend *legFitFunction = new TLegend(0.11,0.5,0.35,0.8);

    // do the fit
    TFitResultPtr frp1 = hToFit[i]->Fit(cb1,"E L I S Q B R","HE",funcRangeLeft,funcRangeRight);
    hToFill->SetBinContent(i+1, frp1->Parameter(2)); // 2 is mu (starts with alpha, which is parameter number 0)       
    hToFill->SetBinError(i+1, frp1->ParError(2)); // 2 is mu (starts with alpha, which is parameter number 0)           

    leg->AddEntry(hToFit[i],"distribution","l");
    legFitFunction->SetHeader("fit functions:");
    legFitFunction->AddEntry(cb1,"Crystal Ball","l");
  
    leg->Draw();
    leg->SetMargin(0.3); 
    leg->SetBorderSize(0);

    legFitFunction->Draw();
    legFitFunction->SetMargin(0.3); 
    legFitFunction->SetBorderSize(0);

    cdistr->Update();
    // box for fit with Crystal Ball
    TPaveStats *stat = (TPaveStats*)(hToFit[i]->FindObject("stats"));
    if(stat) {
      // stat->SetTextColor(kBlue);
      // stat1->SetTextColor(kGreen);
      float width = stat->GetX2NDC() - stat->GetX1NDC();
      // make stat box bigger
      stat->SetX1NDC(stat->GetX1NDC() - 0.5 * width);
      if (FIT_2SIDE_CB) stat->SetY1NDC(stat->GetY2NDC() - 1.5 * (stat->GetY2NDC() - stat->GetY1NDC()));     
      stat->Draw();
    }

    // string whichEnergy = "";
    // if (USE_RAWE) whichEnergy = "raw E";
    // else whichEnergy = "corrected E";
    // if ((sampleName != "DATA") && (hNameID != "EoverP")) whichEnergy = "true E";

    if (whatRatio == "mfOverWgt") hToFit[i]->GetXaxis()->SetTitle("E(multifit)/E(weight)");
    else if (whatRatio == "wgtOverMf") hToFit[i]->GetXaxis()->SetTitle("E(weight)/E(multfit)");

    if (whatBin == "Ecorr") {
      hToFit[i]->SetTitle(Form("%1.0f < E_{corr} [GeV] < %1.0f",energybinEdges[i],energybinEdges[i+1]));
      cdistr->SaveAs(Form("%s%s_Ecorr%1.0fTo%1.0f.pdf",dirName.c_str(),hNameID.c_str(),energybinEdges[i],energybinEdges[i+1]));
      cdistr->SaveAs(Form("%s%s_Ecorr%1.0fTo%1.0f.png",dirName.c_str(),hNameID.c_str(),energybinEdges[i],energybinEdges[i+1]));
    } else if (whatBin == "Eraw") {
      hToFit[i]->SetTitle(Form("%1.0f < E_{raw} [GeV] < %1.0f",energybinEdges[i],energybinEdges[i+1]));
      cdistr->SaveAs(Form("%s%s_Eraw%1.0fTo%1.0f.pdf",dirName.c_str(),hNameID.c_str(),energybinEdges[i],energybinEdges[i+1]));
      cdistr->SaveAs(Form("%s%s_Eraw%1.0fTo%1.0f.png",dirName.c_str(),hNameID.c_str(),energybinEdges[i],energybinEdges[i+1]));
    } else if (whatBin == "Ptrack") {
      hToFit[i]->SetTitle(Form("%1.0f < P_{track} [GeV] < %1.0f",energybinEdges[i],energybinEdges[i+1]));
      cdistr->SaveAs(Form("%s%s_Ptrack%1.0fTo%1.0f.pdf",dirName.c_str(),hNameID.c_str(),energybinEdges[i],energybinEdges[i+1]));
      cdistr->SaveAs(Form("%s%s_Ptrack%1.0fTo%1.0f.png",dirName.c_str(),hNameID.c_str(),energybinEdges[i],energybinEdges[i+1]));
    }
    
    //cout << "CHECK 1" << endl;
    delete cdistr;
    delete leg;
    //cout << "CHECK 2" << endl;
    
    cout << endl;
    
  }

}



//=====================================================================                                                             

void plotGain(TH1F* h1Gain, TH1F* h2Gain, const string dirName) {

  TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2()                                        

  gStyle->SetOptStat(10);

  TCanvas *c = new TCanvas("c","",700,700);
  TLegend *leg = new TLegend(0.45,0.65,0.9,0.87);
  h1Gain->SetStats(0);
  h2Gain->SetStats(0);
 
  h1Gain->SetLineColor(kBlue);
  h2Gain->SetLineColor(kRed);

  h1Gain->Scale(1./h1Gain->Integral());
  h2Gain->Scale(1./h2Gain->Integral());

  h1Gain->Draw("HE");
  h2Gain->Draw("HE SAME");

  Double_t max = h1Gain->GetBinContent(h1Gain->GetMaximumBin());
  if (h2Gain->GetBinContent(h2Gain->GetMaximumBin()) > max) max = h2Gain->GetBinContent(h2Gain->GetMaximumBin());
  h1Gain->SetMaximum(1.4);
  h1Gain->SetMinimum(-0.05);

  TPaveText *pt = new TPaveText(0.15,0.62,0.4,0.87,"NDC");
  pt->AddText("0 : gain 12");
  pt->AddText("1 : gain  6");
  pt->AddText("2 : gain  1");
  pt->Draw();
  pt->SetMargin(0.3); 
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);

  h1Gain->SetTitle("450 < E_{corr} < 650 GeV");
  h1Gain->GetXaxis()->SetTitle("gain");
  h1Gain->GetYaxis()->SetTitle("a.u.");
  h1Gain->GetYaxis()->SetTitleSize(0.045);
  h1Gain->GetYaxis()->SetTitleOffset(1.2);

  leg->AddEntry(h1Gain,"E_{multifit}/E_{weight} < 0.95","l");  
  leg->AddEntry(h2Gain,"E_{multifit}/E_{weight} > 0.95","l");  
  leg->Draw();
  leg->SetMargin(0.3); 
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  c->Update();
  c->SaveAs((dirName+"gain_E450to650.pdf").c_str());
  c->SaveAs((dirName+"gain_E450to650.png").c_str());

  delete c;
  delete leg;
  
}

//=====================================================================                                                             

void plotEtWithGain(TH1F* g1, TH1F* g6, TH1F* g12, const string dirName, const string EorEt = "Et") {

  TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2()                                        

  gStyle->SetOptStat(10);

  TCanvas *c = new TCanvas("c","",700,700);
  TLegend *leg = new TLegend(0.6,0.59,0.9,0.89);
 
  g1->Scale(1./g1->Integral());
  g6->Scale(1./g6->Integral());
  g12->Scale(1./g12->Integral());

  Double_t max = g1->GetBinContent(g1->GetMaximumBin());
  if (g6->GetBinContent(g6->GetMaximumBin()) > max) max = g6->GetBinContent(g6->GetMaximumBin());
  if (g12->GetBinContent(g12->GetMaximumBin()) > max) max = g12->GetBinContent(g12->GetMaximumBin());
  g1->SetMaximum(max * 1.05);

  g1->SetLineColor(kRed);
  g6->SetLineColor(kBlue);
  g12->SetLineColor(kBlack);
  g1->Draw("HE");
  g6->Draw("HE SAME");
  g12->Draw("HE SAME");

  g1->SetStats(0);
  g6->SetStats(0);
  g12->SetStats(0);
 
  if (EorEt == "Et") g1->GetXaxis()->SetTitle("E_{T} [GeV]");
  else if (EorEt == "E") g1->GetXaxis()->SetTitle("E [GeV]");
  g1->GetYaxis()->SetTitle("a.u.");
  g1->GetYaxis()->SetTitleSize(0.045);
  g1->GetYaxis()->SetTitleOffset(1.2);

  leg->SetHeader(Form("EB: |#eta| < %2.2f",ELE_ETA_MAX));
  leg->AddEntry(g1,"gain1","l");  
  leg->AddEntry(g6,"gain 6","l");  
  leg->AddEntry(g12,"gain 12","l");  
  leg->Draw();
  leg->SetMargin(0.3); 
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  c->Update();
  if (EorEt == "Et") {
    c->SaveAs((dirName+"Et_gain.pdf").c_str());
    c->SaveAs((dirName+"Et_gain.png").c_str());
  } else {
    c->SaveAs((dirName+"E_gain.pdf").c_str());
    c->SaveAs((dirName+"E_gain.png").c_str());
  }

  delete c;
  delete leg;

}


//=====================================================================                                                             

void compareWeightAndMultifit() {

  TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2()                                        

  TChain* twgt = new TChain("selected");
  TChain* tmf = new TChain("selected");
  buildChain(tmf,  "DATA"  );
  buildChain(twgt, "DATA_w");

  // just need to build index for the other tree (tmf) that will be attached as friend to twgt 
  //twgt->BuildIndex("runNumber","eventNumber");
  tmf->BuildIndex("runNumber","eventNumber");  

  Int_t           runNumber_wgt, runNumber_mf ;
  ULong64_t       eventNumber_wgt, eventNumber_mf;
  Float_t         PtEle[3];
  Float_t         etaEle[3];
  Float_t         R9Ele[3];
  Float_t         energySCEle_must_regrCorr_ele[3];
  Float_t         rawEnergySCEle_must[3];
  Float_t         pGsfEle[3];

  Float_t         PtEle_mf[3];
  Float_t         etaEle_mf[3];
  Float_t         R9Ele_mf[3];
  Float_t         energySCEle_must_regrCorr_ele_mf[3];
  Float_t         rawEnergySCEle_must_mf[3];
  Float_t         pGsfEle_mf[3];
  UChar_t         gainEle_mf[3];

  twgt->SetBranchAddress("runNumber",&runNumber_wgt);
  twgt->SetBranchAddress("eventNumber",&eventNumber_wgt);
  twgt->SetBranchAddress("PtEle",PtEle);
  twgt->SetBranchAddress("etaEle",etaEle);
  twgt->SetBranchAddress("R9Ele",R9Ele);
  twgt->SetBranchAddress("energySCEle_must_regrCorr_ele",energySCEle_must_regrCorr_ele);
  twgt->SetBranchAddress("rawEnergySCEle_must",rawEnergySCEle_must);
  if (USE_P_AT_VTX) twgt->SetBranchAddress("pAtVtxGsfEle",pGsfEle);
  else twgt->SetBranchAddress("pModeGsfEle",pGsfEle);

  tmf->SetBranchAddress("runNumber",&runNumber_mf);
  tmf->SetBranchAddress("eventNumber",&eventNumber_mf);
  tmf->SetBranchAddress("PtEle",PtEle_mf);
  tmf->SetBranchAddress("etaEle",etaEle_mf);
  tmf->SetBranchAddress("R9Ele",R9Ele_mf);
  tmf->SetBranchAddress("energySCEle_must_regrCorr_ele",energySCEle_must_regrCorr_ele_mf);
  tmf->SetBranchAddress("rawEnergySCEle_must",rawEnergySCEle_must_mf);
  if (USE_P_AT_VTX) tmf->SetBranchAddress("pAtVtxGsfEle",pGsfEle_mf);
  else tmf->SetBranchAddress("pModeGsfEle",pGsfEle_mf);
  tmf->SetBranchAddress("gainEle",gainEle_mf);

  twgt->AddFriend(tmf);

  vector<Float_t> energybinEdges;
  // this binning looks ok                    
  energybinEdges.push_back(25.0);
  energybinEdges.push_back(50.0);
  energybinEdges.push_back(75.0);
  energybinEdges.push_back(100.0);
  energybinEdges.push_back(125.0);
  energybinEdges.push_back(175.0);
  energybinEdges.push_back(225.0);
  energybinEdges.push_back(275.0);
  energybinEdges.push_back(350.0);
  energybinEdges.push_back(450.0);
  energybinEdges.push_back(650.0);
  energybinEdges.push_back(900.0);
  
  Int_t nEnergyBins = energybinEdges.size() - 1;

  TH1F* hModeEcorrMfOverWgtInBin = new TH1F("hModeEcorrMfOverWgtInBin","",nEnergyBins,energybinEdges.data());
  TH1F* hEcorrMfOverWgt = new TH1F("hEcorrMfOverWgt","",80,0.81,1.21);
  vector<TH1F*> hEcorrMfOverWgt_pTrackBin(nEnergyBins,NULL);
  for (Int_t i = 0; i < nEnergyBins; i++) {
    hEcorrMfOverWgt_pTrackBin[i] = new TH1F(Form("hEcorrMfOverWgt_pTrackBin%1.0fTo%1.0f",energybinEdges[i],energybinEdges[i+1]),"",80,0.81,1.21);
  }

  TH1F* hModeEMfOverWgtInBin = new TH1F("hModeEMfOverWgtInBin","",nEnergyBins,energybinEdges.data());
  TH1F* hEMfOverWgt = new TH1F("hEMfOverWgt","",80,0.81,1.21);
  vector<TH1F*> hEMfOverWgt_EBin(nEnergyBins,NULL);
  for (Int_t i = 0; i < nEnergyBins; i++) {
    hEMfOverWgt_EBin[i] = new TH1F(Form("hEMfOverWgt_EBin%1.0fTo%1.0f",energybinEdges[i],energybinEdges[i+1]),"",80,0.81,1.21);
  }

  TH1F *h1Gain = new TH1F("h1Gain","",3,-0.5,2.5);
  TH1F *h2Gain = new TH1F("h2Gain","",3,-0.5,2.5);

  TH1F *hEtGain1 = new TH1F("hEtGain1","",90,0,900);
  TH1F *hEtGain6 = new TH1F("hEtGain6","",90,0,900);
  TH1F *hEtGain12 = new TH1F("hEtGain12","",90,0,900);

  TH1F *hEGain1 = new TH1F("hEGain1","",90,0,900);
  TH1F *hEGain6 = new TH1F("hEGain6","",90,0,900);
  TH1F *hEGain12 = new TH1F("hEGain12","",90,0,900);

  Long64_t count = 0;
  Long64_t countCut = 0;
  Double_t EwgtOverEmf = 0.0;
  Double_t EmfOverEwgt = 0.0;
  Double_t energyToUse = 0.0;
  Double_t energyToUse_mf = 0.0;

  // tree with weight has more events
  Long64_t nentries = twgt->GetEntries();

  for (Long64_t i=0;i<nentries;i++) {
    twgt->GetEntry(i);

    if (i%500000 == 0) cout << i << endl;

    // select only events in both trees
    if ( !(runNumber_wgt == runNumber_mf && eventNumber_wgt == eventNumber_mf)) continue;
    //else if (i < 100) cout << runNumber_wgt << "   " << runNumber_mf << "\t" << eventNumber_wgt << "   " << eventNumber_mf << endl;
    count++;
    if (PtEle[0] < 40 || fabs(etaEle[0]) > ELE_ETA_MAX || R9Ele[0] < 0.94) continue;
    if (PtEle_mf[0] < 40 || fabs(etaEle_mf[0]) > ELE_ETA_MAX || R9Ele_mf[0] < 0.94) continue;
    countCut++;

    // if (i%50000 == 0) cout << energySCEle_must_regrCorr_ele[0]/energySCEle_must_regrCorr_ele_mf[0] << endl;
    if (USE_RAWE) {
      energyToUse = rawEnergySCEle_must[0];
      energyToUse_mf = rawEnergySCEle_must_mf[0];
    } else {
      energyToUse = energySCEle_must_regrCorr_ele[0];
      energyToUse_mf = energySCEle_must_regrCorr_ele_mf[0];
    }
    EwgtOverEmf = energyToUse/energyToUse_mf;
    EmfOverEwgt = energyToUse_mf/energyToUse;

    hEcorrMfOverWgt->Fill(EmfOverEwgt);
    hEMfOverWgt->Fill(EmfOverEwgt);

    Int_t binP = getBinNumber(pGsfEle[0],energybinEdges);  // this function returns negative value if bin not found       
    if (binP >= 0) {
      hEcorrMfOverWgt_pTrackBin[binP]->Fill(EmfOverEwgt);
    } else if (binP == -1) {
      // fill last bin with overflows to gain in statistics                                                                                                           
      hEcorrMfOverWgt_pTrackBin[nEnergyBins-1]->Fill(EmfOverEwgt);
    }

    Int_t binE;
    binE = getBinNumber(energyToUse_mf,energybinEdges);  // this function returns negative value if bin not found       
   
    if (binE >= 0) {
      hEMfOverWgt_EBin[binE]->Fill(EmfOverEwgt);
    } else if (binE == -1) {
      // fill last bin with overflows to gain in statistics                                                                                                           
      hEMfOverWgt_EBin[nEnergyBins-1]->Fill(EmfOverEwgt);
    }

    if (energyToUse_mf >= 450 && energyToUse_mf <= 650) {

      if (EmfOverEwgt < 0.95) {
	h1Gain->Fill(gainEle_mf[0]);
      } else {
	h2Gain->Fill(gainEle_mf[0]);
      }
      
    }

    Double_t Et = energyToUse_mf / cosh(etaEle_mf[0]);
    if (gainEle_mf[0] == 0) {
      hEtGain12->Fill(Et); 
      hEGain12->Fill(energyToUse_mf);
    } else if (gainEle_mf[0] == 1) {
      hEtGain6->Fill(Et);
      hEGain6->Fill(energyToUse_mf);
    } else if (gainEle_mf[0] == 2) {
      hEtGain1->Fill(Et);
      hEGain1->Fill(energyToUse_mf);
    }

  }   // loop end
  ////////////////////////////
  ////////////////////////////

  //void plotDistribution(vector<TH1F*> &hToFit, TH1F * hToFill, const string &dirName, const string &hNameID, const vector<Float_t> &energybinEdges) {
  // for (UInt_t i = 0; i < hEcorrMfOverWgt_pTrackBin.size(); i++) {
  //   cout << hEcorrMfOverWgt_pTrackBin[i]->GetEntries() << endl;
  // }

  if (USE_RAWE) {

    plotDistribution(hEcorrMfOverWgt_pTrackBin, 
		     hModeEcorrMfOverWgtInBin, 
		     finalDirOutName,
		     "ErawMfOverWgt",
		     energybinEdges,
		     "Ptrack",
		     "mfOverWgt");

    plotDistribution(hEMfOverWgt_EBin, 
		     hModeEMfOverWgtInBin, 
		     finalDirOutName,
		     "ErawMfOverWgt",
		     energybinEdges,
		     "Eraw",
		     "mfOverWgt");

  } else {

    plotDistribution(hEcorrMfOverWgt_pTrackBin, 
		     hModeEcorrMfOverWgtInBin, 
		     finalDirOutName,
		     "EcorrMfOverWgt",
		     energybinEdges,
		     "Ptrack",
		     "mfOverWgt");

    plotDistribution(hEMfOverWgt_EBin, 
		     hModeEMfOverWgtInBin, 
		     finalDirOutName,
		     "EcorrMfOverWgt",
		     energybinEdges,
		     "Ecorr",
		     "mfOverWgt");

  }

  TCanvas *c1 = new TCanvas("c1","",700,700);
  hEcorrMfOverWgt->SetStats(0);
  hEcorrMfOverWgt->SetLineColor(kBlue);
  hEcorrMfOverWgt->Draw("HE");
  if (USE_RAWE) hEcorrMfOverWgt->GetXaxis()->SetTitle("E_{raw}(multifit) / E_{raw}(weight)");
  else hEcorrMfOverWgt->GetXaxis()->SetTitle(" E_{corr}(multifit) / E_{corr}(weight)");
  //  hEcorrMfOverWgt->GetXaxis()->SetTitleOffset(0.8);
  //  hEcorrMfOverWgt->GetXaxis()->SetTitleSize(0.04);
  hEcorrMfOverWgt->GetYaxis()->SetTitle("Events");
  hEcorrMfOverWgt->GetYaxis()->SetTitleSize(0.04);
  hEcorrMfOverWgt->GetYaxis()->SetTitleOffset(1.2);
  hEcorrMfOverWgt->GetYaxis()->CenterTitle();
  if (USE_RAWE) {
    c1->SaveAs((finalDirOutName + "ErawMfOverWgt_pBin.pdf").c_str());
    c1->SaveAs((finalDirOutName + "ErawMfOverWgt_pBin.png").c_str());
  } else {
    c1->SaveAs((finalDirOutName + "EcorrMfOverWgt_pBin.pdf").c_str());
    c1->SaveAs((finalDirOutName + "EcorrMfOverWgt_pBin.png").c_str());
  }
  delete c1;

  TCanvas *c2 = new TCanvas("c2","",700,700);
  hModeEcorrMfOverWgtInBin->SetStats(0);
  hModeEcorrMfOverWgtInBin->SetLineColor(kBlue);
  hModeEcorrMfOverWgtInBin->Draw("HE");
  hModeEcorrMfOverWgtInBin->GetXaxis()->SetTitle("p_{track} [GeV]");
  //  hModeEcorrMfOverWgtInBin->GetXaxis()->SetTitleOffset(0.8);
  //  hModeEcorrMfOverWgtInBin->GetXaxis()->SetTitleSize(0.04);
  hModeEcorrMfOverWgtInBin->GetYaxis()->SetTitle("Peak( E(multifit) / E(weight) )");
  if (FIT_2SIDE_CB) hModeEcorrMfOverWgtInBin->SetTitle("fit with double Crystal Ball");
  else hModeEcorrMfOverWgtInBin->SetTitle("fit with Crystal Ball");
  hModeEcorrMfOverWgtInBin->GetYaxis()->SetTitleSize(0.04);
  hModeEcorrMfOverWgtInBin->GetYaxis()->SetTitleOffset(1.2);
  hModeEcorrMfOverWgtInBin->GetYaxis()->CenterTitle();
  if (USE_RAWE) {
    c2->SaveAs((finalDirOutName + "peakErawMfOverWgt_pBin.pdf").c_str());
    c2->SaveAs((finalDirOutName + "peakErawMfOverWgt_pBin.png").c_str());
  } else {
    c2->SaveAs((finalDirOutName + "peakEcorrMfOverWgt_pBin.pdf").c_str());
    c2->SaveAs((finalDirOutName + "peakEcorrMfOverWgt_pBin.png").c_str());
  }
  delete c2;

  //===============================

  TCanvas *c3 = new TCanvas("c3","",700,700);
  hEMfOverWgt->SetStats(0);
  hEMfOverWgt->SetLineColor(kBlue);
  hEMfOverWgt->Draw("HE");
  if (USE_RAWE) hEMfOverWgt->GetXaxis()->SetTitle("E_{raw}(multifit) / E_{raw}(weight)");
  else hEMfOverWgt->GetXaxis()->SetTitle("E_{corr}(multifit) / E_{corr}(weight)");
  //  hEMfOverWgt->GetXaxis()->SetTitleOffset(0.8);
  //  hEMfOverWgt->GetXaxis()->SetTitleSize(0.04);
  hEMfOverWgt->GetYaxis()->SetTitle("Events");
  hEMfOverWgt->GetYaxis()->SetTitleSize(0.04);
  hEMfOverWgt->GetYaxis()->SetTitleOffset(1.2);
  hEMfOverWgt->GetYaxis()->CenterTitle();
  if (USE_RAWE) {
    c3->SaveAs((finalDirOutName + "ErawMfOverWgt.pdf").c_str());
    c3->SaveAs((finalDirOutName + "ErawMfOverWgt.png").c_str());
  } else {
    c3->SaveAs((finalDirOutName + "EcorrMfOverWgt.pdf").c_str());
    c3->SaveAs((finalDirOutName + "EcorrMfOverWgt.png").c_str());
  }
  delete c3;

  TCanvas *c4 = new TCanvas("c4","",700,700);
  hModeEMfOverWgtInBin->SetStats(0);
  hModeEMfOverWgtInBin->SetLineColor(kBlue);
  hModeEMfOverWgtInBin->Draw("HE");
  if (USE_RAWE) hModeEMfOverWgtInBin->GetXaxis()->SetTitle("E_{raw} [GeV]");
  else hModeEMfOverWgtInBin->GetXaxis()->SetTitle("E_{corr} [GeV]");
  //  hModeEMfOverWgtInBin->GetXaxis()->SetTitleOffset(0.8);
  //  hModeEMfOverWgtInBin->GetXaxis()->SetTitleSize(0.04);
  hModeEMfOverWgtInBin->GetYaxis()->SetTitle("Peak( E(multifit) / E(weight) )");
  if (FIT_2SIDE_CB) hModeEMfOverWgtInBin->SetTitle("fit with double Crystal Ball");
  else hModeEMfOverWgtInBin->SetTitle("fit with Crystal Ball");
  hModeEMfOverWgtInBin->GetYaxis()->SetTitleSize(0.04);
  hModeEMfOverWgtInBin->GetYaxis()->SetTitleOffset(1.2);
  hModeEMfOverWgtInBin->GetYaxis()->CenterTitle();
  if (USE_RAWE) {
    c4->SaveAs((finalDirOutName + "peakErawMfOverWgt.pdf").c_str());
    c4->SaveAs((finalDirOutName + "peakErawMfOverWgt.png").c_str());
  } else {
    c4->SaveAs((finalDirOutName + "peakEcorrMfOverWgt.pdf").c_str());
    c4->SaveAs((finalDirOutName + "peakEcorrMfOverWgt.png").c_str());
  }
  delete c4;

  plotGain(h1Gain,h2Gain,finalDirOutName);
  plotEtWithGain(hEtGain1,hEtGain6,hEtGain12,finalDirOutName);
  plotEtWithGain(hEGain1,hEGain6,hEGain12,finalDirOutName,"E");

  cout << "nentries = " << nentries << "\tcount = " << count << "\tcountCut = " << countCut << endl;

}
