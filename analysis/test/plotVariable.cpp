#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <iomanip>

#include "TMath.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h" 
#include "TGraphErrors.h" 
#include "TStyle.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"


#include "TDRStyle.h"
#include "utilities.hh"

int main( int argc, char** argv )
{//main

  std::string varName = "dijet_M";
  //std::string path = "/data/shared/Long_Exercise_Hinvisible/Skims/";
  std::string path = "../../data/Skims8TeV/";

  double lumiData = 20000;//in pb-1

  double xmin  = 1000;
  double xmax = 4000;
  double nBins = 100;
  bool dology = true;

  std::string title = ";M_{jj} (GeV); events/30 GeV";

  const unsigned nCat = 4;
  std::string bkgcat[nCat] = {"Top","VV","W+jets","Z+jets"};

  double dataSF[nCat] = {0.84,1,0.7,0.67*0.6/0.0334};

  std::vector<std::string> files;
  files.push_back("MC_TTJets");
  files.push_back("MC_T-tW");
  files.push_back("MC_Tbar-tW");
  files.push_back("MC_SingleT-s-powheg-tauola");
  files.push_back("MC_SingleTBar-s-powheg-tauola");
  files.push_back("MC_SingleT-t-powheg-tauola");
  files.push_back("MC_SingleTBar-t-powheg-tauola");
  files.push_back("MC_WW-pythia6-tauola");
  files.push_back("MC_WZ-pythia6-tauola");
  files.push_back("MC_ZZ-pythia6-tauola");
  files.push_back("MC_W1JetsToLNu_enu");
  files.push_back("MC_W2JetsToLNu_enu");
  files.push_back("MC_W3JetsToLNu_enu");
  files.push_back("MC_W4JetsToLNu_enu");
  files.push_back("MC_WJetsToLNu-v1_enu");
  files.push_back("MC_WJetsToLNu-v2_enu");
  files.push_back("MC_W1JetsToLNu_munu");
  files.push_back("MC_W2JetsToLNu_munu");
  files.push_back("MC_W3JetsToLNu_munu");
  files.push_back("MC_W4JetsToLNu_munu");
  files.push_back("MC_WJetsToLNu-v1_munu");
  files.push_back("MC_WJetsToLNu-v2_munu");
  files.push_back("MC_W1JetsToLNu_taunu");
  files.push_back("MC_W2JetsToLNu_taunu");
  files.push_back("MC_W3JetsToLNu_taunu");
  files.push_back("MC_W4JetsToLNu_taunu");
  files.push_back("MC_WJetsToLNu-v1_taunu");
  files.push_back("MC_WJetsToLNu-v2_taunu");
  files.push_back("MC_EWK-W2jminus_enu");
  files.push_back("MC_EWK-W2jplus_enu");
  files.push_back("MC_EWK-W2jminus_munu");
  files.push_back("MC_EWK-W2jplus_munu");
  files.push_back("MC_EWK-W2jminus_taunu");
  files.push_back("MC_EWK-W2jplus_taunu");
  //use DY ignoring muons for Znunu
  files.push_back("MC_DYJetsToLL");
  //files.push_back("MC_DYJetsToLL_PtZ-100-madgraph");
  files.push_back("MC_DY1JetsToLL");
  files.push_back("MC_DY2JetsToLL");
  files.push_back("MC_DY3JetsToLL");
  files.push_back("MC_DY4JetsToLL");
  //files.push_back("MC_ZJetsToNuNu_100_HT_200");
  //files.push_back("MC_ZJetsToNuNu_200_HT_400");
  //files.push_back("MC_ZJetsToNuNu_400_HT_inf");
  //files.push_back("MC_ZJetsToNuNu_50_HT_100");
  //files.push_back("MC_EWK-Z2j");
  files.push_back("MC_EWK-Z2j_iglep");

  const unsigned nMC = files.size();

  TFile *fMC[nMC];

  int color[nMC];
  int style[nMC];
  std::string label[nMC];


  for (unsigned iF(0); iF<nMC;++iF){//loop on files
    if (iF<7) {color[iF] = 5; style[iF] = 1001; label[iF] = bkgcat[0];}
    else if (iF<10) {color[iF] = 4; style[iF] = 3004; label[iF] = bkgcat[1];}
    else if (iF<34) {color[iF] = 2; style[iF] = 3005; label[iF] = bkgcat[2];}
    else {color[iF] = 3; style[iF] = 3006; label[iF] = bkgcat[3];}
  };

  TCanvas *myc = new TCanvas("myc","myc",1);

  TH1F *hist[nMC];
  TH1F *bkg[4]={0,0,0,0};

  for (unsigned iF(0); iF<nMC;++iF){//loop on files

    if (!testFile(path+files[iF]+".root",fMC[iF])) return 1;
    
    fMC[iF]->cd();
    TTree *tree = (TTree*)gDirectory->Get("LightTree");
    if (!tree) return 1;

    bool isDY = files[iF].find("MC_DY")!=files[iF].npos;
    //double weight = getNormalisationFactor(lumiData,files[iF]);

    std::ostringstream lselection,lname;
    std::string signal;
    if (!isDY) signal = "(nvetomuons==0 && nvetoelectrons==0)";
    else signal = "(nselmuons==2)";//dummy sel
    //lselection << signal << "*" << weight;
    //for 8TeV trees
    lselection << signal ;
    if (!isDY) lselection << "*total_weight_lepveto";
    else lselection << "*total_weight_leptight";
    lname << "hist" << iF;
    hist[iF] = new TH1F(lname.str().c_str(),title.c_str(),nBins,xmin,xmax);
    hist[iF]->Sumw2();
    lname.str("");
    lname << varName << ">>hist" << iF;
    tree->Draw(lname.str().c_str(),lselection.str().c_str(),"");

    hist[iF]->SetLineColor(color[iF]);
    hist[iF]->SetFillColor(color[iF]);
    hist[iF]->SetFillStyle(style[iF]);

    for (unsigned iC(0);iC<nCat;++iC){
      if (label[iF] != bkgcat[iC]) continue;
      if (!bkg[iC]) bkg[iC] = (TH1F*)hist[iF]->Clone(bkgcat[iC].c_str());
      else bkg[iC]->Add(hist[iF]);
    }

  }//loop on files

  gStyle->SetOptStat(0);

  myc->cd();
  gPad->SetLogy(dology);
  TLegend *leg = new TLegend(0.7,0.7,0.94,0.94);
  leg->SetFillColor(10);
  //myc->Divide(2,2);
  
  std::cout << " -- Yields: " << std::endl;
  for (unsigned iC(0);iC<nCat;++iC){
    bkg[iC]->Scale(dataSF[iC]);
    std::cout <<  bkgcat[iC] << " " << bkg[iC]->Integral() << std::endl;
  }

  //stack
  bkg[3]->Add(bkg[0]);
  bkg[3]->Add(bkg[1]);
  bkg[3]->Add(bkg[2]);
  bkg[2]->Add(bkg[0]);
  bkg[2]->Add(bkg[1]);
  bkg[1]->Add(bkg[0]);

  for (unsigned iC(0);iC<nCat;++iC){
    unsigned newid = nCat-iC-1;
    bkg[newid]->Draw(iC==0? "hist" : "histsame");
    leg->AddEntry(bkg[newid],bkgcat[newid].c_str(),"F");
  }
  leg->Draw("same");
  myc->Update();
  std::ostringstream lsave;
  lsave.str("");
  lsave << "PLOTS/" << varName;
  if (dology) lsave << "_log";
  lsave << ".pdf";
  myc->Print(lsave.str().c_str());

  return 0;

}//main
