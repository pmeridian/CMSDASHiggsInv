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
  std::string path = "/data/shared/meridian/Long_Exercise_Hinvisible/Skims/";
  //std::string path = "../../data/Skims8TeV/";

  double lumiData = 10000;//in pb-1

  double xmin  = 1000;
  double xmax = 4000;
  double nBins = 100;
  bool dology = true;

  std::string title = ";M_{jj} (GeV); events/30 GeV";

  const unsigned nCat = 4;
  std::string bkgcat[nCat] = {"Top","VV","W+jets","Z+jets"};

  double dataSF[nCat] = {0.84,1,0.7,0.67};//*0.6/0.0334};

  std::vector<std::string> files;
  //8 TeV
  /*
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
  */
  files.push_back("TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola");
  //files.push_back("TT_Tune4C_13TeV-pythia8-tauola");
  files.push_back("T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola");
  files.push_back("Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola");
  files.push_back("ZZTo4L_Tune4C_13TeV-powheg-pythia8");
  files.push_back("WJetsToLNu_13TeV-madgraph-pythia8-tauola");
  files.push_back("WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauol");
  files.push_back("WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola");
  files.push_back("WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola");
  files.push_back("WJetsToLNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola");
  files.push_back("ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola");
  files.push_back("ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola");
  files.push_back("ZJetsToNuNu_HT-400to600_Tune4C_13TeV-madgraph-tauola");
  files.push_back("ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola");
  //files.push_back("DYJetsToLL_M-50_13TeV-madgraph-pythia8");
  //files.push_back("DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola");
  //files.push_back("DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola");
  //files.push_back("DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola");
  //files.push_back("DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola");

  std::string sigfile = "VBF_HToInv_M-125_13TeV_powheg-pythia6";

  const unsigned nMC = files.size();

  TFile *fMC[nMC];

  int color[nMC+1];
  int style[nMC+1];
  std::string label[nMC];


  for (unsigned iF(0); iF<nMC;++iF){//loop on files
    if (iF<3) {color[iF] = 5; style[iF] = 1001; label[iF] = bkgcat[0];}
    else if (iF<4) {color[iF] = 4; style[iF] = 3004; label[iF] = bkgcat[1];}
    else if (iF<9) {color[iF] = 2; style[iF] = 3005; label[iF] = bkgcat[2];}
    else {color[iF] = 3; style[iF] = 3006; label[iF] = bkgcat[3];}
  };

  color[nMC] = 7;
  style[nMC] = 1001;

  TCanvas *myc = new TCanvas("myc","myc",1);

  TH1F *hist[nMC+1];
  TH1F *bkg[4]={0,0,0,0};

  TFile *fsig = 0;
  if (!testFile(path+sigfile+".root",fsig)) return 1;
  
  for (unsigned iF(0); iF<nMC+1;++iF){//loop on files
    if (iF==nMC) fsig->cd();
    else {
      if (!testFile(path+files[iF]+".root",fMC[iF])) return 1;
    
      fMC[iF]->cd();
    }
    TTree *tree = (TTree*)gDirectory->Get("LightTree");
    if (!tree) return 1;

    bool isDY = iF<nMC && files[iF].find("DY")!=files[iF].npos;
    
    double weight = getNormalisationFactor(lumiData,iF<nMC?files[iF]:sigfile);

    std::ostringstream lselection,lname;
    std::string signal;
    if (!isDY) signal = "(l1met>=70 && nvetomuons==0 && nvetoelectrons==0)";
    else signal = "(l1met>=70 && nselmuons==2)";//dummy sel
    lselection << signal << "*" << weight;
    //for 8TeV trees
    //lselection << signal ;
    //if (!isDY) lselection << "*total_weight_lepveto";
    //else lselection << "*total_weight_leptight";
    lname << "hist" << iF;
    hist[iF] = new TH1F(lname.str().c_str(),title.c_str(),nBins,xmin,xmax);
    hist[iF]->Sumw2();
    lname.str("");
    lname << varName << ">>hist" << iF;
    tree->Draw(lname.str().c_str(),lselection.str().c_str(),"");

    hist[iF]->SetLineColor(color[iF]);
    hist[iF]->SetFillColor(color[iF]);
    hist[iF]->SetFillStyle(style[iF]);
    if (iF<nMC){
      for (unsigned iC(0);iC<nCat;++iC){
	if (label[iF] != bkgcat[iC]) continue;
	if (!bkg[iC]) bkg[iC] = (TH1F*)hist[iF]->Clone(bkgcat[iC].c_str());
	else bkg[iC]->Add(hist[iF]);
      }
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
  std::cout << " signal " << hist[nMC]->Integral() << std::endl;  

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
  hist[nMC]->Draw("same");
  leg->AddEntry(hist[nMC],"VBF H, m_{H}=125 GeV","F");
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
