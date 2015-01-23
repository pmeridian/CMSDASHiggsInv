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
#include "TLatex.h"

#include "TDRStyle.h"
#include "utilities.hh"

enum SELECTION{
  SIGNAL=1,
  MuMu=2,
  LNu=3,
  Top=4
};

int main( int argc, char** argv )
{//main

  if (argc<8) {
    std::cout << " Usage: " << argv[0] << " "
	      << "<variable> <hist title> "
	      << "<nbin> <binmin> <binmax> "
	      << "<selection: 1=signal, 2=mumu, 3=lnu, 4=top>"
	      << "<dology> "
	      << std::endl;
    return 1;
  }

  std::string varName = argv[1];
  std::string title = argv[2];
  double nBins = atoi(argv[3]);
  double xmin  = atof(argv[4]);
  double xmax = atof(argv[5]);
  SELECTION sel = SELECTION(atoi(argv[6]));
  bool dology = atoi(argv[7]);;

  //std::string path = "/data/shared/meridian/Long_Exercise_Hinvisible/Skims/";
  std::string path = "./skims/";


  double lumiData = 10000;//in pb-1


  const unsigned nCat = 4;
  std::string bkgcat[nCat] = {"Top","VV","W+jets","Z+jets"};

  double dataSF[nCat] = {0.84,1,0.7,0.67};//*0.2/0.0334};

  std::vector<std::string> files;
  files.push_back("TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola");
  //files.push_back("TT_Tune4C_13TeV-pythia8-tauola");
  files.push_back("T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola");
  files.push_back("Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola");
  files.push_back("ZZTo4L_Tune4C_13TeV-powheg-pythia8");
  //files.push_back("WJetsToLNu_13TeV-madgraph-pythia8-tauola");
  files.push_back("WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauol");
  files.push_back("WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola");
  files.push_back("WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola");
  files.push_back("WJetsToLNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola");
  if (sel==SELECTION::SIGNAL){
    files.push_back("ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola");
    files.push_back("ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola");
    files.push_back("ZJetsToNuNu_HT-400to600_Tune4C_13TeV-madgraph-tauola");
    files.push_back("ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola");
  }
  else {
    files.push_back("DYJetsToLL_M-50_13TeV-madgraph-pythia8");
    files.push_back("DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola");
    files.push_back("DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola");
    files.push_back("DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola");
    files.push_back("DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola");
  }

  std::string sigfile = "VBF_HToInv_M-125_13TeV_powheg-pythia6";

  const unsigned nMC = files.size();

  TFile *fMC[nMC];

  int color[nMC+1];
  int style[nMC+1];
  std::string label[nMC];


  for (unsigned iF(0); iF<nMC;++iF){//loop on files
    if (iF<3) {color[iF] = 5; style[iF] = 1001; label[iF] = bkgcat[0];}
    else if (iF<4) {color[iF] = 3; style[iF] = 3004; label[iF] = bkgcat[1];}
    else if (iF<8) {color[iF] = 2; style[iF] = 3005; label[iF] = bkgcat[2];}
    else {color[iF] = 7; style[iF] = 3006; label[iF] = bkgcat[3];}
  };

  color[nMC] = 4;
  style[nMC] = 1001;

  TCanvas *myc = new TCanvas("myc","myc",1);
  SetTdrStyle();


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

    double weight = getNormalisationFactor(lumiData,iF<nMC?files[iF]:sigfile);

    std::ostringstream lselection,lname;
    std::string signal;
    if (sel==SELECTION::SIGNAL) signal = "(l1met>=60 && nvetomuons==0 && nvetoelectrons==0 && dijet_M>1120 && dijet_deta>4.1 && metnomu_significance>4.9 && alljetsmetnomu_mindphi>2.2 && passtrigger>0.5)";
    else if (sel==SELECTION::MuMu) signal = "(l1met>=60 && nselmuons==2 && dijet_M>1120 && dijet_deta>4.1 && metnomu_significance>4.9 && alljetsmetnomu_mindphi>2.2)";
    else if (sel==SELECTION::LNu) signal = "(l1met>=60 && ((nselelectrons==1  && nvetomuons==0) || (nselmuons==1  && nvetoelectrons==0) || (ntaus==1 && nvetoelectrons==0 && nvetomuons==0))&& dijet_M>1120 && dijet_deta>4.1 && metnomu_significance>4.9 && alljetsmetnomu_mindphi>2.2)";
    else if (sel==SELECTION::Top) signal = "(l1met>=60 && nselelectrons==1  && nselmuons==1 && dijet_M>1120 && dijet_deta>4.1 && metnomu_significance>4.9)";
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
    if (iF<nMC){
      hist[iF]->SetFillColor(color[iF]);
      hist[iF]->SetFillStyle(style[iF]);
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
  TLegend *leg = new TLegend(0.6,0.5,0.89,0.89);
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
  if (sel == SELECTION::SIGNAL){
    hist[nMC]->SetLineWidth(2);
    hist[nMC]->Draw("histsame");
    leg->AddEntry(hist[nMC],"VBF H, m_{H}=125 GeV","L");
  }
  leg->Draw("same");
  myc->Update();

  TLatex lat;
  char buf[200];
  lat.DrawLatexNDC(0.1,0.95,"CMSDAS@BARI 2015, CMS simulation");
  sprintf(buf,"#sqrt{s}=13 TeV, L=%3.0f fb^{-1}",lumiData/1000.);
  lat.DrawLatexNDC(0.15,0.85,buf);

  if (sel == SELECTION::SIGNAL)
    lat.DrawLatexNDC(0.15,0.75,"Veto e + veto #mu signal region");
  if (sel == SELECTION::MuMu)
    lat.DrawLatexNDC(0.15,0.75,"#mu+#mu Z control region");
  if (sel == SELECTION::LNu)
    lat.DrawLatexNDC(0.15,0.75,"(e,#mu,#tau) W control region");
  if (sel == SELECTION::Top)
    lat.DrawLatexNDC(0.15,0.75,"e+#mu top control region");

  std::ostringstream lsave;
  lsave.str("");
  lsave << "PLOTS/" << varName;
  if (dology) lsave << "_log";
  lsave << ".png";
  myc->Print(lsave.str().c_str());

  return 0;

}//main
