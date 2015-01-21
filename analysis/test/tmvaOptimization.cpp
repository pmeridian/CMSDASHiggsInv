#include <cstdlib>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include <TTree.h>
#include <TChain.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TAxis.h>
#include <TString.h>
#include "TObjString.h"
#include "TSystem.h"


#include "TMVAGui.C"
#include "utilities.hh"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Config.h"
#endif


using namespace std;

int main(int argc, char* argv[]) {

  if (argc != 3) {
    std::cout << " Usage: " 
              << " <path to input files>"
	      << " <path to output MVA file>"
	      << std::endl;
    return 1;
  }

  std::string inPath(argv[1]);
  std::string outPath(argv[2]);

  std::cout << " -- Input path: " << inPath << std::endl;
  std::cout << " -- Output file: " << outPath << std::endl;

  // This loads the library
  TMVA::Tools::Instance();
  // --- Here the preparation phase begins
  

  TString outfileName( outPath );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  
  // Create the factory object. Later you can choose the methods
  // whose performance you'd like to investigate. The factory is 
  // the only TMVA object you have to interact with
  //
  // The first argument is the base of the name of all the
  // weightfiles in the directory weight/
  //
  // The second argument is the output file for the training results
  // All TMVA output can be suppressed by removing the "!" (not) in
  // front of the "Silent" argument in the option string

   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
   
   factory->AddVariable("jet1_pt","jet1_pt","", 'D');	
   factory->AddVariable("jet2_pt","jet2_pt","", 'D');
   factory->AddVariable("dijet_M","dijet_M","", 'D');
   factory->AddVariable("dijet_deta","dijet_deta","", 'D');
   factory->AddVariable("metnomu_significance","metnomu_significance","",'D');
   factory->AddVariable("alljetsmetnomu_mindphi","alljetsmetnomu_mindphi","",'D');

   //   factory->AddSpectator( "nPhot_presel", "nPhot_presel", "", 'F' );

   std::vector<std::string> backgrounds;

   backgrounds.push_back("WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauol");
   backgrounds.push_back("WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola");
   backgrounds.push_back("WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola");
   backgrounds.push_back("WJetsToLNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola");

   backgrounds.push_back("ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola");
   backgrounds.push_back("ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola");
   backgrounds.push_back("ZJetsToNuNu_HT-400to600_Tune4C_13TeV-madgraph-tauola");
   backgrounds.push_back("ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola");

   std::vector<std::string> signals;
   signals.push_back("VBF_HToInv_M-125_13TeV_powheg-pythia6");
   double lumiData = 10000;//in pb-1

   for (int i=0; i<signals.size(); i++){
     float weight = getNormalisationFactor(lumiData,signals[i]);

     TFile* f=TFile::Open(Form("%s/%s.root",inPath.c_str(),signals[i].c_str()));
     TTree* sig=(TTree*) f->Get("lightTree/LightTree");
     if (!sig)
       {
	 std::cout << "====> ERROR: Sig tree " << signals[i] << " cannot be found" << std::endl;
	 continue;
       }
	 
     factory->AddSignalTree    ( sig, weight);
   } 

   for (int i=0; i<backgrounds.size(); i++){
     float weight = getNormalisationFactor(lumiData,backgrounds[i]);

     TFile* f=TFile::Open(Form("%s/%s.root",inPath.c_str(),backgrounds[i].c_str()));
     TTree* bkg=(TTree*) f->Get("lightTree/LightTree");
     if (!bkg)
       {
	 std::cout << "====> ERROR: Bkg tree " << backgrounds[i] << " cannot be found" << std::endl;
	 continue;
       }
	 
     factory->AddBackgroundTree    ( bkg, weight);
   } 


   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts;
   TCut mycutb;
   
   mycuts = "passtrigger==1 && nvetomuons==0 && nvetoelectrons==0 && ntaus==0 && metnomuons>90 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && dijet_M>0 && jet1_eta*jet2_eta<0"; 
   mycutb = "passtrigger==1 && nvetomuons==0 && nvetoelectrons==0 && ntaus==0 && metnomuons>90 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && dijet_M>0 && jet1_eta*jet2_eta<0";
 
   factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

   // ---- Book MVA methods
   // factory->BookMethod( TMVA::Types::kCuts, "Cuts",
   // 			//			"!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp[0]=FSmart:VarProp[1]=FSmart:VarProp[2]=FSmart:VarProp[3]=FSmart:VarProp[4]=FSmart:VarProp[5]=FSmart" );
    factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
    			"H:!V:FitMethod=GA:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95:VarProp[0]=FSmart:VarProp[1]=FSmart:VarProp[2]=FSmart:VarProp[3]=FSmart:VarProp[4]=FSmart:VarProp[5]=FSmart" );
   // factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
   // 			"!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );
   //factory->BookMethod( TMVA::Types::kBDT, "BDT","!H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs
   // Train MVAs using the set of training events
   factory->TrainAllMethods();
   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();
   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();
   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVAGui( outfileName );
}

