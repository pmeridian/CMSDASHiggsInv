#include <cstdlib>
#include <iostream>
#include <iomanip>
#include<fstream>
#include<sstream>
#include <map>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "LightTreeReader.hh"
#include "utilities.hh"


int main( int argc, char** argv )
{//main

  //////////////////////////////////
  //// PARAMETERS //////////////////
  //////////////////////////////////


  if (argc < 5) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEvts to process (0=all)>"
	      << " <path to input files>"
	      << " <path to output file>"
	      << " <name of dataset>"
	      << " <optional: debug (default=0)>"
	      << std::endl;
    return 1;
  }

  unsigned pNevts  = atoi(argv[1]);
  
  std::string inPath = argv[2];  
  std::string outPath = argv[3];
  std::string fileName = argv[4];

  unsigned debug = 0;
  if (argc>5) debug = atoi(argv[5]);

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << inPath << std::endl
	    << " -- Output file path: " << outPath << std::endl
	    << " -- Dataset name: " << fileName << std::endl;

  //////////////////////////////////
  //// INPUTS & OUTPUTS
  //////////////////////////////////


  std::ostringstream input;
  input << inPath << "/" << fileName << ".root";

  TFile *file = 0;
  if (!testFile(input.str(),file)) return 1;

  std::ostringstream output;
  output << outPath << "/" << fileName << ".root";

  TFile *outFile = 0;
  if (!testFile(output.str(),outFile,true)) return 1;

  TTree *outTree = 0;

  //////////////////////////////////
  //// GET TREE AND EVENT LOOP
  //////////////////////////////////

  LightTreeReader reader;

  reader.tree_ = (TTree*)file->Get("LightTree");
  if (!reader.tree_){
    std::cout << " -- Error, input tree cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  std::cout << " Tree opened successfully. " << std::endl;
  
  reader.SetBranchAddresses();

  const unsigned nEvts = ((pNevts > reader.tree_->GetEntries() || pNevts==0) ? static_cast<unsigned>(reader.tree_->GetEntries()) : pNevts) ;
  
  std::cout << " -- Processing = " << nEvts  << " events out of " << reader.tree_->GetEntries() << std::endl;
  
  //clone tree without copying the entries
  outTree = (TTree*)reader.tree_->CloneTree(0);

  unsigned nPass = 0;

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%1000 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    
    reader.tree_->GetEntry(ievt);

    if (debug) std::cout << " -- Mjj = " << reader.dijet_M_ << std::endl;

    //apply selection
    bool jetmetdphicut=reader.alljetsmetnomu_mindphi_>2.3;
    bool basesel= reader.jet1_eta_*reader.jet2_eta_<0 && reader.jet1_eta_<4.7 && reader.jet2_eta_<4.7 && reader.dijet_M_>=1200 && reader.jet1_pt_>50 && reader.dijet_deta_>3.6 && reader.jet2_pt_>45 && reader.metnomuons_>90 && reader.metnomu_significance_>4.0;
    bool pass = jetmetdphicut&&basesel;

    //fill output tree
    if (!pass) continue;

    outTree->Fill();
    nPass++;

  }//loop on entries

  if (debug) outTree->Print();
  outTree->AutoSave();

  std::cout << " -- Pass " << nPass << " tot " << nEvts << std::endl;

  return 0;
}//main
