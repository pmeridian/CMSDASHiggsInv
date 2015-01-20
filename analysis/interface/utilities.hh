#include <string>
#include "TFile.h"
#include "SimpleParamParser.h"



bool testFile(std::string input, TFile* & file, const bool recreate=false){
  if (recreate) file = TFile::Open(input.c_str(),"RECREATE");
  else file = TFile::Open(input.c_str());
  
  if (!file) {
    std::cout << " -- Error, file " << input.c_str() << " cannot be opened. Skipping..." << std::endl;
    return false;
  }
  else std::cout << " -- file " << file->GetName() << " successfully opened." << std::endl;
  return true;
};


double getNormalisationFactor(const double lumi, std::string dataset){
  //Get lumi xs and events from params file
  SimpleParamParser parser;
  std::string input_params_ = "data/Params13TeV.dat";
  //std::cout << "** Parsing parameter file... **" << input_params_ << std::endl;
  parser.ParseFile(input_params_);
  //std::cout<<"parsed"<<std::endl;
  double xs=parser.GetParam<double>("XS_"+dataset);
  //std::cout<<"got xs"<<std::endl;
  double events=parser.GetParam<double>("EVT_"+dataset);
  //std::cout<<"got events"<<std::endl;

  //std::cout<<"XS is: "<<xs<<"pb"<<std::endl;
  //std::cout<<"EVT is: "<<events<<std::endl;
  //std::cout<<"LUMI is: "<<lumi<<"pb^-1"<<std::endl;
  
  double lumixsweight=xs*lumi/events;
  std::cout<<"-- dataset " << dataset << ", LUMIXSWEIGHT is: "<<lumixsweight<<std::endl;
  return lumixsweight;
};
