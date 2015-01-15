#ifndef UserCode_CMSDASHiggsInv_LightTreeReader_h
#define UserCode_CMSDASHiggsInv_LightTreeReader_h

#include <string>

#include "TTree.h"

class LightTreeReader{

public:
  LightTreeReader();
  ~LightTreeReader();

  void SetBranchAddresses();

  TTree *tree_;
    
  unsigned run_;
  unsigned lumi_;
  unsigned event_;
  double weight_nolep_;
  double total_weight_lepveto_;
  double total_weight_leptight_;
  double jet1_pt_;
  double jet2_pt_;
  double jet1_E_;
  double jet2_E_;
  double jet1_eta_;
  double jet2_eta_;
  double jet1_csv_;
  double jet2_csv_;
  double jet3_csv_;
  double dijet_M_;
  double dijet_deta_;
  double dijet_sumeta_;
  double dijet_dphi_;
  double met_;
  double met_x_;
  double met_y_;
  double met_significance_;
  double sumet_;
  double ht_;
  double mht_;
  double sqrt_ht_;
  double unclustered_et_;
  double jet1met_dphi_;
  double jet2met_dphi_;
  double jetmet_mindphi_;
  double jetunclet_mindphi_;
  double metunclet_dphi_;
  double dijetmet_scalarSum_pt_;
  double dijetmet_vectorialSum_pt_;
  double dijetmet_ptfraction_;
  double jet1met_scalarprod_;
  double jet2met_scalarprod_;
  double cjvjetpt_;
  double jet3pt_;
  unsigned n_jets_cjv_30_;
  unsigned n_jets_cjv_20EB_30EE_;
  double passtrigger_;
  double passparkedtrigger1_;
  double passparkedtrigger2_;
  double l1met_;
  double metnomuons_;
  
  int nvetomuons_;
  int nselmuons_;
  int nvetoelectrons_;
  int nselelectrons_;
  int ntaus_;
  double m_mumu_;
  double m_mumu_gen_;
  double mu1_pt_;
  double mu1_eta_;
  double ele1_pt_;
  double ele1_eta_;
  
private:

  
};

#endif
