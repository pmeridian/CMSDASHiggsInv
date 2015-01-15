#include "LightTreeReader.hh"
#include <iomanip>
#include <iostream>

LightTreeReader::LightTreeReader(){
  tree_ = 0;
  
  run_=-1;
  lumi_=-1;
  event_=-1;
  weight_nolep_=1;
  total_weight_lepveto_ = 1;
  total_weight_leptight_ = 1;
  jet1_pt_ = 0;
  jet2_pt_ = 0;
  jet1_E_ = 0;
  jet2_E_ = 0;
  jet1_eta_ = 0;
  jet2_eta_ = 0;
  jet1_csv_ = 0;
  jet2_csv_ = 0;
  jet3_csv_ = 0;
  dijet_M_ = 0;
  dijet_deta_ = 0;
  dijet_sumeta_ = 0;
  dijet_dphi_ = 0;
  met_ = 0;
  met_x_ = 0;
  met_y_ = 0;
  met_significance_ = 0;
  sumet_ = 0;
  ht_ = 0;
  mht_ = 0;
  sqrt_ht_ = 0;
  unclustered_et_ = 0;
  jet1met_dphi_ = 0;
  jet2met_dphi_ = 0;
  jetmet_mindphi_ = 0;
  jetunclet_mindphi_ = 0;
  metunclet_dphi_ = 0;
  dijetmet_scalarSum_pt_ = 0;
  dijetmet_vectorialSum_pt_ = 0;
  dijetmet_ptfraction_ = 0;
  jet1met_scalarprod_ = 0;
  jet2met_scalarprod_ = 0;
  n_jets_cjv_30_ = 0;
  n_jets_cjv_20EB_30EE_ = 0;
  cjvjetpt_=-1;
  jet3pt_=-1;
  passtrigger_ = -1;
  passparkedtrigger1_ = -1;
  passparkedtrigger2_ = -1;
  l1met_ = 0;
  metnomuons_ =0;
  nvetomuons_=0;
  nselmuons_=0;
  nvetoelectrons_=0;
  nselelectrons_=0;
  ntaus_=0;
  m_mumu_=-1;
  m_mumu_gen_=-1;
  mu1_pt_=-1;
  mu1_eta_=-1;
  ele1_pt_=-1;
  ele1_eta_=-1;
  
}

LightTreeReader::~LightTreeReader(){
  ;
}

void LightTreeReader::SetBranchAddresses(){
  if (!tree_) {
    std::cout << " -- ERROR! Tree pointer is 0." << std::endl;
    exit(1);
  }
  //tree_->SetBranchAddress("run",&run_);
  //tree_->SetBranchAddress("lumi",&lumi_);
  //tree_->SetBranchAddress("event",&event_);
  tree_->SetBranchAddress("weight_nolep",&weight_nolep_);
  tree_->SetBranchAddress("total_weight_lepveto",&total_weight_lepveto_);
  tree_->SetBranchAddress("total_weight_leptight",&total_weight_leptight_);
  tree_->SetBranchAddress("jet1_pt",&jet1_pt_);
  tree_->SetBranchAddress("jet2_pt",&jet2_pt_);
  tree_->SetBranchAddress("jet1_E",&jet1_E_);
  tree_->SetBranchAddress("jet2_E",&jet2_E_);
  tree_->SetBranchAddress("jet1_eta",&jet1_eta_);
  tree_->SetBranchAddress("jet2_eta",&jet2_eta_);
  tree_->SetBranchAddress("jet1_csv",&jet1_csv_);
  tree_->SetBranchAddress("jet2_csv",&jet2_csv_);
  tree_->SetBranchAddress("jet3_csv",&jet3_csv_);
  tree_->SetBranchAddress("dijet_M",&dijet_M_);
  tree_->SetBranchAddress("dijet_deta",&dijet_deta_);
  tree_->SetBranchAddress("dijet_sumeta",&dijet_sumeta_);
  tree_->SetBranchAddress("dijet_dphi",&dijet_dphi_);
  tree_->SetBranchAddress("met",&met_);
  tree_->SetBranchAddress("met_x",&met_x_);
  tree_->SetBranchAddress("met_y",&met_y_);
  tree_->SetBranchAddress("met_significance",&met_significance_);
  tree_->SetBranchAddress("sumet",&sumet_);
  tree_->SetBranchAddress("ht",&ht_);
  tree_->SetBranchAddress("mht",&mht_);
  tree_->SetBranchAddress("sqrt_ht",&sqrt_ht_);
  tree_->SetBranchAddress("unclustered_et",&unclustered_et_);
  tree_->SetBranchAddress("jet1met_dphi",&jet1met_dphi_);
  tree_->SetBranchAddress("jet2met_dphi",&jet2met_dphi_);
  tree_->SetBranchAddress("jetmet_mindphi",&jetmet_mindphi_);
  tree_->SetBranchAddress("jetunclet_mindphi",&jetunclet_mindphi_);
  tree_->SetBranchAddress("metunclet_dphi",&metunclet_dphi_);
  tree_->SetBranchAddress("dijetmet_scalarSum_pt",&dijetmet_scalarSum_pt_);
  tree_->SetBranchAddress("dijetmet_vectorialSum_pt",&dijetmet_vectorialSum_pt_);
  tree_->SetBranchAddress("dijetmet_ptfraction",&dijetmet_ptfraction_);
  tree_->SetBranchAddress("jet1met_scalarprod",&jet1met_scalarprod_);
  tree_->SetBranchAddress("jet2met_scalarprod",&jet2met_scalarprod_);
  tree_->SetBranchAddress("n_jets_cjv_30",&n_jets_cjv_30_);
  tree_->SetBranchAddress("n_jets_cjv_20EB_30EE",&n_jets_cjv_20EB_30EE_);
  tree_->SetBranchAddress("cjvjetpt",&cjvjetpt_);
  tree_->SetBranchAddress("jet3pt",&jet3pt_);
  tree_->SetBranchAddress("passtrigger",&passtrigger_);
  tree_->SetBranchAddress("passparkedtrigger1",&passparkedtrigger1_);
  tree_->SetBranchAddress("passparkedtrigger2",&passparkedtrigger2_);
  tree_->SetBranchAddress("l1met",&l1met_);
  tree_->SetBranchAddress("metnomuons",&metnomuons_);

  tree_->SetBranchAddress("nvetomuons",&nvetomuons_);
  tree_->SetBranchAddress("nselmuons",&nselmuons_);
  tree_->SetBranchAddress("nvetoelectrons",&nvetoelectrons_);
  tree_->SetBranchAddress("nselelectrons",&nselelectrons_);
  tree_->SetBranchAddress("ntaus",&ntaus_);
  tree_->SetBranchAddress("m_mumu",&m_mumu_);
  tree_->SetBranchAddress("m_mumu_gen",&m_mumu_gen_);
  tree_->SetBranchAddress("mu1_pt",&mu1_pt_);
  tree_->SetBranchAddress("mu1_eta",&mu1_eta_);
  tree_->SetBranchAddress("ele1_pt",&ele1_pt_);
  tree_->SetBranchAddress("ele1_eta",&ele1_eta_);
}
