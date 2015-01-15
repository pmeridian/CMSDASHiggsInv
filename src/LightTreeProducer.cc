#include "UserCode/CMSDASHiggsInv/interface/LightTreeProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"


using namespace edm;
using namespace reco;



LightTreeProducer::LightTreeProducer(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales")))
{ 
  outputTree_ = 0;
  
  run_=-1;
  lumi_=-1;
  event_=-1;
  weight_nolep_=1;
  total_weight_lepveto_ = 1;
  total_weight_leptight_ = 1;
  puweight_up_scale_=1;
  puweight_down_scale_=1;
  jet1_pt_ = 0;
  jet2_pt_ = 0;
  jet3_pt_=-1;
  jet1_E_ = 0;
  jet2_E_ = 0;
  jet3_E_ = 0;
  jet1_eta_ = 0;
  jet2_eta_ = 0;
  jet3_eta_ = 0;
  jet1_phi_ = 0;
  jet2_phi_ = 0;
  jet3_phi_ = 0;
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
  metnomu_x_ = 0;
  metnomu_y_ = 0;
  met_significance_ = 0;
  metnomu_significance_ = 0;
  sumet_ = 0;
  ht_ = 0;
  ht30_ = 0;
  mht_ = 0;
  sqrt_ht_ = 0;
  unclustered_et_ = 0;
  jet1met_dphi_ = 0;
  jet2met_dphi_ = 0;
  jet1metnomu_dphi_ = 0;
  jet2metnomu_dphi_ = 0;
  jetmet_mindphi_ = 0;
  jetmetnomu_mindphi_ = 0;
  alljetsmet_mindphi_ = 0;
  alljetsmetnomu_mindphi_ = 0;
  jetunclet_mindphi_ = 0;
  metunclet_dphi_ = 0;
  metnomuunclet_dphi_ = 0;
  dijetmet_scalarSum_pt_ = 0;
  dijetmet_vectorialSum_pt_ = 0;
  dijetmet_ptfraction_ = 0;
  jet1met_scalarprod_ = 0;
  jet2met_scalarprod_ = 0;
  dijetmetnomu_scalarSum_pt_ = 0;
  dijetmetnomu_vectorialSum_pt_ = 0;
  dijetmetnomu_ptfraction_ = 0;
  jet1metnomu_scalarprod_ = 0;
  jet2metnomu_scalarprod_ = 0;
  n_jets_cjv_30_ = 0;
  n_jets_cjv_20EB_30EE_ = 0;
  cjvjetpt_=-1;
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
  mu1_eta_=-10000;
  mu1_phi_=-1;
  mu2_pt_=-1;
  mu2_eta_=-10000;
  mu2_phi_=-1;
  ele1_pt_=-1;
  ele1_eta_=-10000;
  ele1_phi_=-1;
  tau1_pt_=-1;
  tau1_eta_=-10000;
  tau1_phi_=-1;
  lep_mt_=-1;
  n_vertices_=-1;

  
}
  
  
LightTreeProducer::~LightTreeProducer()
{
}


//
// member functions
//

void LightTreeProducer::endJob()
{
}

void LightTreeProducer::beginJob()
{
  std::cout << "--------------------------------------------- " << std::endl
	    << "------ PreAnalysis Info for LightTree ------" << std::endl
	    << "--------------------------------------------- " << std::endl;
  outputTree_=fs_->make<TTree>("LightTree","Tree containing LightTreeAna input variables");
  outputTree_->Branch("run",&run_);
  outputTree_->Branch("lumi",&lumi_);
  outputTree_->Branch("event",&event_);
  outputTree_->Branch("weight_nolep",&weight_nolep_);
  outputTree_->Branch("total_weight_lepveto",&total_weight_lepveto_);
  outputTree_->Branch("total_weight_leptight",&total_weight_leptight_);
  outputTree_->Branch("puweight_up_scale",&puweight_up_scale_);
  outputTree_->Branch("puweight_down_scale",&puweight_down_scale_);
  outputTree_->Branch("jet1_pt",&jet1_pt_);
  outputTree_->Branch("jet2_pt",&jet2_pt_);
  outputTree_->Branch("jet3_pt",&jet3_pt_);
  outputTree_->Branch("jet1_E",&jet1_E_);
  outputTree_->Branch("jet2_E",&jet2_E_);
  outputTree_->Branch("jet3_E",&jet3_E_);
  outputTree_->Branch("jet1_eta",&jet1_eta_);
  outputTree_->Branch("jet2_eta",&jet2_eta_);
  outputTree_->Branch("jet3_eta",&jet3_eta_);
  outputTree_->Branch("jet1_phi",&jet1_phi_);
  outputTree_->Branch("jet2_phi",&jet2_phi_);
  outputTree_->Branch("jet3_phi",&jet3_phi_);
  outputTree_->Branch("jet1_csv",&jet1_csv_);
  outputTree_->Branch("jet2_csv",&jet2_csv_);
  outputTree_->Branch("jet3_csv",&jet3_csv_);
  outputTree_->Branch("dijet_M",&dijet_M_);
  outputTree_->Branch("dijet_deta",&dijet_deta_);
  outputTree_->Branch("dijet_sumeta",&dijet_sumeta_);
  outputTree_->Branch("dijet_dphi",&dijet_dphi_);
  outputTree_->Branch("met",&met_);
  outputTree_->Branch("met_x",&met_x_);
  outputTree_->Branch("met_y",&met_y_);
  outputTree_->Branch("metnomu_x",&metnomu_x_);
  outputTree_->Branch("metnomu_y",&metnomu_y_);
  outputTree_->Branch("met_significance",&met_significance_);
  outputTree_->Branch("metnomu_significance",&metnomu_significance_);
  outputTree_->Branch("sumet",&sumet_);
  outputTree_->Branch("ht",&ht_);
  outputTree_->Branch("ht30",&ht30_);
  outputTree_->Branch("mht",&mht_);
  outputTree_->Branch("sqrt_ht",&sqrt_ht_);
  outputTree_->Branch("unclustered_et",&unclustered_et_);
  outputTree_->Branch("jet1met_dphi",&jet1met_dphi_);
  outputTree_->Branch("jet2met_dphi",&jet2met_dphi_);
  outputTree_->Branch("jet1metnomu_dphi",&jet1metnomu_dphi_);
  outputTree_->Branch("jet2metnomu_dphi",&jet2metnomu_dphi_);
  outputTree_->Branch("jetmet_mindphi",&jetmet_mindphi_);
  outputTree_->Branch("jetmetnomu_mindphi",&jetmetnomu_mindphi_);
  outputTree_->Branch("alljetsmet_mindphi",&alljetsmet_mindphi_);
  outputTree_->Branch("alljetsmetnomu_mindphi",&alljetsmetnomu_mindphi_);
  outputTree_->Branch("jetunclet_mindphi",&jetunclet_mindphi_);
  outputTree_->Branch("metunclet_dphi",&metunclet_dphi_);
  outputTree_->Branch("metnomuunclet_dphi",&metnomuunclet_dphi_);
  outputTree_->Branch("dijetmet_scalarSum_pt",&dijetmet_scalarSum_pt_);
  outputTree_->Branch("dijetmet_vectorialSum_pt",&dijetmet_vectorialSum_pt_);
  outputTree_->Branch("dijetmet_ptfraction",&dijetmet_ptfraction_);
  outputTree_->Branch("jet1met_scalarprod",&jet1met_scalarprod_);
  outputTree_->Branch("jet2met_scalarprod",&jet2met_scalarprod_);
  outputTree_->Branch("dijetmetnomu_scalarSum_pt",&dijetmetnomu_scalarSum_pt_);
  outputTree_->Branch("dijetmetnomu_vectorialSum_pt",&dijetmetnomu_vectorialSum_pt_);
  outputTree_->Branch("dijetmetnomu_ptfraction",&dijetmetnomu_ptfraction_);
  outputTree_->Branch("jet1metnomu_scalarprod",&jet1metnomu_scalarprod_);
  outputTree_->Branch("jet2metnomu_scalarprod",&jet2metnomu_scalarprod_);
  outputTree_->Branch("n_jets_cjv_30",&n_jets_cjv_30_);
  outputTree_->Branch("n_jets_cjv_20EB_30EE",&n_jets_cjv_20EB_30EE_);
  outputTree_->Branch("cjvjetpt",&cjvjetpt_);
  outputTree_->Branch("passtrigger",&passtrigger_);
  outputTree_->Branch("passparkedtrigger1",&passparkedtrigger1_);
  outputTree_->Branch("passparkedtrigger2",&passparkedtrigger2_);
  outputTree_->Branch("l1met",&l1met_);
  outputTree_->Branch("metnomuons",&metnomuons_);
  outputTree_->Branch("nvetomuons",&nvetomuons_);
  outputTree_->Branch("nselmuons",&nselmuons_);
  outputTree_->Branch("nvetoelectrons",&nvetoelectrons_);
  outputTree_->Branch("nselelectrons",&nselelectrons_);
  outputTree_->Branch("ntaus",&ntaus_);
  outputTree_->Branch("m_mumu",&m_mumu_);
  outputTree_->Branch("m_mumu_gen",&m_mumu_gen_);
  outputTree_->Branch("mu1_pt",&mu1_pt_);
  outputTree_->Branch("mu1_eta",&mu1_eta_);
  outputTree_->Branch("mu1_phi",&mu1_phi_);
  outputTree_->Branch("mu2_pt",&mu2_pt_);
  outputTree_->Branch("mu2_eta",&mu2_eta_);
  outputTree_->Branch("mu2_phi",&mu2_phi_);
  outputTree_->Branch("ele1_pt",&ele1_pt_);
  outputTree_->Branch("ele1_eta",&ele1_eta_);
  outputTree_->Branch("ele1_phi",&ele1_phi_);
  outputTree_->Branch("tau1_pt",&tau1_pt_);
  outputTree_->Branch("tau1_eta",&tau1_eta_);
  outputTree_->Branch("tau1_phi",&tau1_phi_);
  outputTree_->Branch("lep_mt",&lep_mt_);
  outputTree_->Branch("n_vertices",&n_vertices_);
}

// ------------ method called to for each event  ------------
void
LightTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  std::cout << "\n === TRIGGER PATHS === " << std::endl;
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    std::cout << "Trigger " << names.triggerName(i) << 
      ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
      ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") 
	      << std::endl;
  }
  std::cout << "\n === TRIGGER OBJECTS === " << std::endl;
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
    obj.unpackPathNames(names);
    std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
    // Print trigger object collection and type
    std::cout << "\t   Collection: " << obj.collection() << std::endl;
    std::cout << "\t   Type IDs:   ";
    for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
    std::cout << std::endl;
    // Print associated trigger filters
    std::cout << "\t   Filters:    ";
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
    std::cout << std::endl;
    std::vector<std::string> pathNamesAll  = obj.pathNames(false);
    std::vector<std::string> pathNamesLast = obj.pathNames(true);
    // Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
    // definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
    // means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
    std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
    for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
      bool isBoth = obj.hasPathName( pathNamesAll[h], true, true ); 
      bool isL3   = obj.hasPathName( pathNamesAll[h], false, true ); 
      bool isLF   = obj.hasPathName( pathNamesAll[h], true, false ); 
      bool isNone = obj.hasPathName( pathNamesAll[h], false, false ); 
      std::cout << "   " << pathNamesAll[h];
      if (isBoth) std::cout << "(L,3)";
      if (isL3 && !isBoth) std::cout << "(*,3)";
      if (isLF && !isBoth) std::cout << "(L,*)";
      if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;


  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &PV = vertices->front();

  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  for (const pat::Muon &mu : *muons) {
    if (mu.pt() < 5 || !mu.isLooseMuon()) continue;
    printf("muon with pt %4.1f, dz(PV) %+5.3f, POG loose id %d, tight id %d\n",
	   mu.pt(), mu.muonBestTrack()->dz(PV.position()), mu.isLooseMuon(), mu.isTightMuon(PV));
  }

  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);
  for (const pat::Electron &el : *electrons) {
    if (el.pt() < 5) continue;
    // printf("elec with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes), lost hits %d, pass conv veto %d\n",
    // 	   el.pt(), el.superCluster()->eta(), el.sigmaIetaIeta(), el.full5x5_sigmaIetaIeta(), el.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(), el.passConversionVeto());
  }

  edm::Handle<pat::PhotonCollection> photons;
  iEvent.getByToken(photonToken_, photons);
  for (const pat::Photon &pho : *photons) {
    if (pho.pt() < 20 or pho.chargedHadronIso()/pho.pt() > 0.3) continue;
    printf("phot with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes)\n",
	   pho.pt(), pho.superCluster()->eta(), pho.sigmaIetaIeta(), pho.full5x5_sigmaIetaIeta());
  }


  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(tauToken_, taus);
  for (const pat::Tau &tau : *taus) {
    if (tau.pt() < 20) continue;
    printf("tau  with pt %4.1f, dxy signif %.1f, ID(byMediumCombinedIsolationDeltaBetaCorr3Hits) %.1f, lead candidate pt %.1f, pdgId %d \n",
	   tau.pt(), tau.dxy_Sig(), tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"), tau.leadCand()->pt(), tau.leadCand()->pdgId());
  }


  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);
  int ijet = 0;
  for (const pat::Jet &j : *jets) {
    if (j.pt() < 20) continue;
    printf("jet  with pt %5.1f (raw pt %5.1f), eta %+4.2f, btag CSV %.3f, CISV %.3f, pileup mva disc %+.2f\n",
	   j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), std::max(0.f,j.bDiscriminator("combinedSecondaryVertexBJetTags")), std::max(0.f,j.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags")), j.userFloat("pileupJetId:fullDiscriminant"));
    if ((++ijet) == 1) { // for the first jet, let's print the leading constituents
      std::vector<reco::CandidatePtr> daus(j.daughterPtrVector());
      std::sort(daus.begin(), daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); }); // the joys of C++11
      for (unsigned int i2 = 0, n = daus.size(); i2 < n && i2 <= 3; ++i2) {
	const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*daus[i2]);
	printf("         constituent %3d: pt %6.2f, dz(pv) %+.3f, pdgId %+3d\n", i2,cand.pt(),cand.dz(PV.position()),cand.pdgId());
      }
    }
  }


  edm::Handle<pat::JetCollection> fatjets;
  iEvent.getByToken(fatjetToken_, fatjets);
  for (const pat::Jet &j : *fatjets) {
    printf("AK8j with pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f ungroomed, %5.1f pruned, %5.1f trimmed, %5.1f filtered. CMS TopTagger %.1f\n",
	   j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), j.mass(), j.userFloat("ak8PFJetsCHSPrunedLinks"), j.userFloat("ak8PFJetsCHSTrimmedLinks"), j.userFloat("ak8PFJetsCHSFilteredLinks"), j.userFloat("cmsTopTagPFJetsCHSLinksAK8"));
  }
 
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  const pat::MET &met = mets->front();
  printf("MET: pt %5.1f, phi %+4.2f, sumEt (%.1f). genMET %.1f. MET with JES up/down: %.1f/%.1f\n",
	 met.pt(), met.phi(), met.sumEt(),
	 met.genMET()->pt(), met.shiftedPt(pat::MET::JetEnUp), met.shiftedPt(pat::MET::JetEnDown));

  printf("\n");
  
  outputTree_->Fill();
} // pruneKids

//define this as a plug-in
DEFINE_FWK_MODULE(LightTreeProducer);
