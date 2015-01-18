#include "UserCode/CMSDASHiggsInv/interface/LightTreeProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Candidate/interface/CompositeCandidate.h"


#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"

#include <boost/foreach.hpp>

#include "DataFormats/Math/interface/deltaR.h"
using namespace edm;
using namespace reco;

bool LightTreeProducer::MinDRToCollection(reco::Candidate const* cand,std::vector<const reco::Candidate*>& coll, double cut) {
  BOOST_FOREACH(const reco::Candidate* cand2, coll) {
    if ( deltaR(cand->p4(), cand2->p4()) < cut) return false;
  }
  return true;
}

bool LightTreeProducer::pu_id_mva_loose(const pat::Jet& j) {
  // Pt2030_Loose   = cms.vdouble(-0.80,-0.85,-0.84,-0.85),                                                                                                                    
  // Pt3050_Loose   = cms.vdouble(-0.80,-0.74,-0.68,-0.77)                                                                                                                     
  // #4 Eta Categories  0-2.5 2.5-2.75 2.75-3.0 3.0-5.0                                                                                                                        
  double abs_eta = fabs(j.eta());
  double pt = j.pt();
  double pu_id_mva_value_ = j.userFloat("pileupJetId:fullDiscriminant");
  if (pt > 20. && pt <= 30) {
    if (abs_eta < 2.5) {
      return (pu_id_mva_value_ > -0.80);
    } else if (abs_eta < 2.75) {
      return (pu_id_mva_value_ > -0.85);
    } else if (abs_eta < 3.0) {
      return (pu_id_mva_value_ > -0.84);
    } else if (abs_eta < 5.0) {
      return (pu_id_mva_value_ > -0.85);
    } else return true;
  } else if (pt > 30.) {
    if (abs_eta < 2.5) {
      return (pu_id_mva_value_ > -0.80);
    } else if (abs_eta < 2.75) {
      return (pu_id_mva_value_ > -0.74);
    } else if (abs_eta < 3.0) {
      return (pu_id_mva_value_ > -0.68);
    } else if (abs_eta < 5.0) {
      return (pu_id_mva_value_ > -0.77);
    } else return true;
  } else return true;
}

LightTreeProducer::LightTreeProducer(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
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

  //objects collections

  //leptons
  std::vector<const pat::Muon*> vetomuons;
  std::vector<const pat::Muon*> selmuons;
  std::vector<const pat::Electron*> vetoelectrons;
  std::vector<const pat::Electron*> selelectrons;
  std::vector<const pat::Tau*> seltaus;
  
  //jets
  std::vector<const Candidate*> const l1met;
  std::vector<const pat::Jet*> allJets;
  std::vector<const pat::Jet*> selJets;
  std::vector<CompositeCandidate> dijet_vec;



  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &PV = vertices->front();

  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  for (const pat::Muon &mu : *muons) {
    if (mu.pt() < 10 ) 
      continue;
    if ( fabs(mu.eta()) > 2.5 )
      continue;

    //VetoMuon Selection
    if (!mu.isPFMuon())
      continue;
    if (! (mu.isGlobalMuon() || mu.isTrackerMuon() ) )
      continue;
    if ( ((mu.chargedHadronIso()+mu.neutralHadronIso()+mu.photonIso()-0.5*mu.puChargedHadronIso())/mu.pt() ) > 0.2 )
      continue;
    vetomuons.push_back(&mu);

    //TightMuon selection
    if (!mu.isTightMuon(PV))
      continue;
    if ( ((mu.chargedHadronIso()+mu.neutralHadronIso()+mu.photonIso()-0.5*mu.puChargedHadronIso())/mu.pt() ) > 0.12 )
      continue;
    selmuons.push_back(&mu);
  }
  printf("NMUONS %d %d\n",int(selmuons.size()),int(vetomuons.size()));

  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);
  for (const pat::Electron &el : *electrons) {
    if (el.pt()<10)
      continue;
    if ( fabs(el.eta()) > 2.5 )
      continue;

    //Remove overlaps with muons    
    bool muOverlap=false;
    BOOST_FOREACH(const pat::Muon* mu, selmuons) {
      if ( deltaR(el.p4(), mu->p4()) < 0.3) 
	muOverlap=true;
    }
    if (muOverlap)
      continue;


    if(el.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-veto")>0.5)
      vetoelectrons.push_back(&el);
    if(el.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-medium")>0.5)
      selelectrons.push_back(&el);
  }
  printf("NELE %d %d\n",int(selelectrons.size()),int(vetoelectrons.size()));


  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(tauToken_, taus);
  for (const pat::Tau &tau : *taus) 
    {
      if (tau.pt() < 20) continue;
      if (fabs(tau.eta()) > 2.3) continue;
      if (tau.tauID("decayModeFinding")<0.5)
	continue;
      if (tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")<0.5)
	continue;
      if (tau.tauID("againstMuonLoose3")<0.5)
	continue;
      if (tau.tauID("againstElectronLooseMVA5")<0.5)
	continue;
      seltaus.push_back(&tau);
    }
  printf("NTAU %d\n",int(seltaus.size()));

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);
  // int ijet = 0;
  for (const pat::Jet &j : *jets) {
    if (j.pt() < 20) continue;
    if (fabs(j.eta()) > 4.7) continue;

    allJets.push_back(&j);
    //Loose JetID
    if (j.numberOfDaughters()<=1)
      continue;
    if (j.neutralHadronEnergyFraction()>0.99)
      continue;
    if (j.neutralEmEnergyFraction()>0.99)
      continue;
    if(!(j.chargedEmEnergyFraction()< 0.99 || fabs(j.eta()) >= 2.4))
      continue;
    if(!(j.chargedHadronEnergyFraction() > 0. || fabs(j.eta()) >= 2.4))
      continue;
    if(!(j.chargedMultiplicity() > 0 || fabs(j.eta()) >= 2.4))
      continue;

    //Remove overlaps with leptons
    bool muOverlap=false;
    bool eleOverlap=false;
    bool tauOverlap=false;
    BOOST_FOREACH(const pat::Muon* mu, vetomuons) {
      if ( deltaR(j.p4(), mu->p4()) < 0.5) 
	muOverlap=true;
    }
    BOOST_FOREACH(const pat::Electron* ele, vetoelectrons) {
      if ( deltaR(j.p4(), ele->p4()) < 0.5) 
	eleOverlap=true;
    }
    BOOST_FOREACH(const pat::Tau* tau, seltaus) {
      if ( deltaR(j.p4(), tau->p4()) < 0.5) 
	eleOverlap=true;
    }

    if (muOverlap || eleOverlap || tauOverlap)
      continue;

    //PU JetID    
    if ( ! pu_id_mva_loose(j) )
      continue;

    selJets.push_back(&j);
  }

  printf("NJETS %d %d\n",int(allJets.size()),int(selJets.size()));

  std::sort(vetomuons.begin(),vetomuons.end(),RefGreaterByPt<pat::Muon>());
  std::sort(selmuons.begin(),selmuons.end(),RefGreaterByPt<pat::Muon>());
  std::sort(vetoelectrons.begin(),vetoelectrons.end(),RefGreaterByPt<pat::Electron>());
  std::sort(selelectrons.begin(),selelectrons.end(),RefGreaterByPt<pat::Electron>());
  std::sort(seltaus.begin(),seltaus.end(),RefGreaterByPt<pat::Tau>());
  std::sort(selJets.begin(),selJets.end(),RefGreaterByPt<pat::Jet>());
  std::sort(allJets.begin(),allJets.end(),RefGreaterByPt<pat::Jet>());
  
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  const pat::MET &met = mets->front();

  Candidate::LorentzVector metnomuons=met.p4();
  BOOST_FOREACH(const pat::Muon* mu, selmuons) {
    metnomuons+=mu->p4();
  }
  printf("MET: pt %5.1f, phi %+4.2f METnomuons: pt %5.1f, phi %+4.2f ",
	 met.pt(), met.phi(), metnomuons.pt(), metnomuons.phi());
  printf("\n");
  
  /*
  EventInfo const* eventInfo = event->GetPtr<EventInfo>("eventInfo");
  run_= eventInfo->run();
  lumi_= eventInfo->lumi_block();
  event_= eventInfo->event();
  n_vertices_=eventInfo->good_vertices();
  
  if (is_data_) {
    
    TriggerPathPtrVec const& triggerPathPtrVec =
      event->GetPtrVec<TriggerPath>("triggerPathPtrVec","triggerPaths");
    //EventInfo const* eventInfo = event->GetPtr<EventInfo>("eventInfo"); //Can be used in future, but commented out to remove compiler warnings      
    //unsigned run = eventInfo->run(); //Can be used in future, but commented out to remove compiler warnings                                         
    passtrigger_=-1;
    passparkedtrigger1_=-1;
    passparkedtrigger2_=-1;
    for (unsigned i = 0; i < triggerPathPtrVec.size(); ++i) {
      std::string name = triggerPathPtrVec[i]->name();
      triggerPathPtrVec[i]->prescale();
      if (name.find(trigger_path_) != name.npos) passtrigger_ = 1;
      if (name.find("HLT_DiJet35_MJJ700_AllJets_DEta3p5_VBF") != name.npos) passparkedtrigger1_ = 1;
      if (name.find("HLT_DiJet30_MJJ700_AllJets_DEta3p5_VBF") != name.npos) passparkedtrigger2_ = 1;
    }
    if(dotrigskim_){
	  if(!(passtrigger_==1||passparkedtrigger1_==1||passparkedtrigger2_==1)){
	    return 1;
	  }
	}
      }
      //for MC                                                                                                                                       
      else {
	passtrigger_=-1;
	passparkedtrigger1_=-1;
	passparkedtrigger2_=-1;
	std::vector<TriggerObject *> const& objs = event->GetPtrVec<TriggerObject>(trig_obj_label_);
	if (objs.size() > 0) passtrigger_=1;
      } // do obj match                                                                            


    double wt = eventInfo->total_weight();
    double vetowt=1;
    double tightwt=1;
    double pileupwt=1;
    double pileupwtup=1;
    double pileupwtdown=1;
    if(!is_data_){
      vetowt= eventInfo->weight("idisoVeto");
      tightwt = eventInfo->weight("idisoTight");
      pileupwt=eventInfo->weight("pileup");
      pileupwtup=eventInfo->weight("pileup_up");
      pileupwtdown=eventInfo->weight("pileup_down");
    }
    puweight_up_scale_=pileupwtup/pileupwt;
    puweight_down_scale_=pileupwtdown/pileupwt;


    if(!ignoreLeptons_) nvetomuons_=vetomuons.size();
    else nvetomuons_=0;
    nselmuons_=selmuons.size();
    nvetoelectrons_=vetoelectrons.size();
    nselelectrons_=selelectrons.size();
    ntaus_=seltaus.size();

    std::sort(selmuons.begin(), selmuons.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));
    std::sort(selelectrons.begin(), selelectrons.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));

    //Get MT
    if(nselmuons_==1&&nselelectrons_==0&&ntaus_==0){//If 1 muon and no electrons use muon to make mt
      lep_mt_=sqrt(2*selmuons[0]->pt()*met->pt()*(1-cos(selmuons[0]->phi()-met->phi())));
    }
    else if(nselelectrons_==1&&nselmuons_==0&&ntaus_==0){//If 1 electron and no muons use electron to make mt
      lep_mt_=sqrt(2*selelectrons[0]->pt()*met->pt()*(1-cos(selelectrons[0]->phi()-met->phi())));      
    }
    else if(ntaus_==1&&nselelectrons_==0&&nselmuons_==0){//If 1 electron and no muons use electron to make mt
      lep_mt_=sqrt(2*seltaus[0]->pt()*met->pt()*(1-cos(seltaus[0]->phi()-met->phi())));      
    }
    else{//If otherwise set mt to dummy value
      lep_mt_=-2;
    }

    if(nselmuons_>=1){
      mu1_pt_=selmuons[0]->pt();
      mu1_eta_=selmuons[0]->eta();
      mu1_phi_=selmuons[0]->phi();
      if(nselmuons_>=2){
	mu2_pt_=selmuons[1]->pt();
	mu2_eta_=selmuons[1]->eta();
	mu2_phi_=selmuons[1]->phi();
      }
      else{
	mu2_pt_=-1;
	mu2_eta_=9999999;
	mu2_phi_=9999999;
      }
    }
    else{
      mu1_pt_=-1;
      mu1_eta_=9999999;
      mu1_phi_=9999999;
    }
    if(nselelectrons_>=1){
      ele1_pt_=selelectrons[0]->pt();
      ele1_eta_=selelectrons[0]->eta();
      ele1_phi_=selelectrons[0]->phi();
    }
    else{
      ele1_pt_=-1;
      ele1_eta_=9999999;
      ele1_phi_=9999999;
    }

    if(ntaus_>=1){
      tau1_pt_=seltaus[0]->pt();
      tau1_eta_=seltaus[0]->eta();
      tau1_phi_=seltaus[0]->phi();
    }
    else{
      tau1_pt_=-1;
      tau1_eta_=9999999;
      tau1_phi_=9999999;
    }

    if(nselmuons_==2){
      m_mumu_=((selmuons.at(0)->vector())+(selmuons.at(1)->vector())).M();
    }
    else m_mumu_=-1;

    //Get gen z mass
    int ngenmuplus=0;
    int ngenmuminus=0;
    m_mumu_gen_=-1;
    if(!is_data_){
      std::vector<GenParticle*> const& parts = event->GetPtrVec<GenParticle>("genParticles");
      GenParticle* lepplus = 0;
      GenParticle* lepminus = 0;
      
      for (unsigned i = 0; i < parts.size(); ++i) {
	if (parts[i]->status() != 3) continue;
	
	int id = parts[i]->pdgid();

	if (id == static_cast<int>(13)) {
	  lepminus = parts[i];
	  ngenmuminus++;
	}
	if (id == static_cast<int>(-13)) {
	lepplus = parts[i];
	ngenmuplus++;
	}  
      }//loop on genparticles                                                                                                                                  
      
      if (ngenmuminus==1&&ngenmuplus==1) {
	m_mumu_gen_ = (lepplus->vector()+lepminus->vector()).M();
      }
    }

    if (dijet_vec.size() != 0) {
      
      CompositeCandidate const* dijet = dijet_vec.at(0);

      Candidate const* jet1 = dijet->GetCandidate("jet1");
      Candidate const* jet2 = dijet->GetCandidate("jet2");
      ROOT::Math::PtEtaPhiEVector jet1vec = jet1->vector();
      ROOT::Math::PtEtaPhiEVector jet2vec = jet2->vector();
      ROOT::Math::PtEtaPhiEVector metvec = met->vector();
      ROOT::Math::PtEtaPhiEVector metnomuvec = metnomuons->vector();

      weight_nolep_ = wt;
      total_weight_lepveto_ =wt*vetowt;
      total_weight_leptight_=wt*tightwt;
      
      jet1_pt_ = jet1->pt();
      jet2_pt_ = jet2->pt();
      jet1_E_ = jet1vec.E();
      jet2_E_ = jet2vec.E();
      jet1_eta_ = jet1->eta();
      jet2_eta_ = jet2->eta();
      jet1_phi_ = jet1->phi();
      jet2_phi_ = jet2->phi();
      dijet_M_ = dijet->M();
      dijet_deta_ = fabs(jet1->eta() - jet2->eta());
      dijet_sumeta_ = jet1->eta() + jet2->eta();
      dijet_dphi_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(jet1vec,jet2vec));
      met_ = met->pt();
      met_x_ = metvec.Px();
      met_y_ = metvec.Py();
      met_significance_ = met->et_sig();
      sumet_ = met->sum_et();
      if(l1met.size()==1){
	l1met_ = l1met[0]->energy();
      }
      metnomuons_ = metnomuons->pt();
      metnomu_x_ = metnomuvec.Px();
      metnomu_y_ = metnomuvec.Py();
      metnomu_significance_ = met_significance_/met_*metnomuons_;

      double ht =0;
      double ht30 =0;
      ROOT::Math::PtEtaPhiEVector mhtVec(0,0,0,0);
      for(unsigned i =0; i<jets.size();++i){
	ht+=jets[i]->vector().Et();
	if(jets[i]->pt()>30)	ht30+=jets[i]->vector().Et();
	mhtVec += jets[i]->vector();
      }
      ROOT::Math::PtEtaPhiEVector unclVec = mhtVec + metvec;

      ht_ = ht;
      ht30_=ht30;
      mht_ = mhtVec.Et();
      sqrt_ht_ = sqrt(ht);
      unclustered_et_ = unclVec.Et();

      double dphi1 = fabs(ROOT::Math::VectorUtil::DeltaPhi(jet1vec,metvec));
      double dphi2 = fabs(ROOT::Math::VectorUtil::DeltaPhi(jet2vec,metvec));
      double nomudphi1 = fabs(ROOT::Math::VectorUtil::DeltaPhi(jet1vec,metnomuvec));
      double nomudphi2 = fabs(ROOT::Math::VectorUtil::DeltaPhi(jet2vec,metnomuvec));
      jet1met_dphi_ = dphi1;
      jet2met_dphi_ = dphi2;
      jet1metnomu_dphi_ = nomudphi1;
      jet2metnomu_dphi_ = nomudphi2;
      jetmet_mindphi_ = std::min(dphi1,dphi2);
      jetmetnomu_mindphi_ = std::min(nomudphi1,nomudphi2);
      

      dijetmet_scalarSum_pt_ = jet1->pt()+jet2->pt()+met->pt();
      dijetmet_vectorialSum_pt_ = (jet1vec+jet2vec+metvec).Pt();
      dijetmet_ptfraction_ = dijet->pt()/(dijet->pt()+met->pt());
      dijetmetnomu_scalarSum_pt_ = jet1->pt()+jet2->pt()+metnomuons->pt();
      dijetmetnomu_vectorialSum_pt_ = (jet1vec+jet2vec+metnomuvec).Pt();
      dijetmetnomu_ptfraction_ = dijet->pt()/(dijet->pt()+metnomuons->pt());

      jet1met_scalarprod_ = (jet1vec.Px()*met_x_+jet1vec.Py()*met_y_)/met_;
      jet2met_scalarprod_ = (jet2vec.Px()*met_x_+jet2vec.Py()*met_y_)/met_;

      jet1metnomu_scalarprod_ = (jet1vec.Px()*metnomu_x_+jet1vec.Py()*metnomu_y_)/met_;
      jet2metnomu_scalarprod_ = (jet2vec.Px()*metnomu_x_+jet2vec.Py()*metnomu_y_)/met_;

      jetunclet_mindphi_ = std::min(fabs(ROOT::Math::VectorUtil::DeltaPhi(jet1vec,unclVec)),
				    fabs(ROOT::Math::VectorUtil::DeltaPhi(jet2vec,unclVec)));
      metunclet_dphi_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(unclVec,metvec));
      metnomuunclet_dphi_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(unclVec,metnomuvec));

      double eta_high = (jet1->eta() > jet2->eta()) ? jet1->eta() : jet2->eta();
      double eta_low = (jet1->eta() > jet2->eta()) ? jet2->eta() : jet1->eta();
      n_jets_cjv_30_ = 0;
      n_jets_cjv_20EB_30EE_ = 0;
      jet3_pt_=-1;
      jet3_eta_=-10000;
      jet3_phi_=-10000;
      cjvjetpt_=-1;
      alljetsmetnomu_mindphi_=jetmetnomu_mindphi_;
      alljetsmet_mindphi_=jetmet_mindphi_;
      if (jets.size() > 2) {
	for (unsigned i = 0; i < jets.size(); ++i) {
	  if(jets[i]->id()==jet1->id()){
	    jet1_csv_=jets[i]->GetBDiscriminator("combinedSecondaryVertexBJetTags");
	  }
	  if(jets[i]->id()==jet2->id()){
	    jet2_csv_=jets[i]->GetBDiscriminator("combinedSecondaryVertexBJetTags");	  
	  }

	  bool isInCentralGap = fabs(jets[i]->eta())<4.7 && jets[i]->eta() > eta_low && jets[i]->eta() < eta_high;
	  double tmppt=jets[i]->pt();
	  if(isInCentralGap&&(tmppt>cjvjetpt_)){
	    cjvjetpt_=tmppt;
	  }
	  if(tmppt>jet3_pt_&&(jets[i]->id()!=jet1->id())&&(jets[i]->id()!=jet2->id())){
	    jet3_csv_=jets[i]->GetBDiscriminator("combinedSecondaryVertexBJetTags");
	    jet3_pt_=tmppt;
	    jet3_eta_=jets[i]->eta();
	    jet3_phi_=jets[i]->phi();
	  }
	  if (jets[i]->pt() > 30.0 && isInCentralGap){
	    ++n_jets_cjv_30_;
	  }
	  if ( ((jets[i]->eta()<2.4 && jets[i]->pt() > 20.0) ||
		(jets[i]->eta()>=2.4 && jets[i]->pt() > 30.0)) && 
	       isInCentralGap){
	    ++n_jets_cjv_20EB_30EE_;
	  }
	  if(jets[i]->pt()>30.0){
	    double thisjetmetnomudphi = fabs(ROOT::Math::VectorUtil::DeltaPhi(jets[i]->vector(),metnomuvec));
	    if(thisjetmetnomudphi<alljetsmetnomu_mindphi_)alljetsmetnomu_mindphi_=thisjetmetnomudphi;
	    double thisjetmetdphi = fabs(ROOT::Math::VectorUtil::DeltaPhi(jets[i]->vector(),metvec));
	    if(thisjetmetdphi<alljetsmet_mindphi_)alljetsmet_mindphi_=thisjetmetdphi;
	  }
	}
      }
      else{
	for (unsigned i = 0; i < jets.size(); ++i) {
	  if(jets[i]->id()==jet1->id()){
	    jet1_csv_=jets[i]->GetBDiscriminator("combinedSecondaryVertexBJetTags");
	  }
	  if(jets[i]->id()==jet2->id()){
	    jet2_csv_=jets[i]->GetBDiscriminator("combinedSecondaryVertexBJetTags");	  
	  }
	}
      }
      static unsigned processed = 0;
      //IF PASSES CUTS FILL TREE
      if(!ignoreLeptons_){
	if(!do_promptskim_){
	  if (metnomu_significance_ > 3.0 &&  dijet_deta_>3.6){
	    //if (dijet_M_>1000 &&  dijet_deta_>3.6 && metnomuons_>100 && jet1_pt_>50){//for prompt presel
	    outputTree_->Fill();
	    ++processed;
	  }
	}
	else{
	  if (passtrigger_==1&&dijet_deta_>3.6&&metnomuons_>90&&jet1_pt_>50){
	    //if (dijet_M_>1000 &&  dijet_deta_>3.6 && metnomuons_>100 && jet1_pt_>50){//for prompt presel
	    outputTree_->Fill();
	    ++processed;
	  }
	}
      }
      else{
	if(!do_promptskim_){
	  if (metnomu_significance_ > 3.0 &&  dijet_deta_>3.6 &&m_mumu_gen_>80&&m_mumu_gen_<100){
	    //if (dijet_M_>1000 &&  dijet_deta_>3.6 && metnomuons_>100 && jet1_pt_>50){//for prompt presel
	    outputTree_->Fill();
	    ++processed;
	  }
	}
	else{
	  if (passtrigger_==1&&dijet_deta_>3.6&&metnomuons_>90&&jet1_pt_>50){
	    //if (dijet_M_>1000 &&  dijet_deta_>3.6 && metnomuons_>100 && jet1_pt_>50){//for prompt presel
	    outputTree_->Fill();
	    ++processed;
	  }
	}
	
      }
      if (processed == 500) outputTree_->OptimizeBaskets();
    }

    return 0;
  */
  outputTree_->Fill();
} // pruneKids

//define this as a plug-in
DEFINE_FWK_MODULE(LightTreeProducer);
