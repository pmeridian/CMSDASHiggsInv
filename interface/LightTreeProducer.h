// system include files
#include <memory>
#include <iostream>
#include <string>
#include <vector>

// user include files
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <map>
#include <set>

#include "TTree.h"

using namespace edm;
using namespace std;
using namespace reco;


class LightTreeProducer : public edm::EDAnalyzer 
{
   public:
      explicit LightTreeProducer(const edm::ParameterSet&);
      ~LightTreeProducer();
      int hltSkim_;

   private:
      
      virtual void beginJob();
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;


      template<typename T>
	struct RefGreaterByPt {
	  typedef T first_argument_type;
	  typedef T second_argument_type;
	  bool operator()( const T*  t1, const T* t2 ) const {
	    return t1->pt() > t2->pt();
	  }
	};


      bool MinDRToCollection(reco::Candidate const* cand,std::vector<const reco::Candidate*>& coll, double cut);
      bool pu_id_mva_loose(const pat::Jet& j);
      // ----------member data ---------------------------
      edm::Service<TFileService> fs_;      

      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      edm::EDGetTokenT<pat::TauCollection> tauToken_;
      edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
      edm::EDGetTokenT<pat::JetCollection> jetToken_;
      edm::EDGetTokenT<pat::METCollection> metToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
      edm::EDGetTokenT<std::vector<l1extra::L1EtMissParticle> > l1MetToken_;

      TTree *outputTree_;
      
      unsigned run_;
      unsigned lumi_;
      unsigned event_;
      double weight_nolep_;
      double total_weight_lepveto_;
      double total_weight_leptight_;
      double puweight_up_scale_;
      double puweight_down_scale_;
      double jet1_pt_;
      double jet2_pt_;
      double jet3_pt_;
      double jet1_E_;
      double jet2_E_;
      double jet3_E_;
      double jet1_eta_;
      double jet2_eta_;
      double jet3_eta_;
      double jet1_phi_;
      double jet2_phi_;
      double jet3_phi_;
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
      double metnomu_x_;
      double metnomu_y_;
      double met_significance_;
      double metnomu_significance_;
      double sumet_;
      double ht_;
      double ht30_;
      double mht_;
      double sqrt_ht_;
      double unclustered_et_;
      double jet1met_dphi_;
      double jet2met_dphi_;
      double jet1metnomu_dphi_;
      double jet2metnomu_dphi_;
      double jetmet_mindphi_;
      double jetmetnomu_mindphi_;
      double alljetsmet_mindphi_;
      double alljetsmetnomu_mindphi_;
      double jetunclet_mindphi_;
      double metunclet_dphi_;
      double metnomuunclet_dphi_;
      double dijetmet_scalarSum_pt_;
      double dijetmet_vectorialSum_pt_;
      double dijetmet_ptfraction_;
      double jet1met_scalarprod_;
      double jet2met_scalarprod_;
      double dijetmetnomu_scalarSum_pt_;
      double dijetmetnomu_vectorialSum_pt_;
      double dijetmetnomu_ptfraction_;
      double jet1metnomu_scalarprod_;
      double jet2metnomu_scalarprod_;
      double cjvjetpt_;
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
      int njets_;
      
      double m_mumu_;
      double m_mumu_gen_;
      double mu1_pt_;
      double mu1_eta_;
      double mu1_phi_;
      double mu2_pt_;
      double mu2_eta_;
      double mu2_phi_;
      double ele1_pt_;
      double ele1_eta_;
      double ele1_phi_;
      double tau1_pt_;
      double tau1_eta_;
      double tau1_phi_;
      double lep_mt_;
      int n_vertices_;
};

