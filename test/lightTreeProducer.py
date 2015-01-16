import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        '/store/mc/Phys14DR/VBF_HToInv_M-125_13TeV_powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/6421CDDA-3A70-E411-8520-00A0D1EEE5CC.root'
        '/store/mc/Phys14DR/DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0A9C5C2C-EE6F-E411-9A3E-E0CB4E29C50A.root'
    )
)

process.TFileService = cms.Service("TFileService",                                                                                                                    fileName = cms.string("hinvLightTree.root"),
                                   closeFileFast = cms.untracked.bool(True)                                                                                                                          
                                   )                                                                                                                                                                    
process.lightTree = cms.EDAnalyzer("LightTreeProducer",
                                   vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                   muons = cms.InputTag("slimmedMuons"),
                                   electrons = cms.InputTag("slimmedElectrons"),
                                   taus = cms.InputTag("slimmedTaus"),
                                   photons = cms.InputTag("slimmedPhotons"),
                                   jets = cms.InputTag("slimmedJets"),
                                   mets = cms.InputTag("slimmedMETs"),
                                   bits = cms.InputTag("TriggerResults","","HLT"),
                                   prescales = cms.InputTag("patTrigger"),
                                   objects = cms.InputTag("selectedPatTrigger"),
)


process.p = cms.Path(process.lightTree)
