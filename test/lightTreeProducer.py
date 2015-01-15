import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/cmst3/user/gpetrucc/miniAOD/v1/TT_Tune4C_13TeV-pythia8-tauola_PU_S14_PAT.root'
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
                                   fatjets = cms.InputTag("slimmedJetsAK8"),
                                   mets = cms.InputTag("slimmedMETs"),
                                   bits = cms.InputTag("TriggerResults","","HLT"),
                                   prescales = cms.InputTag("patTrigger"),
                                   objects = cms.InputTag("selectedPatTrigger"),
)


process.p = cms.Path(process.lightTree)
