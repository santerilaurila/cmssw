import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.PAT_cff import *
from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import patMuons
from HLTrigger.HLTfilters.triggerResultsFilter_cfi import *


## Trigger requirements #FIXME
singleMuonHLTTrigger = cms.EDFilter("TriggerResultsFilter",
    hltResults = cms.InputTag("TriggerResults","","HLT"),
    l1tResults = cms.InputTag(""),
    throw = cms.bool(False),
#    triggerConditions = cms.vstring("HLT_IsoMu27_*") # unprescaled, simple approach: ~37% of SingleMu events pass
    triggerConditions = cms.vstring("HLT_IsoMu24_* OR HLT_IsoTkMu24_v* OR HLT_Mu50_v*") # unprescaled, following HIG-14-023 leptonic analysis, ~48% of SingleMu events pass
#    triggerConditions = cms.vstring("HLT_Mu27_*") # avoids biases but was unprescaled at some point
#    triggerConditions = cms.vstring("HLT_Mu50*") # lowest unprescaled non-isolated single muon trigger
#    triggerConditions = cms.vstring("HLT_Mu45_eta2p1_*") # was prescaled at some point, might lead to bias around eta=2.1
#    triggerConditions = cms.vstring("HLT_Mu27_* OR HLT_Mu45_eta2p1*")
#    triggerConditions = cms.vstring("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v* OR HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*") # for ditau embedding
)


## Muon and jet selection
skimForSingleMuonEmbedding = cms.EDFilter("SingleMuForEmbeddingSelector",
    VertexCollection = cms.InputTag("offlinePrimaryVertices"),
    MuonCollection = cms.InputTag("slimmedMuons"),
    MuonPtCut = cms.double(45.0),
    MuonEtaCut = cms.double(2.5),
    MuonVetoPtCut = cms.double(20.0),
    MuonVetoEtaCut = cms.double(2.5),
    JetCollection = cms.InputTag("slimmedJets"),
    JetPtCut = cms.double(25.0),
    JetEtaCut = cms.double(5.0),
    NJetsCut = cms.int32(2),
    METCollection = cms.InputTag("slimmedMETs"),
    METCut = cms.double(0.0)
)

selectedMuonForSingleMuonEmbedding = cms.EDProducer("SingleMuForEmbeddingProducer",
    VertexCollection = cms.InputTag("offlinePrimaryVertices"),
    MuonCollection = cms.InputTag("slimmedMuons"),
    MuonPtCut = cms.double(45.0),
    MuonEtaCut = cms.double(2.5),
)    

# Define the complete sequence
makePatMuonsForSingleTauEmbedding = cms.Sequence(
    singleMuonHLTTrigger
    + skimForSingleMuonEmbedding
    + selectedMuonForSingleMuonEmbedding
)
