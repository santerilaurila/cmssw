import FWCore.ParameterSet.Config as cms


externalLHEProducer = cms.EDProducer("SingleTauEmbeddingLHEProducer",
    src = cms.InputTag("selectedMuonForSingleMuonEmbedding","",""),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices","","SELECT"),
    switchToMuonEmbedding = cms.bool(False),
    rotate180 = cms.bool(False),
    mirror = cms.bool(False),
    studyFSRmode = cms.untracked.string("reco"),
    lhe_outputfilename = cms.untracked.string("step2a_output.lhe")
)

makeexternalLHEProducer = cms.Sequence( externalLHEProducer)
