import FWCore.ParameterSet.Config as cms

vertexMerger = cms.EDProducer("VertexMerger",
       primaryVertices = cms.InputTag("offlinePrimaryVertices"),
       secondaryVertices = cms.InputTag("inclusiveVertexFinder"),
       maxFraction = cms.double(0.7),
       minSignificance = cms.double(2),
       doMerging = cms.bool(True)
)


