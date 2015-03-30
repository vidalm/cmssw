import FWCore.ParameterSet.Config as cms

vertexMerger = cms.EDProducer("VertexMerger",
       secondaryVertices = cms.InputTag("inclusiveVertexFinder"),
       doMerging = cms.bool(True),
       maxFraction = cms.double(0.7),
       minSignificance = cms.double(2)
#       minSignificance = cms.double(0.5)
#        minSignificance = cms.double(5)
#       minSignificance = cms.double(0) 
)


