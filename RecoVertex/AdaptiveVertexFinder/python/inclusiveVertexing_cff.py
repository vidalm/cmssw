import FWCore.ParameterSet.Config as cms

from RecoVertex.AdaptiveVertexFinder.inclusiveVertexFinder_cfi import *
from RecoVertex.AdaptiveVertexFinder.vertexMerger_cfi import *
from RecoVertex.AdaptiveVertexFinder.trackVertexArbitrator_cfi import *

inclusiveSecondaryVertices = vertexMerger.clone()
inclusiveSecondaryVertices.secondaryVertices = cms.InputTag("trackVertexArbitrator")
inclusiveSecondaryVertices.maxFraction = 0.7
inclusiveSecondaryVertices.minSignificance = 2.
inclusiveSecondaryVertices.doMerging= True

inclusiveVertexing = cms.Sequence(inclusiveVertexFinder*vertexMerger*trackVertexArbitrator*inclusiveSecondaryVertices)

from RecoVertex.AdaptiveVertexFinder.inclusiveCandidateVertexFinder_cfi import *
from RecoVertex.AdaptiveVertexFinder.candidateVertexMerger_cfi import *
from RecoVertex.AdaptiveVertexFinder.candidateVertexArbitrator_cfi import *

inclusiveCandidateSecondaryVertices = candidateVertexMerger.clone()
inclusiveCandidateSecondaryVertices.secondaryVertices = cms.InputTag("candidateVertexArbitrator")
inclusiveCandidateSecondaryVertices.maxFraction = 0.7
inclusiveCandidateSecondaryVertices.minSignificance = 2.
inclusiveCandidateSecondaryVertices.doMerging= True

inclusiveVertexing = cms.Sequence(inclusiveVertexFinder*vertexMerger*trackVertexArbitrator*inclusiveSecondaryVertices)
inclusiveCandidateVertexing = cms.Sequence(inclusiveCandidateVertexFinder*candidateVertexMerger*candidateVertexArbitrator*inclusiveCandidateSecondaryVertices)
