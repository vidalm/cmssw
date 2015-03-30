#Basic example configration file to run the b-tagging validation sequence.
import FWCore.ParameterSet.Config as cms
process = cms.Process("validation")

"""
start customization
"""

#Enter here the Global tags
#tag =  'POSTLS172_V3::All'
tag = 'MCRUN2_73_V7::All'
#Data or MC?
runOnMC    = True
#Flavour plots for MC: "all" = plots for all jets ; "dusg" = plots for d, u, s, dus, g independently ; not mandatory and any combinations are possible 
#b, c, light (dusg), non-identified (NI), PU jets plots are always produced
flavPlots = "allbcldusg"
#Check if jets originate from PU? option recommended (only for MC)
PUid = True
#List of taggers and taginfo to be considered (see example in: DQMOffline/RecoB/python/bTagCommon_cff.py)
from DQMOffline.RecoB.bTagCommon_cff import *
tagConfig = cms.VPSet(
        cms.PSet(
             bTagSimpleSVAnalysisBlock,
             #label = cms.InputTag("pfSimpleSecondaryVertexHighEffBJetTags"),
             label = cms.InputTag("simpleInclusiveSecondaryVertexHighEffBJetTags"),
             folder = cms.string("SISVHE")
        ),
)


# Rerunning the iclusive vertexing 

process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")
process.load("RecoVertex.AdaptiveVertexFinder.trackVertexArbitrator_cfi")
process.load("RecoVertex.AdaptiveVertexFinder.vertexMerger_cfi")
process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexFinder_cfi")

process.MyInclusiveVertex = cms.Sequence(process.inclusiveVertexFinder * process.vertexMerger * process.trackVertexArbitrator * process.inclusiveSecondaryVertices)

# My tagger not default in reco
# Running the simpleInclusiveSecondaryVertexHighEffBJetTags

process.load("RecoBTag.SecondaryVertex.simpleInclusiveSecondaryVertexHighEffBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.inclusiveSecondaryVertexFinderFilteredTagInfos_cfi")
process.load("RecoBTag.SecondaryVertex.bToCharmDecayVertexMerger_cfi")
process.load("RecoBTag.SecondaryVertex.secondaryVertex_cff")
#process.load("RecoBTag.SecondaryVertex.bVertexFilter_cfi")
process.MyTagger = cms.Sequence(process.MyInclusiveVertex * process.inclusiveSecondaryVerticesFiltered * process.bToCharmDecayVertexMerged * process.inclusiveSecondaryVertexFinderFilteredTagInfos * process.simpleInclusiveSecondaryVertexHighEffBJetTags)

"""
end customization
"""

###prints###
print "is it MC ? : ", runOnMC
print "Global Tag : ", tag
############

process.load("DQMServices.Components.DQMEnvironment_cfi")
process.load("DQMServices.Core.DQM_cfg")

#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

if runOnMC:
    #Simple 

    #for MC jet flavour
    process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")
    process.AK4byRef.jets = cms.InputTag("ak4PFJetsCHS")
    process.flavourSeq = cms.Sequence(
        process.myPartons *
        process.AK4Flavour
    )
    #Validation sequence
    process.load("Validation.RecoB.bTagAnalysis_cfi")
    process.bTagValidation.jetMCSrc = 'AK4byValAlgo'
    process.bTagValidation.tagConfig = tagConfig
    process.bTagHarvestMC.tagConfig = tagConfig
    process.bTagValidation.flavPlots = flavPlots
    process.bTagHarvestMC.flavPlots = flavPlots
    process.bTagValidation.doPUid = cms.bool(PUid)
    process.ak4GenJetsForPUid = cms.EDFilter("GenJetSelector",
                                             src = cms.InputTag("ak4GenJets"),
                                             cut = cms.string('pt > 8.'),
                                             filter = cms.bool(False)
                                             )
    process.load("PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi")
    process.patJetGenJetMatch.matched = cms.InputTag("ak4GenJetsForPUid")
    process.patJetGenJetMatch.maxDeltaR = cms.double(0.25)
    process.patJetGenJetMatch.resolveAmbiguities = cms.bool(True)
else :
    process.load("DQMOffline.RecoB.bTagAnalysisData_cfi")
    process.bTagAnalysis.tagConfig = tagConfig
    process.bTagHarvest.tagConfig = tagConfig

# load the full reconstraction configuration, to make sure we're getting all needed dependencies
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")

process.GlobalTag.globaltag = tag

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

if runOnMC:
    process.dqmSeq = cms.Sequence(process.ak4GenJetsForPUid * process.patJetGenJetMatch * process.flavourSeq * process.bTagValidation * process.bTagHarvestMC * process.dqmSaver)
else:
    process.dqmSeq = cms.Sequence(process.bTagAnalysis * process.bTagHarvest * process.dqmSaver)

process.plots = cms.Path(process.MyTagger * process.dqmSeq)
    
process.dqmEnv.subSystemFolder = 'BTAG'
process.dqmSaver.producer = 'DQM'
process.dqmSaver.workflow = '/POG/BTAG/BJET'
process.dqmSaver.convention = 'Offline'
process.dqmSaver.saveByRun = cms.untracked.int32(-1)
process.dqmSaver.saveAtJobEnd =cms.untracked.bool(True) 
process.dqmSaver.forceRunNumber = cms.untracked.int32(1)
process.PoolSource.fileNames = [
    '/store/relval/CMSSW_7_4_0_pre5/RelValTTbar_13/GEN-SIM-RECO/MCRUN2_73_V7-v1/00000/48F1F053-EE9D-E411-B4AB-0025905938D4.root',
    '/store/relval/CMSSW_7_4_0_pre5/RelValTTbar_13/GEN-SIM-RECO/MCRUN2_73_V7-v1/00000/A6FF14A5-F39D-E411-B7EE-0025905964C0.root',
    '/store/relval/CMSSW_7_4_0_pre5/RelValTTbar_13/GEN-SIM-RECO/MCRUN2_73_V7-v1/00000/B6CCCDA5-F39D-E411-889E-0025905A60E0.root' 
]

