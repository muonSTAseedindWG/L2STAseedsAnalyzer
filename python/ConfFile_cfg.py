import FWCore.ParameterSet.Config as cms

process = cms.Process("ALZ")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = 'PRE_ST62_V8::All'


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
                                      #'file:/sps/cms/hbrun/CMSSW_6_2_0_patch1_L2seedingDev/src/files/oldGeometry_RECO.root'
                                      'file:/sps/cms/hbrun/CMSSW_6_2_0_patch1_L2seedingDev/src/files/NewGeometry_RECO.root' #attention Pt 100 
    )
)

process.load("SimMuon.MCTruth.MuonAssociatorByHits_cfi")
process.muonAssociatorByHits.tracksTag = cms.InputTag("standAloneMuons")
process.muonAssociatorByHits.UseTracker = cms.bool(False)
process.muonAssociatorByHits.UseMuon = cms.bool(True)
process.muonAssociatorByHits.PurityCut_muon = cms.double(0.01)
process.muonAssociatorByHits.EfficiencyCut_muon = cms.double(0.01)

import SimMuon.MCTruth.MuonAssociatorByHits_cfi

process.muonAssociatorByHitsL2seeds = process.muonAssociatorByHits.clone() #don't use it, it is not working ;)
process.muonAssociatorByHitsL2seeds.tracksTag = cms.InputTag("ancientMuonSeed")
process.muonAssociatorByHitsL2seeds.UseTracker = cms.bool(False)
process.muonAssociatorByHitsL2seeds.UseMuon = cms.bool(True)
process.muonAssociatorByHitsL2seeds.PurityCut_muon = cms.double(0.01)
process.muonAssociatorByHitsL2seeds.EfficiencyCut_muon = cms.double(0.01)

process.runL2seed = cms.EDAnalyzer('L2seedsAnalyzer',
                              isMC                    = cms.bool(True),
                              muonProducer 		= cms.VInputTag(cms.InputTag("muons")),
                              primaryVertexInputTag   	= cms.InputTag("offlinePrimaryVertices"),
                              StandAloneTrackCollectionLabel = cms.untracked.string("standAloneMuons"),
                              trackingParticlesCollection = cms.InputTag("mix","MergedTrackTruth"),
                              standAloneAssociator = cms.InputTag("muonAssociatorByHits"),
                              L2seedsCollection = cms.InputTag("ancientMuonSeed"),
                              MuonRecHitBuilder = cms.string("MuonRecHitBuilder"),
                              outputFile = cms.string("muonSeedTree.root")
)


process.p = cms.Path(process.muonAssociatorByHits+process.runL2seed)





