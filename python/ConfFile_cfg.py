import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/sps/cms/hbrun/CMSSW_6_2_0_patch1_L2seedingDev/src/files/oldGeometry_RECO.root'
    )
)

process.demo = cms.EDAnalyzer('L2seedsAnalyzer',
                              muonProducer 		= cms.VInputTag(cms.InputTag("muons")),
                              primaryVertexInputTag   	= cms.InputTag("offlinePrimaryVertices"),
                              outputFile = cms.string("muonSeedTree.root")
)


process.p = cms.Path(process.demo)
