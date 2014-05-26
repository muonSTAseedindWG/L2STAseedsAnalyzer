# Auto generated configuration file
# using: 
# Revision: 1.20 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step2 --filein file:step1.root --fileout file:/tmp/hbrun/MUO-Fall13dr-00013.root --mc --eventcontent RECOSIM --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1,Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM-RECO --conditions POSTLS162_V2::All --step RAW2DIGI,L1Reco,RECO --magField 38T_PostLS1 --geometry Extended2015 --python_filename MUO-Fall13dr-00013_2_cfg.py --no_exec -n -1
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
        '/store/user/hbrun/seedTreesTSGsample/subsetByType/events_STAZ2_fail.root'

    )
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('step2 nevts:-1'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('file:/tmp/hbrun/MUO-Fall13dr-00013.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    )
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS162_V2::All', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
#process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

process.load("SimMuon.MCTruth.MuonAssociatorByHitsESProducer_NoSimHits_cfi") 
process.runL2seed = cms.EDAnalyzer('L2seedsAnalyzer',
                              isMC                      = cms.bool(True),
                              selectJpsiOnly            = cms.bool(False),
                              muonProducer 		= cms.VInputTag(cms.InputTag("muons")),
                              primaryVertexInputTag   	= cms.InputTag("offlinePrimaryVertices"),
                              StandAloneTrackCollectionLabel = cms.untracked.string("standAloneMuons"),
                              trackingParticlesCollection = cms.InputTag("mix","MergedTrackTruth"),
                              standAloneAssociator = cms.InputTag("muonAssociatorByHits"),
                              L2seedsCollection = cms.InputTag("ancientMuonSeed"),
                              L2seedTrackCollection = cms.InputTag("seedsInTracks"),
                              L2associator = cms.InputTag("muonAssociatorByHitsL2seeds"),
                              MuonRecHitBuilder = cms.string("MuonRecHitBuilder"),
			      associatorLabel = cms.string("muonAssociatorByHits_NoSimHits"),
                              cscSegmentCollection = cms.InputTag("cscSegments"),
                              dtSegmentCollection = cms.InputTag("dt4DSegments"),
                              outputFile = cms.string("muonSeedTree.root")
)



process.load("SimMuon.MCTruth.MuonTrackProducer_cfi")

process.seedsInTracks= cms.EDProducer('SeedToTrackProducer',
                                         L2seedsCollection = cms.InputTag("ancientMuonSeed")
                                         )


import SimGeneral.MixingModule.trackingTruthProducer_cfi
process.mergedtruthNoSimHits = process.trackingParticles.clone(
                                                            simHitCollections = cms.PSet(
                                                                                         muon = cms.VInputTag(),
                                                                                         tracker = cms.VInputTag(),
                                                                                         pixel = cms.VInputTag()
                                                                                         )
                                                            )

process.mix.digitizers = cms.PSet( mergedtruth = process.mergedtruthNoSimHits )
process.mix.mixObjects = cms.PSet()
del process.simCastorDigis
del process.simEcalUnsuppressedDigis
del process.simHcalUnsuppressedDigis
del process.simSiPixelDigis
del process.simSiStripDigis


process.runMix = cms.Path(process.mix)
process.p = cms.Path(process.seedsInTracks)
process.pAna = cms.EndPath(process.runL2seed)
 

# Schedule definition
process.schedule = cms.Schedule(process.runMix,process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.p, process.pAna)#, process.RECOSIMoutput_step)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# End of customisation functions
