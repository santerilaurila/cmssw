# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: RECO -s RAW2DIGI,L1Reco,RECO,PAT --runUnscheduled --data --scenario pp --conditions 80X_dataRun2_2016SeptRepro_v7 --era Run2_2016 --runUnscheduled --eventcontent RAWRECO,MINIAOD --datatier RAWRECO,MINIAOD --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run2_2016,RecoTracker/Configuration/customizeMinPtForHitRecoveryInGluedDet.customizeHitRecoveryInGluedDetOn,TauAnalysis/MCEmbeddingTools/customisers.customisoptions,TauAnalysis/MCEmbeddingTools/customisers.customiseSelecting_HplusFullyHadronic_Reselect --filein /store/data/Run2016H/DoubleMuon/RAW/v1/000/283/830/00000/E2A70B44-E49A-E611-B10B-02163E0135B9.root --fileout RAWskimmed.root -n 60 --no_exec --python_filename=selection.py
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECO',eras.Run2_2016)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(100)
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'root://cmsxrootd.fnal.gov//store/data/Run2016H/SingleMuon/RAW/v1/000/281/616/00000/92E42FCE-4283-E611-BABE-02163E012A07.root',
        'root://cmsxrootd.fnal.gov//store/data/Run2016H/SingleMuon/RAW/v1/000/281/616/00000/A236B0D8-2D83-E611-B2D4-02163E013620.root',
        'root://cmsxrootd.fnal.gov//store/data/Run2016H/SingleMuon/RAW/v1/000/281/616/00000/AA3A43D1-2C83-E611-A47D-FA163E09D143.root',
        'root://cmsxrootd.fnal.gov//store/data/Run2016H/SingleMuon/RAW/v1/000/281/616/00000/AAD906ED-3083-E611-B6DA-02163E01259E.root'),
#        'root://cmsxrootd.fnal.gov//store/data/Run2016G/SingleMuon/RAW/v1/000/278/820/00000/1A510AD3-7662-E611-836E-02163E0146AB.root',),
#        'root://cmsxrootd.fnal.gov//store/data/Run2016G/SingleMuon/RAW/v1/000/278/820/00000/1A63ACF1-7162-E611-A283-02163E0145AC.root',
#        'root://cmsxrootd.fnal.gov//store/data/Run2016G/SingleMuon/RAW/v1/000/278/820/00000/1A8787A2-7662-E611-A262-02163E011D9F.root',
#        'root://cmsxrootd.fnal.gov//store/data/Run2016G/SingleMuon/RAW/v1/000/278/820/00000/1A98CA2B-7662-E611-97A0-02163E01224A.root',
#        'root://cmsxrootd.fnal.gov//store/data/Run2016G/SingleMuon/RAW/v1/000/278/820/00000/1AAB8E45-7162-E611-9327-02163E011DE2.root',
#        'root://cmsxrootd.fnal.gov//store/data/Run2016G/SingleMuon/RAW/v1/000/278/820/00000/1AC15B65-7362-E611-AC89-02163E01443F.root'),
#    fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/data/Run2016H/DoubleMuon/RAW/v1/000/283/830/00000/E2A70B44-E49A-E611-B10B-02163E0135B9.root'),
#    fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/data/Run2016H/DoubleMuon/RAW/v1/000/283/830/00000/E2A70B44-E49A-E611-B10B-02163E0135B9.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('RECO nevts:60'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RAWRECOoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('RAWRECO'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('step1_output_RAWskimmed.root'),
    outputCommands = process.RAWRECOEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.MINIAODoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAOD'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('step1_output_RAWskimmed_inMINIAOD.root'),
    outputCommands = process.MINIAODEventContent.outputCommands,
    overrideInputFileSplitLevels = cms.untracked.bool(True)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
process.Flag_METFilters = cms.Path(process.metFilters)
process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseIsoFilter)
process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseFilter)
process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWRECOoutput_step = cms.EndPath(process.RAWRECOoutput)
process.MINIAODoutput_step = cms.EndPath(process.MINIAODoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.endjob_step,process.RAWRECOoutput_step,process.MINIAODoutput_step)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.RecoTLR
from Configuration.DataProcessing.RecoTLR import customisePostEra_Run2_2016 

#call to customisation function customisePostEra_Run2_2016 imported from Configuration.DataProcessing.RecoTLR
process = customisePostEra_Run2_2016(process)

# Automatic addition of the customisation function from RecoTracker.Configuration.customizeMinPtForHitRecoveryInGluedDet
from RecoTracker.Configuration.customizeMinPtForHitRecoveryInGluedDet import customizeHitRecoveryInGluedDetOn 

#call to customisation function customizeHitRecoveryInGluedDetOn imported from RecoTracker.Configuration.customizeMinPtForHitRecoveryInGluedDet
process = customizeHitRecoveryInGluedDetOn(process)

# Automatic addition of the customisation function from TauAnalysis.MCEmbeddingTools.customisers_singletau
from TauAnalysis.MCEmbeddingTools.customisers_singletau import customisoptions,customiseSelecting_Reselect 

#call to customisation function customisoptions imported from TauAnalysis.MCEmbeddingTools.customisers_singletau
process = customisoptions(process)

#call to customisation function customiseSelecting_Reselect imported from TauAnalysis.MCEmbeddingTools.customisers_singletau
process = customiseSelecting_Reselect(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)
process.load('Configuration.StandardSequences.PAT_cff')
from FWCore.ParameterSet.Utilities import cleanUnscheduled
process=cleanUnscheduled(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllData 

#call to customisation function miniAOD_customizeAllData imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllData(process)

# End of customisation functions
