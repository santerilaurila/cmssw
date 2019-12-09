# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: LHEprodandCLEAN --filein file:RAWskimmed.root --fileout file:lhe_and_cleaned.root --runUnscheduled --data --era Run2_2016 --scenario pp --conditions 80X_dataRun2_2016SeptRepro_v7 --eventcontent RAWRECO,MINIAOD --datatier RAWRECO,MINIAOD --step RAW2DIGI,RECO,PAT --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run2_2016,RecoTracker/Configuration/customizeMinPtForHitRecoveryInGluedDet.customizeHitRecoveryInGluedDetOn,TauAnalysis/MCEmbeddingTools/customisers.customisoptions,TauAnalysis/MCEmbeddingTools/customisers.customiseLHEandCleaning_Reselect --no_exec -n -1 --python_filename lheprodandcleaning.py
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
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:step1_output_RAWskimmed.root'),
#    fileNames = cms.untracked.vstring('file:/eos/user/s/slaurila/data/embedding_plot_inputs/step1_output_RAWskimmed_3688.root'),
#    fileNames = cms.untracked.vstring('file:multicrab_step1_selection_TEST/multicrab_step1_selection_v8026p1_20171220T1640/crab_SingleMuon_Run2016H_v1_281010_281202/results/step1_output_RAWskimmed_292.root',
#                                      'file:multicrab_step1_selection_TEST/multicrab_step1_selection_v8026p1_20171220T1640/crab_SingleMuon_Run2016H_v1_281010_281202/results/step1_output_RAWskimmed_310.root',    
#                                      'file:multicrab_step1_selection_TEST/multicrab_step1_selection_v8026p1_20171220T1640/crab_SingleMuon_Run2016H_v1_281010_281202/results/step1_output_RAWskimmed_382.root',    
#                                      'file:multicrab_step1_selection_TEST/multicrab_step1_selection_v8026p1_20171220T1640/crab_SingleMuon_Run2016H_v1_281010_281202/results/step1_output_RAWskimmed_349.root',    
#                                      'file:multicrab_step1_selection_TEST/multicrab_step1_selection_v8026p1_20171220T1640/crab_SingleMuon_Run2016H_v1_281010_281202/results/step1_output_RAWskimmed_391.root',    
#                                      'file:multicrab_step1_selection_TEST/multicrab_step1_selection_v8026p1_20171220T1640/crab_SingleMuon_Run2016H_v1_281010_281202/results/step1_output_RAWskimmed_385.root',    
#                                      'file:multicrab_step1_selection_TEST/multicrab_step1_selection_v8026p1_20171220T1640/crab_SingleMuon_Run2016H_v1_281010_281202/results/step1_output_RAWskimmed_404.root',    
#                                      'file:multicrab_step1_selection_TEST/multicrab_step1_selection_v8026p1_20171220T1640/crab_SingleMuon_Run2016H_v1_281010_281202/results/step1_output_RAWskimmed_351.root',    
#                                      'file:multicrab_step1_selection_TEST/multicrab_step1_selection_v8026p1_20171220T1640/crab_SingleMuon_Run2016H_v1_281010_281202/results/step1_output_RAWskimmed_387.root',    
#                                      'file:multicrab_step1_selection_TEST/multicrab_step1_selection_v8026p1_20171220T1640/crab_SingleMuon_Run2016H_v1_281010_281202/results/step1_output_RAWskimmed_415.root',    
#                                      'file:multicrab_step1_selection_TEST/multicrab_step1_selection_v8026p1_20171220T1640/crab_SingleMuon_Run2016H_v1_281010_281202/results/step1_output_RAWskimmed_396.root',    
#                                      'file:multicrab_step1_selection_TEST/multicrab_step1_selection_v8026p1_20171220T1640/crab_SingleMuon_Run2016H_v1_281010_281202/results/step1_output_RAWskimmed_478.root',    
#                                      'file:multicrab_step1_selection_TEST/multicrab_step1_selection_v8026p1_20171220T1640/crab_SingleMuon_Run2016H_v1_281010_281202/results/step1_output_RAWskimmed_386.root',    
#                                      'file:multicrab_step1_selection_TEST/multicrab_step1_selection_v8026p1_20171220T1640/crab_SingleMuon_Run2016H_v1_281010_281202/results/step1_output_RAWskimmed_422.root',    
#                                      'file:multicrab_step1_selection_TEST/multicrab_step1_selection_v8026p1_20171220T1640/crab_SingleMuon_Run2016H_v1_281010_281202/results/step1_output_RAWskimmed_473.root',    
#                                      'file:multicrab_step1_selection_TEST/multicrab_step1_selection_v8026p1_20171220T1640/crab_SingleMuon_Run2016H_v1_281010_281202/results/step1_output_RAWskimmed_417.root',    
#                                      'file:multicrab_step1_selection_TEST/multicrab_step1_selection_v8026p1_20171220T1640/crab_SingleMuon_Run2016H_v1_281010_281202/results/step1_output_RAWskimmed_286.root',    
#                                      'file:multicrab_step1_selection_TEST/multicrab_step1_selection_v8026p1_20171220T1640/crab_SingleMuon_Run2016H_v1_281010_281202/results/step1_output_RAWskimmed_423.root',    
#                                      'file:multicrab_step1_selection_TEST/multicrab_step1_selection_v8026p1_20171220T1640/crab_SingleMuon_Run2016H_v1_281010_281202/results/step1_output_RAWskimmed_406.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True)

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('LHEprodandCLEAN nevts:-1'),
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
    fileName = cms.untracked.string('file:step2_output_lhe_and_cleaned.root'),
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
    fileName = cms.untracked.string('file:step2_output_lhe_and_cleaned_inMINIAOD.root'),
    outputCommands = process.MINIAODEventContent.outputCommands,
    overrideInputFileSplitLevels = cms.untracked.bool(True)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
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
process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.endjob_step,process.RAWRECOoutput_step,process.MINIAODoutput_step)

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
from TauAnalysis.MCEmbeddingTools.customisers_singletau import customisoptions,customiseLHEandCleaning_Reselect 

#call to customisation function customisoptions imported from TauAnalysis.MCEmbeddingTools.customisers_singletau
process = customisoptions(process)

#call to customisation function customiseLHEandCleaning_Reselect imported from TauAnalysis.MCEmbeddingTools.customisers_singletau
process = customiseLHEandCleaning_Reselect(process)

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
