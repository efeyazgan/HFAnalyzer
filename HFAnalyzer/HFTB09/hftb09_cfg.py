import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.tbUnpacker = cms.EDFilter("HcalTBObjectUnpacker",
   IncludeUnmatchedHits = cms.untracked.bool(False),
   HcalTDCFED = cms.untracked.int32(8),
   HcalQADCFED = cms.untracked.int32(8),
   HcalSlowDataFED = cms.untracked.int32(3),
   HcalTriggerFED = cms.untracked.int32(1),
#    ConfigurationFile = cms.untracked.string('configQADCTDC.txt')
   ConfigurationFile = cms.untracked.string('configQADCTDC_WCCalibration_Francesco.txt')
)

process.source = cms.Source("HcalTBSource",
#process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:/tmp/efe/HTB_105861.root'
	'file:/outdata2/HTB_106016.root'
    ),
    streams   = cms.untracked.vstring('HCAL_Trigger','HCAL_SlowData','HCAL_QADCTDC','HCAL_DCC021')                      
)

process.demo = cms.EDAnalyzer('HFTB09'
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('hf_tb09.root')
)


process.p = cms.Path(process.tbUnpacker
                     *process.demo
		)
