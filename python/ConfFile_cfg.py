import FWCore.ParameterSet.Config as cms

process = cms.Process("HiggsAnalysis")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False),
        SkipEvent = cms.untracked.vstring('ProductNotFound')
        )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

from UserCode.HGCanalysis.storeTools_cff import fillFromStore

files =  [f for f in fillFromStore('/afs/cern.ch/work/p/phansen/public/hgcal/CMSSW/hToGammaGammaPU_SLHC21/') ]

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(files),
        )
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.analysis = cms.EDAnalyzer('HiggsAnalyzer',
        genSource = cms.untracked.string("genParticles"),
        jetSource = cms.untracked.string("ak5PFJets"),
        linCorr = cms.untracked.double(1.0)
)

process.filter = cms.EDFilter("MCBarrelEndcapFilter",
       minInvMass = cms.untracked.double(127.0),
       maxInvMass = cms.untracked.double(133.0),
       verbose = cms.untracked.bool(True),
        )

process.TFileService = cms.Service("TFileService", fileName = cms.string('Higgs-analysis.root'))


#process.p = cms.Path(process.filter * process.analysis)
process.p = cms.Path(process.analysis)
