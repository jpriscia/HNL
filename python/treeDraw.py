import FWCore.ParameterSet.Config as cms

process = cms.Process("AnalysisProc")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource", fileNames =  cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/5AE90EBB-54BF-E611-88AC-5065F381E271.root'))
#process.source = cms.Source("PoolSource", fileNames =  cms.untracked.vstring('file:/afs/cern.ch/user/j/jpriscia/heavyNeutrino_215.root') )
#process.source = cms.Source("PoolSource", fileNames =  cms.untracked.vstring('file:heavyNeutrino_235.root'))
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.maxEvents = cms.untracked.PSet(
input = cms.untracked.int32(100)
)

process.printDecay = cms.EDAnalyzer("ParticleDecayDrawer",
    src = cms.InputTag("prunedGenParticles"),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(True)
  )

process.printTree = cms.EDAnalyzer("ParticleListDrawer",
  maxEventsToPrint = cms.untracked.int32(100),
  printVertex = cms.untracked.bool(False),
  printOnlyHardInteraction = cms.untracked.bool(False), 
  src = cms.InputTag("prunedGenParticles")
)
process.printTree2 = cms.EDAnalyzer("ParticleTreeDrawer",
                                   src = cms.InputTag("prunedGenParticles"),#packedGenParticles                                                         
                                   printP4 = cms.untracked.bool(False),
                                   printPtEtaPhi = cms.untracked.bool(False),
                                   printVertex = cms.untracked.bool(False),
                                   printStatus = cms.untracked.bool(False),
                                   printIndex = cms.untracked.bool(False),
                                   #status = cms.untracked.vint32(1)
                                   )
process.p = cms.Path(process.printTree)
