import FWCore.ParameterSet.Config as cms

process = cms.Process("AnalysisProc")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource", fileNames =  cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_119.root'))
#process.source = cms.Source("PoolSource", fileNames =  cms.untracked.vstring('file:/afs/cern.ch/user/j/jpriscia/heavyNeutrino_215.root') )
#process.source = cms.Source("PoolSource", fileNames =  cms.untracked.vstring('file:heavyNeutrino_235.root'))
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.maxEvents = cms.untracked.PSet(
input = cms.untracked.int32(5)
)

process.printDecay = cms.EDAnalyzer("ParticleDecayDrawer",
    src = cms.InputTag("genParticles"),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(True)
  )

process.printTree = cms.EDAnalyzer("ParticleListDrawer",
  maxEventsToPrint = cms.untracked.int32(1),
  printVertex = cms.untracked.bool(False),
  #printOnlyHardInteraction = cms.untracked.bool(False), 
  src = cms.InputTag("genParticles")
)
process.printTree2 = cms.EDAnalyzer("ParticleTreeDrawer",
                                   src = cms.InputTag("genParticles"),                                                                 
                                   printP4 = cms.untracked.bool(True),
                                   printPtEtaPhi = cms.untracked.bool(True),
                                   printVertex = cms.untracked.bool(True),
                                   printStatus = cms.untracked.bool(True),
                                   printIndex = cms.untracked.bool(True),
                                   #status = cms.untracked.vint32( 3 )
                                   )
process.p = cms.Path(process.printTree2)
