import sys, os
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from pdb import set_trace


process = cms.Process("HNLSecondaryVertex")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")


process.options   = cms.untracked.PSet(
      wantSummary = cms.untracked.bool(True),
      SkipEvent = cms.untracked.vstring('ProductNotFound'),
)

process.MessageLogger.cerr.FwkReport.reportEvery = 50
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

print '------------'
print 'inputFiles'
print '------------'

process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

#put input file


#lines = open("files_4GeV.txt").read().splitlines()
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
#process.source = cms.Source("PoolSource",
#                             fileNames = cms.untracked.vstring(lines[options.jobNumber]),
#                            secondaryFileNames = cms.untracked.vstring()
#                            )

process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck")

readFiles.extend(['root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-10_V-0.00244948974278_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_189.root'])


process.inclusiveVertexFinder  = cms.EDProducer("InclusiveVertexFinder",  #????
       beamSpot = cms.InputTag("offlineBeamSpot"),
       primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
       tracks = cms.InputTag("unpackedTracksAndVertices"),
       minHits = cms.uint32(6),#8
       maximumLongitudinalImpactParameter = cms.double(100000000000.),#0.3
       minPt = cms.double(0.4), #was 0.8 (does it make sense?)
       maxNTracks = cms.uint32(40),#30

       clusterizer = cms.PSet(
           seedMax3DIPSignificance = cms.double(9999.),  #between PV and track
           seedMax3DIPValue = cms.double(9999.),         #between PV and track 
           seedMin3DIPSignificance = cms.double(0),    #between PV and track 1.2
           seedMin3DIPValue = cms.double(0),         #between PV and track 0.005 
           clusterMaxDistance = cms.double(0.4), #500um  #among two tracks of the cluster 0.05   
           clusterMaxSignificance = cms.double(4.5), #4.5 sigma
           distanceRatio = cms.double(20), # was cluster scale = 1 / density factor =0.05  (20)
           clusterMinAngleCosine = cms.double(0.5), # only forward decays
       ),

       vertexMinAngleCosine = cms.double(0.95), # scalar prod direction of tracks and flight dir  0.95
       vertexMinDLen2DSig = cms.double(2.5), #2.5 sigma  provo ad alzarlo
       vertexMinDLenSig = cms.double(0.5), #0.5 sigma
       fitterSigmacut =  cms.double(3),
       fitterTini = cms.double(256),
       fitterRatio = cms.double(0.25),
       useDirectVertexFitter = cms.bool(True),
       useVertexReco  = cms.bool(True),
       vertexReco = cms.PSet(
               finder = cms.string('avr'),
               primcut = cms.double(1.0),
               seccut = cms.double(3),
               smoothing = cms.bool(True)
       )


)

#0.04 dr between tracks in ivf and PV

process.HNLSecondaryVertex = cms.Path(process.unpackedTracksAndVertices+process.inclusiveVertexFinder) 

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')
   
########################################################################


#make the pool output
process.Out = cms.OutputModule("PoolOutputModule",
     outputCommands = cms.untracked.vstring(
         "keep *",
         ),
    fileName = cms.untracked.string('/eos/user/j/jpriscia/HNL_IVF/ivf10GeV_47.root'),
    SelectEvents = cms.untracked.PSet(
       SelectEvents = cms.vstring('*')
    ),
)

########################################################################

process.tsk = cms.Task()
for mod in process.producers_().itervalues():
    process.tsk.add(mod)
for mod in process.filters_().itervalues():
    process.tsk.add(mod)

#schedule the sequence
process.endPath1 = cms.EndPath(process.Out)
process.schedule = cms.Schedule(process.HNLSecondaryVertex, process.endPath1)
