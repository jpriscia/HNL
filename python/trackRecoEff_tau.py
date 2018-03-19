import ROOT
import rootpy
import rootpy.plotting as plt
from rootpy import stl
from rootpy.io import root_open
from rootpy.tree import Tree, TreeModel, FloatCol, IntCol
import pprint
from DataFormats.FWLite import Events, Handle
from ROOT import gSystem
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaR2, bestMatch
from math import floor

from pdb import set_trace


def bestMatchedParticles( object, matchCollection):
    '''Return the best match to object in matchCollection, which is the closest object in delta R'''
    deltaR2Min = float('+inf')
    bm = None
    matchList=[]
    for match in matchCollection:
        dR2 = deltaR2( object.eta(), object.phi(),
                       match.eta(), match.phi() )
        if dR2<0.3:
            matchList.append((dR2,match))
    matchList.sort(key=lambda tup: tup[0])
    return matchList

def getDaughters(p):
    childrens = []
    for i in range(0, p.numberOfDaughters()):
        child = p.daughter(i)
        childrens.append(child)
    return childrens


def isAncestor(a, p):
    if a == p :
        return True
    for i in xrange(0,p.numberOfMothers()):
        if isAncestor(a,p.mother(i)):
            return True
        return False

#lines = open("fileList_mu2GeV.txt").read().splitlines()

#events = Events(lines)
########### for test uncomment this ############
events = Events(['root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/163E57C9-7ABE-E611-A73A-0025905B857E.root',
                 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/FA69C320-28C0-E611-8F7A-02163E013A88.root',
                 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/5285C43E-41BF-E611-8034-001C23C0F1F9.root'
                 ])
#event = events.__iter__().next()

handleReco, labelReco  = Handle('vector<pat::PackedCandidate>'),'packedPFCandidates'
handlePruned, labelPruned  = Handle ('vector<reco::GenParticle>'), 'prunedGenParticles'

handleJet, labelJet = Handle('vector<pat::Jet>'), 'slimmedJets'
handlePacked, labelPacked  = Handle ("std::vector<pat::PackedGenParticle> "), "packedGenParticles"

#######OutputTTree###############
outfile = root_open('test.root', 'w')

class GenTree(TreeModel):

        p_tauGen = FloatCol()
        pt_tauGen = FloatCol()
        phi_tauGen = FloatCol()
        eta_tauGen = FloatCol()
        aperture_Gen = FloatCol()
        p_Gen =  stl.vector(float)
        pt_Gen = stl.vector(float)
        px_Gen = stl.vector(float)
        py_Gen = stl.vector(float)
        pz_Gen = stl.vector(float)
        eta_Gen  = stl.vector(float)
        phi_Gen  = stl.vector(float)
        theta_Gen  = stl.vector(float)
        mass_Gen  = stl.vector(float)
        e_Gen  = stl.vector(float)

        p_Reco =  stl.vector(float)
        pt_Reco = stl.vector(float)
        px_Reco  = stl.vector(float)
        py_Reco  = stl.vector(float)
        pz_Reco  = stl.vector(float)

        vx_Reco  = stl.vector(float)
        vy_Reco  = stl.vector(float)
        vz_Reco  = stl.vector(float)

        eta_Reco  = stl.vector(float)
        phi_Reco  = stl.vector(float)
        theta_Reco  = stl.vector(float)

        idx = stl.vector(int)
        pdgId = stl.vector(int)
        dr = stl.vector(float)
        dx_Gen = stl.vector(float)
        dy_Gen = stl.vector(float)
        dz_Gen = stl.vector(float)
        # 'dxy', 'dxyError', 'dz', 'dzError'
        dxy_RecoPV = stl.vector(float)
        dxyErr_RecoPV = stl.vector(float)
        dz_RecoPV = stl.vector(float)
        dzErr_RecoPV = stl.vector(float)
        kaonCharged = IntCol()
        pionCharged = IntCol()
        kaonNeutrals = IntCol()
        pionNeutrals = IntCol()
        muons  = IntCol()
        electrons  = IntCol()
        neutrals  = IntCol()
        charged   = IntCol()
        pt_Jet = FloatCol()
        eta_Jet = FloatCol()
        phi_Jet = FloatCol()
        mass_Jet = FloatCol()
        et_Jet = FloatCol()
        deposits_Jet = stl.vector(float)


tree = Tree('gentree', model=GenTree)

count=1

for evt in events:
    #print 'new event'
    if count%100==0: print 'processing event ', count
    count+=1

    evt.getByLabel(labelPruned, handlePruned)
    evt.getByLabel(labelReco, handleReco)
    evt.getByLabel(labelPacked, handlePacked)    
    evt.getByLabel(labelJet, handleJet)

    # get the product                                                                                                      
    pruned = handlePruned.product()
    recoCand = handleReco.product()
    packedCand = handlePacked.product()
    recoJets = handleJet.product()

    #get the list of maj neutrino - len shoud be always one
    
    taus = [pp for pp in pruned if ((abs(pp.pdgId()) == 15) and (pp.isLastCopy()))]
    #print taus
    if len(taus) ==0 : 
        continue
    elif len(taus)==1 :
        taus = taus[0]
    elif len(taus) > 1:
        taus.sort(key=lambda x: x.pt())
        taus=taus[-1]

    final_particles = []
    final_particles_pi0 = []
    #inter_particles = []
    primaryVtx = None
    
    primaryVtx = (taus.vx(),taus.vy(),taus.vz())
    #print "majorana neutrinos mass: ", maj_n.mass() 

    #loop on the gen and look at final particles. Pi0 is considered as final particle
    for p in [i for i in pruned]:
        if ((p.isLastCopy()) and (p != taus) and (p not in final_particles)):
            if isAncestor(taus,p):
                if p.mother(0).pdgId()==111 and (p.mother(0).isLastCopy()):
                    if p.mother(0) not in final_particles :
                        final_particles.append(p.mother(0))
                        final_particles_pi0.append(p)
                else:
                    final_particles.append(p)

    abs_final_part = [ abs(x.pdgId()) for x in final_particles]
    #print abs_final_part

    if 13 in abs_final_part and 14 in abs_final_part:
        continue
    if 11 in abs_final_part and 12 in abs_final_part:
        continue

    tree.p_tauGen =   taus.p()
    tree.pt_tauGen =  taus.pt()
    tree.phi_tauGen = taus.phi()
    tree.eta_tauGen = taus.eta()
   
    #save some infos of the generated particles
    for x in final_particles:
        tree.p_Gen.push_back(x.p())
        tree.pt_Gen.push_back(x.pt())
        tree.px_Gen.push_back(x.px())
        tree.py_Gen.push_back(x.py())
        tree.pz_Gen.push_back(x.pz())
        tree.pdgId.push_back(x.pdgId())
        tree.dx_Gen.push_back(primaryVtx[0]-x.vx())
        tree.dy_Gen.push_back(primaryVtx[1]-x.vy())
        tree.dz_Gen.push_back(primaryVtx[2]-x.vz())
        tree.eta_Gen.push_back(x.eta())
        tree.phi_Gen.push_back(x.phi())
        tree.theta_Gen.push_back(x.theta())
        tree.mass_Gen.push_back(x.mass())
        tree.e_Gen.push_back(x.energy())

    #this gives the aperture
    tree.aperture_Gen = max(deltaR(taus,i) for i in final_particles)

    charged_gen = [x for x in final_particles if (abs(x.pdgId())==321 or abs(x.pdgId())==211 or abs(x.pdgId())==13 or abs(x.pdgId())==11)]

    #dizionario key: gen particle - value: lista prime n tracce con dr <0.3
    genMatchDict = {}

    for gen in final_particles:
        if abs(gen.pdgId())==321 or abs(gen.pdgId())==211 or abs(gen.pdgId())==13 or abs(gen.pdgId())==11:
            matchedList = bestMatchedParticles(gen, recoCand)
            if len(matchedList)>len(final_particles): matchedList=matchedList[:len(final_particles)]
            if matchedList: genMatchDict[gen] = matchedList
    matchedReco = set()

    while len(genMatchDict) != 0:
        #set_trace()
        minEl = min(genMatchDict.items(), key=lambda x: x[1][0][0])
        if minEl[1][0][1] not in matchedReco:
            tree.idx.push_back(final_particles.index(minEl[0]))
            tree.dr.push_back(minEl[1][0][0])
            tree.dxy_RecoPV.push_back(minEl[1][0][1].dxy())
            tree.dz_RecoPV.push_back(minEl[1][0][1].dz())
            if minEl[1][0][1].hasTrackDetails():
                tree.dxyErr_RecoPV.push_back(minEl[1][0][1].dxyError())
                tree.dzErr_RecoPV.push_back(minEl[1][0][1].dzError())
            else:
                tree.dxyErr_RecoPV.push_back(-1000)
                tree.dzErr_RecoPV.push_back(-1000)
            tree.p_Reco.push_back(minEl[1][0][1].p())
            tree.pt_Reco.push_back(minEl[1][0][1].pt())
            tree.px_Reco.push_back(minEl[1][0][1].px())
            tree.py_Reco.push_back(minEl[1][0][1].py())
            tree.pz_Reco.push_back(minEl[1][0][1].pz())

            tree.vx_Reco.push_back(minEl[1][0][1].vx())
            tree.vy_Reco.push_back(minEl[1][0][1].vy())
            tree.vz_Reco.push_back(minEl[1][0][1].vz())

            tree.eta_Reco.push_back(minEl[1][0][1].eta())
            tree.phi_Reco.push_back(minEl[1][0][1].phi())
            tree.theta_Reco.push_back(minEl[1][0][1].theta())

            matchedReco.add(minEl[1][0][1])
            genMatchDict.pop(minEl[0])
        else:
            minEl[1].pop(0)
            if not minEl[1]:
                genMatchDict.pop(minEl[0])


    #### Fill tree ####
    tree.kaonCharged = abs_final_part.count(321)
    tree.pionCharged = abs_final_part.count(211)
    tree.kaonNeutrals = abs_final_part.count(130)+abs_final_part.count(310)
    tree.pionNeutrals = abs_final_part.count(111)
    tree.muons = abs_final_part.count(13)
    tree.electrons = abs_final_part.count(11)
    tree.neutrals = abs_final_part.count(130)+abs_final_part.count(310)+abs_final_part.count(111)+abs_final_part.count(16)+abs_final_part.count(14)+abs_final_part.count(12)
    tree.charged = abs_final_part.count(321)+abs_final_part.count(211)+abs_final_part.count(13)+abs_final_part.count(11)

    ####### Let's match the jets #######

    jets = [jet for jet in recoJets if jet.pt()>20]
    matchedJet, _ = bestMatch(taus, jets)

    #pt eta phi mass btagging discriminator #isolamento
    #B tagging discriminator: look for displaced vertexes and traks in a jet
    if matchedJet:
        tree.pt_Jet = matchedJet.pt()
        tree.eta_Jet= matchedJet.eta()
        tree.phi_Jet= matchedJet.phi()
        tree.mass_Jet= matchedJet.mass()
        tree.et_Jet= matchedJet.et()
       
        jet_particles =  matchedJet.getJetConstituents()
        deposits = [0]*9
        for prt in jet_particles:
            dr = deltaR(prt, matchedJet)
            idx = min(int(floor(dr/0.05)), 8)
            deposits[idx] += prt.et()
        for i in deposits:
            tree.deposits_Jet.push_back(i)
       #set_trace()
    else:
        tree.pt_Jet = -1000
        tree.eta_Jet= -1000
        tree.phi_Jet= -1000 
        tree.mass_Jet = -1000
        tree.et_Jet= -1000  
        tree.deposits_Jet.push_back(-1000)

    tree.fill(reset=True)

tree.write()
outfile.close()

    
    

