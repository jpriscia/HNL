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


print 'entra'


def bestMatchedJet( object, matchCollection):
    '''Return the best match to object in matchCollection, which is the closest object in delta R'''
    deltaR2Min = float('+inf')
    bm = None
    matchList=[]
    for match in matchCollection:
        matchDone=False
        for charged in match.getPFConstituents() :
            if charged.charge()==0: continue
            if not charged.bestTrack(): continue
            track = charged.bestTrack()
            if abs(track.dxy(object.vertex()))<50 and abs(track.dz(object.vertex()))<50: matchDone=True  #was 0.5 for both dxy and dz
           
        if matchDone:
            dR2 = deltaR2( object.eta(), object.phi(),
                       match.eta(), match.phi() )
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

lines = open("fileList_mu4GeV_TEST.txt").read().splitlines()

events = Events(lines)
########### for test uncomment this ############
#events = Events(['root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_208.root'])
#event = events.__iter__().next()


#########outputfile###########
outfile = root_open('jetMatch4GeV_test.root', 'w')

class GenTree(TreeModel):

        p_NuGen = FloatCol()
        pt_NuGen = FloatCol()
        phi_NuGen = FloatCol()
        eta_NuGen = FloatCol()
        theta_NuGen = FloatCol()
        px_NuGen = FloatCol()
        py_NuGen = FloatCol()
        pz_NuGen = FloatCol()
        vx_NuGen = FloatCol()
        vy_NuGen = FloatCol()
        vz_NuGen = FloatCol()

        deltaR_HNL = FloatCol()
        #jet_HNL_p = FloatCol()
        jet_HNL_pt = FloatCol()
        jet_HNL_eta = FloatCol()
        jet_HNL_phi = FloatCol()

        deltaR_muon = FloatCol()
        #jet_muon_p = FloatCol()
        jet_muon_pt = FloatCol()
        jet_muon_eta = FloatCol()
        jet_muon_phi = FloatCol()
        
        frac_HNLjet = FloatCol()
        frac_muonjet = FloatCol()
        numbOfDaughters = IntCol()

tree = Tree('gentree', model=GenTree)



handleReco, labelReco  = Handle('vector<reco::Track>'),'generalTracks'
handlePruned, labelPruned  = Handle ('std::vector<reco::GenParticle>'), 'genParticles'

handleJet, labelJet = Handle('vector<reco::PFJet>'), 'ak4PFJetsCHS'

handlePileup, labelPileup = Handle('vector<PileupSummaryInfo>'), 'addPileupInfo'
#handlePFCand, labelPFCand = Handle('vector<pat::PackedCandidate>'), 'packedPFCandidates'

count=1

for evt in events:
    if count%100==0: print 'processing event ', count
    count+=1

    evt.getByLabel(labelPruned, handlePruned)
    evt.getByLabel(labelReco, handleReco)
    evt.getByLabel(labelJet, handleJet)
    evt.getByLabel(labelPileup,handlePileup)
    #evt.getByLabel(labelPFCand, handlePFCand)

    # get the product                                                                                                      
    pruned = handlePruned.product()
    recoCand = handleReco.product()
    recoJets = handleJet.product()
    pileup = handlePileup.product()
    #pfCand = handlePFCand.product()

    #for i in pileup:
    #    print i.getTrueNumInteractions()
    #    print i.getPU_NumInteractions()
    #print ' '
    #print [i for i in dir(pileup[0])]
    #set_trace()
    maj_neutrinos = [pp for pp in pruned if ((abs(pp.pdgId()) == 9900012 or abs(pp.pdgId()) == 9900014 or abs(pp.pdgId()) == 9900016) and pp.isLastCopy())]
    #pfCandidatesLeptons = [pp for pp in pfCand if (abs(pp.pdgId()) == 13 or abs(pp.pdgId()) == 11)]
    #set_trace()
    #get the list of maj neutrino - len shoud be always one

    #print [i.pdgId() for i in pfCandidatesLeptons]
    #continue 
    if len(maj_neutrinos) != 1:
        set_trace()
        print "we have multiple majorana neutrinos"

    final_particles = []
    final_particles_pi0 = []
    inter_particles = []
    primaryVtx = None
    dMeson = 0
    cQuark = 0
    for maj_n in maj_neutrinos:
        tree.p_NuGen =   maj_n.p()
        tree.px_NuGen = maj_n.px()
        tree.py_NuGen = maj_n.py()
        tree.pz_NuGen = maj_n.pz()
        tree.pt_NuGen =  maj_n.pt()
        tree.phi_NuGen = maj_n.phi()
        tree.eta_NuGen = maj_n.eta()
        tree.theta_NuGen = maj_n.theta()
        tree.vx_NuGen = maj_n.vx()
        tree.vy_NuGen = maj_n.vy()
        tree.vz_NuGen = maj_n.vz()

        primaryVtx = (maj_n.vx(),maj_n.vy(),maj_n.vz())

        for p in pruned:
            if (p.status()==1 and p.isLastCopy() and (p != maj_n) and (p not in final_particles)):
                if isAncestor(maj_n,p):
                    if p.mother(0).pdgId()==111 and (p.mother(0).isLastCopy()):
                        if p.mother(0) not in final_particles :
                            final_particles.append(p.mother(0))
                            final_particles_pi0.append(p)
                    else:
                        final_particles.append(p)


    tree.numbOfDaughters =  len(final_particles)
    muon_gen = [(x,x.pt()) for x in final_particles if (abs(x.pdgId())==13)]
    muon_gen.sort(key=lambda tup: tup[1], reverse=True) 


    ####### Let's match the jets #######

    jets = [jet for jet in recoJets if jet.pt() > 20]
    #print len(jets)
    if len(jets)!=0:
        jet_HNL, deltaR_HNL = bestMatch(maj_neutrinos[0], jets)

        tree.deltaR_HNL = deltaR_HNL
    #tree.jet_HNL_p = jet_HNL.p()
        tree.jet_HNL_pt = jet_HNL.pt()
        tree.jet_HNL_eta = jet_HNL.eta()
        tree.jet_HNL_phi = jet_HNL.phi()

    #let's check the fraction of the HNL daughters within a cone of .4 from the jet axes
        daughter_deltaR_HNL = []
        for part in final_particles: 
            daughter_deltaR_HNL.append(deltaR(part,jet_HNL))
        tree.frac_HNLjet = float(sum(i < 0.4 for i in daughter_deltaR_HNL))/float(len(daughter_deltaR_HNL))
        
    
        matchJetlist =  bestMatchedJet(muon_gen[0][0], jets)
        if len(matchJetlist)>0:
            deltaR_muon =  matchJetlist[0][0]
            jet_muon =  matchJetlist[0][1]
            tree.deltaR_muon = deltaR_muon
       
    #tree.jet_muon_p = jet_muon.p()
            tree.jet_muon_pt = jet_muon.pt()
            tree.jet_muon_eta = jet_muon.eta()
            tree.jet_muon_phi = jet_muon.phi()

    #same stuff for the jet closer to the muon
            daughter_deltaR_muon= []
            for part in final_particles:
                daughter_deltaR_muon.append(deltaR(part,jet_muon))
                tree.frac_muonjet = float(sum(i < 0.4 for i in daughter_deltaR_muon))/float(len(daughter_deltaR_muon))

        else:
            tree.frac_muonjet = -100
            tree.deltaR_muon= -100
            tree.jet_muon_pt= -100
            tree.jet_muon_eta = -100
            tree.jet_muon_phi = -100
    else:
        tree.frac_HNLjet = -100
        tree.frac_muonjet = -100
        tree.deltaR_HNL = -100
        tree.jet_HNL_pt = -100
        tree.jet_HNL_eta = -100
        tree.jet_HNL_phi = -100
        tree.deltaR_muon= -100
        tree.jet_muon_pt= -100
        tree.jet_muon_eta = -100
        tree.jet_muon_phi = -100

    tree.fill(reset=True)

tree.write()
outfile.close()

    
    

