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
from math import sqrt
import rootpy.plotting.root2matplotlib as rplt
from rootpy.plotting.style import get_style, set_style
import matplotlib.pyplot as mplt
from pylab import arange, show, cm
from pdb import set_trace


def bestMatchedParticles( object, matchCollection):
    '''Return the best match to object in matchCollection, which is the closest object in delta R'''
    deltaR2Min = float('+inf')
    bm = None
    matchList=[]
    for match in matchCollection:
        if abs(match.dxy(object.vertex()))<0.5 and abs(match.dz(object.vertex()))<0.5:
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


events = Events('/afs/cern.ch/user/j/jpriscia/CMSSW_HNL_17/src/HNL/HNL/src/secondVertex/prova_miniAOD.root')


#collections
handlePruned, labelPruned  = Handle ('vector<reco::GenParticle>'), 'prunedGenParticles'
handleJet, labelJet = Handle('vector<reco::GenJet>'), 'slimmedGenJetsAK8'
handleIVF, labelIVF = Handle('vector<reco::Vertex>'),'inclusiveVertexFinder'
handlePFCand, labelPFCand = Handle('vector<pat::PackedCandidate>'),'packedPFCandidates'


#######OutputTTree###############
outfile = root_open('miniAOD_IVF.root', 'w')

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
        vx_SecGen = FloatCol()
        vy_SecGen = FloatCol()
        vz_SecGen = FloatCol()
        vx_SecReco = FloatCol()
        vy_SecReco = FloatCol()
        vz_SecReco = FloatCol()
        dx_SecReco = FloatCol()
        dy_SecReco = FloatCol()
        dz_SecReco = FloatCol()
        d3D_SecReco = FloatCol()

tree = Tree('gentree', model=GenTree)

count=0

for evt in events:

    #initialize
    tree.p_NuGen = -1000.
    tree.pt_NuGen = -1000.
    tree.phi_NuGen = -1000.
    tree.eta_NuGen = -1000.
    tree.theta_NuGen = -1000.
    tree.px_NuGen = -1000.
    tree.py_NuGen = -1000.
    tree.pz_NuGen = -1000.
    tree.vx_NuGen = -1000.
    tree.vy_NuGen = -1000.
    tree.vz_NuGen = -1000.
    tree.vx_SecGen = -1000.
    tree.vy_SecGen = -1000.
    tree.vz_SecGen = -1000.
    tree.vx_SecReco = -1000.
    tree.vy_SecReco = -1000.
    tree.vz_SecReco = -1000.
    tree.dx_SecReco = -1000.
    tree.dy_SecReco = -1000.
    tree.dz_SecReco = -1000.
    tree.d3D_SecReco = -1000.




    if count%100==0: print 'processing event ', count
    count+=1

    evt.getByLabel(labelPruned, handlePruned)
    evt.getByLabel(labelJet, handleJet)
    evt.getByLabel(labelIVF, handleIVF)
    evt.getByLabel(labelPFCand, handlePFCand)

    # get the product                                                                                                      
    pruned = handlePruned.product()
    recoJets = handleJet.product()
    ivf = handleIVF.product()
    pfCands = handlePFCand.product()

    #set_trace()

    # list of mak neutrinos
    maj_neutrinos = [pp for pp in pruned if ((abs(pp.pdgId()) == 9900012 or abs(pp.pdgId()) == 9900014 or abs(pp.pdgId()) == 9900016) and pp.isLastCopy())]

    if len(maj_neutrinos) != 1:
        set_trace()
        print "we have multiple majorana neutrinos"

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


        #loop on the gen and look at final particles. Pi0 is considered as final particle
        final_particles = []
        final_particles_pi0 = []
        for p in pruned:
            if (p.status()==1 and p.isLastCopy() and (p != maj_n) and (p not in final_particles)):
                if isAncestor(maj_n,p):
                    if p.mother(0).pdgId()==111 and (p.mother(0).isLastCopy()):
                        if p.mother(0) not in final_particles :
                            final_particles.append(p.mother(0))
                            final_particles_pi0.append(p)
                    else:
                        final_particles.append(p)

        abs_final_part = [ abs(x.pdgId()) for x in final_particles]

        #take muons
        muons_gen = [x for x in final_particles if abs(x.pdgId())==13]

        # lambda to order in momentum. The one with higher p should be 2nd muon
        muons_gen.sort(key=lambda x: -x.pt())

        if abs(muons_gen[0].vz())>200 or sqrt(muons_gen[0].vx()*muons_gen[0].vx()+muons_gen[0].vy()*muons_gen[0].vy())>60 or abs(muons_gen[0].eta())>2.5:
            tree.fill(reset=True)
            continue 

        secVtxGen = (muons_gen[0].vx(),muons_gen[0].vy(),muons_gen[0].vz())
        tree.vx_SecGen = muons_gen[0].vx()
        tree.vy_SecGen = muons_gen[0].vy()
        tree.vz_SecGen = muons_gen[0].vz()

        #loop on PFCandidates to see how many do have a prompt muon
        muonsPrompt=[]
        muonsDetached=[]
        for cand in pfCands:
            if (cand.isGlobalMuon() or cand.isMuon() or cand.isTrackerMuon() or cand.isCaloMuon()) and cand.dxy()<0.1:
                #print 'muon found'
                muonsPrompt.append(cand)
            if (cand.isGlobalMuon() or cand.isMuon() or cand.isTrackerMuon() or cand.isCaloMuon()) and cand.dxy()>0.1:
                muonsDetached.append(cand)
        #print muonsDetached
        #if none, I don't even bother to check the secondary vertex
        if len(muonsPrompt)==0 or len(muonsDetached)==0:
            tree.fill(reset=True)
            continue
         
        # ivf vertex position
        verteces_good = []
        for vertex in ivf:
            for trk_id in range(vertex.tracksSize()): 
                trk = vertex.trackRefAt(trk_id).get()
                #print trk.algo()
                if trk.dxy()>0.1:
                    for muon_det in muonsDetached:
                        if deltaR(trk,muon_det)<0.05:
                            verteces_good.append(vertex)
                            break   #if I have the detached muon, it is a good vertex!

        
        if len(verteces_good)==0: 
            tree.fill(reset=True)
            continue
        
        #I match to the generated vertex and take the best
        bestVtx = None
        delta_3D_max = 1000000000
        for vtx in verteces_good:
                deltaX = abs(secVtxGen[0]-vertex.x())
                deltaY = abs(secVtxGen[1]-vertex.y())
                deltaZ = abs(secVtxGen[2]-vertex.z())
                delta_3D = sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ)
                
                if delta_3D<delta_3D_max:
                    bestVtx = vtx
                    delta_3D_max=delta_3D
                    

        tree.vx_SecReco = bestVtx.x()
        tree.vy_SecReco = bestVtx.y()
        tree.vz_SecReco = bestVtx.z()
        tree.dx_SecReco = abs(secVtxGen[0]-bestVtx.x())
        tree.dy_SecReco = abs(secVtxGen[0]-bestVtx.y())
        tree.dz_SecReco = abs(secVtxGen[0]-bestVtx.z())
        tree.d3D_SecReco = sqrt(abs(secVtxGen[0]-bestVtx.x())*abs(secVtxGen[0]-bestVtx.x())+abs(secVtxGen[0]-bestVtx.y())*abs(secVtxGen[0]-bestVtx.y())+abs(secVtxGen[0]-bestVtx.z())*abs(secVtxGen[0]-bestVtx.z()))
        
        tree.fill(reset=True)

tree.write()
outfile.close()

    
    

