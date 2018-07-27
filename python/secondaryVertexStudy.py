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
import os, sys
from os.path import join

# Open a file
path = "/eos/user/j/jpriscia/HNL_IVF/SAMPLE_10GEV/"
list_input = [join(path,f) for f in os.listdir( path )]


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


#events = Events(['/eos/user/j/jpriscia/HNL_IVF/SAMPLE_7GEV/ivf7GeV_35.root'])

events = Events(list_input)

#collections
handlePruned, labelPruned  = Handle ('vector<reco::GenParticle>'), 'prunedGenParticles'
handleJet, labelJet = Handle('vector<reco::GenJet>'), 'slimmedGenJetsAK8'
handleIVF, labelIVF = Handle('vector<reco::Vertex>'),'inclusiveVertexFinder'
handlePFCand, labelPFCand = Handle('vector<pat::PackedCandidate>'),'packedPFCandidates'
handleMuon, labelMuon = Handle('vector<pat::Muon>'),'slimmedMuons'
handlePV, labelPV = Handle('vector<reco::Vertex>'),'offlineSlimmedPrimaryVertices'

#######OutputTTree###############
outfile = root_open('miniAOD_10GeV_IVF.root', 'w')

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
        vx_PVReco = FloatCol()
        vy_PVReco = FloatCol()
        vz_PVReco = FloatCol()
        vx_SecReco = FloatCol()
        vy_SecReco = FloatCol()
        vz_SecReco = FloatCol()
        dx_SecReco = FloatCol()
        dy_SecReco = FloatCol()
        dz_SecReco = FloatCol()
        d3D_SecReco = FloatCol()
        vx_SecReco_T = FloatCol()
        vy_SecReco_T = FloatCol()
        vz_SecReco_T = FloatCol()
        dx_SecReco_T = FloatCol()
        dy_SecReco_T = FloatCol()
        dz_SecReco_T = FloatCol()
        d3D_SecReco_T = FloatCol()
        num_looseMuons = IntCol()
        num_tightMuons = IntCol()
        num_promptMuons = IntCol()
        num_muonsLoose_noPrompt = IntCol()
        num_muonsTight_noPrompt = IntCol()
        dist_allVtx = FloatCol()
        muon_vtx_dxy = FloatCol()
        muon_vtx_dz = FloatCol()
        muonT_vtx_dxy = FloatCol()
        muonT_vtx_dz = FloatCol()
        muonPrompt_pt = FloatCol()
        muonPrompt_px = FloatCol()
        muonPrompt_py = FloatCol()
        muonPrompt_pz = FloatCol()
        muonPrompt_eta = FloatCol()

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
    tree.vx_SecReco_T = -1000.
    tree.vy_SecReco_T = -1000.
    tree.vz_SecReco_T = -1000.
    tree.dx_SecReco_T = -1000.
    tree.dy_SecReco_T = -1000.
    tree.dz_SecReco_T = -1000.
    tree.d3D_SecReco_T = -1000.
    tree.num_tightMuons = -1000.
    tree.num_looseMuons = -1000.
    tree.num_promptMuons = -1000.
    tree.num_muonsLoose_noPrompt = -1000.
    tree.num_muonsTight_noPrompt = -1000.
    tree.dist_allVtx = -1000.
    tree.muon_vtx_dxy =  -1000.
    tree.muon_vtx_dz =   -1000.
    tree.muonT_vtx_dxy = -1000.
    tree.muonT_vtx_dz =  -1000.
    tree.vx_PVReco = -1000.
    tree.vy_PVReco = -1000.
    tree.vz_PVReco = -1000.
    tree.muonPrompt_pt =  -1000.
    tree.muonPrompt_px =  -1000.
    tree.muonPrompt_py =  -1000.
    tree.muonPrompt_pz =  -1000.
    tree.muonPrompt_eta = -1000.

    if count%100==0: print 'processing event ', count
    count+=1

    evt.getByLabel(labelPruned, handlePruned)
    evt.getByLabel(labelJet, handleJet)
    evt.getByLabel(labelIVF, handleIVF)
    evt.getByLabel(labelPFCand, handlePFCand)
    evt.getByLabel(labelMuon,handleMuon)
    evt.getByLabel(labelPV,handlePV)

    # get the product                                                                                                      
    pruned = handlePruned.product()
    recoJets = handleJet.product()
    ivf = handleIVF.product()
    pfCands = handlePFCand.product()
    muonCands = handleMuon.product()
    vtxs = handlePV.product()

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

        #5 GeV  pt cut on the second muon at generator level  LooseMuon

         
        #trigger infos
        
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

        #EVENTS IN ACCEPTANCE
        if abs(muons_gen[0].vz())>200 or sqrt(muons_gen[0].vx()*muons_gen[0].vx()+muons_gen[0].vy()*muons_gen[0].vy())>60 or abs(muons_gen[0].eta())>2.4 or muons_gen[0].pt()<5:
            tree.fill(reset=True)
            continue 

        secVtxGen = (muons_gen[0].vx(),muons_gen[0].vy(),muons_gen[0].vz())
        tree.vx_SecGen = muons_gen[0].vx()
        tree.vy_SecGen = muons_gen[0].vy()
        tree.vz_SecGen = muons_gen[0].vz()

        bestPV = None

        for vtx in vtxs:
            if not(vtx.isFake() and vtx.ndof() > 4.and  fabs(vtx.z()) < 24. and vtx.rho() < 2):
                bestPV = vtx
                tree.vx_PVReco = bestPV.x()
                tree.vy_PVReco = bestPV.y()
                tree.vz_PVReco = bestPV.z()
                break 
        if not bestPV: 
            tree.fill(reset=True)
            continue

        
        #loop on PFCandidates to see how many do have a prompt muon
        muonsPrompt=[]
        muonsLoose=[]
        muonsTight=[]
        #set_trace()
        for muon in muonCands:
            if muon.isTightMuon(bestPV) and muon.pt()>5 and abs(muon.eta())<2.4:
                muonsTight.append(muon) 
            if muon.isTightMuon(bestPV) and muon.pt()>24 and abs(muon.eta())<2.4 and muon.muonBestTrack().dxy(bestPV.position())<0.05 and muon.muonBestTrack().dz(bestPV.position())<0.1:
                muonsPrompt.append(muon)
            if muon.isLooseMuon() and muon.pt()>5 and abs(muon.eta())<2.4:
                muonsLoose.append(muon)
            

        tree.num_promptMuons = len(muonsPrompt)
        tree.num_looseMuons = len(muonsLoose)
        tree.num_tightMuons = len(muonsTight)
        
        muonsPrompt.sort(key=lambda x: -x.pt())


        if len(muonsPrompt)>0:
            tree.muonPrompt_pt =  muonsPrompt[0].pt()
            tree.muonPrompt_px =  muonsPrompt[0].px()
            tree.muonPrompt_py =  muonsPrompt[0].py()
            tree.muonPrompt_pz =  muonsPrompt[0].pz()
            tree.muonPrompt_eta =  muonsPrompt[0].eta()

            muonsLoose_noPrompt = [i for i in muonsLoose if i is not muonsPrompt[0]]
            muonsTight_noPrompt = [i for i in muonsTight if i is not muonsPrompt[0]]
        else:
            muonsLoose_noPrompt = muonsLoose
            muonsTight_noPrompt = muonsTight
        tree.num_muonsLoose_noPrompt = len(muonsLoose_noPrompt)
        tree.num_muonsTight_noPrompt = len(muonsTight_noPrompt)

        if len(muonsPrompt)==0: 
            tree.fill(reset=True)
            continue

        #small test to check how many verteces we have -- to comment
        bestVtx = None    
        delta_3D_All=100000000000.
        for vtx in ivf:
                deltaX = abs(secVtxGen[0]-vtx.x())
                deltaY = abs(secVtxGen[1]-vtx.y())
                deltaZ = abs(secVtxGen[2]-vtx.z())
                delta_3D = sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ)

                if delta_3D<delta_3D_All:
                    bestVtx = vtx
                    delta_3D_All=delta_3D
        
        if bestVtx:
            tree.dist_allVtx = delta_3D_All

        
        #tree.fill(reset=True)
        #continue

        #for vertex in ivf:
        #    dist_allVtx = sqrt(abs(secVtxGen[0]-vertex.x())*abs(secVtxGen[0]-vertex.x())+abs(secVtxGen[0]-vertex.y())*abs(secVtxGen[0]-vertex.y())+abs(secVtxGen[0]-vertex.z())*abs(secVtxGen[0]-vertex.z()))
        #    tree.dist_allVtx = dist_allVtx
        #    tree.fill(reset=True)
        #continue    
        #################

         
        # ivf vertex position
        verteces_good = []
        verteces_good_tight = []
        for vertex in ivf:
            for trk_id in range(vertex.tracksSize()): 
                trk = vertex.trackRefAt(trk_id).get()
                for muon_det in muonsLoose_noPrompt:
                    if deltaR(trk,muon_det)<0.4:
                        verteces_good.append((vertex,muon_det,deltaR))

                for muon_t in muonsTight_noPrompt:

                    if deltaR(trk,muon_t)<0.4:
                        verteces_good_tight.append((vertex,muon_t,deltaR))
    
        if len(verteces_good)>0:
            verteces_good.sort(key=lambda tup: tup[2])
            bestRecoVtx = verteces_good[0][0]
            tree.muon_vtx_dxy = verteces_good[0][1].muonBestTrack().dxy(bestRecoVtx.position())
            tree.muon_vtx_dz  = verteces_good[0][1].muonBestTrack().dz(bestRecoVtx.position())

            deltaX = abs(secVtxGen[0]-bestRecoVtx.x())
            deltaY = abs(secVtxGen[1]-bestRecoVtx.y())
            deltaZ = abs(secVtxGen[2]-bestRecoVtx.z())
            delta_3D = sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ)

            tree.vx_SecReco = bestRecoVtx.x()
            tree.vy_SecReco = bestRecoVtx.y()
            tree.vz_SecReco = bestRecoVtx.z()
            tree.dx_SecReco = deltaX
            tree.dy_SecReco = deltaY
            tree.dz_SecReco = deltaZ
            tree.d3D_SecReco = delta_3D
        

        if len(verteces_good_tight)>0:

            verteces_good_tight.sort(key=lambda tup: tup[2])
            bestRecoVtxTight = verteces_good_tight[0][0]
            tree.muonT_vtx_dxy = verteces_good_tight[0][1].muonBestTrack().dxy(bestRecoVtxTight.position())
            tree.muonT_vtx_dz  = verteces_good_tight[0][1].muonBestTrack().dz(bestRecoVtxTight.position())

            deltaX_T = abs(secVtxGen[0]-bestRecoVtxTight.x())
            deltaY_T = abs(secVtxGen[1]-bestRecoVtxTight.y())
            deltaZ_T = abs(secVtxGen[2]-bestRecoVtxTight.z())
            delta_3D_T = sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ)

            tree.vx_SecReco_T = bestRecoVtxTight.x()
            tree.vy_SecReco_T = bestRecoVtxTight.y()
            tree.vz_SecReco_T = bestRecoVtxTight.z()
            tree.dx_SecReco_T = deltaX_T
            tree.dy_SecReco_T = deltaY_T
            tree.dz_SecReco_T = deltaZ_T
            tree.d3D_SecReco_T = delta_3D_T
        


        tree.fill(reset=True)

tree.write()
outfile.close()

    
    

