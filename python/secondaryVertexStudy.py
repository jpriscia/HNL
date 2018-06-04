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


print 'entra'

vtxPlot_IVF = plt.Hist1D(500,-100,100, name="deltaZ_IVF", title="deltaZ_IVF")
vtxPlot_SV = plt.Hist1D(500,-5,5, name="deltaZ_SV", title="deltaZ_SV")
vtxPlot_CSV = plt.Hist1D(500,-5,5, name="deltaZ_CSV", title="deltaZ_CSV")
vtxPlot_CSVVSL = plt.Hist1D(500,-5,5, name="deltaZ_CSVVSL", title="deltaZ_CSVVSL")

vrtxMatch=plt.Hist2D(500,-100,100,500,0,100,name='matchedVtx',title='matchedVtx')
vrtxNoMatch=plt.Hist2D(500,-100,100,500,0,100,name='noMatchedVtx',title='noMatchedVtx')

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

#lines = open("/afs/cern.ch/user/j/jpriscia/CMSSW_HNL_17/src/HNL/HNL/src/secondVertex/prova_AOD.root").read().splitlines()

events = Events('/afs/cern.ch/user/j/jpriscia/CMSSW_HNL_17/src/HNL/HNL/src/secondVertex/prova_AOD.root')

handleReco, labelReco  = Handle('vector<reco::Track>'),'generalTracks'
handlePruned, labelPruned  = Handle ('std::vector<reco::GenParticle>'), 'genParticles'
handleJet, labelJet = Handle('vector<reco::PFJet>'), 'ak4PFJetsCHS'

handleSV, labelSV = Handle('vector<reco::Vertex>'), 'inclusiveSecondaryVertices'
handleCSV, labelCSV = Handle('vector<reco::VertexCompositePtrCandidate>'), 'inclusiveCandidateSecondaryVertices'
handleCSVvsL, labelCSVvsL = Handle('vector<reco::VertexCompositePtrCandidate>'), 'inclusiveCandidateSecondaryVerticesCvsL'
handleIVF, labelIVF = Handle('vector<reco::Vertex>'),'inclusiveVertexFinder'


#######OutputTTree###############
outfile = root_open('prova.root', 'w')

class GenTree(TreeModel):
        
        p_NuGen = FloatCol()
        pt_NuGen = FloatCol()
        phi_NuGen = FloatCol()
        eta_NuGen = FloatCol()
        theta_NuGen = FloatCol()
        aperture_Gen = FloatCol()
        px_NuGen = FloatCol()
        py_NuGen = FloatCol()
        pz_NuGen = FloatCol()
        vx_NuGen = FloatCol()
        vy_NuGen = FloatCol()
        vz_NuGen = FloatCol()
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
        vx_Gen = stl.vector(float)
        vy_Gen = stl.vector(float)
        vz_Gen = stl.vector(float)

        p_Reco =  stl.vector(float)
        pt_Reco = stl.vector(float)
        px_Reco  = stl.vector(float)
        py_Reco  = stl.vector(float)
        pz_Reco  = stl.vector(float)

        vx_Reco  = stl.vector(float)
        vy_Reco  = stl.vector(float)
        vz_Reco  = stl.vector(float)

        eta_Reco  = stl.vector(float)
        etaErr_Reco  = stl.vector(float)
        phi_Reco  = stl.vector(float)
        phiErr_Reco  = stl.vector(float)
        theta_Reco  = stl.vector(float)
        thetaErr_Reco  = stl.vector(float)

        idx = stl.vector(int)
        pdgId = stl.vector(int)
        dr = stl.vector(float)
        dx_Gen = stl.vector(float)
        dy_Gen = stl.vector(float)
        dz_Gen = stl.vector(float)
        algo = stl.vector(int)
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
        dMeson = IntCol()
        cQuark = IntCol()

        pt_Jet = FloatCol()
        eta_Jet = FloatCol()
        phi_Jet = FloatCol()
        mass_Jet = FloatCol()
        et_Jet = FloatCol()
        deposits_Jet = stl.vector(float)


tree = Tree('gentree', model=GenTree)

count=0
csvvsl_count=0
csv_count=0
sv_count=0
ivf_count=0
ivf_all=0
ivf_good=0
maj_in_accept=0

for evt in events:
    if count%100==0: print 'processing event ', count
    count+=1

    evt.getByLabel(labelPruned, handlePruned)
    evt.getByLabel(labelReco, handleReco)
    evt.getByLabel(labelJet, handleJet)
    evt.getByLabel(labelSV, handleSV)  #vector<reco::Vertex> "inclusiveSecondaryVertices"
    evt.getByLabel(labelCSV, handleCSV) #vector<reco::VertexCompositePtrCandidate>    "inclusiveCandidateSecondaryVertices" 
    evt.getByLabel(labelCSVvsL,handleCSVvsL)  #vector<reco::VertexCompositePtrCandidate>    "inclusiveCandidateSecondaryVerticesCvsL"
    evt.getByLabel(labelIVF, handleIVF)

    # get the product                                                                                                      
    pruned = handlePruned.product()
    recoCand = handleReco.product()
    recoJets = handleJet.product()
    sv = handleSV.product()
    csv = handleCSV.product()
    csvvsl = handleCSVvsL.product()
    ivf = handleIVF.product()

    #if count ==10: set_trace()
    #get the list of maj neutrino - len shoud be always one
    maj_neutrinos = [pp for pp in pruned if ((abs(pp.pdgId()) == 9900012 or abs(pp.pdgId()) == 9900014 or abs(pp.pdgId()) == 9900016) and pp.isLastCopy())]

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
        #print "majorana neutrinos mass: ", maj_n.mass() 

        #loop on the gen and look at final particles. Pi0 is considered as final particle
        for p in pruned:
            #print p.vertex().x(),p.vertex().y(),p.vertex().z()
            #print p.vx(), p.vy(), p.vz()
            #print ''
            if (p.status()==1 and p.isLastCopy() and (p != maj_n) and (p not in final_particles)):
                if isAncestor(maj_n,p):
                    if p.mother(0).pdgId()==111 and (p.mother(0).isLastCopy()):
                        if p.mother(0) not in final_particles :
                            final_particles.append(p.mother(0))
                            final_particles_pi0.append(p)
                    else:
                        final_particles.append(p)

            #I want to check if there is a D meson as well...
            if (p.status()!=1 and p.isLastCopy() and (p != maj_n) and (p not in inter_particles)):
                if isAncestor(maj_n,p):
                    inter_particles.append(p)
                    if abs(p.pdgId())==431:
                        dMeson+=1
                    if abs(p.pdgId())==4:
                        cQuark+=1
                        #print cQuark
                
        abs_final_part = [ abs(x.pdgId()) for x in final_particles]

        vertex_pos = []
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
            tree.vx_Gen.push_back(x.vx())
            tree.vy_Gen.push_back(x.vy())
            tree.vz_Gen.push_back(x.vz())
            tree.eta_Gen.push_back(x.eta())
            tree.phi_Gen.push_back(x.phi())
            tree.theta_Gen.push_back(x.theta())
            tree.mass_Gen.push_back(x.mass())
            tree.e_Gen.push_back(x.energy())


            
        for x in final_particles:
            if abs(x.pdgId())==13 :
                vertex_pos.append(x.vx())
                vertex_pos.append(x.vy())
                vertex_pos.append(x.vz())
                vertex_pos.append(x.eta())
        
        #if abs(vertex_pos[2])<200 and abs(sqrt(vertex_pos[0]*vertex_pos[0]+vertex_pos[1]*vertex_pos[1]))<68  :
            #print "maj neutrino decay: ", vertex_pos
        delta_csv = []
        delta_csvvsl = []
        delta_sv = []
        delta_ivf = []

        if abs(vertex_pos[2])<200 and abs(sqrt(vertex_pos[0]*vertex_pos[0]+vertex_pos[1]*vertex_pos[1]))<68 and abs(vertex_pos[3])<2.5 :
            maj_in_accept+=1
            #print "majorana vertex:", vertex_pos[0], vertex_pos[1], vertex_pos[2]
            for vertex in ivf:
                deltaX = abs(vertex_pos[0]-vertex.x())
                deltaY = abs(vertex_pos[1]-vertex.y())
                deltaZ = abs(vertex_pos[2]-vertex.z())
                delta_ivf.append((deltaX,deltaY,deltaZ,sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ)))
                ivf_all+=1
            for vertex in csv: 
                deltaX = abs(vertex_pos[0]-vertex.vx()) 
                deltaY = abs(vertex_pos[1]-vertex.vy())
                deltaZ = abs(vertex_pos[2]-vertex.vz())
                delta_csv.append((deltaX,deltaY,deltaZ,sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ)))

            for vertex in csvvsl:
                deltaX = abs(vertex_pos[0]-vertex.vx())
                deltaY = abs(vertex_pos[1]-vertex.vy())
                deltaZ = abs(vertex_pos[2]-vertex.vz())
                delta_csvvsl.append((deltaX,deltaY,deltaZ,sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ)))
                
            for vertex in sv:
                deltaX = abs(vertex_pos[0]-vertex.x())
                deltaY = abs(vertex_pos[1]-vertex.y())
                deltaZ = abs(vertex_pos[2]-vertex.z())
                delta_sv.append((deltaX,deltaY,deltaZ,sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ)))

            delta_ivf.sort(key=lambda tup: tup[3])
            delta_csv.sort(key=lambda tup: tup[3])
            delta_csvvsl.sort(key=lambda tup: tup[3])
            delta_sv.sort(key=lambda tup: tup[3])
            
            if len(delta_ivf)>0:
                vtxPlot_IVF.fill(delta_ivf[0][3])
                ivf_count+=1
                if delta_ivf[0][3]<1:
                    vrtxMatch.fill(vertex_pos[2],sqrt(vertex_pos[0]*vertex_pos[0]+vertex_pos[1]*vertex_pos[1]))
                    ivf_good+=1
                if delta_ivf[0][3]>1: vrtxNoMatch.fill(vertex_pos[2],sqrt(vertex_pos[0]*vertex_pos[0]+vertex_pos[1]*vertex_pos[1]))

            if len(delta_sv)>0:
                vtxPlot_SV.fill(delta_sv[0][3])
                #if delta_sv[0][3]<1: print "SV", delta_sv[0][0], delta_sv[0][1], delta_sv[0][2]
                sv_count+=1
            if len(delta_csv)>0:
                vtxPlot_CSV.fill(delta_csv[0][3])
                #if delta_csv[0][3]<1: print "CSV", delta_csv[0][0], delta_csv[0][1], delta_csv[0][2]
                csv_count+=1
            if len(delta_csvvsl)>0:
                vtxPlot_CSVVSL.fill(delta_csvvsl[0][3])
                #if delta_csvvsl[0][3]<1: print "CSVVSL", delta_csvvsl[0][0], delta_csvvsl[0][1], delta_csvvsl[0][2]
                csvvsl_count+=1

        tree.fill(reset=True)
       
print "maj_in_accept: ", maj_in_accept

#print sv_count, csv_count, csvvsl_count
print "event with sec vertex : ", ivf_count
print "secondary vertex all: ", ivf_all
print "good secondary vertex: ", ivf_good
cmap=cm.get_cmap('afmhot')
cmap.set_under('w')
cmap.set_bad('gray')

fig_Match, ax_Match = mplt.subplots(figsize=(15, 6))
fig_Match.patch.set_facecolor('white')
vrtxMatch.SetLineColor(2)
vrtxNoMatch.SetLineColor(1)
rplt.hist(vrtxMatch,axes=ax_Match)
rplt.hist(vrtxNoMatch,axes=ax_Match)
mplt.show()

fig_IVF, ax_IVF = mplt.subplots(figsize=(15, 6))
rplt.hist(vtxPlot_IVF,axes=ax_IVF,colorbar=True,cmap=cmap,vmin=.001)
mplt.title('IVF')
ax_IVF.set_xlabel('z [cm]')
mplt.show()

fig_SV, ax_SV = mplt.subplots(figsize=(15, 6))
rplt.hist(vtxPlot_SV,axes=ax_SV,colorbar=True,cmap=cmap,vmin=.001)
mplt.title('SV')
ax_SV.set_xlabel('z [cm]')
mplt.show()

fig_CSV, ax_CSV = mplt.subplots(figsize=(15, 6))
rplt.hist(vtxPlot_CSV,axes=ax_CSV,colorbar=True,cmap=cmap,vmin=.001)
mplt.title('CSV')
ax_CSV.set_xlabel('z [cm]')
mplt.show()


fig_CSVVSL, ax_CSVVSL = mplt.subplots(figsize=(15, 6))
rplt.hist(vtxPlot_CSVVSL,axes=ax_CSVVSL,colorbar=True,cmap=cmap,vmin=.001)
mplt.title('CSVVSL')
ax_CSVVSL.set_xlabel('z [cm]')
mplt.show()



tree.write()
outfile.close()

    
    

