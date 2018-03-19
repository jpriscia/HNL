import ROOT
import rootpy
import rootpy.plotting as plt
from rootpy import stl
from rootpy.io import root_open
from rootpy.tree import Tree, TreeModel, FloatCol, IntCol
import pprint
from DataFormats.FWLite import Events, Handle
from ROOT import gSystem

from pdb import set_trace

inputfile = root_open('tree_2GeVlljj.root')
print inputfile.name
outfile = root_open('outputPlot2GeV.root','recreate')
outfile.cd()

#HPS
#shrinking cone
#apertura vs pt neutrinone 
#prendere quantile al 90% e farne un fit
#eff jet reconstruction vs pt neutrinone
#risoluzione in pt and mass del jet vs pt neutrinone
#risoluzione =  gen-reco/gen

if '2GeV' in inputfile.name: massNeutrino = 2.
elif '4GeV' in inputfile.name: massNeutrino = 4.
else: 'define mass Neutrino', set_trace()


h_apVSpt_dM = plt.Hist2D(100, 0, 100, 100, 0, 3, name="aperture2D_dM", title="aperture vs nu_pt (with D Meson) -- gen level")
h_apVSpt_nodM = plt.Hist2D(100, 0, 100, 100, 0, 3, name="aperture2D_nodM", title="aperture vs nu_pt (no D Meson) -- gen level")
h_apVSpt = plt.Hist2D(100, 0, 100, 100, 0, 3, name="aperture2D", title="aperture vs nu_pt -- gen level")

h_trackEff = plt.Hist(300, 0, 0.3, name="trackEff", title="track matching")
h_algo = plt.Hist(16, 0, 16, name="algo", title="algo")

h_neutr = plt.Hist(5,0,150)
h_neutr_jet = plt.Hist(5,0,150,name="jetEff", title="jet Reco Efficiency")

h_01pt = plt.Hist(5,0,150,name="drEff01", title="jet efficiency cone 0.1")
h_015pt = plt.Hist(5,0,150,name="drEff015", title="jet efficiency cone 0.15")
h_02pt = plt.Hist(5,0,150,name="drEff02", title="jet efficiency cone 0.2")
h_025pt = plt.Hist(5,0,150,name="drEff025", title="jet efficiency cone 0.25")
h_03pt = plt.Hist(5,0,150,name="drEff3", title="jet efficiency cone 0.3")
h_totpt = plt.Hist(5,0,150,name="drEfftot", title="jet efficiency cone tot")

profile_mass = plt.Profile(350,0,350,name="massProfile", title="massProfile",option="s")
profile_pt = plt.Profile(350,0,350,name="ptProfile", title="ptProfile",option="s")


for entry in inputfile.gentree:
    h_apVSpt.fill(entry.pt_NuGen, entry.aperture_Gen)
    for val in entry.dr:
        h_trackEff.fill(val)
    for val in entry.algo:
        h_algo.fill(val)
    if entry.dMeson !=0:
        h_apVSpt_dM.fill(entry.pt_NuGen, entry.aperture_Gen)

    if entry.dMeson == 0:
        h_apVSpt_nodM.fill(entry.pt_NuGen, entry.aperture_Gen)

    h_neutr.fill(entry.pt_NuGen)
    if entry.pt_Jet>-1:  
        h_neutr_jet.fill(entry.pt_NuGen)
        profile_mass.Fill(entry.pt_NuGen,(massNeutrino-entry.mass_Jet)/massNeutrino,1)
        profile_pt.Fill(entry.pt_NuGen,(entry.pt_NuGen-entry.pt_Jet)/entry.pt_NuGen,1)
        
        #if entry.aperture_Gen
        
        #dep_list = [i for i in entry.deposits_Jet]
        #print dep_list
        #for i in range(len(dep_list)-1):
        #    print sum(dep_list[i+1:])/entry.et_Jet

from rootpy.plotting import Graph
graphQuant = Graph(h_apVSpt.GetNbinsX())
graphMedian = Graph(h_apVSpt.GetNbinsX())
graphQuant_dM = Graph(h_apVSpt_dM.GetNbinsX())
graphMedian_dM = Graph(h_apVSpt_dM.GetNbinsX())
graphQuant_nodM = Graph(h_apVSpt_nodM.GetNbinsX())
graphMedian_nodM = Graph(h_apVSpt_nodM.GetNbinsX())

def quantile(histo2D,graphQ,graphM):
    for ix in range(1,histo2D.GetNbinsX()):
        values = [(histo2D[ix,iy].y.center, histo2D[ix,iy].value) for iy in range(1,histo2D.GetNbinsY())]
        integral = sum(i for _, i in values)
        #set_trace()
    #somma consecutivamente gli elementi di una lista e gli applica consecutivamente una funzione a due argomenti
    #il cui primo e` il risultato della stessa operazione sull'elemento precedente
        accumulation = reduce(
            lambda c, x: c + [c[-1] + x], 
            [i for _, i in values], [0])[1:]
        points = zip([i for i,_ in values], accumulation)
        median, _ = min(points, key=lambda x: abs(x[1] - integral/2.))
        quant90, _ = min(points, key=lambda x: abs(x[1] - integral*.9))
        graphM.SetPoint(ix-1, histo2D[ix,0].x.center,median)
        graphQ.SetPoint(ix-1, histo2D[ix,0].x.center,quant90)
   
quantile(h_apVSpt,graphQuant,graphMedian)
quantile(h_apVSpt,graphQuant_dM,graphMedian_dM)
quantile(h_apVSpt,graphQuant_nodM,graphMedian_nodM)

canv_apVSpt = plt.Canvas(800,600)
h_apVSpt.xaxis.title = 'nu_pt'
h_apVSpt.yaxis.title = 'aperture'
h_apVSpt.Draw('COLZ')
graphQuant.SetMarkerStyle(7)
graphMedian.SetMarkerStyle(7)
graphQuant.SetMarkerColor(2)
graphMedian.SetMarkerColor(2)
graphQuant.Draw('PE SAME')
graphMedian.Draw('PE SAME')
canv_apVSpt.Draw()

canv_apVSpt_nodM = plt.Canvas(800,600)
h_apVSpt_nodM.xaxis.title = 'nu_pt'
h_apVSpt_nodM.yaxis.title = 'aperture'
h_apVSpt_nodM.Draw('COLZ')
graphQuant_nodM.SetMarkerStyle(7)
graphMedian_nodM.SetMarkerStyle(7)
graphQuant_nodM.SetMarkerColor(2)
graphMedian_nodM.SetMarkerColor(2)
graphQuant_nodM.Draw('PE SAME')
graphMedian_nodM.Draw('PE SAME')
canv_apVSpt_nodM.Draw()

canv_apVSpt_dM = plt.Canvas(800,600)
h_apVSpt_dM.xaxis.title = 'nu_pt'
h_apVSpt_dM.yaxis.title = 'aperture'
h_apVSpt_dM.Draw('COLZ')
graphQuant_dM.SetMarkerStyle(7)
graphMedian_dM.SetMarkerStyle(7)
graphQuant_dM.SetMarkerColor(2)
graphMedian_dM.SetMarkerColor(2)
graphQuant_dM.Draw('PE SAME')
graphMedian_dM.Draw('PE SAME')
canv_apVSpt_dM.Draw()

canvas_trackEff = plt.Canvas(800,600)
h_trackEff.xaxis.title = 'dr'
h_trackEff.yaxis.title = 'nTracks'
h_trackEff.SetMarkerStyle(7)
h_trackEff.Draw()
canvas_trackEff.Draw()

canvas_algo = plt.Canvas(800,600)
h_algo.xaxis.title = 'trackType'
h_algo.yaxis.title = 'nTracks'
#h_algo.SetMarkerStyle(7)
h_algo.Draw()
canvas_algo.Draw()

canvas_jetRecoEff = plt.Canvas(800,600)
h_neutr_jet.Sumw2()
h_neutr_jet.SetStats(0)
h_neutr_jet.Divide(h_neutr)
h_neutr_jet.SetMarkerStyle(7)
h_neutr_jet.xaxis.title = 'pt'
h_neutr_jet.yaxis.title = 'jet Reco Efficiency'
h_neutr_jet.Draw("ep")  
canvas_jetRecoEff.Draw()



canvas_profileMass = plt.Canvas(800,600)
profile_mass.SetMarkerStyle(7)
profile_mass.xaxis.title = 'pt'
profile_mass.yaxis.title = 'mass resolution'
profile_mass.Draw('s')
canvas_profileMass.Draw()

canvas_profilePt = plt.Canvas(800,600)
profile_pt.SetMarkerStyle(7)
profile_pt.xaxis.title = 'pt'
profile_pt.yaxis.title = 'pt resolution'
profile_pt.Draw('s')
canvas_profilePt.Draw()

outfile.WriteTObject(canvas_algo, 'canvas_algo')
outfile.WriteTObject(canvas_profileMass, 'canvas_profileMass')
outfile.WriteTObject(canvas_profilePt, 'canvas_profilePt')
outfile.WriteTObject(canvas_trackEff, 'canvas_trackEff')
outfile.WriteTObject(canvas_jetRecoEff, 'canvas_jetRecoEff')
outfile.WriteTObject(canv_apVSpt, 'canv_apVSpt')
outfile.WriteTObject(canv_apVSpt_dM, 'canv_apVSpt_dM')
outfile.WriteTObject(canv_apVSpt_nodM, 'canv_apVSpt_nodM')
outfile.WriteTObject(graphQuant, 'graphQuant')
outfile.WriteTObject(graphMedian     , 'graphMedian'     )
outfile.WriteTObject(graphQuant_nodM , 'graphQuant_nodM' )
outfile.WriteTObject(graphMedian_nodM, 'graphMedian_nodM')
outfile.WriteTObject(graphQuant_dM   , 'graphQuant_dM'   )
outfile.WriteTObject(graphMedian_dM  , 'graphMedian_dM'  )
outfile.Write()



