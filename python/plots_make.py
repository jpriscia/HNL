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

inputfile = root_open('tree_lljj.root')

outfile = root_open('outputPlot.root','recreate')
outfile.cd()

#HPS
#shrinking cone
#apertura vs pt neutrinone 
#prendere quantile al 90% e farne un fit

h_apVSpt_dM = plt.Hist2D(100, 0, 100, 100, 0, 3, name="aperture2D_dM", title="aperture vs nu_pt (with D Meson) -- gen level")
h_apVSpt_nodM = plt.Hist2D(100, 0, 100, 100, 0, 3, name="aperture2D_nodM", title="aperture vs nu_pt (no D Meson) -- gen level")
h_apVSpt = plt.Hist2D(100, 0, 100, 100, 0, 3, name="aperture2D", title="aperture vs nu_pt -- gen level")

h_trackEff = plt.Hist(300, 0, 0.3, name="trackEff", title="track matching")

for entry in inputfile.gentree:
    h_apVSpt.fill(entry.pt_NuGen, entry.aperture_Gen)
    for val in entry.dr:
        h_trackEff.fill(val)
 
    if entry.dMeson !=0:
        h_apVSpt_dM.fill(entry.pt_NuGen, entry.aperture_Gen)

    if entry.dMeson == 0:
        h_apVSpt_nodM.fill(entry.pt_NuGen, entry.aperture_Gen)

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
        #set_trace()
        points = zip([i for i,_ in values], accumulation)
        #set_trace()
        median, _ = min(points, key=lambda x: abs(x[1] - integral/2.))
        #set_trace()
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
graphQuant.Draw('L SAME')
graphMedian.Draw('L SAME')
canv_apVSpt.Draw()

canv_apVSpt_nodM = plt.Canvas(800,600)
h_apVSpt_nodM.xaxis.title = 'nu_pt'
h_apVSpt_nodM.yaxis.title = 'aperture'
h_apVSpt_nodM.Draw('COLZ')
graphQuant_nodM.Draw('L SAME')
graphMedian_nodM.Draw('L SAME')
canv_apVSpt_nodM.Draw()

canv_apVSpt_dM = plt.Canvas(800,600)
h_apVSpt_dM.xaxis.title = 'nu_pt'
h_apVSpt_dM.yaxis.title = 'aperture'
h_apVSpt_dM.Draw('COLZ')
graphQuant_dM.Draw('LP SAME')
graphMedian_dM.Draw('LP SAME')
canv_apVSpt_dM.Draw()

h_trackEff.xaxis.title = 'dr'
h_trackEff.yaxis.title = 'nTracks'

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



