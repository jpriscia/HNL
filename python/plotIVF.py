import ROOT
from ROOT import TColor
import rootpy
import rootpy.plotting as plt
from rootpy import stl
from rootpy.io import root_open
from rootpy.tree import Tree, TreeModel, FloatCol, IntCol
import pprint
from DataFormats.FWLite import Events, Handle
from ROOT import gSystem
from math import sqrt
from pdb import set_trace
from rootpy.vector import LorentzVector
from rootpy.plotting.style import get_style, set_style
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as mplt
from root_pandas import read_root
import pandas as pd
from pdb import set_trace
import numpy as np
import os

def eff(num, den):

    if den ==0:
        return(0.,1.)
    else: 
        efficiency = num/den

        if efficiency ==1:
            err_efficiency = 1
        else: 
            err_efficiency = sqrt(efficiency*(1-efficiency)/den)
    return (efficiency, err_efficiency)

style = get_style('ATLAS')
tsize=18
style.SetHistLineWidth(20);
style.SetLabelSize(tsize,"x");
style.SetTitleSize(tsize,"x");
style.SetLabelSize(tsize,"y");
style.SetTitleSize(tsize,"y");

set_style(style)

########################Define Histos##############################3
h_nPMuon_4GeV = plt.Hist(5, 0, 5, name="h_nPMuon_4GeV", title="h_nPMuon_4GeV",legendstyle='lep')
h_nTMuon_4GeV = plt.Hist(5, 0, 5, name="h_nTMuon_4GeV", title="h_nMuon_4GeV",legendstyle='lep')
h_nLMuon_4GeV = plt.Hist(5, 0, 5, name="h_nLMuon_4GeV", title="h_nLMuon_4GeV",legendstyle='lep')

h_nPMuon_7GeV = plt.Hist(5, 0, 5, name="h_nPMuon_7GeV", title="h_nPMuon_7GeV",legendstyle='lep')
h_nTMuon_7GeV = plt.Hist(5, 0, 5, name="h_nTMuon_7GeV", title="h_nTMuon_7GeV",legendstyle='lep')
h_nLMuon_7GeV = plt.Hist(5, 0, 5, name="h_nLMuon_7GeV", title="h_LPMuon_7GeV",legendstyle='lep')

h_nPMuon_10GeV = plt.Hist(5, 0, 5, name="h_nPMuon_10GeV", title="h_nPMuon_10GeV",legendstyle='lep')
h_nTMuon_10GeV = plt.Hist(5, 0, 5, name="h_nTMuon_10GeV", title="h_nTMuon_10GeV",legendstyle='lep')
h_nLMuon_10GeV = plt.Hist(5, 0, 5, name="h_nLMuon_10GeV", title="h_nLMuon_10GeV",legendstyle='lep')


h_3D_4GeV =  plt.Hist(100, 0, 5, name="h_3D_4GeV", title="h_3D_4GeV",legendstyle='lep')
h_3D_7GeV =  plt.Hist(100, 0, 5, name="h_3D_7GeV", title="h_3D_7GeV",legendstyle='lep')
h_3D_10GeV = plt.Hist(100, 0, 5, name="h_3D_10GeV", title="h_3D_10GeV",legendstyle='lep')

h_effBin_4GeV =  plt.Hist(5, 0, 100, name="h_effBin_4GeV", title="h_effBin_4GeV",legendstyle='lep')
h_effBin_7GeV =  plt.Hist(5, 0, 100, name="h_effBin_7GeV", title="h_effBin_7GeV",legendstyle='lep')
h_effBin_10GeV = plt.Hist(5, 0, 100, name="h_effBin_10GeV", title="h_effBin_10GeV",legendstyle='lep')

inputfile_4GeV = 'miniAOD_4GeV_IVF.root'
inputfile_7GeV = 'miniAOD_7GeV_IVF.root'
inputfile_10GeV = 'miniAOD_10GeV_IVF.root'

outfile = root_open('outputPlotIVF.root','recreate')
outfile.cd()


for inputF in [inputfile_4GeV,inputfile_7GeV,inputfile_10GeV]:
    inputfile = root_open(inputF)
    print '#######################################'
    print '#######################################'
    print 'input file is: ', inputfile

    cut_ivfLooseNotPoint = 'd3D_SecReco>0 && num_muonsLoose_noPrompt==1 && num_promptMuons>0 && sqrt((vz_SecReco-vz_PVReco)*(vz_SecReco-vz_PVReco)+(vy_SecReco-vy_PVReco)*(vy_SecReco-vy_PVReco)+(vx_SecReco-vx_PVReco)*(vx_SecReco-vx_PVReco))>0.05'
    cut_L_xyz_gen = 'sqrt((vx_SecGen-vx_NuGen)*(vx_SecGen-vx_NuGen)+(vy_SecGen-vy_NuGen)*(vy_SecGen-vy_NuGen)+(vz_SecGen-vz_NuGen)*(vz_SecGen-vz_NuGen))'


    #cut_ivfTightNotPoint = 'd3D_SecReco_T>0 && num_muonsTight_noPrompt==1 && sqrt((vz_SecReco_T-vz_PVReco)*(vz_SecReco_T-vz_PVReco)+(vy_SecReco_T-vy_PVReco)*(vy_SecReco_T-vy_PVReco)+(vx_SecReco_T-vx_PVReco)*(vx_SecReco_T-vx_PVReco))>0.05'
    
    eff_looseNotPoint_bin=[]
    
    bins=[0,20,40,60,80,100]
    for i in range(5):

        print 'entra'
        num = float(inputfile.gentree.GetEntries(cut_ivfLooseNotPoint+'&&'+ cut_L_xyz_gen +'>'+str(bins[i]) + '&&' + cut_L_xyz_gen +'<'+str(bins[i+1]) ))
        den = float(inputfile.gentree.GetEntries('num_muonsLoose_noPrompt==1 && num_promptMuons>0'+'&&'+ cut_L_xyz_gen +'>'+str(bins[i]) + '&&' + cut_L_xyz_gen +'<'+str(bins[i+1])))
        print num, den
        eff_looseNotPoint_bin.append(eff(num,den))
        print eff_looseNotPoint_bin[i]
        if inputF==inputfile_4GeV:   
            h_effBin_4GeV[i].value = eff_looseNotPoint_bin[i][0]   
            h_effBin_4GeV[i].error = eff_looseNotPoint_bin[i][1]
        if inputF==inputfile_7GeV:   
            h_effBin_7GeV[i].value = eff_looseNotPoint_bin[i][0] 
            h_effBin_7GeV[i].error = eff_looseNotPoint_bin[i][1]
        if inputF==inputfile_10GeV:   
            h_effBin_10GeV[i].value = eff_looseNotPoint_bin[i][0]
            h_effBin_10GeV[i].error = eff_looseNotPoint_bin[i][1]

    all_ev = float(inputfile.gentree.GetEntries())
    in_acc =  float(inputfile.gentree.GetEntries('vz_PVReco>-1000'))
    promptMuon =  float(inputfile.gentree.GetEntries('num_promptMuons>0'))
    num_twoTight =  float(inputfile.gentree.GetEntries('num_muonsTight_noPrompt==1 && num_promptMuons>0'))
    num_LooseTight =  float(inputfile.gentree.GetEntries('num_muonsLoose_noPrompt==1 && num_promptMuons>0'))
    num_ivfLoose =  float(inputfile.gentree.GetEntries('d3D_SecReco>0 && num_muonsLoose_noPrompt==1 &&  && d3D_SecReco<0.2'))
    num_ivfTight =  float(inputfile.gentree.GetEntries('d3D_SecReco_T>0 && num_muonsTight_noPrompt==1 &&  && d3D_SecReco_T<0.2'))
    num_ivfTightNotPoint = float(inputfile.gentree.GetEntries('d3D_SecReco_T>0 && num_muonsTight_noPrompt==1 && sqrt((vz_SecReco_T-vz_PVReco)*(vz_SecReco_T-vz_PVReco)+(vy_SecReco_T-vy_PVReco)*(vy_SecReco_T-vy_PVReco)+(vx_SecReco_T-vx_PVReco)*(vx_SecReco_T-vx_PVReco))>0.05 && d3D_SecReco_T<0.2'))
    num_ivfLooseNotPoint = float(inputfile.gentree.GetEntries('d3D_SecReco>0 && num_muonsLoose_noPrompt==1 && sqrt((vz_SecReco-vz_PVReco)*(vz_SecReco-vz_PVReco)+(vy_SecReco-vy_PVReco)*(vy_SecReco-vy_PVReco)+(vx_SecReco-vx_PVReco)*(vx_SecReco-vx_PVReco))>0.05  && d3D_SecReco<0.2'))

    print 'acceptance: ', eff(in_acc,all_ev)
    print '1 prompt muon efficiency is: ', eff(promptMuon,in_acc)
    print 'one loose one prompt muon eff is: ', eff(num_LooseTight,promptMuon)
    print 'two tight muons (one prompt) efficiency is: ', eff(num_twoTight,promptMuon)
    print 'number of verteces loose: ',  eff(num_ivfLoose,num_LooseTight)
    print 'number of verteces tight: ', eff(num_ivfTight,num_twoTight)
    print 'number of verteces loose with PV cut: ',  eff(num_ivfLooseNotPoint,num_LooseTight)
    print 'number of verteces tight with PV cut: ', eff(num_ivfTightNotPoint,num_twoTight)

    for entry in inputfile.gentree:
        L_z_gen = abs(entry.vz_SecGen-entry.vz_NuGen)
        L_xy_gen = sqrt((entry.vx_SecGen-entry.vx_NuGen)*(entry.vx_SecGen-entry.vx_NuGen)+(entry.vy_SecGen-entry.vy_NuGen)*(entry.vy_SecGen-entry.vy_NuGen))
        L_xyz_gen = sqrt((entry.vx_SecGen-entry.vx_NuGen)*(entry.vx_SecGen-entry.vx_NuGen)+(entry.vy_SecGen-entry.vy_NuGen)*(entry.vy_SecGen-entry.vy_NuGen)+(entry.vz_SecGen-entry.vz_NuGen)*(entry.vz_SecGen-entry.vz_NuGen))

        if entry.vz_PVReco > -1000.: 
        #plt.hist(df_PV['num_promptMuons'], bins=5, normed=True, range=(0,5), histtype='',label='number of prompt muons') 
            if inputF==inputfile_4GeV:  h_nPMuon_4GeV.fill(entry.num_promptMuons)
            if inputF==inputfile_7GeV:  h_nPMuon_7GeV.fill(entry.num_promptMuons)
            if inputF==inputfile_10GeV:  h_nPMuon_10GeV.fill(entry.num_promptMuons)

        if entry.num_promptMuons >0:
            if inputF==inputfile_4GeV:  
                h_nLMuon_4GeV.fill(entry.num_muonsLoose_noPrompt)
                h_nTMuon_4GeV.fill(entry.num_muonsTight_noPrompt)
            if inputF==inputfile_7GeV: 
                h_nLMuon_7GeV.fill(entry.num_muonsLoose_noPrompt)
                h_nTMuon_7GeV.fill(entry.num_muonsTight_noPrompt)
            if inputF==inputfile_10GeV: 
                h_nLMuon_10GeV.fill(entry.num_muonsLoose_noPrompt)
                h_nTMuon_10GeV.fill(entry.num_muonsTight_noPrompt)


        if entry.num_muonsTight_noPrompt==1 and entry.num_promptMuons>0 and entry.d3D_SecReco_T>0:

            deltaZ_T = entry.vz_SecReco_T-entry.vz_PVReco
            deltaX_T = entry.vx_SecReco_T-entry.vx_PVReco
            deltaY_T = entry.vz_SecReco_T-entry.vz_PVReco
            L_z_recoT = abs(deltaZ_T)
            L_xy_recoT = sqrt(deltaX_T*deltaX_T+deltaY_T*deltaY_T)
            L_yxz_recoT = sqrt(deltaX_T*deltaX_T+deltaY_T*deltaY_T+deltaZ_T*deltaZ_T)


        if entry.num_muonsLoose_noPrompt==1 and entry.num_promptMuons>0 and entry.d3D_SecReco>0:

            deltaZ = entry.vz_SecReco-entry.vz_PVReco
            deltaX = entry.vx_SecReco-entry.vx_PVReco
            deltaY = entry.vz_SecReco-entry.vz_PVReco
            L_z_reco = abs(deltaZ)
            L_xy_reco = sqrt(deltaX*deltaX+deltaY*deltaY)
            L_yxz_reco = sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ)

            if inputF==inputfile_4GeV:  h_3D_4GeV.fill(entry.d3D_SecReco)
            if inputF==inputfile_7GeV:  h_3D_7GeV.fill(entry.d3D_SecReco)
            if inputF==inputfile_10GeV: h_3D_10GeV.fill(entry.d3D_SecReco)


canvas_3Ddist=plt.Canvas(800,600) 
h_3D_4GeV.markercolor = 'crimson'                                                                                                                                          
h_3D_7GeV.markercolor = 'darkblue' 
h_3D_10GeV.markercolor = 'green'
h_3D_4GeV.Scale(1/h_3D_4GeV.Integral()) 
h_3D_7GeV.Scale(1/h_3D_7GeV.Integral())
h_3D_10GeV.Scale(1/h_3D_10GeV.Integral())
h_3D_4GeV.xaxis.title = 'd_xyz'                                                                                                                                        
h_3D_4GeV.yaxis.title = 'A.U.'  
h_3D_4GeV.Draw()
h_3D_7GeV.Draw('SAME')
h_3D_10GeV.Draw('SAME')
legend = plt.Legend([h_3D_4GeV,h_3D_7GeV,h_3D_10GeV], leftmargin=0.45, margin=0.3,textsize=20)    
legend.Draw()
canvas_3Ddist.SetLogy()
canvas_3Ddist.Draw()


canvas_effBin=plt.Canvas(800,600)
h_effBin_4GeV.markercolor = 'crimson'
h_effBin_7GeV.markercolor = 'darkblue'
h_effBin_10GeV.markercolor = 'green'
h_effBin_4GeV.xaxis.title = 'd_xyz'
h_effBin_4GeV.yaxis.title = 'A.U.'
h_effBin_4GeV.SetMaximum(1.2)
#h_effBin_4GeV.SetMarkerStyle(7)
h_effBin_4GeV.Draw()
#h_effBin_7GeV.Draw('SAME')
#h_effBin_10GeV.Draw('SAME')
#legend = plt.Legend([h_effBin_4GeV,h_effBin_7GeV,h_effBin_10GeV], leftmargin=0.45, margin=0.3,textsize=20)
#legend.Draw()
canvas_effBin.Draw()



canvas_muon4GeV=plt.Canvas(800,600)
canvas_muon4GeV.Divide(3,1)
canvas_muon4GeV.cd(1)
h_nPMuon_4GeV.xaxis.title = 'N prompt muons' 
h_nPMuon_4GeV.xaxis.title = 'N events'
h_nPMuon_4GeV.Draw()
canvas_muon4GeV.cd(2)
h_nLMuon_4GeV.xaxis.title = 'N loose muons'
h_nLMuon_4GeV.xaxis.title = 'N events'
h_nLMuon_4GeV.Draw()
canvas_muon4GeV.cd(3)
h_nTMuon_4GeV.xaxis.title = 'N tight muons'
h_nTMuon_4GeV.xaxis.title = 'N events'
h_nTMuon_4GeV.Draw()

canvas_muon7GeV=plt.Canvas(800,600)
canvas_muon7GeV.Divide(3,1)
canvas_muon7GeV.cd(1)
h_nPMuon_7GeV.xaxis.title = 'N prompt muons'
h_nPMuon_7GeV.xaxis.title = 'N events'
h_nPMuon_7GeV.Draw()
canvas_muon7GeV.cd(2)
h_nLMuon_7GeV.xaxis.title = 'N loose muons'
h_nLMuon_7GeV.xaxis.title = 'N events'
h_nLMuon_7GeV.Draw()
canvas_muon7GeV.cd(3)
h_nTMuon_7GeV.xaxis.title = 'N tight muons'
h_nTMuon_7GeV.xaxis.title = 'N events'
h_nTMuon_7GeV.Draw()


canvas_muon10GeV=plt.Canvas(800,600)
canvas_muon10GeV.Divide(3,1)
canvas_muon10GeV.cd(1)
h_nPMuon_10GeV.xaxis.title = 'N prompt muons'
h_nPMuon_10GeV.xaxis.title = 'N events'
h_nPMuon_10GeV.Draw()
canvas_muon10GeV.cd(2)
h_nLMuon_10GeV.xaxis.title = 'N loose muons'
h_nLMuon_10GeV.xaxis.title = 'N events'
h_nLMuon_10GeV.Draw()
canvas_muon10GeV.cd(3)
h_nTMuon_10GeV.xaxis.title = 'N tight muons'
h_nTMuon_10GeV.xaxis.title = 'N events'
h_nTMuon_10GeV.Draw()




outfile.WriteTObject(canvas_3Ddist  , 'canvas_3Ddist'  )
outfile.WriteTObject(canvas_effBin  , 'canvas_effBin'  )
outfile.WriteTObject(canvas_muon4GeV  , 'canvas_muon4GeV'  )
outfile.WriteTObject(canvas_muon7GeV  , 'canvas_muon7GeV'  )
outfile.WriteTObject(canvas_muon10GeV  , 'canvas_muon10GeV'  )
outfile.Write()


