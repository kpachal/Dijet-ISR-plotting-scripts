
#/usr/bin/env python

import os
import sys
import ROOT
from art.morisot import Morisot
from itertools import repeat


names = { '0.00':  '#sigma_{G}/m_{G} = Res.',
          '0.05': '#sigma_{G}/m_{G} = 0.05',
          '0.07': '#sigma_{G}/m_{G} = 0.07',
          '0.10': '#sigma_{G}/m_{G} = 0.10',
          '0.15': '#sigma_{G}/m_{G} = 0.15'}

# Initialize painter
myPainter = Morisot()
# Internal
myPainter.setColourPalette("Tropical")
myPainter.setLabelType(2)

def makeBand(graph1, graph2):
  points = []
  for i in range(graph1.GetN()):
    points += [(i,graph1.GetX()[i],graph1.GetY()[i])]
  for i in range(graph2.GetN()-1,-1,-1):
    points += [(i,graph2.GetX()[i],graph2.GetY()[i])]
  graph_band = ROOT.TGraph();
  for i in range (len(points)): graph_band.SetPoint(i,points[i][1],points[i][2])
  return graph_band

def GetCenterAndSigmaDeviations(inputs) :
  inputs = sorted(inputs)
  statVals = []
  quantiles = [0.02275,0.1587,0.5,0.8413,0.9772]
  for q in quantiles:
    wantEvents = len(inputs)*q
    statVals.append(inputs[int(wantEvents)])
  return statVals
  
def getGaussianLimits( inputfileform, ratios, dataset, luminosity, cutstring, makePlot = True, outfolder = "./plots/"):

  basicInputFiles = {}
  for r in ratios: basicInputFiles[r] = []

  import glob
  file_list = glob.glob(inputfileform)
  print "searching for files ",inputfileform
  print "found",file_list
  for f in file_list:
    print f
    ratio_str = f.split('.')[-2].split('_')[-1]
    print ratio_str
    if ratio_str in basicInputFiles.keys(): basicInputFiles[ratio_str].append(f)

  myPainter.cutstring = cutstring
  #myPainter.setLabelType(1)

  minMassVal = {}
  values = {}
  mass_list = []
  allobserved = []
  allexpected = []
  allexpected1Sigma = []
  allexpected2Sigma = []
  
  results = {}
  
  for r in ratios:
    values[r] = {}
    print basicInputFiles[r]
    for f in basicInputFiles[r]:
      file = ROOT.TFile.Open(f)
      CLs_str = "CLsPerMass_widthToMass{0}".format(r)
      if not file or not file.Get(CLs_str): continue
      CLs = file.Get(CLs_str) 
      masses = file.Get("massesUsed")

      if masses == None: continue

      for i,mass in enumerate(masses) :
        
          mass_list += [mass]
          if mass not in values[r]:
            values[r][mass] = {'obs': [], 'exp': [], 'PEs': [] }
          values[r][mass]['obs'].append(CLs[i]/luminosity)
          
          if not doExpected :
            continue
          
          PE_tree = file.Get("ensemble_tree_{0}_{1}".format(mass,r))
          PE_CLs = []
          for event in PE_tree:
              PE_CLs.append( event.GetBranch("95quantile_marginalized_0").GetListOfLeaves().At(0).GetValue() )
          expCLs = GetCenterAndSigmaDeviations(PE_CLs)
          values[r][mass]['exp'].append(expCLs[2]/luminosity)
          values[r][mass]['PEs'] += [e/luminosity for e in PE_CLs]        

    mass_list = sorted(list(set(mass_list)))
    thisobserved = ROOT.TGraph()
    thisexpected = ROOT.TGraph()
    thisexpected_plus1  = ROOT.TGraph()
    thisexpected_minus1 = ROOT.TGraph()
    thisexpected_plus2  = ROOT.TGraph()
    thisexpected_minus2 = ROOT.TGraph()
    for m in mass_list :
      if m not in values[r]: continue
      thisobserved.SetPoint(       thisobserved.GetN(),m,values[r][m]['obs'][0])
      if not doExpected :
        continue
      
      expCLs = GetCenterAndSigmaDeviations(values[r][m]['PEs'])
      print r, m, values[r][m]['obs'][0], values[r][m]['exp'][0], len(values[r][m]['PEs'])
      thisexpected_minus2.SetPoint(thisexpected_minus2.GetN(),m,expCLs[0])
      thisexpected_minus1.SetPoint(thisexpected_minus1.GetN(),m,expCLs[1])
      thisexpected.SetPoint(       thisexpected.GetN(),m,expCLs[2])
      thisexpected_plus1.SetPoint( thisexpected_plus1.GetN(),m,expCLs[3])
      thisexpected_plus2.SetPoint( thisexpected_plus2.GetN(),m,expCLs[4])

    allobserved.append(thisobserved)
    allexpected.append(thisexpected)
    
    thisexpected1Sigma = makeBand(thisexpected_minus1,thisexpected_plus1)
    thisexpected2Sigma = makeBand(thisexpected_minus2,thisexpected_plus2)
    
    if r == 0:
      allexpected1Sigma.append(thisexpected1Sigma)
      allexpected2Sigma.append(thisexpected2Sigma)
      #limtxt = ' #scale[0.8]{(95% CL U.L. #pm 1-2#sigma)}'
      #if limtxt not in names['%1.2f'%r]:
      #  names['%1.2f'%r] += limtxt
    else:
        allexpected1Sigma.append(ROOT.TGraph())
        allexpected2Sigma.append(ROOT.TGraph()) 
    #allexpected1Sigma.append(thisexpected1Sigma)
    #allexpected2Sigma.append(thisexpected2Sigma)
    results[r] = {'obs':thisobserved, 'exp': thisexpected,'exp1':  thisexpected1Sigma,'exp2': thisexpected2Sigma}
    
  if makePlot:
    #print ratios
    #print [names['%1.2f'%r] for r in ratios]
    outname = outfolder+"GenericGaussians_"+channel

    if len(ratios) == 1:
      outname += "_" + str(int(ratios[0]*100))
    myPainter.drawSeveralObservedExpectedLimits(allobserved,allexpected,allexpected1Sigma,allexpected2Sigma,[names['%1.2f'%r] for r in ratios],outname,"m_{G} [GeV]",\
     "#sigma #times #it{A} #times BR [pb]",luminosity,13,400,1850,2E-2,200,[],ATLASLabelLocation="BottomL",cutLocation="Left", doLegendLocation="Left" if len(ratios) == 1 else "Center")
  
  return results

if __name__ == "__main__":

  #==========================
  #   User configurables
  #==========================

  makeCombinedPlot = False
  
  luminosity = {"incSing" : 79826, "incComp" : 76595, "2bSing" : 79826, "2bComp" : 76595}
  doExpected = False

  indir = "/cluster/warehouse/kpachal/DijetISR/Resolved2017/LimitSetting/BayesianFramework/results/gaussian_limits/"
  dirsToCompare = ["katetest_noSysts","katetest_scalesOnly","katetest_JESOnly","katetest_allSignal"]
  channelsToCheck = ["incSing","incComp","2bSing","2bComp"]
  names = ["No syst","Only scale systs","Only JES","All signal systs","All systs"]
  ratios = ["100"]
  
  yrange = {"incSing" : [3E-3,1E-1],
            "incComp" : [3E-3,1E-1],
            "2bSing" : [3E-4,1E-2],
            "2bComp" : [3E-4,1E-2]}
  
  for channel in channelsToCheck :
    for ratio in ratios :
  
      graph_list = []
  
      for dir in dirsToCompare :
        infileform = indir + dir + "/GenericGaussians*{0}*".format(channel)
        interpreted_files = getGaussianLimits( infileform, [ratio], channel, luminosity[channel], "|y*| < 0.75", False)
        graph_list.append(interpreted_files[ratio]["obs"])

      # Now make plot
      plotname = "compare_validation/GenericGaussians_compare_{0}_{1}".format(channel,ratio)
      myPainter.drawSeveralObservedLimits(graph_list,names,plotname,"m_{G}","#sigma x A x #epsilon",luminosity[channel],13,"automatic","automatic",yrange[channel][0],yrange[channel][1],extraLegendLines = [], doLogY=True,doLogX=False,doRectangular=False,doLegendLocation="Right",ATLASLabelLocation="BottomL",isTomBeingDumb=False,addHorizontalLines=[],pairNeighbouringLines=False,cutLocation="Right")

  #==========================
  
 
