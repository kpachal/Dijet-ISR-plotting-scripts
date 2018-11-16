#!/usr/bin/env python

import ROOT
import os,sys
#import pprint
from art.morisot import Morisot
from array import array
from pprint import pprint
from decimal import Decimal

SignalTitles = { "gSp1": "Z' (0.10)",
                 "gSp2": "Z' (0.20)",
                 "gSp3": "Z' (0.30)",
                 "gSp4": "Z' (0.40)",
                }

def rreplace(s, old, new, occurrence):
  li = s.rsplit(old, occurrence)
  return new.join(li)

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

# give it a file template which just needs the mass values to be filled in.
def getModelLimits(filetemplate, coupling, these_massvals, luminosity, limitsDictOut, \
  xname = "M_{Z'} [GeV]", yname = "#sigma #times #it{A} #times #epsilon #times BR [pb]" , label = "", makePlot = True):

  # Initialize painter
  myPainter = Morisot()
  myPainter.setColourPalette("ATLAS")
  myPainter.setEPS(True)
  myPainter.setLabelType(2) # Sets label type i.e. Internal, Work in progress etc.
                            # See below for label explanation
  # 0 Just ATLAS
  # 1 "Preliminary"
  # 2 "Internal"
  # 3 "Simulation Preliminary"
  # 4 "Simulation Internal"
  # 5 "Simulation"
  # 6 "Work in Progress"

  thisobserved = ROOT.TGraph()
  thisexpected = ROOT.TGraph()
  thisexpected_plus1  = ROOT.TGraph()
  thisexpected_minus1 = ROOT.TGraph()
  thisexpected_plus2  = ROOT.TGraph()
  thisexpected_minus2 = ROOT.TGraph()
  thistheory = ROOT.TGraph()
  for mass in these_massvals:

    import glob
    file_list = glob.glob(filetemplate.format(mass))
    print ("file searching format: ", filetemplate.format(mass))
    if len(file_list) == 0:
        print("no files found")
        continue
    allCLs = []
    PE_CLs = []
    for f in file_list:
      file = ROOT.TFile.Open(f)
      if not file or not file.Get("CLOfRealLikelihood"): continue
      CL = file.Get("CLOfRealLikelihood")[0]
      PE_tree = file.Get("ensemble_test")

      if not PE_tree or not CL:
          print("no pe")
          continue

      for event in PE_tree:
          PE_CLs.append( event.GetBranch("95quantile_marginalized_0").GetListOfLeaves().At(0).GetValue() )
      allCLs.append(CL)
    if len(allCLs) == 0:
        print("allCLs has a size of 0")
        continue
    expCLs = GetCenterAndSigmaDeviations(PE_CLs)
    print mass, allCLs[0], expCLs[2], len(PE_CLs)
    m = interpret_mass(mass)
    obsCL = allCLs[0]/luminosity
    print "obsCL: ", obsCL
    expCLs = [e/luminosity for e in expCLs]
    thisobserved.SetPoint(thisobserved.GetN(),m,obsCL)
    thisexpected_minus2.SetPoint(thisexpected_minus2.GetN(),m,expCLs[0])
    thisexpected_minus1.SetPoint(thisexpected_minus1.GetN(),m,expCLs[1])
    thisexpected.SetPoint(       thisexpected.GetN(),m,expCLs[2])
    thisexpected_plus1.SetPoint( thisexpected_plus1.GetN(),m,expCLs[3])
    thisexpected_plus2.SetPoint( thisexpected_plus2.GetN(),m,expCLs[4])

    print "mass:", m
    signal_info = {}
    # y hack
    #signal_acc         = signal_info['acc']
    #signal_thxsec      = signal_info['theory']
    signal_info['exp'] = expCLs[2]
    signal_info['obs'] = obsCL
    signal_info['exp+1'] = expCLs[3]
    signal_info['exp+2'] = expCLs[4]
    signal_info['exp-1'] = expCLs[1]
    signal_info['exp-2'] = expCLs[0]
    signal_info['nPEs'] = len(PE_CLs)
    if coupling not in limitsDictOut: limitsDictOut[coupling] = {}
    limitsDictOut[coupling]['%1.2f'%m] = signal_info

    #thistheory.SetPoint(thistheory.GetN(),m,signal_acc*signal_thxsec)

  if thisobserved.GetN() == 0:
    print "No limits found for %s"%signal
    print("f: ", f)
    #return limitsDictOut

  thisexpected1 = makeBand(thisexpected_minus1,thisexpected_plus1)
  thisexpected2 = makeBand(thisexpected_minus2,thisexpected_plus2)
  outputName = folderextension+"/Limits_"+label
  print "this expected:"
  thisexpected1.Print()

  xlow  = 'automatic'# (int(masses[signal][0]) - 100)/1000.
  xhigh = 'automatic'#(int(masses[signal][-1]) + 100)/1000.
  placeholder=ROOT.TGraph()
  placeholder.SetPoint(0,10,1e6)
  placeholder.SetPoint(1,1000, 1e6)

  if makePlot:
    ROOT.gROOT.SetBatch(True)
    myPainter.drawLimitSettingPlotObservedExpected(thisobserved,thisexpected,thisexpected1, thisexpected2, thistheory,SignalTitles[coupling],\
     outputName, xname,yname,luminosity,13,xlow,xhigh,2E-4,2E-2,False)

    thisobserved.Print()

  return limitsDictOut

def printPEs(limitsDict):
  for c,massdict in limitsDict.iteritems():
    print "Coupling = ", c
    mlist = []
    for m, siginfo in massdict.iteritems(): mlist.append( (m, siginfo['nPEs']))
    mlist = sorted(mlist, key=lambda e: e[0])
    print [m for m in mlist]
    print ["%g"%(float(m[0])*1000) for m in mlist if m[1] < 500]


def writeLimitsDict(folderextension, dataset,limitsDictOut):

  limitsDictFileName = folderextension+'LimitsDict_'+dataset+'.py'
  print 'Writing signal limit info to', limitsDictFileName
  limitsDictFile=open(limitsDictFileName, 'w')
  import pprint
  pp = pprint.PrettyPrinter(indent=4,stream=limitsDictFile)
  limitsDictFile.close()

def interpret_mass(mass) :

  # Please change this next iteration of the analysis
  # It only works for our signal grid.
  string = mass.replace("p",".")
  value = int(eval(string)*1000.0)
  print "Converted",mass,"to",value
  return value

if __name__ == "__main__":
  #====================================================================================
  # User Options
  #====================================================================================


  # Options
  folderextension = "./plots/"
  plotextension = "fullSysts"
  # make plots folder i.e. make folder extension
  if not os.path.exists(folderextension):
      os.makedirs(folderextension)

  #===============================================================================================

  # What is a file name template?
  individualLimitFiles = "setLimitsOneMassPoint_{0}_ZPrime_mZ{1}_{2}_seed*.root"

  # Dictionary to define everything I want to run on
  resultsDict = {
    "single_trigger_nbtag2" : {
      "indir" : "/home/kpachal/project/kpachal/Datasets_DijetISR/limit_results_zprime/2018.10.8/",
      "couplingMassMap" : {
        "gSp1" : ["p25","p35","p45","p55"],
        "gSp2" : ["p25","p35","p45","p55","p75"],
        "gSp3" : ["p25","p3","p35","p4","p45","p5","p55","p75","p95"]
      },
      "luminosity" : 79826
    },
    "compound_trigger_nbtag2" : {
      "indir" : "/home/kpachal/project/kpachal/Datasets_DijetISR/limit_results_zprime/2018.10.8/",
      "couplingMassMap" : {
        "gSp2" : ["p35","p45","p55","p75","p95"],
      },
      "luminosity" : 76595
    },
  }

  #====================================================================================

  # Loop over signals in dictionary
  limitsDictOut = {}
  signalDictOut = {}

  for spectrum in resultsDict.keys() :
  
    spectrum_info = resultsDict[spectrum]
    luminosity = spectrum_info["luminosity"]

    for coupling in spectrum_info["couplingMassMap"].keys() :
    
      masses_strings = spectrum_info["couplingMassMap"][coupling]
      input_dir = spectrum_info["indir"]
      myFileTemplate = input_dir + individualLimitFiles.format(spectrum,"{0}",coupling)
      label = spectrum + "_" + "coupling_"+coupling
  
      limitsDictOut = getModelLimits(myFileTemplate, coupling, masses_strings, luminosity, limitsDictOut, label=label, makePlot = True)

  writeLimitsDict(folderextension, dataset, limitsDictOut)

