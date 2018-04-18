import ROOT
from art.morisot import Morisot
from array import array
import sys,os
import glob

# Want to plot signal peaks overlapping with one another

# What are the signals?
massCouplingDict = {
"trijet" : {
            "p25": {"p1": {}, "p2": {}, "p3": {}, "p4": {}},
            "p35": {"p1": {}, "p2": {}, "p3": {}, "p4": {}},
            "p45": {"p1": {}, "p2": {}, "p3": {}, "p4": {}},
            "p55": {"p1": {}, "p2": {}, "p3": {}, "p4": {}}
           },
"dijetgamma" : {
#            "p" : []
           }
}

# What are our signal selections?
selectionDict = {
"trijet" : {
            "HLT_j380" : "j430_2j25", 
            "HLT_j225_gsc400_boffperf_split" : "j430_2j25"
           },
"dijetgamma" : {}
}

# Plot ranges differ between channels
rangeDict = {
"trijet" : {
            "xmin" : 200,
            "xmax" : 700,
            "ymin" : 1e-3,
            "ymax" : 1},
"dijetgamma" : {
            "xmin" : 100,
            "xmax" : 1200,
            "ymin" : "automatic",
            "ymax" : "automatic"}
}

# For labelling channels, etc
extraLinesDict = {
"trijet" : "X + j",
"dijetgamma" : "X + #gamma"
}

# Initialize painter
myPainter = Morisot()
#myPainter.setColourPalette("Tropical")
myPainter.setLabelType(2) # Sets label type i.e. Internal, Work in progress etc.
                          # See below for label explanation
# 0 Just ATLAS    
# 1 "Preliminary"
# 2 "Internal"
# 3 "Simulation Preliminary"
# 4 "Simulation Internal"
# 5 "Simulation"
# 6 "Work in Progress"

# Pick them up
inputDirs = "/afs/cern.ch/work/k/kpachal/DijetISR/Resolved2017/Local_Samples/Hists/OUT_MGPy8EG_N30LO_A14N23LO_dmA_jjj_Jet350_mR{0}_mD10_gS{1}_gD1_outTree/"

for channel in massCouplingDict.keys() :
  channelDict = massCouplingDict[channel]

  for mass in channelDict.keys() :
    massDict = channelDict[mass]

    for coupling in massDict.keys() :

      inputDir = inputDirs.format(mass,coupling)
      histfile = glob.glob(inputDir+"/hist-*.root")
      if len(histfile) > 1 :
        print "Expected only one file!" 
        print "Found:"
        print histfile
        sys.exit()

      openFile = ROOT.TFile.Open(histfile[0])

      # Different dirs.
      selection = selectionDict[channel]
      for trigger in selection.keys() :
        cuts = selection[trigger]
        dirName = "{0}_trig{1}_{2}".format(channel,trigger,cuts)
        hist = openFile.Get(dirName+"/Zprime_mjj_var")
        # Retrieve and normalise histogram
        hist.SetDirectory(0)
        hist.Scale(1.0/hist.Integral())
        massCouplingDict[channel][mass][coupling][trigger] = hist

      openFile.Close()


# Make the plots, divided by coupling, channel, and trigger, but not mass
for channel in massCouplingDict.keys() :

  couplingList = {}
  for mass in massCouplingDict[channel].keys() :
    couplings = massCouplingDict[channel][mass].keys() 
    for coupling in couplings :
      if not coupling in couplingList.keys() : 
        couplingList[coupling] = [mass]
      else :
        couplingList[coupling].append(mass)

  for coupling in couplingList.keys() :
    for trigger in selectionDict[channel].keys() :
      histList = []
      names = [] 

      # Put the histograms in the to-plot list, in order
      for availableMass in sorted(couplingList[coupling]) :
        histList.append(massCouplingDict[channel][availableMass][coupling][trigger])
        massString = availableMass.replace("p","0.")
        massGeV = eval(massString)*1e3
        massString = "{0}".format(int(massGeV))
        names.append("m_{0} = {1} GeV".format("{Z'}",massString))

      plotName = "signalshapes/signalShapes_{0}_gq{1}_trig_{2}".format(channel,coupling,trigger)

      # Collect plotting ranges
      xmin = rangeDict[channel]["xmin"]
      xmax = rangeDict[channel]["xmax"]
      ymin = rangeDict[channel]["ymin"]
      ymax = rangeDict[channel]["ymax"]

      extraLines = [extraLinesDict[channel]]

      myPainter.drawManyOverlaidHistograms(histList,names,"m_{jj} [GeV]","Norm. Events",plotName,histList[0].FindBin(xmin),histList[0].FindBin(xmax),ymin,ymax,extraLegendLines=extraLines,doLogX=False,doLogY=True,doErrors=False,doRectangular=False,doLegend=True,doLegendLow=False,doLegendLocation="Right",doLegendOutsidePlot=False,doATLASLabel="Low")
