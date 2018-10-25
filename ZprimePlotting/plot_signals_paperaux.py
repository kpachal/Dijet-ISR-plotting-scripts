import ROOT
from art.morisot import Morisot
from array import array
import sys,os
import glob

# Want to plot signal peaks overlapping with one another

# What are the signals?
sigsToPlot = {
  "single_trigger" : {
    "couplings" : ["p1","p2"],
    "masses" : ["p25","p35","p45"],
    "filter" : "100"
  },
  "compound_trigger" : {
    "couplings" : ["p2"],
    "masses" : ["p45","p55","p75","p95"],
    "filter" : "50"
  }
}


# Plot ranges differ between channels
rangeDict = {
"single_trigger" : {
            "xmin" : 169,
            "xmax" : 600,
            "ymin" : "automatic",
            "ymax" : "automatic"},
"compound_trigger" : {
            "xmin" : 350,
            "xmax" : 1200,
            "ymin" : "automatic",
            "ymax" : "automatic"}
}

labelDict = {
"single_trigger" : ["Single photon trigger,", "inclusive selection"],
"compound_trigger" : ["Single photon trigger,", "2 #it{b}-tags"],
}

# Initialize painter
myPainter = Morisot()
myPainter.setColourPalette("notSynthwave")
myPainter.setLabelType(4) # Sets label type i.e. Internal, Work in progress etc.
                          # See below for label explanation
# 0 Just ATLAS    
# 1 "Preliminary"
# 2 "Internal"
# 3 "Simulation Preliminary"
# 4 "Simulation Internal"
# 5 "Simulation"
# 6 "Work in Progress"

# Pick them up
inputDir = "/home/kpachal/project/kpachal/Datasets_DijetISR/limit_setting_sigsamples"

for channel in sigsToPlot.keys() :
  channelDict = sigsToPlot[channel]

  for coupling in channelDict["couplings"] :
  
    massList = channelDict["masses"]
    histList = []
    namesList = []

    for mass in massList :

      histfile = inputDir+"/Zprime_{0}_inclusive_MGPy8EG_N30LO_A14N23LO_dmA_jja_Ph{1}_mR{2}_mD10_gS{3}_gD1.root".format(channel,channelDict["filter"],mass,coupling)

      openFile = ROOT.TFile.Open(histfile)

      hist = openFile.Get("m_jj_Nominal")

      # Retrieve (and normalise?) histogram
      hist.SetDirectory(0)
      hist.SetName("m_jj_Nominal_{0}_{1}".format(mass,coupling))
      # Normalise to 1 fb
      hist.Scale(1000.0)
      #hist.Scale(1.0/hist.Integral())
      openFile.Close()
      histList.append(hist)
      
      # Make name for this
      massVal = eval(mass.replace("p","."))*1000.0
      massString = "{0}".format(int(massVal))
      namesList.append("m_{0} = {1} GeV".format("{Z'}",massString))

    # Make one plot per coupling and trigger
    extraLegendLine = list(reversed(labelDict[channel]))
    plotName = "signalshapes/signalShapes_{0}_gq{1}".format(channel,coupling)

    # Collect plotting ranges
    xmin = rangeDict[channel]["xmin"]
    xmax = rangeDict[channel]["xmax"]
    ymin = rangeDict[channel]["ymin"]
    ymax = rangeDict[channel]["ymax"]
    print xmin,xmax,histList[0].FindBin(xmin),histList[0].FindBin(xmax)

    myPainter.drawManyOverlaidHistograms(histList,namesList,"m_{jj} [GeV]","Events/fb",plotName,histList[0].FindBin(xmin),histList[0].FindBin(xmax),ymin,ymax,extraLegendLines=extraLegendLine,doLogX=False,doLogY=False,doErrors=False,doRectangular=False,doLegend=True,doLegendLow=False,doLegendLocation="Right",doLegendOutsidePlot=False,doATLASLabel="Left")

    myPainter.drawManyOverlaidHistograms(histList,namesList,"m_{jj} [GeV]","Events/fb",plotName+"_logx",histList[0].FindBin(xmin),histList[0].FindBin(xmax),ymin,ymax,extraLegendLines=extraLegendLine,doLogX=False,doLogY=True,doErrors=False,doRectangular=False,doLegend=True,doLegendLow=False,doLegendLocation="Right",doLegendOutsidePlot=False,doATLASLabel="Low")
