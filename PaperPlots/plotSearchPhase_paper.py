#!/usr/bin/env python

import os
import ROOT
from art.morisot import Morisot
from analysisScripts.searchphase import searchFileData
import analysisScripts.generalfunctions
from array import array
import sys
import math


# Initialize painter
myPainter = Morisot()
myPainter.setColourPalette("notSynthwave")
myPainter.setLabelType(2)
myPainter.setEPS(True)
myPainter.dodrawUsersText = True

CME = 13000

doIndividualPlots = False

class multiPlotter :

  def __init__(self,lower_mass_data,higher_mass_data,permitWindow=False) :

    self.data1 = lower_mass_data
    self.data2 = higher_mass_data

  def makeDoubleSearchPhasePlot(self,low1,low2,high,folder,ext,channelLabel="",dataPointsOption=0,fancinessOption=0,includeTriggers = False) :

    firstBin = self.data1.basicData.FindBin(low1)
    lastBin = self.data1.basicData.FindBin(high)
    
    extraLines = []
    if channelLabel :
      extraLines = [channelLabel]
    
    if includeTriggers :
      datastring1 = "single-photon"
      datastring2 = "combined"
    else :
      datastring1 = ""
      datastring2 = ""

    # Fancy figure 1 (two spectra)
    myPainter.drawDataAndFitsOverSignificanceHists_TwoSpectra(\
        self.data1.basicData, self.data1.basicBkgFromFit,self.data1.residualHist,\
        self.data2.basicData, self.data2.basicBkgFromFit,self.data2.residualHist,\
        self.data1.sigHist, self.data2.sigHist, self.data1.sigScale, self.data2.sigScale,\
         'm_{jj} [GeV]','Events','Significance','{0}/figure1'.format(folder)+ext,\
        self.data1.luminosity, self.data2.luminosity, datastring1, datastring2, 13, low1, high,\
        firstBin, lastBin, True, self.data1.bumpLowEdge, self.data1.bumpHighEdge,\
        self.data2.bumpLowEdge, self.data2.bumpHighEdge, extraLegendLines=extraLines, \
        pval = self.data1.bumpHunterPVal, chi2pval = self.data1.chi2PVal,\
        pval2 = self.data2.bumpHunterPVal, chi2pval2 = self.data2.chi2PVal, doRectangular=False,
        dataPointsOption=dataPointsOption,fancinessOption=fancinessOption,extraSpace=("tag" in channelLabel))

# Pass lumi in pb
def ProcessSignalFile(sigfile,lumi,scale,binLow,binHigh,compare_data) :

  openSigFile = ROOT.TFile.Open(sigfile,"READ")
  sig = openSigFile.Get("m_jj_Nominal")
  sig.SetDirectory(0)
  openSigFile.Close()
  sig.Scale(lumi)
    
  for ibin in xrange(1+sig.GetNbinsX()):
    if ibin < binLow or ibin > binHigh:
      sig.SetBinContent(ibin,0)
    if sig.GetBinContent(ibin) > 0:
      sig.SetBinContent(ibin,scale*sig.GetBinContent(ibin)+compare_data.basicData.GetBinContent(ibin))
    if compare_data.basicData.GetBinContent(ibin) > 0:
      # This was 0.0022 for second, initially.
      if sig.GetBinContent(ibin)/compare_data.basicData.GetBinContent(ibin) < 0.0001:
        sig.SetBinContent(ibin,0)

  return sig

#####################################################
# User controlled
folder = "plots/"
filePath = "/home/kpachal/project/kpachal/Datasets_DijetISR/search_official_results/"
signalFileTemplate = "/home/kpachal/project/kpachal/Datasets_DijetISR/signal_inputs/Zprime_{0}_{1}_MGPy8EG_N30LO_A14N23LO_dmA_jja_Ph100_mRp{2}_mD10_gSp{3}_gD1.root"

#make plots folder i.e. make folder extension
if not os.path.exists(folder):
  os.makedirs(folder)

lumiDict = { "single_trigger" : 79826,
             "compound_trigger" : 76595
            }

signalProcessingInfo = {
  "single_trigger" : {
    "inclusive" : {"mass" : "25",
                   "scale" : 150,
                   "binLow" : 3,
                   "binHigh" : 15},
    "nbtag2"      : {"mass" : "25",
                   "scale" : 15,
                   "binLow" : 3,
                   "binHigh" : 15}
  },
  "compound_trigger" : {
    "inclusive" : { "mass" : "55",
                    "scale" : 150,
                    "binLow" : 6,
                    "binHigh" : 25},
    "nbtag2"      : { "mass" : "55",
                    "scale" : 15,
                    "binLow" : 6,
                    "binHigh" : 25}
  }
}

channelLabels = {
"inclusive" : "Flavour inclusive",
"nbtag2" : "2 b-tags"
}

pValCutoff = 0.05
low1 = 169
low2 = 335
high = 1200

for channel in ["inclusive","nbtag2"] :

  channel_data = {}

  for trigger in ["single_trigger", "compound_trigger"] :
  
    fileTemplate = filePath + "SearchPhase_dijetgamma_{0}_{1}.root".format(trigger,channel)
    
    # Get search file data
    data = searchFileData(fileTemplate)

    # Get processing info for signal
    sig_data = signalProcessingInfo[trigger][channel]
    
    # Append luminosity info
    lumi = lumiDict[trigger]
    data.luminosity = lumi
    
    # Append scale info
    scale = sig_data["scale"]
    sig_binLow = sig_data["binLow"]
    sig_binHigh = sig_data["binHigh"]
    data.sigScale = scale
    
    # Process signal file and append its info
    signalFile = signalFileTemplate.format(trigger,channel,sig_data["mass"],"1")
    sigHist = ProcessSignalFile(signalFile,lumi,scale,sig_binLow,sig_binHigh,data)
    data.sigHist = sigHist
  
    channel_data[trigger] = data

  extension = "_"+channel
  label = channelLabels[channel]

  # Change dataPointsOption for more variations

  # Baseline
  paperPlotter = multiPlotter(channel_data["single_trigger"],channel_data["compound_trigger"],True)
  paperPlotter.makeDoubleSearchPhasePlot(low1,low2,high,folder,extension, label, dataPointsOption=2)

  extension_variations = extension+"_fancinessOption1"
  paperPlotter.makeDoubleSearchPhasePlot(low1,low2,high,folder,extension_variations, label, dataPointsOption=2, fancinessOption=1)

  extension_variations = extension+"_fancinessOption2"
  paperPlotter.makeDoubleSearchPhasePlot(low1,low2,high,folder,extension_variations, label, dataPointsOption=2, fancinessOption=2)

  # Including trigger names
  extension_variations = extension+"_withTrigLabels"
  paperPlotter = multiPlotter(channel_data["single_trigger"],channel_data["compound_trigger"],True)
  paperPlotter.makeDoubleSearchPhasePlot(low1,low2,high,folder,extension_variations, label, dataPointsOption=2,includeTriggers = True)

  extension_variations = extension+"_withTrigLabels_fancinessOption1"
  paperPlotter.makeDoubleSearchPhasePlot(low1,low2,high,folder,extension_variations, label, dataPointsOption=2, fancinessOption=1,includeTriggers = True)

  extension_variations = extension+"_withTrigLabels_fancinessOption2"
  paperPlotter.makeDoubleSearchPhasePlot(low1,low2,high,folder,extension_variations, label, dataPointsOption=2, fancinessOption=2,includeTriggers = True)

  # Including trigger names and y* cuts
  extension_variations = extension+"_withTrigLabelsAndYstar"
  paperPlotter = multiPlotter(channel_data["single_trigger"],channel_data["compound_trigger"],True)
  paperPlotter.makeDoubleSearchPhasePlot(low1,low2,high,folder,extension_variations, label+", |y*| < 0.75", dataPointsOption=2,includeTriggers = True)

  extension_variations = extension+"_withTrigLabelsAndYstar_fancinessOption1"
  paperPlotter.makeDoubleSearchPhasePlot(low1,low2,high,folder,extension_variations, label+", |y*| < 0.75", dataPointsOption=2, fancinessOption=1,includeTriggers = True)

  extension_variations = extension+"_withTrigLabelsAndYstar_fancinessOption2"
  paperPlotter.makeDoubleSearchPhasePlot(low1,low2,high,folder,extension_variations, label+", |y*| < 0.75", dataPointsOption=2, fancinessOption=2,includeTriggers = True)

print "Done."
