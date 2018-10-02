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

CME = 13000

doIndividualPlots = False

class multiPlotter :

  def __init__(self,lower_mass_data,higher_mass_data,permitWindow=False) :

    self.data1 = lower_mass_data
    self.data2 = higher_mass_data

  def makeDoubleSearchPhasePlot(self,low1,low2,high,folder,ext,funcName=[],lumi1=0,lumi2=0,sigfile1="",sigfile2="",scale1=1,scale2=2) :

    firstBin = self.basicData.FindBin(low1)
    lastBin = self.basicData.FindBin(high)
    
    # Fancy figure 1 (two spectra)
    myPainter.drawDataAndFitsOverSignificanceHists_TwoSpectra(self.basicData,self.basicBkgFrom4ParamFit,self.residualHist,self.basicData2,self.basicBkgFrom4ParamFit2,self.residualHist2,sig1,sig2,scale1,scale2,\
         'm_{jj} [GeV]','Events','Significance','{0}/figure1'.format(folder)+ext,\
         lumi1,lumi2,"","",13,low1,high,firstBin,lastBin,True,self.bumpLowEdge,self.bumpHighEdge,self.bumpLowEdge2,self.bumpHighEdge2,extraLegendLines=funcName,pval = self.bumpHunterPVal,chi2pval = self.chi2PVal,pval2 = self.bumpHunterPVal2,chi2pval2 = self.chi2PVal2,doRectangular=False)

    # Two separate plots for showing individual uncertainties
    # Materials are always [single, compound] in this order.
    for trigger, info in ["single","compound"] :

      thisData =

        myPainter.drawDataWithFitAsHistogramAndResidualPaper(newbasicdata,newbasicBkgFrom4ParamFit,luminosity,13,"m_{jj} [TeV]","Events",["Data","Fit","Statistical uncertainty on fit","Function choice"],"{0}/compareFitQualityAndFitChoice_Asymm_WithRatioPaper".format(folderextension),drawError = True, [[newNomPlus1,newNomMinus1],[placeHolderNom,newValueNewFuncErrDirected]],[altFitRatio,MinusNomRatio,PlusNomRatio],firstBin,lastBin-1,True,True,True,False,False) # changed from lastBin+2 to lastBin+15E3 to match FancyFigure

# Pass lumi in pb
def ProcessSignalFile(sigfile,lumi,scale,binLow,binHigh) :

  openSigFile = ROOT.TFile.Open(sigfile,"READ")
  sig = openSigFile.Get("m_jj_Nominal")
  sig.SetDirectory(0)
  openSigFile.Close()
  sig1.Scale(lumi)
    
  for ibin in xrange(1+sig1.GetNbinsX()):
    if ibin < binLow or ibin > binHigh:
      sig1.SetBinContent(ibin,0)
    if self.basicData.GetBinContent(ibin) > 0:
      # This was 0.0022 for second, initially.
      if sig.GetBinContent(ibin)/self.basicData.GetBinContent(ibin) < 0.0001:
        sig.SetBinContent(ibin,0)
    if sig.GetBinContent(ibin) > 0:
      sig.SetBinContent(ibin,scale*sig1.GetBinContent(ibin)+self.basicData.GetBinContent(ibin))

  return sig

#####################################################
# User controlled
folder = "plots/"
filePath = "/cluster/warehouse/kpachal/DijetISR/Resolved2017/LimitSetting/BayesianFramework/results/search_official_results/"
signalFileTemplate = "/cluster/warehouse/kpachal/DijetISR/Resolved2017/LimitSetting/inputs/signal_hists/Zprime_{0}_{1}_MGPy8EG_N30LO_A14N23LO_dmA_jja_Ph100_mRp{2}_mD10_gSp{3}_gD1.root"

#make plots folder i.e. make folder extension
if not os.path.exists(folder):
  os.makedirs(folder)

lumiDict = { "single_trigger" : 79826,
             "compound_trigger" : 76595
            }

signalProcessingInfo = {
  "single_trigger" : {
    "inclusive" : {"mass" : "25",
                   "scale" : 10,
                   "binLow" : 3,
                   "binHigh" : 15},
    "2tag"      : {"mass" : "25",
                   "scale" : 1,
                   "binLow" : 3,
                   "binHigh" : 15}
  },
  "compound_trigger" : {
    "inclusive" : { "mass" : "55",
                    "scale" : 10,
                    "binLow" : 6,
                    "binHigh" : 25},
    "2tag"      : { "mass" : "55",
                    "scale" : 1,
                    "binLow" : 6,
                    "binHigh" : 25}
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
    
    # Process signal file and append its info
    signalFile = signalFileTemplate.format(trigger,channel,signalMasses[trigger],"1")
    sigHist = ProcessSignalFile(signalFile)
    data.sigHist = sigHist
    
    # Append luminosity info
    lumi = lumiDict[trigger]
    data.luminosity = lumi
  
    channel_data[trigger] = data

  extension = "_"+channel

  theseData = paperPlotter(fileTemplate1,fileTemplate2,True)
  theseData.makeSearchPhasePlots(low1,low2,high,folder,extension,lumi1=lumi1,lumi2=lumi2,sigfile1=signalFile1,sigfile2=signalFile2,scale1=scale,scale2=scale)

print "Done."
