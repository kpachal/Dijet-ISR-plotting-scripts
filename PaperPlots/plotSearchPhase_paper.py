#!/usr/bin/env python

import os
import ROOT
from art.morisot import Morisot
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

class searchFileData :

  def __init__(self,filename,filename2,permitWindow=False) :

    self.permitWindow = permitWindow

    searchInputFile = ROOT.TFile.Open(filename,"READ")
    searchInputFile2 = ROOT.TFile.Open(filename2,"READ")

    # Retrieve search phase inputs
    self.basicData = searchInputFile.Get("basicData")
    self.basicData2 = searchInputFile2.Get("basicData")
    self.basicData.SetDirectory(0)
    self.basicData2.SetDirectory(0)
    self.basicBkgFrom4ParamFit = searchInputFile.Get("basicBkgFrom4ParamFit")
    self.basicBkgFrom4ParamFit.SetDirectory(0)
    self.basicBkgFrom4ParamFit2 = searchInputFile2.Get("basicBkgFrom4ParamFit")
    self.basicBkgFrom4ParamFit2.SetDirectory(0)
    self.residualHist = searchInputFile.Get("residualHist")
    self.residualHist.SetDirectory(0)
    self.residualHist2 = searchInputFile2.Get("residualHist")
    self.residualHist2.SetDirectory(0)
    self.relativeDiffHist = searchInputFile.Get("relativeDiffHist")
    self.relativeDiffHist.SetDirectory(0)
    self.sigOfDiffHist = searchInputFile.Get("sigOfDiffHist")
    self.sigOfDiffHist.SetDirectory(0)
    self.logLikelihoodPseudoStatHist = searchInputFile.Get("logLikelihoodStatHistNullCase")
    self.logLikelihoodPseudoStatHist.SetDirectory(0)
    self.chi2PseudoStatHist = searchInputFile.Get("chi2StatHistNullCase")
    self.chi2PseudoStatHist.SetDirectory(0)
    self.bumpHunterStatHist = searchInputFile.Get("bumpHunterStatHistNullCase")
    self.bumpHunterStatHist.SetDirectory(0)
    self.bumpHunterTomographyPlot = searchInputFile.Get('bumpHunterTomographyFromPseudoexperiments')

    bumpHunterStatOfFitToData = searchInputFile.Get('bumpHunterStatOfFitToData')
    bumpHunterStatOfFitToData2 = searchInputFile2.Get('bumpHunterStatOfFitToData')
    logLOfFitToDataVec = searchInputFile.Get('logLOfFitToData')
    chi2OfFitToDataVec = searchInputFile.Get('chi2OfFitToData')
    statOfFitToData = searchInputFile.Get('bumpHunterPLowHigh')
    chi2OfFitToDataVec2 = searchInputFile2.Get('chi2OfFitToData')    
    statOfFitToData2 = searchInputFile2.Get('bumpHunterPLowHigh')
    self.logLOfFitToData = logLOfFitToDataVec[0]
    self.logLPVal = logLOfFitToDataVec[1]
    self.chi2OfFitToData = chi2OfFitToDataVec[0]
    self.chi2PVal = chi2OfFitToDataVec[1]
    self.bumpHunterStatFitToData = statOfFitToData[0]
    self.bumpHunterPVal = bumpHunterStatOfFitToData[1]
    self.bumpLowEdge = statOfFitToData[1]
    self.bumpHighEdge = statOfFitToData[2]
    self.chi2OfFitToData2 = chi2OfFitToDataVec2[0]
    self.chi2PVal2 = chi2OfFitToDataVec2[1]
    self.bumpHunterStatFitToData2 = statOfFitToData2[0]
    self.bumpHunterPVal2 = bumpHunterStatOfFitToData2[1]
    self.bumpLowEdge2 = statOfFitToData2[1]
    self.bumpHighEdge2 = statOfFitToData2[2]
   
    print "Chi2 Value: ", float(self.chi2OfFitToData)
    print "Chi2 p Value: ", self.chi2PVal
    print "BH Value: ", float(self.bumpHunterStatFitToData)
    print "BH p Value: ", self.bumpHunterPVal
    print "BH significance: ", GetZVal(self.bumpHunterPVal,False)

    self.NDF = searchInputFile.Get('NDF')[0]

    excludeWindowNums = searchInputFile.Get('excludeWindowNums')
    self.excludeWindow = int(excludeWindowNums[0]+0.5)
    self.bottomWindowEdge = excludeWindowNums[1]
    self.topWindowEdge = excludeWindowNums[2]

    if (self.excludeWindow and self.permitWindow) :
      statsOfRemainingSpectrum = searchInputFile.Get("BHLogLAndChi2OfRemainderAfterWindow")
      self.BHPValRemainder = statsOfRemainingSpectrum[0]
      self.LogLPValRemainder = statsOfRemainingSpectrum[1]
      self.Chi2PValRemainder = self.calculateRemainingChi2() #statsOfRemainingSpectrum[2]
     

    searchInputFile.Close()

  def getPValErrs(self) :

    # (DeltaX/X)^2 = (1/DeltaX)^2 = 1/X: set errors
    nRightBH = self.bumpHunterStatHist.Integral(self.bumpHunterStatHist.FindBin(self.bumpHunterStatFitToData),self.bumpHunterStatHist.GetNbinsX())
    nLeftBH = self.bumpHunterStatHist.Integral() - nRightBH
    if nRightBH > 0 and nLeftBH > 0 : deltaPvalBH = self.bumpHunterPVal * math.sqrt(1/nRightBH + 1/nLeftBH)
    else : deltaPvalBH = 0
    nRightChi2 = self.chi2PseudoStatHist.Integral(self.chi2PseudoStatHist.FindBin(self.chi2OfFitToData),self.chi2PseudoStatHist.GetNbinsX())
    nLeftChi2 = self.chi2PseudoStatHist.Integral() - nRightChi2
    if nRightChi2 > 0 and nLeftChi2 > 0 : deltaPvalChi2 = self.chi2PVal * math.sqrt(1/nRightChi2 + 1/nLeftChi2)
    else : deltaPvalChi2 = 0
    nRightLogL = self.logLikelihoodPseudoStatHist.Integral(self.logLikelihoodPseudoStatHist.FindBin(self.logLOfFitToData),self.logLikelihoodPseudoStatHist.GetNbinsX())
    nLeftLogL = self.logLikelihoodPseudoStatHist.Integral() - nRightLogL
    if nRightLogL > 0 and nLeftLogL > 0 : deltaPvalLogL = self.logLPVal * math.sqrt(1/nRightLogL + 1/nLeftLogL)
    else : deltaPvalLogL = 0

    return deltaPvalBH,deltaPvalChi2,deltaPvalLogL

  def calculateRemainingChi2(self) :

    firstBin = 0
    for bin in range(1,self.basicBkgFrom4ParamFit.GetNbinsX()+2) :
      firstBin = bin
      if self.basicBkgFrom4ParamFit.GetBinContent(bin) > 0 :
        break
    lastBin = 0
    for bin in range(self.basicBkgFrom4ParamFit.GetNbinsX()+1,0,-1) :
      lastBin = bin
      if self.basicBkgFrom4ParamFit.GetBinContent(bin) > 0 :
        break
    firstWindowBin = 0
    lastWindowBin = 0
    if self.excludeWindow :
      for bin in range(1,self.basicBkgFrom4ParamFit.GetNbinsX()+2) :
        if math.fabs(self.basicBkgFrom4ParamFit.GetBinLowEdge(bin) - self.bottomWindowEdge) < 0.1 :
          firstWindowBin = bin
        if math.fabs(self.basicBkgFrom4ParamFit.GetBinLowEdge(bin)+self.basicBkgFrom4ParamFit.GetBinWidth(bin) - self.topWindowEdge) < 0.1 :
          lastWindowBin = bin

    answer = 0
    for bin in range(firstBin,lastBin+1) :

      if self.excludeWindow and bin >= firstWindowBin and bin <= lastWindowBin : continue

      d = self.basicData.GetBinContent(bin)
      if (d==0) : continue
      b = self.basicBkgFrom4ParamFit.GetBinContent(bin)
      deltaB = self.basicBkgFrom4ParamFit.GetBinError(bin)

      term = (d - b) / math.sqrt(b+deltaB*deltaB)
      answer = answer + (term*term)

    nRightChi2 = self.chi2PseudoStatHist.Integral(self.chi2PseudoStatHist.FindBin(answer),self.chi2PseudoStatHist.GetNbinsX())
    nTotal = self.chi2PseudoStatHist.Integral()
    return float(nRightChi2)/float(nTotal)


  def makeSearchPhasePlots(self,low1,low2,high,folder,ext,funcName=[],lumi1=0,lumi2=0,sigfile1="",sigfile2="",scale1=1,scale2=2) :

    openSigFile1 = ROOT.TFile.Open(sigfile1,"READ")
    sig1 = openSigFile1.Get("m_jj_Nominal")
    sig1.SetDirectory(0)
    openSigFile1.Close()
    sig1.Scale(lumi1/1000.0)
    print "Signal 1:"
    sig1.Print("all")
    
    for ibin in xrange(1+sig1.GetNbinsX()):
      if ibin < 3 or ibin > 15:
        sig1.SetBinContent(ibin,0)
      if self.basicData.GetBinContent(ibin) > 0:
        if sig1.GetBinContent(ibin)/self.basicData.GetBinContent(ibin) < 0.0001:
          sig1.SetBinContent(ibin,0)
      if sig1.GetBinContent(ibin) > 0:
        sig1.SetBinContent(ibin,scale1*sig1.GetBinContent(ibin)+self.basicData.GetBinContent(ibin))

    openSigFile2 = ROOT.TFile.Open(sigfile2,"READ")
    sig2 = openSigFile2.Get("m_jj_Nominal")
    sig2.SetDirectory(0)
    openSigFile2.Close()
    sig2.Scale(lumi2/1000.0)
    for ibin in xrange(1+sig2.GetNbinsX()):
      if ibin < 6 or ibin > 25:
        sig2.SetBinContent(ibin,0)
      if self.basicData.GetBinContent(ibin) > 0:
        if sig2.GetBinContent(ibin)/self.basicData.GetBinContent(ibin) < 0.0022:
          sig2.SetBinContent(ibin,0)
      if sig2.GetBinContent(ibin) > 0:
        sig2.SetBinContent(ibin,scale2*sig2.GetBinContent(ibin)+self.basicData2.GetBinContent(ibin))


    firstBin = self.basicData.FindBin(low1)
    lastBin = self.basicData.FindBin(high)
    
    print "start"
    myPainter.drawDataAndFitsOverSignificanceHists_TwoSpectra(self.basicData,self.basicBkgFrom4ParamFit,self.residualHist,self.basicData2,self.basicBkgFrom4ParamFit2,self.residualHist2,sig1,sig2,scale1,scale2,\
         'm_{jj} [GeV]','Events','Significance','{0}/figure1'.format(folder)+ext,\
         lumi1,lumi2,"","",13,low1,high,firstBin,lastBin,True,self.bumpLowEdge,self.bumpHighEdge,self.bumpLowEdge2,self.bumpHighEdge2,extraLegendLines=funcName,pval = self.bumpHunterPVal,chi2pval = self.chi2PVal,pval2 = self.bumpHunterPVal2,chi2pval2 = self.chi2PVal2,doRectangular=False)
    print "end"

def GetZVal (p, excess) :
  #the function normal_quantile converts a p-value into a significance,
  #i.e. the number of standard deviations corresponding to the right-tail of 
  #a Gaussian
  if excess :
    zval = ROOT.Math.normal_quantile(1-p,1);
  else :
    zval = ROOT.Math.normal_quantile(p,1);

  return zval


def MakeHistoFromStats(statistics) :

  nentries = len(statistics)
  nBins = int(float(nentries)/10.0)

  maxVal = max(statistics)
  minVal = min(statistics)
  axisrange = maxVal - minVal;

  thismin = minVal-0.05*axisrange;
  thismax = maxVal+0.05*axisrange;

  statPlot = ROOT.TH1D("statPlot","",nBins,thismin,thismax)
  for val in range(len(statistics)) :
    statPlot.Fill(statistics[val])

  return statPlot

# User controlled
folder = "plots/"
filePath = "/cluster/warehouse/kpachal/DijetISR/Resolved2017/LimitSetting/BayesianFramework/results/search_official_results/"
signalFileTemplate = "/cluster/warehouse/kpachal/DijetISR/Resolved2017/LimitSetting/inputs/signal_hists/Zprime_{0}_{1}_MGPy8EG_N30LO_A14N23LO_dmA_jja_Ph100_mRp{2}_mD10_gSp{3}_gD1.root"

#make plots folder i.e. make folder extension
if not os.path.exists(folder):
  os.makedirs(folder)

lumiList = [
            79826, # single trigger
            76595 # compound trigger
            ]

for channel in ["inclusive","nbtag2"] :

  fileTemplate1 = filePath+"SearchPhase_dijetgamma_single_trigger_{0}.root".format(channel)
  fileTemplate2 = filePath+"SearchPhase_dijetgamma_compound_trigger_{0}.root".format(channel)

  signalFile1 = signalFileTemplate.format("single_trigger",channel,"25","1")
  signalFile2 = signalFileTemplate.format("compound_trigger_2016",channel,"55","1")

  pValCutoff = 0.05
  low1 = 169
  low2 = 335
  high = 1200
  
  lumi1 = lumiList[0]
  lumi2 = lumiList[1]

  extension = "_"+channel

  if "inclusive" in channel :
    scale = 10
  else :
    scale = 1

  theseData = searchFileData(fileTemplate1,fileTemplate2,True)
  theseData.makeSearchPhasePlots(low1,low2,high,folder,extension,lumi1=lumi1,lumi2=lumi2,sigfile1=signalFile1,sigfile2=signalFile2,scale1=scale,scale2=scale)

print "Done."
