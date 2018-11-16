import glob
import ROOT

# Tools for making plots
from art.morisot import Morisot
from analysisScripts.searchphase import searchFileData
import analysisScripts.generalfunctions
import math

# Make a painter
# Initialize painter
myPainter = Morisot()
myPainter.setColourPalette("Tropical")
myPainter.setLabelType(2)
myPainter.setEPS(True)
myPainter.cutstring="|y*| < 0.75"

# Somewhere to put the plots
outDir = "plots_compareFuncChoiceUncertainties/"

# Where are all the files I need to look at?
filesToPlot = "/home/kpachal/project/kpachal/Datasets_DijetISR/search_compare_uncertainties/SearchPhase_dijetgamma_{0}_{1}.root"

spectra = [
           "compound_trigger_inclusive",
           "compound_trigger_nbtag2",
           "single_trigger_inclusive",
           "single_trigger_nbtag2"
           ]

unc_options = [
               "standard", # second best function
               "worstFunc", # worst function still passing requirements
               "changeWHW", # same function, 1 bin smaller SWIFT WHW
               "statistical" # Fit uncertainty NOT from function choice
               ]

# Luminosities in ipb
lumiDict = {"single_trigger" : 79826,
            "compound_trigger" : 76595}

# A few details for plotting
minXDict = {"single_trigger" : 169,
            "compound_trigger" : 335}

yRangeDict = {"nbtag2" : 0.15,
              "inclusive" : 0.005}

yRangeDict_withStat = {"nbtag2" : 0.2,
              "inclusive" : 0.01}

maxXForFit = 1200

for spectrum in spectra :

  plot_list = []
  name_list = []

  dataDict = {}

  for option in unc_options :

    file = filesToPlot.format(spectrum,option)

    # test file is ok
    testOpen = ROOT.TFile.Open(file,"READ")
    bit = testOpen.TestBit(ROOT.TFile.kRecovered)
    testOpen.Close()
    if bit :
      print "File",file,"not closed, skipping..."
      continue

    theseData = searchFileData(file,False)
    dataDict[option] = theseData

    lookup_channel = "_".join(spectrum.split("_")[:-1])
    lookup_selection = spectrum.split("_")[-1]

  luminosity = lumiDict[lookup_channel]
  minX = minXDict[lookup_channel]
  yrange = yRangeDict[lookup_selection]
  yrange_withStat = yRangeDict_withStat[lookup_selection]

  # Errors to plot
  # Before, there were two
  errors_standard = [dataDict["standard"].asymmFitUncertainty,dataDict["standard"].asymmFitUncertainty]
  errors_worstFunc = [dataDict["worstFunc"].asymmFitUncertainty,dataDict["worstFunc"].asymmFitUncertainty]
  errors_changeWHW = [dataDict["changeWHW"].asymmFitUncertainty,dataDict["changeWHW"].asymmFitUncertainty]

  data = dataDict["standard"].basicData
  fit = dataDict["standard"].basicBkgFromFit
  firstBin = data.FindBin(dataDict["standard"].fitLow)-1
  lastBin = data.FindBin(dataDict["standard"].fitHigh)

  # Residuals are ratios taken from the errors
  ratio_standard = fit.Clone()
  ratio_standard.SetName("ratio_standard_{0}".format(spectrum))
  ratio_worstFunc = fit.Clone()
  ratio_worstFunc.SetName("ratio_worstFunc_{0}".format(spectrum))
  ratio_changeWHW = fit.Clone()
  ratio_changeWHW.SetName("ratio_changeWHW_{0}".format(spectrum))

  ratio_standard_newDir = fit.Clone()
  ratio_standard_newDir.SetName("ratio_standard_newDir_{0}".format(spectrum))
  ratio_worstFunc_newDir = fit.Clone()
  ratio_worstFunc_newDir.SetName("ratio_worstFunc_newDir_{0}".format(spectrum))
  ratio_changeWHW_newDir = fit.Clone()
  ratio_changeWHW_newDir.SetName("ratio_changeWHW_newDir_{0}".format(spectrum))

  ratio_statistical_plus1 = fit.Clone()
  ratio_statistical_plus1.SetName("ratio_statistical_plus1_{0}".format(spectrum))
  ratio_statistical_minus1 = fit.Clone()
  ratio_statistical_minus1.SetName("ratio_statistical_minus1_{0}".format(spectrum))

  ratio_symmetric_plus1 = fit.Clone()
  ratio_symmetric_plus1.SetName("ratio_symmetric_plus1_{0}".format(spectrum))
  ratio_symmetric_minus1 = fit.Clone()
  ratio_symmetric_minus1.SetName("ratio_symmetric_minus1_{0}".format(spectrum))

  # Plus and minus stat uncertainty (equal to sqrt(n) for n taken from fit)
  stat_plus1 = fit.Clone()
  stat_plus1.SetName("statunc_plus1_{0}".format(spectrum))
  stat_minus1 = fit.Clone()
  stat_minus1.SetName("statunc_minus1_{0}".format(spectrum))

  for bin in range(0,fit.GetNbinsX()+1) :
    ratio_standard.SetBinError(bin,0)
    ratio_worstFunc.SetBinError(bin,0)
    ratio_changeWHW.SetBinError(bin,0)
    ratio_standard_newDir.SetBinError(bin,0)
    ratio_worstFunc_newDir.SetBinError(bin,0)
    ratio_changeWHW_newDir.SetBinError(bin,0)
    ratio_statistical_plus1.SetBinError(bin,0)
    ratio_statistical_minus1.SetBinError(bin,0)
    stat_plus1.SetBinError(bin,0)
    stat_minus1.SetBinError(bin,0)

    if fit.GetBinContent(bin) == 0 :
      ratio_standard.SetBinContent(bin,0)
      ratio_worstFunc.SetBinContent(bin,0)
      ratio_changeWHW.SetBinContent(bin,0)
      ratio_standard_newDir.SetBinContent(bin,0)
      ratio_worstFunc_newDir.SetBinContent(bin,0)
      ratio_changeWHW_newDir.SetBinContent(bin,0)
      ratio_statistical_plus1.SetBinContent(bin,0)
      ratio_statistical_minus1.SetBinContent(bin,0)
      stat_plus1.SetBinContent(bin,0)
      stat_minus1.SetBinContent(bin,0)
    else :
      ratio_standard.SetBinContent(bin,(dataDict["standard"].asymmFitUncertainty.GetBinContent(bin)-fit.GetBinContent(bin))/fit.GetBinContent(bin))
      ratio_worstFunc.SetBinContent(bin,(dataDict["worstFunc"].asymmFitUncertainty.GetBinContent(bin)-fit.GetBinContent(bin))/fit.GetBinContent(bin))
      ratio_changeWHW.SetBinContent(bin,(dataDict["changeWHW"].asymmFitUncertainty.GetBinContent(bin)-fit.GetBinContent(bin))/fit.GetBinContent(bin))
      ratio_standard_newDir.SetBinContent(bin,(dataDict["standard"].asymmFitUncertainty_averageDir.GetBinContent(bin)-fit.GetBinContent(bin))/fit.GetBinContent(bin))
      ratio_worstFunc_newDir.SetBinContent(bin,(dataDict["worstFunc"].asymmFitUncertainty_averageDir.GetBinContent(bin)-fit.GetBinContent(bin))/fit.GetBinContent(bin))
      ratio_changeWHW_newDir.SetBinContent(bin,(dataDict["changeWHW"].asymmFitUncertainty_averageDir.GetBinContent(bin)-fit.GetBinContent(bin))/fit.GetBinContent(bin))
      
      ratio_statistical_plus1.SetBinContent(bin,dataDict["statistical"].basicBkgFromFit.GetBinError(bin)/fit.GetBinContent(bin))
      ratio_statistical_minus1.SetBinContent(bin,-1.0*dataDict["statistical"].basicBkgFromFit.GetBinError(bin)/fit.GetBinContent(bin))
      stat_plus1.SetBinContent(bin,1.0/math.sqrt(fit.GetBinContent(bin)))
      stat_minus1.SetBinContent(bin,-1.0/math.sqrt(fit.GetBinContent(bin)))

  smallRootFile = "{0}/relevantPlots_{1}.root".format(outDir,spectrum)
  outfile = ROOT.TFile.Open(smallRootFile,"RECREATE")
  outfile.cd()
  data.Write("data")
  fit.Write("fit")
  dataDict["standard"].asymmFitUncertainty.Write("fitfuncerr_standard")
  dataDict["worstFunc"].asymmFitUncertainty.Write("fitfuncerr_worstfunc")
  dataDict["changeWHW"].asymmFitUncertainty.Write("fitfuncerr_changeWHW")
  ratio_standard.Write("ratio_standard")
  ratio_worstFunc.Write("ratio_worstfunc")
  ratio_changeWHW.Write("ratio_changeWHW")
  outfile.Close()
  
  myPainter.drawManyOverlaidHistograms([ratio_standard,ratio_worstFunc,ratio_changeWHW], ["Nominal func. choice","Worst acceptable","Smaller SWIFT WHW"], "m_{jj} [GeV]", "Events", "{0}/compareRatios_{1}".format(outDir,spectrum),firstBin,lastBin,-1.0*yrange,yrange,extraLegendLines = [],doLogX=False,doLogY=False,doErrors=False,doRectangular=False,doLegend=True,doLegendLow=True,doLegendLocation="Left",doLegendOutsidePlot=False,doATLASLabel="Low",pairNeighbouringLines=False,dotLines = [False,False,False],addHorizontalLines=[])

  myPainter.drawManyOverlaidHistograms([ratio_standard,ratio_worstFunc,ratio_changeWHW,[stat_plus1,stat_minus1]], ["Nominal func. choice","Worst acceptable","Smaller SWIFT WHW", "Data stat. unc."], "m_{jj} [GeV]", "Events", "{0}/compareRatios_withStat_{1}".format(outDir,spectrum),firstBin,lastBin,-1.0*yrange_withStat,yrange_withStat,extraLegendLines = [],doLogX=False,doLogY=False,doErrors=False,doRectangular=False,doLegend=True,doLegendLow=True,doLegendLocation="Left",doLegendOutsidePlot=False,doATLASLabel="Low",pairNeighbouringLines=False,dotLines = [False,False,False,True],addHorizontalLines=[])

  myPainter.drawManyOverlaidHistograms([ratio_standard,ratio_worstFunc,ratio_changeWHW,[ratio_statistical_plus1,ratio_statistical_minus1]], ["Nominal func. choice","Worst acceptable","Smaller SWIFT WHW", "Fit uncertainty"], "m_{jj} [GeV]", "Events", "{0}/compareRatios_withFitUnc_{1}".format(outDir,spectrum),firstBin,lastBin,-1.0*yrange_withStat,yrange_withStat,extraLegendLines = [],doLogX=False,doLogY=False,doErrors=False,doRectangular=False,doLegend=True,doLegendLow=True,doLegendLocation="Left",doLegendOutsidePlot=False,doATLASLabel="Low",pairNeighbouringLines=False,dotLines = [False,False,False,True],addHorizontalLines=[])
  
  myPainter.drawManyOverlaidHistograms([ratio_standard_newDir,ratio_worstFunc_newDir,ratio_changeWHW_newDir,[ratio_statistical_plus1,ratio_statistical_minus1]], ["Nominal func. choice","Worst acceptable","Smaller SWIFT WHW", "Fit uncertainty"], "m_{jj} [GeV]", "Events", "{0}/compareRatios_averageDirection_withFitUnc_{1}".format(outDir,spectrum),firstBin,lastBin,-1.0*yrange_withStat,yrange_withStat,extraLegendLines = [],doLogX=False,doLogY=False,doErrors=False,doRectangular=False,doLegend=True,doLegendLow=True,doLegendLocation="Left",doLegendOutsidePlot=False,doATLASLabel="Low",pairNeighbouringLines=False,dotLines = [False,False,False,True],addHorizontalLines=[])
  
  print "done"

