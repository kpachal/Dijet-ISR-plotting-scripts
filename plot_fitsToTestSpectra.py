import glob
import ROOT

# Tools for making plots
from art.morisot import Morisot
from analysisScripts.searchphase import searchFileData
import analysisScripts.generalfunctions

# Settings: plot MC OR data
plotMC = True
# Tag: if you only want to make things matching a pattern
tag = ""#"_MODIFIED"

# Make a painter
# Initialize painter
myPainter = Morisot()
myPainter.setColourPalette("Teals")
myPainter.setLabelType(2)
myPainter.setEPS(True)
myPainter.cutstring="|y*| < 0.75"

# Somewhere to put the plots
if plotMC :
  outDir = "plots_fitsToMC16aMC16d/"
  # Where are all the files I need to look at?
  filesToPlot = "/cluster/warehouse/kpachal/DijetISR/Resolved2017/LimitSetting/BayesianFramework/results/testFitsMC16aMC16d/SearchPhase_*{0}*.root".format(tag)
  # Luminosities in ipb
  lumiDict = {"dijetgamma_single_trigger" : 79826,
            "dijetgamma_compound_trigger" : 76595,
            "trijet" : 79521  }
else :
  outDir = "plots_fitsTo15fbData/"
  filesToPlot = "/cluster/warehouse/kpachal/DijetISR/Resolved2017/LimitSetting/BayesianFramework/results/testFits15ifb/SearchPhase_*{0}*.root".format(tag)
  # Luminosities in ipb
  lumiDict = {"dijetgamma" : 15000,
            "trijet" : 15000  }


# A few details for plotting
minXDict = {"dijetgamma_single_trigger" : 169,
            "dijetgamma_compound_trigger" : 300,
            "trijet" : 300}



maxXForFit = 1200

for file in glob.glob(filesToPlot) :

  print "Beginning file",file

  # test file is ok
  testOpen = ROOT.TFile.Open(file,"READ")
  bit = testOpen.TestBit(ROOT.TFile.kRecovered) 
  testOpen.Close()
  if bit :
    print "File",file,"not closed, skipping..."
    continue

  filename = file.split("/")[-1]
  tokens = filename.replace(".root","").split("_")
  lookup_channel = tokens[1]
  if "trijet" in lookup_channel :
    channel = lookup_channel
    nextval = 2
  else :
    channel = "{0}_{1}_{2}".format(lookup_channel,tokens[2],tokens[3])
    nextval = 4
  tags = tokens[nextval]
  function = tokens[nextval+1]

  if "whw" in filename :
    whw = tokens[nextval+2]
    whw_val = whw.replace("whw","")
    extraLegendLines = ["SWIFT WHW {0}".format(whw_val)]
    extension = whw
  else :
    extraLegendLines = ["Global fit"]
    extension = "global"
  plotNameTag = ""
  if len(tokens) > nextval+3 :
    plotNameTag = "_"+tokens[nextval+3]
  theseData = searchFileData(file,False)
  if not theseData :
    continue
  luminosity = lumiDict[lookup_channel]
  #extraLegendLines.append("BH p-val. = {0}".format(theseData.bumpHunterPVal))
  # Options: painter object, low x, high x, luminosity, folder, plot extension, extra legend lines
  extension = "_{0}_{1}_{2}_{3}{4}".format(channel,tags,function,extension,plotNameTag)
  theseData.makeSearchPhasePlots(myPainter,minXDict[channel],maxXForFit,luminosity,outDir,extension,extraLegendLines)

