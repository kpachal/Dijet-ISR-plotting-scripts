import glob
import ROOT

# Tools for making plots
from art.morisot import Morisot
from analysisScripts.searchphase import searchFileData
import analysisScripts.generalfunctions

# Tag: if you only want to make things matching a pattern
tag = ""

# Make a painter
# Initialize painter
myPainter = Morisot()
myPainter.setColourPalette("Teals")
myPainter.setLabelType(2)
myPainter.setEPS(True)
myPainter.cutstring="|y*| < 0.75"

# Somewhere to put the plots
outDir = "plots_unblindedData/"
# Where are all the files I need to look at?
filesToPlot = "/cluster/warehouse/kpachal/DijetISR/Resolved2017/LimitSetting/inputs/data_outputs/SearchPhase_dijetgamma_*.root".format(tag)
# Luminosities in ipb
lumiDict = {"dijetgamma_single_trigger" : 79826,
            "dijetgamma_compound_trigger" : 76595,
            "trijet" : 79521  }

# A few details for plotting
minXDict = {"dijetgamma_single_trigger" : 169,
            "dijetgamma_compound_trigger" : 335}

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
  channel = "{0}_{1}_{2}".format(tokens[1],tokens[2],tokens[3])
  extratag = tokens[4]

  theseData = searchFileData(file,False)
  if not theseData :
    continue
  luminosity = lumiDict[channel]

  # Options: painter object, low x, high x, luminosity, folder, plot extension, extra legend lines
  extension = "_{0}_{1}".format(channel,extratag)
  theseData.makeSearchPhasePlots(myPainter,minXDict[channel],maxXForFit,luminosity,outDir,extension)

