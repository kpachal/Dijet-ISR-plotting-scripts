import glob
import ROOT

# Tools for making plots
from art.morisot import Morisot
from analysisScripts.searchphase import searchFileData
import analysisScripts.generalfunctions

# Make a painter
# Initialize painter
myPainter = Morisot()
myPainter.setColourPalette("Teals")
myPainter.setLabelType(2)
myPainter.setEPS(True)
myPainter.cutstring="|y*| < 0.75"

# Somewhere to put the plots
outDir = "plots_fitUncIn15ifb/"
# Where are all the files I need to look at?
filesToPlot = "/home/kpachal/project/kpachal/Datasets_DijetISR/search_results_15if_systexamples/SearchResult_*.root"
print "looking for:",filesToPlot

# Luminosities in ipb
lumiDict = {"single_trigger_inclusive" : 15000,
            "compound_trigger_inclusive" : 15000}

# A few details for plotting
minXDict = {"single_trigger_inclusive" : 169,
            "compound_trigger_inclusive" : 335}

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
  channel = "{0}_{1}_{2}".format(tokens[-3],tokens[-2],tokens[-1])

  theseData = searchFileData(file,False)
  if not theseData :
    continue
  luminosity = lumiDict[channel]

  # Options: painter object, low x, high x, luminosity, folder, plot extension, extra legend lines
  theseData.makeSearchPhasePlots(myPainter,minXDict[channel],maxXForFit,luminosity,outDir,channel)

