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

# Somewhere to put the plots
outDir = "plots_unblindedData/"
# Where are all the files I need to look at?
# In the below are the files that were originally used for INT plots after unblinding
#filesToPlot = "/cluster/warehouse/kpachal/DijetISR/Resolved2017/LimitSetting/inputs/data_outputs/SearchPhase_dijetgamma_*.root".format(tag)
# Re-calculated with new definition of function choice uncertainty
filesToPlot = "/home/kpachal/project/kpachal/Datasets_DijetISR/search_official_results/SearchPhase_dijetgamma_*.root".format(tag)

# Luminosities in ipb
lumiDict = {"dijetgamma_single_trigger" : 79826,
            "dijetgamma_compound_trigger" : 76595,
            "trijet" : 79521  }

# A few details for plotting
minXDict = {"dijetgamma_single_trigger" : 169,
            "dijetgamma_compound_trigger" : 335}

maxXForFit = 1200

channelLabels = {
"dijetgamma_single_trigger_inclusive" : "#splitline{Single-photon trigger}{Flavour inclusive}",
"dijetgamma_single_trigger_nbtag2" : "#splitline{Single-photon trigger}{2 b-tags}",
"dijetgamma_compound_trigger_inclusive" : "#splitline{Combined trigger}{Flavour inclusive}",
"dijetgamma_compound_trigger_nbtag2" : "#splitline{Combined trigger}{2 b-tags}"
}

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
  print tokens
  channel = "{0}_{1}_{2}".format(tokens[1],tokens[2],tokens[3])
  selection = tokens[4]

  theseData = searchFileData(file,False)
  if not theseData :
    continue
  luminosity = lumiDict[channel]

  label = channelLabels[channel+"_"+selection]
  myPainter.cutstring=label

  extension = "_{0}_{1}".format(channel,selection)

  # Options: painter object, low x, high x, luminosity, folder, plot extension, extra legend lines
  theseData.makeSearchPhasePlots(myPainter,minXDict[channel],maxXForFit,luminosity,outDir,extension)

