import glob

# Tools for making plots
from art.morisot import Morisot
from analysisScripts.searchphase import searchFileData
import analysisScripts.generalfunctions

# Settings: plot MC OR data
plotMC = False

# Make a painter
# Initialize painter
myPainter = Morisot()
myPainter.setColourPalette("Teals")
myPainter.setLabelType(2)
myPainter.setEPS(True)

# Somewhere to put the tables
if plotMC :
  outDir = "tables_fitsToMC16aMC16d/"
  # Where are all the files I need to look at?
  filesToPlot = "/cluster/warehouse/kpachal/DijetISR/Resolved2017/LimitSetting/BayesianFramework/results/testFitsMC16aMC16d/SearchPhase_*.root"
  # Luminosities in ipb
  lumiDict = {"dijetgamma" : 79826,
            "trijet" : 79521  }
else :
  outDir = "tables_fitsTo15fbData/"
  filesToPlot = "/cluster/warehouse/kpachal/DijetISR/Resolved2017/LimitSetting/BayesianFramework/results/testFits15ifb/SearchPhase_*.root"
  # Luminosities in ipb
  lumiDict = {"dijetgamma" : 15000,
            "trijet" : 15000  }


# A few details for plotting
minXDict = {"dijetgamma_single_trigger" : 169,
            "dijetgamma_compound_trigger" : 300,
            "trijet" : 300}



maxXForFit = 1200

separateFiles = {"channel" : ["dijetgamma_single_trigger","dijetgamma_compound_trigger","trijet"], "ntag" : ["inclusive","nbtag2"], "function" : ["fourpar","fivepar","UA2"]}

# Function to sort files in the order I want
def f(item):
  if "global" in item : return 0, item
  whwtoken = item.split("_")[-1].replace("whw","").replace(".root","")
  if len(whwtoken) < 2 : return 1, item
  else : return 2, item

for channel in separateFiles["channel"] :
  for ntag in separateFiles["ntag"] :
  
    # Make a text file for this table. Format:
    # function & fit detail & BH pval & chi2 pval
    outfilename = outDir + "tab_{0}_{1}.txt".format(channel,ntag)

    for function in separateFiles["function"] :

      fileglob = glob.glob(filesToPlot.replace("*","{0}_{1}_{2}_*".format(channel,ntag,function)))
      filesToUse = sorted(fileglob,key=f)
      
      print "Looking for files with names like",filesToPlot.replace("*","{0}_{1}_{2}_*")
      print "Got:",filesToUse

      valsSoFar = {}
      for filename in filesToUse :
        print filename

        token = filename.split("_")[-1].replace("whw","").replace(".root","")
        print token

        if "whw" in token :
          whw_val = eval(token)
          extraLegendLines = ["SWIFT WHW {0}".format(whw_val)]
          extension = whw
        else :
          extraLegendLines = ["Global fit"]
          extension = "global"
        theseData = searchFileData(filename,False)
        if not theseData :
          continue

        # Get chi2 pval
        chi2 = theseData.chi2PVal
        BH = theseData.bumpHunterPVal

        print extension, chi2, BH

#        theseData.makeSearchPhasePlots(myPainter,minXDict[channel],maxXForFit,luminosity,outDir,extension,extraLegendLines)

