import ROOT
from analysisScripts.searchphase import searchFileData
import analysisScripts.generalfunctions

filePath = "/home/kpachal/project/kpachal/Datasets_DijetISR/search_official_results/"
low1 = 169
low2 = 335
high = 1200

channelLabels = {
"inclusive" : "Flavour inclusive",
"nbtag2" : "2 b-tags"
}
triggerLabels = {
"single_trigger" : "Single-photon trigger",
"compound_trigger" : "Combined trigger"
}

print "========================================="
print "Binning:"

# example hist: just get single inclusive and grab the data hist.
testFile = ROOT.TFile.Open(filePath+"SearchPhase_dijetgamma_single_trigger_inclusive.root","READ")
binningHist = testFile.Get("basicData")
binningHist.SetDirectory(0)
testFile.Close()

bins = []
bins_combined = []
for bin in range (1,binningHist.GetNbinsX()+1) :
  if binningHist.GetBinLowEdge(bin) < low1 or binningHist.GetBinLowEdge(bin)+binningHist.GetBinWidth(bin) > high :
    continue
  print bin, "\t", binningHist.GetBinLowEdge(bin), binningHist.GetBinLowEdge(bin)+binningHist.GetBinWidth(bin)
  bins.append(bin)
  if binningHist.GetBinLowEdge(bin) < low2 :
    continue
  bins_combined.append(bin)

for channel in ["inclusive","nbtag2"] :

  for trigger in ["single_trigger", "compound_trigger"] :

    print "========================================="
    print channelLabels[channel],",",triggerLabels[trigger]

    fileTemplate = filePath + "SearchPhase_dijetgamma_{0}_{1}.root".format(trigger,channel)
    # Get search file data
    data = searchFileData(fileTemplate)

    print "==------------------------------------=="
    print "Data"

    for bin in bins :
      if "compound" in trigger and bin not in bins_combined :
        print bin, "\t 0.0"
        continue

      print bin,"\t",data.basicData.GetBinContent(bin)

    print "==------------------------------------=="
    print "Background prediction"

    for bin in bins :
      if "compound" in trigger and bin not in bins_combined :
        print bin, "\t 0.0"
        continue

      print bin, "\t",data.basicBkgFromFit.GetBinContent(bin)

    print "==------------------------------------=="
    print "Fit uncertainty (stat)"

    for bin in bins :
      if "compound" in trigger and bin not in bins_combined :
        print bin, "\t 0.0"
        continue

      print bin, "\t",data.nominalPlus1Stat.GetBinContent(bin)-data.basicBkgFromFit.GetBinContent(bin)

    print "==------------------------------------=="
    print "Fit uncertainty (function choice)"

    for bin in bins :
      if "compound" in trigger and bin not in bins_combined :
        print bin, "\t 0.0"
        continue

      print bin, "\t",data.asymmFitUncertainty.GetBinContent(bin)-data.basicBkgFromFit.GetBinContent(bin)

    print "\n"

