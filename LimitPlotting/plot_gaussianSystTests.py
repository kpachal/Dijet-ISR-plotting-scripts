import ROOT

import os
import sys
import ROOT
from art.morisot import Morisot
from itertools import repeat

# Initialize painter
myPainter = Morisot()
# Internal
myPainter.setColourPalette("notSynthwave")
myPainter.setLabelType(2)
myPainter.dodrawUsersText = True

luminosity_dict = {
"2bComp" : 76.595*1000,
"2bSing" : 79.826*1000,
}

names = { '0.00':  '#sigma_{G}/m_{G} = Res.',
          '0.05': '#sigma_{G}/m_{G} = 0.05',
          '0.07': '#sigma_{G}/m_{G} = 0.07',
          '0.10': '#sigma_{G}/m_{G} = 0.10',
          '0.15': '#sigma_{G}/m_{G} = 0.15'}

xAxisLimits = [180,1070]

incYlims = [1E-3,5]
twobYlims = [1E-4,1E-1]

# For checking what's in a file
def GetKeyNames(self,dir=""):
  self.cd(dir)
  return [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]
ROOT.TFile.GetKeyNames = GetKeyNames


# Make sure you have already checked compatibility of x points
# in the graphs before you do this
# Does graph 1 divided by graph 2
def ratio_graphs(graph1, graph2) :

  ratioGraph = ROOT.TGraphErrors()
  for point in range(graph1.GetN()) :

    g1_x = ROOT.Double(0)
    g1_y = ROOT.Double(0)
    graph1.GetPoint(point,g1_x,g1_y)

    foundMatchingPoint = False
    for point2 in range(graph2.GetN()) :
      g2_x = ROOT.Double(0)
      g2_y = ROOT.Double(0)
      graph2.GetPoint(point2,g2_x,g2_y)
      if g2_x == g1_x :
        foundMatchingPoint = True
        break
  
    if not foundMatchingPoint : continue

    ratio = g1_y/g2_y
    
    #ratioerr_x = graph1.GetErrorX(point)

    #erry_g1 = graph1.GetErrorY(point)
    #erry_g2 = graph2.GetErrorY(point)
    #ratioerr_y = math.sqrt(math.pow(erry_g1/g1_y,2) + math.pow(erry_g2/g2_y,2))

    ratioGraph.SetPoint(ratioGraph.GetN(),g1_x,ratio)
    #ratioGraph.SetPointError(ratioGraph.GetN(),ratioerr_x,ratioerr_y)

  return ratioGraph

def makeBand(graph1, graph2):
  points = []
  for i in range(graph1.GetN()):
    points += [(i,graph1.GetX()[i],graph1.GetY()[i])]
  for i in range(graph2.GetN()-1,-1,-1):
    points += [(i,graph2.GetX()[i],graph2.GetY()[i])]
  graph_band = ROOT.TGraph();
  for i in range (len(points)): graph_band.SetPoint(i,points[i][1],points[i][2])
  return graph_band

def GetCenterAndSigmaDeviations(inputs) :
  inputs = sorted(inputs)
  statVals = []
  quantiles = [0.02275,0.1587,0.5,0.8413,0.9772]
  for q in quantiles:
    wantEvents = len(inputs)*q
    statVals.append(inputs[int(wantEvents)])
  return statVals

def getGaussianLimits( inputfileform, ratios, dataset, luminosity, cutstring, makePlot = True, outfolder = "./plots/",doEntireRange= True):

  if not os.path.isdir(outfolder):
    os.makedirs(outfolder)

  outname = outfolder+"GenericGaussians_"+dataset

  myPainter.cutstring = cutstring

  if doEntireRange:  
    outname += "_entireRange"

  basicInputFiles = {}
  for r in ratios: basicInputFiles[r] = []

  import glob
  file_list = glob.glob(inputfileform)
  
  if file_list == []:
    print "WARNING: empty file list for form", inputfileform

  for f in file_list:
    ratio_str = f.split('.')[-2].split('_')[-1]
    if ratio_str == 'resolutionwidth': ratio = 0.0
    else: ratio = float(ratio_str)/1000.
    if ratio in basicInputFiles: basicInputFiles[ratio].append(f)

  minMassVal = {}
  values = {}
  mass_list = []
  allobserved = []
  allexpected = []
  allexpected1Sigma = []
  allexpected2Sigma = []
  
  results = {}
  
  for r in ratios:
    values[r] = {}
    for f in basicInputFiles[r]:

      file = ROOT.TFile.Open(f)
      print "opened file",f
      file_keys = file.GetKeyNames()
      
      CLs_str = "CLsPerMass_widthToMass%d"%(r*1000)
      if not file or not file.Get(CLs_str): continue
      CLs = file.Get(CLs_str) 
      masses = file.Get("massesUsed")

      if masses == None: continue

      for i,mass in enumerate(masses) :

          #if doEntireRange is set, we want to plot each trigger channel over the entire mass range. 
          #If not, e.g. when making the combined plot, we split at the place where they are to be joined:
          if not doEntireRange and "Sing" in dataset and mass > channelSplit: continue

          if not doEntireRange and "Comp" in dataset and mass < channelSplit: continue
          #we never want to plot any limits above 1200 anyway
          if mass > 1200: continue
         
          mass_list += [mass]
          
          doExp = False
          treename = "ensemble_tree_%d_%d"%(mass,r*1000)
          if treename in file_keys :
            print "Getting",treename
            doExp = True
            PE_tree = file.Get(treename)
            PE_CLs = []
            for event in PE_tree:
                PE_CLs.append( event.GetBranch("95quantile_marginalized_0").GetListOfLeaves().At(0).GetValue() )
            expCLs = GetCenterAndSigmaDeviations(PE_CLs)

            #print mass, CLs[i]/luminosity, [e/luminosity for e in expCLs]
          if mass not in values[r]:
            values[r][mass] = {'obs': [], 'exp': [], 'PEs': [] }
          values[r][mass]['obs'].append(CLs[i]/luminosity)
          if doExp :
            values[r][mass]['exp'].append(expCLs[2]/luminosity)
            values[r][mass]['PEs'] += [e/luminosity for e in PE_CLs]
          else :
            values[r][mass]['exp'].append(None)
            values[r][mass]['PEs'] += []

    mass_list = sorted(list(set(mass_list)))
    thisobserved = ROOT.TGraph()
    thisexpected = ROOT.TGraph()
    thisexpected_plus1  = ROOT.TGraph()
    thisexpected_minus1 = ROOT.TGraph()
    thisexpected_plus2  = ROOT.TGraph()
    thisexpected_minus2 = ROOT.TGraph()
    for m in mass_list :
      if m not in values[r]: continue
      thisobserved.SetPoint(       thisobserved.GetN(),m,values[r][m]['obs'][0])
      if values[r][m]['PEs'] :
        expCLs = GetCenterAndSigmaDeviations(values[r][m]['PEs'])
        thisexpected_minus2.SetPoint(thisexpected_minus2.GetN(),m,expCLs[0])
        thisexpected_minus1.SetPoint(thisexpected_minus1.GetN(),m,expCLs[1])
        thisexpected.SetPoint(       thisexpected.GetN(),m,expCLs[2])
        thisexpected_plus1.SetPoint( thisexpected_plus1.GetN(),m,expCLs[3])
        thisexpected_plus2.SetPoint( thisexpected_plus2.GetN(),m,expCLs[4])

    allobserved.append(thisobserved)
    allexpected.append(thisexpected)

    thisexpected1Sigma = makeBand(thisexpected_minus1,thisexpected_plus1)
    thisexpected2Sigma = makeBand(thisexpected_minus2,thisexpected_plus2)

    #if there is just one ratio, we want bands
    if len(ratios) == 1:
      allexpected1Sigma.append(thisexpected1Sigma)
      allexpected2Sigma.append(thisexpected2Sigma)
    #if there are more than one, but this is the first, we want bands
    elif r == 0:
      allexpected1Sigma.append(thisexpected1Sigma)
      allexpected2Sigma.append(thisexpected2Sigma)
    #if there are more than one and this is not the first, we want no bands
    else:
      allexpected1Sigma.append(ROOT.TGraph())
      allexpected2Sigma.append(ROOT.TGraph())

    results[r] = {'obs':thisobserved, 'exp': thisexpected,'exp1':  thisexpected1Sigma,'exp2': thisexpected2Sigma}

  if len(ratios) == 1:
    outname += "_" + str(int(ratios[0]*100))

  #set channel-specific axis limits
  if "inc" in dataset:
    yAxisLimits = incYlims
  elif "2b" in dataset: 
    yAxisLimits = twobYlims
  else: 
   print  "Unknown dataset", dataset

  xAxisLimits = [180,1070]


    #definition: drawSeveralObservedExpectedLimits(self,observedlist,expectedlist,expected1list,expected2list,signallegendlist,name,nameX,nameY,luminosity,CME,xmin,xmax,ymin,ymax,extraLegendLines = [], doLogY=True,doLogX=False,doRectangular=False,doLegendLocation="Left",ATLASLabelLocation="BottomL",addHorizontalLines=[],cutLocation="Right",labels=[]) :

  print allobserved, allexpected

  myPainter.drawSeveralObservedExpectedLimits(allobserved,allexpected,allexpected1Sigma,allexpected2Sigma,[names['%1.2f'%r] for r in ratios],outname,"m_{G} [GeV]","#sigma #times #it{A} #times BR [pb]",luminosity,13,xAxisLimits[0],xAxisLimits[1],yAxisLimits[0],yAxisLimits[1],[],ATLASLabelLocation="TopL",cutLocation="Left", doLegendLocation="Right")
  
  return results

#===================================================================================

def makePlots(ratios, compare_list, names_list, makeComparisons=True, destination = "./plots/"):


  #make sure the destination dir exists
  if not os.path.isdir(destination):
    os.makedirs(destination)

  # single plots
  dataset = compare_list[0].split("/")[-1]
  dataset = dataset.split("_")[1]
  print dataset

  comparison_plots = {}

  luminosity = luminosity_dict[dataset]
  cutstring = "|y*| < 0.75"
  for filetype in compare_list :
    thisname = names_list[compare_list.index(filetype)]
    destination_single = destination+thisname+"/"
    if not os.path.isdir(destination_single):
      os.makedirs(destination_single)
    plots = getGaussianLimits(filetype, ratios, dataset, luminosity, cutstring, makePlot = True, outfolder = destination_single,doEntireRange = True)
    for ratio in plots.keys() :
      if not ratio in comparison_plots.keys() :
        comparison_plots[ratio] = {"graphs": [], "names": []}
      comparison_plots[ratio]["graphs"].append(plots[ratio]["obs"])
      comparison_plots[ratio]["names"].append(thisname)

# definition
#  def drawSeveralObservedLimits(observedlist,signallegendlist,name,nameX,nameY,luminosity,CME,xmin,xmax,ymin,ymax,extraLegendLines = [], doLogY=True,doLogX=False,doRectangular=False,doLegendLocation="Right",ATLASLabelLocation="BottomL",isTomBeingDumb=False,addHorizontalLines=[],pairNeighbouringLines=False,cutLocation="Right") :

  #set channel-specific axis limits
  if "inc" in dataset:
    yAxisLimits = incYlims
  elif "2b" in dataset: 
    yAxisLimits = twobYlims
  else: 
   print  "Unknown dataset", dataset

  outRoot = ROOT.TFile.Open(destination+"full_graphs.root","RECREATE")
  outRoot.cd()

  ratios_list_all = []
  ratios_names_list_all = []
  for ratio in comparison_plots.keys() :
    graphs_list = comparison_plots[ratio]["graphs"]
    names_list = comparison_plots[ratio]["names"]
    
    write_ratio = "{0}".format(ratio)
    write_ratio = write_ratio.replace(".","p")
    # Save
    for graph in graphs_list :
      graph.Write("graph_{0}_{1}".format(names_list[graphs_list.index(graph)],write_ratio))
    
    ratios_list = []
    ratios_names_list = []
    for index in range(1,len(graphs_list)) :
      print "looking at index",index
      print graphs_list[index]
      print graphs_list[0]
      ratio_wrt_first = ratio_graphs(graphs_list[index],graphs_list[0])
      ratios_list.append(ratio_wrt_first)
      ratios_names_list.append("{0}/{1}, {2}".format(names_list[index],names_list[0],ratio))
    ratios_list_all += ratios_list
    ratios_names_list_all += ratios_names_list

    # Compare values
    outname = destination+"/comparisons/raw_comparison_{0}".format(ratio)
    myPainter.drawSeveralObservedLimits(graphs_list,names_list,outname,"m_{G} [GeV]","#sigma #times #it{A} #times BR [pb]",luminosity,13,xAxisLimits[0],xAxisLimits[1],yAxisLimits[0],yAxisLimits[1],[],ATLASLabelLocation="TopL",cutLocation="Left", doLegendLocation="Right")

    # Make separate plot of ratios
    outname = destination+"/comparisons/ratios_{0}".format(ratio)
    myPainter.drawSeveralObservedLimits(ratios_list,ratios_names_list,outname,"m_{G} [GeV]","#sigma #times #it{A} #times BR [pb]",luminosity,13,xAxisLimits[0],xAxisLimits[1],0.95,1.05,[],doLogY=False,doLogX=False,addHorizontalLines=[1.0],ATLASLabelLocation="TopL",cutLocation="Left", doLegendLocation="Right")

  outRoot.Close()

  # Make separate plot of ratios (all)
  outname = destination+"/comparisons/ratios_all".format(ratio)
  myPainter.drawSeveralObservedLimits(ratios_list_all,ratios_names_list_all,outname,"m_{G} [GeV]","#sigma #times #it{A} #times BR [pb]",luminosity,13,xAxisLimits[0],xAxisLimits[1],0.96,1.1,[],doLogY=False,doLogX=False,addHorizontalLines=[1.0],ATLASLabelLocation="TopL",cutLocation="Left", doLegendLocation="Right")


if __name__ == "__main__":

  #==========================
  #   User configurables
  #==========================

  makeCombinedPlots = False

  dest = "./compare_systs/"

  indir = "/cluster/warehouse/kpachal/DijetISR/Resolved2017/LimitSetting/BayesianFramework/results/gaussian_limits/"
  compareForms = [
      indir + "/2018.08.31_total_noSyst/GenericGaussians_2bSing_*.root",
      indir + "/2018.08.31_total_JES/GenericGaussians_2bSing_*.root",
  ]
  compareNames = ["noSyst","JES"]

  makePlots([0.0, 0.07, 0.10, 0.15], compareForms, compareNames, makeComparisons=True, destination = dest)

  #==========================
  
