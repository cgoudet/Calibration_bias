import ROOT
import csv
from ROOT import *

path="/sps/atlas/a/aguerguichon/Calibration/Bias/"
rootFile= path+"RootFiles/BiasBinFit.root"
openedFile= TFile.Open(rootFile)

rootFileName="rootFileName="+rootFile+"\n"
configFile = open("CompareBiasBin.boost", "w")
configFile.write("inputType=0 \n")
 
for iBin in range(0, 6):
    objName = "histBias"+str(iBin)
    if openedFile.Get(objName):
        configFile.write( rootFileName )
        configFile.write( "objName="+objName+"\n" )
        configFile.write( "legend= Bin "+str(iBin)+"\n")
                    
                    
configFile.write("legendPos= 0.65 0.9 \n" )            
configFile.write("latex=Input: 0.007\n")
configFile.write("latexOpt= 0.3 0.85 \n")    
configFile.write("xTitle=Statistics (10^{3}) \n")
configFile.write("yTitle=Bias \n")
configFile.write("plotDirectory="+path+"Plots/ \n")
configFile.write("extendUp=0.2 \n")
configFile.write("rangeUserX= 0 1100 \n")
configFile.write("rangeUserY= -0.006 0.007\n")
#configFile.write("drawStyle=3 \n")
configFile.close()
