import ROOT
import csv
from ROOT import *

path="/sps/atlas/a/aguerguichon/Calibration/Bias/"
rootFile= path+"RootFiles/BiasConf.root"
openedFile= TFile.Open(rootFile)

rootFileName="rootFileName="+rootFile+"\n"
configFile = open("CompareBiasConf_hist.boost", "w")
configFile.write("inputType=0 \n")
 
configurations=[ "(0, 0)", "(1, 0)", "(1, 1)", "(2, 0)", "(2, 1)", "(2, 2)", "(3, 2)", "(3, 3)", "(4, 3)", "(4, 4)", "(5, 3)", "(5, 4)", "(5, 5)"  ]

for iConf in range(0, 13):
    objName = "histBias"+str(iConf)
    if openedFile.Get(objName):
        configFile.write( rootFileName )
        configFile.write( "objName="+objName+"\n" )
        configFile.write( "legend= Conf"+configurations[iConf] +"\n")
                    
                    
configFile.write("legendPos= 0.65 0.9 \n" )            
configFile.write("latex=Input: 0.007\n")
configFile.write("latexOpt= 0.3 0.85 \n")    
configFile.write("xTitle=Statistics (10^{3}) \n")
configFile.write("yTitle=Bias \n")
configFile.write("plotDirectory="+path+"Plots/ \n")
configFile.write("extendUp=0.2 \n")
configFile.write("rangeUserX= 0 1100")
#configFile.write("drawStyle=3 \n")
configFile.close()
