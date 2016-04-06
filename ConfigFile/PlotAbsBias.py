import ROOT
import csv
from ROOT import *

path="/sps/atlas/a/aguerguichon/Calibration/Bias/"
rootFile= path+"RootFiles/BiasBinFit.root"
openedFile= TFile.Open(rootFile)

rootFileName="rootFileName="+rootFile+"\n"

statFile=open(path+"Stats/BiasBinFit.csv", "rb")

reader= csv.reader(statFile)

for iBin in range(0, 6):
    configFile = open("AbsBiasBin_"+str(iBin)+".boost", "w") 
    configFile.write("inputType=0 \n")
    for row in reader:                
        objName = "inputC_7000_statTree_2774685_indepTemplates_0_indepDistorded_1_bootstrap_0_iBin_" +str(iBin)
        if openedFile.Get(objName):
            if row[0]==objName:
                configFile.write( rootFileName )
                configFile.write( "objName="+objName+"\n" )
                configFile.write( "legend= 2.7M, Mean="+str(row[9])+", RMS="+row[10]+"\n" )
                    
                                    
    configFile.write("legendPos= 0.45 0.9 \n" )
    configFile.write("normalize=1 \n")            
    configFile.write("latex=Bin "+ str(iBin)+"\n")
    configFile.write("latexOpt= 0.15 0.85 \n")
    configFile.write("xTitle=c^{meas}-c^{input} \n")
    configFile.write("yTitle=Number of events \n")
    configFile.write("plotDirectory="+path+"Plots/ \n")
    configFile.write("extendUp=0.4 \n")
    configFile.close()
   # statFile.seek(0)

statFile.close()
