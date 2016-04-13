import ROOT
import csv
from ROOT import *

path="/sps/atlas/a/aguerguichon/Calibration/Bias/"
rootFile= path+"RootFiles/BiasConfFit.root"
openedFile= TFile.Open(rootFile)

rootFileName="rootFileName="+rootFile+"\n"

statFile=open(path+"Stats/BiasConfFit.csv", "rb")

reader= csv.reader(statFile)

for iConf in range(0, 6):
    for jConf in range (0, iConf+1):
        configFile = open("CompareStatConf_"+str(iConf)+"_"+str(jConf)+".boost", "w") 
        configFile.write("inputType=0 \n")
        for row in reader:
            objName = "inputC_7000_statTree_100000_indepTemplates_1_indepDistorded_1_bootstrap_1_iConf_" +str(iConf)+ "_jConf_"+str(jConf)
            if openedFile.Get(objName):
                if row[0]==objName:
                    configFile.write( rootFileName )
                    configFile.write( "objName="+objName+"\n" )
                    configFile.write( "legend= 100k, Mean="+str(row[10])+", RMS="+row[11]+"\n" )
                    
            objName = "inputC_7000_statTree_1000000_indepTemplates_1_indepDistorded_1_bootstrap_1_iConf_" +str(iConf)+ "_jConf_"+str(jConf)
            if openedFile.Get(objName):
                if row[0]==objName:
                    configFile.write( rootFileName )
                    configFile.write( "objName="+objName+"\n" )
                    configFile.write( "legend= 1M, Mean="+str(row[10])+", RMS="+row[11]+"\n" )
                    
            objName = "inputC_7000_statTree_2774685_indepTemplates_1_indepDistorded_1_bootstrap_1_iConf_" +str(iConf)+ "_jConf_"+str(jConf)
            if openedFile.Get(objName):
                if row[0]==objName:
                    configFile.write( rootFileName )
                    configFile.write( "objName="+objName+"\n" )
                    configFile.write( "legend= 2.7M, Mean="+str(row[10])+", RMS="+row[11]+"\n" )
                    
                     
        configFile.write("legendPos= 0.45 0.9 \n" )
        configFile.write("normalize=1 \n")            
        configFile.write("latex=Configuration ("+str(iConf)+", "+str(jConf)+")\n")
        configFile.write("latexOpt= 0.15 0.85 \n")
        configFile.write("latex=Input: 0.007\n")
        configFile.write("latexOpt= 0.15 0.8 \n")
        configFile.write("xTitle=c^{meas}-c^{input} \n")
        configFile.write("yTitle=Number of events \n")
        configFile.write("plotDirectory="+path+"Plots/ \n")
        configFile.write("extendUp=0.4\n")
        configFile.close()
        statFile.seek(0)

statFile.close()
