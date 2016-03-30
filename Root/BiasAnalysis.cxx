#include <iostream>
#include <istream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <boost/program_options.hpp>
#include <boost/multi_array.hpp>

#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TString.h"

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooArgSet.h"

#include "Calibration_bias/BiasAnalysis.h"
#include "PlotFunctions/DrawPlot.h"
#include "PlotFunctions/SideFunctionsTpp.h"
#include "PlotFunctions/SideFunctions.h"
#include "PlotFunctions/MapBranches.h"

using namespace std;
using namespace RooFit;
namespace po = boost::program_options;

//====================================================
//Read configuration options from a configuration file
BiasAnalysis::BiasAnalysis(string configFileName)
{  
  po::options_description configOptions("Configuration options");
  configOptions.add_options()
    ("help", "Display this help message")
    ("variablesBias", po::value<vector<string>>(&m_variablesBias)->multitoken(),"Variables to study bias")
    /*accepted variables for ConfugurationsCTree: 
       - double: sigma, errSigma, inputC, dataRMS, nOptim;
       - unsigned int: iConf, jConf, statConf, statTree, indepDistorded, indepTemplates, runNumber, nBins, bootstrap, fitMethod;
    accepted variables for scalesTree:
    - double: sigma, errSigma, inputC, dataRMS, nOptim;
    - unsigned int: statTree, indepDistorded, indepTemplates, runNumber, nBins, bootstrap, fitMethod, inversionMethod;

    */
    ("variablesStats", po::value<vector<unsigned int>>(&m_variablesStats)->multitoken(),"Variables for the csv file")
    /*0: mean, rms, errMean
     */
    ("methodStats", po::value<unsigned int>(&m_methodStats), "")
    /*0: compute mean and rms
      1: h->GetMean(), h->GetRMS()
      2: get mean given by the Gaussian fit of the histogram
      3: RooFit gauss
    */
    ("selectTree", po::value<string>(&m_inTreeName), "")
    /*ConfigurationsCTree or scalesTree
     */
    ;
    
  po::variables_map vm;
  ifstream ifs( configFileName, ifstream::in );
  po::store(po::parse_config_file(ifs, configOptions), vm);
  po::notify(vm);  

  m_nHist=0;
}

BiasAnalysis::~BiasAnalysis(){}

//===================================================
// Read files, link tree branches to local variables and sorting variables according to the ones selected.
// Fill the map m_mapHist with unique histogramms for each combination of all possible values of each variable.

void BiasAnalysis::SelectVariables(vector <string> dataFiles)
{
  map <string, unsigned int> mapUInt;
  map <string, double> mapDouble, mapMean, mapXMin, mapXMax;

  TTree *inTree;

  unsigned int nEntries;
  double bias;

  string histName;
  TString value;
  

  //1st loop over files: get min & max of each histogram, fill maps needed to compute mean an rms
  for (unsigned int iFile=0; iFile <dataFiles.size(); iFile++)
    { 
      TFile *f= TFile::Open(dataFiles[iFile].c_str());
      if (f == 0) { cout<<"Error: cannot open "<<dataFiles[iFile]<<" file\n"<<endl; return;}

      inTree = (TTree*) f->Get(m_inTreeName.c_str());  
      
      MapBranches mapBranches; 
      mapBranches.LinkTreeBranches(inTree);
      nEntries= inTree->GetEntries();
    
      for (unsigned int iEntry=0; iEntry<nEntries; iEntry++)
      	{
      	  inTree->GetEntry(iEntry);
	  mapDouble=mapBranches.GetMapDouble();
      	  mapUInt=mapBranches.GetMapUnsigned();
	  histName ="";
	  bias = mapDouble.at("sigma")-mapDouble.at("inputC");
	  
	  if (mapUInt.at("nBins")!=6) continue;

	  for (unsigned int iVar =0; iVar < m_variablesBias.size(); iVar++)
      	    {
	      if (mapUInt.count(m_variablesBias[iVar])>0) value = TString::Format("%d",mapUInt.find(m_variablesBias[iVar])->second);
	      if (mapDouble.count(m_variablesBias[iVar])>0) value = TString::Format("%d",(int) floor ( (mapDouble.find(m_variablesBias[iVar])->second)*1e6) );
      
	      if (iVar == m_variablesBias.size()-1) 
		{
		  histName+= m_variablesBias[iVar]+"_"+value;
		  if (m_mapNEff.count(histName) == 0)
		    {
		      mapXMin.insert(pair<string, double>(histName, bias));
		      mapXMax.insert(pair<string, double>(histName, bias));
		      
		      m_mapSumX.insert(pair<string, double>(histName, bias));
		      // m_mapSumXSquare.insert(pair<string, double>(histName, bias*bias));
		      m_mapNEff.insert(pair<string, unsigned int>(histName, 1));
		      m_mapHistPosition.insert(pair<string, unsigned int>(histName, m_nHist));
		      m_nHist++;
		    }

		  if (m_mapNEff.count(histName) > 0)
		    {		     	      
		      if(bias<mapXMin[histName]) mapXMin[histName]=bias;
		      if(bias>mapXMax[histName]) mapXMax[histName]=bias;      
		      m_mapSumX[histName]+=bias;
		      //m_mapSumXSquare[histName]+=bias*bias;
		      m_mapNEff[histName]+=1;
		    }
		}
	      histName += m_variablesBias[iVar]+"_"+value+"_";
	      
	    }//end iVar (1st loop)

	}//end iEntry (1st loop)

    }//end iFile
  
  delete inTree;

  //2nd loop over files: fill m_mapHist     
  for (unsigned int iFile=0; iFile <dataFiles.size(); iFile++)
    { 
      TFile *f= TFile::Open(dataFiles[iFile].c_str());
      if (f == 0) { cout<<"Error: cannot open "<<dataFiles[iFile]<<" file\n"<<endl; return;}

      inTree = (TTree*) f->Get(m_inTreeName.c_str());  
      
      MapBranches mapBranches; 
      mapBranches.LinkTreeBranches(inTree);
      nEntries= inTree->GetEntries();

      for (unsigned int iEntry=0; iEntry<nEntries; iEntry++)
	{
	  inTree->GetEntry(iEntry);
	  mapDouble=mapBranches.GetMapDouble();
	  mapUInt=mapBranches.GetMapUnsigned();
	  histName ="";
	  bias = mapDouble.at("sigma")-mapDouble.at("inputC");

	  if (mapUInt.at("nBins")!=6) continue;

	  for (unsigned int iVar =0; iVar < m_variablesBias.size(); iVar++)
	    {
	      if (mapUInt.count(m_variablesBias[iVar])>0) value = TString::Format("%d",mapUInt.find(m_variablesBias[iVar])->second);
	      if (mapDouble.count(m_variablesBias[iVar])>0) value = TString::Format("%d",(int) floor ( (mapDouble.find(m_variablesBias[iVar])->second)*1e6) );
	      
	      if (iVar == m_variablesBias.size()-1) 
		{
		  histName+= m_variablesBias[iVar]+"_"+value; 
		  if(m_mapHist.count(histName)==0)
		    {
		      m_mapHist.insert(pair<string, TH1D*> (histName, new TH1D(histName.c_str(), "", 100, mapXMin[histName], mapXMax[histName])));

		      m_mapHist[histName]->Fill(bias);
		      
		      mapMean.insert(pair <string, double> (histName, m_mapSumX[histName]/m_mapNEff[histName]));
		      m_mapSumXM.insert(pair<string, double> (histName, pow (bias-mapMean[histName], 2)) ); 
		    }

		  if (m_mapHist.count(histName) > 0)
		    {
		      m_mapHist[histName]->Fill(bias);
		      m_mapSumXM[histName]+= pow (bias-mapMean[histName], 2);
		    }

		}

	      histName += m_variablesBias[iVar]+"_"+value+"_";
	      
	    }//end iVar (2nd loop)
	}//end iEntry (2nd loop)
      
    }//end iFile
  
  delete inTree;
  //cout <<"End of selection"<<endl;
  return;
}



//=====================================================
//For each histogram, fill the 2D multi_array with:
// - 1st dim: histogram
// - 2nd dim: mean (for a given method), mean error, rms...
//Fill a csv file with those values.

void BiasAnalysis::MeasureBias(string outFileName)
{
  m_histStats.resize(extents[m_nHist][m_variablesStats.size()+2]);
  unsigned int iHist=0;
  unsigned int nBins;
  double mean, errMean, rms, xMin, xMax;
  string histName;

  ofstream outputFile(outFileName, ios::out);
  if (outputFile == 0) {cout<<"Error while opening outputFile"<<endl; return ;}

  outputFile <<"Histogram"<<","<<"Mean"<<","<<"RMS"<<","<<"Error mean"<<endl;
  

  map <string, unsigned int>::iterator it=m_mapHistPosition.begin();
  while(it != m_mapHistPosition.end())
    {
      histName= it->first;
      iHist= it->second;
      for (unsigned int iVar=0; iVar<m_variablesStats.size(); iVar++)
  	{
	  switch (m_methodStats)
	    {

	    case 0://compute mean and rms
	      {
		if (m_variablesStats[iVar] == 0) 
		  {
		    mean= m_mapSumX[histName]/m_mapNEff[histName];
		    rms= sqrt( m_mapSumXM[histName]/m_mapNEff[histName] );
		    errMean= rms/sqrt(m_mapNEff[histName]);
		  }
		break;
	      }

	    case 1://get mean and rms from histogram
	      {
		if (m_variablesStats[iVar]== 0)
		  {
		    mean= m_mapHist[histName]->GetMean();
		    rms= m_mapHist[histName]->GetRMS();
		    errMean= rms/sqrt(m_mapNEff[histName]);
		  } 
		break;
	      }

	    case 2://get from gaussian fit
	      {
		nBins = m_mapHist[histName]->GetNbinsX();
		
		xMin = m_mapHist[histName]->GetXaxis()->GetBinCenter(2);
		xMax = m_mapHist[histName]->GetXaxis()->GetBinCenter(nBins-1);
		//cout<<xMax<<endl;

		//cout <<xMin <<"xMax last bin"<< xMax<<endl;
		TF1 *f1 = new TF1("f1", "gaus", xMin, xMax);
		m_mapHist[histName]->Fit("f1","R");
		if (m_variablesStats[iVar]== 0)
		  {
		    mean = m_mapHist[histName]->GetFunction("f1")->GetParameter(1);
		    rms = m_mapHist[histName]->GetFunction("f1")->GetParameter(2);
		    errMean = m_mapHist[histName]->GetFunction("f1")->GetParError(1); 
		  }
		delete f1; f1=0;
		break;
	      }
	    }//end switch

	  m_histStats[iHist][0]=mean;
	  m_histStats[iHist][1]=rms;
	  m_histStats[iHist][2]=errMean;
	  
  	}//end iVar
      
	  outputFile << histName <<","<<mean<<","<<rms<<","<<errMean<<endl;
      
      it++;
    }//end iteration over histograms (while loop)

  cout<<"End of measure"<<endl;
  return;
}





//==================================================
//Draw plots and save them into a pdf file

void BiasAnalysis::MakePlots(string latexFileName)
{
  //Prepare latex file to store plots 
  fstream stream;
  stream.open( latexFileName.c_str(), fstream::out | fstream::trunc );
  WriteLatexHeader( stream, "Test", "Antinea Guerguichon" );

  //Draw plots
  vector <string> vectHistNames;
  vector <string> vectStatNames;
  vectStatNames.push_back("Mean");
  vectStatNames.push_back("RMS");
  vectStatNames.push_back("Error mean");

  vector <string> vectOptDraw;//Options to draw the histograms (cf DrawPlot.cxx

  string histName;
  unsigned int iHist=0;
  TString statVal, statPos;
  string statLatex;

  map <string, TH1D*>::iterator it=m_mapHist.begin();
  while(it != m_mapHist.end())
    {
      histName = it->first;
      vectHistNames.push_back(histName);

      for (unsigned int i=0; i<m_variablesStats.size()+2; i++)
      	{
      	  statVal= TString::Format("%f", m_histStats[iHist][i]);
	  statLatex= "latex="+ vectStatNames[i]+ ": " +statVal;
      	  vectOptDraw.push_back(statLatex.c_str());
	  statPos= TString::Format("%f", 0.9-i*0.05);
	  statLatex= "latexOpt= 0.7 "+ statPos;
	  vectOptDraw.push_back(statLatex.c_str());
      	}
      
      vectOptDraw.push_back("yTitle=#Events");
      vectOptDraw.push_back("xTitle=C^{meas}-C^{input}");
                  
      DrawPlot({it->second}, histName, {vectOptDraw});
      
      vectOptDraw.clear();
      
      iHist++;
      it++;
    }

  //  Store plots into the file
  //stream << "\\section{Method to get stats: "<<m_methodStats<<"}"<< endl;
  WriteLatexMinipage( stream, vectHistNames, 2, true );
  vectHistNames.clear();
  stream << "\\end{document}" << endl;
  string commandLine = "pdflatex  -interaction=batchmode " + latexFileName;
  system( commandLine.c_str() );
  system( commandLine.c_str() );
  system( commandLine.c_str() );

  // commandLine = "rm " + m_variablesBias[0]+ "*";
  // system( commandLine.c_str() );

  cout<<"Plots drawn and stored into a pdf file"<<endl;
  return;
}
