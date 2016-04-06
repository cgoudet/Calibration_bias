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
#include "RooGaussian.h"

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
    - unsigned int: iBin, statTree, indepDistorded, indepTemplates, runNumber, nBins, bootstrap, fitMethod, inversionMethod;
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
    ("checkDistri", po::value<unsigned int>(&m_checkDistri)->default_value(0),"")
    /*0= bias study
      1= errSigma distribution
     */
    
    ;
    
  po::variables_map vm;
  ifstream ifs( configFileName, ifstream::in );
  po::store(po::parse_config_file(ifs, configOptions), vm);
  po::notify(vm);  

  m_nHist=0;
}


//==================================================
//Clear attributes
BiasAnalysis::~BiasAnalysis()
{
  for (unsigned int iHist=0; iHist<m_nHist; iHist++)
    {
      for(unsigned int iVar=0; iVar<m_variablesStats.size()+2; iVar++)
  	{
  	  m_histStats[iHist][iVar]=-999;
  	}
    }

  m_mapHist.clear();
  m_mapHistPosition.clear();
  m_mapSumX.clear();
  m_mapSumXM.clear();
  m_mapNEff.clear();

  m_variablesBias.clear();
  m_variablesStats.clear();

  m_inTreeName.clear();
  cout<<"Cleaning ok"<<endl;

}



//===================================================
// Read files, link tree branches to local variables and sorting variables according to the ones selected.
// Fill the map m_mapHist with unique histogramms for each combination of all possible values of each variable.

void BiasAnalysis::SelectVariables(vector <string> dataFiles)
{
  map <string, unsigned int> mapUInt;
  map <string, double> mapDouble, mapMean;
  map <string, RooArgSet*> mapArgSet;
   
  TTree *inTree;
  TFile *inFile;
  unsigned int nEntries;
  double bias;

  string histName, rooName, errSigma;
  TString value, rooNum;

  TH1::AddDirectory(kFALSE);

  //1st loop over files: get min & max of each histogram, fill maps needed to compute mean an rms
  for (unsigned int iFile=0; iFile <dataFiles.size(); iFile++)
    { 
      inFile= TFile::Open(dataFiles[iFile].c_str());
      if (inFile == 0) { cout<<"Error: cannot open "<<dataFiles[iFile]<<" file\n"<<endl; return;}

      inTree = (TTree*) inFile->Get(m_inTreeName.c_str());  
      
      MapBranches mapBranches; 
      mapBranches.LinkTreeBranches(inTree);
      nEntries= inTree->GetEntries();
    
      for (unsigned int iEntry=0; iEntry<nEntries; iEntry++)
      	{
      	  inTree->GetEntry(iEntry);
	  mapDouble=mapBranches.GetMapDouble();
      	  mapUInt=mapBranches.GetMapUnsigned();
	  histName ="";

	  switch (m_checkDistri)
	    {
	    case 0:
	      { 
		bias = mapDouble.at("sigma")-mapDouble.at("inputC");
		break;
	      }
	    case 1:
	      {
		bias = mapDouble.at("errSigma");
		break;
	      }
	    default:
	      bias =  mapDouble.at("sigma")-mapDouble.at("inputC");
	    }
	
	  //	  if (mapUInt.at("nBins")!=6 || mapUInt.at("indepTemplates")==0 || mapUInt.at("bootstrap")==0 || mapUInt.at("indepDistorded")==0) continue;
	  if (mapUInt.at("nBins")!=6) continue;

	  for (unsigned int iVar =0; iVar < m_variablesBias.size(); iVar++)
      	    {
	      if (mapUInt.count(m_variablesBias[iVar])>0) value = TString::Format("%d",mapUInt.find(m_variablesBias[iVar])->second);
	      if (mapDouble.count(m_variablesBias[iVar])>0) value = TString::Format("%d",(int) floor ( (mapDouble.find(m_variablesBias[iVar])->second)*1e6) );
	      if (mapUInt.at("statTree")>=2E6)
	      {
		if (m_variablesBias[iVar] == "statTree") value = TString::Format("%d", 2774685);	      
	      }
      
	      if (iVar == m_variablesBias.size()-1) 
		{
		  histName+= m_variablesBias[iVar]+"_"+value;
		  if (m_mapNEff.count(histName) == 0)
		    {
		      m_mapXMin.insert(pair<string, double>(histName, bias));
		      m_mapXMax.insert(pair<string, double>(histName, bias));
		      
		      m_mapSumX.insert(pair<string, double>(histName, bias));
		      // m_mapSumXSquare.insert(pair<string, double>(histName, bias*bias));
		      m_mapNEff.insert(pair<string, unsigned int>(histName, 1));
		      m_mapHistPosition.insert(pair<string, unsigned int>(histName, m_nHist));
		           		      
		      m_nHist++;
		    }

		  if (m_mapNEff.count(histName) > 0)
		    {		     	      
		      if(bias<m_mapXMin[histName]) m_mapXMin[histName]=bias;
		      if(bias>m_mapXMax[histName]) m_mapXMax[histName]=bias;      
		      m_mapSumX[histName]+=bias;
		      //m_mapSumXSquare[histName]+=bias*bias;
		      m_mapNEff[histName]+=1;
		    }
		}
	      histName += m_variablesBias[iVar]+"_"+value+"_";
	      
	    }//end iVar (1st loop)

	}//end iEntry (1st loop)
    
    }//end iFile

  //2nd loop over files: fill m_mapHist     
  for (unsigned int iFile=0; iFile <dataFiles.size(); iFile++)
    { 
      inFile= TFile::Open(dataFiles[iFile].c_str());
      if (inFile == 0) { cout<<"Error: cannot open "<<dataFiles[iFile]<<" file\n"<<endl; return;}

      inTree = (TTree*) inFile->Get(m_inTreeName.c_str());  
      
      MapBranches mapBranches; 
      mapBranches.LinkTreeBranches(inTree);
      nEntries= inTree->GetEntries();
      
      for (unsigned int iEntry=0; iEntry<nEntries; iEntry++)
	{
	  inTree->GetEntry(iEntry);
	  mapDouble=mapBranches.GetMapDouble();
	  mapUInt=mapBranches.GetMapUnsigned();
	  histName ="";
	  

	  switch (m_checkDistri)
	    {
	    case 0:
	      { 
		bias = mapDouble.at("sigma")-mapDouble.at("inputC");
		break;
	      }
	    case 1:
	      {
		bias = mapDouble.at("errSigma");
		break;
	      }
	    default:
	      bias =  mapDouble.at("sigma")-mapDouble.at("inputC");
	    }

	  // if (mapUInt.at("nBins")!=6 || mapUInt.at("indepTemplates")==0 || mapUInt.at("bootstrap")==0 || mapUInt.at("indepDistorded")==0) continue;
	  if (mapUInt.at("nBins")!=6) continue;

	  for (unsigned int iVar =0; iVar < m_variablesBias.size(); iVar++)
	    {
	      if (mapUInt.count(m_variablesBias[iVar])>0) value = TString::Format("%d",mapUInt.find(m_variablesBias[iVar])->second);
	      if (mapDouble.count(m_variablesBias[iVar])>0) value = TString::Format("%d",(int) floor ( (mapDouble.find(m_variablesBias[iVar])->second)*1e6) );
	      if (mapUInt.at("statTree")>=2E6)
	      {
		if (m_variablesBias[iVar] == "statTree") value = TString::Format("%d", 2774685);	      
	      }
	      if (iVar == m_variablesBias.size()-1) 
		{
		  histName+= m_variablesBias[iVar]+"_"+value; 
		  
		  if(m_mapHist.count(histName)==0)
		    {
		      if (m_mapXMin[histName] == m_mapXMax[histName]) m_mapHist.insert(pair<string, TH1D*> (histName, new TH1D(histName.c_str(), "", 100, -0.1, 0.1)));
		      else m_mapHist.insert(pair<string, TH1D*> (histName, new TH1D(histName.c_str(), "", 100, m_mapXMin[histName], m_mapXMax[histName])));
		      m_mapHist[histName]->Sumw2();
		      m_mapHist[histName]->Fill(bias);
		      
		      mapMean.insert(pair <string, double> (histName, m_mapSumX[histName]/m_mapNEff[histName]));
		      m_mapSumXM.insert(pair<string, double> (histName, pow (bias-mapMean[histName], 2)) );
		      m_histNames.push_back(histName);

		      //Fill m_mapDataSet to fit using RooFit
		      rooName= "bias_"+histName;
		      m_mapBias.insert(pair<string, RooRealVar*>(histName, new RooRealVar(rooName.c_str(), "C^{meas}-C^{input}", m_mapXMin[histName], m_mapXMax[histName])));
		      m_mapBias[histName]->setVal(bias);
		      
		      rooName= "set_"+histName;
		      mapArgSet.insert(pair<string, RooArgSet*>(histName, new RooArgSet(rooName.c_str())));
		      mapArgSet[histName]->add(*m_mapBias[histName]);
		      
		      rooName= "data_"+histName;
		      m_mapDataSet.insert(pair<string, RooDataSet*> (histName, new RooDataSet(rooName.c_str(), "data", *mapArgSet[histName])));
		      m_mapDataSet[histName]->add(*mapArgSet[histName]);
		    }

		  if (m_mapHist.count(histName) > 0)
		    {
		      m_mapHist[histName]->Fill(bias);
		      m_mapSumXM[histName]+= pow (bias-mapMean[histName], 2);
		      m_mapBias[histName]->setVal(bias);
		      mapArgSet[histName]->add(*m_mapBias[histName]);
		      m_mapDataSet[histName]->add(*mapArgSet[histName]);
		    }
		}

	      histName += m_variablesBias[iVar]+"_"+value+"_";
	      
	    }//end iVar (2nd loop)
	}//end iEntry (2nd loop)
    }//end iFile


  inFile->Close(); //close file and delete tree
  delete inFile;
  return;
}



//=====================================================
//For each histogram, fill the 2D multi_array with:
// - 1st dim: histogram
// - 2nd dim: mean (for a given method), mean error, rms...
//Fill a csv file with those values.

void BiasAnalysis::MeasureBias(string outFileName, string outRootFileName)
{
  m_histStats.resize(extents[m_nHist][m_variablesStats.size()+2]);
  unsigned int iHist=0;
  unsigned int nBins, skip;
  double mean=0.;
  double errMean=0.;
  double rms=0; 
  double xMin=0.;
  double xMax=0.;
  string histName;
  char *token;

  TFile *outRootFile = new TFile(outRootFileName.c_str(), "RECREATE"); 

  ofstream outputFile(outFileName, ios::out);
  if (outputFile == 0) {cout<<"Error while opening outputFile"<<endl; return ;}
 

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
		    m_mapHist[histName]->Write();
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
		    m_mapHist[histName]->Write();
		  } 
		break;
	      }

	    case 2://get from gaussian fit
	      {
		cout << histName<<endl;
		nBins = m_mapHist[histName]->GetNbinsX();
		
		xMin = m_mapHist[histName]->GetXaxis()->GetBinCenter(2);
		xMax = m_mapHist[histName]->GetXaxis()->GetBinCenter(nBins-1);
		cout <<histName<<" " <<xMin << " "<< xMax<<endl;
		TF1 *f0 = new TF1("f0", "gaus", xMin, xMax);
		m_mapHist[histName]->Fit("f0","R");
		mean = m_mapHist[histName]->GetFunction("f0")->GetParameter(1);
		rms = m_mapHist[histName]->GetFunction("f0")->GetParameter(2);
		if (mean-1.5*rms>=xMin) xMin= mean-1.5*rms;
		if (mean+1.5*rms<=xMax && mean+1.5*rms>xMin) xMax= mean+1.5*rms;
		cout <<histName<<" " <<xMin << " "<< xMax<<endl;
		TF1 *f1= new TF1("f1", "gaus", xMin, xMax);
		m_mapHist[histName]->Fit("f1","R");
		
		if (m_variablesStats[iVar]== 0)
		  {
		    mean = m_mapHist[histName]->GetFunction("f1")->GetParameter(1);
		    rms = m_mapHist[histName]->GetFunction("f1")->GetParameter(2);
		    errMean = m_mapHist[histName]->GetFunction("f1")->GetParError(1); 
		  }
		m_mapHist[histName]->Write();
		delete f1; f1=0;
		delete f0; f0=0;
		break;
	      }
	    case 3://gauss RooFit
	      {
		RooRealVar* meanFit = new RooRealVar("meanFit", "mean", -1, 1);
		RooRealVar* sigmaFit= new RooRealVar("sigmaFit", "sigma",0, 1);
		m_mapGauss.insert(pair<string, RooGaussian*>(histName, new RooGaussian("gauss", "gauss", *m_mapBias[histName], *meanFit, *sigmaFit) ));

		m_mapGauss[histName]->fitTo(*m_mapDataSet[histName], Range( m_mapXMin[histName]+(m_mapXMax[histName]-m_mapXMin[histName])*0.02, m_mapXMax[histName] ));
		rms= sigmaFit->getValV();
		mean= meanFit->getValV();
		m_mapGauss[histName]->fitTo(*m_mapDataSet[histName], Range(mean-1.5*rms, mean+1.5*rms));
		m_mapHist[histName]->Write();
		delete meanFit;
		delete sigmaFit;
		break;
	      }
	    }//end switch

	  m_histStats[iHist][0]=mean;
	  m_histStats[iHist][1]=rms;
	  m_histStats[iHist][2]=errMean;
	  //m_histStats[iHist][3]=m_mapNEff[histName];
  	}//end iVar

      //writing the cvs file
      if (iHist==0) 
	{
	  outputFile<<"Histogram name"<<","<<"Histogram index"<<","<<"Number of entries"<<",";
	  for (unsigned int iVarBias=0; iVarBias<m_variablesBias.size(); iVarBias++)
	    {
	      outputFile<<m_variablesBias[iVarBias]<<",";
	    }
	  outputFile<<"Mean"<<","<<"RMS"<<","<<"Error mean"<<"\n";
	}

      outputFile<<histName<<","<<iHist<<","<<m_mapNEff[histName]<<",";
      token = strtok((char*)histName.c_str(), "_");
      skip=1;
      while(token !=NULL)
      	{
	  if (skip>m_nHist) break;
      	  if (skip%2==0) outputFile<<token<<",";
	  skip++;
	  token=strtok(NULL, "_");
      	}
           
      outputFile<<mean<<","<<rms<<","<<errMean<<"\n";
      
      //next histogram
      it++;
    }//end iteration over histograms (while loop)

  cout<<"End of measure"<<endl;
  outRootFile->Close();
  delete outRootFile;
  return;
}





//==================================================
//Draw plots and save them into a pdf file

void BiasAnalysis::MakePlots(string path, string latexFileName)
{
  //Prepare latex file to store plots 
  fstream stream;
  string latexTitle = "Bias study";
  stream.open( (path+latexFileName).c_str(), fstream::out | fstream::trunc );
  WriteLatexHeader( stream, latexTitle , "Antinea Guerguichon" );

  //Draw plots
  //vector <string> vectHistNames;
  vector <string> vectStatNames;
  if (m_methodStats == 1) vectStatNames.push_back("Mean hist");
  if (m_methodStats == 2) vectStatNames.push_back("Mean fit");
  if (m_methodStats!= 1 && m_methodStats !=2) vectStatNames.push_back("Mean");
  vectStatNames.push_back("RMS");
  vectStatNames.push_back("Error mean");

  vector <string> vectOptDraw;//Options to draw the histograms (cf DrawPlot.cxx

  string histName;
  unsigned int iHist=0;
  TString statVal, statPos;
  string legLatex;

  for (unsigned int i=0; i<m_histNames.size(); i++)
    {
      histName = m_histNames[i];
      m_histNames[i]=path+m_histNames[i];

      for (unsigned int i=0; i<m_variablesStats.size()+2; i++)
      	{
      	  statVal= TString::Format("%f", m_histStats[iHist][i]);
	  legLatex= "latex="+ vectStatNames[i]+ ": " +statVal;
      	  vectOptDraw.push_back(legLatex.c_str());
	  statPos= TString::Format("%f", 0.85-i*0.05);
	  legLatex= "latexOpt= 0.7 "+ statPos;
	  vectOptDraw.push_back(legLatex.c_str());
      	}
      
      vectOptDraw.push_back("yTitle=#Events");
      legLatex= "latex="+histName;
      vectOptDraw.push_back(legLatex.c_str());
      vectOptDraw.push_back("latexOpt= 0.1 0.9");
      vectOptDraw.push_back("extendUp= 0.4");
      
      if (m_methodStats == 3)  DrawPlot(m_mapBias[histName], {m_mapDataSet[histName], m_mapGauss[histName]}, path+histName,{vectOptDraw} );
      
      switch (m_checkDistri)
	{
	case 0:
	  {
	    vectOptDraw.push_back("xTitle=C^{meas}-C^{input}");
	    break;
	  }
	case 1:
	  {
	    vectOptDraw.push_back("xTitle=C^{meas} error");
	    break;
	  }
	default:
	  vectOptDraw.push_back("xTitle=C^{meas}-C^{input}");
	}


      
      if (m_methodStats!=3) DrawPlot({m_mapHist[histName]}, path+histName, {vectOptDraw});       
      vectOptDraw.clear();
      
      iHist++;
    }

  //  Store plots into the file
  //stream << "\\section{Method to get stats: "<<m_methodStats<<"}"<< endl;
  
  stream << "Tree: "<< m_inTreeName<<"\\newline  "<<endl;
  stream << "\\indent Variables: ";
  for (unsigned int iVar=0; iVar< m_variablesBias.size(); iVar++)
    {
      if (iVar == m_variablesBias.size()-1) stream<<m_variablesBias[iVar] <<"\\newline  ";
      else  stream  << m_variablesBias[iVar] <<", ";
    }

  if (m_checkDistri ==1) stream << "\\indent Check for errSigma distribution with bootstrap=0, indepTemplates=0, indepDistorded=1 \\newline"<<endl;
  WriteLatexMinipage( stream, m_histNames, 2, true );
  stream << "\\end{document}" << endl;
  string commandLine = "pdflatex  -interaction=batchmode " + path+latexFileName;
  system( commandLine.c_str() );
  system( commandLine.c_str() );
  system( commandLine.c_str() );

  commandLine = "rm " + path+ m_variablesBias[0]+ "*";
  system( commandLine.c_str() );

 
  cout<<"Plots drawn and stored into a pdf file"<<endl;
  return;
}
