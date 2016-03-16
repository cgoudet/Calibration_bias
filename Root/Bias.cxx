#include <iostream>
#include <string>
#include <vector>

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"

#include "PlotFunctions/DrawPlot.h"
#include "Calibrations_bias/Bias.h"

using namespace std;


//Constructor: initialisation of the array of histograms
Bias::Bias(unsigned int NBins)
{
  nBinUp = NBins;
  nEff[100][100]=0; //number of entries for a given histogram
}
//===================
double Bias::SetArrayHist()
{
  char *histName = new char[100];                                     
  for (unsigned int i=0; i<nBinUp; i++)                               
    {                                                                 
      for (unsigned int j=0; j<nBinUp; j++)                           
	{                                                             
	  sprintf(histName, "hBiasConf_%d_%d",i,j);                    
	  arHisto[i][j]= new TH1D(histName, "", 100, -0.05, 0.05);    
	}                                                             
    }
}

//====================
double Bias::GetBias (double input, double measure, vector <double> optVariables)
{ 
  return measure-input;
}


//====================
//For each configuration, histogram of #Events vs bias is drawn with the data from each root file (for the tree ConfigurationCTree)
void Bias::MakeHist (string ToyRootFile)
{
  TFile *f= TFile::Open(ToyRootFile);
  if (f == 0) { cout<<"Error: cannot open file "<<ToyRootFile<<"\n"<<endl; return 0;}
     
  TTree *t = (TTree*) f->Get("ConfigurationsCTree");

  const int nEntries= t->GetEntries();
  unsigned int iConf, jConf, statTree, nBins;
  double sigma, errSigma, inputC, dataRMS;

  t->SetBranchAddress("iConf", &iConf);
  t->SetBranchAddress("jConf", &jConf);
  t->SetBranchAddress("nBins", &nBins);
  t->SetBranchAddress("statTree", &statTree);
  t->SetBranchAddress("sigma", &sigma);
  t->SetBranchAddress("errSigma", &errSigma);
  t->SetBranchAddress("inputC", &inputC);
  t->SetBranchAddress("dataRMS", &dataRMS);

  double bias=0.;
  xMin[100][100]=10.;
  xMax[100][100]=-10.;
  statTreeVar[100][100]=0;
  inputCVar[100][100]=0.;
 
  for (int n=0; n<nEntries; n++)
    {
      t->GetEntry(n);
      if (nBins!=nBinUp) continue;
      bias= sigma-inputC;
      statTreeVar[iConf][jConf] = statTree;
      inputCVar[iConf][jConf]= inputC;
      arHisto[iConf][jConf]->Fill(bias);
      nEff[iConf][jConf]++;
      if (bias < xMin[iConf][jConf]) xMin[iConf][jConf]=bias;
      if (bias > xMax[iConf][jConf]) xMax[iConf][jConf]=bias;
    }
  f->Close();
  printf("\n File %s : OK", ToyRootFile);
}


