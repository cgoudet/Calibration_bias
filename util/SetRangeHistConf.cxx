#include <iostream>
#include <istream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <sstream>

#include "TString.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"

#include "Calibration_bias/Analysis.h"
#include "PlotFunctions/DrawPlot.h"
#include "PlotFunctions/SideFunctionsTpp.h"

#include <boost/program_options.hpp>
#include <boost/multi_array.hpp>


using namespace std;

int main(int argc, char *argv[])
{ 
 
  //Set the number of bins
  unsigned int nBinUp = 6;
  
  //Declare an array of histo
  typedef boost::multi_array <TH1D*,2> maHist;
  maHist hist(boost::extents[nBinUp][nBinUp]);

  TString histName;
  for (unsigned int i=0; i<nBinUp; i++)
    {
      for (unsigned int j=0; j<nBinUp; j++)
	{
	  histName = TString::Format("hBiasConf_%d_%d",i,j);
	  hist[i][j]= new TH1D(histName, "", 100, -0.05, 0.05);
	}
    }

  
  //Open file, get tree and branches
  typedef boost::multi_array <unsigned int, 2> maUInt;
  maUInt statTreeVar(boost::extents[nBinUp][nBinUp]);
  maUInt nEff(boost::extents[nBinUp][nBinUp]);

  typedef boost::multi_array <double, 2> maDouble;
  maDouble inputCVar(boost::extents[nBinUp][nBinUp]);
  maDouble xMin(boost::extents[nBinUp][nBinUp]);
  maDouble xMax(boost::extents[nBinUp][nBinUp]);
  
  fill( nEff.origin(), nEff.origin() + nEff.size(), 0 );

  int k=1;

  while(argv[k]!=NULL)
    {
      TFile *f= TFile::Open(argv[k]);
      if (f == 0) { cout<<"Error: cannot open root file\n"<<endl; return 0;}
     
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
      fill( xMin.origin(), xMin.origin() + xMin.size(), 0 );
      fill( xMax.origin(), xMax.origin() + xMax.size(), -10 );
      fill( statTreeVar.origin(), statTreeVar.origin() + statTreeVar.size(), 0 );
      fill( inputCVar.origin(), inputCVar.origin() + inputCVar.size(), 0. );

      //Fill the array of histograms 
      for (int n=0; n<nEntries; n++)
	{
	  t->GetEntry(n);
	  if (nBins!=nBinUp) continue;
	  bias= sigma-inputC;
	  statTreeVar[iConf][jConf] = statTree;
	  inputCVar[iConf][jConf]= inputC;
	  hist[iConf][jConf]->Fill(bias);
	  nEff[iConf][jConf]++;
	  if (bias < xMin[iConf][jConf]) xMin[iConf][jConf]=bias;
	  if (bias > xMax[iConf][jConf]) xMax[iConf][jConf]=bias;
	}
      
      f->Close();
      printf("\n File %s : OK", argv[k]);
      k++;
    }


  //Draw the histograms
  stringstream iSs[10];
  stringstream jSs[10];
  string iStr[10];
  string jStr[10];
  string namePlot;
  
  vector <string> vecOptions;//Options to draw the histograms (cf DrawPlot.cxx)
  vecOptions.push_back("yTitle=#Events");
  vecOptions.push_back("xTitle=C^{meas}-C^{input}");
  
  Analysis myAnalysis;

  string path = "/sps/atlas/a/aguerguichon/Calibration/Bias/";
  
  for (unsigned int i=0; i<nBinUp; i++)
    {
      for (unsigned int j=0; j<=i; j++)
	{  
	  iSs[i]<<i;
	  iSs[i]>>iStr[i];
	  jSs[j]<<j;
	  jSs[j]>>jStr[j];
	  if (hist[i][j]->GetEntries() == 0) continue;
	  hist[i][j]->GetXaxis()->SetRangeUser(xMin[i][j], xMax[i][j]);
	  namePlot = path+"Plots/biasConf_" + iStr[i]+"_" +jStr[j];
	  DrawPlot( {hist[i][j]}, namePlot,{vecOptions});
	  
	}
    }
 
  string outFileName = path + "Stats/Stats.txt";
  myAnalysis.WriteFileStats(hist, nEff, statTreeVar, inputCVar, nBinUp ,outFileName );

  //End of program
  return 0;
}