#include <iostream>
#include <istream>
#include <fstream>
#include <map>
#include <vector>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include <string>
#include <sstream>
#include "PlotFunctions/DrawPlot.h"

using namespace std;

int main(int argc, char *argv[])
{
  
  //Declare an array of histo
  unsigned int nBinUp = 6;
  TH1D* arHisto[100][100];
  char *histName = new char[100];
  for (unsigned int i=0; i<nBinUp; i++)
    {
      for (unsigned int j=0; j<nBinUp; j++)
	{
	  sprintf(histName, "hBiasConf_%d%d",i,j);
	  arHisto[i][j]= new TH1D(histName, "", 100, -0.05, 0.05);
	}
    }

  //File to store mean, RMS and mean error of the histograms
  string nBinUpStr;
  stringstream nBinUpSs;
  nBinUpSs<<nBinUp;
  nBinUpSs>>nBinUpStr;
  ofstream outputFile("/sps/atlas/a/aguerguichon/Calibration/Bias/Stats/stats.txt", ios::out);
  if (outputFile == 0) {cout<<"Error while opening stats.txt"<<endl; return 0;}
  outputFile <<"Conf"<<"\tMean"<<"\t RMS"<<"\tMean error\n";

  //Open file, get tree and branches
 
  int  nEff[100][100];
  double xMin[100][100];
  double xMax[100][100];
  
  int k=1;
  while(argv[k]!=NULL)
    {
      TFile *f= TFile::Open(argv[k]);
      if (f == 0) { cout<<"Error: cannot open root file\n"<<endl; return 0;}
      //printf("\n File %s is opened", argv[k]);

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
      nEff[100][100]=0;
      xMin[100][100]=10;
      xMax[100][100]=-10;
    
      //Fill the array of histograms 
      for (int n=0; n<nEntries; n++)
	{
	  t->GetEntry(n);
	  if (nBins!=nBinUp) continue;
	  bias= sigma-inputC;
	  arHisto[iConf][jConf]->Fill(bias);
	  nEff[iConf][jConf]++;
	  if (bias < xMin[iConf][jConf]) xMin[iConf][jConf]=bias;
	  if (bias > xMax[iConf][jConf]) xMax[iConf][jConf]=bias;
	}

      f->Close();
      k++;
    }


  //Draw the histograms
  stringstream iSs[10];
  stringstream jSs[10];
  string iStr[10];
  string jStr[10];
  string nameHisto;
  
  vector <string> vecOptions;//Options to draw the histograms (cf DrawPlot.cxx)
  vecOptions.push_back("yTitle=#Events");
  vecOptions.push_back("xTitle=C^{meas}-C^{input}");

  for (unsigned int i=0; i<nBinUp; i++)
    {
      for (unsigned int j=0; j<=i; j++)
	{  
	  iSs[i]<<i;
	  iSs[i]>>iStr[i];
	  iSs[i].str("");
	  jSs[j]<<j;
	  jSs[j]>>jStr[j];
	  jSs[j].str("");
	  if (arHisto[i][j]->GetEntries() == 0) continue;
	  arHisto[i][j]->GetXaxis()->SetRangeUser(xMin[i][j], xMax[i][j]);
	  nameHisto = "/sps/atlas/a/aguerguichon/Calibration/Bias/Plots/nBins"+nBinUpStr+"/biasConf_" + iStr[i] +jStr[j];
	  DrawPlot( {arHisto[i][j]}, nameHisto,{vecOptions});
	  
	  //Write stats.txt
	  double hMean = arHisto[i][j]->GetMean();
	  double hRMS = arHisto[i][j]->GetRMS();
	  double hMeanErr = hRMS/nEff[i][j];
	  outputFile <<iStr[i]<<jStr[j]<<"\t"<<hMean<<"\t"<<hRMS<<"\t"<< hMeanErr<<"\n";
	}
    }

  outputFile.close();
 
  //End of program
  return 0;
}
