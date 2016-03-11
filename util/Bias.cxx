#include <iostream>
#include <istream>
#include <fstream>
#include <map>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include <string>
#include <sstream>

using namespace std;

int main()
{

  //Open file, get tree and branches
  TFile *f = TFile::Open ("/sps/atlas/c/cgoudet/Calibration/Bias/Data/TreeToyTemplates_14093222.root");

  if (f == 0) 
    { // print an error message if the file cannot be opened
      cout<<"Error: cannot open root file\n"<<endl;
      return 0;
    }

  //Get the tree
  TTree *t = (TTree*) f->Get("ConfigurationsCTree");

  const int nEntries= t->GetEntries();
  unsigned int iConf, jConf, statTree;
  double sigma, errSigma, inputC, dataRMS;
    
  t->SetBranchAddress("iConf", &iConf);
  t->SetBranchAddress("jConf", &jConf);
  t->SetBranchAddress("statTree", &statTree);
  t->SetBranchAddress("sigma", &sigma);
  t->SetBranchAddress("errSigma", &errSigma);
  t->SetBranchAddress("inputC", &inputC);
  t->SetBranchAddress("dataRMS", &dataRMS);

  /*TBranch *br;
  TIter next( t->GetListOfBranches() );
  string brName='brName';
  double ptr;
  while ( (br = (TBranch *) next()) )
    {
      brName= br->GetName();
      t->SetBranchAddress(brName, &ptr);
      //cout<<ptr<<endl;
      }*/

  double bias;
  int  nEff;
  double xMin;
  double xMax;
  unsigned int nBins=6;
  stringstream iSs[10];
  stringstream jSs[10];
  string iStr[10];
  string jStr[10];
  string nameHisto;
  
  ofstream outputFile("stats.txt", ios::out);
  if (outputFile == 0) {cout<<"Error while opening stats.txt"<<endl; return 0;}
 outputFile <<"Conf"<<"\tMean"<<"\t RMS"<<"\tMean error\n";
  //Create, fill, save the histo for a given configuration
  

  for (unsigned int i=0; i<nBins; i++)
    {
      for (unsigned int j=0; j<=i; j++)
	{
	  TH1D *h= new TH1D("hBias","",100, -0.05, 0.05);
	  xMin=10;
	  xMax=-10;
	  bias=0.;
	  nEff=0;
	  for (int n=0; n<nEntries; n++) 
	    {
	      t->GetEntry(n);
	      if(iConf==i && jConf==j)
		{
		  nEff++;
		  bias= sigma-inputC;
		  if (bias < xMin) xMin=bias;
		  if (bias > xMax) xMax=bias;
		  h->Fill(bias);
		  //cout << iConf << "\t jConf"<< jConf<<endl;
		}
	    }
	  
	  iSs[i]<<i;
	  iSs[i]>>iStr[i];
	  iSs[i].str("");
	  jSs[j]<<j;
	  jSs[j]>>jStr[j];
	  jSs[j].str("");
	  if (h->GetEntries() == 0) continue;
	  TCanvas *c= new TCanvas();
	  h->SetMarkerStyle(2);
	  h->Draw("P");
	  h->GetXaxis()->SetRangeUser(xMin, xMax);
	  h->Draw("P");
	  nameHisto = "~/public/Calibration/biasConf_" + iStr[i] +jStr[j] +".pdf"; 
	  c->Print(nameHisto.c_str());

	  double hMean = h->GetMean();
	  double hRMS = h->GetRMS();
	  double hMeanErr = hRMS/nEff;
	  outputFile <<iStr[i]<<jStr[j]<<"\t"<<hMean<<"\t"<<hRMS<<"\t"<< hMeanErr<<"\n";
	}
    }
  
 
  outputFile.close();

  //End of program
  return 0;
}
