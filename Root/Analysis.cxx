#include <iostream>
#include <string>
#include <istream>
#include <fstream>

#include <TH1.h>
#include <TF1.h>

#include <boost/multi_array.hpp>
#include "Calibration_bias/Analysis.h"

using namespace std; 



//=======================
Analysis::Analysis()
{
}

Analysis::~Analysis()
{
}


//========================
double Analysis::GetMeanForConf(TH1D* h, unsigned int stat, double rms, double input)
{
  return h->GetMean();
} 

//=========================
TH1D* Analysis::FitGauss(TH1D* h, double xmin, double xmax)
{
  TF1 *f1 = new TF1("f1", "gaus", xmin, xmax);
  h->Fit("f1", "R");
  return h;
}

//=========================

void Analysis::WriteFileStats (boost::multi_array <TH1D*, 2> h, boost::multi_array <unsigned int, 2> nEff, boost::multi_array <unsigned int, 2> stat, boost::multi_array <double, 2> input, unsigned int nbin, string outFileName)
{
  Analysis a;

  ofstream outputFile(outFileName, ios::out);
  if (outputFile == 0) {cout<<"Error while opening stats.txt"<<endl; return;}
  
  outputFile <<"Conf"<<"\tMean"<<"\t RMS"<<"\tMean error\n";
  
  double hRMS, hMean, hMeanErr;

  stringstream iSs[10];
  stringstream jSs[10];
  string iStr[10];
  string jStr[10];
  
  for (unsigned int i=0; i<nbin; i++)
    {
      for (unsigned int j=0; j<=i; j++)
	{
	  if (h[i][j]->GetEntries()==0) continue; 
	  iSs[i]<<i;
	  iSs[i]>>iStr[i];
	  jSs[j]<<j;
	  jSs[j]>>jStr[j];
	 	  
	  hRMS = h[i][j]->GetRMS();
	  hMean = a.GetMeanForConf(h[i][j], stat[i][j], hRMS, input[i][j]);
	  hMeanErr = hRMS/sqrt(nEff[i][j]);
	  outputFile <<iStr[i]<<"_"<<jStr[j]<<"\t"<<hMean<<"\t"<<hRMS<<"\t"<< hMeanErr<<"\n";
	}
    }
  outputFile.close();
  cout<<"\n File "<<outFileName<<" has been written.\n"<<endl;
}
