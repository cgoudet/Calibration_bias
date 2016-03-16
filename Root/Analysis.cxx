#include <iostream>
#include <TH1.h>
#include "Calibration_bias/Analysis.h"

using namespace std; 



//=======================
Analysis::Analysis()
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
