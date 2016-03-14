#include <iostream>
#include <TH1.h>
#include "Calibration_bias/Stats.h"

using namespace std; 

Stats::Stats()
{
}

Stats::~Stats()
{
}

double Stats::GetMeanForConf(TH1D* h, unsigned int stat, double rms, double input)
{
  return h->GetMean();
} 
