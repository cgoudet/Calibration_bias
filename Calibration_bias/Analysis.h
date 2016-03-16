#include <iostream>
#include "TH1.h"


Class Analysis
{
public:
  Analysis();
  ~Analysis();
  double GetMeanForConf(TH1D* h, unsigned int stat, double rms, double input);
  TH1D* FitGauss(TH1D* h, double xmin, double xmax);
};
