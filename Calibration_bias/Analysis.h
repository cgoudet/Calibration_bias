#include <iostream>
#include <string>

#include "TH1.h"

#include <boost/multi_array.hpp>

class Analysis
{

 public:
  Analysis();
  ~Analysis();
  double GetMeanForConf(TH1D* h, unsigned int stat, double rms, double input);
  TH1D* FitGauss(TH1D* h, double xmin, double xmax);
  void WriteFileStats (boost::multi_array <TH1D*, 2> h, boost::multi_array <unsigned int, 2> nEff, boost::multi_array <unsigned int, 2> stat, boost::multi_array <double, 2> input, unsigned int nbin, std::string outFileName);

};
