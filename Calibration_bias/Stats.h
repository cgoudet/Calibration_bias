#include <iostream>
#include "TH1.h"

using namespace std;


class Stats
{
public:
  Stats();
  ~Stats();
  double GetMeanForConf(TH1D* h, unsigned int stat, double rms, double input);
};
