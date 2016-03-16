#include <iostream>
#include "PlotFunctions/DrawPlot"
#include <string>
#include <vector>


Class Bias
{

public:
  Bias (unsigned int NBins);
  double GetBias (double input, double measure, vector <double> optVariables);
  void MakeHist (std::string ToyRootFile);
  void SaveHist (std::string path, std::string nameHist);

private:
  unsigned int nBinUp;
  TH1D* arHisto[100][100];
  unsigned int nEff[100][100];
};
