#include <iostream>
#include <istream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include "boost/multi_array.hpp"
#include "TH1.h"

class BiasAnalysis
{

public:
  BiasAnalysis(std::string configFileName);
  ~BiasAnalysis();
  void SelectVariables(std::vector <std::string> dataFiles);
  void MeasureBias(std::string outFileName);
  void MakePlots(std::string latexFileName);


 private:
  std::vector <std::string> m_variablesBias;
  std::vector <unsigned int> m_variablesStats;

  std::map <std::string, TH1D*> m_mapHist;
  std::map <std::string, unsigned int> m_mapHistPosition;
  std::map <std::string, double> m_mapSumX; 
  std::map <std::string, double> m_mapSumXM;
  std::map <std::string, unsigned int> m_mapNEff;

  std::string m_inTreeName;

  
  typedef boost::multi_array<double, 2> maDouble;
  maDouble m_histStats;
  maDouble m_histMinMax;

  unsigned int m_nHist;
  unsigned int m_methodStats;
};
