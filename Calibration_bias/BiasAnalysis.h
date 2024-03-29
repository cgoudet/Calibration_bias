#include <iostream>
#include <istream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include "boost/multi_array.hpp"
#include "TH1.h"
#include "TMatrixD.h"
#include "RooDataSet.h"
#include "RooGaussian.h"

using namespace RooFit;


class BiasAnalysis
{

public:
  BiasAnalysis(std::string configFileName);
  ~BiasAnalysis();
  void SelectVariables(std::vector <std::string> dataFiles);
  void MeasureBias(std::string outFileName, std::string outRootFileName);
  void MakePlots(std::string path, std::string latexFileName, std::string comment=" ");
  void InvertCijMatrix(unsigned int inversionProcedure);

 private:
  std::vector <std::string> m_variablesBias;
  std::vector <unsigned int> m_variablesStats;
  std::vector <std::string> m_varBiasValues;
  std::vector <std::string> m_histNames;
  std::vector <double> m_bias;
  std::vector <double> m_errBias;

  std::map <std::string, TH1D*> m_mapHist;
  std::map <std::string, double> m_mapXMin;
  std::map <std::string, double> m_mapXMax;
  std::map <std::string, unsigned int> m_mapHistPosition;
  std::map <std::string, double> m_mapSumX; 
  std::map <std::string, double> m_mapSumXM;
  std::map <std::string, unsigned int> m_mapNEff;
  std::map <std::string, RooRealVar*> m_mapBias;
  std::map <std::string, RooDataSet*> m_mapDataSet;
  std::map <std::string, RooGaussian*> m_mapGauss;
  std::map <unsigned int, TMatrixD> m_mapCij;
  std::map <unsigned int, TMatrixD> m_mapErrCij; 

  std::string m_inTreeName;

  unsigned int m_nHist;
  unsigned int m_methodStats;
  unsigned int m_checkDistri;
  
  typedef boost::multi_array<double, 2> maDouble;
  maDouble m_histStats;

};
