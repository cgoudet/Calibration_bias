#include <iostream>
#include <istream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <boost/program_options.hpp>
#include <boost/multi_array.hpp>

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"

#include "Calibration_bias/BiasAnalysis.h"

using namespace std;
namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  
  //=================
  //Define options to read name of files from the command line
  po::options_description desc("Input data files (.root)");
  vector<string> dataFiles;
  
  desc.add_options()
    ("help", "Display this help message")
    ("dataFiles", po::value<vector <string> >(&dataFiles), "Absolute path" )
    ;
                                   
  po::positional_options_description p;
  p.add("dataFiles", -1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).style(po::command_line_style::unix_style ^ po::command_line_style::allow_short).run(), vm);
  po::notify(vm);
  
  if (vm.count("help")) {cout << desc; return 0;}


  //================

  BiasAnalysis BA("Calibration_bias/ConfigFile/Bias.boost");


  BA.SelectVariables(dataFiles);
  BA.MeasureBias("/sps/atlas/a/aguerguichon/Calibration/Bias/Stats/Stats_2.csv");


  vector <string> vectOptDraw;//Options to draw the histograms (cf DrawPlot.cxx)
  vectOptDraw.push_back("yTitle=#Events");
  vectOptDraw.push_back("xTitle=C^{meas}-C^{input}");

  BA.MakePlots("latex_2.tex", vectOptDraw);

  //================
  //End of program
  cout <<"End of programm"<<endl;
  return 0;
}