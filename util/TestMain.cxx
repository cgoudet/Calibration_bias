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

  BiasAnalysis BA("Calibration_bias/ConfigFile/BiasConf.boost");

  string path= "/sps/atlas/a/aguerguichon/Calibration/Bias/";
  string latexFileName= "BiasConf_CheckRepro_26455434";
  
  BA.SelectVariables(dataFiles);
  BA.MeasureBias(path+"Stats/BiasConf_CheckRepro.csv", path+"RootFiles/BiasConf_CheckRepro.root");

  path= "/sps/atlas/a/aguerguichon/Calibration/Bias/Plots/";
  
  BA.MakePlots(path, latexFileName, "Check for reproducibility: set2, TreeToyTemplates 26455434");

  string commandLine = "mv ./"+latexFileName+".pdf "+path;
  system ( commandLine.c_str() );
  commandLine = "rm "+latexFileName+"*";
  system ( commandLine.c_str() );
  

  //  BA.InvertCijMatrix(11);


  
  //================
  //End of program
  cout <<"End of programm"<<endl;
  return 0;}

