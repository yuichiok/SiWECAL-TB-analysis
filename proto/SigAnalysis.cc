#include "TROOT.h"
#include "TFile.h"
#include "protoAnalysis.cc"

void SigAnalysis(TString conf="conf1",TString energy="5.8GeV",TString grid="grid20", double mipcut=0.5){


  TString filename = "/home/irles/cernbox/TB2017/ConvertedData/pass1/Tungsten/"+conf+"/"+grid+"/"+energy+"__build.root";
  //eos/project/s/siw-ecal/TB2017-06/DESY/ConvertedData/pass1/Tungsten
  protoAnalysis ss(filename);
  //ss.SimpleDistributionsShower("_"+conf+"_"+energy+"_"+grid);
  ss.ShowerDistributions("results_showers/zmlt1_6slabs/", conf, energy, grid, mipcut);
  
  gSystem->Exit(0);

}
