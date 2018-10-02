#include "TROOT.h"
#include "TFile.h"
#include "protoAnalysis.cc"

void SigAnalysis(TString energy="80GeV", double mipcut=0.5){


  TString filename = "/home/irles/WorkAreaECAL/2018/TB201809/electrons/80GeV/electrons_80GeV_build.root";
  protoAnalysis ss(filename);
  //ss.SimpleDistributionsShower("_"+conf+"_"+energy+"_"+grid);
  ss.ShowerDistributions("results_showers/", energy, mipcut);
  
  gSystem->Exit(0);

}
