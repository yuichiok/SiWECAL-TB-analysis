#include "TROOT.h"
#include "TFile.h"
#include "protoAnalysis.cc"

void SigAnalysis(TString energy=""){


  TString filename = "/eos/project/s/siw-ecal/TB2018-09/built_files/ECAL_"+energy+"_pions_plus_built.root";
  protoAnalysis ss(filename);
  //ss.SimpleDistributionsShower("_"+conf+"_"+energy+"_"+grid);
  double mipcut=0.5;
  ss.ShowerDistributions("results_showers/", energy, mipcut);
  
  gSystem->Exit(0);

}
