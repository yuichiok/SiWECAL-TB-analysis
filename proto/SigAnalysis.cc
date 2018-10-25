#include "TROOT.h"
#include "TFile.h"
#include "protoAnalysis.cc"

void SigAnalysis(TString energy=""){


  TString filename = "/home/irles/WorkAreaECAL/2018/TB201809/Common_stuff/Common/ECAL//ECAL_"+energy+"_pions_plus_build.root";
  //  filename = "/home/irles/WorkAreaECAL/2018/TB201809/Common_stuff/Common/ECAL/ECAL_200GeV_muons_build.root";
  protoAnalysis ss(filename);
  //ss.SimpleDistributionsShower("_"+conf+"_"+energy+"_"+grid);
  double mipcut=0.5;
  ss.ShowerDistributions("results_showers/pions/", energy, mipcut);
  
  gSystem->Exit(0);

}
