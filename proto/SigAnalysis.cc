#include "TROOT.h"
#include "TFile.h"
#include "protoAnalysis.cc"

void SigAnalysis(TString conf="conf1",TString energy="5.8GeV",TString grid="grid20", double mipcut=0.5){


  TString filename = "/home/irles/cernbox/TB2017/TBdata/Tungsten/"+conf+"/"+grid+"/"+energy+"__build.root";
  protoAnalysis ss(filename);
  ss.SimpleDistributionsShower("_"+conf+"_"+energy+"_"+grid);
  
  gSystem->Exit(0);

}
