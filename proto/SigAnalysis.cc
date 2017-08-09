#include "TROOT.h"
#include "TFile.h"
#include "protoAnalysis.cc"

void SigAnalysis(){
  //, TString filename_out,TString filename_grid) {

  protoAnalysis ss("/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh15/fullscan/3GeV_calibration_build.root"); 
  ss.SimpleMIPAnalysis("");
  
  gSystem->Exit(0);

}
