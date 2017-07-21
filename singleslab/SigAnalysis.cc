#include "TROOT.h"
#include "TFile.h"
#include "singleSlabAnalysis.cc"

void SigAnalysis(TString filename_in = "dif_1_1_1"){
  //, TString filename_out,TString filename_grid) {

  filename_in="/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh15/"+filename_in+".raw.root"
  singleSlabAnalysis ss(filename_in); 
  ss.SignalAnalysis();
  
  gSystem->Exit(0);

}
