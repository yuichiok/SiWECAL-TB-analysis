#include "TROOT.h"
#include "TFile.h"
#include "singleSlabAnalysis.cc"

void SigAnalysis(TString filename_in="/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh15/", TString dif = "dif_1_1_1", TString outputname=""){
  //, TString filename_out,TString filename_grid) {

  filename_in=filename_in+dif+".raw.root";
  singleSlabAnalysis ss(filename_in); 
  ss.SignalAnalysis(dif,outputname,true);
  
  gSystem->Exit(0);

}
