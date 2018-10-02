#include "TROOT.h"
#include "TFile.h"
#include "singleSlabAnalysis.cc"

void SigAnalysis(TString filename_in="/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh15/", TString dif = "dif_1_1_1", TString outputname=""){
  //, TString filename_out,TString filename_grid) {

  singleSlabAnalysis ss(filename_in);
  ss.PedestalAnalysis(dif);
  ss.SignalAnalysis(dif,outputname,true);
  ss.Retriggers(dif,outputname);
  gSystem->Exit(0);

}
