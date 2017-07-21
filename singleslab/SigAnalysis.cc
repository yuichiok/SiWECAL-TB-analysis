#include "TROOT.h"
#include "TFile.h"
#include "singleSlabAnalysis.cc"

void SigAnalysis(){//
  TString filename_in = "grid00_dif_1_1_2.raw.root";
  //, TString filename_out,TString filename_grid) {

  
  singleSlabAnalysis ss(filename_in); 
  //if(filename_grid=="") ss.FindMasked(filename_out);
  ss.SignalAnalysis();
  //ss.BcidCorrelations("dif_1_1_2");
  
  //  gSystem->Exit(0);

}
