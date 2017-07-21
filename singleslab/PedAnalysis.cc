#include "TROOT.h"
#include "TFile.h"
#include "savePedestal.cc"

void PedAnalysis(TString filename_in, TString filename_out,TString filename_grid) {
 
  savePedestal ss(filename_in); 
  //if(filename_grid=="") ss.FindMasked(filename_out);
  ss.PedestalAnalysis(filename_out,filename_grid);
  //ss.BcidCorrelations("dif_1_1_2");
  
  gSystem->Exit(0);

}
