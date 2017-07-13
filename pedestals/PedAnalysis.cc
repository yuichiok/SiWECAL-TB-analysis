#include "TROOT.h"
#include "TFile.h"
#include "savePedestal.cc"

void PedAnalysis(TString filename_in, TString filename_out) {
 
  
  savePedestal ss(filename_in); 
  ss.Loop(filename_out);
  
  gSystem->Exit(0);

}
