#include "TROOT.h"
#include "TFile.h"
#include "singleSlabAnalysis.cc"

void PedAnalysis(TString filename_in, TString filename_out,TString out2) {
 
  singleSlabAnalysis ss(filename_in); 
  ss.PedestalAnalysis(filename_out,out2);

  gSystem->Exit(0);

}
