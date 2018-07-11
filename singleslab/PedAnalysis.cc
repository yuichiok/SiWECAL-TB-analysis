#include "TROOT.h"
#include "TFile.h"
#include "singleSlabAnalysis.cc"

void PedAnalysis(TString filename_in, TString filename_out) {
 
  singleSlabAnalysis ss(filename_in); 
  ss.PedestalAnalysis(filename_out);

  gSystem->Exit(0);

}
