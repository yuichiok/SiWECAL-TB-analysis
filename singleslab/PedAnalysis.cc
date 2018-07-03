#include "TROOT.h"
#include "TFile.h"
#include "singleSlabAnalysis.cc"

void PedAnalysis(TString filename_in, TString filename_out,TString filename_grid) {
 
  singleSlabAnalysis ss(filename_in); 
  ss.PedestalAnalysis(filename_out,filename_grid);

  gSystem->Exit(0);

}
