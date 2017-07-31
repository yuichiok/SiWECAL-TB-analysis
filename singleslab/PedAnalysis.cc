#include "TROOT.h"
#include "TFile.h"
#include "singleSlabAnalysis.cc"

void PedAnalysis(TString filename_in, TString filename_out,TString filename_grid) {
 
  singleSlabAnalysis ss(filename_in); 
  if(filename_grid=="") {
    //ss.FindMasked(filename_out);
    //ss.PedestalAnalysis_bcid(filename_out,filename_grid);
    ss.BcidCorrelations(filename_out);
  } else {
    //    ss.PedestalAnalysis_bcid(filename_out,filename_grid);
    ss.PedestalAnalysis_gridpoints(filename_out,filename_grid);
  }
  gSystem->Exit(0);

}
