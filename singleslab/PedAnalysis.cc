#include "TROOT.h"
#include "TFile.h"
#include "singleSlabAnalysis.cc"

void PedAnalysis(TString filename_in, TString filename_out,TString filename_grid) {
 
  singleSlabAnalysis ss(filename_in);
  ss.PedestalAnalysis(filename_out,filename_grid);
  //if(filename_grid=="") {
    //ss.FindMasked(filename_out);
  //  ss.PedestalAnalysis_bcid(filename_out,filename_grid);
    // ss.BcidCorrelations(filename_out);
  //} else {
    //ss.PedestalAnalysis_bcid(filename_out,filename_grid);
  // ss.PedestalAnalysis_gridpoints(filename_out,filename_grid);
  //}
  //TString pedestal="results_pedestal/gridpoints/Pedestal_"+filename_out+".txt";
  //ss.ReadPedestals(pedestal);
  gSystem->Exit(0);

}
