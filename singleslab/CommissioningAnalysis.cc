#include "TROOT.h"
#include "TFile.h"
#include "singleSlabAnalysis.cc"

void CommissioningAnalysis(TString filename_in, TString output, TString dif){

  filename_in=filename_in+dif+".raw.root";
  singleSlabAnalysis ss(filename_in); 
  ss.PedestalAnalysis(dif,output);
  ss.SignalAnalysis(dif,output,true);
  
  gSystem->Exit(0);

}
