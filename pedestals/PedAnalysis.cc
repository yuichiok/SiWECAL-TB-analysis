#include "TROOT.h"
#include "TFile.h"
#include "savePedestal.cc"

void PedAnalysis(TString filename_in, TString filename_out,TString filename_grid) {
 
  cout<<filename_in <<" " << filename_out << " " <<filename_grid << endl;
  savePedestal ss(filename_in); 
  //ss.FindMasked(filename_out);
  ss.PedestalAnalysis(filename_out,filename_grid);
  
  gSystem->Exit(0);

}
