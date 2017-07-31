#include "TROOT.h"
#include "TFile.h"
#include "singleSlabAnalysis.cc"

void SigAnalysisMagnet(TString filename_in="/home/irles/cernbox/TB2017/TBdata/Magnet/", TString dif = "dif_1_1_1", TString outputname=""){


  //filename_in=filename_in+dif+".raw.root";
  singleSlabAnalysis ss(filename_in); 
  ss.PedestalAnalysis(dif,outputname);
  ss.SignalAnalysis(dif,outputname,false);
  //root -l SigAnalysisMagnet.cc\(\"$file\",\"$dif\",\"$output_file\"\)
  
  gSystem->Exit(0);

}

