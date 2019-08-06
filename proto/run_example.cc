#include "TROOT.h"
#include "TFile.h"
#include "example.cc"

void run_example(){

  TString run="run_32015";
  TString filename = "/home/irles/WorkAreaECAL/2019/TB201906/EvBuilt/"+run+"_build.root";
  example ss(filename);
  ss.SimpleDistributionsTrack("_"+run);
  ss.SimpleEvDisplayTrack("_"+run);
  gSystem->Exit(0);

}
