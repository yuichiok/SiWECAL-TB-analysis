//# Copyright 2020  Adrián Irles IJCLab (CNRS/IN2P3)

#include "TROOT.h"
#include "TFile.h"
#include "DecodedSLBAnalysis.cc"

void SlowControlMonitoring(TString filename_in, TString output=""){


  filename_in=filename_in+".root";
  cout<<" SlowControl Monitoring of file: "<<filename_in<<endl;
  DecodedSLBAnalysis ss(filename_in);

  ss.SlowControlMonitoring(output);
  // gSystem->Exit(0);

}
