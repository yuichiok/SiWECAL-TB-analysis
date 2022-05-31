//# Copyright 2020  Adrián Irles IJCLab (CNRS/IN2P3)

#include "TROOT.h"
#include "TFile.h"
#include "DecodedSLBAnalysis.cc"

int Proto(TString filename_in, TString output="", int monitoring=0){


  filename_in=filename_in+".root";
  cout<<" Display of file: "<<filename_in<<endl;
  DecodedSLBAnalysis ss(filename_in,"siwecaldecoded");

  cout<<"Start NSlabsNAalysis"<<endl;
  int result1=1;
  int result2=1;
  if(monitoring>0) ss.HitMonitoring(output,3,5);
  else result2=ss.NSlabsAnalysis(output,2,9);
  return result1*result2;

}
