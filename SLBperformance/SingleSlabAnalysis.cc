//# Copyright 2020 Adrián Irles IJCLab (CNRS/IN2P3)

#include "TROOT.h"
#include "TFile.h"
#include "DecodedSLBAnalysis.cc"

void SingleSlabAnalysis(TString filename_in, TString output="", int i_slboard=2){


  TString map="../mapping/fev10_chip_channel_x_y_mapping.txt";
  if(i_slboard==8 || i_slboard==12) map="../mapping/fev11_cob_chip_channel_x_y_mapping.txt"; 
  /*  if(i_slboard==0) map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev11_cob_chip_channel_x_y_mapping.txt";
  if(i_slboard==1) map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev10_chip_channel_x_y_mapping.txt"; 
  if(i_slboard==2) map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev11_cob_chip_channel_x_y_mapping.txt";
  */
  gStyle->SetOptFit(1);

  //  TString treename_st="siwecaldecoded";
  
  filename_in=filename_in+".root";
  
  cout<<" Analyze file: "<<filename_in<<endl;
  DecodedSLBAnalysis ss(filename_in);
  ss.ReadMap(map,i_slboard);
  ss.PedestalAnalysis(i_slboard,output,4);
  ss.SignalAnalysis(i_slboard,output,4);
  ss.Retriggers(i_slboard,output,10);

  gSystem->Exit(0);

}
