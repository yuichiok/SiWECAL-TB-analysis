<<<<<<< HEAD
//# Copyright 2020 Adrián Irles
=======
//# Copyright 2020 Adrián Irles IJCLab (CNRS/IN2P3)
>>>>>>> abab2458134d850f9ab68580bedac4a57d9401fb

#include "TROOT.h"
#include "TFile.h"
#include "DecodedSLBAnalysis.cc"

void SingleSlabAnalysis(TString filename_in, TString output="", int i_layer=2){


  TString map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev10_chip_channel_x_y_mapping.txt";
  if(i_layer==1) map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev11_cob_chip_channel_x_y_mapping.txt"; 
  /*  if(i_layer==0) map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev11_cob_chip_channel_x_y_mapping.txt";
  if(i_layer==1) map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev10_chip_channel_x_y_mapping.txt"; 
  if(i_layer==2) map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev11_cob_chip_channel_x_y_mapping.txt";
  */
  gStyle->SetOptFit(1);

  //  TString treename_st="siwecaldecoded";
  
  filename_in=filename_in+".root";
  
  cout<<" Analyze file: "<<filename_in<<endl;
  DecodedSLBAnalysis ss(filename_in);
  ss.ReadMap(map,i_layer);
  ss.PedestalAnalysis(i_layer,output,4);
  ss.SignalAnalysis(i_layer,output,4);
  ss.Retriggers(i_layer,output,10);
<<<<<<< HEAD
=======
  
>>>>>>> abab2458134d850f9ab68580bedac4a57d9401fb
  gSystem->Exit(0);

}
