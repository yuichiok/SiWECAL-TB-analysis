//# Copyright 2020 Adrián Irles IJCLab (CNRS/IN2P3)

#include "TROOT.h"
#include "TFile.h"
#include "DecodedSLBAnalysis.cc"

void SingleSlabAnalysis(TString filename_in, TString output="", int i_slboard=2){


  TString map="../mapping/fev10_chip_channel_x_y_mapping.txt";
  //if(i_slboard==8 || i_slboard==12) map="../mapping/fev11_cob_chip_channel_x_y_mapping.txt"; 
  /*  if(i_slboard==0) map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev11_cob_chip_channel_x_y_mapping.txt";
  if(i_slboard==1) map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev10_chip_channel_x_y_mapping.txt"; 
  if(i_slboard==2) map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev11_cob_chip_channel_x_y_mapping.txt";
  */
  gStyle->SetOptFit(1);

  //  TString treename_st="siwecaldecoded";
  
  filename_in=filename_in+".root";
  
  cout<<" Analyze file: "<<filename_in<<endl;
  DecodedSLBAnalysis ss(filename_in);

  ss.slboard_array_mapping[0]=1;//the first ss.slboard found by the core module has Add. 1
  ss.slboard_array_mapping[1]=2;
  ss.slboard_array_mapping[2]=3;
  ss.slboard_array_mapping[3]=4;
  ss.slboard_array_mapping[4]=5;
  ss.slboard_array_mapping[5]=6;
  ss.slboard_array_mapping[6]=7;
  ss.slboard_array_mapping[7]=10;
  ss.slboard_array_mapping[8]=11;
  ss.slboard_array_mapping[9]=13;
  ss.slboard_array_mapping[10]=14;

  ss.ReadMap(map,i_slboard);
  ss.PedestalAnalysis(i_slboard,output,4,6); //analysis filtering out events with more than 4 triggers in a chip and less than 6 slabs with hit at the same time
  ss.SignalAnalysis(i_slboard,output,4,6);//analysis filtering out events with more than 4 triggers in a chip and less than 6 slabs with hit at the same time
  ss.Retriggers(i_slboard,output,4);

  gSystem->Exit(0);

}
