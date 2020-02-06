//# Copyright 2012-2018  Roman Poeschl, Adrián Irles
//# This file is part of Calicoes.  

#include "TROOT.h"
#include "TFile.h"
#include "singleSlabAnalysis.cc"

void CommissioningAnalysis(TString filename_in, TString output="", TString slboard="_SLB_2"){

  TString map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev10_chip_channel_x_y_mapping.txt";  
  if(slboard=="_SLB_0") map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev11_cob_chip_channel_x_y_mapping.txt";
  if(slboard=="_SLB_1") map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev11_cob_chip_channel_x_y_mapping.txt";
  if(slboard=="_SLB_2") map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev11_cob_chip_channel_x_y_mapping.txt";
  if(slboard=="_SLB_3") map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev11_cob_chip_channel_x_y_mapping.txt";

  gStyle->SetOptFit(1);

  TString treename_st="slboard";
  if(slboard== "_SLB_0" || slboard== "_SLB_1" || slboard== "_SLB_2" || slboard== "_SLB_3" || slboard =="_SLB_USB"){
    filename_in=filename_in+slboard+".root";
    treename_st="slboard";
    slboard=slboard+"_";
  } else {
    filename_in=filename_in+slboard+".raw.root";
    treename_st="fev10";
  }
  cout<<filename_in<<" "<<treename_st<<endl;

  singleSlabAnalysis ss(filename_in,treename_st);
  ss.PedestalAnalysis(slboard,output,map,4);
  ss.SignalAnalysis(slboard,output,map,4);
  ss.Retriggers(slboard,output,map,10);
  //ss.PedestalAnalysis(slboard,output,map,4);
  //ss.SignalAnalysis(slboard,output,map,4);
  gSystem->Exit(0);

}
