//# Copyright 2012-2018  Roman Poeschl, Adrián Irles
//# This file is part of Calicoes.  

#include "TROOT.h"
#include "TFile.h"
#include "DecodedSLBAnalysis.cc"

void Monitoring(TString filename_in, TString output="", int freq=1){



  filename_in=filename_in+".root";
  cout<<" Monitoring of file: "<<filename_in<<endl;
  DecodedSLBAnalysis ss(filename_in);

  TString map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev10_chip_channel_x_y_mapping.txt";

  for(int i_slboard=0; i_slboard<3; i_slboard++) {
    if(i_slboard==0) map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev11_cob_chip_channel_x_y_mapping.txt";
    if(i_slboard==1) map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev10_chip_channel_x_y_mapping.txt"; 
    if(i_slboard==2) map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev11_cob_chip_channel_x_y_mapping.txt";
    ss.ReadMap(map,i_slboard);
  }

  //  ss.n_slboards=3;
  ss.Monitoring(output,freq);
  //  gSystem->Exit(0);

}
