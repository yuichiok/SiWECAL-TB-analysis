//# Copyright 2020  Adrián Irles IJCLab (CNRS/IN2P3)

#include "TROOT.h"
#include "TFile.h"
#include "DecodedSLBAnalysis.cc"

void Monitoring(TString filename_in, TString output="", int freq=1, bool shifter=false){


  filename_in=filename_in+".root";
  cout<<" Monitoring of file: "<<filename_in<<endl;
  DecodedSLBAnalysis ss(filename_in);

  TString map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev10_chip_channel_x_y_mapping.txt";

  for(int i_layer=0; i_layer<15; i_layer++) {
    if(i_layer==8 || i_layer==12) map="../mapping/fev11_cob_chip_channel_x_y_mapping.txt";
    else map="../mapping/fev10_chip_channel_x_y_mapping.txt"; 
    
    ss.ReadMap(map,i_layer);
  }

  //  ss.n_layers=3;
  ss.Monitoring(output,freq,shifter);
  //ss.SynchronizationStudies(output,freq,shifter); 
  if(shifter==false) gSystem->Exit(0);

}
