//# Copyright 2020  Adrián Irles IJCLab (CNRS/IN2P3)

#include "TROOT.h"
#include "TFile.h"
#include "DecodedSLBAnalysis.cc"

void Monitoring(TString filename_in, TString output="", int freq=1, bool shifter=false, bool debug1=false, bool debug2=false){


  filename_in=filename_in+".root";
  cout<<" Monitoring of file: "<<filename_in<<endl;
  DecodedSLBAnalysis ss(filename_in);

  TString map;

  for(int i_layer=0; i_layer<15; i_layer++) {
    map="../mapping/fev10_chip_channel_x_y_mapping.txt"; 
    ss.ReadMap(map,i_layer);
  }

  //  ss.n_layers=3;
  if(debug2==true) ss.HitMapsSimpleTracks(output,10);
  if(debug2==true) ss.Monitoring(output,freq,shifter);
  ss.HitMapsSimpleTracks(output,10);
  //ss.SynchronizationStudies(output,freq,shifter); 
  if(shifter==false) gSystem->Exit(0);

}
