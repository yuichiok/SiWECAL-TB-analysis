//# Copyright 2020  Adrián Irles IJCLab (CNRS/IN2P3)

#include "TROOT.h"
#include "TFile.h"
#include "DecodedSLBAnalysis.cc"

int Proto(TString filename_in, TString output="", TString type="pedestal",TString pedestal_file=""){


  filename_in=filename_in+".root";
  cout<<" Display of file: "<<filename_in<<endl;
  DecodedSLBAnalysis ss(filename_in);

  TString map="../mapping/fev10_chip_channel_x_y_mapping.txt";


  // layer 0 is the slab with slboard 17, etc... setup July 2021
  ss.slboard_array_mapping[0]=18;
  ss.slboard_array_mapping[1]=23;
  ss.slboard_array_mapping[2]=17;
  ss.slboard_array_mapping[3]=22;
  ss.slboard_array_mapping[4]=25;
  ss.slboard_array_mapping[5]=24;
  ss.slboard_array_mapping[6]=31;
  ss.slboard_array_mapping[7]=30;
  ss.slboard_array_mapping[8]=21;
  ss.slboard_array_mapping[9]=20;
  ss.slboard_array_mapping[11]=19;
  ss.slboard_array_mapping[12]=15;
  ss.slboard_array_mapping[13]=14;
  ss.slboard_array_mapping[14]=13;


  for(int i_slboard=0; i_slboard<15; i_slboard++) {
    cout<<filename_in<<" "<<output<<" islboard " <<i_slboard<<endl;
    ss.ReadMap(map,i_slboard);
    if(type=="retriggers") ss.Retriggers(i_slboard,output,4);
  }
  
  cout<<"Start NSlabsNAalysis"<<endl;
  if(type=="pedestal") int result=ss.NSlabsAnalysis(output,"pedestal");
  if(type=="mip")   int result=ss.NSlabsAnalysis(output,"mip",5,pedestal_file);

  //for(int i_slboard=0; i_slboard<1; i_slboard++)  ss.Retriggers(i_slboard,output,4);
  gSystem->Exit(0);

  return 0;

}
