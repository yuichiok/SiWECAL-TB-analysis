//# Copyright 2020  Adrián Irles IJCLab (CNRS/IN2P3)

#include "TROOT.h"
#include "TFile.h"
#include "DecodedSLBAnalysis.cc"

int Proto(TString filename_in, TString output="", TString type="pedestal",TString pedestal_file=""){


  filename_in=filename_in+".root";
  cout<<" Display of file: "<<filename_in<<endl;
  DecodedSLBAnalysis ss(filename_in);

  TString map="/mnt/HardDrive/cernbox_hd/SiWECAL/TB2021/SiWECAL-TB-analysis/mapping/fev10_chip_channel_x_y_mapping.txt";


  // layer 0 is the slab with slboard 17, etc... setup July 2021
  ss.slboard_array_mapping[0]=17;
  ss.slboard_array_mapping[1]=8;
  ss.slboard_array_mapping[2]=10;
  ss.slboard_array_mapping[3]=5;
  ss.slboard_array_mapping[4]=1;
  ss.slboard_array_mapping[5]=13;
  ss.slboard_array_mapping[6]=11;
  ss.slboard_array_mapping[7]=7;
  ss.slboard_array_mapping[8]=14;
  ss.slboard_array_mapping[9]=3;
  ss.slboard_array_mapping[11]=6;
  ss.slboard_array_mapping[12]=9;
  ss.slboard_array_mapping[13]=2;
  ss.slboard_array_mapping[14]=0;


  for(int i_slboard=0; i_slboard<15; i_slboard++) {
    //  if(i_slboard==8 || i_slboard==12) map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis-git/mapping/fev11_cob_chip_channel_x_y_mapping.txt";
    //  else map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis-git/mapping/fev10_chip_channel_x_y_mapping.txt"; 
    //
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
