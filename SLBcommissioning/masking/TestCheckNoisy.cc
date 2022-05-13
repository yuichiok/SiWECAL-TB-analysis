//# Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)

#include "SimpleNoiseChecks.h"

void cosmics_check_test(TString filename_in, TString settings, TString settings_new, int window=150, float trig=0, bool global=false){


  read_configuration_file(settings,false);
  cout<<"Reading configuration file: "<<settings<<endl;
  window = detector.acq_window;
  int delay = detector.acq_delay;

  TString filename=filename_in;
  
  DecodedSLBAnalysis ss_0(filename);
  std::vector<std::array<int,9>> noiselevels_0= ss_0.NoiseLevels(window,delay,false,false);
  for(unsigned i=0; i<noiselevels_0.size(); i++) {

    /* if( detector.slab[0][noiselevels_0.at(i)[0]].add==10) {
      cout<<" XXXX  slboard: "<<detector.slab[0][noiselevels_0.at(i)[0]].add;
      cout_check(noiselevels_0.at(i),trig*6.,999999,trig*6.,trig,trig*10);
      }*/
    
    if(    trig_check(noiselevels_0.at(i),trig*5) == true
           || retrig_check(noiselevels_0.at(i),trig*5) == true
           || retrigall_check(noiselevels_0.at(i),trig*10) == true
           || under_or_over_flowcheck_trig(noiselevels_0.at(i),trig*5) == true
	   || all_check(noiselevels_0.at(i),trig*10) == true
           ) {
      //if ind threshold != 0, then decrease it                                                                                                                                         
      if(debug==true && mask.at(noiselevels_0.at(i)[0]).at(noiselevels_0.at(i)[1]).at(noiselevels_0.at(i)[2])==0) cout_check(noiselevels_0.at(i),trig*2.,999999,trig*2.,trig,trig*2);

      if(detector.slab[0][noiselevels_0.at(i)[0]].asu[0].skiroc[noiselevels_0.at(i)[1]].chn_threshold_adj[noiselevels_0.at(i)[2]]>4) {
        // detector.slab[0][noiselevels_0.at(i)[0]].asu[0].skiroc[noiselevels_0.at(i)[1]].chn_threshold_adj[noiselevels_0.at(i)[2]]=detector.slab[0][noiselevels_0.at(i)[0]].asu[0].skiroc[noiselevels_0.at(i)[1]].chn_threshold_adj[noiselevels_0.at(i)[2]]-5;
        // cout<<" Adjusting:  slboard:"<<detector.slab[0][noiselevels_0.at(i)[0]].add<<" skiroc:"<<noiselevels_0.at(i)[1]<<" chn:"<<noiselevels_0.at(i)[2]<< " - " <<detector.slab[0][noiselevels_0.at(i)[0]].asu[0].skiroc[noiselevels_0.at(i)[1]].chn_threshold_adj[noiselevels_0.at(i)[2]]<<endl;
      } else {
        //if ind threshold =0 , then mask the channel                                                                                                                                   
        if(detector.slab[0][noiselevels_0.at(i)[0]].asu[0].skiroc[noiselevels_0.at(i)[1]].chn_threshold_adj[noiselevels_0.at(i)[2]]==0) {
           mask.at(noiselevels_0.at(i)[0]).at(noiselevels_0.at(i)[1]).at(noiselevels_0.at(i)[2])=1;
	}

      }
    }
  }
  noiselevels_0.clear();
  apply_mask(global);
  write_configuration_file( settings_new);
  cout<<"Writing configuration file: "<< settings_new<<endl;

}


// int TestCheckNoisy(TString filename_in="../../converter_SLB/convertedfiles/Run_ILC_10152020_3_cosmic_it19_Ascii/Run_ILC_10152020_3_cosmic_it19_Ascii.root", TString settingsfile="/mnt/win2/Run_data/Run_ILC_10152020_3_cosmic_it19_Ascii/Run_Settings.txt", TString settings_new="Run_Settings_longcosmicrun_6.txt"){
int TestCheckNoisy(TString filename_in="default", TString settingsfile="default", TString settings_new="default"){

  if(filename_in=="default" || settingsfile=="default" || settings_new=="default"){
    std::cout << "ERROR: files not properly set!" << std::endl;
    return 0;
  }

  //filename_in is the full name of the converted root file
  // setingfile is the full name of the settings file
 
  init(filename_in,true);

  cosmics_check_test(filename_in, settingsfile, settings_new, 100,10,true);
 
  return 0;
}



