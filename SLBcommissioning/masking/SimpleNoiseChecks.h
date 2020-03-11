//# Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)

#include "DecodedSLBAnalysis.cc"
#include "TROOT.h"
#include "TFile.h"
//#include "../conf_struct.h"
#include "../scurves/fithistos.C"

std::vector<std::vector<std::vector<int> > > mask;

void init() {

  for(int i=0; i<15; i++) {
    std::vector<std::vector<int> >temp1;

    for(int j=0; j<16; j++) {
      std::vector<int> temp2;

      for(int k=0; k<64; k++) {
	temp2.push_back(0);
      }
      temp1.push_back(temp2);
    }
    mask.push_back(temp1);
  }
  
}

TString run_settings="../Run_Settings_";
TString scurves="RateVsThresholdScan_02192020_SLBoard_test2";

bool trig_check(std::array<int,8> noiselevels, int threshold) {
	if(noiselevels[4]>threshold*noiselevels[3]) return true;
	else return false;
}

bool adc_check(std::array<int,8> noiselevels, int threshold) {
	if(noiselevels[7]>threshold*noiselevels[3]) return true;
	else return false;
}

bool retrig_check(std::array<int,8> noiselevels, int threshold) {
	if(noiselevels[5]>threshold*noiselevels[3]) return true;
	else return false;
}

bool retrigall_check(std::array<int,8> noiselevels, int threshold) {
  if((noiselevels[5]+noiselevels[6])>threshold*noiselevels[3]) return true;
	else return false;
}

bool all_check(std::array<int,8> noiselevels, int threshold) {
  if(( noiselevels[4]+noiselevels[5]+noiselevels[6])>threshold*noiselevels[3]) return true;
	else return false;
}

void apply_mask() {

  for(int i=0; i<15; i++) {
    for(int j=0; j<16; j++) {
      for(int k=0; k<64; k++) {
	if(detector.slab[0][i].asu[0].skiroc[j].mask[k]==0 && mask.at(i).at(j).at(k)==1) {
	  detector.slab[0][i].asu[0].skiroc[j].mask[k]=1;
	  cout<<" Masking:  slab:"<<i<<" skiroc:"<<j<<" chn:"<<k<<endl;
	}
      }
    }
  }
  
}

void double_check(TString filename_in, int iteration=1, int window=1, float adc=1, float retrig=5, float trig=5){

  TString settings=run_settings+TString::Format("it%i.txt",iteration);
  
  read_configuration_file(settings,false);
  cout<<"Reading configuration file: "<<settings<<endl;

  TString filename=filename_in+"_first_0.root";
  DecodedSLBAnalysis ss_0(filename);
  std::vector<std::array<int,8>> noiselevels_0= ss_0.NoiseLevels(window);
  for(unsigned i=0; i<noiselevels_0.size(); i++) {
    if( adc_check(noiselevels_0.at(i),adc) == true ||  trig_check(noiselevels_0.at(i),trig) == true ||  retrig_check(noiselevels_0.at(i),retrig) == true)
      mask.at(noiselevels_0.at(i)[0]).at(noiselevels_0.at(i)[1]).at(noiselevels_0.at(i)[2])=1;
  }
  
  
  filename=filename_in+"_first_1.root";
  DecodedSLBAnalysis ss_1(filename);
  std::vector<std::array<int,8>> noiselevels_1= ss_1.NoiseLevels(window);
  for(unsigned i=0; i<noiselevels_1.size(); i++) {
    if( mask.at(noiselevels_0.at(i)[0]).at(noiselevels_0.at(i)[1]).at(noiselevels_0.at(i)[2])==1 && (adc_check(noiselevels_1.at(i),adc) == true ||  trig_check(noiselevels_1.at(i),trig) == true ||  retrig_check(noiselevels_1.at(i),retrig) == true) )
      mask.at(noiselevels_0.at(i)[0]).at(noiselevels_0.at(i)[1]).at(noiselevels_0.at(i)[2])=1;
    else mask.at(noiselevels_0.at(i)[0]).at(noiselevels_0.at(i)[1]).at(noiselevels_0.at(i)[2])=0;
  }

  settings=run_settings+TString::Format("it%i.txt",iteration+1);
  apply_mask();
  write_configuration_file(settings);
  cout<<"Writing configuration file: "<<settings<<endl;

}

void single_check(TString filename_in, int iteration=1, int window=1, float adc=1, float retrig=5, float trig=5){
  TString settings=run_settings+TString::Format("it%i.txt",iteration);
  read_configuration_file(settings,false);
  cout<<"Reading configuration file: "<<settings<<endl;

  TString filename=filename_in+TString::Format("_second_%i.root",iteration);
  DecodedSLBAnalysis ss_0(filename);
  std::vector<std::array<int,8>> noiselevels_0= ss_0.NoiseLevels(window);
  for(unsigned i=0; i<noiselevels_0.size(); i++) {
    if( adc_check(noiselevels_0.at(i),adc) == true ||  trig_check(noiselevels_0.at(i),trig) == true ||  retrig_check(noiselevels_0.at(i),retrig) == true)
      mask.at(noiselevels_0.at(i)[0]).at(noiselevels_0.at(i)[1]).at(noiselevels_0.at(i)[2])=1;
  }
  
  settings=run_settings+TString::Format("it%i.txt",iteration+1);
  apply_mask();
  write_configuration_file(settings);
  cout<<"Writing configuration file: "<<settings<<endl;

}


void cosmics_check(TString filename_in, int iteration=1, int window=150, float trig=20, bool thresholds=false){

  if(thresholds==true) fithistos(scurves,6,run_settings,iteration);
  
  TString settings=run_settings+TString::Format("thresholds_it%i.txt",iteration);
  read_configuration_file(settings,false);
  cout<<"Reading configuration file: "<<settings<<endl;

  TString filename=filename_in+TString::Format("_cosmic_%i.root",iteration);
  DecodedSLBAnalysis ss_0(filename);
  std::vector<std::array<int,8>> noiselevels_0= ss_0.NoiseLevels(window);
  for(unsigned i=0; i<noiselevels_0.size(); i++) {
    if(    trig_check(noiselevels_0.at(i),trig) == true ) {
      //if ind threshold != 0, then decrease it

      if(detector.slab[0][noiselevels_0.at(i)[0]].asu[0].skiroc[noiselevels_0.at(i)[1]].chn_threshold_adj[noiselevels_0.at(i)[2]]!=0) {
	detector.slab[0][noiselevels_0.at(i)[0]].asu[0].skiroc[noiselevels_0.at(i)[1]].chn_threshold_adj[noiselevels_0.at(i)[2]]=TMath::Max(0,detector.slab[0][noiselevels_0.at(i)[0]].asu[0].skiroc[noiselevels_0.at(i)[1]].chn_threshold_adj[noiselevels_0.at(i)[2]]-5);
      }

      //if ind threshold =0 , then mask the channel
      if(detector.slab[0][noiselevels_0.at(i)[0]].asu[0].skiroc[noiselevels_0.at(i)[1]].chn_threshold_adj[noiselevels_0.at(i)[2]]==0)
	detector.slab[0][noiselevels_0.at(i)[0]].asu[0].skiroc[noiselevels_0.at(i)[1]].mask[noiselevels_0.at(i)[2]]=1;
	    
		
    }
  }

  settings=run_settings+TString::Format("thresholds_it%i.txt",iteration+1);
  write_configuration_file(settings);
  cout<<"Writing configuration file: "<<settings<<endl;
  
}

