//# Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)

#include "DecodedSLBAnalysis.cc"
#include "TROOT.h"
#include "TFile.h"
//#include "../conf_struct.h"
#include "../scurves/fithistos.C"

std::vector<std::vector<std::vector<int> > > mask;

TString run_settings;
TString scurves;
bool debug;

void init(TString s1="", bool debug_=false) {

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

  run_settings="../"+s1+"/Run_Settings_";
  scurves="RateVsThresholdScan_"+s1+"_SLBoard.txt";
  debug=debug_;
}



bool trig_check(std::array<int,9> noiselevels, float threshold) {
	if(noiselevels[4]>threshold*noiselevels[3]) return true;
	else return false;
}

bool under_or_over_flowcheck_trig(std::array<int,9> noiselevels, float threshold) {
	if(noiselevels[7]>threshold*noiselevels[3]) return true;
	else return false;
}

bool under_or_over_flowcheck_all(std::array<int,9> noiselevels, float threshold) {
        if(noiselevels[8]>threshold*noiselevels[3]) return true;
        else return false;
}


bool retrig_check(std::array<int,9> noiselevels, float threshold) {
	if(noiselevels[5]>threshold*noiselevels[3]) return true;
	else return false;
}

bool retrigall_check(std::array<int,9> noiselevels, float threshold) {
  if((noiselevels[5]+noiselevels[6])>threshold*noiselevels[3]) return true;
	else return false;
}

bool all_check(std::array<int,9> noiselevels, float threshold) {
  if(( noiselevels[4]+noiselevels[5]+noiselevels[6])>threshold*noiselevels[3]) return true;
	else return false;
}

void cout_check(std::array<int,9> noiselevels, float th_underflow=9999999, float th_underflow_all=9999999, float th_retr_first=9999999, float th_trig=9999999, float th_retr_all=9999999) {
  cout<<"Noisy channel candidate from Layer:"<<noiselevels[0]<<" skiroc:"<<noiselevels[1]<<" chn:"<<noiselevels[2]<<
    "  Exp. Cosmics:"<<noiselevels[3]<<" underflows (trig): "<<noiselevels[7]<<"(th:"<<th_underflow<<"xExp.)"<<
    " underflows (all): "<<noiselevels[8]<<"(th:"<<th_underflow<<"xExp.)"<<
    " retriger_start: "<<noiselevels[5]<<"(th:"<<th_retr_first<<"xExp.)"<<
    " retrigger_all: "<<noiselevels[6]<<"(th:"<<th_retr_all<<"xExp.)"<<
    " trigers: "<<noiselevels[4]<<"(th:"<<th_trig<<"xExp.)"<<endl;
  /*noisylevels.at(0)=ilayer;
    noisylevels.at(1)=ichip;
    noisylevels.at(2)=ichn;
    noisylevels.at(3)=int(expected);
    noisylevels.at(4)=trigger[ilayer][ichip][ichn];
    noisylevels.at(5)=retrigger_start[ilayer][ichip][ichn];
    noisylevels.at(6)=retrigger_train[ilayer][ichip][ichn];
    noisylevels.at(7)=under_or_over_flow[ilayer][ichip][ichn];
  */

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

void double_check(TString filename_in, int iteration=1,int repetition=0, int window=1, float underflow=1, float retrig=1, float trig=1){

  filename_in="../../converter_SLB/convertedfiles/Run_ILC_"+filename_in;
  TString settings=run_settings+TString::Format("it%i.txt",iteration);
  
  read_configuration_file(settings,false);
  cout<<"Reading configuration file: "<<settings<<endl;

  TString filename="";
  if(repetition==1) filename=filename_in+"_first1_0.root";
  if(repetition==2) filename=filename_in+"_first2_0.root";
  if(repetition==3) filename=filename_in+"_first3_0.root";

  cout<<"Reading root file: "<<filename<<endl;
  DecodedSLBAnalysis ss_0(filename);
  std::vector<std::array<int,9>> noiselevels_0= ss_0.NoiseLevels(window);
  for(unsigned i=0; i<noiselevels_0.size(); i++) {
    if( under_or_over_flowcheck_all(noiselevels_0.at(i),underflow) == true ||  trig_check(noiselevels_0.at(i),trig) == true ||  retrig_check(noiselevels_0.at(i),retrig) == true) {
      mask.at(noiselevels_0.at(i)[0]).at(noiselevels_0.at(i)[1]).at(noiselevels_0.at(i)[2])=1;
      if(debug==true) cout_check(noiselevels_0.at(i),999999,underflow,retrig,trig);
    }
    
  }
  
  
  if(repetition==1) filename=filename_in+"_first1_1.root";
  if(repetition==2) filename=filename_in+"_first2_1.root";
  if(repetition==3) filename=filename_in+"_first3_1.root";
  
  DecodedSLBAnalysis ss_1(filename);
  std::vector<std::array<int,9>> noiselevels_1= ss_1.NoiseLevels(window);
  for(unsigned i=0; i<noiselevels_1.size(); i++) {
    if( mask.at(noiselevels_0.at(i)[0]).at(noiselevels_0.at(i)[1]).at(noiselevels_0.at(i)[2])==1 && (under_or_over_flowcheck_all(noiselevels_1.at(i),underflow) == true ||  trig_check(noiselevels_1.at(i),trig) == true ||  retrig_check(noiselevels_1.at(i),retrig) == true) ){
      if(debug==true) cout_check(noiselevels_0.at(i),999999,underflow,retrig,trig);
      mask.at(noiselevels_0.at(i)[0]).at(noiselevels_0.at(i)[1]).at(noiselevels_0.at(i)[2])=1;
    }else mask.at(noiselevels_0.at(i)[0]).at(noiselevels_0.at(i)[1]).at(noiselevels_0.at(i)[2])=0;
  }

  settings=run_settings+TString::Format("comm_it%i.txt",iteration+1);
  apply_mask();
  write_configuration_file(settings);
  cout<<"Writing configuration file: "<<settings<<endl;

}

void single_check(TString filename_in, int iteration=1, int window=1, float underflow=1, float retrig=5, float trig=5){

  filename_in="../../converter_SLB/convertedfiles/Run_ILC_"+filename_in;
 
  TString settings=run_settings+TString::Format("it%i.txt",iteration);
  read_configuration_file(settings,false);
  cout<<"Reading configuration file: "<<settings<<endl;

  TString filename=filename_in+TString::Format("_second_%i.root",iteration);
  DecodedSLBAnalysis ss_0(filename);
  std::vector<std::array<int,9>> noiselevels_0= ss_0.NoiseLevels(window);
  for(unsigned i=0; i<noiselevels_0.size(); i++) {
    if( under_or_over_flowcheck_trig(noiselevels_0.at(i),underflow) == true ||  trig_check(noiselevels_0.at(i),trig) == true ||  retrig_check(noiselevels_0.at(i),retrig) == true) {
      if(debug==true && mask.at(noiselevels_0.at(i)[0]).at(noiselevels_0.at(i)[1]).at(noiselevels_0.at(i)[2])==0) cout_check(noiselevels_0.at(i),underflow,99999,retrig,trig);
      mask.at(noiselevels_0.at(i)[0]).at(noiselevels_0.at(i)[1]).at(noiselevels_0.at(i)[2])=1;
    }
  }
  
  settings=run_settings+TString::Format("comm_it%i.txt",iteration+1);
  apply_mask();
  write_configuration_file(settings);
  cout<<"Writing configuration file: "<<settings<<endl;

}


void cosmics_check(TString filename_in, int iteration=1, int window=150, float trig=0){

  
  TString settings=run_settings+TString::Format("it%i.txt",iteration);
  read_configuration_file(settings,false);
  cout<<"Reading configuration file: "<<settings<<endl;

  filename_in="../../converter_SLB/convertedfiles/Run_ILC_"+filename_in;
  TString filename=filename_in+TString::Format("_cosmic_%i.root",iteration);
  DecodedSLBAnalysis ss_0(filename);
  std::vector<std::array<int,9>> noiselevels_0= ss_0.NoiseLevels(window);
  for(unsigned i=0; i<noiselevels_0.size(); i++) {
    if(    trig_check(noiselevels_0.at(i),trig) == true
	   || retrig_check(noiselevels_0.at(i),trig/2.) == true
	   || retrigall_check(noiselevels_0.at(i),trig*5.) == true
	   || under_or_over_flowcheck_trig(noiselevels_0.at(i),trig/5.) == true
	   ) {
      //if ind threshold != 0, then decrease it
      if(debug==true && mask.at(noiselevels_0.at(i)[0]).at(noiselevels_0.at(i)[1]).at(noiselevels_0.at(i)[2])==0) cout_check(noiselevels_0.at(i),trig/5.,999999,trig/2.,trig,trig*5);

      if(detector.slab[0][noiselevels_0.at(i)[0]].asu[0].skiroc[noiselevels_0.at(i)[1]].chn_threshold_adj[noiselevels_0.at(i)[2]]>4) {
	detector.slab[0][noiselevels_0.at(i)[0]].asu[0].skiroc[noiselevels_0.at(i)[1]].chn_threshold_adj[noiselevels_0.at(i)[2]]=detector.slab[0][noiselevels_0.at(i)[0]].asu[0].skiroc[noiselevels_0.at(i)[1]].chn_threshold_adj[noiselevels_0.at(i)[2]]-5;
	cout<<" Adjusting:  slab:"<<noiselevels_0.at(i)[0]<<" skiroc:"<<noiselevels_0.at(i)[1]<<" chn:"<<noiselevels_0.at(i)[2]<< " - " <<detector.slab[0][noiselevels_0.at(i)[0]].asu[0].skiroc[noiselevels_0.at(i)[1]].chn_threshold_adj[noiselevels_0.at(i)[2]]<<endl;
      } else {
	//if ind threshold =0 , then mask the channel
	if(detector.slab[0][noiselevels_0.at(i)[0]].asu[0].skiroc[noiselevels_0.at(i)[1]].chn_threshold_adj[noiselevels_0.at(i)[2]]==0) {
	  detector.slab[0][noiselevels_0.at(i)[0]].asu[0].skiroc[noiselevels_0.at(i)[1]].mask[noiselevels_0.at(i)[2]]=1;
	  cout<<" Masking:  slab:"<<noiselevels_0.at(i)[0]<<" skiroc:"<<noiselevels_0.at(i)[1]<<" chn:"<<noiselevels_0.at(i)[2]<<endl;
	}
		
      }
    }
  }

  settings=run_settings+TString::Format("comm_it%i.txt",iteration+1);
  write_configuration_file(settings);
  cout<<"Writing configuration file: "<<settings<<endl;

}

