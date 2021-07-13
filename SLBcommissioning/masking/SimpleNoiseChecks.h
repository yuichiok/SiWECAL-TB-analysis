//# Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)

#include "DecodedSLBAnalysis.cc"
#include "TROOT.h"
#include "TFile.h"
//#include "../conf_struct.h"
#include "../scurves/fithistos.C"

std::vector<std::vector<std::vector<int> > > mask;
std::vector<std::vector<std::vector<int> > > mask1;
std::vector<std::vector<std::vector<int> > > mask2;
std::vector<std::vector<std::vector<int> > > mask3;

TString run_settings;
TString scurves;
bool debug;
int acqwindow=1;
int acqwindow_cosm=100;

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
    mask1.push_back(temp1);
    mask2.push_back(temp1);
    mask3.push_back(temp1);
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
  cout<<"Noisy channel candidate from SLboard:"<<detector.slab[0][noiselevels[0]].add<<" skiroc:"<<noiselevels[1]<<" chn:"<<noiselevels[2]<<
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

void apply_mask(bool global_threshold=false) {

  for(int i=0; i<15; i++) {
    for(int j=0; j<16; j++) {
      int channels_to_mask=0;
      int channels_notmasked=0;
      int channels_masked=0;
      for(int k=0; k<64; k++) {
	if(detector.slab[0][i].asu[0].skiroc[j].mask[k]==0 && mask.at(i).at(j).at(k)==1) channels_to_mask++;
	if(detector.slab[0][i].asu[0].skiroc[j].mask[k]==0 && mask.at(i).at(j).at(k)==0) channels_notmasked++;
	if(detector.slab[0][i].asu[0].skiroc[j].mask[k]==1) channels_masked++;
      }

      if(global_threshold==true && detector.slab[0][i].asu[0].skiroc[j].threshold_dac<270 && float(channels_to_mask+channels_masked)/float(channels_notmasked+channels_to_mask+channels_masked)>0.2) {
	cout<<" WARNING!! TOO MANY CHANNELS TO BE MASKED ("<<channels_to_mask <<") + already masked: "<<channels_masked <<" in  Slboard:"<<detector.slab[0][i].add<<" skiroc:"<<j<<endl;
      	cout<<" instead we increase the global thershold, from"<<detector.slab[0][i].asu[0].skiroc[j].threshold_dac<<" to: "<<detector.slab[0][i].asu[0].skiroc[j].threshold_dac+5<<endl;
	detector.slab[0][i].asu[0].skiroc[j].threshold_dac+=5;
      } else {
	for(int k=0; k<64; k++) {
	  if(detector.slab[0][i].asu[0].skiroc[j].mask[k]==0 && mask.at(i).at(j).at(k)==1) {
	    detector.slab[0][i].asu[0].skiroc[j].mask[k]=1;
	    cout<<" Masking:  Slboard idx"<<i<< " add:"<< detector.slab[0][i].add <<"  skiroc:"<<j<<" chn:"<<k<<endl;
	  }
	}
      }
      
    }//j
  }//i
  
}

void triple_check(TString filename_in, int iteration=1, int voting=3, int window=1, float underflow_trig=1, float underflow=1, float retrig=1, float trig=1){

  if(debug==true) cout<< "Start Triple Chech with voting "<< voting <<endl;

  
  //if voting =0, no masking
  //if voting =1, masking if is noisy in at least one run (or)
  //if voting =2, masking if is noisy in at least two runs (or)
  //if voting =3, masking if is noisy in the three runs

  filename_in="../../converter_SLB/convertedfiles/Run_ILC_"+filename_in;
  if(debug==true) cout<< "Filename Prefix = "<<filename_in <<endl;
  TString settings=run_settings+TString::Format("it%i.txt",iteration);
  cout<<"Reading configuration file: "<<settings<<endl;
  read_configuration_file(settings,false);

  if(voting!=0) {

    for(int irepetition=0; irepetition<3; irepetition++) {
      TString filename="";
      filename=TString::Format("%s_%i.root",filename_in.Data(),irepetition);
      cout<<"Reading root file: "<<filename<<endl;
      DecodedSLBAnalysis ss(filename);
      std::vector<std::array<int,9>> noiselevels= ss.NoiseLevels(window);
      for(unsigned i=0; i<noiselevels.size(); i++) {
	if( under_or_over_flowcheck_trig(noiselevels.at(i),underflow_trig) == true || under_or_over_flowcheck_all(noiselevels.at(i),underflow) == true ||  trig_check(noiselevels.at(i),trig) == true ||  retrig_check(noiselevels.at(i),retrig) == true) {
	  if(irepetition==0) mask1.at(noiselevels.at(i)[0]).at(noiselevels.at(i)[1]).at(noiselevels.at(i)[2])=1;
	  if(irepetition==1) mask2.at(noiselevels.at(i)[0]).at(noiselevels.at(i)[1]).at(noiselevels.at(i)[2])=1;
	  if(irepetition==2) mask3.at(noiselevels.at(i)[0]).at(noiselevels.at(i)[1]).at(noiselevels.at(i)[2])=1;

	  if(debug==true) cout_check(noiselevels.at(i),underflow_trig,underflow,retrig,trig);
	}
      }
      noiselevels.clear();
    }
        
  }//end if voting!=0
  
  if(voting==1) {
    for(int i=0; i<15; i++) {
      for(int j=0; j<16; j++) {
  	for(int k=0; k<64; k++) {
  	  if( mask1.at(i).at(j).at(k)==1 || mask2.at(i).at(j).at(k)==1 || mask3.at(i).at(j).at(k)==1) mask.at(i).at(j).at(k)=1;
  	}
      }
    }
  }//voting==1
  
  if(voting==2) {
    for(int i=0; i<15; i++) {
      for(int j=0; j<16; j++) {
  	for(int k=0; k<64; k++) {
  	  if( (mask1.at(i).at(j).at(k)==1 && mask2.at(i).at(j).at(k)==1) ||
  	      (mask1.at(i).at(j).at(k)==1 && mask3.at(i).at(j).at(k)==1) ||
  	      (mask2.at(i).at(j).at(k)==1 && mask3.at(i).at(j).at(k)==1) )
  	    mask.at(i).at(j).at(k)=1;
  	}
      }
    }
  }//voting==2
  
  if(voting==3) {
    for(int i=0; i<15; i++) {
      for(int j=0; j<16; j++) {
  	for(int k=0; k<64; k++) {
  	  if( (mask1.at(i).at(j).at(k)==1 && mask2.at(i).at(j).at(k)==1 && mask3.at(i).at(j).at(k)==1))
  	    mask.at(i).at(j).at(k)=1;
  	}
      }
    }
  }//voting==3
  

  settings=run_settings+TString::Format("comm_it%i.txt",iteration+1);
  apply_mask();
  write_configuration_file(settings);
  cout<<"Writing configuration file: "<<settings<<endl;

}

/* void double_check(TString filename_in, int iteration=1,int repetition=0, int window=1, float underflow=1, float retrig=1, float trig=1){ */

/*   filename_in="../../converter_SLB/convertedfiles/Run_ILC_"+filename_in; */
/*   TString settings=run_settings+TString::Format("it%i.txt",iteration); */
  
/*   read_configuration_file(settings,false); */
/*   cout<<"Reading configuration file: "<<settings<<endl; */

/*   TString filename=""; */
/*   if(repetition==1) filename=filename_in+"_first1_0.root"; */
/*   if(repetition==2) filename=filename_in+"_first2_0.root"; */
/*   if(repetition==3) filename=filename_in+"_first3_0.root"; */

/*   cout<<"Reading root file: "<<filename<<endl; */
/*   DecodedSLBAnalysis ss_0(filename); */
/*   std::vector<std::array<int,9>> noiselevels_0= ss_0.NoiseLevels(window); */
/*   for(unsigned i=0; i<noiselevels_0.size(); i++) { */
/*     if( under_or_over_flowcheck_all(noiselevels_0.at(i),underflow) == true ||  trig_check(noiselevels_0.at(i),trig) == true ||  retrig_check(noiselevels_0.at(i),retrig) == true) { */
/*       mask2.at(noiselevels_0.at(i)[0]).at(noiselevels_0.at(i)[1]).at(noiselevels_0.at(i)[2])=1; */
/*       if(debug==true) cout_check(noiselevels_0.at(i),999999,underflow,retrig,trig); */
/*     } */
    
/*   } */
  
  
/*   if(repetition==1) filename=filename_in+"_first1_1.root"; */
/*   if(repetition==2) filename=filename_in+"_first2_1.root"; */
/*   if(repetition==3) filename=filename_in+"_first3_1.root"; */
  
/*   DecodedSLBAnalysis ss_1(filename); */
/*   std::vector<std::array<int,9>> noiselevels_1= ss_1.NoiseLevels(window); */
/*   for(unsigned i=0; i<noiselevels_1.size(); i++) { */
/*     if( mask.at(noiselevels_1.at(i)[0]).at(noiselevels_1.at(i)[1]).at(noiselevels_1.at(i)[2])==1 && (under_or_over_flowcheck_all(noiselevels_1.at(i),underflow) == true ||  trig_check(noiselevels_1.at(i),trig) == true ||  retrig_check(noiselevels_1.at(i),retrig) == true) ){ */
/*       if(debug==true) cout_check(noiselevels_1.at(i),999999,underflow,retrig,trig); */
/*       mask.at(noiselevels_1.at(i)[0]).at(noiselevels_1.at(i)[1]).at(noiselevels_1.at(i)[2])=1; */
/*     }else mask.at(noiselevels_1.at(i)[0]).at(noiselevels_1.at(i)[1]).at(noiselevels_1.at(i)[2])=0; */
/*   } */

/*   settings=run_settings+TString::Format("comm_it%i.txt",iteration+1); */
/*   apply_mask(); */
/*   write_configuration_file(settings); */
/*   cout<<"Writing configuration file: "<<settings<<endl; */

/* } */

/* void single_check(TString filename_in, int iteration=1, int window=1, float underflow=1, float retrig=5, float trig=5){ */

/*   filename_in="../../converter_SLB/convertedfiles/Run_ILC_"+filename_in; */
 
/*   TString settings=run_settings+TString::Format("it%i.txt",iteration); */
/*   read_configuration_file(settings,false); */
/*   cout<<"Reading configuration file: "<<settings<<endl; */

/*   TString filename=filename_in+TString::Format("_second_%i.root",iteration); */
/*   DecodedSLBAnalysis ss_0(filename); */
/*   std::vector<std::array<int,9>> noiselevels_0= ss_0.NoiseLevels(window); */
/*   for(unsigned i=0; i<noiselevels_0.size(); i++) { */
/*     if( under_or_over_flowcheck_trig(noiselevels_0.at(i),underflow) == true ||  trig_check(noiselevels_0.at(i),trig) == true ||  retrig_check(noiselevels_0.at(i),retrig) == true) { */
/*       if(debug==true && mask.at(noiselevels_0.at(i)[0]).at(noiselevels_0.at(i)[1]).at(noiselevels_0.at(i)[2])==0) cout_check(noiselevels_0.at(i),underflow,99999,retrig,trig); */
/*       mask.at(noiselevels_0.at(i)[0]).at(noiselevels_0.at(i)[1]).at(noiselevels_0.at(i)[2])=1; */
/*     } */
/*   } */
  
/*   settings=run_settings+TString::Format("comm_it%i.txt",iteration+1); */
/*   apply_mask(); */
/*   write_configuration_file(settings); */
/*   cout<<"Writing configuration file: "<<settings<<endl; */

/* } */


 void cosmics_check(TString filename_in, int iteration=1, int window=150, float trig=0, bool global=false){

  
  TString settings=run_settings+TString::Format("it%i.txt",iteration);
  read_configuration_file(settings,false);
  cout<<"Reading configuration file: "<<settings<<endl;

  filename_in="../../converter_SLB/convertedfiles/Run_ILC_"+filename_in;
  TString filename=filename_in+TString::Format("_cosmic_it%i.root",iteration);
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
	detector.slab[0][noiselevels_0.at(i)[0]].asu[0].skiroc[noiselevels_0.at(i)[1]].chn_threshold_adj[noiselevels_0.at(i)[2]]-=5;
	cout<<" Adjusting:  slabidx:"<<noiselevels_0.at(i)[0]<<" SLB:"<<detector.slab[0][noiselevels_0.at(i)[0]].add<<" skiroc:"<<noiselevels_0.at(i)[1]<<" chn:"<<noiselevels_0.at(i)[2]<< " - " <<detector.slab[0][noiselevels_0.at(i)[0]].asu[0].skiroc[noiselevels_0.at(i)[1]].chn_threshold_adj[noiselevels_0.at(i)[2]]<<endl;
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
  settings=run_settings+TString::Format("comm_it%i.txt",iteration+1);
  write_configuration_file(settings);
  cout<<"Writing configuration file: "<<settings<<endl;

}

