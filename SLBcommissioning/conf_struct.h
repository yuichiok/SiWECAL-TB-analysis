//# Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)                                                  
// example about structures
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TH1.h"
#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TString.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

//# ASU: 0 
//## ChipIndex: 0 ChipId: 0  FeedbackCap: 3 ThresholdDAC: 230 HoldDelay: 130 FSPeakTime: 2 GainSelectionThreshold: 255 
//### Ch: 0 TrigMask: 0 ChThreshold: 0 PAMask: 0

struct skiroc_t {
  int idx;
  int id;
  int n_channels=64;
  int feedback_cap;
  int threshold_dac;
  int hold_delay;
  int fs_peak_time;
  int gain_selection_threshold;
  int mask[64];
  int chn_threshold_adj[64];
  int preamplifier_mask[64];
} skiroc;

struct asu_t {
  int idx;
  int n_chips=16;
  skiroc_t skiroc[16];
} asu;

struct slab_t {
  int idx;
  int add;
  int ser_num;
  string slb_fpga;
  int nb_asus;
  asu_t asu[5];
} slab;


struct detector_t {
  string sl_soft_v;
  string date_run_unix;
  string date_run_date;
  string date_run_time;
  string usb_sernum;
  string fpga_ver;
  int n_core_daughters;
  int trigger_type;
  int acq_window_source;
  int acq_window;
  int acq_delay;
  int ext_sig_level;
  int ext_sig_edge;
  string core_daughter_fpga_ver[2];
  int core_daughter_n_slabs[2];
  slab_t slab[2][15];
} detector;


void read_configuration_file(TString filename="Run_Settings.txt", bool debug=true) {

  //===== Daughter: 0 SlabIdx: 0 SlabAdd: 1 SL_Board_SerNum: -1 FPGA_Version: V2.3.7  Nb_Of_Connected_ASUs: 1 =====

    
  std::ifstream reading_file(filename);
  if(!reading_file){
    std::cout<<" dameyo - damedame"<<std::endl;
  }

  TString tmpst;
  //first line
  // == SETTINGS FILE SAVED WITH ILC SL-SOFTWARE VERSION: V2.16  == DATE OF RUN: UnixTime = 1582821640.563 date = 2020.2.27 time = 17h.40m.40s  ===
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> detector.sl_soft_v >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> detector.date_run_unix >> tmpst>> tmpst >> detector.date_run_date >> tmpst >> tmpst >> detector.date_run_time >> tmpst ;
  if(debug) cout<<detector.sl_soft_v<<" "<<detector.date_run_unix<<" "<< detector.date_run_date <<" "<< detector.date_run_time<<endl;
  //second line
  //----------------------
  //== SYSTEM_TYPE: SL_COREMODULE_INTERFACE USB_SerNum: 2.41A FPGA_Version: V2.1.11  NB_Of_Core_Daughters: 1 ==
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> detector.usb_sernum >> tmpst >> detector.fpga_ver >> tmpst >> detector.n_core_daughters;
  if(debug) cout<<" COREMODULE: USB_SerNum " <<detector.usb_sernum<<", FPGA Ver "<<detector.fpga_ver<<", NB_Of_Core_Daughters "<< detector.n_core_daughters<<endl;

  //third line
  //----------------------
  //== TriggerType: 1  ('0' = FORCE_TRIGGER, '1' = SELF_TRIGGER) ACQWindowSource: 0 ('0' = AUTO, '1'= SOFT_CMD, '2' = EXT_SIG) ACQWindow: 150 (ms) DelayBetweenCycle: 50 (ms) ExtSigLevel: 0 ('0' = TTL, '1' = NIM, '-1' = NA) ExtSigEdge: 0 ('0' = RISING_EDGE, '1' = FALLING_EDGE, '-1' = NA) ==
  reading_file >> tmpst >> tmpst >> tmpst >> detector.trigger_type >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> detector.acq_window_source >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> detector.acq_window >> tmpst >> tmpst >> detector.acq_delay >> tmpst >> tmpst >> detector.ext_sig_level >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> detector.ext_sig_edge >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst ;
  if(debug) cout<<" Trigger Type: "<<detector.trigger_type << " AcqWindowSource ('0' = AUTO, '1'= SOFT_CMD, '2' = EXT_SIG): " <<detector.acq_window_source << " AcqWindow: "<<detector.acq_window << "ms DelayBetweenCycle: "<<detector.acq_delay<< "ms ExtSigLevel ('0' = TTL, '1' = NIM, '-1' = NA): " << detector.ext_sig_level << " ExtSigEdge ('0' = RISING_EDGE, '1' = FALLING_EDGE, '-1' = NA): " << detector.ext_sig_edge <<endl;
  
  //----------------------
  //=== CORE_DAUGHTER: 0 FPGA_Version: V1.3.1 Nb_Of_Connected_Slabs: 6 ===
  for(int i=0; i<detector.n_core_daughters; i++) {
    reading_file >> tmpst >> tmpst  >> tmpst  >> tmpst >> detector.core_daughter_fpga_ver[i] >> tmpst  >> detector.core_daughter_n_slabs[i] >> tmpst;
    if(debug) cout<<" CORE_DAUGHTER: " << i <<" FPGA_Version: " << detector.core_daughter_fpga_ver[i] << " Nb_Of_Connected_Slabs: "<<detector.core_daughter_n_slabs[i]<< endl;
  }

  //----------------------
  //===== Daughter: 0 SlabIdx: 0 SlabAdd: 1 SL_Board_SerNum: -1 FPGA_Version: V2.3.7  Nb_Of_Connected_ASUs: 1 =====
  for(int i=0; i<detector.n_core_daughters; i++) {
    for(int j=0; j<detector.core_daughter_n_slabs[i]; j++) {
      reading_file >> tmpst >> tmpst >> tmpst >>  tmpst >> detector.slab[i][j].idx >> tmpst >> detector.slab[i][j].add >> tmpst >> detector.slab[i][j].ser_num >> tmpst >> detector.slab[i][j].slb_fpga >> tmpst >> detector.slab[i][j].nb_asus>>tmpst ;
      if(debug) cout<< " Daughter: " << i << " SlabIdx: "<< detector.slab[i][j].idx<< " SlabAdd: "<<detector.slab[i][j].add << " SL_Board_SerNum: " << detector.slab[i][j].ser_num << " FPGA_Version: "<<detector.slab[i][j].slb_fpga << " Nb_ASUS: "<<detector.slab[i][j].nb_asus<<endl;

      for(int k=0; k<detector.slab[i][j].nb_asus; k++) {
	reading_file >> tmpst >> tmpst >> detector.slab[i][j].asu[k].idx;
	if(debug) cout<< " ASU: "<< detector.slab[i][j].asu[k].idx<<" "<<endl;

	for(int ichip=0; ichip<detector.slab[i][j].asu[k].n_chips; ichip++) {
	  reading_file >> tmpst >> tmpst >> detector.slab[i][j].asu[k].skiroc[ichip].idx >> tmpst >> detector.slab[i][j].asu[k].skiroc[ichip].id >> tmpst >> detector.slab[i][j].asu[k].skiroc[ichip].feedback_cap >>  tmpst >> detector.slab[i][j].asu[k].skiroc[ichip].threshold_dac >> tmpst >> detector.slab[i][j].asu[k].skiroc[ichip].hold_delay >> tmpst >> detector.slab[i][j].asu[k].skiroc[ichip].fs_peak_time >> tmpst >>detector.slab[i][j].asu[k].skiroc[ichip].gain_selection_threshold;
	  if(debug) cout<< " ChipIdx: "<< detector.slab[i][j].asu[k].skiroc[ichip].idx << " ChipId: "<< detector.slab[i][j].asu[k].skiroc[ichip].id<< " FeedbackCap: "<< detector.slab[i][j].asu[k].skiroc[ichip].feedback_cap << "  ThresholdDAC: "<< detector.slab[i][j].asu[k].skiroc[ichip].threshold_dac<< " HoldDelay: "<< detector.slab[i][j].asu[k].skiroc[ichip].hold_delay<< " FSPeakTime: " << detector.slab[i][j].asu[k].skiroc[ichip].fs_peak_time <<" GainSelectionThreshold: "<< detector.slab[i][j].asu[k].skiroc[ichip].gain_selection_threshold<<endl;

	  for(int ichn=0; ichn<detector.slab[i][j].asu[k].skiroc[ichip].n_channels; ichn++) {
	    TString tmpst1, tmpst2, tmpst3, tmpst4, tmpst5, tmpst6;
	    reading_file >> tmpst1 >>tmpst2 >> tmpst3 >> tmpst4 >> detector.slab[i][j].asu[k].skiroc[ichip].mask[ichn] >> tmpst5 >> detector.slab[i][j].asu[k].skiroc[ichip].chn_threshold_adj[ichn] >> tmpst6 >> detector.slab[i][j].asu[k].skiroc[ichip].preamplifier_mask[ichn];
	    if(debug) cout<< " Chn: " << ichn << " TrigMask: "<< detector.slab[i][j].asu[k].skiroc[ichip].mask[ichn] << " ChThreshold: "<< detector.slab[i][j].asu[k].skiroc[ichip].chn_threshold_adj[ichn]<< " PAMask: "<<detector.slab[i][j].asu[k].skiroc[ichip].preamplifier_mask[ichn]<<endl;
	  }
	}
      }
    }
  }

}

void write_configuration_file(TString filename="Run_Settings_2.txt") {

  //===== Daughter: 0 SlabIdx: 0 SlabAdd: 1 SL_Board_SerNum: -1 FPGA_Version: V2.3.7  Nb_Of_Connected_ASUs: 1 =====
ofstream fout(filename,ios::out);

  fout<<"== SETTINGS FILE SAVED WITH ILC SL-SOFTWARE VERSION: "<<detector.sl_soft_v<<"  == DATE OF RUN: UnixTime = "<<detector.date_run_unix<<" date = "<<detector.date_run_date<<" time = "<<detector.date_run_time<<"  ==="<<endl;
  fout<<"== SYSTEM_TYPE: SL_COREMODULE_INTERFACE USB_SerNum: "<<detector.usb_sernum<<" FPGA_Version: "<<detector.fpga_ver<<"  NB_Of_Core_Daughters: "<<detector.n_core_daughters<<" =="<<endl;
  fout<<"== TriggerType: "<<detector.trigger_type<<"  ('0' = FORCE_TRIGGER, '1' = SELF_TRIGGER) ACQWindowSource: "<<detector.acq_window_source <<" ('0' = AUTO, '1'= SOFT_CMD, '2' = EXT_SIG) ACQWindow: "<< detector.acq_window <<" (ms) DelayBetweenCycle: "<< detector.acq_delay <<" (ms) ExtSigLevel: "<< detector.ext_sig_level  <<" ('0' = TTL, '1' = NIM, '-1' = NA) ExtSigEdge: "<< detector.ext_sig_edge <<" ('0' = RISING_EDGE, '1' = FALLING_EDGE, '-1' = NA) =="<<endl;

  //----------------------
  //=== CORE_DAUGHTER: 0 FPGA_Version: V1.3.1 Nb_Of_Connected_Slabs: 6 ===
  for(int i=0; i<detector.n_core_daughters; i++) {
    fout<<"=== CORE_DAUGHTER: "<<i<<" FPGA_Version: "<<detector.core_daughter_fpga_ver[i]<<" Nb_Of_Connected_Slabs: "<<detector.core_daughter_n_slabs[i]<<" ==="<<endl;
    
    for(int j=0; j<detector.core_daughter_n_slabs[i]; j++) {
      fout<<"===== Daughter: "<<j<<" SlabIdx: "<< detector.slab[i][j].idx <<" SlabAdd: "<< detector.slab[i][j].add <<" SL_Board_SerNum: "<<detector.slab[i][j].ser_num<<" FPGA_Version: "<< detector.slab[i][j].slb_fpga <<"  Nb_Of_Connected_ASUs: "<< detector.slab[i][j].nb_asus <<" ====="<<endl;

      for(int k=0; k<detector.slab[i][j].nb_asus; k++) {
	fout<<"# ASU: "<<detector.slab[i][j].asu[k].idx<<endl;
	for(int ichip=0; ichip<detector.slab[i][j].asu[k].n_chips; ichip++) {
	  fout<<"## ChipIndex: "<< detector.slab[i][j].asu[k].skiroc[ichip].idx<<" ChipId: "<< detector.slab[i][j].asu[k].skiroc[ichip].id<<"  FeedbackCap: "<<detector.slab[i][j].asu[k].skiroc[ichip].feedback_cap<<" ThresholdDAC: "<<detector.slab[i][j].asu[k].skiroc[ichip].threshold_dac <<" HoldDelay: "<<  detector.slab[i][j].asu[k].skiroc[ichip].hold_delay <<" FSPeakTime: "<< detector.slab[i][j].asu[k].skiroc[ichip].fs_peak_time <<" GainSelectionThreshold: "<< detector.slab[i][j].asu[k].skiroc[ichip].gain_selection_threshold <<" "<<endl;
	  
	  for(int ichn=0; ichn<detector.slab[i][j].asu[k].skiroc[ichip].n_channels; ichn++) {
	    fout<<"### Ch: "<<ichn<<" TrigMask: "<<detector.slab[i][j].asu[k].skiroc[ichip].mask[ichn]<<" ChThreshold: "<< detector.slab[i][j].asu[k].skiroc[ichip].chn_threshold_adj[ichn] <<" PAMask: "<< detector.slab[i][j].asu[k].skiroc[ichip].preamplifier_mask[ichn] <<" "<<endl;
	  }
	}
      }
    }
  }

}
