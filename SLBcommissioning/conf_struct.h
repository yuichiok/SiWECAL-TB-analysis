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

Float_t map_pointX[16][64];
Float_t map_pointY[16][64];

void ReadMap(TString filename) 
{

  std::ifstream reading_file(filename);
  if(!reading_file){
    cout<<" dameyo - damedame"<<endl;
  }
  
  for(int i=0; i<16; i++) {
    for(int j=0; j<64; j++) {
      map_pointX[i][j] = -1000.;
      map_pointY[i][j] = -1000.;
    }
  }

  Int_t tmp_chip = 0,tmp_channel = 0;
  Float_t tmp_x0 = 0 ,tmp_y0 = 0 , tmp_x = 0 , tmp_y = 0 ;
  TString tmpst;
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst ;
  while(reading_file){
    reading_file >> tmp_chip >> tmp_x0 >> tmp_y0 >> tmp_channel >> tmp_x >> tmp_y ;
    map_pointX[tmp_chip][tmp_channel] = -tmp_x ;
    map_pointY[tmp_chip][tmp_channel] = -tmp_y ;
  }

}

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
  int default_ch_threshold=0;
  int default_trig_mask=0;
  int default_pa_mask=0;
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
  int external_clock;
  int n_core_daughters;
  int trigger_type;
  int acq_window_source;
  int acq_window;
  int acq_delay;
  int start_acq_delay;
  int ext_sig_level;
  int ext_sig_edge;
  string core_daughter_fpga_ver[2];
  int core_daughter_n_slabs[2];
  slab_t slab[2][15];
  bool isCore;
} detector;

void read_configuration_file(TString filename="Run_Settings.txt", bool debug=true) {

  //===== Daughter: 0 SlabIdx: 0 SlabAdd: 1 SL_Board_SerNum: -1 FPGA_Version: V2.3.7  Nb_Of_Connected_ASUs: 1 =====

    
  std::ifstream reading_file(filename);
  if(!reading_file){
    std::cout<<" dameyo - damedame"<<std::endl;
  }


  //CRP 11/12/20 initialise bool to check whether we have a core card or a direct interface
  detector.isCore=false;
  TString tmpst;
  std::string line;
  //first line
  //CRP 7/12/20 For test purposes and later usage we store the beginning of the line (see usual c++ references for details)
  streampos begin,end;
  begin = reading_file.tellg();
  std::cout << "Begin 1: " << begin << std::endl;
  // == SETTINGS FILE SAVED WITH ILC SL-SOFTWARE VERSION: V2.16  == DATE OF RUN: UnixTime = 1582821640.563 date = 2020.2.27 time = 17h.40m.40s  ===
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> detector.sl_soft_v >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> detector.date_run_unix >> tmpst>> tmpst >> detector.date_run_date >> tmpst >> tmpst >> detector.date_run_time >> tmpst ;
  //CRP 7/12/20 Playing with positions in the file
  //CRP The action should bring us back to the beginning of the line (if successful we can analyse the string in parallel to the existing code) 
  //reading_file.clear();  
  //reading_file.seekg(begin);  
  //std::cout << "tellg() returned "
  //          << std::setw(5) << reading_file.tellg()
  //          << " after " << line << "\n";
  // extracing the one character from current location
  //char c = reading_file.get();
  // printing the character
  //cout << "The character 1: " << c << endl;
  reading_file.clear();  
  reading_file.seekg(begin);  
  getline (reading_file,line);
  std::cout << "Before printing string" << std::endl;
  std::cout << line << std::endl;
  //getline (reading_file,line);
  //std::cout << "CRP Second line: " << line << std::endl; 
  //std::cout << "begin 1: " <<  begin << std::endl;
  std::size_t found = line.find("ILC SL-SOFTWARE VERSION:");    
  //std::cout << "found: " << found << std::endl;
  std::size_t found1 = line.find_first_of(":",found);
  //std::cout << "found1: " << found1 << std::endl;
  std::size_t found2 = line.find_first_not_of(' ',found1+1);
  //std::cout << "found2: " << found2 << std::endl;
  //std::cout << "found2: " << found2 << std::endl;
  //std::cout << "The first non white space is " << line[found2] << std::endl;
  std::size_t found3 = line.find_first_of(' ',found2+1);
  std::string version = line.substr(found2,found3-found2);
  std::cout << "version: " << version << std::endl;
  found = line.find("UnixTime");
  found1 = line.find_first_of("123456789",found);
  found2 = line.find_first_of(' ',found1);
  std::string unixtime = line.substr(found1,found2-found1);;
  std::cout << "unixtime: " << unixtime << std::endl;
 
  std::cout << "After printing string" << std::endl;
  //CRP
  if(debug) cout<<detector.sl_soft_v<<" "<<detector.date_run_unix<<" "<< detector.date_run_date <<" "<< detector.date_run_time<<endl;
  //second line
  //----------------------
  //CRP 7/12/20 For test purposes and later usage we store the beginning of the line (see usual c++ references for details)
  begin = reading_file.tellg();
  //std::cout << "begin 2: " <<  begin << std::endl;
  //reading_file.seekg(begin);  
  // extracing a set of characters from current location
  //char A[6];
  //reading_file.read(A, 5);
  //A[5] = 0;
  // printing the character stream 
  //cout << "The character 2+: " << A << endl;
  //== SYSTEM_TYPE: SL_COREMODULE_INTERFACE USB_SerNum: 2.41A FPGA_Version: V2.1.11  NB_Of_Core_Daughters: 1 EXT_CLOCK: 0 ( 1 = YES, 0= NO) ==
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> detector.usb_sernum >> tmpst >> detector.fpga_ver >> tmpst >> detector.n_core_daughters >> tmpst >> detector.external_clock >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst;
  if(debug) cout<<" COREMODULE: USB_SerNum " <<detector.usb_sernum<<", FPGA Ver "<<detector.fpga_ver<<", NB_Of_Core_Daughters "<< detector.n_core_daughters<<", EXT_CLOCK "<<detector.external_clock<< "( 1 = YES, 0= NO)"<<endl;
  //CRP Analysing the second line
  //return to beginning of line
  //bool iscore(false);
  reading_file.clear();
  reading_file.seekg(begin);
  getline (reading_file,line);
  std::cout << "CRP Second line: " << line << std::endl; 
  
  //CRP 14/12/20 Check whether we have a CORE Module or a direct interface 
  found = line.find("SL_COREMODULE_INTERFACE");
  //of we find a direct interface we retrieve here the number of connected ASUs
  int nconnasu(0);
  if (found!=std::string::npos) detector.isCore=true; 
  //CRP 14/12/20 If we don't find a CORE Module we suppose it's the direct interface (if there is a third
  //interface this part of the code will have to be adjusted
  else {
    found = line.find("Nb_Of_Connected_ASUs:");
    if (found!=std::string::npos) { 
      found1 = line.find_first_of("123456789",found);
      found2 = line.find_first_of(' ',found1);
      //std::string nconnasu_str = line.substr(found1,found2-found1);;
      nconnasu = std::stoi(line.substr(found1,found2-found1));
      std::cout << "CRP Number of connected ASUs: " << nconnasu << std::endl;
    } else std::cout << "Attention: DIRECT INTERFACE BUT NUMBER OF CONNECTED ASUS NOT GIVEN" << std::endl; 
  }

  std::cout << "CRP iscore ?: " << detector.isCore << std::endl; 
  
  //third line
  //----------------------
  //=== TriggerType: 0  ('0' = FORCE_TRIGGER, '1' = SELF_TRIGGER) ACQWindowSource: 0 ('0' = AUTO, '1'= SOFT_CMD, '2' = EXT_SIG) ACQWindow: 10 (ms) DelayBetweenCycle: 100 (ms) DelayForStartAcq: 0 (5Mhz_Clock period) ExtSigLevel: 0 ('0' = TTL, '1' = NIM, '-1' = NA) ExtSigEdge: 0 ('0' = RISING_EDGE, '1' = FALLING_EDGE, '-1' = NA) ==
  reading_file >> tmpst >> tmpst >> detector.trigger_type >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> detector.acq_window_source >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> detector.acq_window >> tmpst >> tmpst >> detector.acq_delay >> tmpst >> tmpst >> detector.start_acq_delay >> tmpst  >> tmpst >> tmpst >>detector.ext_sig_level >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> detector.ext_sig_edge >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst ;
  if(debug) cout<<" Trigger Type: "<<detector.trigger_type << " AcqWindowSource ('0' = AUTO, '1'= SOFT_CMD, '2' = EXT_SIG): " <<detector.acq_window_source << " AcqWindow: "<<detector.acq_window << "ms DelayBetweenCycle: "<<detector.acq_delay<< " (ms) DelayForStartAcq:"<<detector.start_acq_delay<<" (5Mhz_Clock period) ExtSigLevel ('0' = TTL, '1' = NIM, '-1' = NA): " << detector.ext_sig_level << " ExtSigEdge ('0' = RISING_EDGE, '1' = FALLING_EDGE, '-1' = NA): " << detector.ext_sig_edge <<endl;

  //----------------------
  //CRP 14/12/20 This line is not present in case of a direct interface (at least not at the moment)
  if(detector.isCore) {
    //=== CORE_DAUGHTER: 0 FPGA_Version: V1.3.1 Nb_Of_Connected_Slabs: 6 ===
    for(int i=0; i<detector.n_core_daughters; i++) {
      int icore=-1;
      reading_file >> tmpst >> tmpst  >> icore  >> tmpst >> detector.core_daughter_fpga_ver[i] >> tmpst  >> detector.core_daughter_n_slabs[i] >> tmpst;
      if(debug) cout<<" CORE_DAUGHTER: " << icore <<" FPGA_Version: " << detector.core_daughter_fpga_ver[i] << " Nb_Of_Connected_Slabs: "<<detector.core_daughter_n_slabs[i]<< endl;
    }
  }
  //----------------------
  //CRP 14/12/20 Now we make a terrible hack
  //IF we have the direct interface we set the number of core daughters to 1 and and the number of slabs to 1
  if(!detector.isCore) {
    detector.n_core_daughters=1;
    detector.core_daughter_n_slabs[0]=1;
    detector.slab[0][0].nb_asus=nconnasu;
  }
  //CRP
  //===== Daughter: 0 SlabIdx: 0 SlabAdd: 1 SL_Board_SerNum: -1 FPGA_Version: V2.3.7  Nb_Of_Connected_ASUs: 1 =====
  for(int i=0; i<detector.n_core_daughters; i++) {
    for(int j=0; j<detector.core_daughter_n_slabs[i]; j++) {
      //CRP in case of direct interface we have to prevent the program from reading the next line 
      if(detector.isCore) {
	reading_file >> tmpst >> tmpst >> tmpst >>  tmpst >> detector.slab[i][j].idx >> tmpst >> detector.slab[i][j].add >> tmpst >> detector.slab[i][j].ser_num >> tmpst >> detector.slab[i][j].slb_fpga >> tmpst >> detector.slab[i][j].nb_asus>>tmpst ;
	if(debug) cout<< " Daughter: " << i << " SlabIdx: "<< detector.slab[i][j].idx<< " SlabAdd: "<<detector.slab[i][j].add << " SL_Board_SerNum: " << detector.slab[i][j].ser_num << " FPGA_Version: "<<detector.slab[i][j].slb_fpga << " Nb_ASUS: "<<detector.slab[i][j].nb_asus<<endl;
      } 
      //CRP
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
	  }//End of loop over channels
	}//End of loop over chips of an ASU
      }//End of loop over ASUs
    }//End of loop over slabs connected to a core daughter
  }//End of loop over core daughters  
}

void mask_full_chip(int idaughter, int islab, int iasu, int ichip) {
  for(int i=0; i<detector.slab[idaughter][islab].asu[iasu].skiroc[ichip].n_channels; i++) detector.slab[idaughter][islab].asu[iasu].skiroc[ichip].mask[i]=1;
}


void enable_trig_row(int idaughter, int islab, int iasu, int ichip, int irow) {
  //row 0 --> chns 0-7
  //row 2--> chns 8-15
  //etc
  cout<<irow<<endl;
  for(int i=irow*8; i<(irow*8+8); i++) {
    cout<<idaughter<<" "<<islab<<" "<<iasu<<" "<< ichip<<" "<<i <<endl;
    detector.slab[idaughter][islab].asu[iasu].skiroc[ichip].mask[i]=0;
  }
}

TString trigger_type_string(int i) {

  TString result="SELF_TRIGGER";
  if(i==0) result="FORCE_TRIGGER";
  if(i==1) result ="SELF_TRIGGER";

  return result;
}

TString external_clock_string(int i) {

  TString result="NO";
  if(i==1) result="YES";
  return result;
}

TString acq_window_source_string(int i) {
  
  //('0' = AUTO, '1'= SOFT_CMD, '2' = EXT_SIG)
  TString result="AUTO";
  if(i==0) result="AUTO";
  if(i==1) result ="SOFT_CMD";
  if(i==2) result ="EXT_SIG";

  return result;
}

TString ext_sig_level_string(int i) {
  
  //  ('0' = TTL, '1' = NIM, '-1' = NA)
  TString result="TTL";
  //  if(i==0) result="TTL";
  //if(i==1) result ="NIM";
  //if(i==-1) result ="NA";

  return result;
}


TString ext_sig_edge_string(int i) {
  
  //  ('0' = RISING_EDGE, '1' = FALLING_EDGE, '-1' = NA)
  TString result="POS";
  // if(i==0) result="RISING_EDGE";
  //if(i==1) result ="FALLING_EDGE";
  //if(i==-1) result ="NA";
  
  return result;
}


void write_configuration_file(TString filename="Run_Settings_2.txt") {

  //===== Daughter: 0 SlabIdx: 0 SlabAdd: 1 SL_Board_SerNum: -1 FPGA_Version: V2.3.7  Nb_Of_Connected_ASUs: 1 =====
ofstream fout(filename,ios::out);

//fout<<"== SETTINGS FILE SAVED WITH ILC SL-SOFTWARE VERSION: "<<detector.sl_soft_v<<"  == DATE OF RUN: UnixTime = "<<detector.date_run_unix<<" date = "<<detector.date_run_date<<" time = "<<detector.date_run_time<<"  ==="<<endl;
//fout<<"== SYSTEM_TYPE: SL_COREMODULE_INTERFACE USB_SerNum: "<<detector.usb_sernum<<" FPGA_Version: "<<detector.fpga_ver<<"  NB_Of_Core_Daughters: "<<detector.n_core_daughters<<" =="<<endl;
//fout<<"== TriggerType: "<<detector.trigger_type<<"  ('0' = FORCE_TRIGGER, '1' = SELF_TRIGGER) ACQWindowSource: "<<detector.acq_window_source <<" ('0' = AUTO, '1'= SOFT_CMD, '2' = EXT_SIG) ACQWindow: "<< detector.acq_window <<" (ms) DelayBetweenCycle: "<< detector.acq_delay <<" (ms) ExtSigLevel: "<< detector.ext_sig_level  <<" ('0' = TTL, '1' = NIM, '-1' = NA) ExtSigEdge: "<< detector.ext_sig_edge <<" ('0' = RISING_EDGE, '1' = FALLING_EDGE, '-1' = NA) =="<<endl;

 fout<<"# AcqParams TrigType ' "<<trigger_type_string(detector.trigger_type)<<" ' AcqWindowSource ' "<<acq_window_source_string(detector.acq_window_source) <<" ' ACQWindow "<< detector.acq_window <<" DelayBetweenCycles "<< detector.acq_delay <<" DelayForStartAcq 0 ExtSigLvel ' "<< ext_sig_level_string(detector.ext_sig_level)  <<" ' ExtSigEdge ' "<< ext_sig_edge_string(detector.ext_sig_edge) <<" '"<<endl;

 //CRP 14/12/20 This is a terrible hack!!!  
 //It is to make use of the current organisation of the code
 //This clearly has to be rewritten soon the more since the code that fills the detector struct also tacitly supposed only
 //one CORE Daughter and will thus collapse as sonn as a second CORE Daughter will be added!!!!
 if(!detector.isCore) {
   std::cout << "No CORE Module present" << std::endl;
   detector.n_core_daughters = 1;
 }
 //
 //----------------------
  //=== CORE_DAUGHTER: 0 FPGA_Version: V1.3.1 Nb_Of_Connected_Slabs: 6 ===
   for(int i=0; i<detector.n_core_daughters; i++) {
     //  fout<<"=== CORE_DAUGHTER: "<<i<<" FPGA_Version: "<<detector.core_daughter_fpga_ver[i]<<" Nb_Of_Connected_Slabs: "<<detector.core_daughter_n_slabs[i]<<" ==="<<endl;
     std::cout << "In writing loop" << std::endl;
     for(int j=0; j<detector.core_daughter_n_slabs[i]; j++) {
     std::cout << "In writing loop 2" << std::endl;
       for(int k=0; k<detector.slab[i][j].nb_asus; k++) {
     std::cout << "In writing loop 3" << std::endl;
	 //	fout<<"# SlabSerNum "<<detector.slab[i][j].ser_num<<" SlabAdd "<< detector.slab[i][j].add <<" Asu "<<detector.slab[i][j].asu[k].idx<<endl;
	 fout<<"# SlabSerNum 0.0 SlabAdd "<< detector.slab[i][j].add <<" Asu "<<detector.slab[i][j].asu[k].idx<<endl;
	 for(int ichip=0; ichip<detector.slab[i][j].asu[k].n_chips; ichip++) {
	   fout<<"## ChipIndex "<< detector.slab[i][j].asu[k].skiroc[ichip].idx<<" FeedbackCap "<<detector.slab[i][j].asu[k].skiroc[ichip].feedback_cap<<" ThresholdDAC "<<detector.slab[i][j].asu[k].skiroc[ichip].threshold_dac <<" HoldDelay "<<  detector.slab[i][j].asu[k].skiroc[ichip].hold_delay <<" FSPeakTime "<< detector.slab[i][j].asu[k].skiroc[ichip].fs_peak_time <<" GainSelectionThreshold "<< detector.slab[i][j].asu[k].skiroc[ichip].gain_selection_threshold <<" DefaultChThreshold "<<detector.slab[i][j].asu[k].skiroc[ichip].default_ch_threshold <<" DefaultTrigMask "<<detector.slab[i][j].asu[k].skiroc[ichip].default_trig_mask<<" DefaultPAMask "<< detector.slab[i][j].asu[k].skiroc[ichip].default_pa_mask <<endl;
	  
	   for(int ichn=0; ichn<detector.slab[i][j].asu[k].skiroc[ichip].n_channels; ichn++) {
	     fout<<"### Ch "<<ichn<<" TrigMask "<<detector.slab[i][j].asu[k].skiroc[ichip].mask[ichn]<<" ChThreshold "<< detector.slab[i][j].asu[k].skiroc[ichip].chn_threshold_adj[ichn] <<" PAMask "<< detector.slab[i][j].asu[k].skiroc[ichip].preamplifier_mask[ichn] <<" "<<endl;
	   }
	 }
       }
     }//loop over slabs
   }//loop over core daughters

}
