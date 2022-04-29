//# Copyright 2020 Adrián Irles IJCLab (CNRS/IN2P3)

#ifndef SLBdecoded2ROOT_CC
#define SLBdecoded2ROOT_CC

#include <iomanip>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <bitset>
#include <math.h>
#include <vector>
#include <string>
#include <algorithm>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <sstream>
#include <string>


// .x SLBdecoded2ROOT.C+

using std::cout;
using std::endl;

#define NB_OF_SKIROCS_PER_ASU 16
#define NB_OF_CHANNELS_IN_SKIROC 64
#define NB_OF_SCAS_IN_SKIROC 15
#define SINGLE_SKIROC_EVENT_SIZE (129*2) 
#define SLBDEPTH 15
#define NB_CORE_DAUGHTERS 1 
#define NEGDATA_THR 11
#define BCIDTHRES 3

class SLBdecoded2ROOT {

public:
  SLBdecoded2ROOT(){
  };
  ~SLBdecoded2ROOT(){

  };
  
  void ReadFile(TString inputFileName, bool overwrite=false, TString outFileName = "default",bool zerosupression=false);

protected:


  int R2Rstate;

  void Initialisation();
  void treeInit(bool);
  //  int  GetTree(TString rootfilename);
  void GetBadBCID();

  TFile* fout;
  TTree* tree;

  TFile *finroot;
  TTree* slboardread;

  int _bcid[SLBDEPTH][NB_OF_SKIROCS_PER_ASU][NB_OF_SCAS_IN_SKIROC];
  int _corrected_bcid[SLBDEPTH][NB_OF_SKIROCS_PER_ASU][NB_OF_SCAS_IN_SKIROC];
  int _badbcid[SLBDEPTH][NB_OF_SKIROCS_PER_ASU][NB_OF_SCAS_IN_SKIROC];
  int _nhits[SLBDEPTH][NB_OF_SKIROCS_PER_ASU][NB_OF_SCAS_IN_SKIROC];
  int _adc_low[SLBDEPTH][NB_OF_SKIROCS_PER_ASU][NB_OF_SCAS_IN_SKIROC][NB_OF_CHANNELS_IN_SKIROC];
  int _adc_high[SLBDEPTH][NB_OF_SKIROCS_PER_ASU][NB_OF_SCAS_IN_SKIROC][NB_OF_CHANNELS_IN_SKIROC];
  int _autogainbit_low[SLBDEPTH][NB_OF_SKIROCS_PER_ASU][NB_OF_SCAS_IN_SKIROC][NB_OF_CHANNELS_IN_SKIROC];
  int _autogainbit_high[SLBDEPTH][NB_OF_SKIROCS_PER_ASU][NB_OF_SCAS_IN_SKIROC][NB_OF_CHANNELS_IN_SKIROC];
  int _hitbit_low[SLBDEPTH][NB_OF_SKIROCS_PER_ASU][NB_OF_SCAS_IN_SKIROC][NB_OF_CHANNELS_IN_SKIROC];
  int _hitbit_high[SLBDEPTH][NB_OF_SKIROCS_PER_ASU][NB_OF_SCAS_IN_SKIROC][NB_OF_CHANNELS_IN_SKIROC];
  //  int _event;
  //int numbcid;
  int _numCol[SLBDEPTH][NB_OF_SKIROCS_PER_ASU];
  int _chipId[SLBDEPTH][NB_OF_SKIROCS_PER_ASU];
  int _slot[SLBDEPTH];
  int _slboard_id[SLBDEPTH];
  int _n_slboards;
  int _acqNumber;

  //slowcontrol
  double _startACQ[SLBDEPTH];
  int _rawTSD[SLBDEPTH];
  float _TSD[SLBDEPTH];
  int _rawAVDD0[SLBDEPTH];
  int _rawAVDD1[SLBDEPTH];
  float _AVDD0[SLBDEPTH];
  float _AVDD1[SLBDEPTH];


  //  InfoChip * info;
};

//******************************************************************************************************************

void SLBdecoded2ROOT::Initialisation() {
  fout->cd(); R2Rstate=-1;

  tree = new TTree("siwecaldecoded","siwecaldecoded");

  //    tree->Branch("event",&_event,"event/I");
  tree->Branch("acqNumber",&_acqNumber,"acqNumber/I");
  tree->Branch("n_slboards",&_n_slboards,"n_slboards/I");

  TString name;
  name= TString::Format("slot[%i]/I",SLBDEPTH);
  tree->Branch("slot",_slot,name);
  
  name= TString::Format("slboard_id[%i]/I",SLBDEPTH);
  tree->Branch("slboard_id",_slboard_id,name);
  
  name= TString::Format("chipid[%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU);
  tree->Branch("chipid",_chipId,name);

  name= TString::Format("nColumns[%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU);
  tree->Branch("nColumns",_numCol,name);

  name= TString::Format("startACQ[%i]/F",SLBDEPTH);
  tree->Branch("startACQ",_startACQ,name);

  name= TString::Format("rawTSD[%i]/I",SLBDEPTH);
  tree->Branch("rawTSD",_rawTSD,name);

  name= TString::Format("TSD[%i]/F",SLBDEPTH);
  tree->Branch("TSD",_TSD,name);
  
  name= TString::Format("rawAVDD0[%i]/I",SLBDEPTH);
  tree->Branch("rawAVDD0",_rawAVDD0,name);

  name= TString::Format("rawAVDD1[%i]/I",SLBDEPTH);
  tree->Branch("rawAVDD1",_rawAVDD1,name);

  name= TString::Format("AVDD0[%i]/F",SLBDEPTH);
  tree->Branch("AVDD0",_AVDD0,name);

  name= TString::Format("AVDD1[%i]/F",SLBDEPTH);
  tree->Branch("AVDD1",_AVDD1,name);
  
  name= TString::Format("bcid[%i][%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU,NB_OF_SCAS_IN_SKIROC);
  tree->Branch("bcid",_bcid,name);

  name= TString::Format("corrected_bcid[%i][%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU,NB_OF_SCAS_IN_SKIROC);
  tree->Branch("corrected_bcid",_corrected_bcid,name);

  name= TString::Format("badbcid[%i][%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU,NB_OF_SCAS_IN_SKIROC);
  tree->Branch("badbcid",_badbcid,name);

  name= TString::Format("nhits[%i][%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU,NB_OF_SCAS_IN_SKIROC);
  tree->Branch("nhits",_nhits,name);

  name= TString::Format("adc_low[%i][%i][%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU,NB_OF_SCAS_IN_SKIROC,NB_OF_CHANNELS_IN_SKIROC);
  tree->Branch("adc_low",_adc_low,name);

  name= TString::Format("adc_high[%i][%i][%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU,NB_OF_SCAS_IN_SKIROC,NB_OF_CHANNELS_IN_SKIROC);
  tree->Branch("adc_high",_adc_high,name);

  name= TString::Format("autogainbit_low[%i][%i][%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU,NB_OF_SCAS_IN_SKIROC,NB_OF_CHANNELS_IN_SKIROC);
  tree->Branch("autogainbit_low",_autogainbit_low,name);

  name= TString::Format("autogainbit_high[%i][%i][%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU,NB_OF_SCAS_IN_SKIROC,NB_OF_CHANNELS_IN_SKIROC);
  tree->Branch("autogainbit_high",_autogainbit_high,name);
    
  name= TString::Format("hitbit_low[%i][%i][%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU,NB_OF_SCAS_IN_SKIROC,NB_OF_CHANNELS_IN_SKIROC);
  tree->Branch("hitbit_low",_hitbit_low,name);

  name= TString::Format("hitbit_high[%i][%i][%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU,NB_OF_SCAS_IN_SKIROC,NB_OF_CHANNELS_IN_SKIROC);
  tree->Branch("hitbit_high",_hitbit_high,name);

  return;

}






//******************************************************************************************************************

void SLBdecoded2ROOT::treeInit(bool zerosupression=false) { //init data for a single SPILL ?

  for (int isl=0; isl<SLBDEPTH; isl++) {
    for (int k=0; k<NB_OF_SKIROCS_PER_ASU; k++) {
      for (int i=0; i<NB_OF_SCAS_IN_SKIROC; i++) {
	_bcid[isl][k][i]=-999;
	_badbcid[isl][k][i]=-999;
	_corrected_bcid[isl][k][i]=-999;
	_nhits[isl][k][i]=-999;
	for (int j=0; j<NB_OF_CHANNELS_IN_SKIROC; j++) {
	  _adc_low[isl][k][i][j]=-999;
	  _adc_high[isl][k][i][j]=-999;
	  _autogainbit_low[isl][k][i][j]=-999;
	  _autogainbit_high[isl][k][i][j]=-999;
	  _hitbit_low[isl][k][i][j]=-999;
	  _hitbit_high[isl][k][i][j]=-999;

	}
      }
      
      _chipId[isl][k]=-999;
      _numCol[isl][k]=0;
      _startACQ[isl]=-1;
      _rawTSD[isl]=-1;
      _TSD[isl]=-1;
      _rawAVDD0[isl]=-1;
      _rawAVDD1[isl]=-1;
      _AVDD0[isl]=-1;
      _AVDD1[isl]=-1;
    }
    _slot[isl]=-1;
    _slboard_id[isl]=-1;
  }
  _n_slboards=-1;
  _acqNumber=-1;

  return;
}


//******************************************************************************************************************


//******************************************************************************************************************

void SLBdecoded2ROOT::ReadFile(TString inputFileName, bool overwrite, TString outFileName, bool zerosupression) {

  //  event=0;
  _acqNumber=0;
  
  if(outFileName == "default"){
    outFileName = TString::Format("%s.root",inputFileName.Data());
    //cout<<outFileName<<endl;
  }
  
  if(!overwrite){
    fout = new TFile(outFileName,"create");
    if(!fout->IsOpen()){
      return;
    }
  }
  else {
    fout = new TFile(outFileName,"recreate");
  }
  

  std::ifstream reading_file(inputFileName);
 
  if(!reading_file){
    cout<<" ERROR  ----------------------------------- No file: "<<inputFileName<<endl;
    return 0;
  } else {
    cout<<" Read File "<<inputFileName<<endl;
    //cout<<" slb= " <<slboard_index<<endl;
  }

  cout<<" ##### ##### ##### ##### "<<endl;
  cout<<" #####  WARNING  ######"<<endl;
  cout<<" # --> this code is deprecated. The converter is not dealing properly with interleaved data (cycleIDs that arrived to the aggregator interleaved from different SLBoards" <<endl;
  cout<<" # PLEASE USE THE BINARY FORMAT to record your data and the SLBraw2ROOT.cc macro "<<endl;
  cout<<" # A. Irles 29/04/2022 "<<endl;
  cout<<" #####  WARNING  ######"<<endl;
  cout<<" ##### ##### ##### ##### "<<endl;
  

  //FTDI
  int cycleID=-1;
  int tmp_n_slboards=0;
  TString tmpst;

  std::string strheader;
  std::getline(reading_file, strheader);
  std::cout << "header 1: " << strheader << std::endl;
  std::getline(reading_file, strheader);
  std::cout << "header 2: " << strheader << std::endl;
  std::cout << "FTDI Interface? " << strheader.find("SL_DIRECT_INTERFACE") << std::endl;
  std::cout << "COREMODULE " << strheader.find("COREMODULE") << std::endl;
  std::cout << "std::string::npos " << std::string::npos   << std::endl;
  TString readout_type="COREMODULE";
  //CRP 24/11/20 Correct analysis of std::string  
  //FTDI readout DIRECT Interface, only one slab
  std::size_t found = strheader.find("SL_DIRECT_INTERFACE");
  if (found!=std::string::npos) {
    //if(strheader.find("SL_DIRECT_INTERFACE")!=0) {
    std::cout << "FTDI found" << std::endl;
    readout_type="FTDI";
    tmp_n_slboards=1;
    std::getline(reading_file, strheader);
    std::getline(reading_file, strheader);
    std::getline(reading_file, strheader);
    std::getline(reading_file, strheader);
  } else {
    reading_file >> tmpst >> tmpst >> tmpst  >> tmpst >>  tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst  >> tmpst >> tmp_n_slboards >> tmpst;
    cout<<" NB OF CONNECTED SLABs = " << tmp_n_slboards <<endl;
    std::getline(reading_file, strheader);

    for(int islboard=0; islboard<tmp_n_slboards; islboard++) {
      std::getline(reading_file, strheader);
      std::cout<<islboard<<" "<<strheader<< std::endl;
    }

    std::getline(reading_file, strheader);
    std::getline(reading_file, strheader);
    std::getline(reading_file, strheader);
    std::getline(reading_file, strheader);

  }

  
  // Initialisation((const int)tmp_n_slboards);
  Initialisation();

  
  // int bcid_cycle[20][16][15];

 

  while (reading_file) {
    // output the line
    TString tmpst1;
    TString tmpst;
    int size=0;
    int chip=-1;
    int slabidx=-1;
    int slabadd=-1;
    double start_acq=-1;
    int raw_tsd=-1;
    float _tsd=-1;
    int raw_avdd0=-1;
    int raw_avdd1=-1;
    float _avdd0=-1;
    float _avdd1=-1;

    reading_file >> tmpst1 >> tmpst >> size >> tmpst >> chip >> tmpst >>  tmpst >> tmpst >> slabidx >> tmpst >> slabadd >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst  >> tmpst >> tmpst >> cycleID >> tmpst >> start_acq >> tmpst >> raw_tsd >> tmpst >> raw_avdd0 >>  tmpst >> raw_avdd1 >> tmpst >> _tsd >> tmpst >> _avdd0 >>  tmpst >> _avdd1;
    //    if(slabidx>20)
    //std::cout << " -----" <<tmpst1<<" "<<size<< " "<<chip<<" "<<slabidx<<" "<<slabadd<<" "<<cycleID <<" "<<acqNumber<<" "<<start_acq<<" "<<tmpst1<<std::endl;
    //  if(slabidx<0) break;
    if( readout_type=="COREMODULE" && slabidx==-1) break;//ftdi
    if( readout_type=="FTDI"  && slabidx!=-1) break;
    if( readout_type=="FTDI") {
      slabidx=0;
      slabadd=0;
    }

    if(_acqNumber==0) treeInit(zerosupression);
    if(_acqNumber>0 && _acqNumber!=cycleID) {
      GetBadBCID();
      tree->Fill();
      treeInit(zerosupression);
    }

    

    //save variables
    _acqNumber=cycleID;
    _n_slboards=tmp_n_slboards;
    //    event++;
    int previousBCID=-1000;
    int loopBCID=0;

    _startACQ[slabadd]=start_acq;
    _rawTSD[slabadd]=raw_tsd;
    _TSD[slabadd]=_tsd;
    _rawAVDD0[slabadd]=raw_avdd0;
    _rawAVDD1[slabadd]=raw_avdd1;
    _AVDD0[slabadd]=raw_avdd0;
    _AVDD1[slabadd]=raw_avdd1;

    _slot[slabadd]=-1;//we don't know this information... the DQ only provides addresses.
    _slboard_id[slabadd]=slabadd;

    
    for(int i=0; i<size; i++) {
      //##0 BCID 1627 SCA 1 #Hits 49
      int bcid_tmp=-1;
      int sca=-1;
      int nhits_tmp=-1;
      reading_file >> tmpst >> tmpst >> bcid_tmp >> tmpst >> sca >>  tmpst  >>  nhits_tmp ;
      //cout<<bcid_tmp<<" "<<sca<<" "<<nhits_tmp<<" "<<endl;
      sca=size-(sca+1);
      _bcid[slabadd][chip][sca]=bcid_tmp;
      _nhits[slabadd][chip][sca]=nhits_tmp;
      _numCol[slabadd][chip]++;
      
      if(_bcid[slabadd][chip][sca] > 0 && _bcid[slabadd][chip][sca]-previousBCID < 0) loopBCID++;
      if(_bcid[slabadd][chip][sca] > 0 ) _corrected_bcid[slabadd][chip][sca] = _bcid[slabadd][chip][sca]+loopBCID*4096;
      previousBCID=bcid_tmp;
      
      if(chip>-1 && chip<16) {
	_chipId[slabadd][chip]=chip;
      } else {
	cout<<"Wrong chipID = "<<chip<<endl;
	break;
      }
      //      }
      int nchn=0;
      if(zerosupression==false) nchn=NB_OF_CHANNELS_IN_SKIROC;
      else nchn = nhits_tmp;
      //Ch 0 LG 322 1 0 HG 408 1 0
      for(int ichn=0; ichn<nchn; ichn++) {
	int chn=-1;
	int lg=-1;
	int hg=-1;
	int lg_hbit=-1;
	int hg_hbit=-1;
	int lg_gbit=-1;
	int hg_gbit=-1;
	reading_file >> tmpst >> chn >> tmpst >> lg >>  lg_hbit >> lg_gbit >>  tmpst >> hg >> hg_hbit >> hg_gbit;
	//cout << chn << " "<< hg <<" " <<hg_bit<<endl;
	//	if(slabid==slboard_index) {
	_adc_low[slabadd][chip][sca][ichn]=lg;
	_autogainbit_low[slabadd][chip][sca][ichn]=lg_gbit;
	_hitbit_low[slabadd][chip][sca][ichn]=lg_hbit;
	_adc_high[slabadd][chip][sca][ichn]=hg;
	_autogainbit_high[slabadd][chip][sca][ichn]=hg_gbit;
	_hitbit_high[slabadd][chip][sca][ichn]=hg_hbit;

	//}
      }//end channel  
    }//end events of the chi

  }
  
  fout->cd(); 
  fout->Write(0);
  //  for(int ichip=0; ichip<16; ichip++) {
  //   bcid_diff[ichip]->Write();
  //    bcid_correl[ichip]->Write();
  //  }
  fout->Close();
    
  return;
}

 
void SLBdecoded2ROOT::GetBadBCID() {

  //add tags
  //  int count_negdata[SLBDEPTH];
  //  for(int i=0; i<SLBDEPTH; i++) count_negdata[i]=0;

  for(int i=0; i<SLBDEPTH; i++) {
    for (int k=0; k<NB_OF_SKIROCS_PER_ASU; k++) {
      //only for valid chips in this spill
      if (_chipId[i][k]>=0) {
	for (int ibc=0; ibc<_numCol[i][k]; ibc++) {

	  // if sca+1 is filled with consec bcid, and sca+2 also, then _badbcid[sca]==3 && _badbcid[sca+1]==3 (retriggering)
	  // if sca+1 is not filled with consec bcid,  _badbcid==0
	  // if sca and sca+1 have consecutive bcids but are not retriggers, then we consider 3 types of EMPTY events
	  //     case A: empty events after a event with triggers
	  //     case B: empty events before a event with triggers --> TO BE UNDERSTOOD but it seems that the good one is the one with the trigger (the second)
	  //     case C&D: both bcids have triggers but, we do nothing
	    
	  if(ibc==0) {
	    _badbcid[i][k][ibc]=0;
	    int _corr_bcid=_corrected_bcid[i][k][ibc];
	    int _corr_bcid1=0;
	    int _corr_bcid2=0;

	    if(_corrected_bcid[i][k][ibc+1]>0 && _corrected_bcid[i][k][ibc]>0 && (_corrected_bcid[i][k][ibc+1]-_corrected_bcid[i][k][ibc])>0) 
	      _corr_bcid1=_corrected_bcid[i][k][ibc+1];

	    if(_corrected_bcid[i][k][ibc+2]>0 && (_corrected_bcid[i][k][ibc+2]-_corrected_bcid[i][k][ibc+1])>0) 
	      _corr_bcid2=_corrected_bcid[i][k][ibc+2];

	    if(_corr_bcid2>0) {
	      //empty events
	      if( ( _corr_bcid2-_corr_bcid1) >(BCIDTHRES - 1) && (_corr_bcid1-_corr_bcid) ==1) {
		//case A: empty events after a event with triggers
		if(_nhits[i][k][ibc]>0 && _nhits[i][k][ibc+1]==0) {
		  _badbcid[i][k][ibc]=0;
		  _badbcid[i][k][ibc+1]=2;//this one is to not be used in any case
		}
		//case B: empty events before a event with triggers --> TO BE UNDERSTOOD but it seems that the good one is the one with trigger (the second bcid)
		if(_nhits[i][k][ibc]==0 && _nhits[i][k][ibc+1]>0) {
		  _badbcid[i][k][ibc]=1; //empty event
		  _badbcid[i][k][ibc+1]=0;//good one
		}
		//case C&D: both bcids have triggers but, we do nothing
	      }
	      // pure retriggers
	      if( ( _corr_bcid2-_corr_bcid1) < BCIDTHRES && (_corr_bcid1-_corr_bcid) < BCIDTHRES) {
		_badbcid[i][k][ibc]=3;
		_badbcid[i][k][ibc+1]=3;
		_badbcid[i][k][ibc+2]=3;
	      }
	    }

	    // pure retriggers
	    if( _corr_bcid1 > 0 && (_corr_bcid1-_corr_bcid) > 1 && (_corr_bcid1-_corr_bcid) <BCIDTHRES) {
	      _badbcid[i][k][ibc]=3;
	      _badbcid[i][k][ibc+1]=3;
	    }  
	  } //ibc==0 if 

	  if(ibc>0 && _badbcid[i][k][ibc]<0 && _corrected_bcid[i][k][ibc] >0 &&  (_corrected_bcid[i][k][ibc]-_corrected_bcid[i][k][ibc-1])>0 ) {
	    _badbcid[i][k][ibc]=0;
	    int _corr_bcid=_corrected_bcid[i][k][ibc];
	    int _corr_bcidminus=_corrected_bcid[i][k][ibc-1];

	    if(_corrected_bcid[i][k][ibc+1]>0 && (_corrected_bcid[i][k][ibc+1]-_corrected_bcid[i][k][ibc])>0) {
	      int _corr_bcid1=_corrected_bcid[i][k][ibc+1];

	      if(_corrected_bcid[i][k][ibc+2]>0 && (_corrected_bcid[i][k][ibc+2]-_corrected_bcid[i][k][ibc+1])>0) {
		int _corr_bcid2=_corrected_bcid[i][k][ibc+2];
		if( ( _corr_bcid2-_corr_bcid1) < BCIDTHRES && (_corr_bcid1-_corr_bcid) < BCIDTHRES) {
		  _badbcid[i][k][ibc]=3;
		  _badbcid[i][k][ibc+1]=3;
		  _badbcid[i][k][ibc+2]=3;
		}
		if( (_corr_bcid1-_corr_bcid) < BCIDTHRES && (_corr_bcid-_corr_bcidminus) < BCIDTHRES ) _badbcid[i][k][ibc]=3;

		if( _badbcid[i][k][ibc]!=3 && ( _corr_bcid2-_corr_bcid1) >(BCIDTHRES - 1) && (_corr_bcid1-_corr_bcid) ==1) {
		  //case A: empty events after a event with triggers
		  if(_nhits[i][k][ibc]>0 && _nhits[i][k][ibc+1]==0) {
		    _badbcid[i][k][ibc]=0;
		    _badbcid[i][k][ibc+1]=2;//this one is to not be used in any case
		  }
		  //case B: empty events before a event with triggers --> TO BE UNDERSTOOD but it seems that the good one is the one with the trgger (the second bcid)
		  if(_nhits[i][k][ibc]==0 && _nhits[i][k][ibc+1]>0) {
		    _badbcid[i][k][ibc]=1; //empty event before the real event
		    _badbcid[i][k][ibc+1]=0;//good event
		  }
		  //case C&D: both bcids have triggers: we do nothing
		}
		if( _badbcid[i][k][ibc]!=3 && ( _corr_bcid2-_corr_bcid1) >(BCIDTHRES - 1) && (_corr_bcid1-_corr_bcid) > 1 && (_corr_bcid1-_corr_bcid) <BCIDTHRES) {
		  _badbcid[i][k][ibc]=3;
		  _badbcid[i][k][ibc+1]=3;
		}
		if( (_corr_bcid-_corr_bcidminus) < BCIDTHRES ) _badbcid[i][k][ibc]=3;

		//if( _badbcid[i][k][ibc-1]==1 && (_corr_bcid1-_corr_bcid) > (BCIDTHRES - 1)) _badbcid[i][k][ibc]=2;
	      } else {
		if( (_corr_bcid1-_corr_bcid) < BCIDTHRES ) _badbcid[i][k][ibc]=3;
		if( (_corr_bcid-_corr_bcidminus) < BCIDTHRES ) _badbcid[i][k][ibc]=3;
	      } //ibc+2 if
	    } else {
	      if( (_corr_bcid-_corr_bcidminus) < BCIDTHRES ) _badbcid[i][k][ibc]=3;
	    }//ibc+1 if
	  } //ibc>0 if 

      
	    //tag zero (under/over flow) data
	    //	count_negdata=0;
	
	    //	for (int ichan=0; ichan<NB_OF_CHANNELS_IN_SKIROC; ichan++) {
	    //  if  (adc_high[i][k][ibc][ichan] < NEGDATA_THR) count_negdata++;
	    //}//ichan

	    //if (count_negdata>0) {_badbcid[i][k][ibc]+=32;}

	}//ibc

      }//chipId
    }//k
  }//i   slboard
  
}

#endif

