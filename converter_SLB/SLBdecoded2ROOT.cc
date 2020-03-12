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

class SLBdecoded2ROOT {

public:
  SLBdecoded2ROOT(){
  };
  ~SLBdecoded2ROOT(){

  };
  
  void ReadFile(TString inputFileName, bool overwrite=false, TString outFileName = "default",bool zerosupression=false);

protected:

  enum {
    SLBDEPTH=15,
    MEMDEPTH=15,
    NCHANNELS=64,
    NCHIP=16,
    NEGDATA_THR=11, //event with data below are tagged badbcid+=32
    BCIDTHRES=15
  };

  int R2Rstate;

  void Initialisation();
  void treeInit(bool);
  int  GetTree(TString rootfilename);
  void GetBadBCID();

  TFile* fout;
  TTree* tree;

  TFile *finroot;
  TTree* slboardread;

  int bcid[SLBDEPTH][NCHIP][MEMDEPTH];
  int corrected_bcid[SLBDEPTH][NCHIP][MEMDEPTH];
  int badbcid[SLBDEPTH][NCHIP][MEMDEPTH];
  int nhits[SLBDEPTH][NCHIP][MEMDEPTH];
  int charge_low[SLBDEPTH][NCHIP][MEMDEPTH][NCHANNELS];
  int charge_high[SLBDEPTH][NCHIP][MEMDEPTH][NCHANNELS];
  int gain_hit_low[SLBDEPTH][NCHIP][MEMDEPTH][NCHANNELS];
  int gain_hit_high[SLBDEPTH][NCHIP][MEMDEPTH][NCHANNELS];
  int event;
  //int numbcid;
  int numCol[SLBDEPTH][NCHIP];
  int chipID[SLBDEPTH][NCHIP];
  int slot[SLBDEPTH];
  int slboard_id[SLBDEPTH];
  int n_slboards;
  int acqNumber;


  //  InfoChip * info;
};

//******************************************************************************************************************

void SLBdecoded2ROOT::Initialisation() {
  fout->cd(); R2Rstate=-1;

  tree = new TTree("siwecaldecoded","siwecaldecoded");

  tree->Branch("event",&event,"event/I");
  tree->Branch("acqNumber",&acqNumber,"acqNumber/I");
  tree->Branch("n_slboards",&n_slboards,"n_slboards/I");

  TString name;
  name= TString::Format("slot[%i]/I",SLBDEPTH);
  tree->Branch("slot",slot,name);
  
  name= TString::Format("slboard_id[%i]/I",SLBDEPTH);
  tree->Branch("slboard_id",slboard_id,name);
  
  name= TString::Format("chipid[%i][%i]/I",SLBDEPTH,NCHIP);
  tree->Branch("chipid",chipID,name);

  name= TString::Format("nColumns[%i][%i]/I",SLBDEPTH,NCHIP);
  tree->Branch("nColumns",numCol,name);

  name= TString::Format("bcid[%i][%i][%i]/I",SLBDEPTH,NCHIP,MEMDEPTH);
  tree->Branch("bcid",bcid,name);

  name= TString::Format("corrected_bcid[%i][%i][%i]/I",SLBDEPTH,NCHIP,MEMDEPTH);
  tree->Branch("corrected_bcid",corrected_bcid,name);

  name= TString::Format("badbcid[%i][%i][%i]/I",SLBDEPTH,NCHIP,MEMDEPTH);
  tree->Branch("badbcid",badbcid,name);

  name= TString::Format("nhits[%i][%i][%i]/I",SLBDEPTH,NCHIP,MEMDEPTH);
  tree->Branch("nhits",nhits,name);

  name= TString::Format("lowGain[%i][%i][%i][%i]/I",SLBDEPTH,NCHIP,MEMDEPTH,NCHANNELS);
  tree->Branch("charge_lowGain",charge_low,name);

  name= TString::Format("highGain[%i][%i][%i][%i]/I",SLBDEPTH,NCHIP,MEMDEPTH,NCHANNELS);
  tree->Branch("charge_hiGain",charge_high,name);

  name= TString::Format("gain_hit_low[%i][%i][%i][%i]/I",SLBDEPTH,NCHIP,MEMDEPTH,NCHANNELS);
  tree->Branch("gain_hit_low",gain_hit_low,name);

  name= TString::Format("gain_hit_high[%i][%i][%i][%i]/I",SLBDEPTH,NCHIP,MEMDEPTH,NCHANNELS);
  tree->Branch("gain_hit_high",gain_hit_high,name);

  return;
}






//******************************************************************************************************************

void SLBdecoded2ROOT::treeInit(bool zerosupression=false) { //init data for a single SPILL ?

  for (int isl=0; isl<SLBDEPTH; isl++) {
    for (int k=0; k<NCHIP; k++) {
      for (int i=0; i<MEMDEPTH; i++) {
	bcid[isl][k][i]=-999;
	badbcid[isl][k][i]=-999;
	corrected_bcid[isl][k][i]=-999;
	nhits[isl][k][i]=-999;
	for (int j=0; j<NCHANNELS; j++) {
	  charge_low[isl][k][i][j]=-999;
	  charge_high[isl][k][i][j]=-999;
	  gain_hit_low[isl][k][i][j]=-999;
	  gain_hit_high[isl][k][i][j]=-999;
	}
      }
      
      chipID[isl][k]=-999;
      numCol[isl][k]=0; 
    }
    slot[isl]=-1;
    slboard_id[isl]=-1;
  }
  n_slboards=-1;
  acqNumber=-1;

  return;
}


//******************************************************************************************************************


//******************************************************************************************************************

void SLBdecoded2ROOT::ReadFile(TString inputFileName, bool overwrite, TString outFileName, bool zerosupression) {

  event=0;
  acqNumber=0;
  
  if(outFileName == "default"){
    outFileName = TString::Format("%s.root",inputFileName.Data());
    cout<<outFileName<<endl;
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
  
  Initialisation();
  ifstream fin;
  unsigned short int dataResult=0;

  std::ifstream reading_file(inputFileName);
 
  if(!reading_file){
    cout<<" ERROR  ----------------------------------- No file: "<<inputFileName<<endl;
    return 0;
  } else {
    cout<<" Read File "<<inputFileName<<endl;
    //cout<<" slb= " <<slboard_index<<endl;
  }


  /*
  TH1F *bcid_diff[16];
  TH2F *bcid_correl[16];
  TH1F *bcid_4coinc[16];
  TH1F *bcid_5coinc[16];

  for(int ichip=0; ichip<16; ichip++) {
    bcid_diff[ichip]=new TH1F(TString::Format("bcid_diff_chip%i",ichip),TString::Format("bcid_diff_chip%i",ichip),4001,-2000.5,2000.5);
    bcid_correl[ichip]=new TH2F(TString::Format("bcid_correl_chip%i",ichip),TString::Format("bcid_correl_chip%i",ichip),4096,-0.5,4095.5,4001,-2000.5,2000.5);
  }
  */

  
  int cycleID=-1;
  std::string strheader;
  std::getline(reading_file, strheader);
  std::getline(reading_file, strheader);
  TString tmpst;
  int tmp_n_slboards=0;
  reading_file >> tmpst >> tmpst >> tmpst  >> tmpst >>  tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst  >> tmpst >> tmp_n_slboards >> tmpst;
  cout<<" NB OF CONNECTED SLABs = " << tmp_n_slboards <<endl;

  for(int islboard=0; islboard<tmp_n_slboards; islboard++) 
    std::getline(reading_file, strheader);

  std::getline(reading_file, strheader);
  std::getline(reading_file, strheader);
  std::getline(reading_file, strheader);
  std::getline(reading_file, strheader);
  std::getline(reading_file, strheader);

  std::string str;

  // int bcid_cycle[20][16][15];

  int mapping_z[15];
  //2020 03 12
  mapping_z[0]=8;
  mapping_z[1]=9;
  mapping_z[2]=10;
  mapping_z[3]=5;
  mapping_z[4]=1;
  mapping_z[5]=13;
  mapping_z[6]=11;
  mapping_z[7]=7;

  int mapping_slot[15];
  //2020 03 12
  mapping_slot[0]=1;
  mapping_slot[1]=2;
  mapping_slot[2]=3;
  mapping_slot[3]=4;
  mapping_slot[4]=5;
  mapping_slot[5]=6;
  mapping_slot[6]=7;
  mapping_slot[7]=8;

  
  while (reading_file) {
    // output the line
    TString tmpst1;
    TString tmpst;
    int size=0;
    int chip=-1;
    int slabidx=-1;
    int slabadd=-1;
    reading_file >> tmpst1 >> tmpst >> size >> tmpst >> chip >> tmpst >>  tmpst >> tmpst >> slabidx >> tmpst >> slabadd >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst  >> tmpst >> tmpst >> cycleID >> tmpst >> tmpst >> tmpst >> tmpst;
    
    //    if(slabidx>20)
    std::cout << " -----" <<tmpst1<<" "<<size<< " "<<chip<<" "<<slabidx<<" "<<slabadd<<" "<<cycleID <<" "<<acqNumber<< std::endl;
    if(slabidx<0) break;
    if(acqNumber==0) treeInit(zerosupression);
     if(acqNumber>0 && acqNumber!=cycleID) {
      GetBadBCID();
      tree->Fill();
      treeInit(zerosupression);
    }

    

    //save variables
    acqNumber=cycleID;
    n_slboards=tmp_n_slboards;
    event++;
    int previousBCID=-1000;
    int loopBCID=0;

    slot[mapping_z[slabidx]]=mapping_slot[slabidx];
    slboard_id[mapping_z[slabidx]]=slabadd;
    
    for(int i=0; i<size; i++) {
      //##0 BCID 1627 SCA 1 #Hits 49
      int bcid_tmp=-1;
      int sca=-1;
      int nhits_tmp=-1;
      reading_file >> tmpst >> tmpst >> bcid_tmp >> tmpst >> sca >>  tmpst  >>  nhits_tmp ;
      //      bcid_cycle[slabid][chip][sca]=bcid_tmp;
      // cout<<bcid_tmp<<endl;
      //      if(slabid==slboard_index) {
      sca=size-(sca+1);
      bcid[mapping_z[slabidx]][chip][sca]=bcid_tmp;
      nhits[mapping_z[slabidx]][chip][sca]=nhits_tmp;
      numCol[mapping_z[slabidx]][chip]++;
      
      if(bcid[mapping_z[slabidx]][chip][sca] > 0 && bcid[mapping_z[slabidx]][chip][sca]-previousBCID < 0) loopBCID++;
      if(bcid[mapping_z[slabidx]][chip][sca] > 0 ) corrected_bcid[mapping_z[slabidx]][chip][sca] = bcid[mapping_z[slabidx]][chip][sca]+loopBCID*4096;
      previousBCID=bcid_tmp;
      
      if(chip>-1 && chip<16) {
	chipID[mapping_z[slabidx]][chip]=chip;
      } else {
	cout<<"Wrong chipID = "<<chip<<endl;
	break;
      }
      //      }
      int nchn=0;
      if(zerosupression==false) nchn=NCHANNELS;
      else nchn = nhits_tmp;
      //Ch 0 LG 322 1 0 HG 408 1 0
      for(int ichn=0; ichn<nchn; ichn++) {
	int chn=-1;
	int lg=-1;
	int hg=-1;
	int lg_bit=-1;
	int hg_bit=-1;
	reading_file >> tmpst >> chn >> tmpst >> lg >>  lg_bit >> tmpst >>  tmpst >> hg >> hg_bit >> tmpst;
	//	cout << chn << " "<< hg <<" " <<hg_bit<<endl;
	//	if(slabid==slboard_index) {
	  charge_low[mapping_z[slabidx]][chip][sca][chn]=lg;
	  gain_hit_low[mapping_z[slabidx]][chip][sca][chn]=lg_bit;
	  charge_high[mapping_z[slabidx]][chip][sca][chn]=hg;
	  gain_hit_high[mapping_z[slabidx]][chip][sca][chn]=hg_bit;
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
  for (int k=0; k<NCHIP; k++) {
    //only for valid chips in this spill
    if (chipID[i][k]>=0) {
      for (int ibc=0; ibc<numCol[i][k]; ibc++) {

	// if sca+1 is filled with consec bcid, but sca+2 not, then badbcid[sca]==1 && badbcid[sca+1]==2 (bcid+1 issue, events are not bad, just the next sca trigger is empty)
	// if sca+1 is filled with consec bcid, and sca+2 also, then badbcid[sca]==3 && badbcid[sca+1]==3 (retriggering)
	// if sca+1 is not filled with consec bcid,  badbcid==0

	if(ibc==0) {
	  badbcid[i][k][ibc]=0;
	  int corr_bcid=corrected_bcid[i][k][ibc];
	  int corr_bcid1=0;
	  int corr_bcid2=0;

	  if(corrected_bcid[i][k][ibc+1]>0 && corrected_bcid[i][k][ibc]>0 && (corrected_bcid[i][k][ibc+1]-corrected_bcid[i][k][ibc])>0) 
	    corr_bcid1=corrected_bcid[i][k][ibc+1];

	  if(corrected_bcid[i][k][ibc+2]>0 && (corrected_bcid[i][k][ibc+2]-corrected_bcid[i][k][ibc+1])>0) 
	    corr_bcid2=corrected_bcid[i][k][ibc+2];

	  if(corr_bcid2>0) {
	    //empty events
	    if( ( corr_bcid2-corr_bcid1) >(BCIDTHRES - 1) && (corr_bcid1-corr_bcid) ==1) {
	      badbcid[i][k][ibc]=1;
	      badbcid[i][k][ibc+1]=2;
	    }
	    // pure retriggers
	    if( ( corr_bcid2-corr_bcid1) < BCIDTHRES && (corr_bcid1-corr_bcid) < BCIDTHRES) {
	      badbcid[i][k][ibc]=3;
	      badbcid[i][k][ibc+1]=3;
	      badbcid[i][k][ibc+2]=3;
	    }
    	  }

	  if( corr_bcid1 > 0 && (corr_bcid1-corr_bcid) > 1 && (corr_bcid1-corr_bcid) <BCIDTHRES) {
	    badbcid[i][k][ibc]=3;
	    badbcid[i][k][ibc+1]=3;
	  }  
	} //ibc==0 if 

	if(ibc>0 && badbcid[i][k][ibc]<0 && corrected_bcid[i][k][ibc] >0 &&  (corrected_bcid[i][k][ibc]-corrected_bcid[i][k][ibc-1])>0 ) {
	  badbcid[i][k][ibc]=0;
	  int corr_bcid=corrected_bcid[i][k][ibc];
	  int corr_bcidminus=corrected_bcid[i][k][ibc-1];

	  if(corrected_bcid[i][k][ibc+1]>0 && (corrected_bcid[i][k][ibc+1]-corrected_bcid[i][k][ibc])>0) {
	    int corr_bcid1=corrected_bcid[i][k][ibc+1];

	    if(corrected_bcid[i][k][ibc+2]>0 && (corrected_bcid[i][k][ibc+2]-corrected_bcid[i][k][ibc+1])>0) {
	      int corr_bcid2=corrected_bcid[i][k][ibc+2];
	      if( ( corr_bcid2-corr_bcid1) < BCIDTHRES && (corr_bcid1-corr_bcid) < BCIDTHRES) {
		badbcid[i][k][ibc]=3;
		badbcid[i][k][ibc+1]=3;
		badbcid[i][k][ibc+2]=3;
	      }
	      if( (corr_bcid1-corr_bcid) < BCIDTHRES && (corr_bcid-corr_bcidminus) < BCIDTHRES ) badbcid[i][k][ibc]=3;

	      if( badbcid[i][k][ibc]!=3 && ( corr_bcid2-corr_bcid1) >(BCIDTHRES - 1) && (corr_bcid1-corr_bcid) ==1) {
		badbcid[i][k][ibc]=1;
		badbcid[i][k][ibc+1]=2;
	      }
	      if( badbcid[i][k][ibc]!=3 && ( corr_bcid2-corr_bcid1) >(BCIDTHRES - 1) && (corr_bcid1-corr_bcid) > 1 && (corr_bcid1-corr_bcid) <BCIDTHRES) {
		badbcid[i][k][ibc]=3;
		badbcid[i][k][ibc+1]=3;
	      }
	      if( (corr_bcid-corr_bcidminus) < BCIDTHRES ) badbcid[i][k][ibc]=3;

	      //if( badbcid[i][k][ibc-1]==1 && (corr_bcid1-corr_bcid) > (BCIDTHRES - 1)) badbcid[i][k][ibc]=2;
	    } else {
	      if( (corr_bcid1-corr_bcid) < BCIDTHRES ) badbcid[i][k][ibc]=3;
	      if( (corr_bcid-corr_bcidminus) < BCIDTHRES ) badbcid[i][k][ibc]=3;
	    } //ibc+2 if
	  } else {
	    if( (corr_bcid-corr_bcidminus) < BCIDTHRES ) badbcid[i][k][ibc]=3;
	  }//ibc+1 if
	} //ibc>0 if 

      
	//tag zero (under/over flow) data
	//	count_negdata=0;
	
	//	for (int ichan=0; ichan<NCHANNELS; ichan++) {
	//  if  (charge_high[i][k][ibc][ichan] < NEGDATA_THR) count_negdata++;
	//}//ichan

	//if (count_negdata>0) {badbcid[i][k][ibc]+=32;}

      }//ibc

    }//chipID
  }//k
  }//i   slboard
  
}

//******************************************************************************************************************

////////////////////////  GETTREE   ////////////////////////
int SLBdecoded2ROOT::GetTree (TString rootfilename) { //from raw2root 1st pass
  //open pre-precessed data
  finroot = new TFile(rootfilename,"read");
  if (! finroot ) {return 0;}
  slboardread = (TTree*)finroot->Get("siwecaldecoded"); //rawdata

  //get data from file
  slboardread->SetBranchAddress("gain_hit_high",gain_hit_high );
  slboardread->SetBranchAddress("gain_hit_low",gain_hit_low );
  slboardread->SetBranchAddress("charge_hiGain",charge_high );
  slboardread->SetBranchAddress("charge_lowGain",charge_low );
  slboardread->SetBranchAddress("bcid", bcid);
  slboardread->SetBranchAddress("badbcid", badbcid);
  slboardread->SetBranchAddress("nhits", nhits);
  slboardread->SetBranchAddress("event", &event);
  slboardread->SetBranchAddress("acqNumber", &acqNumber);
  slboardread->SetBranchAddress("corrected_bcid",corrected_bcid);
  slboardread->SetBranchAddress("nColumns",numCol);
  slboardread->SetBranchAddress("chipid",chipID);
  slboardread->SetBranchAddress("slot",slot);
  slboardread->SetBranchAddress("slboard_id",slboard_id);
  slboardread->SetBranchAddress("n_slboards",&n_slboards);
  R2Rstate=1;
  return 1;

}//method GetTree


#endif

