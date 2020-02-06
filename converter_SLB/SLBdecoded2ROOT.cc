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
#include <fstream>
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
  
  void ReadFile(TString inputFileName, bool overwrite=false, TString outFileName = "default", int slboard_index=2, int nslboards=4,bool zerosupression=false);

protected:

  enum {
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

  int bcid[NCHIP][MEMDEPTH];
  int corrected_bcid[NCHIP][MEMDEPTH];
  int badbcid[NCHIP][MEMDEPTH];
  int nhits[NCHIP][MEMDEPTH];
  int charge_low[NCHIP][MEMDEPTH][NCHANNELS];
  int charge_high[NCHIP][MEMDEPTH][NCHANNELS];
  int gain_hit_low[NCHIP][MEMDEPTH][NCHANNELS];
  int gain_hit_high[NCHIP][MEMDEPTH][NCHANNELS];
  int event;
  //int numbcid;
  int numCol[NCHIP];
  int chipID[NCHIP];
  int acqNumber;


  //  InfoChip * info;
};

//******************************************************************************************************************

void SLBdecoded2ROOT::Initialisation() {
  fout->cd(); R2Rstate=-1;

  tree = new TTree("slboard","slboard");

  tree->Branch("event",&event,"event/I");
  tree->Branch("acqNumber",&acqNumber,"acqNumber/I");

  TString name;

  name= TString::Format("chipid[%i]/I",NCHIP);
  tree->Branch("chipid",chipID,name);

  name= TString::Format("nColumns[%i]/I",NCHIP);
  tree->Branch("nColumns",numCol,name);

  name= TString::Format("bcid[%i][%i]/I",NCHIP,MEMDEPTH);
  tree->Branch("bcid",bcid,name);

  name= TString::Format("corrected_bcid[%i][%i]/I",NCHIP,MEMDEPTH);
  tree->Branch("corrected_bcid",corrected_bcid,name);

  name= TString::Format("badbcid[%i][%i]/I",NCHIP,MEMDEPTH);
  tree->Branch("badbcid",badbcid,name);

  name= TString::Format("nhits[%i][%i]/I",NCHIP,MEMDEPTH);
  tree->Branch("nhits",nhits,name);

  name= TString::Format("lowGain[%i][%i][%i]/I",NCHIP,MEMDEPTH,NCHANNELS);
  tree->Branch("charge_lowGain",charge_low,name);

  name= TString::Format("highGain[%i][%i][%i]/I",NCHIP,MEMDEPTH,NCHANNELS);
  tree->Branch("charge_hiGain",charge_high,name);

  name= TString::Format("gain_hit_low[%i][%i][%i]/I",NCHIP,MEMDEPTH,NCHANNELS);
  tree->Branch("gain_hit_low",gain_hit_low,name);

  name= TString::Format("gain_hit_high[%i][%i][%i]/I",NCHIP,MEMDEPTH,NCHANNELS);
  tree->Branch("gain_hit_high",gain_hit_high,name);

  return;
}






//******************************************************************************************************************

void SLBdecoded2ROOT::treeInit(bool zerosupression=false) { //init data for a single SPILL ?

  if(zerosupression == true) {
    for (int k=0; k<NCHIP; k++) {
      for (int i=0; i<MEMDEPTH; i++) {
	bcid[k][i]=1;
	badbcid[k][i]=1;
	corrected_bcid[k][i]=1;
	nhits[k][i]=1;
	for (int j=0; j<NCHANNELS; j++) {
	  charge_low[k][i][j]=1;
	  charge_high[k][i][j]=1;
	  gain_hit_low[k][i][j]=0;
	  gain_hit_high[k][i][j]=0;
	}
      }
      chipID[k]=-1;
      numCol[k]=1;
    }
  } else {
      for (int k=0; k<NCHIP; k++) {
      for (int i=0; i<MEMDEPTH; i++) {
	bcid[k][i]=-999;
	badbcid[k][i]=-999;
	corrected_bcid[k][i]=-999;
	nhits[k][i]=-999;
	for (int j=0; j<NCHANNELS; j++) {
	  charge_low[k][i][j]=-999;
	  charge_high[k][i][j]=-999;
	  gain_hit_low[k][i][j]=-999;
	  gain_hit_high[k][i][j]=-999;
	}
      }
      
      chipID[k]=-999;
      numCol[k]=0; 
    }
  }


  return;
}


//******************************************************************************************************************


//******************************************************************************************************************

void SLBdecoded2ROOT::ReadFile(TString inputFileName, bool overwrite, TString outFileName, int slboard_index, int nslboards, bool zerosupression) {

  event=0;
  acqNumber=0;
  
  if(outFileName == "default"){
    outFileName = TString::Format("%s_SLB_%i.root",inputFileName.Data(),slboard_index);
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
    cout<<" slb= " <<slboard_index<<endl;
  }


  TH1F *bcid_diff[16];
  TH2F *bcid_correl[16];
  TH1F *bcid_4coinc[16];
  TH1F *bcid_5coinc[16];

  for(int ichip=0; ichip<16; ichip++) {
    bcid_diff[ichip]=new TH1F(TString::Format("bcid_diff_chip%i",ichip),TString::Format("bcid_diff_chip%i",ichip),4001,-2000.5,2000.5);
    bcid_correl[ichip]=new TH2F(TString::Format("bcid_correl_chip%i",ichip),TString::Format("bcid_correl_chip%i",ichip),4096,-0.5,4095.5,4001,-2000.5,2000.5);
  }


  
  int cycleID=-1;
  std::string strheader;
  std::getline(reading_file, strheader);
  std::getline(reading_file, strheader);
  std::getline(reading_file, strheader);
  std::getline(reading_file, strheader);
  std::getline(reading_file, strheader);
  std::getline(reading_file, strheader);
  //std::getline(reading_file, strheader);
  std::getline(reading_file, strheader);
  if(nslboards==1) {
    std::getline(reading_file, strheader);
  }
  if(nslboards==2){
    std::getline(reading_file, strheader);
    std::getline(reading_file, strheader);
  }
  if(nslboards==3){
    std::getline(reading_file, strheader);
    std::getline(reading_file, strheader);
    std::getline(reading_file, strheader);
  }

  std::string str;

  int bcid_cycle[20][16][15];
  
  while (reading_file) {
    // output the line
    TString tmpst;
    int size=0;
    int chip=-1;
    int slabid=-1;
    reading_file >> tmpst >> tmpst >> size >> tmpst >> chip >> tmpst >>  tmpst >> tmpst >> tmpst >> tmpst >> slabid >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst  >> tmpst >> tmpst >> cycleID >> tmpst >> tmpst >> tmpst >> tmpst;

    if(acqNumber!=cycleID && acqNumber>0 ) {
	
      for(int i=0; i<16;i++) 
	for(int k=0;k<20; k++) 
	  for(int j=0;j<15; j++) 
	    for(int k2=k+1;k2<20; k2++) 
	      for(int j2=0;j2<15; j2++) 
		if(bcid_cycle[k][i][j]>0 && bcid_cycle[k2][i][j2]>0 ) {
		  bcid_diff[i]->Fill(bcid_cycle[k][i][j]-bcid_cycle[k2][i][j2]);
		  bcid_correl[i]->Fill(bcid_cycle[k][i][j],bcid_cycle[k][i][j]-bcid_cycle[k2][i][j2]);
		}
		
    }
  
    if(acqNumber!=cycleID) {
      for(int k=0; k<20;k++)
	for(int i=0; i<16;i++)
	  for(int j=0; j<15;j++) bcid_cycle[k][i][j]=-1;
    }
    //std::cout << " -----" <<size<< " "<<chip<<" "<<slabid<<" "<<cycleID <<" "<<acqNumber<< std::endl;

    if(acqNumber==0) treeInit(zerosupression);
     if(acqNumber>0 && acqNumber!=cycleID) {
      GetBadBCID();
      // if(slabid==slboard_index) {
	tree->Fill();
	//	std::cout << " Filling trees, cycleId=" <<cycleID << std::endl;

	//	BuildTimeEvents(outFileName);
	//  }
      treeInit(zerosupression);
    }


    //save variables
    acqNumber=cycleID;
    event++;
    int previousBCID=-1000;
    int loopBCID=0;
    
    for(int i=0; i<size; i++) {
      //##0 BCID 1627 SCA 1 #Hits 49
      int bcid_tmp=-1;
      int sca=-1;
      int nhits_tmp=-1;
      reading_file >> tmpst >> tmpst >> bcid_tmp >> tmpst >> sca >>  tmpst  >>  nhits_tmp ;
      bcid_cycle[slabid][chip][sca]=bcid_tmp;
      // cout<<bcid_tmp<<endl;
      if(slabid==slboard_index) {
	sca=size-(sca+1);
	bcid[chip][sca]=bcid_tmp;
	nhits[chip][sca]=nhits_tmp;
	numCol[chip]++;
      
	if(bcid[chip][sca] > 0 && bcid[chip][sca]-previousBCID < 0) loopBCID++;
	if(bcid[chip][sca] > 0 ) corrected_bcid[chip][sca] = bcid[chip][sca]+loopBCID*4096;
	previousBCID=bcid_tmp;
	
	if(chip>-1 && chip<16) {
	  chipID[chip]=chip;
	} else {
	  cout<<"Wrong chipID = "<<chip<<endl;
	  break;
	}
      }
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
	if(slabid==slboard_index) {
	  charge_low[chip][sca][chn]=lg;
	  gain_hit_low[chip][sca][chn]=lg_bit;
	  charge_high[chip][sca][chn]=hg;
	  gain_hit_high[chip][sca][chn]=hg_bit;
	}
      }//end channel  
    }//end events of the chi

  }
  
  fout->cd(); 
  fout->Write(0);
  for(int ichip=0; ichip<16; ichip++) {
    bcid_diff[ichip]->Write();
    bcid_correl[ichip]->Write();
  }
  fout->Close();
    
  return;
}

void SLBdecoded2ROOT::GetBadBCID() {

  //add tags
  int  count_negdata=0;

  for (int k=0; k<NCHIP; k++) {
    //only for valid chips in this spill
    if (chipID[k]>=0) {
      for (int ibc=0; ibc<numCol[k]; ibc++) {

	// if sca+1 is filled with consec bcid, but sca+2 not, then badbcid[sca]==1 && badbcid[sca+1]==2 (bcid+1 issue, events are not bad, just the next sca trigger is empty)
	// if sca+1 is filled with consec bcid, and sca+2 also, then badbcid[sca]==3 && badbcid[sca+1]==3 (retriggering)
	// if sca+1 is not filled with consec bcid,  badbcid==0

	if(ibc==0) {
	  badbcid[k][ibc]=0;
	  int corr_bcid=corrected_bcid[k][ibc];
	  int corr_bcid1=0;
	  int corr_bcid2=0;

	  if(corrected_bcid[k][ibc+1]>0 && corrected_bcid[k][ibc]>0 && (corrected_bcid[k][ibc+1]-corrected_bcid[k][ibc])>0) 
	    corr_bcid1=corrected_bcid[k][ibc+1];

	  if(corrected_bcid[k][ibc+2]>0 && (corrected_bcid[k][ibc+2]-corrected_bcid[k][ibc+1])>0) 
	    corr_bcid2=corrected_bcid[k][ibc+2];

	  if(corr_bcid2>0) {
	    //empty events
	    if( ( corr_bcid2-corr_bcid1) >(BCIDTHRES - 1) && (corr_bcid1-corr_bcid) ==1) {
	      badbcid[k][ibc]=1;
	      badbcid[k][ibc+1]=2;
	    }
	    // pure retriggers
	    if( ( corr_bcid2-corr_bcid1) < BCIDTHRES && (corr_bcid1-corr_bcid) < BCIDTHRES) {
	      badbcid[k][ibc]=3;
	      badbcid[k][ibc+1]=3;
	      badbcid[k][ibc+2]=3;
	    }
    	  }

	  if( corr_bcid1 > 0 && (corr_bcid1-corr_bcid) > 1 && (corr_bcid1-corr_bcid) <BCIDTHRES) {
	    badbcid[k][ibc]=3;
	    badbcid[k][ibc+1]=3;
	  }  
	} //ibc==0 if 

	if(ibc>0 && badbcid[k][ibc]<0 && corrected_bcid[k][ibc] >0 &&  (corrected_bcid[k][ibc]-corrected_bcid[k][ibc-1])>0 ) {
	  badbcid[k][ibc]=0;
	  int corr_bcid=corrected_bcid[k][ibc];
	  int corr_bcidminus=corrected_bcid[k][ibc-1];

	  if(corrected_bcid[k][ibc+1]>0 && (corrected_bcid[k][ibc+1]-corrected_bcid[k][ibc])>0) {
	    int corr_bcid1=corrected_bcid[k][ibc+1];

	    if(corrected_bcid[k][ibc+2]>0 && (corrected_bcid[k][ibc+2]-corrected_bcid[k][ibc+1])>0) {
	      int corr_bcid2=corrected_bcid[k][ibc+2];
	      if( ( corr_bcid2-corr_bcid1) < BCIDTHRES && (corr_bcid1-corr_bcid) < BCIDTHRES) {
		badbcid[k][ibc]=3;
		badbcid[k][ibc+1]=3;
		badbcid[k][ibc+2]=3;
	      }
	      if( (corr_bcid1-corr_bcid) < BCIDTHRES && (corr_bcid-corr_bcidminus) < BCIDTHRES ) badbcid[k][ibc]=3;

	      if( badbcid[k][ibc]!=3 && ( corr_bcid2-corr_bcid1) >(BCIDTHRES - 1) && (corr_bcid1-corr_bcid) ==1) {
		badbcid[k][ibc]=1;
		badbcid[k][ibc+1]=2;
	      }
	      if( badbcid[k][ibc]!=3 && ( corr_bcid2-corr_bcid1) >(BCIDTHRES - 1) && (corr_bcid1-corr_bcid) > 1 && (corr_bcid1-corr_bcid) <BCIDTHRES) {
		badbcid[k][ibc]=3;
		badbcid[k][ibc+1]=3;
	      }
	      if( (corr_bcid-corr_bcidminus) < BCIDTHRES ) badbcid[k][ibc]=3;

	      //if( badbcid[k][ibc-1]==1 && (corr_bcid1-corr_bcid) > (BCIDTHRES - 1)) badbcid[k][ibc]=2;
	    } else {
	      if( (corr_bcid1-corr_bcid) < BCIDTHRES ) badbcid[k][ibc]=3;
	      if( (corr_bcid-corr_bcidminus) < BCIDTHRES ) badbcid[k][ibc]=3;
	    } //ibc+2 if
	  } else {
	    if( (corr_bcid-corr_bcidminus) < BCIDTHRES ) badbcid[k][ibc]=3;
	  }//ibc+1 if
	} //ibc>0 if 

      
	//tag zero (under/over flow) data
	count_negdata=0;
	
	for (int ichan=0; ichan<NCHANNELS; ichan++) {
	  if  (charge_high[k][ibc][ichan] < NEGDATA_THR) count_negdata++;
	}//ichan

	if (count_negdata>0) {badbcid[k][ibc]+=32;}

      }//ibc

    }//chipID
  }//k
  
}

//******************************************************************************************************************

////////////////////////  GETTREE   ////////////////////////
int SLBdecoded2ROOT::GetTree (TString rootfilename) { //from raw2root 1st pass
  //open pre-precessed data
  finroot = new TFile(rootfilename,"read");
  if (! finroot ) {return 0;}
  slboardread = (TTree*)finroot->Get("slboard"); //rawdata

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


  R2Rstate=1;
  return 1;

}//method GetTree


#endif

