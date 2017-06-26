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
#include "InfoChip.C"
#include "TChain.h"




class PedestalObj{
public:
  PedestalObj(){
    info = new InfoChip();
    mDebug=false;
  }
  ~PedestalObj(){}

  void SetDebug(){mDebug=true;}

  void Execute(TString pathName /*, int dif */, TString outFileName="results/pedestals/MIPPedestal.root"){
    TString filename = pathName;

    filename += "/*_merge.root";

    TFile *fout = new TFile(outFileName, "recreate");

    for(int dif=0;dif<NSLABS;dif++){
      for(int i = 0; i < NCHIP; i++){
	TString name = "Pedestal_dif";
	name+=dif;
	name+="_chip";
	name += i+1;
	dPedestal[dif][i] = fout->mkdir(name);
      }
    }


    /*    for(int dif=0;dif<NSLABS;dif++){
	  for(int i = 0; i<NCHIP; i++){
	  fout->cd();
	  TString name = "Pedestal_Mean_dif";
	  name+=dif;
	  name+="_Chip_";
	  name += i+1;
	  TString title = "Mean of the Pedestal for column 0 of chip ";
	  title += i+1;
	  hPedestal_Mean[dif][i] =  new TH1D(name,title,64,0,64);
	  hPedestal_Mean[dif][i]->GetXaxis()->SetTitle("Channels");
	  hPedestal_Mean[dif][i]->GetYaxis()->SetTitle("Mean (uADC)");

	  name = "Pedestal_RMS_dif";
	  name+=dif;
	  name+="_Chip_";
	  name += i+1;
	  title = "RMS of the Pedestal for column 0 of chip ";
	  title += i+1;
	  hPedestal_RMS[dif][i] =  new TH1D(name,title,64,0,64);
	  hPedestal_RMS[dif][i]->GetXaxis()->SetTitle("Channels");
	  hPedestal_RMS[dif][i]->GetYaxis()->SetTitle("RMS (uADC)");

	  dPedestal[dif][i]->cd();
	  for(int k = 0; k<NCHANNELS; k++){
	  for(int j = 0; j < MEMDEPTH; j++){
	  name = "Dif";
	  name+=dif;
	  name += "_Chan_";
	  name += k;
	  name += "_Col_";
	  name += j;
	  hPedestal[dif][i][k][j] = new TH1D(name,name,4096,0,4096);
	  hPedestal[dif][i][k][j]->GetXaxis()->SetTitle("Pedestal");
	  hPedestal[dif][i][k][j]->GetYaxis()->SetTitle("#");
	  }
	  }
	  }
	  }
    */

    fout->cd();
    TH2D * htime1[NCHANNELS];
    TH2D * htime2[NCHANNELS];
    TH2D * htime3[NCHANNELS];
    TH2D * htime4[NCHANNELS];

    for(int i = 0; i < NCHANNELS; i++){
      TString name;
      name = "chip1_Chan_";
      name += i;
      htime1[i] = new TH2D(name,"htime",500,0,4000000,90,260,350);
      name = "chip2_Chan_";
      name += i;
      htime2[i] = new TH2D(name,"htime",500,0,4000000,90,260,350);
      name = "chip3_Chan_";
      name += i;
      htime3[i] = new TH2D(name,"htime",500,0,4000000,90,260,350);
      name = "chip4_Chan_";
      name += i;
      htime4[i] = new TH2D(name,"htime",500,0,4000000,90,260,350);
    }

   // TH2D * hfiles = new TH2D("hfiles","hfiles",500,0,4000000,30,0,30);

    /*    TH2 *h2Mean[MEMDEPTH];
	  TH2 *h2RMS[MEMDEPTH];

	  TH1D * hMeanPedestals[MEMDEPTH]; 

	  for(int j = 0; j < MEMDEPTH; j++){
	  TString name = "h2Mean_";
	  name += j;
	  h2Mean[j] = new TH2D(name,"Mean of the Pedestal for 1 column",18,0.5,18.5,18,0.5,18.5);
	  h2Mean[j]->GetXaxis()->SetTitle("X");
	  h2Mean[j]->GetYaxis()->SetTitle("Y");

	  name = "h2RMS_";
	  name += j;
	  h2RMS[j] = new TH2D(name,"RMS of the Pedestal for 1 column",18,0.5,18.5,18,0.5,18.5);
	  h2RMS[j]->GetXaxis()->SetTitle("X");
	  h2RMS[j]->GetYaxis()->SetTitle("Y");

	  name = "hMeanPedAllChannels_";
	  name += j;
	  hMeanPedestals[j] = new TH1D(name,"Mean of the Pedestals",200,200,400);
	  hMeanPedestals[j]->GetXaxis()->SetTitle("Mean value (uADC)");
	  }
    */


    //   TFile *f;
    //    TTree *tree=0 ;

    TChain * tree;
    tree = new TChain("fev10");
    tree->Add(filename);

    //f = new TFile(filename,"read");
    //if(f->IsOpen()){

    //tree = (TTree*)f->Get("fev8");


    tree->SetBranchAddress("chipid", chipID_in);
    tree->SetBranchAddress("acqNumber", &acqNumber_in);

    tree->SetBranchAddress("bcid"           , bcid_in);
    tree->SetBranchAddress("nhits"           , nhits_in);
    tree->SetBranchAddress("badbcid"        , badbcid_in);
    tree->SetBranchAddress("nColumns"       , numCol_in);
    tree->SetBranchAddress("gain_hit_high"  , gain_hit_high_in );
    tree->SetBranchAddress("charge_hiGain"  , charge_high_in );
    tree->SetBranchAddress("gain_hit_low"   , gain_hit_low_in );
    tree->SetBranchAddress("charge_lowGain" , charge_low_in );

    //}






    double zpos[NSLABS] = {30.3,60.3};//0.3,30.3,60.3,90.3,120.3,135.3};
    for(int iz = 0;iz< NSLABS;iz++) zpos[iz] -= 0.3/2. ;

    int MAXBCID = MEMDEPTH*NCHIP*NSLABS;

    int pointgr = 0;

    int listOfBCID[MAXBCID];
    int listOfbadBCID[MAXBCID];
    int orderedListOfBCID[MAXBCID];

    int listOfIsSuspect[MAXBCID];

    for(unsigned i = 0 ;i<tree->GetEntries();i++){

      tree->GetEntry(i);

      bool empty = true;
      for(int slab = 0;slab<NSLABS;slab++){
	for (int chip=0; chip<NCHIP; chip++) {
	  if(numCol_in[slab][chip]>0) empty = false;
	}
      }
      if(!empty){

	int c=-1;

	for(int tt = 0;tt<MAXBCID;tt++){
	  listOfBCID[tt]=-1;
	  orderedListOfBCID[tt]=-1;
	  listOfIsSuspect[tt]=0;
	  listOfbadBCID[tt]=-1;
	}
	int nBCID = 0;
	int nbadBCID = 0;

	for(int slab = 0;slab<NSLABS;slab++){
	  for (int chip=0; chip<NCHIP; chip++) {
	    for (int j=0; j<MEMDEPTH; j++) {

	      int currentBCID = bcid_in[slab][chip][j];
	      if(currentBCID>=0){
		if(badbcid_in[slab][chip][j]){
		  listOfbadBCID[nbadBCID] = currentBCID;
		  nbadBCID++;
		  continue;
		}

		if(j>0) if(fabs(bcid_in[slab][chip][j]-bcid_in[slab][chip][j-1])<6) continue;

		//bool isbad = false;
		//for(int ibad = 0;ibad<nbadBCID;ibad++){
		//  if(listOfbadBCID[ibad]==currentBCID)isbad=true;
		//}
		//if(isbad) continue;

		bool newBCID = true;
		for(int k = 0;k<nBCID+1;k++){
		  if(listOfBCID[k] == currentBCID) newBCID = false;
		}
		if(newBCID){
		  listOfBCID[nBCID]=currentBCID;

		  nBCID++;
		}

	      }
	    }
	  }
	}


	int min = 1000000;
	int pos_max = 0;
	int pos = 0;

	for(int k1=0;k1<nBCID;k1++){

	  for(int k=0;k<nBCID;k++){
	    if(listOfBCID[k]>=0 && listOfBCID[k]<min){
	      min = listOfBCID[k];
	      pos_max = k;
	    }
	  }
	  if(min < 100000){
	    orderedListOfBCID[pos] = min;
	    listOfBCID[pos_max]=-1;
	    //cout<<orderedListOfBCID[pos]<<"  ";
	    min = 10000000;
	    pos++;
	  }
	}

	int start = 0;

	for(int k1=start;k1<nBCID;k1++){
	  listOfIsSuspect[k1]=0;
	  if(k1>0 && (orderedListOfBCID[k1]-orderedListOfBCID[k1-1]<6)) listOfIsSuspect[k1]=1;
	  if(acqNumber_in==30) cout<<"k1 = "<<k1<<"  "<<orderedListOfBCID[k1]<<"  "<<listOfIsSuspect[k1]<<endl;
	}

	bool next = false;
	eventType = 0;
	//if(nBCID>2)  if(orderedListOfBCID[1]==orderedListOfBCID[0]+1)start=1;
	if(mDebug) cout<<"--------------- event "<<i<<" / "<< tree->GetEntries()   <<" --------------- "<<endl;
	for(int k1=start;k1<nBCID;k1++){
	  next = false;
	  eventType = 0;
	  //suspectBCID = 0;
	  treeInit();

	  //if(k1>0) cout<<k1<<"  "<<orderedListOfBCID[k1]<<"  "<<orderedListOfBCID[k1-1]<<endl;
	  //if(listOfIsSuspect[k1]==1) suspectBCID=1;


	  //if(orderedListOfBCID[k1+1]==orderedListOfBCID[k1]+1) next = true;
	  if(mDebug){
	    if(orderedListOfBCID[k1+2]==orderedListOfBCID[k1]+2) 
	      cout<<"WARNING  :  too many consecutive BCID (+2) : "<< orderedListOfBCID[k1] <<endl;
	    if(orderedListOfBCID[k1+3]==orderedListOfBCID[k1]+3) 
	      cout<<"WARNING  :  too many consecutive BCID (+3) : "<< orderedListOfBCID[k1] <<endl;
	  }

	  //bcid = orderedListOfBCID[k1];
	  pointgr++;

	  for(int slab = 0;slab<NSLABS;slab++){
	    for (int chip=0; chip<NCHIP; chip++) {
	      c = info->GetASUChipNumberFromChipID(chipID_in[slab][chip])-1;
	      for (int j=0; j<MEMDEPTH; j++) {
		if( (orderedListOfBCID[k1]==bcid_in[slab][chip][j] || orderedListOfBCID[k1]==bcid_in[slab][chip][j]-1 )  && !badbcid_in[slab][chip][j]){
                  if( nhits_in[slab][chip][j]>=5  ){//case with an event with high number of hits
		    eventType=2;
		  }
		


		}
	      }
	    }
	  }

	  if(eventType!=2){
	    for(int slab = 0;slab<NSLABS;slab++){
	      for (int chip=0; chip<NCHIP; chip++) {
		c = info->GetASUChipNumberFromChipID(chipID_in[slab][chip])-1;
		for (int j=0; j<MEMDEPTH; j++) {


		  if( (orderedListOfBCID[k1]==bcid_in[slab][chip][j] || orderedListOfBCID[k1]==bcid_in[slab][chip][j]-1 )  && !badbcid_in[slab][chip][j]){
		    if( nhits_in[slab][chip][j]>=5  ){//case with an event with high number of hits
		      eventType=2;
		    }
		    else {
		      if(eventType==0) eventType=1;
		    }

		    if(orderedListOfBCID[k1]!=bcid_in[slab][chip][j]) next = true;
		    //else next = false;

		    if(mDebug){
		      cout<<"slab = "<<slab<<"  bcid = "<<bcid_in[slab][chip][j]<<"  chip = "<< chip  <<"  j = "<< j <<"   nhits = "<<nhits_in[slab][chip][j];//    <<endl;
		      int jj = 0;
		      while(!gain_hit_high_in[slab][chip][j][jj]%2) jj++;
		      cout<<"   hit = "<<jj<<endl;
		    }


		    //if(orderedListOfBCID[memory]!=bcid[slab][chip][j]) cout<<"     + slab "<<slab<<"   bcid = "<<bcid[slab][chip][j]<<endl;
		    for(int k=0;k<NCHANNELS;k++){

		      if ( gain_hit_high_in[slab][chip][j][k]%2 != 1 ){
			if(slab==0 &&  j==0 && c==0) htime1[k]->Fill(pointgr,charge_high_in[slab][chipID_in[slab][chip]][j][k] );
		        if(slab==0 &&  j==0 && c==1) htime2[k]->Fill(pointgr,charge_high_in[slab][chipID_in[slab][chip]][j][k] );
			if(slab==0 &&  j==0 && c==2) htime3[k]->Fill(pointgr,charge_high_in[slab][chipID_in[slab][chip]][j][k] );
			if(slab==0 &&  j==0 && c==3) htime4[k]->Fill(pointgr,charge_high_in[slab][chipID_in[slab][chip]][j][k] );

		      }
		    }
		  }
		}
	      }
	    }
	  }

	  if(mDebug) cout<<"remove bcid = "<< orderedListOfBCID[k1]<<"  "<<k1<<endl;
	  orderedListOfBCID[k1]=-1;
	  if(next){
	    k1++;//=2;
	    if(mDebug) cout<<"remove bcid = "<< orderedListOfBCID[k1]<<"  "<<k1<<endl;
	    orderedListOfBCID[k1]=-1;
	  }
	}
      }
    }

    fout->cd();
    fout->Write(0);
    fout->Close();
  }

  void treeInit() {

  }


protected:
  //TFile *fout;

  enum {
    MEMDEPTH=15,
    NCHANNELS=64,
    NCHIP=16,
    MIN_BCID=5,
    NSLABS=2
  };


  TDirectory *dPedestal[NSLABS][NCHIP];
  TH1D *hPedestal[NSLABS][NCHIP][NCHANNELS][MEMDEPTH];
  TH1D *hPedestal_Mean[NSLABS][NCHIP];
  TH1D *hPedestal_RMS[NSLABS][NCHIP];
  bool nohitnear;

  int bcid_in[NSLABS][NCHIP][MEMDEPTH];
  int badbcid_in[NSLABS][NCHIP][MEMDEPTH];
  int charge_low_in[NSLABS][NCHIP][MEMDEPTH][NCHANNELS];
  int charge_high_in[NSLABS][NCHIP][MEMDEPTH][NCHANNELS];
  int gain_hit_low_in[NSLABS][NCHIP][MEMDEPTH][NCHANNELS];
  int gain_hit_high_in[NSLABS][NCHIP][MEMDEPTH][NCHANNELS];
  int numCol_in[NSLABS][NCHIP];
  int chipID_in[NSLABS][NCHIP];
  int acqNumber_in;
  int corrected_bcid_in[NSLABS][NCHIP][MEMDEPTH];
  int nhits_in[NSLABS][NCHIP][MEMDEPTH];

  bool mDebug;

  int eventType;

  InfoChip * info;

};
