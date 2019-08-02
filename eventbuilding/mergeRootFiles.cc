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

class mergeRootFiles{
public:
  mergeRootFiles(){

  }
  ~mergeRootFiles(){}

  void Merge(TString run, TString indir1, TString indir2, TString outdir, TString mode){


    TString outfilename= outdir + "/run_" + run + "_merge.root";
    cout << "Merging into output file: " << outfilename << endl;
    
    TFile *f[NSLABS];
    TTree *tree[NSLABS] ;


    int entries[NSLABS];
    int maxEntries = 0;
    int entry[NSLABS];

    for(int i = 0;i<NSLABS;i++){
      entries[i] = 0;
      entry[i]=0;

      TString filename="";
     
      if(i == 0) filename = indir1+"/run_"+run + "_dif_1_1_1.raw.root";
      if(i == 1) filename = indir1+"/run_"+run + "_dif_1_1_2.raw.root";
      if(i == 2) filename = indir1+"/run_"+run + "_dif_1_1_3.raw.root";
      if(i == 3) filename = indir1+"/run_"+run + "_dif_1_1_4.raw.root";
      if(i == 4) filename = indir1+"/run_"+run + "_dif_1_1_5.raw.root";
      if(i == 5) filename = indir2+"/run_"+run + "_SLB_2.root";
      if(i == 6) filename = indir2+"/run_"+run + "_SLB_1.root";
      if(i == 7) filename = indir2+"/run_"+run + "_SLB_3.root";
      if(i == 8) filename = indir2+"/run_"+run + "_SLB_0.root";
	    
    

      cout << "merging " << filename << endl;
      f[i] = new TFile(filename,"read");
      if(f[i]->IsOpen()){
	      
	if(i>4) tree[i] = (TTree*)f[i]->Get("slboard");
	else tree[i] = (TTree*)f[i]->Get("fev10");
	      

	tree[i]->SetBranchAddress("chipid", chipID_in[i]);
	//int acq = -999;
	tree[i]->SetBranchAddress("acqNumber", &acqNumber_in[i]);
	//acqNumber_in[i]= acq;

	tree[i]->SetBranchAddress("corrected_bcid" , corrected_bcid_in[i]);
	tree[i]->SetBranchAddress("bcid"           , bcid_in[i]);
	tree[i]->SetBranchAddress("nhits"          , nhits_in[i]);
	tree[i]->SetBranchAddress("badbcid"        , badbcid_in[i]);
	tree[i]->SetBranchAddress("nColumns"       , numCol_in[i]);
	tree[i]->SetBranchAddress("gain_hit_high"  , gain_hit_high_in[i] );
	tree[i]->SetBranchAddress("charge_hiGain"  , charge_high_in[i] );
	tree[i]->SetBranchAddress("gain_hit_low"   , gain_hit_low_in[i] );
	tree[i]->SetBranchAddress("charge_lowGain" , charge_low_in[i] );

	entries[i] = tree[i]->GetEntries();
	if(entries[i]>maxEntries)maxEntries = entries[i];
      }
    }


    fout = new TFile(outfilename,"recreate");
    treeout = new TTree("ecal_raw","ecal_raw");




    TString name;

    treeout->Branch("acqNumber",&acqNumber,"acqNumber/I");

    name= TString::Format("chipid[%i][%i]/I",NSLABS,NCHIP);
    treeout->Branch("chipid",chipID,name);
	
    name= TString::Format("nColumns[%i][%i]/I",NSLABS,NCHIP);
    treeout->Branch("nColumns",numCol,name);
	
    name= TString::Format("bcid[%i][%i][%i]/I",NSLABS,NCHIP,MEMDEPTH);
    treeout->Branch("bcid",bcid,name);
	
    name= TString::Format("corrected_bcid[%i][%i][%i]/I",NSLABS,NCHIP,MEMDEPTH);
    treeout->Branch("corrected_bcid",corrected_bcid,name);
	
    name= TString::Format("badbcid[%i][%i][%i]/I",NSLABS,NCHIP,MEMDEPTH);
    treeout->Branch("badbcid",badbcid,name);
  
    name= TString::Format("nhits[%i][%i][%i]/I",NSLABS,NCHIP,MEMDEPTH);
    treeout->Branch("nhits",nhits,name);

    name= TString::Format("lowGain[%i][%i][%i][%i]/I",NSLABS,NCHIP,MEMDEPTH,NCHANNELS);
    treeout->Branch("charge_lowGain",charge_low,name);

    name= TString::Format("highGain[%i][%i][%i][%i]/I",NSLABS,NCHIP,MEMDEPTH,NCHANNELS);
    treeout->Branch("charge_hiGain",charge_high,name);

    name= TString::Format("tdc[%i][%i][%i][%i]/I",NSLABS,NCHIP,MEMDEPTH,NCHANNELS);
    treeout->Branch("tdc",tdc,name);

    name= TString::Format("gain_hit_low[%i][%i][%i][%i]/I",NSLABS,NCHIP,MEMDEPTH,NCHANNELS);
    treeout->Branch("gain_hit_low",gain_hit_low,name);

    name= TString::Format("gain_hit_high[%i][%i][%i][%i]/I",NSLABS,NCHIP,MEMDEPTH,NCHANNELS);
    treeout->Branch("gain_hit_high",gain_hit_high,name);
  
    //  currentSpill[i]=0;

    int currentSpill=0;
    // for(int i=0;i<maxEntries;i++){
    bool continueScan = true;
    int maxEvents=2000000000;
    int evt = 0;
    bool isSlabEnd[NSLABS];
    for(int i=0;i<NSLABS;i++) isSlabEnd[i]=false;


    while(continueScan && evt<maxEvents){
      treeInit();
      currentSpill=2000000000;//INT_MAX;
      for(int i = 0;i<NSLABS;i++){
	if(f[i]->IsOpen()){
	  tree[i]->GetEntry(entry[i]);
	  //cout<<"-> current = "<<entry[i]<<"    acq = "<<acqNumber_in[i]<<endl;
	  if(currentSpill>acqNumber_in[i] && !isSlabEnd[i])currentSpill=acqNumber_in[i];
	}
      }



      for(int slab = 0;slab<NSLABS;slab++){
	if(currentSpill==acqNumber_in[slab] && entries[slab]>entry[slab]){
	  entry[slab]++;
	  //cout<<"-> current = "<<slab<<"    acq = "<<acqNumber_in[slab]<<endl;


	  for (int k=0; k<NCHIP; k++) {

	    
	    //correct DIF bcid for the offset
	    for (int iraw=0; iraw<MEMDEPTH; iraw++) {
	      if(corrected_bcid_in[slab][k][iraw]<0 || bcid_in[slab][k][iraw]<0) continue;
	      if(slab<5) {
		bcid_in[slab][k][iraw]=bcid_in[slab][k][iraw]-2492;
		corrected_bcid_in[slab][k][iraw]=corrected_bcid_in[slab][k][iraw]-2492;
	      }
	    }

	    int i =0;
	    int loopBCID = 0;
	    for (int iraw=0; iraw<MEMDEPTH; iraw++) {

	      if(bcid_in[slab][k][i]<0) continue;

	      //calculate overrunning bcids
	      // this is supposed to be don in the RAW2ROOT or SLBdecodedROOT macros, in the corrected_bcid
	      if(iraw==0) i = iraw;
	      else if( (bcid_in[slab][k][iraw] - bcid_in[slab][k][iraw-1])>0) i = iraw;
	      else if( (bcid_in[slab][k][iraw] - bcid_in[slab][k][iraw-1])<0){
		i = iraw;
		loopBCID++;
	      } else continue;

	      bcid[slab][k][i] = bcid_in[slab][k][i]+loopBCID*4096;
	      badbcid[slab][k][i]=badbcid_in[slab][k][i];
	      corrected_bcid[slab][k][i]=corrected_bcid_in[slab][k][i];
	      nhits[slab][k][i]=nhits_in[slab][k][i];

	      for (int j=0; j<NCHANNELS; j++) {
		if(mode=="TDC" ) {
		  if(slab <5 ) {
		    charge_high[slab][k][i][j]=charge_low_in[slab][k][i][j];
		    charge_low[slab][k][i][j]=-1;
		    tdc[slab][k][i][j]=charge_high_in[slab][k][i][j];
		    gain_hit_low[slab][k][i][j]=gain_hit_high_in[slab][k][i][j];
		    gain_hit_high[slab][k][i][j]=gain_hit_low_in[slab][k][i][j];
		  } else {
		    charge_high[slab][k][i][j]=charge_high_in[slab][k][i][j];
		    charge_low[slab][k][i][j]=charge_low_in[slab][k][i][j];
		    tdc[slab][k][i][j]=-1;
		    gain_hit_low[slab][k][i][j]=gain_hit_low_in[slab][k][i][j];
		    gain_hit_high[slab][k][i][j]=gain_hit_high_in[slab][k][i][j];
		  }
		}
		if(mode=="HighLow" ) {
		  charge_high[slab][k][i][j]=charge_high_in[slab][k][i][j];
		  charge_low[slab][k][i][j]=charge_low_in[slab][k][i][j];
		  tdc[slab][k][i][j]=-1;
		  gain_hit_low[slab][k][i][j]=gain_hit_low_in[slab][k][i][j];
		  gain_hit_high[slab][k][i][j]=gain_hit_high_in[slab][k][i][j];
		}
			      
	      }
	    }
	    chipID[slab][k]=chipID_in[slab][k];
	    numCol[slab][k]=numCol_in[slab][k];

	  }
	  acqNumber=currentSpill;

	}
	else {
	  //cout<<"rerore"<<endl;
	}



      }
      treeout->Fill();
      //cout<<"SPILL = "<<currentSpill<<endl;


      evt++;
      continueScan = false;
      for(int i = 0;i<NSLABS;i++){
	//cout<<"   slab = "<<i<<"   tot = "<<entries[i]<<"   current = "<<entry[i]  <<endl;
	isSlabEnd[i]=true;
	if(entries[i]>entry[i]){
	  continueScan = true;
	  isSlabEnd[i]=false;
	}
      }

    }


    fout->cd();
    fout->Write(0);
    fout->Close();

  }


  void treeInit() {
    for (int slab=0; slab<NSLABS; slab++) {
      // acqNumber_in[slab]=-9999;
      for (int k=0; k<NCHIP; k++) {
	for (int i=0; i<MEMDEPTH; i++) {
	  bcid[slab][k][i]= -9999;
	  badbcid[slab][k][i]= -999;
	  //  bcid_in[slab][k][i]= -9999;
	  //  badbcid_in[slab][k][i]= 0;
	  corrected_bcid[slab][k][i]= -9999;
	  nhits[slab][k][i]= -9999;
	  // corrected_bcid_in[slab][k][i]= -9999;
	  nhits_in[slab][k][i]= -9999;
	  for (int j=0; j<NCHANNELS; j++) {
	    charge_low[slab][k][i][j]= -999;
	    charge_high[slab][k][i][j]= -999;
	    tdc[slab][k][i][j]= -999;
	    gain_hit_low[slab][k][i][j]= -999;
	    gain_hit_high[slab][k][i][j]= -999;
	    // charge_low_in[slab][k][i][j]= -9999;
	    //charge_high_in[slab][k][i][j]= -9999;
	    //gain_hit_low_in[slab][k][i][j]= -9999;
	    //gain_hit_high_in[slab][k][i][j]= -9999;
	  }
	}

	chipID[slab][k]= -999;
	numCol[slab][k]= -999;

	//  chipID_in[slab][k]= -9999;
	//  numCol_in[slab][k]= -9999;

      }
    }
    acqNumber= -9999;
  }

protected:
  TFile *fout;
  TTree *treeout;

  //int currentSpill[10];


  enum {
    MEMDEPTH=15,
    NCHANNELS=64,
    NCHIP=16,
    NSLABS=9
  };


  int bcid[NSLABS][NCHIP][MEMDEPTH];
  int badbcid[NSLABS][NCHIP][MEMDEPTH];
  int charge_low[NSLABS][NCHIP][MEMDEPTH][NCHANNELS];
  int charge_high[NSLABS][NCHIP][MEMDEPTH][NCHANNELS];
  int tdc[NSLABS][NCHIP][MEMDEPTH][NCHANNELS];
  int gain_hit_low[NSLABS][NCHIP][MEMDEPTH][NCHANNELS];
  int gain_hit_high[NSLABS][NCHIP][MEMDEPTH][NCHANNELS];
  int numCol[NSLABS][NCHIP];
  int chipID[NSLABS][NCHIP];
  int acqNumber;
  int corrected_bcid[NSLABS][NCHIP][MEMDEPTH];
  int nhits[NSLABS][NCHIP][MEMDEPTH];

  int bcid_in[NSLABS][NCHIP][MEMDEPTH];
  int badbcid_in[NSLABS][NCHIP][MEMDEPTH];
  int charge_low_in[NSLABS][NCHIP][MEMDEPTH][NCHANNELS];
  int charge_high_in[NSLABS][NCHIP][MEMDEPTH][NCHANNELS];
  int gain_hit_low_in[NSLABS][NCHIP][MEMDEPTH][NCHANNELS];
  int gain_hit_high_in[NSLABS][NCHIP][MEMDEPTH][NCHANNELS];
  int numCol_in[NSLABS][NCHIP];
  int chipID_in[NSLABS][NCHIP];
  int acqNumber_in[NSLABS];
  int corrected_bcid_in[NSLABS][NCHIP][MEMDEPTH];
  int nhits_in[NSLABS][NCHIP][MEMDEPTH];

};
