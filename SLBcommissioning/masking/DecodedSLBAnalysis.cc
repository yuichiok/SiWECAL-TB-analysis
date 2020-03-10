//# Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)
#define DecodedSLBAnalysis_cxx
#include "DecodedSLBAnalysis.h"
#include <TPaveStats.h>
#include <TH2.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <TFitResult.h>
#include <TSpectrum.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

std::vector<std::array<int,8>>  DecodedSLBAnalysis::NoiseLevels(int acqwindow=150)
{

  
  int retrigger_start[15][16][64];
  int retrigger_train[15][16][64];
  int trigger[15][16][64];
  int adc4[15][16][64];

  int n_SLB=0;


  for(int ilayer=0; ilayer<15; ilayer++) {
    for(int ichip=0; ichip<16; ichip++) {     
      for(int ichn=0; ichn<64; ichn++) {
	retrigger_start[ilayer][ichip][ichn]=0;
	retrigger_train[ilayer][ichip][ichn]=0;
	trigger[ilayer][ichip][ichn]=0;
	adc4[ilayer][ichip][ichn]=0;
     }
    }
  }
  
  Long64_t nentries = fChain->GetEntriesFast();

  float expected= 0.25*nentries*float(acqwindow)/60.;

  // -----------------------------------------------------------------------------------------------------   
  // SCA analysis
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(jentry==0) n_SLB=n_slboards;

    if ( jentry > 1000 && jentry % 1000 ==0 ) std::cout << "Progress: " << 100.*jentry/nentries <<" %"<<endl;

    for(int ilayer=0; ilayer<n_slboards; ilayer++) {
      
      for(int ichip=0; ichip<16; ichip++) {

	bool first_retrig=false;
      
	for(int isca=0; isca<15; isca++) {
	  bool retrig=false;
	  bool burst=false;
	  if(bcid[ilayer][ichip][isca]<0) continue;
	
	  if(badbcid[ilayer][ichip][isca]>2 ) {
	    retrig=true;
	    if(first_retrig==false) {
	      first_retrig=true;
	      for(int ichn=0;ichn<64;ichn++) {
		if(gain_hit_high[ilayer][ichip][isca][ichn]==1) retrigger_start[ilayer][ichip][ichn]++;
		if(gain_hit_high[ilayer][ichip][isca][ichn]==1 && (charge_hiGain[ilayer][ichip][isca][ichn]<100 || charge_hiGain[ilayer][ichip][isca][ichn]>1000) ) adc4[ilayer][ichip][ichn]++;
	      }
	    }
	    for(int ichn=0;ichn<64;ichn++) {
	      if(gain_hit_high[ilayer][ichip][isca][ichn]==1) retrigger_train[ilayer][ichip][ichn]++;
	    }
	  }

	  if(badbcid[ilayer][ichip][isca]==0) {
	    for(int ichn=0;ichn<64;ichn++) {
	      if(gain_hit_high[ilayer][ichip][isca][ichn]==1 ) trigger[ilayer][ichip][ichn]++;
	    }
	  }
	    

	}//isca
      }//ichip
    }//ilayer
	  
      
  }

  std::vector<std::array<int,8>>  result;
  
  for(int ilayer=0; ilayer<n_SLB; ilayer++) {
    for(int ichip=0; ichip<16; ichip++) {     
      for(int ichn=0; ichn<64; ichn++) {
	std::array<int,8> temp;
	temp[0]=ilayer;
	temp[1]=ichip;
	temp[2]=ichn;
	temp[3]=int(expected);
	temp[4]=trigger[ilayer][ichip][ichn];
	temp[5]=retrigger_start[ilayer][ichip][ichn];
	temp[6]=retrigger_train[ilayer][ichip][ichn];
	temp[7]=adc4[ilayer][ichip][ichn];
	result.push_back(temp);
      }
    }
  }
	

  return result;
}




