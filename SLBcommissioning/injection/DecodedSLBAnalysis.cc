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

std::vector<std::array<float,6>> DecodedSLBAnalysis::HoldscanAnalysis(int sca_=0, int readentries=-1)
{

  float y[15][16][64];
  float ey[15][16][64];
  TH1F* adc[15][16][64];

  int board_id[15];
  int board_slot[15];

    
  for(int i=0; i<15; i++) {
    board_id[i]=-1;
    board_slot[i]=-1;
    for(int j=0; j<16; j++) {
      for(int k=0; k<64; k++) {
	y[i][j][k]=0;
	ey[i][j][k]=0;
	adc[i][j][k]= new TH1F(TString::Format("adc_slab%i_chip%i_chn%i_sca%i",i,j,k,sca_),TString::Format("adc_slab%i_chip%i_chn%i_sca%i",i,j,k,sca_),4096,-0.5,4095.5);
      }
    }
  }



  // --------------
  Long64_t nentries = fChain->GetEntriesFast();
  if(readentries!=-1) nentries=Long64_t(readentries);
  //-----------------
  cout<<"Hold scan analysis "<<endl;
  cout<<"Total number of entries: "<< nentries<<endl;

  // -----------------------------------------------------------------------------------------------------
  // Signal readout
  Long64_t nbytes = 0, nb = 0;
  int n_SLB=0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if ( jentry > 1000 && jentry % 1000 ==0 ) std::cout << "Progress: " << 100.*jentry/nentries <<" %"<<endl;

    if(jentry==0) n_SLB=n_slboards;

    //if(jentry % 1000 !=0 ) continue;


    for(int islboard=0; islboard<n_slboards; islboard++) {
      	if(board_id[islboard]==-1) board_id[islboard]=slboard_id[islboard];
	if(board_slot[islboard]==-1) board_slot[islboard]=slot[islboard];
      //channels starting retriggers
      for(int ichip=0; ichip<16; ichip++) {
	for(int isca=0; isca<15; isca++) {
	  if(isca==sca_ || sca_==15) {
	    for(int ichn=0; ichn<64; ichn++) {
	      if(gain_hit_high[islboard][ichip][isca][ichn]==1 &&
		 ( fabs(bcid[islboard][ichip][isca]-22)<5
		   || fabs(bcid[islboard][ichip][isca]-38)<5
		   || fabs(bcid[islboard][ichip][isca]-54)<5) && charge_hiGain[islboard][ichip][isca][ichn]>200 &&  charge_hiGain[islboard][ichip][isca][ichn]<500)  {
		adc[islboard][ichip][ichn]->Fill(charge_hiGain[islboard][ichip][isca][ichn]);
	      }
	    }
	  }
	}//for isca
      }
    }


  }

  std::vector<std::array<float,6>> result;
  
  for(int islboard=0; islboard<n_slboards; islboard++) {
    for(int ichip=0; ichip<16; ichip++) {
      for(int ichn=0; ichn<64; ichn++) {
	y[islboard][ichip][ichn]=adc[islboard][ichip][ichn]->GetMean();
	ey[islboard][ichip][ichn]=adc[islboard][ichip][ichn]->GetRMS();
	std::array<float,6> temparray;
	temparray[0]=islboard;
	temparray[1]=board_id[islboard];
	temparray[2]=ichip;
	temparray[3]=ichn;
	temparray[4]=y[islboard][ichip][ichn];
	temparray[5]=ey[islboard][ichip][ichn];
	result.push_back(temparray);
      }
    }
  }

  return result;

}


std::vector<std::array<float,9>> DecodedSLBAnalysis::InjectionAnalysis(int readentries=-1)
{

    cout<<"START"<< endl;

  float y[15][16][64][3];
  float ey[15][16][64][3];
  TH1F* high_gain[15][16][64][3];
  TH1F* low_gain[15][16][64][3];

  int board_id[15];
  int board_slot[15];

    
  for(int i=0; i<15; i++) {
    board_id[i]=-1;
    board_slot[i]=-1;
    for(int j=0; j<16; j++) {
      for(int k=0; k<64; k++) {
	for(int l=0; l<3; l++) {
	  y[i][j][k][l]=0;
	  ey[i][j][k][l]=0;
	  high_gain[i][j][k][l]= new TH1F(TString::Format("high_gain_slab%i_chip%i_sca%i_chn%i",i,j,l,k),TString::Format("high_gain_slab%i_chip%i_sca%i_chn%i",i,j,l,k),4096,-0.5,4095.5);
	  low_gain[i][j][k][l]= new TH1F(TString::Format("low_gain_slab%i_chip%i_sca%i_chn%i",i,j,l,k),TString::Format("low_gain_slab%i_chip%i_sca%i_chn%i",i,j,l,k),4096,-0.5,4095.5);
	}
      }
    }
  }


  // --------------
  Long64_t nentries = fChain->GetEntriesFast();
  if(readentries!=-1) nentries=Long64_t(readentries);
  //-----------------
  cout<<"Injection analysis "<<endl;
  cout<<"Total number of entries: "<< nentries<<endl;

  // -----------------------------------------------------------------------------------------------------
  // Signal readout
  Long64_t nbytes = 0, nb = 0;
  int n_SLB=0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if ( jentry > 1000 && jentry % 1000 ==0 ) std::cout << "Progress: " << 100.*jentry/nentries <<" %"<<endl;

    if(jentry==0) n_SLB=n_slboards;

    //if(jentry % 1000 !=0 ) continue;

    for(int islboard=0; islboard<n_slboards; islboard++) {

      if(board_id[islboard]==-1) {
	board_id[islboard]=slboard_id[islboard];
      }
      
      //channels starting retriggers
      for(int ichip=0; ichip<16; ichip++) {
	for(int isca=0; isca<15; isca++) {
	  if(isca>2) continue;
	 for(int ichn=0; ichn<64; ichn++) {
	   if( fabs(bcid[islboard][ichip][isca]-22)<5
		|| fabs(bcid[islboard][ichip][isca]-38)<5
	       || fabs(bcid[islboard][ichip][isca]-54)<5) {

	     if( gain_hit_high[islboard][ichip][isca][ichn]==1 && charge_hiGain[islboard][ichip][isca][ichn]>200)  {
	       high_gain[islboard][ichip][ichn][isca]->Fill(charge_hiGain[islboard][ichip][isca][ichn]);
	     }
	     if( gain_hit_low[islboard][ichip][isca][ichn]==1 && charge_lowGain[islboard][ichip][isca][ichn]>200)  {
	       low_gain[islboard][ichip][ichn][isca]->Fill(charge_lowGain[islboard][ichip][isca][ichn]);
	     }
	   }
	 }
	}
      }
    }
    
  }

  std::vector<std::array<float,9>> result;
  
  for(int islboard=0; islboard<n_slboards; islboard++) {
    for(int ichip=0; ichip<16; ichip++) {
      for(int ichn=0; ichn<64; ichn++) {
	for(int isca=0; isca<3; isca++) {
	  
	  std::array<float,9> temparray;
	  temparray[0]=islboard;
	  temparray[1]=board_id[islboard];

	  temparray[2]=ichip;
	  temparray[3]=ichn;
	  temparray[4]=isca;
	  temparray[5]=high_gain[islboard][ichip][ichn][isca]->GetMean();
	  temparray[6]=high_gain[islboard][ichip][ichn][isca]->GetRMS();
	  temparray[7]=low_gain[islboard][ichip][ichn][isca]->GetMean();
	  temparray[8]=low_gain[islboard][ichip][ichn][isca]->GetRMS();
	  result.push_back(temparray);
	}
      }
    }
  }
  return result;

}

