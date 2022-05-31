#include "TFile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TH1.h"
#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TString.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "TLatex.h"
#include "../include/utils.h"

//#include "../../style/Style.C"
//#include "../../style/Labels.C"

using namespace std;

void readwritepedestals(){

  float high[15][16][64][15]={0.};
  float low[15][16][64][15]={0.};

  float low_av[15][16]={0.};
  float high_av[15][16]={0.};
  float n_av[15][16]={0.};

  ReadPedestalsProtoCovariance("Pedestal_method2_Pedestal_run_050306_injection_merged_highgain.txt");

  for(int layer=0;layer<15;layer++){
    for(int i=0;i<16;i++){
      n_av[layer][i]=0;
      for(int j=0; j<64; j++) {
	for(int k=0; k<15; k++) {
	  high[layer][i][j][k]=ped_mean_slboard.at(layer).at(i).at(j).at(k);
	  if(high[layer][i][j][k]>0) {
	    n_av[layer][i]++;
	    high_av[layer][i]+=high[layer][i][j][k];
	  }
	}
      }
      high_av[layer][i]/=n_av[layer][i];
    }
  }

  ReadPedestalsProtoCovariance("Pedestal_method2_Pedestal_run_050306_injection_merged_lowgain.txt");

  for(int layer=0;layer<15;layer++){
    for(int i=0;i<16;i++){
      n_av[layer][i]=0;
      for(int j=0; j<64; j++) {
	for(int k=0; k<15; k++) {
	  low[layer][i][j][k]=ped_mean_slboard.at(layer).at(i).at(j).at(k);
	  if(low[layer][i][j][k]>0) {
	    n_av[layer][i]++;
	    low_av[layer][i]+=low[layer][i][j][k];
	  }
	}
      }
      low_av[layer][i]/=n_av[layer][i];
    }
  }
  


  for(int layer=0; layer<15; layer++) {
    ofstream fout_ped(TString::Format("Ped_Calib_Core0_SlabAdd%i_Asu0.txt",layer).Data(),ios::out);
    fout_ped<<"Skiroc Ch SCA Gain PedestalValue [float in ADC count]"<<endl;
    for(int i=0;i<16;i++){
      for(int j=0; j<64; j++) {
	fout_ped<<TString::Format("Skiroc %i Ch %i ",i,j)<<endl; 
	fout_ped<<"Gain 0 "<<endl;
	for(int k=0; k<15; k++) {
	  if(low[layer][i][j][k]>0) fout_ped<<setprecision(6)<<low[layer][i][j][k]<<" ";
	  else fout_ped<<setprecision(6)<<low_av[layer][i]<<" ";
	}
	fout_ped<<endl;
	fout_ped<<"Gain 1 "<<endl;
	for(int k=0; k<15; k++) {
	  if(high[layer][i][j][k]>0) fout_ped<<setprecision(6)<<high[layer][i][j][k]<<" ";
	  else fout_ped<<setprecision(6)<<high_av[layer][i]<<" ";
	}
	fout_ped<<endl;
      }
    }
  }

}

