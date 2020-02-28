//# Copyright 2020  Adrián Irles

#include "TROOT.h"
#include "TFile.h"
#include "DecodedSLBAnalysis.cc"
#include "TGraphErrors.h"
#include "TLegend.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"

void PlotLinearity4Channels(TString run="calib_02272020_pa1.2fb_trig230_chn0.8.16.24_calib_small_Ascii", int nslboards=6, int channel1=0, int channel2=8, int channel3=16, int channel4=24){
  
  cout<<" Graphs file: "<<run<<endl;
 
  TGraphErrors* injection_high[15][16][64][3];
  //  TGraphErrors* injection_high_sca1[15][16][64][3];
  TGraphErrors* injection_high_sca2[15][16][64][3];

  TGraphErrors* injection_low[15][16][64][3];
  //  TGraphErrors* injection_low_sca1[15][16][64];
  //TGraphErrors* injection_low_sca2[15][16][64];


  TFile *file_read = new TFile("results/graphs_"+run+".root" , "READ");
  file_read->cd();
  
  for(int i=0; i<nslboards; i++) {
    for(int j=0; j<16; j++) {
      for(int k=0; k<64; k++) {
	
	injection_high[i][j][k][0]=(TGraphErrors*)file_read->Get(TString::Format("high_gain_layer%i_chip%i_chn%i_sca0",i,j,k));
	injection_high[i][j][k][1]=(TGraphErrors*)file_read->Get(TString::Format("high_gain_layer%i_chip%i_chn%i_sca1",i,j,k));
	injection_high[i][j][k][2]=(TGraphErrors*)file_read->Get(TString::Format("high_gain_layer%i_chip%i_chn%i_sca2",i,j,k));

	injection_low[i][j][k][0]=(TGraphErrors*)file_read->Get(TString::Format("low_gain_layer%i_chip%i_chn%i_sca0",i,j,k));
	injection_low[i][j][k][1]=(TGraphErrors*)file_read->Get(TString::Format("low_gain_layer%i_chip%i_chn%i_sca1",i,j,k));
	injection_low[i][j][k][2]=(TGraphErrors*)file_read->Get(TString::Format("low_gain_layer%i_chip%i_chn%i_sca2",i,j,k));
      }
    }
  }

  file_read->Close();
  
  TFile *file_write = new TFile("results/canvas_"+run+".root" , "RECREATE");
  file_write->cd();
  for(int i=0; i<nslboards; i++) {
    for(int j=0; j<16; j++) {
      for(int isca=0; isca<3; isca++) {

	TCanvas * canvas = new TCanvas(TString::Format("canvas_layer%i_chip%i_sca%i",i,j,isca),TString::Format("canvas_layer%i_chip%i_sca%i",i,j,isca),1200,1000);
	canvas->Divide(2,2);
	for(int ichan=0; ichan<4; ichan++) {
	  int chan=0;
	  if(ichan==0) chan=channel1;
	  if(ichan==1) chan=channel2;
	  if(ichan==2) chan=channel3;
	  if(ichan==3) chan=channel4;
	  
	  TLegend *leg = new TLegend(0.1,0.5,0.25,0.85);
	  canvas->cd(1+ichan);
	  injection_high[i][j][chan][isca]->SetLineColor(1);
	  injection_high[i][j][chan][isca]->SetLineStyle(1);
	  injection_high[i][j][chan][isca]->SetMarkerColor(1);
	  injection_high[i][j][chan][isca]->SetMarkerStyle(21);
	  injection_high[i][j][chan][isca]->SetLineWidth(2);
	  injection_high[i][j][chan][isca]->GetYaxis()->SetRangeUser(100,500);

	  injection_high[i][j][chan][isca]->SetTitle(TString::Format("Layer %i, Chip %i, channel: %i, SCA %i",i,j,chan,isca));
	  injection_high[i][j][chan][isca]->Draw("alp");
	  leg->AddEntry(injection_high[i][j][chan][isca],"High Gain","lep");

	  injection_low[i][j][chan][isca]->SetLineColor(2);
	  injection_low[i][j][chan][isca]->SetLineStyle(2);
	  injection_low[i][j][chan][isca]->SetMarkerColor(2);
	  injection_low[i][j][chan][isca]->SetMarkerStyle(22);
	  injection_low[i][j][chan][isca]->SetLineWidth(2);
	  injection_low[i][j][chan][isca]->Draw("lp");
	  leg->AddEntry(injection_low[i][j][chan][isca],"Low Gain","lep");
	  leg->Draw();
	}
	canvas->Write();
	delete	canvas;
      }
    }
  }



				 
  
}

