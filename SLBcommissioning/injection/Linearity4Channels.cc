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

void Linearity4Channels(TString filename_in="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis_TB2020/converter_SLB/convertedfiles/",
		       TString run="calib_02272020_pa1.2fb_trig230_chn0.8.16.24_calib_small_Ascii", int nslboards=6, int channel1=0, int channel2=8, int channel3=16, int channel4=24){
  
  cout<<" Injection folder: "<<filename_in<<endl;

  //slabs, chips, channels, sca, size --> yhigh, eyhigh, ylow, eylow
  std::vector<std::vector<std::array<float,4>>> y;
 
  float x[100];
  float ex[100];

  for(int l=0; l<100; l++) {
    x[l]=0;
    ex[l]=0;
  }
  
  for(int i=0; i<15; i++) {
    for(int j=0; j<16; j++) {
      for(int k=0; k<64; k++) {
	for(int isca=0; isca<3; isca++) {
	std:vector<std::array<float,4>> tempvec;
	  for(int l=0; l<100; l++) {
	    std::array<float,4> temp;
	    temp[0]=0;
	    temp[1]=0;
	    temp[2]=0;
	    temp[3]=0;
	    tempvec.push_back(temp);
	  }
	  y.push_back(tempvec);
	}
      }
    }
  }

  
  
  int count=0;
  for(int i=0; i<15; i++) {
    float pt=0.0;
    if(i==0) pt=0.0;
    if(i==1) pt=0.1;
    if(i==2) pt=0.2;
    if(i==3) pt=0.3;
    if(i==4) pt=0.4;
    if(i==5) pt=0.5;
    if(i==6) pt=0.65;
    if(i==7) pt=0.7;
    if(i==8) pt=0.8;
    if(i==9) pt=0.9;
    if(i==10) pt=1.0;
    if(i==11) pt=1.15;
    if(i==12) pt=1.3;
    if(i==13) pt=1.5;
    if(i==14) pt=1.7;
    
    TString filename=filename_in+TString::Format("%s/injection_small_%.1fV.root",run.Data(),pt);
    if(i==6 || i==11) filename=filename_in+TString::Format("%s/injection_small_%.2fV.root",run.Data(),pt);

    cout<<" i: "<<i<<endl;
    cout<<" Injection file: "<<filename_in<<endl;

    DecodedSLBAnalysis ss(filename);
    std::vector<std::array<float,8>> injecvalues=ss.InjectionAnalysis(-1);
    
    x[count]=(3.3-pt);
    
    for(int j=0; j<injecvalues.size(); j++) {
      y.at(j).at(i)[0]=injecvalues.at(j)[4]; //yhigh
      y.at(j).at(i)[1]=injecvalues.at(j)[5]; //eyhigh
      
      y.at(j).at(i)[2]=injecvalues.at(j)[6]; //ylow
      y.at(j).at(i)[3]=injecvalues.at(j)[7]; //eylow
    }
    count++;
    
  }


  
  TGraphErrors* injection_high_sca0[15][16][64];
  TGraphErrors* injection_high_sca1[15][16][64];
  TGraphErrors* injection_high_sca2[15][16][64];

  TGraphErrors* injection_low_sca0[6][16][64];
  TGraphErrors* injection_low_sca1[6][16][64];
  TGraphErrors* injection_low_sca2[6][16][64];


  int counter=0;
  for(int i=0; i<nslboards; i++) {
    for(int j=0; j<16; j++) {
      for(int k=0; k<64; k++) {
	for(int isca=0; isca<3; isca++) {

	  float y0[100], ey0[100];
	  float y1[100], ey1[100];

	  for(int ivec=0; ivec<count; ivec++) {
	    y0[ivec]=y.at(counter).at(ivec)[0];
	    ey0[ivec]=y.at(counter).at(ivec)[1];
	    y1[ivec]=y.at(counter).at(ivec)[2];
	    ey1[ivec]=y.at(counter).at(ivec)[3];
	  }

	  counter++;

	  if(isca==0) {
	    injection_high_sca0[i][j][k]= new TGraphErrors(count,x,y0,ex,ey0);
	    injection_low_sca0[i][j][k]= new TGraphErrors(count,x,y1,ex,ey1);
	  }
	  if(isca==1) {
	    injection_high_sca1[i][j][k]= new TGraphErrors(count,x,y0,ex,ey0);
	    injection_low_sca1[i][j][k]= new TGraphErrors(count,x,y1,ex,ey1);
	  }
	  if(isca==2) {
	    injection_high_sca2[i][j][k]= new TGraphErrors(count,x,y0,ex,ey0);
	    injection_low_sca2[i][j][k]= new TGraphErrors(count,x,y1,ex,ey1);
	  }	  

	}
      }
    }
  }

  TFile *file_summary = new TFile("results/graphs_"+run+".root" , "RECREATE");
  file_summary->cd();
  
  for(int i=0; i<nslboards; i++) {
    for(int j=0; j<16; j++) {
      for(int k=0; k<64; k++) {
	injection_high_sca0[i][j][k]->SetName(TString::Format("high_gain_layer%i_chip%i_chn%i_sca0",i,j,k));
	injection_high_sca1[i][j][k]->SetName(TString::Format("high_gain_layer%i_chip%i_chn%i_sca1",i,j,k));
	injection_high_sca2[i][j][k]->SetName(TString::Format("high_gain_layer%i_chip%i_chn%i_sca2",i,j,k));
	injection_low_sca0[i][j][k]->SetName(TString::Format("low_gain_layer%i_chip%i_chn%i_sca0",i,j,k));
	injection_low_sca1[i][j][k]->SetName(TString::Format("low_gain_layer%i_chip%i_chn%i_sca1",i,j,k));
	injection_low_sca2[i][j][k]->SetName(TString::Format("low_gain_layer%i_chip%i_chn%i_sca2",i,j,k));

	injection_high_sca0[i][j][k]->Write();
	injection_high_sca1[i][j][k]->Write();
	injection_high_sca2[i][j][k]->Write();
	injection_low_sca0[i][j][k]->Write();
	injection_low_sca1[i][j][k]->Write();
	injection_low_sca2[i][j][k]->Write();
      }
    }
  }
  
  /*   for(int i=0; i<nslboards; i++) {

      TCanvas * canvas = new TCanvas(TString::Format("canvas_%i",i),TString::Format("canvas_%i",i),1200,1000);
      canvas->Divide(4,3);
      for(int ichan=0; ichan<4; ichan++) {
	int chan=0;
	if(ichan==0) chan=channel1;
	if(ichan==1) chan=channel2;
	if(ichan==2) chan=channel3;
	if(ichan==3) chan=channel4;

	TLegend *leg = new TLegend(0.1,0.5,0.25,0.85);

	//sca0
	for(int j=0; j<16; j++) {
	  canvas->cd(1+ichan);
	  if(j<9) {
	    injection_high_sca0[i][j][chan]->SetLineColor(j+1);
	    injection_high_sca0[i][j][chan]->SetLineStyle(j+1);
	    injection_high_sca0[i][j][chan]->SetMarkerColor(j+1);
	    injection_high_sca0[i][j][chan]->SetMarkerStyle(j+1);
	  } else {
	    injection_high_sca0[i][j][chan]->SetLineColor(j-9+2);
	    injection_high_sca0[i][j][chan]->SetLineStyle(j-9+1);
	    injection_high_sca0[i][j][chan]->SetMarkerColor(j-9+2);
	    injection_high_sca0[i][j][chan]->SetMarkerStyle(j-9+1);
	  }
	  leg->AddEntry(injection_high_sca0[i][j][chan],TString::Format("ASIC %i, SCA %i",j,0),"lep");
	  if(j==0) {
	    injection_high_sca0[i][j][chan]->SetTitle(TString::Format("Layer %i, channel: %i",i,chan));
	    injection_high_sca0[i][j][chan]->GetYaxis()->SetRangeUser(0,500);
	    injection_high_sca0[i][j][chan]->Draw("alp");
	  }
	  else injection_high_sca0[i][j][chan]->Draw("lp");
	}
	leg->Draw();

	//sca1
	for(int j=0; j<16; j++) {
	  canvas->cd(5+ichan);
	  if(j<9) {
	    injection_high_sca1[i][j][chan]->SetLineColor(j+1);
	    injection_high_sca1[i][j][chan]->SetLineStyle(j+1);
	    injection_high_sca1[i][j][chan]->SetMarkerColor(j+1);
	    injection_high_sca1[i][j][chan]->SetMarkerStyle(j+1);
	  } else {
	    injection_high_sca1[i][j][chan]->SetLineColor(j-9+2);
	    injection_high_sca1[i][j][chan]->SetLineStyle(j-9+1);
	    injection_high_sca1[i][j][chan]->SetMarkerColor(j-9+2);
	    injection_high_sca1[i][j][chan]->SetMarkerStyle(j-9+1);
	  }
	  if(j==0) {
	    injection_high_sca1[i][j][chan]->SetTitle(TString::Format("Layer %i, channel: %i",i,chan));
	    injection_high_sca1[i][j][chan]->GetYaxis()->SetRangeUser(0,500);
	    injection_high_sca1[i][j][chan]->Draw("alp");
	  }
	  else injection_high_sca1[i][j][chan]->Draw("lp");
	}
	leg->Draw();

	//sca2
	for(int j=0; j<16; j++) {
	  canvas->cd(9+ichan);
	  if(j<9) {
	    injection_high_sca2[i][j][chan]->SetLineColor(j+1);
	    injection_high_sca2[i][j][chan]->SetLineStyle(j+1);
	    injection_high_sca2[i][j][chan]->SetMarkerColor(j+1);
	    injection_high_sca2[i][j][chan]->SetMarkerStyle(j+1);
	  } else {
	    injection_high_sca2[i][j][chan]->SetLineColor(j-9+2);
	    injection_high_sca2[i][j][chan]->SetLineStyle(j-9+1);
	    injection_high_sca2[i][j][chan]->SetMarkerColor(j-9+2);
	    injection_high_sca2[i][j][chan]->SetMarkerStyle(j-9+1);
	  }
	  if(j==0) {
	    injection_high_sca2[i][j][chan]->SetTitle(TString::Format("Layer %i, channel: %i",i,chan));
	    injection_high_sca2[i][j][chan]->GetYaxis()->SetRangeUser(0,500);
	    injection_high_sca2[i][j][chan]->Draw("alp");
	  }
	  else injection_high_sca2[i][j][chan]->Draw("lp");
	}
	leg->Draw();
      }
    canvas->Print(TString::Format("results/%s_layer%i.root",run.Data(),i));
  }

  */
				 
  
}

