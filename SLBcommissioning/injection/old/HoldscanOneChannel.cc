//# Copyright 2020  Adrián Irles

#include "TROOT.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "DecodedSLBAnalysis.cc"

void HoldscanOneChannel(TString filename_in="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis_TB2020/converter_SLB/convertedfiles/", TString run="Run_ILC_trig230_chn0_DAC1.15V_Ascii", int nslboards=8, int channel=0){
  
  cout<<" Holdscan file: "<<filename_in<<endl;

  float y[15][16][64][30];
  float x[30];
  float ex[30];
  float ey[15][16][64][30];

  for(int l=0; l<30; l++) {
    x[l]=0;
    ex[l]=0;
    for(int i=0; i<15; i++) {
      for(int j=0; j<16; j++) {
  	for(int k=0; k<64; k++) {
  	  y[i][j][k][l]=0;
  	  ey[i][j][k][l]=0;
  	}
      }
    }
  }

  int count=0;
  for(int i=50; i<210; i=i+10) {
    cout<<" i: "<<i<<endl;

    TString filename=filename_in+TString::Format("%s/holdscan_hold%i.root",run.Data(),i);
    cout<<" Holdscan file: "<<filename_in<<endl;
    DecodedSLBAnalysis ss(filename);
    std::vector<std::array<float,5>> holdvalues =ss.HoldscanAnalysis(-1);
    
    x[count]=i;
    
    for(int i=0; i<holdvalues.size(); i++) {
      y[int(holdvalues.at(i)[0])][int(holdvalues.at(i)[1])][int(holdvalues.at(i)[2])][count]=holdvalues.at(i)[3];
      ey[int(holdvalues.at(i)[0])][int(holdvalues.at(i)[1])][int(holdvalues.at(i)[2])][count]=holdvalues.at(i)[4];
    }
    count++;
    
  }

  TGraphErrors* holdscan[15][16][64];

  for(int i=0; i<nslboards; i++) {
    for(int j=0; j<16; j++) {
      for(int k=0; k<64; k++) {

  	holdscan[i][j][k]= new TGraphErrors(count,x,y[i][j][k],ex,ey[i][j][k]);
      }
    }
  }

  
  for(int i=0; i<nslboards; i++) {
    TCanvas * canvas = new TCanvas(TString::Format("canvas_%i",i),TString::Format("canvas_%i",i),1200,1200);
    TLegend *leg = new TLegend(0.1,0.5,0.25,0.85);
    for(int j=0; j<16; j++) {
      if(j<9) {
	holdscan[i][j][channel]->SetLineColor(j+1);
	holdscan[i][j][channel]->SetLineStyle(j+1);
	holdscan[i][j][channel]->SetMarkerColor(j+1);
	holdscan[i][j][channel]->SetMarkerStyle(j+1);
      } else {
	holdscan[i][j][channel]->SetLineColor(j-9+2);
	holdscan[i][j][channel]->SetLineStyle(j-9+1);
	holdscan[i][j][channel]->SetMarkerColor(j-9+2);
	holdscan[i][j][channel]->SetMarkerStyle(j-9+1);
      }
      leg->AddEntry(holdscan[i][j][channel],TString::Format("ASIC %i",j),"lep");
      if(j==0) {
	holdscan[i][j][channel]->SetTitle(TString::Format("Layer %i",i));
	holdscan[i][j][channel]->GetYaxis()->SetRangeUser(0,500);
	holdscan[i][j][channel]->Draw("alp");
      }
      else holdscan[i][j][channel]->Draw("lp");
    }
    leg->Draw();
    canvas->Print(TString::Format("results/%s_layer%i.root",run.Data(),i));
  }



  
}

