//# Copyright 2020  Adrián Irles

#include "TROOT.h"
#include "TFile.h"
#include "DecodedSLBAnalysis.cc"
#include "TGraphErrors.h"
#include "TLegend.h"

//void HoldscanPlots8Channels(TString run="07282021_dac1.2V_small", int nslboards=15, int channel1=0, int channel2=8, int channel3=16, int channel4=24, int channel5=32, int channel6=40, int channel7=48, int channel8=56){
void HoldscanPlots8Channels(TString run="08312021_dac1.2V_small", int nslboards=15, int channel1=0, int channel2=8, int channel3=16, int channel4=24, int channel5=32, int channel6=40, int channel7=48, int channel8=56){
  

  TFile *file_read = new TFile("results/graphs_holdscan_"+run+"_sca15.root" , "READ");
  file_read->cd();

  TGraphErrors* holdscan[15][16][64];

  TH1F *h1 = (TH1F*)file_read->Get("layer_slboard_relation");

  
  for(int i=0; i<nslboards; i++) {
    for(int j=0; j<16; j++) {
      for(int k=0; k<64; k++) {
	holdscan[i][j][k]=(TGraphErrors*)file_read->Get(TString::Format("holdscan_layer%i_slboard%i_chip%i_chn%i",i,int(h1->GetBinContent(i+1)),j,k));
      }
    }
  }

  
  for(int i=0; i<nslboards; i++) {
    
    TCanvas * canvas = new TCanvas(TString::Format("canvas_%i",i),TString::Format("canvas_%i",i),1200,600);
    canvas->Divide(2,4);

    for(int ichan=0; ichan<8; ichan++) {
      int chan=0;
      if(ichan==0) chan=channel1;
      if(ichan==1) chan=channel2;
      if(ichan==2) chan=channel3;
      if(ichan==3) chan=channel4;
      if(ichan==4) chan=channel5;
      if(ichan==5) chan=channel6;
      if(ichan==6) chan=channel7;
      if(ichan==7) chan=channel8;

      TLegend *leg = new TLegend(0.1,0.5,0.25,0.85);
     
      for(int j=0; j<16; j++) {
	
	canvas->cd(1+ichan);
	
	if(j<9) {
	  holdscan[i][j][chan]->SetLineColor(j+1);
	  holdscan[i][j][chan]->SetLineStyle(j+1);
	  holdscan[i][j][chan]->SetMarkerColor(j+1);
	  holdscan[i][j][chan]->SetMarkerStyle(j+1);
	} else {
	  holdscan[i][j][chan]->SetLineColor(j-9+2);
	  holdscan[i][j][chan]->SetLineStyle(j-9+1);
	  holdscan[i][j][chan]->SetMarkerColor(j-9+2);
	  holdscan[i][j][chan]->SetMarkerStyle(j-9+1);
	}
	leg->AddEntry(holdscan[i][j][chan],TString::Format("ASIC %i",j),"lep");
	if(j==0) {
	  holdscan[i][j][chan]->SetTitle(TString::Format("Layer %i, SLBoard %i, channel: %i",i,int(h1->GetBinContent(i+1)),chan));
	  holdscan[i][j][chan]->GetYaxis()->SetRangeUser(0,500);
	  holdscan[i][j][chan]->Draw("alp");
	}
	else holdscan[i][j][chan]->Draw("lp");
      }
      leg->Draw();
    }
    canvas->Print(TString::Format("results/%s_layer%i_slboard%i.root",run.Data(),i,int(h1->GetBinContent(i+1))));
  }
  
				 
  
}

