#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "../../style/Style.C"
#include "../../style/Labels.C"

using namespace std;


int Noise_LS(){


  SetIrlesStyle();
  gStyle->SetOptFit(0111); 
  gStyle->SetOptStat(1110);

  gStyle->SetOptTitle(1);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.2);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleStyle(0);
  gStyle->SetMarkerSize(0.8);

  TString angle = "angle45";
  // int asu = 8; // asu number from 1 (closest to the dif) to 8.
  int channel[16]= {4,11,15,18,20,23,27,30,33,36,37,40,42,45,48,49};
  
  for(int iasu=8; iasu<9; iasu++) {
    //open file
    TString filename =TString::Format("../results_pedestal/Pedestal_dif_1_1_1_ASU%i_%s.root",iasu,angle.Data());
    TFile *f = new TFile(filename);
  
    TGraphErrors *pedestal_channel[16];
    TGraphErrors *width_channel[16];

    int ichip=10+16*(iasu-1);

    for (int ichn=0; ichn<16; ichn++) {
      double x[15], y[15], ey[15];
      int nsca=0;
      
      for(int isca=0; isca<15; isca ++ ) {
	x[isca]=isca;
	y[isca]=0;
	ey[isca]=0;
		
	TH1F *hpedestal_good  = (TH1F*)f->Get(TString::Format("ped_chip%i_chn%i_sca%i",ichip,channel[ichn],isca));
	
	if(hpedestal_good->GetEntries()>50) {
	  y[isca]=hpedestal_good->GetMean();
	  ey[isca]=hpedestal_good->GetRMS();
	  nsca++;
	}
      }  // end sca

      pedestal_channel[ichn]=new TGraphErrors(nsca,x,y,0,ey);
      width_channel[ichn]=new TGraphErrors(nsca,x,ey,0,0);

    }//ichn

    TLegend *leg1 = new TLegend(0.2,0.2,0.5,0.5);
    TLegend *leg2 = new TLegend(0.55,0.2,0.85,0.5);

    TCanvas *canvas = new TCanvas(TString::Format("canvas_asu%i_%s",iasu,angle.Data()),TString::Format("canvas_asu%i_%s",iasu,angle.Data()),1400,800);
    canvas->Divide(2,1);
    canvas->cd(1);
    for(int ichn=0;ichn<16;ichn++) {
     if(ichn<8) {
	pedestal_channel[ichn]->SetMarkerColor(1+ichn);
	pedestal_channel[ichn]->SetMarkerStyle(20+ichn);
	pedestal_channel[ichn]->SetLineColor(1+ichn);
	pedestal_channel[ichn]->SetLineStyle(1);
	pedestal_channel[ichn]->SetLineWidth(2);
      } else {
	pedestal_channel[ichn]->SetMarkerColor(1+ichn-8);
	pedestal_channel[ichn]->SetMarkerStyle(20+ichn-8);
	pedestal_channel[ichn]->SetLineColor(1+ichn-8);
	pedestal_channel[ichn]->SetLineStyle(2);
	pedestal_channel[ichn]->SetLineWidth(2);
      }
       if(ichn==0) {
	pedestal_channel[ichn]->SetTitle(TString::Format("ASU%i %s",iasu,angle.Data()));
	pedestal_channel[ichn]->GetXaxis()->SetTitle("# SCA");
	pedestal_channel[ichn]->GetYaxis()->SetTitle("pedestal [ADC]");
	pedestal_channel[ichn]->GetYaxis()->SetRangeUser(200,400);
	pedestal_channel[ichn]->Draw("alp");
      } else pedestal_channel[ichn]->Draw("lp");
      
    }//ichn

    canvas->cd(2);
    for(int ichn=0;ichn<16;ichn++) {
      if(ichn<8) {
	width_channel[ichn]->SetMarkerColor(1+ichn);
	width_channel[ichn]->SetMarkerStyle(20+ichn);
	width_channel[ichn]->SetLineColor(1+ichn);
	width_channel[ichn]->SetLineStyle(1);
	width_channel[ichn]->SetLineWidth(2);
      } else {
	width_channel[ichn]->SetMarkerColor(1+ichn-8);
	width_channel[ichn]->SetMarkerStyle(20+ichn-8);
	width_channel[ichn]->SetLineColor(1+ichn-8);
	width_channel[ichn]->SetLineStyle(2);
	width_channel[ichn]->SetLineWidth(2);
      }
      if(ichn==0) {
	width_channel[ichn]->SetTitle(TString::Format("ASU%i %s",iasu,angle.Data()));
	width_channel[ichn]->GetXaxis()->SetTitle("# SCA");
	width_channel[ichn]->GetYaxis()->SetTitle("noise [ADC]");
	width_channel[ichn]->GetYaxis()->SetRangeUser(-1,7);
	width_channel[ichn]->Draw("alp");
      } else width_channel[ichn]->Draw("lp");
      if(ichn<8) leg1->AddEntry(width_channel[ichn],TString::Format("channel=%i",channel[ichn]),"lp");
      else leg2->AddEntry(width_channel[ichn],TString::Format("channel=%i",channel[ichn]),"lp");
      
    }//ichn
    leg1->Draw();
    leg2->Draw();

    canvas->Print(TString::Format("canvas_asu%i_%s.png",iasu,angle.Data()));
 
  }

   return 0;
 
}
