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

int Noise_LS_2(){


  SetIrlesStyle();
  gStyle->SetOptFit(0111); 
  gStyle->SetOptStat(1110);

  gStyle->SetOptTitle(1);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.2);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleStyle(0);
  gStyle->SetMarkerSize(1.2);

  TString angle = "angle0";
  // int asu = 8; // asu number from 1 (closest to the dif) to 8.
  int channel[16]= {4,11,15,18,20,23,27,30,33,36,37,40,42,45,48,49};

  TGraphErrors *pedestal_channel[16];
  TGraphErrors *width_channel[16];
  
  for (int ichn=0; ichn<16; ichn++) {
    double x[122];
    double y[122], ey[122];

    for(int i=0; i<122; i ++ ) {
	y[i]=0;
	ey[i]=0;
	x[i]=0.;
    }
    int nsca=0;
    
    for(int iasu=1; iasu<9; iasu++) {
      //open file
      TString filename =TString::Format("../results_pedestal/Pedestal_dif_1_1_1_ASU%i_%s.root",iasu,angle.Data());
      TFile *f = new TFile(filename);
           
      int ichip=10+16*(iasu-1);
                  
      for(int isca=0; isca<15; isca ++ ) {
	
	x[nsca]=isca+20*(iasu-1);
	//y[nsca]=0;
	//ey[nsca]=0;
	nsca++;
	//	cout<<nsca<<" "<<isca+20*(iasu-1)<<endl;
	TH1F *hpedestal_good  = (TH1F*)f->Get(TString::Format("ped_chip%i_chn%i_sca%i",ichip,channel[ichn],isca));

	///cout<<TString::Format("ped_chip%i_chn%i_sca%i",ichip,channel[ichn],isca)<<endl;
	if(hpedestal_good) {
	  y[nsca]=hpedestal_good->GetMean();
	  ey[nsca]=hpedestal_good->GetRMS();
	  //cout<<filename<<" "<<TString::Format("ped_chip%i_chn%i_sca%i",ichip,channel[ichn],isca)<<" "<<hpedestal_good->GetEntries()<<" "<<nsca<<" "<<x[nsca]<<" "<<y[nsca]<<endl;

	}
      } // end sca
      // f->Close();

    }//asu

   pedestal_channel[ichn]=new TGraphErrors(nsca,x,y,0,ey);
   width_channel[ichn]=new TGraphErrors(nsca,x,ey,0,0);
  }//chn


    TLegend *leg1 = new TLegend(0.2,0.2,0.5,0.5);
    TLegend *leg2 = new TLegend(0.55,0.2,0.85,0.5);

    TCanvas *canvas = new TCanvas(TString::Format("pedestal_allasu_%s",angle.Data()),TString::Format("canvas_allasu_%s",angle.Data()),1400,800);
    canvas->Divide(2,1);
    canvas->cd(1);
    for(int ichn=0;ichn<16;ichn++) {
      if(ichn> 4) pedestal_channel[ichn]->SetMarkerStyle(16+ichn);
      if(ichn<8) {
	pedestal_channel[ichn]->SetMarkerColor(1+ichn);
	pedestal_channel[ichn]->SetLineColor(1+ichn);
	pedestal_channel[ichn]->SetLineStyle(1);
	pedestal_channel[ichn]->SetLineWidth(2);
      } else {  
	pedestal_channel[ichn]->SetMarkerColor(1+ichn-8);
	pedestal_channel[ichn]->SetLineColor(1+ichn-8);
	pedestal_channel[ichn]->SetLineStyle(2);
	pedestal_channel[ichn]->SetLineWidth(2);
      }
      if(ichn<3) pedestal_channel[ichn]->SetMarkerStyle(ichn+2);
      else  pedestal_channel[ichn]->SetMarkerStyle(16+ichn);

      
       if(ichn==0) {
	pedestal_channel[ichn]->SetTitle(TString::Format("all ASU %s",angle.Data()));
	pedestal_channel[ichn]->GetXaxis()->SetTitle("# SCA + 20x(ASU-1)");
	pedestal_channel[ichn]->GetYaxis()->SetTitle("pedestal [ADC]");
	pedestal_channel[ichn]->GetYaxis()->SetRangeUser(200,400);
	pedestal_channel[ichn]->Draw("ap");
      } else pedestal_channel[ichn]->Draw("p");
      
    }//ichn

    canvas->cd(2);
    for(int ichn=0;ichn<16;ichn++) {
      if(ichn> 4) width_channel[ichn]->SetMarkerStyle(16+ichn);
      if(ichn<8) {
	width_channel[ichn]->SetMarkerColor(1+ichn);
	width_channel[ichn]->SetLineColor(1+ichn);
	width_channel[ichn]->SetLineStyle(1);
	width_channel[ichn]->SetLineWidth(2);
      } else {  
	width_channel[ichn]->SetMarkerColor(1+ichn-8);
	width_channel[ichn]->SetLineColor(1+ichn-8);
	width_channel[ichn]->SetLineStyle(2);
	width_channel[ichn]->SetLineWidth(2);
      }
      if(ichn<3) width_channel[ichn]->SetMarkerStyle(ichn+2);
      else  width_channel[ichn]->SetMarkerStyle(16+ichn);

      if(ichn==0) {
	width_channel[ichn]->SetTitle(TString::Format("all ASU %s",angle.Data()));
	width_channel[ichn]->GetXaxis()->SetTitle("# SCA + 20x(ASU-1)");
	width_channel[ichn]->GetYaxis()->SetTitle("noise [ADC]");
	width_channel[ichn]->GetYaxis()->SetRangeUser(-10,10);
	width_channel[ichn]->Draw("ap");
      } else width_channel[ichn]->Draw("p");
      if(ichn<8) leg1->AddEntry(width_channel[ichn],TString::Format("channel=%i",channel[ichn]),"lp");
      else leg2->AddEntry(width_channel[ichn],TString::Format("channel=%i",channel[ichn]),"lp");
      
    }//ichn
    leg1->Draw();
    leg2->Draw();

    canvas->Print(TString::Format("pedestal_allasu_%s.png",angle.Data()));
 
  

   return 0;
 
}
