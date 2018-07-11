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
#include "../../../style/Style.C"
#include "../../../style/Labels.C"

using namespace std;


int Noise(){


  SetIrlesStyle();
  gStyle->SetOptFit(0111); 
  gStyle->SetOptStat(1110);

  gStyle->SetOptTitle(1);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.2);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleStyle(0);
  gStyle->SetMarkerSize(0.8);


  //open file
  TString filename ="../results_pedestal/Pedestal_dif_1_1_1.root";
  TFile *f = new TFile(filename);
  
  TCanvas *canvas = new TCanvas("canvas","canvas",1400,1400);
  canvas->Divide(4,4);

  // do analysis of chip 0, ichn 15
  Int_t ichip = 6;
  Int_t ichn = 15;
  cout<< "Analysis of chip: "<<ichip<< " channel:"<<ichn<<endl;

  Double_t pedestal[15];
  Double_t noise[15];
  Int_t entries[15]; 
  Int_t nsca=0; //total scas with enough data
  Double_t sca[15];

  cout<< "isca, entries, pedestal, noise " <<endl;

  for(int isca=0; isca<15; isca ++ ) {

    pedestal[isca]=-1;
    noise[isca]=-1;
    entries[isca]=-1;
    sca[isca]=isca;
    
    TH1F *hpedestal_good  = (TH1F*)f->Get(TString::Format("ped_chip%i_chn%i_sca%i",ichip,ichn,isca));
    
    entries[isca]=hpedestal_good->GetEntries();
    if(entries[isca]>5) {
      pedestal[isca]=hpedestal_good->GetMean();
      noise[isca]=hpedestal_good->GetRMS();
      nsca++;
    }

    canvas->cd(1+isca);
    hpedestal_good->GetXaxis()->SetRangeUser(260,340);
    hpedestal_good->Draw();


    
    cout<< isca << " " << entries[isca] << " " << pedestal[isca] << " " << noise[isca] << endl;
  }  // end analysis

         
  canvas->Print(TString::Format("pedestal_histograms_chip%i_chn%i.png",ichip,ichn));

  // create graphs
  TGraph *pedestal_mean = new TGraph(nsca,sca,pedestal);

  
  // draw the results
  TCanvas *canvas_summary = new TCanvas("canvas_summary","canvas_summary",800,600);
  canvas_summary->cd(1);
  pedestal_mean->SetMarkerColor(1);
  pedestal_mean->SetMarkerStyle(22);
  pedestal_mean->SetLineColor(1);
  pedestal_mean->SetLineStyle(1);
  pedestal_mean->SetLineWidth(2);
  pedestal_mean->SetTitle(TString::Format("chip%i chn%i",ichip,ichn));

  pedestal_mean->GetXaxis()->SetTitle("# SCA");
  pedestal_mean->GetYaxis()->SetTitle("pedestal [ADC]");
  pedestal_mean->GetYaxis()->SetRangeUser(200,400);
  pedestal_mean->Draw("alp");
    
  canvas_summary->Print(TString::Format("pedestal_mean_chip%i_chn%i.png",ichip,ichn));

  return 0;

}
