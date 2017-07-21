#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

#include "TFile.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TH1.h"
#include "TH2F.h"
#include "TImage.h"

#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"


using namespace std;

int Analysis_bcid(TString dif="1_1_2", bool ratio=true) {
  gROOT->Reset();

  gStyle->SetOptFit(0); 
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleBorderSize(0);
  //  gStyle->SetTitleX(0.2);
  //gStyle->SetTitleY(0.9);
  //  // gStyle->SetTitleStyle(0);

  gROOT->ForceStyle();

  for(int igrid=0;igrid<4;igrid++) {
    TString grid="";
    if(igrid == 1) grid = "00";
    if(igrid == 2) grid = "41";
    if(igrid == 3) grid = "80";


    for(int ibcid=0; ibcid<5; ibcid++) {
      TString bcid="5";
      if(ibcid==1) bcid="10";
      if(ibcid==2) bcid="15";
      if(ibcid==3) bcid="30";
      if(ibcid==4) bcid="60";

  
      if(grid!="") grid="_grid"+grid+"_bcidTh"+bcid;
      else  grid="_bcidTh"+bcid;

      TCanvas *pedestal = new TCanvas("pedestal","pedestal",1000,1000);
      TCanvas *pedestal_width = new TCanvas("pedestal_width","pedestal_width",1000,1000);
      TCanvas *pedestal_npeaks = new TCanvas("pedestal_npeaks","pedestal_npeaks",1000,1000);
      TCanvas *pedestal_tagged_npeaks = new TCanvas("pedestal_tagged_npeaks","pedestal_tagged_npeaks",1000,1000);
  
      TH2F *tmp_pedestal_total[15];
      TH2F *tmp_pedestal_width_total[15];
      TH2F *tmp_pedestal_npeaks_total[15];
  
      TH2F *tmp_pedestal[15];
      TH2F *tmp_pedestal_width[15];
      TH2F *tmp_pedestal_npeaks[15];
      TH2F *tmp_pedestal_tagged_npeaks[15];
  
      TString file="/home/irles/cernbox/TB2017/TBdata/MIPscan/pedestal/rootfiles/Pedestal_summary_dif_"+dif+"_bcidTh15.root";
      std::cout<<"Opening file with full run: "<<file<<std::endl;
      TFile *f_total = new TFile(file);
  
      for(int isca=0; isca<15;isca++) {
	tmp_pedestal_total[isca]= (TH2F*)f_total->Get(TString::Format("pedestal_map_sca%i",isca));
	tmp_pedestal_width_total[isca]= (TH2F*)f_total->Get(TString::Format("pedestal_width_map_sca%i",isca));
	tmp_pedestal_npeaks_total[isca]= (TH2F*)f_total->Get(TString::Format("pedestal_npeaks_map_sca%i",isca));
      }
      // f_total->Close();
  
      file="/home/irles/cernbox/TB2017/TBdata/MIPscan/pedestal/rootfiles/Pedestal_summary_dif_"+dif+grid+".root";
      std::cout<<"Opening file for single run: "<<file<<std::endl;
      TFile *f = new TFile(file);
  
      for(int isca=0; isca<15;isca++) {
	tmp_pedestal[isca]= (TH2F*)f->Get(TString::Format("pedestal_map_sca%i",isca));
	tmp_pedestal_width[isca]= (TH2F*)f->Get(TString::Format("pedestal_width_map_sca%i",isca));
	tmp_pedestal_npeaks[isca]= (TH2F*)f->Get(TString::Format("pedestal_npeaks_map_sca%i",isca));
	tmp_pedestal_tagged_npeaks[isca]= (TH2F*)f->Get(TString::Format("pedestal_tagged_npeaks_map_sca%i",isca));
      }
  
      //f->Close();
      TString title="map";
      if(ratio==true) {
	title="ratiomap";
	for(int isca=0; isca<15;isca++) {
	  // tmp_pedestal_total[isca]->Draw("colz");
	  tmp_pedestal[isca]->Divide(tmp_pedestal_total[isca]);
	  tmp_pedestal_width[isca]->Divide(tmp_pedestal_width_total[isca]);
	  tmp_pedestal_npeaks[isca]->Divide(tmp_pedestal_npeaks_total[isca]);
	  tmp_pedestal_tagged_npeaks[isca]->Divide(tmp_pedestal_npeaks_total[isca]);
	}
      }
  
      // good pedestal events (not tagged events)
      pedestal->Divide(4,4);
      for(int isca=0; isca<15; isca ++) {
	TString dif0=dif+TString::Format("_sca%i",isca);
	pedestal->cd(isca+1);
	tmp_pedestal[isca]->SetStats(kFALSE);
	tmp_pedestal[isca]->SetTitle("Ped-"+title+" "+grid+" dif="+dif);
	tmp_pedestal[isca]->GetXaxis()->SetTitle("x");
	tmp_pedestal[isca]->GetYaxis()->SetTitle("y");
	if(ratio==true) tmp_pedestal[isca]->GetZaxis()->SetRangeUser(0.98,1.02);
	else tmp_pedestal[isca]->GetZaxis()->SetRangeUser(200,400);
	
	tmp_pedestal[isca]->Draw("colz");
      }
      pedestal->Print("plots/pedestal_"+title+"_dif_"+dif+grid+".png");
      
      // good pedestal_width events (not tagged events)
      pedestal_width->Divide(4,4);
      for(int isca=0; isca<15; isca ++) {
	TString dif0=dif+TString::Format("_sca%i",isca);
	pedestal_width->cd(isca+1);
	tmp_pedestal_width[isca]->SetStats(kFALSE);
	tmp_pedestal_width[isca]->SetTitle("Width-"+title+" "+grid+" dif="+dif);
	tmp_pedestal_width[isca]->GetXaxis()->SetTitle("x");
	tmp_pedestal_width[isca]->GetYaxis()->SetTitle("y");
	if(ratio==true) tmp_pedestal_width[isca]->GetZaxis()->SetRangeUser(0.8,1.2);
	else tmp_pedestal_width[isca]->GetZaxis()->SetRangeUser(0,4.5);
	tmp_pedestal_width[isca]->Draw("colz");
      }
      pedestal_width->Print("plots/pedestal_width_"+title+"_dif_"+dif+grid+".png");

      // good pedestal_npeaks events (not tagged events)
      pedestal_npeaks->Divide(4,4);
      for(int isca=0; isca<15; isca ++) {
	TString dif0=dif+TString::Format("_sca%i",isca);
	pedestal_npeaks->cd(isca+1);
	tmp_pedestal_npeaks[isca]->SetStats(kFALSE);
	tmp_pedestal_npeaks[isca]->SetTitle("NPeaks-"+title+" "+grid+" dif="+dif);
	tmp_pedestal_npeaks[isca]->GetXaxis()->SetTitle("x");
	tmp_pedestal_npeaks[isca]->GetYaxis()->SetTitle("y");
	tmp_pedestal_npeaks[isca]->GetZaxis()->SetRangeUser(0,3);
	
	tmp_pedestal_npeaks[isca]->Draw("colz");
      }
      pedestal_npeaks->Print("plots/pedestal_npeaks_"+title+"_dif_"+dif+grid+".png");

      // good pedestal_tagged_npeaks events (not tagged events)
      pedestal_tagged_npeaks->Divide(4,4);
      for(int isca=0; isca<15; isca ++) {
	TString dif0=dif+TString::Format("_sca%i",isca);
	pedestal_tagged_npeaks->cd(isca+1);
	tmp_pedestal_tagged_npeaks[isca]->SetStats(kFALSE);
	tmp_pedestal_tagged_npeaks[isca]->SetTitle("NPeaks(tagged)-"+title+" "+grid+" dif="+dif);
	tmp_pedestal_tagged_npeaks[isca]->GetXaxis()->SetTitle("x");
	tmp_pedestal_tagged_npeaks[isca]->GetYaxis()->SetTitle("y");
	tmp_pedestal_tagged_npeaks[isca]->GetZaxis()->SetRangeUser(0,3);
	
	tmp_pedestal_tagged_npeaks[isca]->Draw("colz");
      }
      pedestal_tagged_npeaks->Print("plots/pedestal_tagged_npeaks_"+title+"_dif_"+dif+grid+".png");
      
      f_total->Close();
      f->Close();
    }

  }
  
  return 1;
      
}
