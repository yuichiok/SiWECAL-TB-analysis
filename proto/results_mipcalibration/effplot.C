#include "TFile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TH1.h"
#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TString.h"
#include "../../style/Style.C"
#include "../../style/Labels.C"

using namespace std;

void effplot(){

  gROOT->Reset();
  SetIrlesStyle();
  //  gROOT->LoadMacro("Labels.C");
  gStyle->SetOptFit(0); 
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
 
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.1);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleFontSize(0.03);
  gStyle->SetMarkerSize(1.2);
  int bcid_max=-1;
  int nout=-1;

  // int nhits=4;
  //int nout=1000;

  for(int nhits=2; nhits<7; nhits++) {

    for(int ibcid_max=0; ibcid_max<1; ibcid_max++) {
      if(ibcid_max==0) bcid_max=2850;
      if(ibcid_max==1) bcid_max=10000;

      for(int inout=0; inout<3; inout++) {
	if(inout==0) nout=0;
	if(inout==1) nout=3;
	if(inout==2) nout=10;


  	TFile *_file0 = TFile::Open(TString::Format("eff_files/Signal_summary_eff0.15_hitsintrack%i_hitouttrack%i_bcidcut%i.root",nhits,nout,bcid_max));
	cout<<"Open: "<<TString::Format("eff_files/Signal_summary_eff0.15_hitsintrack%i_hitouttrack%i_bcidcut%i.root",nhits,nout,bcid_max)<<endl;
	TH2F* ineff[7];
	ineff[0]= (TH2F*)_file0->Get("mip_ineff_layer0_chip_chn"); 
	ineff[1]= (TH2F*)_file0->Get("mip_ineff_layer1_chip_chn");
	ineff[2]= (TH2F*)_file0->Get("mip_ineff_layer2_chip_chn");
	ineff[3]= (TH2F*)_file0->Get("mip_ineff_layer3_chip_chn");
	ineff[4]= (TH2F*)_file0->Get("mip_ineff_layer4_chip_chn");
	ineff[5]= (TH2F*)_file0->Get("mip_ineff_layer5_chip_chn");
	ineff[6]= (TH2F*)_file0->Get("mip_ineff_layer6_chip_chn");

	// // -------------------------------------

	// TH2F* ineff_full[7];
	// ineff_full[0]= (TH2F*)_file0->Get("mip_ineff_full_layer0_chip_chn"); 
	// ineff_full[1]= (TH2F*)_file0->Get("mip_ineff_full_layer1_chip_chn");
	// ineff_full[2]= (TH2F*)_file0->Get("mip_ineff_full_layer2_chip_chn");
	// ineff_full[3]= (TH2F*)_file0->Get("mip_ineff_full_layer3_chip_chn");
	// ineff_full[4]= (TH2F*)_file0->Get("mip_ineff_full_layer4_chip_chn");
	// ineff_full[5]= (TH2F*)_file0->Get("mip_ineff_full_layer5_chip_chn");
	// ineff_full[6]= (TH2F*)_file0->Get("mip_ineff_full_layer6_chip_chn");

	
	double x[7][1024];
	double y[7][1024];

	for(int ilayer=0; ilayer<7; ilayer++) {
	  int ibin=0;
	  for(int ichip=0; ichip<16; ichip++) {
	    for(int ichn=0; ichn<64; ichn++) {
	      x[ilayer][ibin]=ichip*200+ichn;
	      y[ilayer][ibin]=100.-ineff[ilayer]->GetBinContent(ichip+1,ichn+1);//+ineff_full[ilayer]->GetBinContent(ichip+1,ichn+1));
	      ibin++;
	    }
	  }
	}

	TGraph *g_eff[7];
	for(int ilayer=0; ilayer<7; ilayer++) {
	  g_eff[ilayer]= new TGraph(1024,x[ilayer],y[ilayer]);
	}

	//	-------------------------------------

	TH2F* ineff_full[7];
	ineff_full[0]= (TH2F*)_file0->Get("mip_ineff_full_layer0_chip_chn"); 
	ineff_full[1]= (TH2F*)_file0->Get("mip_ineff_full_layer1_chip_chn");
	ineff_full[2]= (TH2F*)_file0->Get("mip_ineff_full_layer2_chip_chn");
	ineff_full[3]= (TH2F*)_file0->Get("mip_ineff_full_layer3_chip_chn");
	ineff_full[4]= (TH2F*)_file0->Get("mip_ineff_full_layer4_chip_chn");
	ineff_full[5]= (TH2F*)_file0->Get("mip_ineff_full_layer5_chip_chn");
	ineff_full[6]= (TH2F*)_file0->Get("mip_ineff_full_layer6_chip_chn");

	double x_full[7][1024];
	double y_full[7][1024];
	double x_full_mean[7][16];
	double y_full_mean[7][16];

	for(int ilayer=0; ilayer<7; ilayer++) {
	  int ibin=0;
	  for(int ichip=0; ichip<16; ichip++) {
	    int ibin_chip_1=0;
	    int ibin_chip_2=0;
	    double sum_chip_1=0;
	    double sum_chip_2=0;
	    for(int ichn=0; ichn<64; ichn++) {
	      x_full[ilayer][ibin]=ichip*200+ichn;
	      y_full[ilayer][ibin]=100.-ineff_full[ilayer]->GetBinContent(ichip+1,ichn+1);
	      ibin++;
	      if(ineff_full[ilayer]->GetBinContent(ichip+1,ichn+1)>0) {
		ibin_chip_1++;
		sum_chip_1+=ineff_full[ilayer]->GetBinContent(ichip+1,ichn+1);
	      }
	      if(ineff[ilayer]->GetBinContent(ichip+1,ichn+1)>0) {
		ibin_chip_2++;
		sum_chip_2+=ineff[ilayer]->GetBinContent(ichip+1,ichn+1);
	      }
	    }

	    double sum =0.;
	    //if(ibin_chip_1>0 ) sum+=sum_chip_1/ibin_chip_1;
	    if(ibin_chip_2>0 ) sum+=sum_chip_2/ibin_chip_2;

	    y_full_mean[ilayer][ichip]=100.-sum;
	    x_full_mean[ilayer][ichip]=ichip;
	    
	  }
	}

	TGraph *g_eff_full[7];
	for(int ilayer=0; ilayer<7; ilayer++) {
	  g_eff_full[ilayer]= new TGraph(1024,x_full[ilayer],y_full[ilayer]);
	}

	TGraph *g_eff_full_mean[7];
	for(int ilayer=0; ilayer<7; ilayer++) {
	  g_eff_full_mean[ilayer]= new TGraph(16,x_full_mean[ilayer],y_full_mean[ilayer]);
	}
  

	TCanvas *canvas = new TCanvas(TString::Format("efficiency_nhits%i",nhits),TString::Format("efficiency_nhits%i",nhits),2000,600);
	canvas->Divide(3,1);
	// gPad->SetLogy();
	canvas->cd(1);
	TLegend *leg=new TLegend(0.75,0.65,0.9,0.9);//0.75,0.25,0.9,0.5);
	for(int ilayer=0; ilayer<7; ilayer++) {
	  g_eff[ilayer]->SetMarkerStyle(20+ilayer);
	  g_eff[ilayer]->SetMarkerColor(ilayer+2);

	  if(ilayer==0) {
	    g_eff[ilayer]->SetTitle(TString::Format("Hit detection efficiency for tracks with at least %i hits",nhits));
	    g_eff[ilayer]->SetTitle(TString::Format("Hit detection efficiency for tracks with at least %i hits",nhits));
	    g_eff[ilayer]->GetYaxis()->SetRangeUser(90,110);
	    g_eff[ilayer]->GetYaxis()->SetTitle("%");
	    g_eff[ilayer]->GetXaxis()->SetTitle("200*chip+channel");
	    g_eff[ilayer]->Draw("ap");
	  }
	  else g_eff[ilayer]->Draw("p");
	  leg->AddEntry(g_eff[ilayer],TString::Format("Layer %i",ilayer),"p");
	}
	leg->SetFillColor(0);
	leg->SetLineColor(0);
	leg->SetShadowColor(0);
	leg->Draw();   
	IRLESLabel(0.2,0.2,"",kGray+2);

  
	canvas->cd(2);
	TLegend *leg2=new TLegend(0.75,0.65,0.9,0.9);

	for(int ilayer=0; ilayer<7; ilayer++) {
	  g_eff_full[ilayer]->SetMarkerStyle(20+ilayer);
	  g_eff_full[ilayer]->SetMarkerColor(ilayer+2);
	  g_eff_full[ilayer]->SetLineColor(ilayer+2);
	  if(ilayer==0) {
	    g_eff_full[ilayer]->SetTitle(TString::Format("Hit detection efficiency for tracks with at least %i hits (full buffers)",nhits));
	    // g_eff_full[ilayer]->GetYaxis()->SetRangeUser(-0.1,1);
	    g_eff_full[ilayer]->GetYaxis()->SetRangeUser(90,110);
	    g_eff_full[ilayer]->GetYaxis()->SetTitle("%");
	    g_eff_full[ilayer]->GetXaxis()->SetTitle("200*chip+channel");
	    g_eff_full[ilayer]->Draw("ap");
	  }
	  else g_eff_full[ilayer]->Draw("p");

	  leg2->AddEntry(g_eff[ilayer],TString::Format("Layer %i",ilayer),"p");
	}

	leg2->SetFillColor(0);
	leg2->SetLineColor(0);
	leg2->SetShadowColor(0);
	leg2->Draw();

	IRLESLabel(0.2,0.2,"",kGray+2);

	  canvas->cd(3);
	TLegend *leg3=new TLegend(0.75,0.65,0.9,0.9);

	for(int ilayer=0; ilayer<7; ilayer++) {
	  g_eff_full_mean[ilayer]->SetMarkerStyle(20+ilayer);
	  g_eff_full_mean[ilayer]->SetMarkerColor(ilayer+2);
	  g_eff_full_mean[ilayer]->SetLineColor(ilayer+2);
	  if(ilayer==0) {
	    g_eff_full_mean[ilayer]->SetTitle(TString::Format("Average hit detection efficiency for tracks with at least %i hits",nhits));
	    // g_eff_full_mean[ilayer]->GetYaxis()->SetRangeUser(-0.1,1);
	    g_eff_full_mean[ilayer]->GetYaxis()->SetRangeUser(90,110);
	    g_eff_full_mean[ilayer]->GetYaxis()->SetTitle("%");
	    g_eff_full_mean[ilayer]->GetXaxis()->SetTitle("chip");
	    g_eff_full_mean[ilayer]->Draw("ap");
	  }
	  else g_eff_full_mean[ilayer]->Draw("p");

	  leg3->AddEntry(g_eff[ilayer],TString::Format("Layer %i",ilayer),"p");

	}

	leg3->SetFillColor(0);
	leg3->SetLineColor(0);
	leg3->SetShadowColor(0);
	leg3->Draw();
	
	IRLESLabel(0.2,0.2,"",kGray+2);
	
	canvas->Print(TString::Format("eff_files/MIPefficiency_hitsintrack%i_hitouttrack%i_bcidcut%i.eps",nhits,nout,bcid_max));

      }//nout
    }//bcid
  }//nhits

}
  
