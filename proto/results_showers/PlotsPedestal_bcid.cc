#include "TROOT.h"
#include "TFile.h"
#include "../../style/Style.C"
#include "../../style/Labels.C"
#include "TGraphErrors.h"
#include "TMath.h"


//  for(int iz=1; iz<2; iz++) {
//    TString zm="zm0.5";
//     if(iz==1) zm="zm1.0";
//    if(iz==2) zm="zm1.0";

void PlotsPedestal_bcid(){
 
  gROOT->Reset();
  SetIrlesStyle();
  gStyle->SetOptFit(0); 
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.4);
  gStyle->SetTitleY(0.9);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleFontSize(0.03);
  gStyle->SetMarkerSize(1.2);

  TString grid="grid20";
  TString conf="conf2";
  
  int islabs=6;
  TString nslabs_pedestal=TString::Format("pedestal_nslabs%i",islabs);

  TString bcid="bcidmax2850";
  
  TString energy_string;

  for(int ienergy=0; ienergy<6; ienergy++) {
    if(ienergy==0 ) energy_string="1GeV";
    if(ienergy==1 ) energy_string="2GeV";
    if(ienergy==2 ) energy_string="3GeV";
    if(ienergy==3 ) energy_string="4GeV";
    if(ienergy==4 ) energy_string="5GeV";
    if(ienergy==5 ) energy_string="5.8GeV";
  
  TString s_file=nslabs_pedestal+"_"+bcid+"/"+conf+"_"+grid+"_"+energy_string+"_mipcut0.5_showers.root";
  std::cout<<"Opening file: "<<s_file<<std::endl;
  TFile *file = new TFile(s_file);
  
    double max_x=0;
    for(int ilayer=0; ilayer<7; ilayer++) {

      TGraphErrors *pedestal_mean[15];
      for(int isca=0;isca<4; isca++) {
	TString histostring = TString::Format("pedestal_vs_prevbcid_z%i_chip12_sca%i",ilayer,isca);
	cout<<histostring<<endl;
	TH2F *temp= (TH2F*)file->Get(histostring);
	double x[50];
	double y[50];
	double ex[50];
	double ey[50];

	int nbins =0 ;
	for(int iy=1; iy<50; iy++) {
	  TH1D * proj_temp = temp->ProjectionX("proj_temp",iy,iy+1);  
	  x[iy-1]=20*(iy-1);
	  ex[iy-1]=10;
	  if(proj_temp->GetEntries()>50) {
	    y[iy-1]=100*proj_temp->GetMean();
	    ey[iy-1]=100*proj_temp->GetRMS();
	  } else {
	    y[iy-1]=0;
	    ey[iy-1]=0;
	  }
	  if(proj_temp->GetEntries()>0) nbins++;
	}

	pedestal_mean[isca] = new TGraphErrors(49,x,y,ex,ey);
    }
      
  
    
  TCanvas *c_energy = new TCanvas(TString::Format("c_energy_ilayer%i",ilayer),TString::Format("c_energy_ilayer%i",ilayer),800,600);

    
    c_energy->cd(1);
    pedestal_mean[0]->GetXaxis()->SetTitle("BCID-BCID_{previous}");
    pedestal_mean[0]->GetYaxis()->SetTitle("average pedestal deviation [%MIP]");
    pedestal_mean[0]->GetYaxis()->SetRangeUser(-20,20);
    pedestal_mean[0]->SetLineWidth(2);
    pedestal_mean[0]->SetLineColor(1);
    pedestal_mean[0]->SetMarkerStyle(20);
    pedestal_mean[0]->SetMarkerColor(1);
    pedestal_mean[0]->Draw("alp");

    for(int isca=1;isca<4; isca++) {
      pedestal_mean[isca]->SetLineWidth(1);
      pedestal_mean[isca]->SetLineColor(isca+1);
      pedestal_mean[isca]->SetMarkerStyle(isca+21);
      pedestal_mean[isca]->SetMarkerColor(isca+1);
      pedestal_mean[isca]->Draw("lp");
    }
    
    
  TLegend *l_1 = new TLegend(0.75,0.7,0.9,0.9);
  // l_1->SetHeader(TString::Format("SiW-ECAL: wafer 3, W-configuration 2, layer %i",ilayer));

  for(int isca=0;isca<4; isca++) {
    l_1->AddEntry(pedestal_mean[isca],TString::Format("SCA=%i",isca),"lp");
  }

  l_1->SetFillColor(0);
  l_1->SetLineColor(0);
  l_1->SetShadowColor(0);
  l_1->Draw();
  IRLESLabel(0.2,0.88,"",kGray+2);

  c_energy->Print("pedestal_plots/pedestal_deviation_previousbcid_energy"+energy_string+"_nslabs6_"+bcid+"_"+conf+"_"+grid+"_mipcut0.5_layer"+TString::Format("%i.eps",ilayer));

   }

}

}
