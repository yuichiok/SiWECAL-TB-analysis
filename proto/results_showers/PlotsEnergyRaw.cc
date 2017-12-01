#include "TROOT.h"
#include "TFile.h"
#include "../../style/Style.C"
#include "TGraphErrors.h"
#include "TMath.h"

void PlotsEnergyRaw(){
 


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
  TString conf="1";
  TString energy_s="4";
  
  
  TString s_file="zmlt1/conf"+conf+"_"+grid+"_"+energy_s+"GeV_mipcut0.5_showers.root";
  std::cout<<"Opening file: "<<s_file<<std::endl;
  TFile *file = new TFile(s_file);
  
  TH1F *energy = (TH1F*)file->Get("energy");
  
  energy->Rebin(4);
    
  TCanvas *c_energy = new TCanvas("c_energy","c_energy",800,600);
  //c_energy->Divide(2,1);
  c_energy->cd(1);
  // energy->SetTitle("SiW-ECAL (wafer-3), e^{+} beam");
  energy->GetXaxis()->SetTitle("E^{raw}/MIP = #sum_{i=layers} #sum_{j=cells} #omega_{i} E_{i,j}^{raw}");
  energy->GetYaxis()->SetTitle("# entries");
  energy->Draw();
  
  TLegend *l_energy = new TLegend(0.6,0.7,0.92,0.9);
  l_energy->SetHeader("wafer 3, W-configuration "+conf);  
  l_energy->AddEntry(energy,"e^{+} "+energy_s+" GeV","l");
  l_energy->Draw();

  c_energy->Print("shower_plots/raw_"+energy_s+"GeV_"+grid+"_"+conf+".eps");


}
