#include "TROOT.h"
#include "TFile.h"
#include "Style.C"

void Plots(){


  gROOT->Reset();

  //gROOT->LoadMacro("Style.C");
  //gROOT->LoadMacro("Labels.C");
  SetIrlesStyle();
  gStyle->SetOptFit(0); 
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.4);
  gStyle->SetTitleY(0.9);
    //gStyle->SetTitleColor(0);
  gStyle->SetTitleStyle(0);
  //gStyle->SetTitleFont(44);
  gStyle->SetTitleFontSize(0.03);

  TString grid="grid24";
  TString conf="conf1";

  for(int icases=0; icases<4; icases++) {

    if(icases==1) grid="grid20";
    if(icases==2) conf="conf2";
    if(icases==3) conf="conf3";


  
  TString s_file=conf+"_"+grid+"_2GeV_mipcut0.5_showers.root";
  std::cout<<"Opening file: "<<s_file<<std::endl;
  TFile *file = new TFile(s_file);
  
  TH1F *energy = (TH1F*)file->Get("energy");
  TH1F *energy_center = (TH1F*)file->Get("energy_center");

  energy->Rebin(4);
  energy_center->Rebin(4);

  //  file->Close();

  TString s_file_2=conf+"_"+grid+"_4GeV_mipcut0.5_showers.root";

  std::cout<<"Opening file: "<<s_file_2<<std::endl;
  TFile *file_2 = new TFile(s_file_2);
  
  TH1F *energy_2 = (TH1F*)file_2->Get("energy");
  TH1F *energy_center_2 = (TH1F*)file_2->Get("energy_center");

   energy_2->Rebin(4);
   energy_center_2->Rebin(4);

  TCanvas *c_energy = new TCanvas("c_energy","c_energy",800,600);
  //c_energy->Divide(2,1);
  c_energy->cd(1);
  // energy->SetTitle("SiW-ECAL (wafer-3), e^{+} beam");
  energy->GetXaxis()->SetTitle("E^{raw}/MIP = #sum_{i=layers} #sum_{j=cells} #omega_{i} E_{i,j}^{raw}");
  energy->GetYaxis()->SetTitle("# entries");
  energy->Draw();
  
  energy_center->SetLineColor(4);
  energy_center->Draw("same");

  energy_2->SetLineStyle(2);
  energy_2->Draw("same");

  energy_center_2->SetLineColor(4);
  energy_center_2->SetLineStyle(2);
  energy_center_2->Draw("same");

  TLegend *l_energy = new TLegend(0.4,0.7,0.9,0.9);
  l_energy->SetHeader("SiW-ECAL: wafer 4, W-configuration 1");
  if(icases==1)   l_energy->SetHeader("SiW-ECAL: wafer 3, W-configuration 1");
  if(icases==2)   l_energy->SetHeader("SiW-ECAL: wafer 3, W-configuration 2");
  if(icases==3)   l_energy->SetHeader("SiW-ECAL: wafer 3, W-configuration 3");

  l_energy->AddEntry(energy,"e^{+} 2 GeV beam","l");
  l_energy->AddEntry(energy_center,"e^{+} 2 GeV beam, contained showers","l");
  l_energy->AddEntry(energy_2,"e^{+} 4 GeV beam","l");
  l_energy->AddEntry(energy_center_2,"e^{+} 4 GeV, contained showers","l");
  l_energy->SetFillColor(0);
  l_energy->SetLineColor(0);
  l_energy->SetShadowColor(0);
  l_energy->Draw();
  
  c_energy->Print("energy_"+grid+"_"+conf+".eps");

  }
  
}
