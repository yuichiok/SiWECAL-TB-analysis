#include "TROOT.h"
#include "TFile.h"
#include "../../style/Style.C"
#include "TGraphErrors.h"
#include "TMath.h"

void TungstenProfile(){
 


  gROOT->Reset();
  SetIrlesStyle();
  gStyle->SetOptFit(0); 
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.4);
  gStyle->SetTitleY(0.9);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleFontSize(0.03);
  gStyle->SetMarkerSize(1.2);


  TH1F *conf1 = new TH1F("conf1","conf1",7,0.5,7.5);
  TH1F *conf2 = new TH1F("conf2","conf2",7,0.5,7.5);
  TH1F *conf3 = new TH1F("conf3","conf3",7,0.5,7.5);

  conf1->SetBinContent(1,2*0.56);
  conf1->SetBinContent(2,2*0.56);
  conf1->SetBinContent(3,2*0.56);
  conf1->SetBinContent(4,2*0.56);
  conf1->SetBinContent(5,4*0.56);
  conf1->SetBinContent(6,4*0.56);
  conf1->SetBinContent(7,6*0.56);

  conf2->SetBinContent(1,4*0.56);
  conf2->SetBinContent(2,2*0.56);
  conf2->SetBinContent(3,2*0.56);
  conf2->SetBinContent(4,4*0.56);
  conf2->SetBinContent(5,4*0.56);
  conf2->SetBinContent(6,6*0.56);
  conf2->SetBinContent(7,6*0.56);
  
  conf3->SetBinContent(1,6*0.56);
  conf3->SetBinContent(2,2*0.56);
  conf3->SetBinContent(3,4*0.56);
  conf3->SetBinContent(4,4*0.56);
  conf3->SetBinContent(5,6*0.56);
  conf3->SetBinContent(6,6*0.56);
  conf3->SetBinContent(7,6*0.56);
  
    
  TCanvas *c_energy = new TCanvas("c_energy","c_energy",800,600);
  //c_energy->Divide(2,1);
  c_energy->cd(1);
  conf1->SetTitle("W Configuration per layer");
  conf1->GetXaxis()->SetTitle("layer number");
  conf1->GetYaxis()->SetTitle("W thickness [X_{0}]");
  conf1->GetYaxis()->SetRangeUser(0,6);
  conf1->GetXaxis()->SetRangeUser(0,8);
  conf1->SetLineWidth(2);
  conf1->Draw();

  conf2->SetLineColor(2);
  conf2->SetLineStyle(2);
  conf2->SetLineWidth(4);
  conf2->Draw("same");

  conf3->SetLineColor(4);
  conf3->SetLineStyle(4);
  conf3->SetLineWidth(4);
  conf3->Draw("same");

  TLegend *l_energy = new TLegend(0.4,0.7,0.8,0.8);
  l_energy->AddEntry(conf1,"configuration 1 (12.32 X_{0})","l");
  l_energy->AddEntry(conf2,"configuration 2 (15.68 X_{0})","l");
  l_energy->AddEntry(conf3,"configuration 3 (19.04 X_{0})","l");

  l_energy->SetFillColor(0);
  l_energy->SetLineColor(0);
  l_energy->SetShadowColor(0);

  l_energy->Draw();

  c_energy->Print("tungsten_configurations.eps");


  TH1F *acc_conf1 = new TH1F("acc_conf1","acc_conf1",7,0.5,7.5);
  TH1F *acc_conf2 = new TH1F("acc_conf2","acc_conf2",7,0.5,7.5);
  TH1F *acc_conf3 = new TH1F("acc_conf3","acc_conf3",7,0.5,7.5);

  acc_conf1->SetBinContent(1,2*0.56);
  acc_conf1->SetBinContent(2,4*0.56);
  acc_conf1->SetBinContent(3,6*0.56);
  acc_conf1->SetBinContent(4,8*0.56);
  acc_conf1->SetBinContent(5,12*0.56);
  acc_conf1->SetBinContent(6,16*0.56);
  acc_conf1->SetBinContent(7,22*0.56);

  acc_conf2->SetBinContent(1,4*0.56);
  acc_conf2->SetBinContent(2,6*0.56);
  acc_conf2->SetBinContent(3,8*0.56);
  acc_conf2->SetBinContent(4,12*0.56);
  acc_conf2->SetBinContent(5,16*0.56);
  acc_conf2->SetBinContent(6,22*0.56);
  acc_conf2->SetBinContent(7,28*0.56);
  
  acc_conf3->SetBinContent(1,6*0.56);
  acc_conf3->SetBinContent(2,8*0.56);
  acc_conf3->SetBinContent(3,12*0.56);
  acc_conf3->SetBinContent(4,16*0.56);
  acc_conf3->SetBinContent(5,22*0.56);
  acc_conf3->SetBinContent(6,28*0.56);
  acc_conf3->SetBinContent(7,34*0.56);
  
    
  TCanvas *c_energy2 = new TCanvas("c_energy2","c_energy2",800,600);
  //c_energy->Divide(2,1);
  c_energy2->cd(1);
  acc_conf1->SetTitle("Accumulated W Configuration");
  acc_conf1->GetXaxis()->SetTitle("layer number");
  acc_conf1->GetYaxis()->SetTitle("W thickness [X_{0}]");
  acc_conf1->GetYaxis()->SetRangeUser(0,30);
  acc_conf1->GetXaxis()->SetRangeUser(0,8);

  acc_conf1->SetLineWidth(2);
  acc_conf1->Draw();

  acc_conf2->SetLineColor(2);
  acc_conf2->SetLineStyle(2);
  acc_conf2->SetLineWidth(4);
  acc_conf2->Draw("same");

  acc_conf3->SetLineColor(4);
  acc_conf3->SetLineStyle(4);
  acc_conf3->SetLineWidth(4);
  acc_conf3->Draw("same");

  TLegend *l_energy2 = new TLegend(0.2,0.7,0.6,0.8);
  l_energy2->AddEntry(acc_conf1,"configuration 1 (12.32 X_{0})","l");
  l_energy2->AddEntry(acc_conf2,"configuration 2 (15.68 X_{0})","l");
  l_energy2->AddEntry(acc_conf3,"configuration 3 (19.04 X_{0})","l");

  l_energy2->SetFillColor(0);
  l_energy2->SetLineColor(0);
  l_energy2->SetShadowColor(0);

  l_energy2->Draw();

  c_energy2->Print("tungsten_acc_configurations.eps");


}
