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

void mipfitplot(){

  gROOT->Reset();
  SetIrlesStyle();
  gStyle->SetOptFit(0111); 
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.3);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleFontSize(0.03);
  gStyle->SetMarkerSize(1.2);


  TFile *_file0 = TFile::Open("eff_files/Signal_summary_eff0.15_hitsintrack5_hitouttrack3_bcidcut2850.root");

  TH1F * energy_distribution = (TH1F*)_file0->Get("energy_distribution"); 

  TH1F * first_iteration_peak1 = (TH1F*)_file0->Get("first_iteration_peak1"); 
  TH1F * first_iteration_peak2 = (TH1F*)_file0->Get("first_iteration_peak2"); 
  TH1F * first_iteration_peak3 = (TH1F*)_file0->Get("first_iteration_peak3"); 

  TH1F * first_peak_fit_histo = (TH1F*)_file0->Get("first_peak_fit_histo"); 
  TH1F * second_peak_fit_histo = (TH1F*)_file0->Get("second_peak_fit_histo"); 
  TH1F * third_peak_fit_histo = (TH1F*)_file0->Get("third_peak_fit_histo"); 


  TCanvas *plot = new TCanvas("plot","plot",1200,800);
  
  gPad->SetLogy();

  energy_distribution->SetLineWidth(1);
  energy_distribution->SetLineColor(1);
  energy_distribution->SetStats("nemr");
  energy_distribution->GetXaxis()->SetRangeUser(0,6);
  energy_distribution->GetYaxis()->SetRangeUser(1000,500000);
  energy_distribution->SetTitle("Single cell energy distribution for 3 GeV e^{+} beam w/o absorber");
  energy_distribution->Draw("h");

  first_iteration_peak1->SetLineWidth(2);
  first_iteration_peak1->SetLineStyle(1);
  first_iteration_peak1->SetLineColor(2);
  first_iteration_peak1->SetStats(0);
  first_iteration_peak1->Draw("lsame");

  gPad->Update();

  first_peak_fit_histo->Draw("samesaxis");
  gPad->Update();
  TPaveStats *st_histo1 = (TPaveStats*)first_peak_fit_histo->FindObject("stats");
  //gPad->Update();
  st_histo1->SetLineColor(2);
  st_histo1->SetShadowColor(0);
  st_histo1->SetX1NDC(0.5);
  st_histo1->SetX2NDC(0.7);
  st_histo1->SetY1NDC(0.78);
  st_histo1->SetY2NDC(0.9);
  st_histo1->SetTextColor(2);
  gPad->Modified();

 
  //-----------

  first_iteration_peak2->SetLineWidth(2);
  first_iteration_peak2->SetLineStyle(1);
  first_iteration_peak2->SetLineColor(8);
  first_iteration_peak2->SetStats(0);
  first_iteration_peak2->Draw("lsame");

  gPad->Update();

  second_peak_fit_histo->Draw("samesaxis");
  gPad->Update();
  TPaveStats *st_histo2 = (TPaveStats*)second_peak_fit_histo->FindObject("stats");
  //gPad->Update();
  st_histo2->SetLineColor(3);
  st_histo2->SetShadowColor(0);
  st_histo2->SetX1NDC(0.71);
  st_histo2->SetX2NDC(0.91);
  st_histo2->SetY1NDC(0.78);
  st_histo2->SetY2NDC(0.9);
  st_histo2->SetTextColor(8);
  gPad->Modified();

 //-----------

  first_iteration_peak3->SetLineWidth(2);
  first_iteration_peak3->SetLineStyle(1);
  first_iteration_peak3->SetLineColor(4);
  first_iteration_peak3->SetStats(0);
  first_iteration_peak3->Draw("lsame");

  gPad->Update();

  third_peak_fit_histo->Draw("samesaxis");
  gPad->Update();
  TPaveStats *st_histo3 = (TPaveStats*)third_peak_fit_histo->FindObject("stats");
  //gPad->Update();
  st_histo3->SetLineColor(4);
  st_histo3->SetShadowColor(0);
  st_histo3->SetX1NDC(0.71);
  st_histo3->SetX2NDC(0.91);
  st_histo3->SetY1NDC(0.64);
  st_histo3->SetY2NDC(0.76);
  st_histo3->SetTextColor(4);
  gPad->Modified();
  IRLESLabel(0.55,0.55,"",kGray+2);



   
}
// TFile**		Signal_summary.root	
//  TFile*		Signal_summary.root	
//   KEY: TH1F	energy_distribution;1	Single cell hit energy in 3GeV e^{+} beam
//   KEY: TH1F	first_peak_fit_histo;1	First MIP
//   KEY: TH1F	second_peak_fit_histo;1	Second MIP (first subtracted)
//   KEY: TH1F	third_peak_fit_histo;1	Third MIP (first+second subtracted)
//   KEY: TF1	FirstPeak_iter1;1	fpeak1_iter1
//   KEY: TF1	SecondPeak_iter1;1	fpeak2_iter1
//   KEY: TF1	ThirdPeak_iter1;1	fpeak3_iter1
//   KEY: TH1F	first_iteration_peak1;1	first_iteration_peak1
//   KEY: TH1F	first_iteration_peak2;1	first_iteration_peak2
//   KEY: TH1F	first_iteration_peak3;1	first_iteration_peak3
