#include "TFile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TH1.h"
#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TString.h"

using namespace std;

void mipfitplot(){
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0111); 
  gStyle->SetOptTitle(1);
  gStyle->SetTitleBorderSize(0);
  // gStyle->SetTitleX(0.2);
  //gStyle->SetTitleY(0.9);
  //gStyle->SetTitleStyle(0);


  TFile *_file0 = TFile::Open("Signal_summarytest_old.root");

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
  energy_distribution->GetYaxis()->SetRangeUser(100,25000);
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
  st_histo1->SetX1NDC(0.35);
  st_histo1->SetX2NDC(0.6);
  st_histo1->SetY1NDC(0.7);
  st_histo1->SetY2NDC(0.85);
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
  st_histo2->SetX1NDC(0.62);
  st_histo2->SetX2NDC(0.87);
  st_histo2->SetY1NDC(0.7);
  st_histo2->SetY2NDC(0.85);
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
  st_histo3->SetX1NDC(0.62);
  st_histo3->SetX2NDC(0.87);
  st_histo3->SetY1NDC(0.5);
  st_histo3->SetY2NDC(0.65);
  st_histo3->SetTextColor(4);
  gPad->Modified();



   
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
