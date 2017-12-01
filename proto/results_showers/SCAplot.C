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

void SCAplot(){
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0111); 
  gStyle->SetOptTitle(1);
  gStyle->SetTitleBorderSize(0);
  // gStyle->SetTitleX(0.2);
  //gStyle->SetTitleY(0.9);
  //gStyle->SetTitleStyle(0);


  TFile *file_1GeV= TFile::Open("SimpleDist_conf3_1GeV_grid20.root");
  TFile *file_2GeV= TFile::Open("SimpleDist_conf3_2GeV_grid20.root");
  TFile *file_3GeV= TFile::Open("SimpleDist_conf3_3GeV_grid20.root");
  TFile *file_4GeV= TFile::Open("SimpleDist_conf3_4GeV_grid20.root");
  TFile *file_5GeV= TFile::Open("SimpleDist_conf3_5GeV_grid20.root");

  TH1F * SCA_1GeV = (TH1F*)file_1GeV->Get("SCA_shower"); 
  TH1F * SCA_2GeV = (TH1F*)file_2GeV->Get("SCA_shower"); 
  TH1F * SCA_3GeV = (TH1F*)file_3GeV->Get("SCA_shower"); 
  TH1F * SCA_4GeV = (TH1F*)file_4GeV->Get("SCA_shower"); 
  TH1F * SCA_5GeV = (TH1F*)file_5GeV->Get("SCA_shower"); 


  TCanvas *plot = new TCanvas("plot","plot",1200,800);
  
  SCA_1GeV->SetLineWidth(1);
  SCA_1GeV->SetLineColor(1);
  //SCA_1GeV->GetYaxis()->SetRangeUser(0,50000);
  SCA_1GeV->GetYaxis()->SetTitle("#entries / total #entries");
   SCA_1GeV->SetTitle("SCA distribution (e^{+} beam absorber configuration 3)");
  SCA_1GeV->DrawNormalized("h");

  SCA_2GeV->SetLineWidth(2);
  SCA_2GeV->SetLineColor(2);
  SCA_2GeV->DrawNormalized("hsame");

  SCA_3GeV->SetLineWidth(2);
  SCA_3GeV->SetLineColor(3);
  SCA_3GeV->DrawNormalized("hsame");

  SCA_4GeV->SetLineWidth(2);
  SCA_4GeV->SetLineColor(4);
  SCA_4GeV->DrawNormalized("hsame");

  SCA_5GeV->SetLineWidth(2);
  SCA_5GeV->SetLineColor(6);
  SCA_5GeV->DrawNormalized("hsame");

  TLegend *l_energy = new TLegend(0.5,0.65,0.93,0.9);
  
  l_energy->AddEntry(SCA_1GeV,"1 GeV beam","l");
  l_energy->AddEntry(SCA_2GeV,"2 GeV beam","l");
  l_energy->AddEntry(SCA_3GeV,"3 GeV beam","l");
  l_energy->AddEntry(SCA_4GeV,"4 GeV beam","l");
  l_energy->AddEntry(SCA_5GeV,"5 GeV beam","l");

  l_energy->SetFillColor(0);
  l_energy->SetLineColor(0);
  l_energy->SetShadowColor(0);
  l_energy->Draw();
  



   
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
