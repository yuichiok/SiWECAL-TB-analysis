#include "analysis.h"


 
void analysis(TString run, int layer) {

  gROOT->Reset();
  //SetIrlesStyle();
  //  gROOT->LoadMacro("Labels.C");
  gStyle->SetOptFit(0); 
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.2);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetMarkerSize(1.2);

  TFile *file = new TFile(TString::Format("plots/layer_%i_%s.root",layer,run.Data()) , "RECREATE");
  pedanalysis(file,run, layer);
  mipanalysis(file,run, layer);
  retriggeranalysis(file,run, layer);
  file->Close();

}
