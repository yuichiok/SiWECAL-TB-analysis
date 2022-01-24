#include "analysis.h"

 
void analysis(TString run="3GeVMIPscan", TString gain="high", bool pedestal=true, bool mip=true) {

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
  
  if(pedestal==true)  pedanalysis(run,gain);
  else {
    if(mip==true) mipanalysis_summary(run,gain);
    else s_n_analysis_summary(run,gain);
  }

  gSystem->Exit(0);
  return 0;
}
