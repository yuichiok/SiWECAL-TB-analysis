#include "analysis.h"


 
void analysis(TString run, int slboard) {

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

  int cob=0;
  if(slboard==8 || slboard==12) cob=1;
  triggers(run, slboard,cob);

  //  TFile *file = new TFile(TString::Format("plots/slboard_%i_%s.root",slboard,run.Data()) , "RECREATE");  
  //mipanalysis(file,run, slboard,cob);
  //  file->Close();

}
