#include "TROOT.h"
#include "TFile.h"
#include "../../style/Style.C"
#include "TGraphErrors.h"
#include "TLatex.h"

void PlotsBarycenter(){


  gROOT->Reset();
  SetIrlesStyle();
  gStyle->SetOptFit(1); 
  gStyle->SetOptStat(1);
  gStyle->SetOptTitle(0);
  
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.4);
  gStyle->SetTitleY(0.9);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleFontSize(0.03);
  gStyle->SetMarkerSize(1.);

  TString grid="grid24";
  TString conf="conf1";

  for(int icases=0; icases<4; icases++) {

    if(icases>0) grid="grid20";
    if(icases==2) conf="conf2";
    if(icases==3) conf="conf3";

    //energy_1 GeV
    TString s_file="zmlt1/"+conf+"_"+grid+"_1GeV_mipcut0.5_showers.root";
    std::cout<<"Opening file: "<<s_file<<std::endl;
    TFile *file_1 = new TFile(s_file);
    TH1F *xbarycenter = (TH1F*)file_1->Get("x-barycenter");
    TH1F *ybarycenter = (TH1F*)file_1->Get("y-barycenter");
    TH1F *zbarycenter = (TH1F*)file_1->Get("z-barycenter");

    TH2F *xybarycenter = (TH2F*)file_1->Get("Energy xy-barycenter");
    TH2F *xzbarycenter = (TH2F*)file_1->Get("Energy xz-barycenter");
    TH2F *yzbarycenter = (TH2F*)file_1->Get("Energy yz-barycenter");

    gStyle->SetPadRightMargin(0.1);

    
    TCanvas *c_energy = new TCanvas("c_energy","c_energy",1600,800);
    c_energy->Divide(3,2);
    c_energy->cd(1);
    xbarycenter->GetXaxis()->SetTitle("#bar{x}");//#sum_{i=all cells}x_{i}#omega_{i}E_{i}^{raw}}{#sum_{i=all cells}#omega_{i}E_{i}^{raw}}");
    xbarycenter->GetYaxis()->SetTitle("# entries");
    xbarycenter->GetYaxis()->SetRangeUser(0,1.25*xbarycenter->GetMaximum());
    xbarycenter->Draw("");

    c_energy->cd(2);
    ybarycenter->GetXaxis()->SetTitle("#bar{y}");
    ybarycenter->GetYaxis()->SetTitle("# entries");
    ybarycenter->GetYaxis()->SetRangeUser(0,1.25*ybarycenter->GetMaximum());
    ybarycenter->Draw("");

    c_energy->cd(3);
    zbarycenter->GetXaxis()->SetTitle("#bar{z}");
    zbarycenter->GetYaxis()->SetTitle("# entries");
    zbarycenter->GetYaxis()->SetRangeUser(0,1.25*zbarycenter->GetMaximum());
    zbarycenter->Draw("");

    c_energy->cd(4);
    xybarycenter->SetStats(kFALSE);
    xybarycenter->GetXaxis()->SetTitle("#bar{x}");
    xybarycenter->GetYaxis()->SetTitle("#bar{y}");
    xybarycenter->Draw("colz");
  
    c_energy->cd(5);
    xzbarycenter->SetStats(kFALSE);
    xzbarycenter->GetXaxis()->SetTitle("#bar{x}");
    xzbarycenter->GetYaxis()->SetTitle("#bar{z}");
    xzbarycenter->Draw("colz");

    c_energy->cd(6);
    yzbarycenter->SetStats(kFALSE);
    yzbarycenter->GetXaxis()->SetTitle("#bar{y}");
    yzbarycenter->GetYaxis()->SetTitle("#bar{z}");
    yzbarycenter->Draw("colz");
   
  
    c_energy->Print("zmlt1/barycenter_plots_1GeV__"+grid+"_"+conf+".eps");

  }



  


  
}
