#include "TROOT.h"
#include "TFile.h"
#include "../../style/Style.C"
#include "../../style/Labels.C"
#include "TGraphErrors.h"
#include "TMath.h"


//  for(int iz=1; iz<2; iz++) {
//    TString zm="zm0.5";
//     if(iz==1) zm="zm1.0";
//    if(iz==2) zm="zm1.0";

void PlotsEnergyPedestal(){
 
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

  TString grid="grid24";
  TString conf="conf1";

  int islabs=6;
  TString nslabs_pedestal=TString::Format("pedestal_nslabs%i",islabs);
  TString nslabs=TString::Format("nslabs%i",islabs);

  TString bcid="bcidmax2850";
      
  for(int icases=1; icases<3; icases++) {
    
    
    if(icases==1) {
      grid="grid20";
      conf="conf1";
    }
    if(icases==2) conf="conf2";
    if(icases==3) conf="conf3";
    
    
    TString s_file=nslabs+"_"+bcid+"/"+conf+"_"+grid+"_1GeV_mipcut0.5_showers.root";
    std::cout<<"Opening file: "<<s_file<<std::endl;
    TFile *file = new TFile(s_file);
    TH1F *energy_center_zm1_1GeV = (TH1F*)file->Get("energy_center_zm05");
    energy_center_zm1_1GeV->Rebin(4);
    
    s_file=nslabs+"_"+bcid+"/"+conf+"_"+grid+"_2GeV_mipcut0.5_showers.root";
    std::cout<<"Opening file: "<<s_file<<std::endl;
    file = new TFile(s_file);
    TH1F *energy_center_zm1_2GeV = (TH1F*)file->Get("energy_center_zm05");
    energy_center_zm1_2GeV->Rebin(4);

    s_file=nslabs+"_"+bcid+"/"+conf+"_"+grid+"_4GeV_mipcut0.5_showers.root";
    std::cout<<"Opening file: "<<s_file<<std::endl;
    file = new TFile(s_file);
    TH1F *energy_center_zm1_4GeV = (TH1F*)file->Get("energy_center_zm05");
    energy_center_zm1_4GeV->Rebin(4);

    s_file=nslabs+"_"+bcid+"/"+conf+"_"+grid+"_5GeV_mipcut0.5_showers.root";
    std::cout<<"Opening file: "<<s_file<<std::endl;
    file = new TFile(s_file);
    TH1F *energy_center_zm1_5GeV = (TH1F*)file->Get("energy_center_zm05");
    energy_center_zm1_5GeV->Rebin(4);

    s_file=nslabs+"_"+bcid+"/"+conf+"_"+grid+"_5.8GeV_mipcut0.5_showers.root";
    std::cout<<"Opening file: "<<s_file<<std::endl;
    file = new TFile(s_file);
    TH1F *energy_center_zm1_58GeV = (TH1F*)file->Get("energy_center_zm05");
    energy_center_zm1_58GeV->Rebin(4);

    // -----------------------

    s_file=nslabs_pedestal+"_"+bcid+"/"+conf+"_"+grid+"_1GeV_mipcut0.5_showers.root";
    std::cout<<"Opening file: "<<s_file<<std::endl;
    file = new TFile(s_file);
    TH1F *energy_center_zm1_pedestal_1GeV = (TH1F*)file->Get("energy_center_zm05");
    energy_center_zm1_pedestal_1GeV->Rebin(4);
    
    s_file=nslabs_pedestal+"_"+bcid+"/"+conf+"_"+grid+"_2GeV_mipcut0.5_showers.root";
    std::cout<<"Opening file: "<<s_file<<std::endl;
    file = new TFile(s_file);
    TH1F *energy_center_zm1_pedestal_2GeV = (TH1F*)file->Get("energy_center_zm05");
    energy_center_zm1_pedestal_2GeV->Rebin(4);

    s_file=nslabs_pedestal+"_"+bcid+"/"+conf+"_"+grid+"_4GeV_mipcut0.5_showers.root";
    std::cout<<"Opening file: "<<s_file<<std::endl;
    file = new TFile(s_file);
    TH1F *energy_center_zm1_pedestal_4GeV = (TH1F*)file->Get("energy_center_zm05");
    energy_center_zm1_pedestal_4GeV->Rebin(4);

    s_file=nslabs_pedestal+"_"+bcid+"/"+conf+"_"+grid+"_5GeV_mipcut0.5_showers.root";
    std::cout<<"Opening file: "<<s_file<<std::endl;
    file = new TFile(s_file);
    TH1F *energy_center_zm1_pedestal_5GeV = (TH1F*)file->Get("energy_center_zm05");
    energy_center_zm1_pedestal_5GeV->Rebin(4);

    s_file=nslabs_pedestal+"_"+bcid+"/"+conf+"_"+grid+"_5.8GeV_mipcut0.5_showers.root";
    std::cout<<"Opening file: "<<s_file<<std::endl;
    file = new TFile(s_file);
    TH1F *energy_center_zm1_pedestal_58GeV = (TH1F*)file->Get("energy_center_zm05");
    energy_center_zm1_pedestal_58GeV->Rebin(4);
 
    TCanvas *c_energy = new TCanvas("c_energy","c_energy",800,600);
    //c_energy->Divide(2,1);
    c_energy->cd(1);
    gPad->SetLogy();
    energy_center_zm1_1GeV->SetTitle("SiW-ECAL (wafer-3), e^{+} beam");
    energy_center_zm1_1GeV->GetXaxis()->SetTitle("E^{raw}/MIP = #sum_{i=layers} #sum_{j=cells} #omega_{i} E_{i,j}^{raw}");
    energy_center_zm1_1GeV->GetYaxis()->SetTitle("# entries");
    energy_center_zm1_1GeV->GetYaxis()->SetRangeUser(5,5000);
    energy_center_zm1_1GeV->SetLineWidth(1);
    energy_center_zm1_1GeV->SetLineColor(1);
    energy_center_zm1_1GeV->Draw();
    
    energy_center_zm1_2GeV->SetLineColor(2);
    energy_center_zm1_2GeV->SetLineWidth(1);
    energy_center_zm1_2GeV->Draw("same");
    
    energy_center_zm1_4GeV->SetLineColor(3);
    energy_center_zm1_4GeV->SetLineWidth(1);
    energy_center_zm1_4GeV->Draw("same");
    
    energy_center_zm1_5GeV->SetLineColor(4);
    energy_center_zm1_5GeV->SetLineWidth(1);
    energy_center_zm1_5GeV->Draw("same");
    
    energy_center_zm1_58GeV->SetLineColor(6);
    energy_center_zm1_58GeV->SetLineWidth(1);
    energy_center_zm1_58GeV->Draw("same");

    
    energy_center_zm1_pedestal_1GeV->SetLineStyle(2);
    energy_center_zm1_pedestal_1GeV->SetLineColor(1);
    energy_center_zm1_pedestal_1GeV->Draw("same");
    
    energy_center_zm1_pedestal_2GeV->SetLineColor(2);
    energy_center_zm1_pedestal_2GeV->SetLineStyle(2);
    energy_center_zm1_pedestal_2GeV->Draw("same");
    
    energy_center_zm1_pedestal_4GeV->SetLineColor(3);
    energy_center_zm1_pedestal_4GeV->SetLineStyle(2);
    energy_center_zm1_pedestal_4GeV->Draw("same");
    
    energy_center_zm1_pedestal_5GeV->SetLineColor(4);
    energy_center_zm1_pedestal_5GeV->SetLineStyle(2);
    energy_center_zm1_pedestal_5GeV->Draw("same");
    
    energy_center_zm1_pedestal_58GeV->SetLineColor(6);
    energy_center_zm1_pedestal_58GeV->SetLineStyle(2);
    energy_center_zm1_pedestal_58GeV->Draw("same");
    
    TLegend *l_1 = new TLegend(0.45,0.65,0.9,0.85);
    l_1->SetHeader("SiW-ECAL: wafer 4, W-configuration 1");
    if(icases==1)   l_1->SetHeader("SiW-ECAL: wafer 3, W-configuration 1");
    if(icases==2)   l_1->SetHeader("SiW-ECAL: wafer 3, W-configuration 2");
    if(icases==3)   l_1->SetHeader("SiW-ECAL: wafer 3, W-configuration 3");
    
    l_1->AddEntry(energy_center_zm1_1GeV,"e^{+} 1 GeV, precalculated pedestals","l");
    l_1->AddEntry(energy_center_zm1_pedestal_1GeV,"e^{+} 1 GeV, pedestals calculated on the fly","l");
    l_1->AddEntry(energy_center_zm1_2GeV,"e^{+} 2 GeV, precalculated pedestals","l");
    l_1->AddEntry(energy_center_zm1_4GeV,"e^{+} 4 GeV, precalculated pedestals","l");
    l_1->AddEntry(energy_center_zm1_5GeV,"e^{+} 5 GeV, precalculated pedestals","l");
    l_1->AddEntry(energy_center_zm1_58GeV,"e^{+} 5.8 GeV, precalculated pedestals","l");

    // l_1->AddEntry(energy_center_zm1,TString::Format("e^{+} 2 GeV, contained showers:  gaussian fit #chi^{2}/ndf=%0.1f",chi3),"l");
    // l_1->AddEntry(energy_2,TString::Format("e^{+} 4 GeV: gaussian fit #chi^{2}/ndf=%0.1f",chi2),"l");
    // l_1->AddEntry(energy_center_zm1_2,TString::Format("e^{+} 4 GeV, contained showers: gaussian fit #chi^{2}/ndf=%0.1f",chi4),"l");
    l_1->SetFillColor(0);
    l_1->SetLineColor(0);
    l_1->SetShadowColor(0);
    l_1->Draw();
    IRLESLabel(0.2,0.88,"",kGray+2);

    c_energy->Print("pedestal_plots/energy_"+nslabs+"_"+bcid+"_"+conf+"_"+grid+"_"+"_mipcut0.5.eps");


    TCanvas *c_energy_2 = new TCanvas("c_energy_2","c_energy_2",800,600);
    //c_energy->Divide(2,1);
    c_energy_2->cd(1);
    // gPad->SetLogy();

    energy_center_zm1_1GeV->Rebin(5);
    energy_center_zm1_2GeV->Rebin(5);
    energy_center_zm1_4GeV->Rebin(5);
    energy_center_zm1_5GeV->Rebin(5);
    energy_center_zm1_58GeV->Rebin(5);

    energy_center_zm1_pedestal_1GeV->Rebin(5);
    energy_center_zm1_pedestal_2GeV->Rebin(5);
    energy_center_zm1_pedestal_4GeV->Rebin(5);
    energy_center_zm1_pedestal_5GeV->Rebin(5);
    energy_center_zm1_pedestal_58GeV->Rebin(5);
    
    energy_center_zm1_1GeV->Add(energy_center_zm1_pedestal_1GeV,-1);
    energy_center_zm1_2GeV->Add(energy_center_zm1_pedestal_2GeV,-1);
    energy_center_zm1_4GeV->Add(energy_center_zm1_pedestal_4GeV,-1);
    energy_center_zm1_5GeV->Add(energy_center_zm1_pedestal_5GeV,-1);
    energy_center_zm1_58GeV->Add(energy_center_zm1_pedestal_58GeV,-1);

    energy_center_zm1_1GeV->Divide(energy_center_zm1_pedestal_1GeV);
    energy_center_zm1_2GeV->Divide(energy_center_zm1_pedestal_2GeV);
    energy_center_zm1_4GeV->Divide(energy_center_zm1_pedestal_4GeV);
    energy_center_zm1_5GeV->Divide(energy_center_zm1_pedestal_5GeV);
    energy_center_zm1_58GeV->Divide(energy_center_zm1_pedestal_58GeV);

    energy_center_zm1_1GeV->SetTitle("SiW-ECAL (wafer-3), e^{+} beam");
    energy_center_zm1_1GeV->GetXaxis()->SetTitle("E^{raw}/MIP = #sum_{i=layers} #sum_{j=cells} #omega_{i} E_{i,j}^{raw}");
    energy_center_zm1_1GeV->GetYaxis()->SetTitle("deviation [%]");
    energy_center_zm1_1GeV->GetYaxis()->SetRangeUser(-5,5);
    energy_center_zm1_1GeV->SetLineWidth(1);
    energy_center_zm1_1GeV->SetLineColor(1);
    energy_center_zm1_1GeV->Draw();
    
    energy_center_zm1_2GeV->SetLineColor(2);
    energy_center_zm1_2GeV->SetLineWidth(1);
    energy_center_zm1_2GeV->Draw("same");
    
    energy_center_zm1_4GeV->SetLineColor(3);
    energy_center_zm1_4GeV->SetLineWidth(1);
    energy_center_zm1_4GeV->Draw("same");
    
    energy_center_zm1_5GeV->SetLineColor(4);
    energy_center_zm1_5GeV->SetLineWidth(1);
    energy_center_zm1_5GeV->Draw("same");
    
    energy_center_zm1_58GeV->SetLineColor(6);
    energy_center_zm1_58GeV->SetLineWidth(1);
    energy_center_zm1_58GeV->Draw("same");
    
    TLegend *l_2 = new TLegend(0.45,0.65,0.9,0.85);
    l_2->SetHeader("SiW-ECAL: wafer 4, W-configuration 1");
    if(icases==1)   l_2->SetHeader("SiW-ECAL: wafer 3, W-configuration 1");
    if(icases==2)   l_2->SetHeader("SiW-ECAL: wafer 3, W-configuration 2");
    if(icases==3)   l_2->SetHeader("SiW-ECAL: wafer 3, W-configuration 3");
    
    l_2->AddEntry(energy_center_zm1_1GeV,"e^{+} 1 GeV","l");
    l_2->AddEntry(energy_center_zm1_2GeV,"e^{+} 2 GeV","l");
    l_2->AddEntry(energy_center_zm1_4GeV,"e^{+} 4 GeV","l");
    l_2->AddEntry(energy_center_zm1_5GeV,"e^{+} 5 GeV","l");
    l_2->AddEntry(energy_center_zm1_58GeV,"e^{+} 5.8 GeV","l");

    // l_2->AddEntry(energy_center_zm1,TString::Format("e^{+} 2 GeV, contained showers:  gaussian fit #chi^{2}/ndf=%0.1f",chi3),"l");
    // l_2->AddEntry(energy_2,TString::Format("e^{+} 4 GeV: gaussian fit #chi^{2}/ndf=%0.1f",chi2),"l");
    // l_2->AddEntry(energy_center_zm1_2,TString::Format("e^{+} 4 GeV, contained showers: gaussian fit #chi^{2}/ndf=%0.1f",chi4),"l");
    l_2->SetFillColor(0);
    l_2->SetLineColor(0);
    l_2->SetShadowColor(0);
    l_2->Draw();

    IRLESLabel(0.2,0.88,"",kGray+2);

    
    c_energy_2->Print("pedestal_plots/ratio_energy_"+nslabs+"_"+bcid+"_"+conf+"_"+grid+"_"+"_mipcut0.5.eps");
    
  }
         
}
