#include "TROOT.h"
#include "TFile.h"
#include "../../style/Style.C"
#include "TGraphErrors.h"

void PlotsShower(){


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

  for(int icases=2; icases<4; icases++) {

    if(icases>0) grid="grid20";
    if(icases==2) conf="conf2";
    if(icases==3) conf="conf3";

    //energy_1 GeV
    TString s_file="zmlt1/"+conf+"_"+grid+"_1GeV_mipcut0.5_showers.root";
    std::cout<<"Opening file: "<<s_file<<std::endl;
    TFile *file_1 = new TFile(s_file);
    TH1F *energy_profile_1 = (TH1F*)file_1->Get("energy_profile_z");

    //energy_2 GeV
    s_file="zmlt1/"+conf+"_"+grid+"_2GeV_mipcut0.5_showers.root";
    std::cout<<"Opening file: "<<s_file<<std::endl;
    TFile *file_2 = new TFile(s_file);
    TH1F *energy_profile_2 = (TH1F*)file_2->Get("energy_profile_z");
    
    //energy_3 GeV
    s_file="zmlt1/"+conf+"_"+grid+"_3GeV_mipcut0.5_showers.root";
    std::cout<<"Opening file: "<<s_file<<std::endl;
    TFile *file_3 = new TFile(s_file);
    TH1F *energy_profile_3 = (TH1F*)file_3->Get("energy_profile_z");

    //energy_4 GeV
    s_file="zmlt1/"+conf+"_"+grid+"_4GeV_mipcut0.5_showers.root";
    std::cout<<"Opening file: "<<s_file<<std::endl;
    TFile *file_4 = new TFile(s_file);
    TH1F *energy_profile_4 = (TH1F*)file_4->Get("energy_profile_z");

    //energy_5 GeV
    s_file="zmlt1/"+conf+"_"+grid+"_5GeV_mipcut0.5_showers.root";
    std::cout<<"Opening file: "<<s_file<<std::endl;
    TFile *file_5 = new TFile(s_file);
    TH1F *energy_profile_5 = (TH1F*)file_5->Get("energy_profile_z");

    //energy_5.8 GeV
    s_file="zmlt1/"+conf+"_"+grid+"_5.8GeV_mipcut0.5_showers.root";
    std::cout<<"Opening file: "<<s_file<<std::endl;
    TFile *file_6 = new TFile(s_file);
    TH1F *energy_profile_6 = (TH1F*)file_6->Get("energy_profile_z");

    //  file->Close();

    TGraphErrors *g_energy_profile_1 = new TGraphErrors(energy_profile_1);
    Double_t *x_1 = g_energy_profile_1->GetX();
    Double_t *ex_1 = g_energy_profile_1->GetEX();
    for(int i=0; i<g_energy_profile_1->GetN(); i++) {
      x_1[i]+=1.;
      ex_1[i]=0;
    }

         TGraphErrors *g_energy_profile_2 = new TGraphErrors(energy_profile_2);
    Double_t *x_2 = g_energy_profile_2->GetX();
    Double_t *ex_2 = g_energy_profile_2->GetEX();
    for(int i=0; i<g_energy_profile_2->GetN(); i++) {
      x_2[i]+=1.;
      ex_2[i]=0;
    }

       TGraphErrors *g_energy_profile_3 = new TGraphErrors(energy_profile_3);
    Double_t *x_3 = g_energy_profile_3->GetX();
    Double_t *ex_3 = g_energy_profile_3->GetEX();
    for(int i=0; i<g_energy_profile_3->GetN(); i++) {
      x_3[i]+=1.;
      ex_3[i]=0;
    }

        TGraphErrors *g_energy_profile_4 = new TGraphErrors(energy_profile_4);
    Double_t *x_4 = g_energy_profile_4->GetX();
    Double_t *ex_4 = g_energy_profile_4->GetEX();
    for(int i=0; i<g_energy_profile_4->GetN(); i++) {
      x_4[i]+=1.;
      ex_4[i]=0;
    }

           TGraphErrors *g_energy_profile_5 = new TGraphErrors(energy_profile_5);
    Double_t *x_5 = g_energy_profile_5->GetX();
    Double_t *ex_5 = g_energy_profile_5->GetEX();
    for(int i=0; i<g_energy_profile_5->GetN(); i++) {
      x_5[i]+=1.;
      ex_5[i]=0;
    }

        TGraphErrors *g_energy_profile_6 = new TGraphErrors(energy_profile_6);
    Double_t *x_6 = g_energy_profile_6->GetX();
    Double_t *ex_6 = g_energy_profile_6->GetEX();
    for(int i=0; i<g_energy_profile_6->GetN(); i++) {
      x_6[i]+=1.;
      ex_6[i]=0;
    }
    
    TCanvas *c_energy = new TCanvas("c_energy","c_energy",800,600);
    //c_energy->Divide(2,1);
    c_energy->cd(1);
    // energy->SetTitle("SiW-ECAL (wafer-3), e^{+} beam");

    
    g_energy_profile_1->GetYaxis()->SetRangeUser(0,0.5);
    g_energy_profile_1->GetXaxis()->SetTitle("layer number");
    g_energy_profile_1->GetYaxis()->SetTitle("average contribution to  #tilde{E}^{raw}");
    g_energy_profile_1->Draw("alp");

    g_energy_profile_2->SetLineColor(2);
    g_energy_profile_2->SetMarkerColor(2);
    g_energy_profile_2->SetMarkerStyle(21);
    g_energy_profile_2->Draw("lp");

    g_energy_profile_3->SetLineColor(4);
    g_energy_profile_3->SetMarkerColor(4);
    g_energy_profile_3->SetMarkerStyle(22);
    g_energy_profile_3->Draw("lp");
  
    g_energy_profile_4->SetLineColor(6);
    g_energy_profile_4->SetMarkerColor(6);
    g_energy_profile_4->SetMarkerStyle(24);
    g_energy_profile_4->Draw("lp");

    g_energy_profile_5->SetLineColor(8);
    g_energy_profile_5->SetMarkerColor(8);
    g_energy_profile_5->SetMarkerStyle(25);
    g_energy_profile_5->Draw("lp");

    g_energy_profile_6->SetLineColor(5);
    g_energy_profile_6->SetMarkerColor(5);
    g_energy_profile_6->SetMarkerStyle(26);
    g_energy_profile_6->Draw("lp");
  

    TLegend *l_energy = new TLegend(0.5,0.65,0.93,0.9);
    l_energy->SetHeader("SiW-ECAL: wafer 4, W-configuration 1");
    if(icases==1)   l_energy->SetHeader("SiW-ECAL: wafer 3, W-configuration 1");
    if(icases==2)   l_energy->SetHeader("SiW-ECAL: wafer 3, W-configuration 2");
    if(icases==3)   l_energy->SetHeader("SiW-ECAL: wafer 3, W-configuration 3");

    l_energy->AddEntry(g_energy_profile_1,"e^{+} 1 GeV beam","lp");
    l_energy->AddEntry(g_energy_profile_2,"e^{+} 2 GeV beam","lp");
    l_energy->AddEntry(g_energy_profile_3,"e^{+} 3 GeV beam","lp");
    l_energy->AddEntry(g_energy_profile_4,"e^{+} 4 GeV beam","lp");
    l_energy->AddEntry(g_energy_profile_5,"e^{+} 5 GeV beam","lp");
    l_energy->AddEntry(g_energy_profile_6,"e^{+} 5.8 GeV beam","lp");

    l_energy->SetFillColor(0);
    l_energy->SetLineColor(0);
    l_energy->SetShadowColor(0);
    l_energy->Draw();
  
    c_energy->Print("zmlt1/energy_z_profile_"+grid+"_"+conf+".eps");

  }


  
  grid="grid24";
  conf="conf1";
  
  TGraphErrors *g_energy_profile_1;
  TGraphErrors *g_energy_profile_2;
  TGraphErrors *g_energy_profile_3;
  TGraphErrors *g_energy_profile_4;
  

  for(int icases=0; icases<4; icases++) {

    if(icases>0) grid="grid20";
    if(icases==2) conf="conf2";
    if(icases==3) conf="conf3";

    //energy_4 GeV
    TString s_file="zmlt1/"+conf+"_"+grid+"_4GeV_mipcut0.5_showers.root";
    std::cout<<"Opening file: "<<s_file<<std::endl;
    TFile *file = new TFile(s_file);
    TH1F *energy_profile = (TH1F*)file->Get("energy_profile_z");

    if(icases==0) g_energy_profile_1 = new TGraphErrors(energy_profile);
    if(icases==1) g_energy_profile_2 = new TGraphErrors(energy_profile);
    if(icases==2) g_energy_profile_3 = new TGraphErrors(energy_profile);
    if(icases==3) g_energy_profile_4 = new TGraphErrors(energy_profile);

  }


    Double_t *x_1 = g_energy_profile_1->GetX();
    Double_t *ex_1 = g_energy_profile_1->GetEX();
    for(int i=0; i<g_energy_profile_1->GetN(); i++) {
      x_1[i]+=1.;
      ex_1[i]=0;
    }

    Double_t *x_2 = g_energy_profile_2->GetX();
    Double_t *ex_2 = g_energy_profile_2->GetEX();
    for(int i=0; i<g_energy_profile_2->GetN(); i++) {
      x_2[i]+=1.;
      ex_2[i]=0;
    }

    Double_t *x_3 = g_energy_profile_3->GetX();
    Double_t *ex_3 = g_energy_profile_3->GetEX();
    for(int i=0; i<g_energy_profile_3->GetN(); i++) {
      x_3[i]+=1.;
      ex_3[i]=0;
    }

    Double_t *x_4 = g_energy_profile_4->GetX();
    Double_t *ex_4 = g_energy_profile_4->GetEX();
    for(int i=0; i<g_energy_profile_4->GetN(); i++) {
      x_4[i]+=1.;
      ex_4[i]=0;
    }
    

  
    TCanvas *c_energy2 = new TCanvas("c_energy2","c_energy2",800,600);
    //c_energy->Divide(2,1);
    c_energy2->cd(1);
    // energy->SetTitle("SiW-ECAL (wafer-3), e^{+} beam");

    g_energy_profile_4->SetLineColor(8);
    g_energy_profile_4->SetMarkerColor(8);
    g_energy_profile_4->SetMarkerStyle(24);
    g_energy_profile_4->GetYaxis()->SetRangeUser(0,0.5);
    g_energy_profile_4->GetXaxis()->SetTitle("layer number");
    g_energy_profile_4->GetYaxis()->SetTitle("average contribution to  #tilde{E}^{raw}");

    g_energy_profile_4->Draw("alp");
    
    g_energy_profile_1->Draw("lp");

    g_energy_profile_2->SetLineColor(1);
    g_energy_profile_2->SetLineStyle(2);
    g_energy_profile_2->SetMarkerColor(1);
    g_energy_profile_2->SetMarkerStyle(4);
    g_energy_profile_2->Draw("lp");

    g_energy_profile_3->SetLineColor(4);
    g_energy_profile_3->SetMarkerColor(4);
    g_energy_profile_3->SetMarkerStyle(22);
    g_energy_profile_3->Draw("lp");
  
    g_energy_profile_4->SetLineColor(8);
    g_energy_profile_4->SetMarkerColor(8);
    g_energy_profile_4->SetMarkerStyle(24);
    g_energy_profile_4->Draw("lp");


    TLegend *l_energy = new TLegend(0.5,0.7,0.8,0.9);
    l_energy->SetHeader("SiW-ECAL, 4 GeV positrons");
 
    l_energy->AddEntry(g_energy_profile_1,"Wafer 3, W-configuration 1","lp");
    l_energy->AddEntry(g_energy_profile_2,"Wafer 4, W-configuration 1","lp");
    l_energy->AddEntry(g_energy_profile_3,"Wafer 4, W-configuration 2","lp");
    l_energy->AddEntry(g_energy_profile_4,"Wafer 4, W-configuration 3","lp");
    l_energy->SetFillColor(0);
    l_energy->SetLineColor(0);
    l_energy->SetShadowColor(0);
    l_energy->Draw();
  
    c_energy2->Print("zmlt1/4GeV_energy_z_profile.eps");

  


  
}
