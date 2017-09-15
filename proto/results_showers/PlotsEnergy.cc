#include "TROOT.h"
#include "TFile.h"
#include "Style.C"
#include "TGraphErrors.h"
#include "TMath.h"

void PlotsEnergy(){
 


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

  
  for(int icases=0; icases<5; icases++) {

    if(icases==1) grid="grid20";
    if(icases==2) conf="conf2";
    if(icases==3) conf="conf3";

  
    TString s_file="zmlt1_5slabs/"+conf+"_"+grid+"_2GeV_mipcut0.5_showers.root";
    std::cout<<"Opening file: "<<s_file<<std::endl;
    TFile *file = new TFile(s_file);
  
    TH1F *energy = (TH1F*)file->Get("energy");
    TH1F *energy_center = (TH1F*)file->Get("energy_center");

    energy->Rebin(4);
    energy_center->Rebin(4);

    //  file->Close();

    TString s_file_2="zmlt1_5slabs/"+conf+"_"+grid+"_4GeV_mipcut0.5_showers.root";

    std::cout<<"Opening file: "<<s_file_2<<std::endl;
    TFile *file_2 = new TFile(s_file_2);
  
    TH1F *energy_2 = (TH1F*)file_2->Get("energy");
    TH1F *energy_center_2 = (TH1F*)file_2->Get("energy_center");

    energy_2->Rebin(4);
    energy_center_2->Rebin(4);

    TCanvas *c_energy = new TCanvas("c_energy","c_energy",800,600);
    //c_energy->Divide(2,1);
    c_energy->cd(1);
    // energy->SetTitle("SiW-ECAL (wafer-3), e^{+} beam");
    energy->GetXaxis()->SetTitle("E^{raw}/MIP = #sum_{i=layers} #sum_{j=cells} #omega_{i} E_{i,j}^{raw}");
    energy->GetYaxis()->SetTitle("# entries");
    energy->Draw();
  
    energy_center->SetLineColor(4);
    energy_center->Draw("same");

    energy_2->SetLineStyle(2);
    energy_2->Draw("same");

    energy_center_2->SetLineColor(4);
    energy_center_2->SetLineStyle(2);
    energy_center_2->Draw("same");

    double xmin=energy->GetMean()-2*energy->GetRMS();
    double xmax=energy->GetMean()+2*energy->GetRMS();
    TF1 *fit_1 = new TF1("fit_1","gaus",xmin,xmax);
    fit_1->SetLineColor(1);
    fit_1->SetLineWidth(2);
    energy->Fit(fit_1,"MRQ");

    xmin=energy_2->GetMean()-3*energy_2->GetRMS();
    xmax=energy_2->GetMean()+3*energy_2->GetRMS();
    TF1 *fit_2 = new TF1("fit_2","gaus",xmin,xmax);
    energy_2->Fit(fit_2,"MRQN");
    xmin=fit_2->GetParameter(1)-2*fit_2->GetParameter(2);
    xmax=fit_2->GetParameter(1)+2*fit_2->GetParameter(2);
    TF1 *fit_3 = new TF1("fit_3","gaus",xmin,xmax);
    fit_3->SetLineColor(1);
    fit_3->SetLineStyle(2);
    fit_3->SetLineWidth(2);
    energy_2->Fit(fit_3,"MRQ");

    //--------

    xmin=energy_center->GetMean()-2*energy_center->GetRMS();
    xmax=energy_center->GetMean()+2*energy_center->GetRMS();
    TF1 *fit_center1 = new TF1("fit_center1","gaus",xmin,xmax);
    fit_center1->SetLineColor(4);
    fit_center1->SetLineWidth(2);
    energy_center->Fit(fit_center1,"MRQ");

    xmin=energy_center_2->GetMean()-3*energy_center_2->GetRMS();
    xmax=energy_center_2->GetMean()+3*energy_center_2->GetRMS();
    TF1 *fit_center2 = new TF1("fit_center2","gaus",xmin,xmax);
    energy_center_2->Fit(fit_center2,"MRQN");
    xmin=fit_center2->GetParameter(1)-2*fit_center2->GetParameter(2);
    xmax=fit_center2->GetParameter(1)+2*fit_center2->GetParameter(2);
    TF1 *fit_center3 = new TF1("fit_center3","gaus",xmin,xmax);
    fit_center3->SetLineColor(4);
    fit_center3->SetLineStyle(2);
    fit_center3->SetLineWidth(2);
    energy_center_2->Fit(fit_center3,"MRQ");

    TLegend *l_energy = new TLegend(0.4,0.7,0.92,0.9);
    l_energy->SetHeader("SiW-ECAL: wafer 4, W-configuration 1");
    if(icases==1)   l_energy->SetHeader("SiW-ECAL: wafer 3, W-configuration 1");
    if(icases==2)   l_energy->SetHeader("SiW-ECAL: wafer 3, W-configuration 2");
    if(icases==3)   l_energy->SetHeader("SiW-ECAL: wafer 3, W-configuration 3");


    double chi1=fit_1->GetChisquare()/fit_1->GetNDF();
    double chi2=fit_3->GetChisquare()/fit_3->GetNDF();
    double chi3=fit_center1->GetChisquare()/fit_center1->GetNDF();
    double chi4=fit_center3->GetChisquare()/fit_center3->GetNDF();
    
    
    l_energy->AddEntry(energy,TString::Format("e^{+} 2 GeV: gaussian fit #chi^{2}/ndf=%0.1f",chi1),"l");
    l_energy->AddEntry(energy_center,TString::Format("e^{+} 2 GeV, contained showers:  gaussian fit #chi^{2}/ndf=%0.1f",chi3),"l");
    l_energy->AddEntry(energy_2,TString::Format("e^{+} 4 GeV: gaussian fit #chi^{2}/ndf=%0.1f",chi2),"l");
    l_energy->AddEntry(energy_center_2,TString::Format("e^{+} 4 GeV, contained showers: gaussian fit #chi^{2}/ndf=%0.1f",chi4),"l");
    l_energy->SetFillColor(0);
    l_energy->SetLineColor(0);
    l_energy->SetShadowColor(0);
    l_energy->Draw();
  
    c_energy->Print("zmlt1_5slabs/energy_"+grid+"_"+conf+".eps");

  }

  
  /// ####################################################33
  // ENERGY LINEARITY

  grid="grid24";
  conf="conf1";
  
  TGraphErrors * g_grid24_conf1;
  TGraphErrors * g_grid20_conf1;
  TGraphErrors * g_grid20_conf2;
  TGraphErrors * g_grid20_conf3;

  
  for(int icases=0; icases<4; icases++) {

    if(icases==1) grid="grid20";
    if(icases==2) conf="conf2";
    if(icases==3) conf="conf3";
  
    int np=0;
    double x[7], y[7], ex[7],ey[7];
    
    
    for(int ienergy=0; ienergy<6; ienergy++) {
      
      x[ienergy]=0.;
      y[ienergy]=0.;
      ex[ienergy]=0.;
      ey[ienergy]=0.;
      
      if(ienergy==4 && grid=="grid24") continue;
      if(ienergy==2 && grid=="grid20" && conf=="conf1") continue;
      
      
      double ie=0;
      if(ienergy == 0) ie=1;
      if(ienergy == 1) ie=2;
      if(ienergy == 2) ie=3;
      if(ienergy == 3) ie=4;
      if(ienergy == 4) ie=5;
      if(ienergy == 5) ie=5.8;
      
      
      TString s_energy;
      if(ienergy<5) s_energy=TString::Format("%iGeV",int(ie));
      else s_energy = "5.8GeV";
      
      x[np]=ie;
      
      
      TString s_file="zmlt1_5slabs/"+conf+"_"+grid+"_"+s_energy+"_mipcut0.5_showers.root";
      std::cout<<"Opening file: "<<s_file<<std::endl;
      TFile *file = new TFile(s_file);
      
      TH1F *energy_center = (TH1F*)file->Get("energy_center");
      energy_center->Rebin(4);
      
      double xmin=energy_center->GetMean()-3*energy_center->GetRMS();
      double xmax=energy_center->GetMean()+3*energy_center->GetRMS();
      
      TF1 *fit_1 = new TF1("fit_1","gaus",xmin,xmax);
      energy_center->Fit(fit_1,"MRQ");
      
      xmin=fit_1->GetParameter(1)-3*fit_1->GetParameter(2);
      xmax=fit_1->GetParameter(1)+3*fit_1->GetParameter(2);
      
      TF1 *fit_2 = new TF1("fit_2","gaus",xmin,xmax);
      energy_center->Fit(fit_2,"MRQ");
      
      y[np]=fit_2->GetParameter(1);
      ey[np]=fit_2->GetParError(1);

      TAxis *xaxis = energy_center->GetXaxis();
      TAxis *yaxis = energy_center->GetYaxis();
      Int_t binxmin = xaxis->FindBin(xmin);
      Int_t binxmax = xaxis->FindBin(xmax);

      //ey[np]=TMath::Min(fit_2->GetParameter(2)/sqrt(energy_center->Integral(binxmin,binxmax)),ey[np]);
      
      np++;
      
      
    }

    if(icases==0) g_grid24_conf1 = new TGraphErrors(np,x,y,ex,ey);
    if(icases==1) g_grid20_conf1 = new TGraphErrors(np,x,y,ex,ey);
    if(icases==2) g_grid20_conf2 = new TGraphErrors(np,x,y,ex,ey);
    if(icases==3) g_grid20_conf3 = new TGraphErrors(np,x,y,ex,ey);

  }

 
    TCanvas *c_energy_linearity = new TCanvas("c_energy_linearity","c_energy_linearity",800,600);
    //c_energy_linearity->Divide(2,1);
    c_energy_linearity->cd(1);

    g_grid24_conf1->SetLineWidth(2);
    g_grid24_conf1->SetLineStyle(2);
    g_grid24_conf1->SetMarkerStyle(4);
    g_grid24_conf1->GetYaxis()->SetTitle("E^{raw}/MIP");
    g_grid24_conf1->GetXaxis()->SetTitle("E^{beam}/GeV");
    g_grid24_conf1->GetYaxis()->SetRangeUser(0,900);
    g_grid24_conf1->Draw("alp");

    g_grid20_conf1->SetLineWidth(2);
    g_grid20_conf1->SetLineStyle(1);
    g_grid20_conf1->SetMarkerStyle(20);
    g_grid20_conf1->Draw("lp");


    g_grid20_conf2->SetLineWidth(2);
    g_grid20_conf2->SetLineColor(8);
    g_grid20_conf2->SetMarkerStyle(22);
    g_grid20_conf2->SetMarkerColor(8);
    g_grid20_conf2->Draw("lp");

    
    g_grid20_conf3->SetLineWidth(2);
    g_grid20_conf3->SetLineColor(4);
    g_grid20_conf3->SetMarkerStyle(23);
    g_grid20_conf3->SetMarkerColor(4);
    g_grid20_conf3->Draw("lp");
    

    TLegend *l_energy_linearity = new TLegend(0.2,0.7,0.7,0.9);
    l_energy_linearity->SetHeader("SiW-ECAL: reconstructed energy linearity");

    l_energy_linearity->AddEntry(g_grid24_conf1,"Wafer 3, W-configuration 1","lp");
    l_energy_linearity->AddEntry(g_grid20_conf1,"Wafer 4, W-configuration 1","lp");
    l_energy_linearity->AddEntry(g_grid20_conf2,"Wafer 4, W-configuration 2","lp");
    l_energy_linearity->AddEntry(g_grid20_conf3,"Wafer 4, W-configuration 3","lp");
    l_energy_linearity->SetFillColor(0);
    l_energy_linearity->SetLineColor(0);
    l_energy_linearity->SetShadowColor(0);
    l_energy_linearity->Draw();
  
    c_energy_linearity->Print("zmlt1_5slabs/energy_linearity_zcut_1_mipcut0.5.eps");
    c_energy_linearity->Print("zmlt1_5slabs/energy_linearity_zcut_1_mipcut0.5.C");


    //###  #######################################################################
    // -------------------------- Fit


    //24/1 ----------------------------------------------------
    TF1 *fit_pol1_1 = new TF1("fit_pol1_1","pol1",0.5,6);
    g_grid24_conf1->Fit(fit_pol1_1,"RMN");

    Double_t *x_fit =  g_grid24_conf1->GetX();
    Double_t *y_fit =  g_grid24_conf1->GetY();
    Double_t *ex_fit = g_grid24_conf1->GetEX();
    Double_t *ey_fit = g_grid24_conf1->GetEY();

    int n_fit =  g_grid24_conf1->GetN();

    double xfit[10], yfit[10], exfit[10], eyfit[10];

    for(int i=0; i<g_grid24_conf1->GetN(); i++) {
      double e1=  fit_pol1_1->GetParameter(0)+fit_pol1_1->GetParError(0)+ (fit_pol1_1->GetParameter(1) - fit_pol1_1->GetParError(1)) * x_fit[i];
      double e2=  fit_pol1_1->GetParameter(0)+fit_pol1_1->GetParError(0)+ (fit_pol1_1->GetParameter(1) + fit_pol1_1->GetParError(1)) * x_fit[i];
      double e3=  fit_pol1_1->GetParameter(0)-fit_pol1_1->GetParError(0)+ (fit_pol1_1->GetParameter(1) - fit_pol1_1->GetParError(1)) * x_fit[i];
      double e4=  fit_pol1_1->GetParameter(0)-fit_pol1_1->GetParError(0)+ (fit_pol1_1->GetParameter(1) + fit_pol1_1->GetParError(1)) * x_fit[i];

      double error1= TMath::Min(e1,TMath::Min(e2,TMath::Min(e3,e4)));
      double error2= TMath::Max(e1,TMath::Max(e2,TMath::Max(e3,e4)));
      double errorb= (error2 - error1)/2.;
      double b=fit_pol1_1->GetParameter(0)+fit_pol1_1->GetParameter(1) * x_fit[i];

      double errora=ey_fit[i];
      double a=y_fit[i];

      yfit[i]= 0.;      
      eyfit[i] = 0.;
      xfit[i]= 0.;
      exfit[i]=0.;
      
      
      yfit[i]= 100*(a-b)/a;           
      eyfit[i] = yfit[i] * sqrt(TMath::Power(errora/a,2) + TMath::Power(errorb/b,2) ) ;
      xfit[i]=x_fit[i];
      exfit[i]=ex_fit[i];
	    
    }

  
    TGraphErrors* g_dev_grid24_conf1 = new TGraphErrors(n_fit,xfit,yfit,exfit,eyfit);

    
    //20/1 ------------------------------------------------
    TF1 *fit_pol1_2 = new TF1("fit_pol1_2","pol1",0.5,6);
    g_grid20_conf1->Fit(fit_pol1_2,"RMN");

    x_fit =  g_grid20_conf1->GetX();
    y_fit =  g_grid20_conf1->GetY();
    ex_fit = g_grid20_conf1->GetEX();
    ey_fit = g_grid20_conf1->GetEY();

    n_fit =  g_grid20_conf1->GetN();

    // double xfit[10], yfit[10], exfit[10], eyfit[10];

    for(int i=0; i<g_grid20_conf1->GetN(); i++) {
      double e1=  fit_pol1_2->GetParameter(0)+fit_pol1_2->GetParError(0)+ (fit_pol1_2->GetParameter(1) - fit_pol1_2->GetParError(1)) * x_fit[i];
      double e2=  fit_pol1_2->GetParameter(0)+fit_pol1_2->GetParError(0)+ (fit_pol1_2->GetParameter(1) + fit_pol1_2->GetParError(1)) * x_fit[i];
      double e3=  fit_pol1_2->GetParameter(0)-fit_pol1_2->GetParError(0)+ (fit_pol1_2->GetParameter(1) - fit_pol1_2->GetParError(1)) * x_fit[i];
      double e4=  fit_pol1_2->GetParameter(0)-fit_pol1_2->GetParError(0)+ (fit_pol1_2->GetParameter(1) + fit_pol1_2->GetParError(1)) * x_fit[i];

      double error1= TMath::Min(e1,TMath::Min(e2,TMath::Min(e3,e4)));
      double error2= TMath::Max(e1,TMath::Max(e2,TMath::Max(e3,e4)));
      double errorb= (error2 - error1)/2.;
      double b=fit_pol1_2->GetParameter(0)+fit_pol1_2->GetParameter(1) * x_fit[i];

      double errora=ey_fit[i];
      double a=y_fit[i];
      
      yfit[i]= 0.;      
      eyfit[i] = 0.;
      xfit[i]= 0.;
      exfit[i]=0.;
      
      yfit[i]= 100*(a-b)/a;           
      eyfit[i] = yfit[i] * sqrt(TMath::Power(errora/a,2) + TMath::Power(errorb/b,2) ) ;
      xfit[i]=x_fit[i];
      exfit[i]=ex_fit[i];
	    
    }

  
    TGraphErrors* g_dev_grid20_conf1 = new TGraphErrors(n_fit,xfit,yfit,exfit,eyfit);

      //20/2 ------------------------------------------------
    TF1 *fit_pol1_3 = new TF1("fit_pol1_3","pol1",0.5,6);
    g_grid20_conf2->Fit(fit_pol1_3,"RMN");

    x_fit =  g_grid20_conf2->GetX();
    y_fit =  g_grid20_conf2->GetY();
    ex_fit = g_grid20_conf2->GetEX();
    ey_fit = g_grid20_conf2->GetEY();

    n_fit =  g_grid20_conf2->GetN();

    // double xfit[10], yfit[10], exfit[10], eyfit[10];

    for(int i=0; i<g_grid20_conf2->GetN(); i++) {
      double e1=  fit_pol1_3->GetParameter(0)+fit_pol1_3->GetParError(0)+ (fit_pol1_3->GetParameter(1) - fit_pol1_3->GetParError(1)) * x_fit[i];
      double e2=  fit_pol1_3->GetParameter(0)+fit_pol1_3->GetParError(0)+ (fit_pol1_3->GetParameter(1) + fit_pol1_3->GetParError(1)) * x_fit[i];
      double e3=  fit_pol1_3->GetParameter(0)-fit_pol1_3->GetParError(0)+ (fit_pol1_3->GetParameter(1) - fit_pol1_3->GetParError(1)) * x_fit[i];
      double e4=  fit_pol1_3->GetParameter(0)-fit_pol1_3->GetParError(0)+ (fit_pol1_3->GetParameter(1) + fit_pol1_3->GetParError(1)) * x_fit[i];

      double error1= TMath::Min(e1,TMath::Min(e2,TMath::Min(e3,e4)));
      double error2= TMath::Max(e1,TMath::Max(e2,TMath::Max(e3,e4)));
      double errorb= (error2 - error1)/2.;
      double b=fit_pol1_3->GetParameter(0)+fit_pol1_3->GetParameter(1) * x_fit[i];

      double errora=ey_fit[i];
      double a=y_fit[i];
      
      yfit[i]= 0.;      
      eyfit[i] = 0.;
      xfit[i]= 0.;
      exfit[i]=0.;
      
      yfit[i]= 100*(a-b)/a;           
      eyfit[i] = yfit[i] * sqrt(TMath::Power(errora/a,2) + TMath::Power(errorb/b,2) ) ;
      xfit[i]=x_fit[i];
      exfit[i]=ex_fit[i];
	    
    }

  
    TGraphErrors* g_dev_grid20_conf2 = new TGraphErrors(n_fit,xfit,yfit,exfit,eyfit);

      //20/3 ------------------------------------------------
    TF1 *fit_pol1_4 = new TF1("fit_pol1_4","pol1",0.5,6);
    g_grid20_conf3->Fit(fit_pol1_4,"RMN");

    x_fit =  g_grid20_conf3->GetX();
    y_fit =  g_grid20_conf3->GetY();
    ex_fit = g_grid20_conf3->GetEX();
    ey_fit = g_grid20_conf3->GetEY();

    n_fit =  g_grid20_conf3->GetN();

    // double xfit[10], yfit[10], exfit[10], eyfit[10];

    for(int i=0; i<g_grid20_conf3->GetN(); i++) {
      double e1=  fit_pol1_4->GetParameter(0)+fit_pol1_4->GetParError(0)+ (fit_pol1_4->GetParameter(1) - fit_pol1_4->GetParError(1)) * x_fit[i];
      double e2=  fit_pol1_4->GetParameter(0)+fit_pol1_4->GetParError(0)+ (fit_pol1_4->GetParameter(1) + fit_pol1_4->GetParError(1)) * x_fit[i];
      double e3=  fit_pol1_4->GetParameter(0)-fit_pol1_4->GetParError(0)+ (fit_pol1_4->GetParameter(1) - fit_pol1_4->GetParError(1)) * x_fit[i];
      double e4=  fit_pol1_4->GetParameter(0)-fit_pol1_4->GetParError(0)+ (fit_pol1_4->GetParameter(1) + fit_pol1_4->GetParError(1)) * x_fit[i];

      double error1= TMath::Min(e1,TMath::Min(e2,TMath::Min(e3,e4)));
      double error2= TMath::Max(e1,TMath::Max(e2,TMath::Max(e3,e4)));
      double errorb= (error2 - error1)/2.;
      double b=fit_pol1_4->GetParameter(0)+fit_pol1_4->GetParameter(1) * x_fit[i];

      double errora=ey_fit[i];
      double a=y_fit[i];
      
      yfit[i]= 0.;      
      eyfit[i] = 0.;
      xfit[i]= 0.;
      exfit[i]=0.;
      
      yfit[i]= 100*(a-b)/a;           
      eyfit[i] = yfit[i] * sqrt(TMath::Power(errora/a,2) + TMath::Power(errorb/b,2) ) ;
      xfit[i]=x_fit[i];
      exfit[i]=ex_fit[i];
	    
    }

  
    TGraphErrors* g_dev_grid20_conf3 = new TGraphErrors(n_fit,xfit,yfit,exfit,eyfit);

    
    
    TCanvas *c_energy_linearity_2 = new TCanvas("c_energy_linearity_2","c_energy_linearity_2",800,600);
    //c_energy_linearity->Divide(2,1);
    c_energy_linearity_2->cd(1);

    g_dev_grid24_conf1->SetLineWidth(2);
    g_dev_grid24_conf1->SetLineStyle(2);
    g_dev_grid24_conf1->SetMarkerStyle(4);
    g_dev_grid24_conf1->GetYaxis()->SetTitle("deviation from linear fit [%]");
    g_dev_grid24_conf1->GetXaxis()->SetTitle("E^{beam}/GeV");
    g_dev_grid24_conf1->GetYaxis()->SetRangeUser(-10,15);
    g_dev_grid24_conf1->Draw("alp");

    g_dev_grid20_conf1->SetLineWidth(2);
    g_dev_grid20_conf1->SetLineStyle(1);
    g_dev_grid20_conf1->SetMarkerStyle(20);
    g_dev_grid20_conf1->Draw("lp");

    g_dev_grid20_conf2->SetLineWidth(2);
    g_dev_grid20_conf2->SetLineColor(8);
    g_dev_grid20_conf2->SetMarkerStyle(22);
    g_dev_grid20_conf2->SetMarkerColor(8);
    g_dev_grid20_conf2->Draw("lp");

    
    g_dev_grid20_conf3->SetLineWidth(2);
    g_dev_grid20_conf3->SetLineColor(4);
    g_dev_grid20_conf3->SetMarkerStyle(23);
    g_dev_grid20_conf3->SetMarkerColor(4);
    g_dev_grid20_conf3->Draw("lp");
    

    TLegend *l_energy_linearity_2 = new TLegend(0.2,0.7,0.7,0.9);
    //  l_energy_linearity_2->SetHeader("SiW-ECAL: reconstructed energy linearity");

    double chi1=fit_pol1_1->GetChisquare()/fit_pol1_1->GetNDF();
    double chi2=fit_pol1_2->GetChisquare()/fit_pol1_2->GetNDF();
    double chi3=fit_pol1_3->GetChisquare()/fit_pol1_3->GetNDF();
    double chi4=fit_pol1_4->GetChisquare()/fit_pol1_4->GetNDF();

    l_energy_linearity_2->AddEntry(g_dev_grid24_conf1,TString::Format("Wafer 3, W-configuration 1, #chi^{2}/ndf=%.1f",chi1),"lp");
    l_energy_linearity_2->AddEntry(g_dev_grid20_conf1,TString::Format("Wafer 4, W-configuration 1, #chi^{2}/ndf=%.1f",chi2),"lp");
    l_energy_linearity_2->AddEntry(g_dev_grid20_conf2,TString::Format("Wafer 4, W-configuration 2, #chi^{2}/ndf=%.1f",chi3),"lp");
    l_energy_linearity_2->AddEntry(g_dev_grid20_conf3,TString::Format("Wafer 4, W-configuration 3, #chi^{2}/ndf=%.1f",chi4),"lp");

    //    l_energy_linearity_2->AddEntry(g_dev_grid20_conf1,"Wafer 4, W-configuration 1","lp");
    //    l_energy_linearity_2->AddEntry(g_dev_grid20_conf2,"Wafer 4, W-configuration 2","lp");
    //    l_energy_linearity_2->AddEntry(g_dev_grid20_conf3,"Wafer 4, W-configuration 3","lp");
    l_energy_linearity_2->SetFillColor(0);
    l_energy_linearity_2->SetLineColor(0);
    l_energy_linearity_2->SetShadowColor(0);
    l_energy_linearity_2->Draw();

    
    c_energy_linearity_2->Print("zmlt1_5slabs/energy_linearity_dev_zcut_1_mipcut0.5.eps");
    c_energy_linearity_2->Print("zmlt1_5slabs/energy_linearity_dev_zcut_1_mipcut0.5.C");
    

}
