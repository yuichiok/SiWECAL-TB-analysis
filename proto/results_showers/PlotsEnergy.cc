#include "TROOT.h"
#include "TFile.h"
#include "../../style/Style.C"
#include "TGraphErrors.h"
#include "TMath.h"


//  for(int iz=1; iz<2; iz++) {
//    TString zm="zm0.5";
//     if(iz==1) zm="zm1.0";
//    if(iz==2) zm="zm1.0";

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

  for(int islabs=4;islabs<8; islabs++) {
    TString nslabs=TString::Format("pedestal_nslabs%i",islabs);
      
    for(int ibcid=0; ibcid<2; ibcid++) {
      TString bcid="bcidmax2850";
      if(ibcid==1) bcid="bcidmax2850_removedIsolatedHits";
      
      for(int icases=0; icases<4; icases++) {
	if(icases==0) {
	  grid="grid24";
	  conf="conf1";
	}

	if(icases==1) grid="grid20";
	if(icases==2) conf="conf2";
	if(icases==3) conf="conf3";
	
	
	TString s_file=nslabs+"_"+bcid+"/"+conf+"_"+grid+"_2GeV_mipcut1.0_showers.root";
	std::cout<<"Opening file: "<<s_file<<std::endl;
	TFile *file = new TFile(s_file);
	
	TH1F *energy = (TH1F*)file->Get("energy");
	TH1F *energy_center_zm1 = (TH1F*)file->Get("energy_center_zm2");
	
	energy->Rebin(4);
	energy_center_zm1->Rebin(4);
	
	//  file->Close();
	
	TString s_file_2=nslabs+"_"+bcid+"/"+conf+"_"+grid+"_4GeV_mipcut1.0_showers.root";
	
	std::cout<<"Opening file: "<<s_file_2<<std::endl;
	TFile *file_2 = new TFile(s_file_2);
	
	TH1F *energy_2 = (TH1F*)file_2->Get("energy");
	TH1F *energy_center_zm1_2 = (TH1F*)file_2->Get("energy_center_zm2");
	
	energy_2->Rebin(4);
	energy_center_zm1_2->Rebin(4);
	
	TCanvas *c_energy = new TCanvas("c_energy","c_energy",800,600);
	//c_energy->Divide(2,1);
	c_energy->cd(1);
	// energy->SetTitle("SiW-ECAL (wafer-3), e^{+} beam");
	energy->GetXaxis()->SetTitle("E^{raw}/MIP = #sum_{i=layers} #sum_{j=cells} #omega_{i} E_{i,j}^{raw}");
	energy->GetYaxis()->SetTitle("# entries");
	energy->Draw();
	
	energy_center_zm1->SetLineColor(4);
	energy_center_zm1->Draw("same");
	
	energy_2->SetLineStyle(2);
	energy_2->Draw("same");
	
	energy_center_zm1_2->SetLineColor(4);
	energy_center_zm1_2->SetLineStyle(2);
	energy_center_zm1_2->Draw("same");

	double xmin=energy->GetMean()-2*energy->GetRMS();
	double xmax=energy->GetMean()+2*energy->GetRMS();
	double xmin2=TMath::Max(xmin,75.);//TMath::Max(e1,TMath::Max(e2,TMath::Max(e3,e4)));
	TF1 *fit_1 = new TF1("fit_1","gaus",xmin2,xmax);
	energy->Fit(fit_1,"MERQN");
	xmin=fit_1->GetParameter(1)-1.*fit_1->GetParameter(2);
	xmax=fit_1->GetParameter(1)+1.*fit_1->GetParameter(2);
	TF1 *fit_2 = new TF1("fit_2","gaus",xmin,xmax);
	fit_2->SetLineColor(1);
	fit_2->SetLineWidth(2);
	energy->Fit(fit_1,"MERQ");

	energy_2->GetXaxis()->SetRangeUser(125,1000);
	xmin=energy_2->GetMean()-2.*energy_2->GetRMS();
	xmax=energy_2->GetMean()+2.*energy_2->GetRMS();
	TF1 *fit_3 = new TF1("fit_3","gaus",xmin,xmax);
	energy_2->Fit(fit_3,"MERQN");
	xmin=fit_3->GetParameter(1)-1*fit_3->GetParameter(2);
	xmax=fit_3->GetParameter(1)+1*fit_3->GetParameter(2);
	TF1 *fit_4 = new TF1("fit_4","gaus",xmin,xmax);
	fit_4->SetLineColor(1);
	fit_4->SetLineStyle(2);
	fit_4->SetLineWidth(2);
	energy_2->Fit(fit_4,"MERQ");
	energy_2->GetXaxis()->SetRangeUser(0,1000);

	
	
	//--------

	xmin=energy_center_zm1->GetMean()-1*energy_center_zm1->GetRMS();
	xmax=energy_center_zm1->GetMean()+1*energy_center_zm1->GetRMS();
	TF1 *fit_center1 = new TF1("fit_center1","gaus",xmin,xmax);
	fit_center1->SetLineColor(4);
	fit_center1->SetLineWidth(2);
	energy_center_zm1->Fit(fit_center1,"MERQ");

	energy_center_zm1_2->GetXaxis()->SetRangeUser(75,1000);
	xmin=energy_center_zm1_2->GetMean()-2*energy_center_zm1_2->GetRMS();
	xmax=energy_center_zm1_2->GetMean()+2*energy_center_zm1_2->GetRMS();
	energy_center_zm1_2->GetXaxis()->SetRangeUser(xmin,xmax);
	xmin=energy_center_zm1_2->GetMean()-2*energy_center_zm1_2->GetRMS();
	xmax=energy_center_zm1_2->GetMean()+2*energy_center_zm1_2->GetRMS();
	
	TF1 *fit_center2 = new TF1("fit_center2","gaus",xmin,xmax);
	energy_center_zm1_2->Fit(fit_center2,"MERQN");
	xmin=fit_center2->GetParameter(1)-1*fit_center2->GetParameter(2);
	xmax=fit_center2->GetParameter(1)+1*fit_center2->GetParameter(2);
	TF1 *fit_center3 = new TF1("fit_center3","gaus",xmin,xmax);
	fit_center3->SetLineColor(4);
	fit_center3->SetLineStyle(2);
	fit_center3->SetLineWidth(2);
	energy_center_zm1_2->Fit(fit_center3,"MERQ");
	energy_center_zm1_2->GetXaxis()->SetRangeUser(0,1000);

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
	l_energy->AddEntry(energy_center_zm1,TString::Format("e^{+} 2 GeV, contained showers:  gaussian fit #chi^{2}/ndf=%0.1f",chi3),"l");
	l_energy->AddEntry(energy_2,TString::Format("e^{+} 4 GeV: gaussian fit #chi^{2}/ndf=%0.1f",chi2),"l");
	l_energy->AddEntry(energy_center_zm1_2,TString::Format("e^{+} 4 GeV, contained showers: gaussian fit #chi^{2}/ndf=%0.1f",chi4),"l");
	l_energy->SetFillColor(0);
	l_energy->SetLineColor(0);
	l_energy->SetShadowColor(0);
	l_energy->Draw();
	
	c_energy->Print("shower_plots_zm2/energy_"+nslabs+"_"+bcid+"_"+conf+"_"+grid+"_"+"_mipcut1.0.eps");
	
      }
         
      /// ####################################################
      // ENERGY LINEARITY
      
      grid="grid24";
      conf="conf1";
      
      TGraphErrors * g_grid24_conf1;
      TGraphErrors * g_grid20_conf1;
      TGraphErrors * g_grid20_conf2;
      TGraphErrors * g_grid20_conf3;
      
      
      for(int icases=0; icases<4; icases++) {
	if(icases==0) {
	  grid="grid24";
	  conf="conf1";
	}
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
	  
	  
	  TString s_file=nslabs+"_"+bcid+"/"+conf+"_"+grid+"_"+s_energy+"_mipcut1.0_showers.root";
	  std::cout<<"Opening file: "<<s_file<<std::endl;
	  TFile *file = new TFile(s_file);
	  
	  TH1F *energy_center_zm1_fit = (TH1F*)file->Get("energy_center_zm2");
	  //energy_center_zm1_fit->Rebin(2);
	  
	  double xmin=energy_center_zm1_fit->GetMean()-2*energy_center_zm1_fit->GetRMS();
	  double xmax=energy_center_zm1_fit->GetMean()+2*energy_center_zm1_fit->GetRMS();
	  energy_center_zm1_fit->GetXaxis()->SetRangeUser(xmin,xmax);
	  xmin=energy_center_zm1_fit->GetMean()-2*energy_center_zm1_fit->GetRMS();
	  xmax=energy_center_zm1_fit->GetMean()+2*energy_center_zm1_fit->GetRMS();
	  
	  TF1 *fit_center = new TF1("fit_center","gaus",xmin,xmax);
	  energy_center_zm1_fit->Fit(fit_center,"MERQN");
	  xmin=fit_center->GetParameter(1)-1*fit_center->GetParameter(2);
	  xmax=fit_center->GetParameter(1)+1*fit_center->GetParameter(2);
	  TF1 *fit_center2 = new TF1("fit_center2","gaus",xmin,xmax);
	  energy_center_zm1_fit->Fit(fit_center2,"MERQN");
	  energy_center_zm1_fit->GetXaxis()->SetRangeUser(0,1000);
	
  
	  y[np]=fit_center2->GetParameter(1);
	  ey[np]=fit_center2->GetParameter(2)/sqrt(fit_center2->GetHistogram()->GetEntries());//Error(1);
	  
	  TAxis *xaxis = energy_center_zm1_fit->GetXaxis();
	  TAxis *yaxis = energy_center_zm1_fit->GetYaxis();
	  Int_t binxmin = xaxis->FindBin(xmin);
	  Int_t binxmax = xaxis->FindBin(xmax);
	  
	  
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

      c_energy_linearity->Print("shower_plots_zm2/"+nslabs+"_"+bcid+"_"+"energy_linearity_mipcut1.0.eps");
      // c_energy_linearity->Print("shower_plots/"+nslabs+"_"+bcid+"_"+"energy_linearity_mipcut1.0.C");


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
	eyfit[i] = 100*ey_fit[i]/a;// sqrt(TMath::Power(errora/a,2) + TMath::Power(errorb/b,2) ) ;
	xfit[i]=x_fit[i];
	exfit[i]=ex_fit[i];
	    
      }

  
      TGraphErrors* g_dev_grid24_conf1 = new TGraphErrors(n_fit,xfit,yfit,exfit,eyfit);

    
      //20/1 ------------------------------------------------
      TF1 *fit_pol1_2 = new TF1("fit_pol1_2","pol1",0.5,6);
      g_grid20_conf1->Fit(fit_pol1_2,"RMNQE");

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
	eyfit[i] = 100*ey_fit[i]/a;// sqrt(TMath::Power(errora/a,2) + TMath::Power(errorb/b,2) ) ;
	xfit[i]=x_fit[i];
	exfit[i]=ex_fit[i];
	    
      }

  
      TGraphErrors* g_dev_grid20_conf1 = new TGraphErrors(n_fit,xfit,yfit,exfit,eyfit);

      //20/2 ------------------------------------------------
      TF1 *fit_pol1_3 = new TF1("fit_pol1_3","pol1",0.5,6.);
      g_grid20_conf2->Fit(fit_pol1_3,"RMNQEI");

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
	double errorb= (error2 - error1);
	double b=fit_pol1_3->GetParameter(0)+fit_pol1_3->GetParameter(1) * x_fit[i];

	double errora=ey_fit[i];
	double a=y_fit[i];
      
	yfit[i]= 0.;      
	eyfit[i] = 0.;
	xfit[i]= 0.;
	exfit[i]=0.;
      
	eyfit[i] = 100*ey_fit[i]/a;// sqrt(TMath::Power(errora/a,2) + TMath::Power(errorb/b,2) ) ;
	yfit[i]= 100*(a-b)/a;           

	xfit[i]=x_fit[i];
	exfit[i]=ex_fit[i];
	    
      }

  
      TGraphErrors* g_dev_grid20_conf2 = new TGraphErrors(n_fit,xfit,yfit,exfit,eyfit);

      //20/3 ------------------------------------------------

    
      TF1 *fit_pol1_4 = new TF1("fit_pol1_4","pol1",0.5,6);
      g_grid20_conf3->Fit(fit_pol1_4,"RMNE");

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
	eyfit[i] = 100*ey_fit[i]/a;// sqrt(TMath::Power(errora/a,2) + TMath::Power(errorb/b,2) ) ;
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
      g_dev_grid24_conf1->Draw("alpe0");

      g_dev_grid20_conf1->SetLineWidth(2);
      g_dev_grid20_conf1->SetLineStyle(1);
      g_dev_grid20_conf1->SetMarkerStyle(20);
      g_dev_grid20_conf1->Draw("lpe0");

      g_dev_grid20_conf2->SetLineWidth(2);
      g_dev_grid20_conf2->SetLineColor(8);
      g_dev_grid20_conf2->SetMarkerStyle(22);
      g_dev_grid20_conf2->SetMarkerColor(8);
      g_dev_grid20_conf2->Draw("lpe0");

    
      g_dev_grid20_conf3->SetLineWidth(2);
      g_dev_grid20_conf3->SetLineColor(4);
      g_dev_grid20_conf3->SetMarkerStyle(23);
      g_dev_grid20_conf3->SetMarkerColor(4);
      g_dev_grid20_conf3->Draw("lpe0");
    

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

    
      c_energy_linearity_2->Print("shower_plots_zm2/"+nslabs+"_"+bcid+"_"+"energy_linearity_dev_mipcut1.0.eps");
      //  c_energy_linearity_2->Print("shower_plots/"+nslabs+"_"+bcid+"_"+conf+"_"+grid+"_"+"energy_linearity_dev_mipcut0.5.C");

    }
  }

}
