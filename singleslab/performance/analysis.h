#include "TFile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TH1.h"
#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TString.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "TLatex.h"

//#include "../../style/Style.C"
//#include "../../style/Labels.C"

using namespace std;
Float_t map_pointX[16][64];
Float_t map_pointY[16][64];
int nchips=16;

Double_t langaufun(Double_t *x, Double_t *par) {
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.
  
  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  
  
  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 
  
  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  
  step = (xupp-xlow) / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  
  return (par[2] * step * sum * invsq2pi / par[3]);
}



TF1* langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
   //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf

   Int_t i;
   Char_t FunName[100];

   sprintf(FunName,"Fitfcn_%s",his->GetName());

   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MP","Area","GSigma");
   
   for (i=0; i<4; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }

   his->Fit(FunName,"RBQM");   // fit within specified range, use ParLimits, quiet, improve fit results (TMINUIT)

   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<4; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf

   return (ffit);              // return fit function

}

void ReadMap(TString filename="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev10_chip_channel_x_y_mapping.txt")
{

  std::ifstream reading_file(filename);
  if(!reading_file){
    cout<<" dameyo - damedame"<<endl;
  }
  for(int i=0; i<16; i++) {
    for(int j=0; j<64; j++) {
      map_pointX[i][j] = -1000.;
      map_pointY[i][j] = -1000.;
    }
  }

  Int_t tmp_chip = 0,tmp_channel = 0;
  Float_t tmp_x0 = 0 ,tmp_y0 = 0 , tmp_x = 0 , tmp_y = 0 ;
  TString tmpst;
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst ;
  while(reading_file){
    reading_file >> tmp_chip >> tmp_x0 >> tmp_y0 >> tmp_channel >> tmp_x >> tmp_y ;
    map_pointX[tmp_chip][tmp_channel] = -tmp_x ;
    map_pointY[tmp_chip][tmp_channel] = -tmp_y ;
  }

}

void pedanalysis(TFile *file, TString run="Run_ILC_cosmic_test_11222019", int slboard=2){
  
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
  
  // Comparing nbr entries in tag or not tag // GetWidth and Mean
  TString map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev10_chip_channel_x_y_mapping.txt";
  if(slboard==0 || slboard==2)  map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev11_cob_chip_channel_x_y_mapping.txt";
  ReadMap(map);

  
  TFile *_file0 = TFile::Open(TString::Format("../results_pedestal/Pedestal_SLB_%i_%s.root",slboard,run.Data()));

 
  TH2F* PedW[15];
  TH2F* PedM[15];
  TH2F* PedNPeaks[15];
  TH2F* PedRMS[15];

  for(int i=0; i<15; i++) {
    PedW[i]=new TH2F(TString::Format("PedW_sca%i",i),TString::Format("PedW_sca%i",i),32,-90,90,32,-90,90);
    PedM[i]=new TH2F(TString::Format("PedM_sca%i",i),TString::Format("PedM_sca%i",i),32,-90,90,32,-90,90);

    PedRMS[i]=new TH2F(TString::Format("PedRMS_sca%i",i),TString::Format("PedRMS_sca%i",i),32,-90,90,32,-90,90);
    PedNPeaks[i]=new TH2F(TString::Format("PedNPeaks_sca%i",i),TString::Format("PedNPeaks_sca%i",i),32,-90,90,32,-90,90);
  }

  TH2F* PedW_sca;
  TH2F* PedM_sca;
  PedW_sca=new TH2F("PedW_sca","PedW_sca",16,0,15,15,0,14);
  PedM_sca=new TH2F("PedM_sca","PedM_sca",16,0,15,15,0,14);

  for(int n=0; n<15;n++) {
    for(int i=0;i<nchips;i++){

      double nch=0;
      float avmean=0, avwidth=0;
      
      for(int j=0; j<64; j++) {
	//GetGoodEntries
	TH1F *temp=(TH1F*)_file0->Get(TString::Format("ped_chip%i_chn%i_sca%i",i,j,n));
	
	if(temp->GetEntries()> 200 ){ //max_entries/2 ) {
          TSpectrum *s = new TSpectrum();
          int npeaks = s->Search(temp,2,"",0.8); 
          if(npeaks > 0) {
            Double_t *mean_peak=s->GetPositionX();
            Double_t *mean_high=s->GetPositionY();
            double mean_peak_higher=0;
            double mean_high_higher=0;
            int npeak_max=0;
            for(int ipeak=0; ipeak<npeaks; ipeak ++) {
              if(mean_high[ipeak]>mean_high_higher && mean_high[ipeak]>50) {
                mean_high_higher=mean_high[ipeak];
                mean_peak_higher=mean_peak[ipeak];
                npeak_max=ipeak;
              }
            }
	    PedNPeaks[n]->Fill(map_pointX[i][j],map_pointY[i][j],npeaks);
	    
            if(npeaks ==1 ) {
	      nch++;

              Double_t *mean_peak=s->GetPositionX();
              mean_peak[0]=mean_peak_higher;
              
              TF1 *f0 = new TF1("f0","gaus",mean_peak[0]-temp->GetRMS()/2,mean_peak[0]+temp->GetRMS()/2);
              temp->Fit("f0","RQNOC");
              
              TF1 *f1 = new TF1("f1","gaus",f0->GetParameter(1)-f0->GetParameter(2)/2,f0->GetParameter(1)+f0->GetParameter(2)/2);
              temp->Fit("f1","RQME");
	      avmean+=f1->GetParameter(1);
	      avwidth+=f1->GetParameter(2);
	      PedW[n]->Fill(map_pointX[i][j],map_pointY[i][j],f1->GetParameter(2));
	      PedRMS[n]->Fill(map_pointX[i][j],map_pointY[i][j],temp->GetRMS());
	      PedM[n]->Fill(map_pointX[i][j],map_pointY[i][j],f1->GetParameter(1));
            }
	  }
       }

      }
      
      if(nch>20) {
	PedW_sca->Fill(i,n,avwidth/nch);
	PedM_sca->Fill(i,n,avmean/nch);
      }
    }

    
  }

  //  _file0 ->Close();

  TCanvas* canvas= new TCanvas(TString::Format("PedAna_%s_SLB_%i",run.Data(),slboard),TString::Format("PedAna_%s_SLB_%i",run.Data(),slboard),1600,1600);   
  canvas->Divide(4,4);
  
  for(int i=0; i<15; i++) {
    canvas->cd(1+i);
    PedM[i]->GetXaxis()->SetTitle("x");
    PedM[i]->GetYaxis()->SetTitle("y");
    PedM[i]->GetZaxis()->SetRangeUser(200,400);
    PedM[i]->Draw("colz");

    canvas->cd(5+i);
    PedW[i]->GetXaxis()->SetTitle("x");
    PedW[i]->GetYaxis()->SetTitle("y");
    PedW[i]->GetZaxis()->SetRangeUser(2,12);
    PedW[i]->Draw("COLZ");
    
    canvas->cd(9+i);
    PedNPeaks[i]->GetXaxis()->SetTitle("x");
    PedNPeaks[i]->GetYaxis()->SetTitle("y");
    PedNPeaks[i]->GetZaxis()->SetRangeUser(0,5);
    PedNPeaks[i]->Draw("COLZ");

    canvas->cd(13+i);
    PedRMS[i]->GetXaxis()->SetTitle("x");
    PedRMS[i]->GetYaxis()->SetTitle("y");
    PedRMS[i]->GetZaxis()->SetRangeUser(2,12);
    PedRMS[i]->Draw("COLZ");
    
  }
  canvas->Print(TString::Format("plots/PedAna_%s_SLB_%i.eps",run.Data(),slboard));

  TCanvas* canvas2= new TCanvas(TString::Format("PedAna2_%s_SLB_%i",run.Data(),slboard),TString::Format("PedAna2_%s_SLB_%i",run.Data(),slboard),1200,600);   
  canvas2->Divide(2,1);
  
  canvas2->cd(1);
  PedM_sca->GetXaxis()->SetTitle("ASIC");
  PedM_sca->GetYaxis()->SetTitle("SCA");
  PedM_sca->GetZaxis()->SetRangeUser(200,400);
  PedM_sca->Draw("colz");

  canvas2->cd(2);
  PedW_sca->GetXaxis()->SetTitle("ASIC");
  PedW_sca->GetYaxis()->SetTitle("SCA");
  PedW_sca->GetZaxis()->SetRangeUser(2,12);
  PedW_sca->Draw("colz");
   
  canvas2->Print(TString::Format("plots/PedAna2_%s_SLB_%i.eps",run.Data(),slboard));

  file->cd();
  canvas->Write();
  canvas->Write();
  for(int i=0; i<15; i++) {
    PedM[i]->Write();
    PedW[i]->Write();
    PedNPeaks[i]->Write();
    PedRMS[i]->Write();
  }
  PedM_sca->Write();
  PedW_sca->Write();
  

}


void mipanalysis(TFile* file, TString run="Run_ILC_cosmic_test_11222019", int slboard=2){
  
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


  file->cd();
  TDirectory *cdhisto = file->GetDirectory("histograms_mip");
  if(!cdhisto) cdhisto = file->mkdir("histograms_mip");
  cdhisto->cd();

    // Comparing nbr entries in tag or not tag // GetWidth and Mean
  TString map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev10_chip_channel_x_y_mapping.txt";
  if(slboard==0 || slboard==2)  map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev11_cob_chip_channel_x_y_mapping.txt";
  ReadMap(map);

  
  TFile *_file0 = TFile::Open(TString::Format("../results_mipcalibration/Signal_summary_SLB_%i_%s.root",slboard,run.Data()));

 
  TH2F* MIPW =new TH2F("MIPW","MIPW",32,-90,90,32,-90,90);
  TH2F* MIPM=new TH2F("MIPM","MIPM",32,-90,90,32,-90,90);
  TH2F* MIPRMS=new TH2F("MIPRMS","MIPRMS",32,-90,90,32,-90,90);
  TH2F* MIPN=new TH2F("MIPN","MIPN",32,-90,90,32,-90,90);

  TCanvas* canvas_mip[16];
  for(int i=0;i<nchips;i++){
    
    canvas_mip[i]= new TCanvas(TString::Format("MIPs_%s_SLB_%i_chip%i",run.Data(),slboard,i),TString::Format("MIPs_%s_SLB_%i_chip%i",run.Data(),slboard,i),1600,1600);   
    canvas_mip[i]->Divide(8,8);

      for(int j=0; j<64; j++) {
	TCanvas* canvastemp= new TCanvas("temp","temp",100,100);   
	canvastemp->cd();
	//GetGoodEntries
	TH1F *temp=(TH1F*)_file0->Get(TString::Format("histograms/charge__SLB_%i_%s_chip%i_chn%i",slboard,run.Data(),i,j));
	if(temp==NULL) continue;
	temp->Rebin(5);

	MIPN->Fill(map_pointX[i][j],map_pointY[i][j],temp->GetEntries());

	 Double_t fr[2];
	 Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];   
	 pllo[0]=1.0; pllo[1]=20; pllo[2]=100.0; pllo[3]=1;
	 plhi[0]=100.0; plhi[1]=100.0; plhi[2]=100000000.0; plhi[3]=30.0;
	 sv[0]=15.0;
	 Double_t chisqr;
	 Int_t    ndf;

	 
	 if(temp->GetEntries()>100){
	   fr[0]=30;
	   fr[1]=150;//fr[0]+0.5*temp->GetRMS();
	   TF1 *fitsnr_temp=langaufit(temp,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
	
	   double mpv=fitsnr_temp->GetParameter(1);
	   double empv=fitsnr_temp->GetParError(1);
	   double wmpv=fitsnr_temp->GetParameter(0);
	   double chi2ndf=0;
	   if(ndf>0) chi2ndf=chisqr/ndf;

	   MIPM->Fill(map_pointX[i][j],map_pointY[i][j],mpv);
	   MIPRMS->Fill(map_pointX[i][j],map_pointY[i][j],fitsnr_temp->GetParameter(3));
	   MIPW->Fill(map_pointX[i][j],map_pointY[i][j],wmpv);

	   canvas_mip[i]->cd(j+1);
	   temp->GetXaxis()->SetRangeUser(0,200);
	   temp->GetXaxis()->SetLabelSize(0.1);
	   temp->GetYaxis()->SetLabelSize(0.1);
	   float ymax=temp->GetMaximum();
	   temp->GetYaxis()->SetRangeUser(0,ymax*1.5);
	   temp->SetName(TString::Format("mip_chip%i_chn%i",i,j));
	   temp->Draw();
	   cdhisto->cd();
	   temp->Write();
	   TLatex t;
	   t.SetTextSize(0.15);
	   t.SetTextAlign(13);  //align at top
 	   t.DrawLatex(10,ymax*1.5,TString::Format("N=%i",int(temp->GetEntries())));
	   t.DrawLatex(10,ymax*1.15,TString::Format("MPV=%.1f",mpv));
	 } 
      }
      file->cd();
      canvas_mip[i]->Print(TString::Format("plots/MIPs_%s_SLB_%i_chip%i.eps",run.Data(),slboard,i));
      canvas_mip[i]->Write();
      
  }



  //  _file0 ->Close();

  TCanvas* canvas= new TCanvas(TString::Format("MIPAna_%s_SLB_%i",run.Data(),slboard),TString::Format("MIPAna_%s_SLB_%i",run.Data(),slboard),1600,1600);   
  canvas->Divide(2,2);
  
  canvas->cd(1);
    MIPM->GetXaxis()->SetTitle("x");
    MIPM->GetYaxis()->SetTitle("y");
    MIPM->GetZaxis()->SetRangeUser(50,150);
    MIPM->Draw("colz");

    canvas->cd(2);
    MIPW->GetXaxis()->SetTitle("x");
    MIPW->GetYaxis()->SetTitle("y");
    MIPW->GetZaxis()->SetRangeUser(1,30);
    MIPW->Draw("COLZ");
    
    canvas->cd(3);
    gPad->SetLogz();
    MIPN->GetXaxis()->SetTitle("x");
    MIPN->GetYaxis()->SetTitle("y");
    MIPN->GetZaxis()->SetRangeUser(50,10000);
    MIPN->Draw("COLZ");

    canvas->cd(4);
    MIPRMS->GetXaxis()->SetTitle("x");
    MIPRMS->GetYaxis()->SetTitle("y");
    MIPRMS->GetZaxis()->SetRangeUser(5,200);
    MIPRMS->Draw("COLZ");

  canvas->Print(TString::Format("plots/MIPAna_%s_SLB_%i.eps",run.Data(),slboard));

  file->cd();
  MIPM->Write();
  MIPW->Write();
  MIPN->Write();
  MIPRMS->Write();
  canvas->Write();
  
}


void retriggeranalysis(TFile* file, TString run="Run_ILC_cosmic_test_11222019", int slboard=2){
  
  gROOT->Reset();
  //SetIrlesStyle();
  //  gROOT->LoadMacro("Labels.C");
  gStyle->SetOptFit(0); 
  gStyle->SetOptStat(1);
  gStyle->SetOptTitle(1);
  
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.4);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetMarkerSize(1.2);


  file->cd();
  TDirectory *cdhisto = file->GetDirectory("histograms_mip");
  if(!cdhisto) cdhisto = file->mkdir("histograms_mip");
  cdhisto->cd();

    // Comparing nbr entries in tag or not tag // GetWidth and Mean
  TString map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev10_chip_channel_x_y_mapping.txt";
  if(slboard==0 || slboard==2)  map="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis/mapping/fev11_cob_chip_channel_x_y_mapping.txt";
  ReadMap(map);

  
  TFile *_file0 = TFile::Open(TString::Format("../results_retriggers/Retriggers_SLB_%i_%s.root",slboard,run.Data()));

   
  TH2F* f_ret_xy= (TH2F*)_file0->Get(TString::Format("first_retriggering_xy_SLB_%i_%s",slboard,run.Data()));
  TH2F* ret_xy= (TH2F*)_file0->Get(TString::Format("all_retriggering_xy_SLB_%i_%s",slboard,run.Data()));
  TH2F* trig_xy= (TH2F*)_file0->Get(TString::Format("trig_xy_SLB_%i_%s",slboard,run.Data()));

  TH2F* f_ret= (TH2F*)_file0->Get(TString::Format("first_retriggering_SLB_%i_%s",slboard,run.Data()));
  TH2F* ret= (TH2F*)_file0->Get(TString::Format("all_retriggering_SLB_%i_%s",slboard,run.Data()));
  TH2F* trig= (TH2F*)_file0->Get(TString::Format("trig_SLB_%i_%s",slboard,run.Data()));

  TCanvas* canvas= new TCanvas(TString::Format("Ret_%s_SLB_%i",run.Data(),slboard),TString::Format("Ret_%s_SLB_%i",run.Data(),slboard),1200,1600);   
  canvas->Divide(2,3);
  
  canvas->cd(1);
  trig_xy->Draw("colz");
  canvas->cd(2);
  trig->Draw("colz");

  canvas->cd(3);
  f_ret_xy->Draw("colz");
  canvas->cd(4);
  f_ret->Draw("colz");

  canvas->cd(5);
  ret_xy->Draw("colz");
  canvas->cd(6);
  ret->Draw("colz"); 

  canvas->Print(TString::Format("plots/Ret_%s_SLB_%i.eps",run.Data(),slboard));
  
  file->cd();
  trig_xy->Write();
  trig->Write();
  f_ret_xy->Write();
  f_ret->Write();
  ret_xy->Write();
  ret->Write();
  canvas->Write();
  
}
