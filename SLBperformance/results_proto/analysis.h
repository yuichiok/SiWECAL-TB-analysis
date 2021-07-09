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

void ReadMap(TString filename="../../mapping/fev10_chip_channel_x_y_mapping.txt")
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


void triggers( TString run="Run_ILC_cosmic_test_11222019", int layer=2, int cob=0){
  
  gROOT->Reset();
  //SetIrlesStyle();
  //  gROOT->LoadMacro("Labels.C");
  gStyle->SetOptFit(0); 
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.4);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetPadRightMargin(0.2);
  //gStyle->SetPalette(kTemperatureMap);
  ////gStyle->SetPalette(kThermometer);
  gStyle->SetPalette(kLightTemperature);

    // Comparing nbr entries in tag or not tag // GetWidth and Mean
  TString map="../../mapping/fev10_chip_channel_x_y_mapping.txt";
  if(cob==1)  map="../../mapping/fev11_cob_chip_channel_x_y_mapping.txt";
  ReadMap(map);

  TFile *_file0 = TFile::Open(TString::Format("MIPs_15_layers_%s.root",run.Data()));
  
  TH1F *trig=(TH1F*)_file0->Get(TString::Format("layer_%i/trig_layer_%i",layer,layer));
  TH1F *trig_xy=(TH1F*)_file0->Get(TString::Format("layer_%i/trig_xy_layer_%i",layer,layer));
  TH1F *all_retriggering=(TH1F*)_file0->Get(TString::Format("layer_%i/all_retriggering_layer_%i",layer,layer));
  TH1F *all_retriggering_xy=(TH1F*)_file0->Get(TString::Format("layer_%i/all_retriggering_xy_layer_%i",layer,layer));
  

  TCanvas* canvas= new TCanvas(TString::Format("Triggers_layer_%i",layer),TString::Format("Triggers_layer_%i",layer),1400,1200);   
  canvas->Divide(2,2);
  
  canvas->cd(1);
  trig->GetXaxis()->SetTitle("CHIP");
  trig->GetYaxis()->SetTitle("chn");
  float nchn=1024.;
  if(cob==1) nchn=256.;
  float zmin=trig->GetEntries()/nchn-5*trig->GetEntries()/nchn;
  float zmax=trig->GetEntries()/nchn+5*trig->GetEntries()/nchn;
  float n_average=0;
  float average=0;
  for(int i=0; i<32; i++) {
    for(int j=0; j<32; j++) {
      if(trig->GetBinContent(i+1,j+1)<zmax && trig->GetBinContent(i+1,j+1)>zmin)  {
	average+=trig->GetBinContent(i+1,j+1);
	n_average++;
      }
    }
  }
  average/=n_average;
  zmin=average*(1-3);
  zmax=average*(1+3);
  if(zmin<0) zmin =0;
  if(cob==1) zmax=550;

  trig->GetZaxis()->SetRangeUser(zmin,zmax);
  trig->Draw("colz");
  
  canvas->cd(2);
  // gPad->SetLogz();
  trig_xy->GetXaxis()->SetTitle("x");
  trig_xy->GetYaxis()->SetTitle("y");
  trig_xy->GetZaxis()->SetRangeUser(zmin,zmax);
  trig_xy->Draw("COLZ");
  
  canvas->cd(3);
  zmin=all_retriggering->GetEntries()/nchn-5*all_retriggering->GetEntries()/nchn;
  zmax=all_retriggering->GetEntries()/nchn+5*all_retriggering->GetEntries()/nchn;
  n_average=0;
  average=0;
  for(int i=0; i<32; i++) {
    for(int j=0; j<32; j++) {
      if(all_retriggering->GetBinContent(i+1,j+1)<zmax && all_retriggering->GetBinContent(i+1,j+1)>zmin)  {
	average+=all_retriggering->GetBinContent(i+1,j+1);
	n_average++;
      }
    }
  }
  average/=n_average;
  zmin=average*(1-3);
  zmax=average*(1+3);
  if(zmin<0) zmin =0;
  
  
  all_retriggering->GetXaxis()->SetTitle("CHIP");
  all_retriggering->GetYaxis()->SetTitle("chn");
  //all_retriggering->GetZaxis()->SetRangeUser(50,150);
  all_retriggering->GetZaxis()->SetRangeUser(zmin,zmax);
  all_retriggering->Draw("colz");
  
  canvas->cd(4);
  // gPad->SetLogz();

  all_retriggering_xy->GetXaxis()->SetTitle("x");
  all_retriggering_xy->GetYaxis()->SetTitle("y");
  //all_retriggering_xy->GetZaxis()->SetRangeUser(1,30);
  all_retriggering_xy->GetZaxis()->SetRangeUser(zmin,zmax);
  all_retriggering_xy->Draw("COLZ");

  canvas->Print(TString::Format("plots/Triggers_%s_layer_%i.eps",run.Data(),layer));
  
}



void mipanalysis(TFile* file, TString run="Run_ILC_cosmic_test_11222019", int layer=2, int cob=0){
  
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
  gStyle->SetPadRightMargin(0.2);
  gStyle->SetPalette(kLightTemperature);

  file->cd();
  TDirectory *cdhisto = file->GetDirectory("histograms_mip");
  if(!cdhisto) cdhisto = file->mkdir("histograms_mip");
  cdhisto->cd();

    // Comparing nbr entries in tag or not tag // GetWidth and Mean
  TString map="../../mapping/fev10_chip_channel_x_y_mapping.txt";
  if(cob==1)  map="../../mapping/fev11_cob_chip_channel_x_y_mapping.txt";
  ReadMap(map);

  
  TFile *_file0 = TFile::Open(TString::Format("MIPs_15_layers_%s.root",run.Data()));

  cout<< TString::Format("MIPs_15_layers_%s.root",run.Data()) <<endl;
  TH2F* MIPW =new TH2F("MIPW","MIPW",32,-90,90,32,-90,90);
  TH2F* MIPM_xy=new TH2F("MIPM_xy","MIPM_xy",32,-90,90,32,-90,90);
  TH2F* MIPM=new TH2F(TString::Format("MIPM_slbooard_%i",layer),TString::Format("MIPM_slbooard_%i",layer),16,-0.5,15.5,64,-0.5,63.5);
  TH2F* MIPRMS=new TH2F("MIPRMS","MIPRMS",32,-90,90,32,-90,90);
  TH2F* MIPN=new TH2F("MIPN","MIPN",32,-90,90,32,-90,90);

  ofstream fout_mip(TString::Format("MIPs_layer_%i_%s",layer,run.Data()).Data(),ios::out);
  fout_mip<<"#mip results PROTO15"<<endl;
  fout_mip<<"#layer chip channel mpv empv widthmpv chi2ndf nentries"<<endl;

  for(int i=0;i<nchips;i++){
    
    TCanvas* canvas_mip= new TCanvas(TString::Format("MIPs_layer_%i_chip%i",layer,i),TString::Format("MIPs_layer_%i_chip%i",layer,i),1600,1600);   
      canvas_mip->Divide(8,8);

      for(int j=0; j<64; j++) {
	TCanvas* canvastemp= new TCanvas("temp","temp",100,100);   
	canvastemp->cd();
	//GetGoodEntries
	TH1F *temp=(TH1F*)_file0->Get(TString::Format("layer_%i/charge_layer%i_chip%i_chn%i",layer,layer,i,j));
	if(temp==NULL) continue;
	temp->Rebin(4);

	MIPN->Fill(map_pointX[i][j],map_pointY[i][j],temp->GetEntries());

	 Double_t fr[2];
	 Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];   
	 pllo[0]=1.0; pllo[1]=15; pllo[2]=1.0; pllo[3]=1;
	 plhi[0]=100.0; plhi[1]=100.0; plhi[2]=100000000.0; plhi[3]=20.0;
	 sv[0]=15.0;
	 Double_t chisqr;
	 Int_t    ndf;

	 
	 if(temp->GetEntries()>10){
	   fr[0]=15;
	   fr[1]=150;//fr[0]+0.5*temp->GetRMS();
	   TF1 *fitsnr_temp=langaufit(temp,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
	
	   double mpv=fitsnr_temp->GetParameter(1);
	   double empv=fitsnr_temp->GetParError(1);
	   double wmpv=fitsnr_temp->GetParameter(0);
	   double chi2ndf=0;
	   fout_mip<<layer<<" "<<i<<" "<<j<<" "<<mpv<<" "<<empv<<" "<<wmpv<<" "<<chi2ndf<<" "<<temp->GetEntries()<<"\n";

	   if(ndf>0) chi2ndf=chisqr/ndf;

	   MIPM->Fill(i,j,mpv);
	   MIPM_xy->Fill(map_pointX[i][j],map_pointY[i][j],mpv);
	   MIPRMS->Fill(map_pointX[i][j],map_pointY[i][j],fitsnr_temp->GetParameter(3));
	   MIPW->Fill(map_pointX[i][j],map_pointY[i][j],wmpv);

	   canvas_mip->cd(j+1);
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
	 } else {
           fout_mip<<layer<<" "<<i<<" "<<j<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<"\n";
	 }
      }
      file->cd();
      //canvas_mip->Print(TString::Format("plots/MIPs_%s_layer_%i_chip%i.eps",run.Data(),layer,i));
      canvas_mip->Write();
      delete canvas_mip;
  }


  //  _file0 ->Close();

  TCanvas* canvas2= new TCanvas(TString::Format("MIPAna_layer_%i",layer),TString::Format("MIPAna_layer_%i",layer),800,1600);   
  canvas2->Divide(1,2);
  
  canvas2->cd(1);
  MIPM->GetXaxis()->SetTitle("CHIP");
  MIPM->GetYaxis()->SetTitle("chn");
  MIPM->GetZaxis()->SetRangeUser(15,80);
  MIPM->Draw("colz");


  canvas2->cd(2);
  MIPM_xy->GetXaxis()->SetTitle("x");
  MIPM_xy->GetYaxis()->SetTitle("y");
  MIPM_xy->GetZaxis()->SetRangeUser(15,80);
  MIPM_xy->Draw("colz");

  canvas2->Print(TString::Format("plots/MIPAna_%s_layer_%i.eps",run.Data(),layer));

  TCanvas* canvas= new TCanvas(TString::Format("MIPAna2_layer_%i",layer),TString::Format("MIPAna2_layer_%i",layer),1600,1600);   
  canvas->Divide(2,2);
  
  canvas->cd(1);
    MIPM->GetXaxis()->SetTitle("x");
    MIPM->GetYaxis()->SetTitle("y");
    MIPM->GetZaxis()->SetRangeUser(20,80);
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

    // canvas->Print(TString::Format("plots/MIPAna_%s_layer_%i.eps",run.Data(),layer));

  file->cd();
  MIPM->Write();
  MIPM_xy->Write();
  MIPW->Write();
  MIPN->Write();
  MIPRMS->Write();
  canvas2->Write();
  canvas->Write();
  delete  canvas;
  delete  canvas2;
  file->Close();
  
}

