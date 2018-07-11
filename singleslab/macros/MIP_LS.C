#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "../../style/Style.C"
#include "../../style/Labels.C"

using namespace std;
// LANDAU STUFF
//---------------------------------------------------------------------------------------

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



int MIP_LS(){


  SetIrlesStyle();
  gStyle->SetOptFit(0111); 
  gStyle->SetOptStat(1110);

  gStyle->SetOptTitle(1);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.2);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleStyle(0);
  gStyle->SetMarkerSize(1.2);

  TString angle = "angle0";
  int channel[16]= {4,11,15,18,20,23,27,30,33,36,37,40,42,45,48,49};


  TGraphErrors *position_channel[16];
  TGraphErrors *width_channel[16];
  TGraphErrors *entries_channel[16];

  for(int ichn=0; ichn<16; ichn++) {

    double x[8], y[8], ey[8], y2[8], ey2[8], entries[8];
    
    for(int iasu=1; iasu<9; iasu++) {
      
      x[iasu-1]=iasu;
      y[iasu-1]=0;
      ey[iasu-1]=0;
      y2[iasu-1]=0;
      ey2[iasu-1]=0;
      entries[iasu-1]=0;

      //open file
      TString filename =TString::Format("../results_mipcalibration/Signal_summary_dif_1_1_1_ASU%i_%s.root",iasu,angle.Data());

      TFile *f = new TFile(filename);
      cout<<filename<<" "<<x[iasu-1]<<endl;
      
      int ichip=10+16*(iasu-1);

      TH1F *hmip_good  = (TH1F*)f->Get(TString::Format("histograms/charge_dif_1_1_1_chip%i_chn%i",ichip,channel[ichn]));
      if(hmip_good ) {
	// Setting fit range and start values
	Double_t fr[2];
	Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];   
	pllo[0]=0.; pllo[1]=0.0; pllo[2]=10.0; pllo[3]=0.;
	plhi[0]=100.0; plhi[1]=200.0; plhi[2]=100000000.0; plhi[3]=20.0;
	sv[0]=15.0;
	Double_t chisqr;
	Int_t    ndf;
	
	fr[0]=hmip_good->GetMean()-0.8*hmip_good->GetRMS();
	fr[1]=hmip_good->GetMean();
	sv[0]=hmip_good->GetRMS()*0.5;
	sv[1]=hmip_good->GetMean()*0.6;
	sv[2]=hmip_good->Integral("width");
	sv[3]=hmip_good->GetRMS()/5.;
	
	TF1 *fitsnr_temp=langaufit(hmip_good,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
	
	entries[iasu-1]=hmip_good->GetEntries();
	
	if(hmip_good->GetEntries()>50) {
	  y[iasu-1]=fitsnr_temp->GetParameter(1);
	  ey[iasu-1]=fitsnr_temp->GetParError(1);
	  y2[iasu-1]=fitsnr_temp->GetParameter(0);
	  ey2[iasu-1]=fitsnr_temp->GetParError(0);

	}
      }
      f->Close();

    }//iasu

    position_channel[ichn]=new TGraphErrors(8,x,y,0,ey);
    width_channel[ichn]=new TGraphErrors(8,x,y2,0,0);
    entries_channel[ichn]=new TGraphErrors(8,x,entries,0,0);


  }//ichn


  TLegend *leg1 = new TLegend(0.2,0.6,0.5,0.9);
  TLegend *leg2 = new TLegend(0.55,0.6,0.85,0.9);

  TCanvas *canvas = new TCanvas(angle,angle,1400,600);
  canvas->Divide(3,1);
  canvas->cd(1);
  
  for(int ichn=0;ichn<16;ichn++) {
    if(ichn<8) {
	position_channel[ichn]->SetMarkerColor(1+ichn);
	position_channel[ichn]->SetMarkerStyle(20+ichn);
	position_channel[ichn]->SetLineColor(1+ichn);
	position_channel[ichn]->SetLineStyle(1);
	position_channel[ichn]->SetLineWidth(2);
      } else {
	position_channel[ichn]->SetMarkerColor(1+ichn-8);
	position_channel[ichn]->SetMarkerStyle(20+ichn-8);
	position_channel[ichn]->SetLineColor(1+ichn-8);
	position_channel[ichn]->SetLineStyle(2);
	position_channel[ichn]->SetLineWidth(2);
      }
      if(ichn==0) {
	position_channel[ichn]->SetTitle(angle);
	position_channel[ichn]->GetXaxis()->SetTitle("ASU");
	position_channel[ichn]->GetYaxis()->SetTitle("MPV [ADC]");
	position_channel[ichn]->GetYaxis()->SetRangeUser(40,80);
	position_channel[ichn]->Draw("ap");
      } else position_channel[ichn]->Draw("p");
      
  }//ichn

  canvas->cd(2);
  for(int ichn=0;ichn<16;ichn++) {
    if(ichn<8) {
      width_channel[ichn]->SetMarkerColor(1+ichn);
      width_channel[ichn]->SetMarkerStyle(20+ichn);
      width_channel[ichn]->SetLineColor(1+ichn);
      width_channel[ichn]->SetLineStyle(1);
	width_channel[ichn]->SetLineWidth(2);
      } else {
	width_channel[ichn]->SetMarkerColor(1+ichn-8);
	width_channel[ichn]->SetMarkerStyle(20+ichn-8);
	width_channel[ichn]->SetLineColor(1+ichn-8);
	width_channel[ichn]->SetLineStyle(2);
	width_channel[ichn]->SetLineWidth(2);
      }
      if(ichn==0) {
	width_channel[ichn]->SetTitle(angle);
	width_channel[ichn]->GetXaxis()->SetTitle("ASU");
	width_channel[ichn]->GetYaxis()->SetTitle("width landau [ADC]");
	width_channel[ichn]->GetYaxis()->SetRangeUser(0,20);
	width_channel[ichn]->Draw("ap");
      } else width_channel[ichn]->Draw("p");
      if(ichn<8) leg1->AddEntry(width_channel[ichn],TString::Format("channel=%i",channel[ichn]),"lp");
      else leg2->AddEntry(width_channel[ichn],TString::Format("channel=%i",channel[ichn]),"lp");
      
  }//ichn

  leg1->Draw();
  leg2->Draw();
  
  canvas->cd(3);
  for(int ichn=0;ichn<16;ichn++) {
    if(ichn<8) {
      entries_channel[ichn]->SetMarkerColor(1+ichn);
      entries_channel[ichn]->SetMarkerStyle(20+ichn);
      entries_channel[ichn]->SetLineColor(1+ichn);
      entries_channel[ichn]->SetLineStyle(1);
	entries_channel[ichn]->SetLineWidth(2);
      } else {
	entries_channel[ichn]->SetMarkerColor(1+ichn-8);
	entries_channel[ichn]->SetMarkerStyle(20+ichn-8);
	entries_channel[ichn]->SetLineColor(1+ichn-8);
	entries_channel[ichn]->SetLineStyle(2);
	entries_channel[ichn]->SetLineWidth(2);
      }
      if(ichn==0) {
	entries_channel[ichn]->SetTitle(angle);
	entries_channel[ichn]->GetXaxis()->SetTitle("ASU");
	entries_channel[ichn]->GetYaxis()->SetTitle("Ntriggers");
	entries_channel[ichn]->GetYaxis()->SetRangeUser(0,20000);
	entries_channel[ichn]->Draw("ap");
      } else entries_channel[ichn]->Draw("p");
  }
  
    canvas->Print(TString::Format("mip_allasu_canvas_%s.png",angle.Data()));

  return 0;
 
}
