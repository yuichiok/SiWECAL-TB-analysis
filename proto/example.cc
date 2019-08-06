#define example_cxx
#include "example.h"
#include <TPaveStats.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <TFitResult.h>
#include <TSpectrum.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;


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



bool example::IsHit(int ihit) {

  if( (hit_isMasked[ihit]==0) && hit_isHit[ihit]==1) return true;
  return false;
}

bool example::IsPedestal(int ihit) {

  if( (hit_isMasked[ihit]==0) &&  hit_isHit[ihit]==0) return true;
  return false;
  
}

bool example::TrackBasicSelection( int nslabs_selection=7) {


  bool bool_hit_slab[10];
  for(int i=0; i<NSLABS; i++ ) bool_hit_slab[i] = false;
  int nhits_sl=0;
  int nhits_dif=0;

  for(int ihit=0; ihit< nhit_chan; ihit ++) {
    int z=hit_slab[ihit];

    if(IsHit(ihit)==true ) {
      bool_hit_slab[z]=true;
      if(z>4) nhits_sl++;
      else nhits_dif++;
    }
  }
  

  int nslabhitted_fev13 =0;
  int nslabhitted_sl =0;
  int nslabhitted =0;
  for(int i=0; i<NSLABS; i++ ){
    if( bool_hit_slab[i]==true)  nslabhitted++;
  }
  for(int i=0; i<5; i++ ) if( bool_hit_slab[i]==true) nslabhitted_fev13++;
  for(int i=5; i<9; i++ ) if( bool_hit_slab[i]==true) nslabhitted_sl++;

  if( (nslabhitted_fev13<4 || nslabhitted_sl<4)) return false;// || nhits_dif>5 || nhits_sl>4 ) return false;
  return true;

}

bool example::TrackTightSelection( int nslabs_selection=9, int nhits_perslab=3) {


  bool bool_hit_slab[10];
  for(int i=0; i<NSLABS; i++ ) bool_hit_slab[i] = false;

  int n_hit_slab[10];
  for(int i=0; i<NSLABS; i++ ) n_hit_slab[i] = 0;


  for(int ihit=0; ihit< nhit_chan; ihit ++) {
    int z=hit_slab[ihit];
    if(IsHit(ihit)==true ) {
      bool_hit_slab[z]=true;
      n_hit_slab[z]++;
    }
  }
  
  int nslabhitted =0;
  for(int i=0; i<NSLABS; i++ ){
    if( bool_hit_slab[i]==true)  nslabhitted++;
  }
  
 
  if( nslabhitted<9) return false;
  for(int i=0; i<NSLABS; i++ ) if( n_hit_slab[i]>nhits_perslab) return false;
  
  return true;

}



void example::SimpleEvDisplayTrack(TString outputname="")
{

   // --------------
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //-----------------

  //Objects to store the chip/channel results:
  // numbers 
  TH3F* mip_evdisp[500];
  for(int i=0; i<500; i++) mip_evdisp[i]= new TH3F(TString::Format("mip_evdisp_%i",i),TString::Format("mip_evdisp_%i",i),32,-90,90,20,-0.5,19.5,32,-90,90);
 


  int ntracks=0;
  int ntracks_50=0;
  // -----------------------------------------------------------------------------------------------------   
  // Signal readout
  Long64_t nbytes = 0, nb = 0;
  //  nentries=500000;
  for (Long64_t jentry=0; jentry<nentries;jentry++){//nentries-100;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(ntracks>499) break;

    int fivepercent=nentries/20;
    if( (jentry % fivepercent) == 0 ) cout<< "Progress: "<<100.*jentry/nentries<<" %" <<endl;

    if( bcid<50 || (bcid>899 && bcid<930)) continue;
    if(TrackBasicSelection()==false) continue;
    // std::cout<<ntracks<<" "<<jentry<<endl;
    // ntracks_50++;
    // if(ntracks_50>) {

    bool tracks=false;
    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      if(  hit_isMasked[ihit] == 0 && hit_isHit[ihit]==1 ) {
	tracks=true;
	mip_evdisp[ntracks]->Fill(hit_x[ihit],hit_z[ihit],hit_y[ihit],hit_energy[ihit]);
	//	else mip_evdisp[ntracks]->Fill(-hit_x[ihit],hit_z[ihit],-hit_y[ihit],hit_lg[ihit]);
	
      }
    }
    if(tracks==true) ntracks++;
    //   ntracks_50=0;
      //  }

    
  }

  // -----------------------------------------------------------------------------------------------------   
  // Signal analysis
  TFile *signalfile_summary = new TFile("results_mipcalibration/EvDisplay"+outputname+".root" , "RECREATE");
  signalfile_summary->cd();
  for(int i=0; i<500; i++) {
    mip_evdisp[i]->GetXaxis()->SetTitle("X");
    mip_evdisp[i]->GetYaxis()->SetTitle("Z");
    mip_evdisp[i]->GetZaxis()->SetTitle("Y");
    mip_evdisp[i]->Write();
  }

  signalfile_summary->Close();


}



void example::SimpleDistributionsTrack(TString outputname="")
{

   // --------------
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //-----------------

  //Objects to store the chip/channel results:
  // numbers 
  TH1F* mip_histo_all = new TH1F("mip_histo_all","mip_histo_all",4096,0.1,4096.1);
  TH2F* delta_x = new TH2F("delta_x","delta_x",9,-0.5,8.5,33,-92.8125,92.8125);
  TH2F* delta_y = new TH2F("delta_y","delta_y",9,-0.5,8.5,33,-92.8125,92.8125);


  int ntracks=0;
 
  // -----------------------------------------------------------------------------------------------------   
  // Signal readout
  Long64_t nbytes = 0, nb = 0;
  //  nentries=500000;
  for (Long64_t jentry=0; jentry<nentries;jentry++){//nentries-100;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    int fivepercent=nentries/20;
    if( (jentry % fivepercent) == 0 ) cout<< "Progress: "<<100.*jentry/nentries<<" %" <<endl;

    if( bcid<50 || (bcid>899 && bcid<950)) continue;
    if(TrackTightSelection()==false) continue;

    double hit_x_slab5[1000];
    double hit_y_slab5[1000];
    for(int i=0; i<1000; i++) {
      hit_x_slab5[i]=-10000;
      hit_y_slab5[i]=-10000;
    }
    int hit_slab5=0;
    
    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      if(  hit_isMasked[ihit] == 0 && hit_isHit[ihit]==1 ) {
	int z=hit_slab[ihit];
	mip_histo_all->Fill(hit_energy[ihit]);
	if(z==5) {
	  hit_x_slab5[hit_slab5]=hit_x[ihit];
	  hit_y_slab5[hit_slab5]=hit_y[ihit];
	  hit_slab5++;
	}
      }
    }

    for(int i=0; i<hit_slab5; i++) {
      for(int ihit=0; ihit< nhit_chan; ihit ++) {
	int z=hit_slab[ihit];
	if(  hit_isHit[ihit]==1 && z!=5) {
	  delta_x->Fill(z,hit_x_slab5[i]-hit_x[ihit]);
	  delta_y->Fill(z,hit_y_slab5[i]-hit_y[ihit]);
	}
      }
    }
    
  }

  // -----------------------------------------------------------------------------------------------------   
  // Signal analysis
  TFile *signalfile_summary = new TFile("results_mipcalibration/Signal_summary"+outputname+".root" , "RECREATE");
  signalfile_summary->cd();

  mip_histo_all->GetXaxis()->SetTitle("Energy [MIP]");
  mip_histo_all->GetYaxis()->SetTitle("# entries");
  mip_histo_all->SetTitle("Single cell hit energy in 3GeV e^{+} beam");
  mip_histo_all->SetName("energy_distribution");
  mip_histo_all->Write();

  delta_x->Write();
  delta_y->Write();

  /*
  // Setting fit range and start values
  Double_t fr[2];
  Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];   
  pllo[0]=0.001; pllo[1]=50; pllo[2]=100.0; pllo[3]=0.001;
  plhi[0]=30.; plhi[1]=500.0; plhi[2]=100000000.0; plhi[3]=30.0;

  sv[0]=15.0;
  Double_t chisqr;
  Int_t    ndf;

  // -----------------------------------------------------------------------------------
  // iteration 1
  // ***********************************************************************************

  //-----------------------------------
  // Peak 1
  //------------------------------------
  TH1F *hist_peak1_iter1 = (TH1F*)mip_histo_all->Clone();
  TSpectrum *s_iter1 = new TSpectrum(2);
  int npeaks_iter1_0 = s_iter1->Search(hist_peak1_iter1,2,"",0.4); 
  Double_t *mean_peak_iter1_0=s_iter1->GetPositionX();
  fr[0]=mean_peak_iter1_0[1]-0.15*hist_peak1_iter1->GetRMS();
  fr[1]=mean_peak_iter1_0[1]+0.4*hist_peak1_iter1->GetRMS();
  sv[0]=5;//hist_peak1_iter1->GetRMS()/10;
  sv[1]=90;//hist_peak1_iter1->GetMean();
  sv[2]=100000;//hist_peak1_iter1->Integral("width");
  sv[3]=5;//hist_peak1_iter1->GetRMS()/10.;
  
  TF1 *fit_1st_iter_temp=langaufit(hist_peak1_iter1,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  
  hist_peak1_iter1->GetXaxis()->SetTitle("Energy [MIP]");
  hist_peak1_iter1->GetYaxis()->SetTitle("# entries");
  hist_peak1_iter1->SetTitle("First MIP");
  hist_peak1_iter1->SetName("first_peak_fit_histo");
  hist_peak1_iter1->Write();

  double mpv=fit_1st_iter_temp->GetParameter(1);
  double empv=fit_1st_iter_temp->GetParError(1);
  double wmpv=fit_1st_iter_temp->GetParameter(0);
  double wg=fit_1st_iter_temp->GetParameter(3);
  double chi2ndf=0;
  if(ndf>0) chi2ndf=chisqr/ndf;
  double mipentries=mip_histo_all->GetEntries();

  TF1 *fpeak1_iter1 = new TF1("fpeak1_iter1",langaufun,0,5,4);
  fpeak1_iter1->SetParameter(0,fit_1st_iter_temp->GetParameter(0));
  fpeak1_iter1->SetParameter(1,fit_1st_iter_temp->GetParameter(1));
  fpeak1_iter1->SetParameter(2,fit_1st_iter_temp->GetParameter(2));
  fpeak1_iter1->SetParameter(3,fit_1st_iter_temp->GetParameter(3));
 
  cout<<npeaks_iter1_0<<" "<<mean_peak_iter1_0[0]<<endl;
  cout<<" "<<mpv<<" "<<empv<<" "<<wmpv<<" "<<wg<<" "<<chi2ndf<<" "<<mipentries<<endl;

  
  // ------------------------------------------------------
  // Peak 2
  // -----------------------------------------------------
  TH1F *hist_peak2_iter1 = (TH1F*)mip_histo_all->Clone();
  hist_peak2_iter1->Add(fpeak1_iter1,-1);

  TSpectrum *s_iter1_1 = new TSpectrum(1);
  int npeaks_iter1_1 = s_iter1_1->Search(hist_peak2_iter1,0.5,"",0.4); 
  Double_t *mean_peak_iter1_1=s_iter1_1->GetPositionX();
  fr[0]=mean_peak_iter1_1[0]-0.2*mip_histo_all->GetRMS();//-0.65*mip_histo_all->GetRMS();
  fr[1]=mean_peak_iter1_1[0]+0.3*mip_histo_all->GetRMS();//mip_histo_all->GetMean();
  sv[0]=wmpv;
  sv[1]=2*mpv;
  sv[2]=hist_peak2_iter1->Integral("width");
  sv[3]=hist_peak2_iter1->GetRMS()/10.;

  pllo[1]=2*mpv - 15*wmpv;
  plhi[1]=2*mpv + 15*wmpv;

  TF1 *fit_1st_iter_temp2=langaufit(hist_peak2_iter1,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  mpv=fit_1st_iter_temp2->GetParameter(1);
  empv=fit_1st_iter_temp2->GetParError(1);
  wmpv=fit_1st_iter_temp2->GetParameter(0);
  wg=fit_1st_iter_temp2->GetParameter(3);
  chi2ndf=0;
  if(ndf>0) chi2ndf=chisqr/ndf;

  cout<<npeaks_iter1_1<<" "<<mean_peak_iter1_1[0]<<endl;

  hist_peak2_iter1->GetXaxis()->SetTitle("Energy [MIP]");
  hist_peak2_iter1->GetYaxis()->SetTitle("# entries");
  hist_peak2_iter1->SetTitle("Second MIP (first subtracted)");
  hist_peak2_iter1->SetName("second_peak_fit_histo");
  hist_peak2_iter1->Write();

  TF1 *fpeak2_iter1 = new TF1("fpeak2_iter1",langaufun,0,5,4);
  fpeak2_iter1->SetParameter(0,fit_1st_iter_temp2->GetParameter(0));
  fpeak2_iter1->SetParameter(1,fit_1st_iter_temp2->GetParameter(1));
  fpeak2_iter1->SetParameter(2,fit_1st_iter_temp2->GetParameter(2));
  fpeak2_iter1->SetParameter(3,fit_1st_iter_temp2->GetParameter(3));


  // ------------------------------------------------------
  // Peak 2
  // -----------------------------------------------------
  TH1F *hist_peak3_iter1 = (TH1F*)mip_histo_all->Clone();
  hist_peak3_iter1->Add(fpeak1_iter1,-1);
  hist_peak3_iter1->Add(fpeak2_iter1,-1);

  TSpectrum *s_iter1_2 = new TSpectrum(1);
  int npeaks_iter1_2 = s_iter1_2->Search(hist_peak3_iter1,0.5,"",0.8); 
  Double_t *mean_peak_iter1_2=s_iter1_2->GetPositionX();
  fr[0]=mean_peak_iter1_2[0]-0.4*mip_histo_all->GetRMS();//-0.65*mip_histo_all->GetRMS();
  fr[1]=mean_peak_iter1_2[0]+0.4*mip_histo_all->GetRMS();//mip_histo_all->GetMean();
  sv[0]=wmpv;
  sv[1]=3/2*mpv;
  sv[2]=hist_peak3_iter1->Integral("width");
  sv[3]=hist_peak3_iter1->GetRMS()/10.;

  pllo[1]=3/2*mpv - 5*wmpv;
  plhi[1]=3/2*mpv + 5*wmpv;

  TF1 *fit_1st_iter_temp3=langaufit(hist_peak3_iter1,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  mpv=fit_1st_iter_temp3->GetParameter(1);
  empv=fit_1st_iter_temp3->GetParError(1);
  wmpv=fit_1st_iter_temp3->GetParameter(0);
  wg=fit_1st_iter_temp3->GetParameter(3);
  chi2ndf=0;
  if(ndf>0) chi2ndf=chisqr/ndf;

  cout<<npeaks_iter1_2<<" "<<mean_peak_iter1_2[0]<<endl;

  hist_peak3_iter1->GetXaxis()->SetTitle("Energy [MIP]");
  hist_peak3_iter1->GetYaxis()->SetTitle("# entries");
  hist_peak3_iter1->SetTitle("Third MIP (first+second subtracted)");
  hist_peak3_iter1->SetName("third_peak_fit_histo");
  hist_peak3_iter1->Write();

  TF1 *fpeak3_iter1 = new TF1("fpeak3_iter1",langaufun,0,5,4);
  fpeak3_iter1->SetParameter(0,fit_1st_iter_temp3->GetParameter(0));
  fpeak3_iter1->SetParameter(1,fit_1st_iter_temp3->GetParameter(1));
  fpeak3_iter1->SetParameter(2,fit_1st_iter_temp3->GetParameter(2));
  fpeak3_iter1->SetParameter(3,fit_1st_iter_temp3->GetParameter(3));

  // ------------------------------------------------------
  //write TF1s
  fpeak1_iter1->SetName("FirstPeak_iter1");
  fpeak1_iter1->Write();

  fpeak2_iter1->SetName("SecondPeak_iter1");
  fpeak2_iter1->Write();

  fpeak3_iter1->SetName("ThirdPeak_iter1");
  fpeak3_iter1->Write();


  TH1F* first_iteration_peak1 = new TH1F("first_iteration_peak1","first_iteration_peak1",925,0.005,1.855);
  TH1F* first_iteration_peak2 = new TH1F("first_iteration_peak2","first_iteration_peak2",500,1.855,2.855);
  TH1F* first_iteration_peak3 = new TH1F("first_iteration_peak3","first_iteration_peak3",1500,2.855,5.855);
  first_iteration_peak1->Add(fpeak1_iter1);
  first_iteration_peak1->Add(fpeak2_iter1);
  first_iteration_peak1->Add(fpeak3_iter1);

  first_iteration_peak2->Add(fpeak1_iter1);
  first_iteration_peak2->Add(fpeak2_iter1);
  first_iteration_peak2->Add(fpeak3_iter1);

  first_iteration_peak3->Add(fpeak1_iter1);
  first_iteration_peak3->Add(fpeak2_iter1);
  first_iteration_peak3->Add(fpeak3_iter1);

  first_iteration_peak1->Write();
  first_iteration_peak2->Write();
  first_iteration_peak3->Write();*/

  signalfile_summary->cd();

  signalfile_summary->Close();


}



// LANDAU STUFF
//---------------------------------------------------------------------------------------




TF1* example::langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
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

  his->Fit(FunName,"RBQMS","same");   // fit within specified range, use ParLimits, quiet, improve fit results (TMINUIT)

  ffit->GetParameters(fitparams);    // obtain fit parameters
  for (i=0; i<4; i++) {
    fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
  }
  ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
  NDF[0] = ffit->GetNDF();           // obtain ndf

  return (ffit);              // return fit function

}


