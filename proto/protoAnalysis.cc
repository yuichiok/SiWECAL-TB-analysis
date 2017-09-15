#define protoAnalysis_cxx
#include "protoAnalysis.h"
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


void protoAnalysis::SimpleDistributionsShower(TString outputname="")
{

  TH1F* SCA_shower = new TH1F("SCA_shower","SCA_shower",15,-0.5,14.5);
  TH1F* BCID_shower = new TH1F("BCID_shower","BCID_shower",60,25,3025);
  TH2F* SCA_BCID_shower = new TH2F("SCA_BCID_shower","SCA_BCID_shower",15,-0.5,14.5,60,25,3025);
  TH2F* BCID_PREV_shower = new TH2F("BCID_PREV_shower","BCID_PREV_shower",1200,2.5,6002.5,200,2.5,1002.5);

  // --------------
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //-----------------

  double mipcut=0.5;
 // Signal readout
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(bcid<1300) continue;
    if(bcid>2900) continue;

    if(nhit_slab<7) continue;


    bool bool_hit_slab[7];
    bool_hit_slab[0] = false;
    bool_hit_slab[1] = false;
    bool_hit_slab[2] = false;
    bool_hit_slab[3] = false;
    bool_hit_slab[4] = false;
    bool_hit_slab[5] = false;
    bool_hit_slab[6] = false;

    int n_hitstotal =0;  

    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      int z=hit_z[ihit];
      if(z==9) z=6;
      if(hit_energy[ihit]>mipcut && hit_isMasked[ihit]==0 ) {
	bool_hit_slab[z]=true;
      }
    }
 

    int nslabhitted =0;

    for(int i=0; i<7; i++ ){
      if( bool_hit_slab[i]==true) nslabhitted++;
    }

    if( nslabhitted<7) continue;

    bool fill=false;
    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      
      if(hit_energy[ihit]>mipcut && hit_isMasked[ihit]==0 ) {
    	double xi=hit_x[ihit];
    	double yi=hit_y[ihit];
    	double zi=hit_z[ihit];
	
    	double mindist=1000000000.;
    	double dist=-10;
    	for(int jhit=0; jhit< nhit_chan; jhit ++) {
    	  if(hit_energy[jhit]>mipcut && hit_isMasked[jhit]==0 && ihit!=jhit) {
    	    double xj=hit_x[jhit];
    	    double yj=hit_y[jhit];
    	    double zj=hit_z[jhit];
    	    dist= sqrt( (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + 150*150*(zi-zj)*(zi-zj) );
	    if(dist < mindist) mindist = dist;
    	  }
    	}

	if(mindist<7.) {
	  fill=true;
	  SCA_shower->Fill(hit_sca[ihit]);
	  SCA_BCID_shower->Fill(hit_sca[ihit],bcid);
	}
      }

    }

    if(fill==true) {
    
      BCID_PREV_shower->Fill(bcid,bcid-prev_bcid);
      BCID_shower->Fill(bcid);
    }

   

  }
  // -----------------------------------------------------------------------------------------------------   
  // Signal analysis
  TFile *signalfile_summary = new TFile("results_showers/SimpleDist"+outputname+".root" , "RECREATE");
  signalfile_summary->cd();

  SCA_shower->GetXaxis()->SetTitle("SCA");
  SCA_shower->GetYaxis()->SetTitle("# entries");
  SCA_shower->SetTitle("SCA distribution for shower-like events");
  SCA_shower->SetName("SCA_shower");
  SCA_shower->Write();

  BCID_shower->GetXaxis()->SetTitle("bcid");
  BCID_shower->GetYaxis()->SetTitle("# entries");
  BCID_shower->SetTitle("bcid distribution for shower-like events");
  BCID_shower->SetName("BCID_shower");
  BCID_shower->Write();

  SCA_BCID_shower->GetXaxis()->SetTitle("SCA");
  SCA_BCID_shower->GetYaxis()->SetTitle("bcid");
  SCA_BCID_shower->SetTitle("SCA vs bcid distribution for shower-like events");
  SCA_BCID_shower->SetName("SCA_BCID_shower");
  SCA_BCID_shower->Write();
  
  BCID_PREV_shower->GetXaxis()->SetTitle("bcid");
  BCID_PREV_shower->GetYaxis()->SetTitle("bcid-prev_bcid");
  BCID_PREV_shower->SetTitle("bcid vs bcid-prev_bcid distribution for shower-like events");
  BCID_PREV_shower->SetName("BCID_PREV_shower");
  BCID_PREV_shower->Write();

  signalfile_summary->Close();

}


void protoAnalysis::ShowerDistributions(TString folder="", TString configuration="conf1", TString energy_string="3GeV", TString gridpoint ="grid20", double mipcut=0.5)
{

  // --------------
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //-----------------

  Float_t w[10];
  Float_t thickness[7];
  thickness[0]=2.;
  thickness[1]=2.;
  thickness[2]=2.;
  thickness[3]=2.;
  thickness[4]=4.;
  thickness[5]=4.;
  thickness[6]=6.;

  if(configuration=="conf2") {
    thickness[0]=4.;
    thickness[1]=2.;
    thickness[2]=2.;
    thickness[3]=4.;
    thickness[4]=4.;
    thickness[5]=6.;
    thickness[6]=6.;
  }

  if(configuration=="conf3") {
    thickness[0]=6.;
    thickness[1]=2.;
    thickness[2]=4.;
    thickness[3]=4.;
    thickness[4]=6.;
    thickness[5]=6.;
    thickness[6]=6.;
  }

  w[0]=0.56*thickness[0];
  w[1]=0.56*thickness[1];
  w[2]=0.56*thickness[2];
  w[3]=0.56*thickness[3];
  w[4]=0.56*thickness[4];
  w[5]=0.56*thickness[5];
  w[6]=0.56*0;
  w[7]=0.56*0;
  w[8]=0.56*0;
  w[9]=0.56*thickness[6];

  Float_t edge =0.0001;

  Float_t bins[] = { 0, w[0]+edge, w[0]+edge+w[1], w[0]+edge+w[1]+w[2], w[0]+edge+w[1]+w[2]+w[3], w[0]+edge+w[1]+w[2]+w[3]+w[4],  w[0]+edge+w[1]+w[2]+w[3]+w[4]+w[5],  w[0]+edge+w[1]+w[2]+w[3]+w[4]+w[5]+w[9] };
  Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1; 

  Float_t binsx[33];
  binsx[0]=-86.;
  for(int ibins=1;ibins<33;ibins++) binsx[ibins]=binsx[ibins-1]+5.5;
  Int_t  binnumx = 32; 

  Float_t binsz[11];
  binsz[0]=-0.5;
  for(int ibins=1;ibins<11;ibins++) binsz[ibins]=binsz[ibins-1]+1;
  Int_t  binnumz = 10; 

  // energy histograms
  //TH1 histograms
  TH1F* energy = new TH1F("energy","energy",500,1,1001);
  TH1F* energy_center = new TH1F("energy_center","energy_center",500,1,1001);
  TH1F* energy_profile_z = new TH1F("energy_profile_z","energy_profile_z",7,-0.5,6.5);//binnumz,binsz);
  energy_profile_z->Sumw2();

  // x,y,z histograms
  TH2F* energy_xy = new TH2F("energy_xy","energy_xy",binnumx,binsx,binnumx,binsx);
  TH2F* energy_xz = new TH2F("energy_xz","energy_xz",binnumx,binsx,binnumz,binsz);
  TH2F* energy_yz = new TH2F("energy_yz","energy_yz",binnumx,binsx,binnumz,binsz);


  /// Xo thickness 2d plots

  TH1F* energy_profile_X0 = new TH1F("energy_profile_X0","energy_profile_X0",binnum,bins);
  TH2F* energy_xX0 = new TH2F("energy_xX0","energy_xX0",binnumx,binsx,binnum,bins);
  TH2F* energy_yX0 = new TH2F("energy_yX0","energy_yX0",binnumx,binsx,binnum,bins);

  TH1F *xbarycenter = new TH1F("x-barycenter","x-barycenter",binnumx,binsx);
  TH1F *ybarycenter = new TH1F("y-barycenter","y-barycenter",binnumx,binsx);
  TH1F *zbarycenter = new TH1F("z-barycenter","z-barycenter",binnumz,binsz);

  TH2F* xybarycenter = new TH2F("Energy xy-barycenter","Energy xy-barycenter",binnumx,binsx,binnumx,binsx);
  TH2F* xzbarycenter = new TH2F("Energy xz-barycenter","Energy xz-barycenter",binnumx,binsx,binnumz,binsz);
  TH2F* yzbarycenter = new TH2F("Energy yz-barycenter","Energy yz-barycenter",binnumx,binsx,binnumz,binsz);

  TH1F* n_x= new TH1F("n_x","n_x",binnumx,binsx);
  TH1F* n_y= new TH1F("n_y","n_y",binnumx,binsx);
  TH1F* n_z= new TH1F("n_z","n_z",binnumz,binsz);

  TH1F* nhitstotal= new TH1F("nhitstotal","nhitstotal",501,-0.5,500.5);
  TH2F* nhits_vs_energy= new TH2F("nhits_vs_energy","nhits_vs_energy",51,-2.5,252.5,101,-5,1005);
  TH2F* min_dist_vs_energy= new TH2F("min_dist_vs_energy","min_dist_vs_energy",250,0,250,100,-5,995);

  int number_hits=0;
  // -----------------------------------------------------------------------------------------------------   
  // Signal readout
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(bcid<1300) continue;
    if(nhit_slab!=7) continue;


    bool bool_hit_slab[7];
    bool_hit_slab[0] = false;
    bool_hit_slab[1] = false;
    bool_hit_slab[2] = false;
    bool_hit_slab[3] = false;
    bool_hit_slab[4] = false;
    bool_hit_slab[5] = false;
    bool_hit_slab[6] = false;

    int n_hitstotal =0;  

    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      int z=hit_z[ihit];
      if(z==9) z=6;
      if(hit_energy[ihit]>mipcut && hit_isMasked[ihit]==0 ) {
	bool_hit_slab[z]=true;
	n_hitstotal++;
      }
    }

    nhitstotal->Fill(n_hitstotal);
    

    int nslabhitted =0;

    for(int i=0; i<7; i++ ){
      if( bool_hit_slab[i]==true) nslabhitted++;
    }

    if( nslabhitted<7) continue;

    // ---------------------------------------------
    //   the energy and other variables,
    double energy_sum_tmp=0;
    double weight_X0_sum_tmp=0;
    double energy_X0_sum_tmp=0;
    double xm=0, ym=0., zm=0.;

    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      
      if(hit_energy[ihit]>mipcut && hit_isMasked[ihit]==0 ) {
    	  energy_sum_tmp += hit_energy[ihit];
    	  energy_X0_sum_tmp += w[int(hit_z[ihit])] * hit_energy[ihit];
	  weight_X0_sum_tmp += w[int(hit_z[ihit])];

	  xm +=  - hit_x[ihit] * hit_energy[ihit] * w[int(hit_z[ihit])];
	  ym +=  - hit_y[ihit] * hit_energy[ihit] * w[int(hit_z[ihit])];
	  zm +=  hit_z[ihit] * hit_energy[ihit] * w[int(hit_z[ihit])];
      }
    }

    //-------------------------------------------
    //calcualte min distance with other hits 
    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      if(hit_energy[ihit]>mipcut && hit_isMasked[ihit]==0 ) {
	double xi=hit_x[ihit];
	double yi=hit_y[ihit];
	double zi=hit_z[ihit];
	
	double mindist=1000000000.;
	double dist=-10;
	for(int jhit=0; jhit< nhit_chan; jhit ++) {
	  if(hit_energy[jhit]>mipcut && hit_isMasked[jhit]==0 && ihit!=jhit) {
	    double xj=hit_x[jhit];
	    double yj=hit_y[jhit];
	    double zj=hit_z[jhit];
	    dist= sqrt( (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + 150*150*(zi-zj)*(zi-zj) );
	    if(dist < mindist) mindist = dist;
	  }
	}
	min_dist_vs_energy->Fill(mindist,energy_X0_sum_tmp);
      }
    }
    // ---------------------------------------------
    //  recalculate the energy and other variables, with a cut in the distance (avoid isolated cells)
    energy_sum_tmp=0;
    weight_X0_sum_tmp=0;
    energy_X0_sum_tmp=0;
    xm=0, ym=0., zm=0.;

    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      
      if(hit_energy[ihit]>mipcut && hit_isMasked[ihit]==0 ) {
    	double xi=hit_x[ihit];
    	double yi=hit_y[ihit];
    	double zi=hit_z[ihit];
	
    	double mindist=1000000000.;
    	double dist=-10;
    	for(int jhit=0; jhit< nhit_chan; jhit ++) {
    	  if(hit_energy[jhit]>mipcut && hit_isMasked[jhit]==0 && ihit!=jhit) {
    	    double xj=hit_x[jhit];
    	    double yj=hit_y[jhit];
    	    double zj=hit_z[jhit];
    	    dist= sqrt( (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + 150*150*(zi-zj)*(zi-zj) );
    	    if(dist < mindist) mindist = dist;
    	  }
    	}
	if(mindist<7.) {
    	  energy_sum_tmp += hit_energy[ihit];
    	  energy_X0_sum_tmp += w[int(hit_z[ihit])] * hit_energy[ihit];
	  weight_X0_sum_tmp += w[int(hit_z[ihit])];

    	  xm +=  - hit_x[ihit] * hit_energy[ihit] * w[int(hit_z[ihit])];
    	  ym +=  - hit_y[ihit] * hit_energy[ihit] * w[int(hit_z[ihit])];
    	  zm +=  hit_z[ihit] * hit_energy[ihit] * w[int(hit_z[ihit])];
	}
      }
    }

    if(energy_sum_tmp<0.5) continue;
    number_hits++;

    nhits_vs_energy->Fill(n_hitstotal,energy_X0_sum_tmp);


    xm = xm / energy_X0_sum_tmp;
    ym = ym / energy_X0_sum_tmp;
    zm = zm / energy_X0_sum_tmp;

    for(int ihit=0; ihit< nhit_chan; ihit ++) {

      if(hit_energy[ihit]>mipcut && hit_isMasked[ihit]==0 ) {
	double zlayer=hit_z[ihit];
	if(hit_z[ihit]>8) zlayer=6;
	energy_profile_z->Fill(zlayer,hit_energy[ihit]/energy_sum_tmp);
	energy_profile_X0->Fill(hit_x0[ihit],hit_energy[ihit]/energy_sum_tmp);
	
	energy_xy->Fill(-hit_x[ihit],-hit_y[ihit],hit_energy[ihit]/energy_sum_tmp );
	energy_xz->Fill(-hit_x[ihit],hit_z[ihit],hit_energy[ihit]/energy_sum_tmp );
	energy_yz->Fill(-hit_y[ihit],hit_z[ihit],hit_energy[ihit]/energy_sum_tmp );     
	
	energy_xX0->Fill(-hit_x[ihit],hit_x0[ihit],hit_energy[ihit]/energy_sum_tmp);
	energy_yX0->Fill(-hit_y[ihit],hit_x0[ihit],hit_energy[ihit]/energy_sum_tmp);
      }

    }
    
    energy->Fill(energy_X0_sum_tmp);

    xbarycenter->Fill(xm, energy_X0_sum_tmp);
    ybarycenter->Fill(ym, energy_X0_sum_tmp);
    zbarycenter->Fill(zm,energy_X0_sum_tmp);

    xybarycenter->Fill(xm,ym,energy_X0_sum_tmp);
    xzbarycenter->Fill(xm,zm,energy_X0_sum_tmp);
    yzbarycenter->Fill(ym,zm,energy_X0_sum_tmp);

    n_x->Fill(xm);
    n_y->Fill(ym);
    n_z->Fill(zm);


  }


  // First analysis, normalization and barycenter graph filling
  double rawE_x[binnumx];
  double rawE_y[binnumx];
  double rawE_z[binnumz];

  double erawE_x[binnumx];
  double erawE_y[binnumx];
  double erawE_z[binnumz];

  double x[binnumx];
  double y[binnumx];
  double z[binnumz];

  for(int i=0; i<binnumx; i++) {
    rawE_x[i]=0.; rawE_y[i]=0.;
    x[i]=n_x->GetBinCenter(i+1);
    y[i]=n_y->GetBinCenter(i+1);

    if(  n_x->GetBinContent(i+1) >50 ) rawE_x[i]=xbarycenter->GetBinContent(i+1)/n_x->GetBinContent(i+1);
    if(  n_y->GetBinContent(i+1) >50 ) rawE_y[i]=ybarycenter->GetBinContent(i+1)/n_y->GetBinContent(i+1);


    if(  n_x->GetBinContent(i+1) >50 ) erawE_x[i]=rawE_x[i]/sqrt(n_x->GetBinContent(i+1));
    if(  n_y->GetBinContent(i+1) >50 ) erawE_y[i]=rawE_y[i]/sqrt(n_y->GetBinContent(i+1));

  }

  for(int i=0; i<binnumz; i++) {
    rawE_z[i]=0.;
    z[i]=n_z->GetBinCenter(i+1);

    if(  n_z->GetBinContent(i+1) >50 ) rawE_z[i]=zbarycenter->GetBinContent(i+1)/n_z->GetBinContent(i+1);
    if(  n_y->GetBinContent(i+1) >50 ) erawE_z[i]=rawE_z[i]/sqrt(n_z->GetBinContent(i+1));


  }

  TGraphErrors * E_xbarycenter = new TGraphErrors(32,x,rawE_x,0,erawE_x);
  TGraphErrors * E_ybarycenter = new TGraphErrors(32,y,rawE_y,0,erawE_y);
  TGraphErrors * E_zbarycenter = new TGraphErrors(10,z,rawE_z,0,erawE_z);

  energy_profile_z->Scale(1./number_hits);
  energy_profile_X0->Scale(1./number_hits);

  energy_xy->Scale(1./number_hits);
  energy_xz->Scale(1./number_hits);
  energy_yz->Scale(1./number_hits);

  energy_xX0->Scale(1./number_hits);
  energy_yX0->Scale(1./number_hits);

  xbarycenter->Scale(1./number_hits);
  ybarycenter->Scale(1./number_hits);
  zbarycenter->Scale(1./number_hits);
  
  xybarycenter->Scale(1./number_hits);
  xzbarycenter->Scale(1./number_hits);
  yzbarycenter->Scale(1./number_hits);

  E_xbarycenter->SetName("mean_Eraw_vs_xm");
  E_ybarycenter->SetName("mean_Eraw_vs_ym");
  E_zbarycenter->SetName("mean_Eraw_vs_zm");



  // #########################################################################################################
  // -----------------------------------------------------
  // SECOND ANALYSIS, using barycenter of the shower as input for selection
  TF1 *fit_xbarycenter = new TF1("fit_xbarycenter","gaus");
  xbarycenter->Fit(fit_xbarycenter,"M");
  
  double xm_max=fit_xbarycenter->GetParameter(1)+2.*fit_xbarycenter->GetParameter(2);
  double xm_min=fit_xbarycenter->GetParameter(1)-2.*fit_xbarycenter->GetParameter(2);

   TF1 *fit_ybarycenter = new TF1("fit_ybarycenter","gaus");
  ybarycenter->Fit(fit_ybarycenter,"M");
  
  double ym_max=fit_ybarycenter->GetParameter(1)+2.*fit_ybarycenter->GetParameter(2);
  double ym_min=fit_ybarycenter->GetParameter(1)-2.*fit_ybarycenter->GetParameter(2);

  TF1 *fit_zbarycenter = new TF1("fit_zbarycenter","gaus");
  zbarycenter->Fit(fit_zbarycenter,"M");
  
  double zm_max=fit_zbarycenter->GetParameter(1)+1.0*fit_zbarycenter->GetParameter(2);
  double zm_min=fit_zbarycenter->GetParameter(1)-1.0*fit_zbarycenter->GetParameter(2);

  // -----------------------------------------------------------------------------------------------------   

  // Signal readout
  nbytes = 0; nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(bcid<1300) continue;
    if(nhit_slab!=7) continue;


    bool bool_hit_slab[7];
    bool_hit_slab[0] = false;
    bool_hit_slab[1] = false;
    bool_hit_slab[2] = false;
    bool_hit_slab[3] = false;
    bool_hit_slab[4] = false;
    bool_hit_slab[5] = false;
    bool_hit_slab[6] = false;

    int n_hitstotal =0;  

    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      int z=hit_z[ihit];
      if(z==9) z=6;
      if(hit_energy[ihit]>mipcut && hit_isMasked[ihit]==0 ) {
	bool_hit_slab[z]=true;
	n_hitstotal++;
      }
    }

  

    int nslabhitted =0;

    for(int i=0; i<7; i++ ){
      if( bool_hit_slab[i]==true) nslabhitted++;
    }

    if( nslabhitted<7) continue;

    // ---------------------------------------------
    //   the energy and other variables,
    double energy_sum_tmp=0;
    double weight_X0_sum_tmp=0;
    double energy_X0_sum_tmp=0;
    double xm=0, ym=0., zm=0.;

    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      
      if(hit_energy[ihit]>mipcut && hit_isMasked[ihit]==0 ) {
    	  energy_sum_tmp += hit_energy[ihit];
    	  energy_X0_sum_tmp += w[int(hit_z[ihit])] * hit_energy[ihit];
	  weight_X0_sum_tmp += w[int(hit_z[ihit])];

	  xm +=  - hit_x[ihit] * hit_energy[ihit] * w[int(hit_z[ihit])];
	  ym +=  - hit_y[ihit] * hit_energy[ihit] * w[int(hit_z[ihit])];
	  zm +=  hit_z[ihit] * hit_energy[ihit] * w[int(hit_z[ihit])];
      }
    }

    //-------------------------------------------
    //calcualte min distance with other hits 
    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      if(hit_energy[ihit]>mipcut && hit_isMasked[ihit]==0 ) {
	double xi=hit_x[ihit];
	double yi=hit_y[ihit];
	double zi=hit_z[ihit];
	
	double mindist=1000000000.;
	double dist=-10;
	for(int jhit=0; jhit< nhit_chan; jhit ++) {
	  if(hit_energy[jhit]>mipcut && hit_isMasked[jhit]==0 && ihit!=jhit) {
	    double xj=hit_x[jhit];
	    double yj=hit_y[jhit];
	    double zj=hit_z[jhit];
	    dist= sqrt( (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + 150*150*(zi-zj)*(zi-zj) );
	    if(dist < mindist) mindist = dist;
	  }
	}
	min_dist_vs_energy->Fill(mindist,energy_X0_sum_tmp);
      }
    }
    // ---------------------------------------------
    //  recalculate the energy and other variables, with a cut in the distance (avoid isolated cells)
    energy_sum_tmp=0;
    weight_X0_sum_tmp=0;
    energy_X0_sum_tmp=0;
    xm=0, ym=0., zm=0.;

    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      
      if(hit_energy[ihit]>mipcut && hit_isMasked[ihit]==0 ) {
    	double xi=hit_x[ihit];
    	double yi=hit_y[ihit];
    	double zi=hit_z[ihit];
	
    	double mindist=1000000000.;
    	double dist=-10;
    	for(int jhit=0; jhit< nhit_chan; jhit ++) {
    	  if(hit_energy[jhit]>mipcut && hit_isMasked[jhit]==0 && ihit!=jhit) {
    	    double xj=hit_x[jhit];
    	    double yj=hit_y[jhit];
    	    double zj=hit_z[jhit];
    	    dist= sqrt( (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + 150*150*(zi-zj)*(zi-zj) );
    	    if(dist < mindist) mindist = dist;
    	  }
    	}
	if(mindist<7.) {
    	  energy_sum_tmp += hit_energy[ihit];
    	  energy_X0_sum_tmp += w[int(hit_z[ihit])] * hit_energy[ihit];
	  weight_X0_sum_tmp += w[int(hit_z[ihit])];

    	  xm +=  - hit_x[ihit] * hit_energy[ihit] * w[int(hit_z[ihit])];
    	  ym +=  - hit_y[ihit] * hit_energy[ihit] * w[int(hit_z[ihit])];
    	  zm +=  hit_z[ihit] * hit_energy[ihit] * w[int(hit_z[ihit])];
	}
      }
    }

    if(energy_sum_tmp<0.5) continue;

    xm = xm / energy_X0_sum_tmp;
    ym = ym / energy_X0_sum_tmp;
    zm = zm / energy_X0_sum_tmp;


    //    if( xm > xm_min && xm < xm_max && ym > ym_min && ym < ym_max && zm > zm_min && zm < zm_max ) 
    //  energy_center->Fill(energy_X0_sum_tmp);
   
    if( zm > zm_min && zm < zm_max ) 
      energy_center->Fill(energy_X0_sum_tmp);

  }


  // Signal analysis
  TFile *file = new TFile(TString::Format("%s%s_%s_%s_mipcut%0.1f_showers.root",folder.Data(),configuration.Data(),gridpoint.Data(),energy_string.Data(),mipcut) , "RECREATE");
  file->cd();

  // energy->GetXaxis()->SetTitle(TString::Format("E^{raw}/MIP = #sum_{i=layers} #omega_{i} #(#sum_{j=cells (E_{j}/MIP>%0.1f)} E_{j}/MIP)#)",mipcut));
  // energy->GetYaxis()->SetTitle("# entries");
  energy->Write();
  
  //energy_center->GetYaxis()->SetTitle("# entries");
  //energy_center->GetXaxis()->SetTitle(TString::Format("E^{raw}/MIP = #sum_{i=layers} #omega_{i} #(){#sum_{j=cells (E_{j}/MIP>%0.1f)} E_{j}/MIP)}",mipcut));
  energy_center->Write();


  energy_profile_z->Write();
  energy_profile_X0->Write();

  energy_xy->Write();
  energy_xz->Write();
  energy_yz->Write();

  energy_xX0->Write();
  energy_yX0->Write();

  xbarycenter->Write();
  ybarycenter->Write();
  zbarycenter->Write();
  
  xybarycenter->Write();
  xzbarycenter->Write();
  yzbarycenter->Write();

  E_xbarycenter->Write();
  E_ybarycenter->Write();
  E_zbarycenter->Write();

  n_x->Write();
  n_y->Write();
  n_z->Write();
  nhitstotal->Write();
  nhits_vs_energy->Write();
  min_dist_vs_energy->Write();
  file->Close();
  
}

void protoAnalysis::SimpleDistributionsTrack(TString outputname="")
{

  TH1F* SCA_track = new TH1F("SCA_track","SCA_track",15,-0.5,14.5);
  TH1F* BCID_track = new TH1F("BCID_track","BCID_track",60,25,3025);
  TH2F* SCA_BCID_track = new TH2F("SCA_BCID_track","SCA_BCID_track",15,-0.5,14.5,60,25,3025);
  TH2F* BCID_PREV_track = new TH2F("BCID_PREV_track","BCID_PREV_track",1200,2.5,6002.5,200,2.5,1002.5);

  // --------------
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //-----------------


 
  // -----------------------------------------------------------------------------------------------------   
  // Signal readout
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
     
    // for(int minhit=5;minhit<8; minhit++) {
    int minhit=7;
    if( nhit_slab == minhit ) {
	bool track=true;
	
	int nout=0;
	int ngoodhit=0;

	for(int ihit=0; ihit< nhit_chan; ihit ++) {
	  if( fabs(hit_x[ihit]-hit_x[0]) > 0 || fabs(hit_y[ihit]-hit_y[0]) > 0) nout++;
	  if( hit_energy[ihit]> 0.5 > 0) ngoodhit++;

	}

	if(nout > 0 && ngoodhit!=7 ) track=false;
	
	if(track==true && bcid>1260 && bcid<2850) {
	  BCID_PREV_track->Fill(bcid,bcid-prev_bcid);
	  BCID_track->Fill(bcid);
	  for(int ihit=0; ihit< nhit_chan; ihit ++) {
	    SCA_track->Fill(hit_sca[ihit]);
	    SCA_BCID_track->Fill(hit_sca[ihit],bcid);
	  }
	}



    }
  //
    
  }

  // -----------------------------------------------------------------------------------------------------   
  // Signal analysis
  TFile *signalfile_summary = new TFile("results_mipcalibration/SimpleDist"+outputname+".root" , "RECREATE");
  signalfile_summary->cd();

  SCA_track->GetXaxis()->SetTitle("SCA");
  SCA_track->GetYaxis()->SetTitle("# entries");
  SCA_track->SetTitle("SCA distribution for track-like events");
  SCA_track->SetName("SCA_track");
  SCA_track->Write();

  BCID_track->GetXaxis()->SetTitle("bcid");
  BCID_track->GetYaxis()->SetTitle("# entries");
  BCID_track->SetTitle("bcid distribution for track-like events");
  BCID_track->SetName("BCID_track");
  BCID_track->Write();

  SCA_BCID_track->GetXaxis()->SetTitle("SCA");
  SCA_BCID_track->GetYaxis()->SetTitle("bcid");
  SCA_BCID_track->SetTitle("SCA vs bcid distribution for track-like events");
  SCA_BCID_track->SetName("SCA_BCID_track");
  SCA_BCID_track->Write();
  
  BCID_PREV_track->GetXaxis()->SetTitle("bcid");
  BCID_PREV_track->GetYaxis()->SetTitle("bcid-prev_bcid");
  BCID_PREV_track->SetTitle("bcid vs bcid-prev_bcid distribution for track-like events");
  BCID_PREV_track->SetName("BCID_PREV_track");
  BCID_PREV_track->Write();

  signalfile_summary->Close();

}


void protoAnalysis::SimpleMIPAnalysis(TString outputname="")
{

  ofstream fout_mip("results_mipcalibration/MIP"+outputname+".txt",ios::out);
  fout_mip<<"#mip results, all channels together SimpleMIPAnalysis using tracks"<<endl;
  fout_mip<<"#mpv empv widthmpv widthgauss chi2ndf nentries"<<endl;

  // --------------
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //-----------------

  //Objects to store the chip/channel results:
  // numbers 
  std::vector<std::vector<std::vector<TH1F*> > > mip_histo;
  //  TH1F* mip_histo_all = new TH1F("mip_histo_all","mip_histo_all",500,0.05,10.05);
  //TH1F* mip_histo_all2 = new TH1F("mip_histo_all2","mip_histo_all2",500,0.05,10.05);
  //TH1F* mip_histo_all_subtracted1 = new TH1F("mip_histo_all_subtracted1","mip_histo_all_subtracted1",500,0.05,10.05);
  //TH1F* mip_histo_all_subtracted2 = new TH1F("mip_histo_all_subtracted2","mip_histo_all_subtracted2",500,0.05,10.05);

  TH1F* mip_histo_all = new TH1F("mip_histo_all","mip_histo_all",250,0.1,10.1);
  //  TH1F* mip_histo_all2 = new TH1F("mip_histo_all2","mip_histo_all2",250,0.1,10.1);
  //TH1F* mip_histo_all_subtracted1 = new TH1F("mip_histo_all_subtracted1","mip_histo_all_subtracted1",250,0.1,10.1);
  //TH1F* mip_histo_all_subtracted2 = new TH1F("mip_histo_all_subtracted2","mip_histo_all_subtracted2",250,0.1,10.1);


  for(int islab=0; islab<7; islab++) {
    std::vector<std::vector<TH1F*>  >slab_mip_histo;
    for(int ichip=0; ichip<16; ichip++) {
      std::vector<TH1F*>  chip_mip_histo;
      for(int ichn=0;ichn<64;ichn++) {
	TString histo_title=TString::Format("charge_slab%i_chip%i_chn%i",islab,ichip,ichn);
	TH1F *chn_mip = new TH1F(histo_title,histo_title,4197,-100.5,4096.5);
	chip_mip_histo.push_back(chn_mip);
      }      
      slab_mip_histo.push_back(chip_mip_histo);
    }
    mip_histo.push_back(slab_mip_histo);
  }

 
  // -----------------------------------------------------------------------------------------------------   
  // Signal readout
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
     
    for(int minhit=5;minhit<8; minhit++) {
      if( (nhit_chan > (minhit-1) && nhit_chan< 10 ) && nhit_slab == minhit ) {
	bool track=true;
	int nout=0;

	for(int ihit=0; ihit< nhit_chan; ihit ++) {
	  if( fabs(hit_x[ihit]-hit_x[0]) > 0 || fabs(hit_y[ihit]-hit_y[0]) > 0) nout++;
	}

	if(nout> (nhit_chan-nhit_slab) ) track=false;
	
	  if(track==true && bcid>1250 && bcid<2850) {
	    for(int ihit=0; ihit< nhit_chan; ihit ++) {
	      mip_histo_all->Fill(hit_energy[ihit]);
	      //	      mip_histo_all2->Fill(hit_energy[ihit]);
	      //mip_histo_all_subtracted1->Fill(hit_energy[ihit]);
	      //mip_histo_all_subtracted2->Fill(hit_energy[ihit]);
	    }
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

  // Setting fit range and start values
  Double_t fr[2];
  Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];   
  pllo[0]=0.001; pllo[1]=0.5; pllo[2]=10.0; pllo[3]=0.001;
  plhi[0]=1.; plhi[1]=5.0; plhi[2]=100000000.0; plhi[3]=1.0;

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
  TSpectrum *s_iter1 = new TSpectrum(1);
  int npeaks_iter1_0 = s_iter1->Search(hist_peak1_iter1,1,"",0.8); 
  Double_t *mean_peak_iter1_0=s_iter1->GetPositionX();
  fr[0]=mean_peak_iter1_0[0]-0.15*hist_peak1_iter1->GetRMS();
  fr[1]=mean_peak_iter1_0[0]+0.4*hist_peak1_iter1->GetRMS();
  sv[0]=hist_peak1_iter1->GetRMS()/10;
  sv[1]=hist_peak1_iter1->GetMean();
  sv[2]=hist_peak1_iter1->Integral("width");
  sv[3]=hist_peak1_iter1->GetRMS()/10.;
  
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
 
  fout_mip<<" "<<mpv<<" "<<empv<<" "<<wmpv<<" "<<wg<<" "<<chi2ndf<<" "<<mipentries<<endl;
  cout<<npeaks_iter1_0<<" "<<mean_peak_iter1_0[0]<<endl;
  cout<<" "<<mpv<<" "<<empv<<" "<<wmpv<<" "<<wg<<" "<<chi2ndf<<" "<<mipentries<<endl;

  
  // ------------------------------------------------------
  // Peak 2
  // -----------------------------------------------------
  TH1F *hist_peak2_iter1 = (TH1F*)mip_histo_all->Clone();
  hist_peak2_iter1->Add(fpeak1_iter1,-1);

  TSpectrum *s_iter1_1 = new TSpectrum(1);
  int npeaks_iter1_1 = s_iter1_1->Search(hist_peak2_iter1,0.5,"",0.8); 
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

  fout_mip<<" "<<mpv<<" "<<empv<<" "<<wmpv<<" "<<wg<<" "<<chi2ndf<<" "<<mipentries<<endl;
  cout<<npeaks_iter1_1<<" "<<mean_peak_iter1_1[0]<<endl;
  cout<<" "<<mpv<<" "<<empv<<" "<<wmpv<<" "<<wg<<" "<<chi2ndf<<" "<<mipentries<<endl;

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

  fout_mip<<" "<<mpv<<" "<<empv<<" "<<wmpv<<" "<<wg<<" "<<chi2ndf<<" "<<mipentries<<endl;
  cout<<npeaks_iter1_2<<" "<<mean_peak_iter1_2[0]<<endl;
  cout<<" "<<mpv<<" "<<empv<<" "<<wmpv<<" "<<wg<<" "<<chi2ndf<<" "<<mipentries<<endl;

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
  first_iteration_peak3->Write();

  signalfile_summary->cd();

  signalfile_summary->Close();

}



// LANDAU STUFF
//---------------------------------------------------------------------------------------




TF1* protoAnalysis::langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
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

