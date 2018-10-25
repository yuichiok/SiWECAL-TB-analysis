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

void protoAnalysis::ReadMap(TString filename="../fev10_chip_channel_x_y_mapping.txt") 
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


bool protoAnalysis::IsHit(int ihit) {

  if( (hit_isMasked[ihit]==0) && hit_isHit[ihit]==1) return true;
  return false;
}

bool protoAnalysis::IsPedestal(int ihit) {

  if( (hit_isMasked[ihit]==0) &&  hit_isHit[ihit]==0) return true;
  return false;
  
}

bool protoAnalysis::ShowerBasicSelection(double mipcut=0.5, int bcid_max=99999999, int nslabs_selection=3) {

  if(bcid>bcid_max) return false;
  

  bool bool_hit_slab[10];
  bool_hit_slab[0] = false;
  bool_hit_slab[1] = false;
  bool_hit_slab[2] = false;
  bool_hit_slab[3] = false;
  bool_hit_slab[4] = false;
  bool_hit_slab[5] = false;
  bool_hit_slab[6] = false;
  bool_hit_slab[7] = false;
  bool_hit_slab[8] = false;
  bool_hit_slab[9] = false;

  for(int ihit=0; ihit< nhit_chan; ihit ++) {
    int z=hit_z[ihit];

    if(hit_energy[ihit]>mipcut && IsHit(ihit)==true ) {
      bool_hit_slab[z]=true;
    }
  }
  

  int nslabhitted =0;
  
  for(int i=0; i<10; i++ ){
    if( bool_hit_slab[i]==true) nslabhitted++;
  }
  
  if( nslabhitted < nslabs_selection) return false;
  return true;

}


bool protoAnalysis::ShowerAntiSelection(double mipcut=0.5, int bcid_max=99999999, int nslabs_selection=3) {

  if(bcid>bcid_max) return false;
 
  bool bool_hit_slab[10];
  bool_hit_slab[0] = false;
  bool_hit_slab[1] = false;
  bool_hit_slab[2] = false;
  bool_hit_slab[3] = false;
  bool_hit_slab[4] = false;
  bool_hit_slab[5] = false;
  bool_hit_slab[6] = false;
  bool_hit_slab[7] = false;
  bool_hit_slab[8] = false;
  bool_hit_slab[9] = false;

  for(int ihit=0; ihit< nhit_chan; ihit ++) {
    int z=hit_z[ihit];

    if(hit_energy[ihit]>mipcut && IsHit(ihit)==true ) {
      bool_hit_slab[z]=true;
    }
  }

  int nslabhitted =0;
  
  for(int i=0; i<10; i++ ){
    if( bool_hit_slab[i]==true) nslabhitted++;
  }
  
  if( nslabhitted == nslabs_selection) return true;

  return false;

}


double protoAnalysis::MinDistCalc(int ihit, double mipcut) {

  double xorient=1;
  if(hit_z[ihit]>2 && hit_z[ihit]<7) xorient=-1;
  double xi=-hit_x[ihit];
  double yi=-xorient * hit_y[ihit];
  double zi=hit_z[ihit];
  
  double mindist=1000.;
  double dist=-10;
  for(int jhit=0; jhit< nhit_chan; jhit ++) {
    if(hit_energy[jhit]>mipcut && IsHit(ihit)==true && ihit!=jhit) {
      double xorient=1;
      if(hit_z[ihit]>2 && hit_z[ihit]<7) xorient=-1;
      double xj=-hit_x[jhit];
      double yj=-xorient * hit_y[jhit];
      double zj=hit_z[jhit];
      dist= sqrt( (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + 150*150*(zi-zj)*(zi-zj) );
      if(dist < mindist) mindist = dist;
    }
  }

  return mindist;
}


void protoAnalysis::ShowerDistributions(TString folder="", TString energy_string="40GeV", double mipcut=0.5, int bcid_max=99999999, int nslabs_selection=4, bool RMIsolatedHits=false)
{


  // --------------
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //-----------------

  Float_t w[10];
  Float_t thickness[10];
  thickness[0]=6.3;
  thickness[1]=6.3;
  thickness[2]=6.3;
  thickness[3]=6.3;
  thickness[4]=6.3;
  thickness[5]=6.3;
  thickness[6]=6.3;
  thickness[7]=6.3;
  thickness[8]=6.3;
  thickness[9]=6.3;

  w[0]=(1/3.5)*thickness[0];
  w[1]=(1/3.5)*thickness[1];
  w[2]=(1/3.5)*thickness[2];
  w[3]=(1/3.5)*thickness[3];
  w[4]=(1/3.5)*thickness[4];
  w[5]=(1/3.5)*thickness[5];
  w[6]=(1/3.5)*thickness[6];
  w[7]=(1/3.5)*thickness[7];
  w[8]=(1/3.5)*thickness[8];
  w[9]=(1/3.5)*thickness[9];
  
  
  Float_t edge =0.0001;

  Float_t bins[] = { 0, w[0]+edge, w[0]+edge+w[1], w[0]+edge+w[1]+w[2], w[0]+edge+w[1]+w[2]+w[3], w[0]+edge+w[1]+w[2]+w[3]+w[4],  w[0]+edge+w[1]+w[2]+w[3]+w[4]+w[5],  w[0]+edge+w[1]+w[2]+w[3]+w[4]+w[5]+w[6], w[0]+edge+w[1]+w[2]+w[3]+w[4]+w[5]+w[6]+w[7], w[0]+edge+w[1]+w[2]+w[3]+w[4]+w[5]+w[6]+w[7]+w[8], w[0]+edge+w[1]+w[2]+w[3]+w[4]+w[5]+w[6]+w[7]+w[8]+w[9]};
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
  TH1F* energy = new TH1F("energy","energy",100,10,2010);
  TH1F* energy_center_05 = new TH1F("energy_center_zm05","energy_center_zm05",100,20,2010);
  TH1F* energy_center_1 = new TH1F("energy_center_zm1","energy_center_zm1",100,20,2010);
  TH1F* energy_center_2 = new TH1F("energy_center_zm2","energy_center_zm2",100,20,2010);
  
  TH1F* energy_profile_z = new TH1F("energy_profile_z","energy_profile_z",11,-0.5,10.5);//binnumz,binsz);
  TH1F* energy_profile_z_center_05 = new TH1F("energy_profile_z_center_05","energy_profile_z_center_05",11,-0.5,10.5);//binnumz,binsz);
  TH1F* energy_profile_z_center_1 = new TH1F("energy_profile_z_center_1","energy_profile_z_center_1",11,-0.5,10.5);//binnumz,binsz);
  TH1F* energy_profile_z_center_2 = new TH1F("energy_profile_z_center_2","energy_profile_z_center_2",11,-0.5,10.5);//binnumz,binsz);

  TH1F* time_dist = new TH1F("time_dist","time_dist",2500,0,25000);

  energy_profile_z->Sumw2();
  energy_profile_z_center_05->Sumw2();
  energy_profile_z_center_1->Sumw2();
  energy_profile_z_center_2->Sumw2();

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
  TH1F* nhitstotal_center_05 = new TH1F("nhitstotal_center_zm05","nhitstotal_center_zm05",100,20,2010);
  TH1F* nhitstotal_center_1 = new TH1F("nhitstotal_center_zm1","nhitstotal_center_zm1",100,20,2010);
  TH1F* nhitstotal_center_2 = new TH1F("nhitstotal_center_zm2","nhitstotal_center_zm2",100,20,2010);

  
  Float_t binsnhit[501];
  binsnhit[0]=-0.5;
  for(int ibins=1;ibins<501;ibins++) binsnhit[ibins]=binsnhit[ibins-1]+1;
  Int_t  binnhit = 500;
  
  Float_t binsz2[101];
  binsz2[0]=0;
  for(int ibins=1;ibins<101;ibins++) binsz2[ibins]=binsz2[ibins-1]+0.2;
  Int_t  binnumz2 = 100;
  
  TH2F* nhits_vs_energy= new TH2F("nhits_vs_energy","nhits_vs_energy",51,-2.5,252.5,101,-5,1005);
  TH2F* zbarycenter_vs_nhits= new TH2F("zbarycenter_vs_nhits","zbarycenter_vs_nhits",binnumz2,binsz2,binnhit,binsnhit);
  TH2F* min_dist_vs_energy= new TH2F("min_dist_vs_energy","min_dist_vs_energy",250,0,250,100,-5,995);

  int number_hits=0;
  // -----------------------------------------------------------------------------------------------------
  Long64_t nbytes = 0, nb = 0;

 
  // Signal readout
  nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(ShowerBasicSelection(mipcut,bcid_max,nslabs_selection)==false) continue;

    time_dist->Fill(bcid*0.2);

    // ---------------------------------------------
    //   the energy and other variables,
    double energy_sum_tmp=0;
    double weight_X0_sum_tmp=0, energy_X0_sum_tmp_dist=0;
    double energy_X0_sum_tmp=0;
    double xm=0, ym=0., zm=0.;

    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      if(hit_energy[ihit]>mipcut && IsHit(ihit)==true) {
	energy_X0_sum_tmp_dist += w[int(hit_z[ihit])] * (hit_energy[ihit]);
      }
    }
    //-------------------------------------------
    //calcualte min distance with other hits
    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      if(hit_energy[ihit]>mipcut && IsHit(ihit)==true ) {
	double mindist=MinDistCalc(ihit,mipcut);
	min_dist_vs_energy->Fill(mindist,energy_X0_sum_tmp_dist);
      }
    }

    int n_hitstotal=0;

    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      if(hit_energy[ihit]>mipcut && IsHit(ihit)==true) {
	double mindist=MinDistCalc(ihit,mipcut);
	double mindist_th=10000.;
	if(RMIsolatedHits==true) mindist_th=7.0;

	if(mindist< mindist_th) {
	  n_hitstotal++;

	  energy_sum_tmp += (hit_energy[ihit]);
	  energy_X0_sum_tmp += w[int(hit_z[ihit])] * (hit_energy[ihit]);
	  weight_X0_sum_tmp += w[int(hit_z[ihit])];
	  
	  double xorient=1;
	  if(hit_z[ihit]>2 && hit_z[ihit]<7) xorient=-1;

	  xm +=  - hit_x[ihit] * (hit_energy[ihit]) * w[int(hit_z[ihit])];
	  ym +=  - xorient * hit_y[ihit] * (hit_energy[ihit]) * w[int(hit_z[ihit])];
	  zm +=    hit_z[ihit] * (hit_energy[ihit]) * w[int(hit_z[ihit])];
	}
      }
    }

    if(energy_sum_tmp<0.5) continue;
    number_hits++;

    nhitstotal->Fill(n_hitstotal);
    nhits_vs_energy->Fill(n_hitstotal,energy_X0_sum_tmp);


    xm = xm / energy_X0_sum_tmp;
    ym = ym / energy_X0_sum_tmp;
    zm = zm / energy_X0_sum_tmp;

    zbarycenter_vs_nhits->Fill(zm,n_hitstotal);

    for(int ihit=0; ihit< nhit_chan; ihit ++) {

      if(hit_energy[ihit]>mipcut && IsHit(ihit)==true ) {
	double zlayer=hit_z[ihit];
	energy_profile_z->Fill(zlayer,hit_energy[ihit]/energy_sum_tmp);
	energy_profile_X0->Fill(hit_x0[ihit],hit_energy[ihit]/energy_sum_tmp);
	
	double xorient=1;
	if(hit_z[ihit]>2 && hit_z[ihit]<7) xorient=-1;

	energy_xy->Fill(-hit_x[ihit],-xorient * hit_y[ihit],hit_energy[ihit]/energy_sum_tmp );
	energy_xz->Fill(-hit_x[ihit],hit_z[ihit],hit_energy[ihit]/energy_sum_tmp );
	energy_yz->Fill(-xorient * hit_y[ihit],hit_z[ihit],hit_energy[ihit]/energy_sum_tmp );     
	
	energy_xX0->Fill(-hit_x[ihit],hit_x0[ihit],hit_energy[ihit]/energy_sum_tmp);
	energy_yX0->Fill(-xorient * hit_y[ihit],hit_x0[ihit],hit_energy[ihit]/energy_sum_tmp);
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
  /*
  energy_profile_z->Scale(1./number_hits);
  energy_profile_X0->Scale(1./number_hits);

  energy_xy->Scale(1./number_hits);
  energy_xz->Scale(1./number_hits);
  energy_yz->Scale(1./number_hits);

  energy_xX0->Scale(1./number_hits);
  energy_yX0->Scale(1./number_hits);
  */
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
  
  double xm_max=fit_xbarycenter->GetParameter(1)+1.*fit_xbarycenter->GetParameter(2);
  double xm_min=fit_xbarycenter->GetParameter(1)-1.*fit_xbarycenter->GetParameter(2);

  TF1 *fit_ybarycenter = new TF1("fit_ybarycenter","gaus");
  ybarycenter->Fit(fit_ybarycenter,"M");
  
  double ym_max=fit_ybarycenter->GetParameter(1)+1.*fit_ybarycenter->GetParameter(2);
  double ym_min=fit_ybarycenter->GetParameter(1)-1.*fit_ybarycenter->GetParameter(2);

  TF1 *fit_zbarycenter = new TF1("fit_zbarycenter","gaus");
  zbarycenter->Fit(fit_zbarycenter,"M");

  double zm_max_05=fit_zbarycenter->GetParameter(1)+0.5*fit_zbarycenter->GetParameter(2);
  double zm_min_05=fit_zbarycenter->GetParameter(1)-0.5*fit_zbarycenter->GetParameter(2);

  double zm_max_1=fit_zbarycenter->GetParameter(1)+1.*fit_zbarycenter->GetParameter(2);
  double zm_min_1=fit_zbarycenter->GetParameter(1)-1.*fit_zbarycenter->GetParameter(2);

  double zm_max_2=fit_zbarycenter->GetParameter(1)+2.*fit_zbarycenter->GetParameter(2);
  double zm_min_2=fit_zbarycenter->GetParameter(1)-2.*fit_zbarycenter->GetParameter(2);

  
  // -----------------------------------------------------------------------------------------------------   
  int number_events_zm_05=0;
  int number_events_zm_1=0;
  int number_events_zm_2=0;

  // Signal readout
  nbytes = 0; nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(ShowerBasicSelection(mipcut,bcid_max,nslabs_selection)==false) continue;

    // ---------------------------------------------
    //   the energy and other variables,
    double energy_sum_tmp=0;
    double weight_X0_sum_tmp=0;
    double energy_X0_sum_tmp=0;

    double xm=0, ym=0, zm=0;

    int n_hitstotal=0;

    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      
      if(hit_energy[ihit]>mipcut && IsHit(ihit)==true) {

	double mindist=MinDistCalc(ihit,mipcut);
	double mindist_th=10000.;
	if(RMIsolatedHits==true) mindist_th=7.0;

	if(mindist< mindist_th && IsHit(ihit)==true && hit_energy[ihit]>mipcut ) {
	  n_hitstotal++;
	  energy_sum_tmp += (hit_energy[ihit]);
	  energy_X0_sum_tmp += w[int(hit_z[ihit])] * (hit_energy[ihit]);
	  weight_X0_sum_tmp += w[int(hit_z[ihit])];
	  double xorient=1;
	  if(hit_z[ihit]>2 && hit_z[ihit]<7) xorient=-1;

	  xm +=  - hit_x[ihit] * hit_energy[ihit] * w[int(hit_z[ihit])];
	  ym +=  - xorient * hit_y[ihit] * hit_energy[ihit] * w[int(hit_z[ihit])];
	  zm +=  hit_z[ihit] * hit_energy[ihit] * w[int(hit_z[ihit])];
	}
      }
    }

    if(energy_sum_tmp<0.5) continue;
    
    xm = xm / energy_X0_sum_tmp;
    ym = ym / energy_X0_sum_tmp;
    zm = zm / energy_X0_sum_tmp;

   
    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      if(hit_energy[ihit]>mipcut && IsHit(ihit)==true ) {
	double zlayer=hit_z[ihit];
	if( zm > zm_min_05 && zm < zm_max_05 )	energy_profile_z_center_05->Fill(zlayer,hit_energy[ihit]/energy_sum_tmp);
	if( zm > zm_min_1 && zm < zm_max_1 )	energy_profile_z_center_1->Fill(zlayer,hit_energy[ihit]/energy_sum_tmp);
	if( zm > zm_min_2 && zm < zm_max_2 )	energy_profile_z_center_2->Fill(zlayer,hit_energy[ihit]/energy_sum_tmp);
      }
    }
   
    if( zm > zm_min_05 && zm < zm_max_05 ) {
      energy_center_05->Fill(energy_X0_sum_tmp);
      nhitstotal_center_05->Fill(n_hitstotal);
      number_events_zm_05++;
    }

    if( zm > zm_min_1 && zm < zm_max_1 ) {
      energy_center_1->Fill(energy_X0_sum_tmp);
      nhitstotal_center_1->Fill(n_hitstotal);
      number_events_zm_1++;
    }

    if( zm > zm_min_2 && zm < zm_max_2 ) {
      energy_center_2->Fill(energy_X0_sum_tmp);
      nhitstotal_center_2->Fill(n_hitstotal);
      number_events_zm_2++;
    }
    
  }

  energy_profile_z_center_05->Scale(1./number_events_zm_05);
  energy_profile_z_center_1->Scale(1./number_events_zm_1);
  energy_profile_z_center_2->Scale(1./number_events_zm_2);

  // Signal analysis
  TFile *file = new TFile(TString::Format("%s/%s_mipcut%0.1f_showers.root",folder.Data(),energy_string.Data(),mipcut) , "RECREATE");
  file->cd();

  time_dist->GetXaxis()->SetTitle("#mu s");
  time_dist->Write();
  
  energy->Write();
  energy_center_05->Write();
  energy_center_1->Write();
  energy_center_2->Write();
  
  energy_profile_z->Write();
  energy_profile_z_center_05->Write();
  energy_profile_z_center_1->Write();
  energy_profile_z_center_2->Write();

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
  nhitstotal_center_05->Write();
  nhitstotal_center_1->Write();
  nhitstotal_center_2->Write();
  
  nhits_vs_energy->Write();
  zbarycenter_vs_nhits->Write();
  min_dist_vs_energy->Write();
  file->Close();
  
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

