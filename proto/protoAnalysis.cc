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


void protoAnalysis::ReadCalibrated(TString filename) 
{

  for(int i=0; i<7; i++) {
    for(int j=0; j<16; j++) {
      for(int k=0; k<64; k++) {
	calibrated[i][j][k] = 1;
      }
    }
  }

  cout<<filename<<endl;
  

  for(int iz=0; iz<7; iz++) {
    TString dif="dif_1_1_1.txt";
    if(iz==1) dif="dif_1_1_2.txt";
    if(iz==2) dif="dif_1_1_3.txt";
    if(iz==3) dif="dif_1_1_4.txt";
    if(iz==4) dif="dif_1_1_5.txt";
    if(iz==5) dif="dif_1_2_1.txt";
    if(iz==6) dif="dif_1_2_2.txt";

    TString filename_f=filename+dif;
    cout<<filename_f<<endl;
    std::ifstream reading_file(filename_f);
    if(!reading_file){
      cout<<" dameyo - damedame"<<endl;
    } else cout<<filename_f<<endl;
    
    Int_t tmp_chip = 0,tmp_channel = 0;
    Double_t tmp_value = 0;
    Double_t tmp_mip_value = 0;

    Int_t tmp_calibrated = 0;
    TString tmpst;
    //  #mip results dif_1_1_1
    //#chip channel mpv empv widthmpv chi2ndf nentries

    reading_file >> tmpst >> tmpst >> tmpst  ;
    reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst  ;
    
    cout<<"Read Calibrated: "<<filename<<endl;
    while(reading_file){
      reading_file >> tmp_chip >> tmp_channel >> tmp_mip_value >> tmp_value >> tmp_value >> tmp_value >> tmp_value ;
      if(tmp_mip_value>0) calibrated[iz][tmp_chip][tmp_channel] = 0;
    }
  

    double ncalibrated=0.;
    for(int i=0; i<16; i++) {
      for(int j=0; j<64; j++) {
	if(calibrated[iz][i][j]==1) ncalibrated++;
      }
    }
    ncalibrated=100.*ncalibrated/1024;
    cout<< 100-ncalibrated<<"% of channels are calibrated in layer "<<iz<<endl;
  }
  
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

  int z=hit_z[ihit];
  if(z==9) z=6;
  if( (hit_isMasked[ihit]==0 && calibrated[z][hit_chip[ihit]][hit_chan[ihit]]==0) && hit_isHit[ihit]==1) return true;
  return false;
}

bool protoAnalysis::IsPedestal(int ihit) {

  int z=hit_z[ihit];
  if(z==9) z=6;
  if( (hit_isMasked[ihit]==0 && calibrated[z][hit_chip[ihit]][hit_chan[ihit]]==0) &&  hit_isHit[ihit]==0) return true;
  return false;
  
}

bool protoAnalysis::ShowerBasicSelection(double mipcut=0.5, int bcid_max=2900, int nslabs_selection=6) {

  if(bcid<1250)  return false;
  if(bcid>bcid_max) return false;
  

  bool bool_hit_slab[7];
  bool_hit_slab[0] = false;
  bool_hit_slab[1] = false;
  bool_hit_slab[2] = false;
  bool_hit_slab[3] = false;
  bool_hit_slab[4] = false;
  bool_hit_slab[5] = false;
  bool_hit_slab[6] = false;
  
  for(int ihit=0; ihit< nhit_chan; ihit ++) {
    int z=hit_z[ihit];
    if(z==9) z=6;
    if(hit_energy[ihit]>mipcut && IsHit(ihit)==true ) {
      bool_hit_slab[z]=true;
    }
  }
  

  int nslabhitted =0;
  
  for(int i=0; i<7; i++ ){
    if( bool_hit_slab[i]==true) nslabhitted++;
  }
  
  if( nslabhitted < nslabs_selection) return false;

  return true;

}


bool protoAnalysis::ShowerAntiSelection(double mipcut=0.5, int bcid_max=2900, int nslabs_selection=6) {

  if(bcid<1250)  return false;
  if(bcid>bcid_max) return false;
  

  bool bool_hit_slab[7];
  bool_hit_slab[0] = false;
  bool_hit_slab[1] = false;
  bool_hit_slab[2] = false;
  bool_hit_slab[3] = false;
  bool_hit_slab[4] = false;
  bool_hit_slab[5] = false;
  bool_hit_slab[6] = false;
  
  for(int ihit=0; ihit< nhit_chan; ihit ++) {
    int z=hit_z[ihit];
    if(z==9) z=6;
    if( hit_energy[ihit]>mipcut && IsHit(ihit)==true ) {
      bool_hit_slab[z]=true;
    }
  }
  

  int nslabhitted =0;
  
  for(int i=0; i<7; i++ ){
    if( bool_hit_slab[i]==true) nslabhitted++;
  }
  
  if( nslabhitted == nslabs_selection) return true;

  return false;

}


double protoAnalysis::MinDistCalc(int ihit, double mipcut) {

  double xi=hit_x[ihit];
  double yi=hit_y[ihit];
  double zi=hit_z[ihit];
  
  double mindist=1000.;
  double dist=-10;
  for(int jhit=0; jhit< nhit_chan; jhit ++) {
    if(hit_energy[jhit]>mipcut && IsHit(ihit)==true && ihit!=jhit) {
      double xj=hit_x[jhit];
      double yj=hit_y[jhit];
      double zj=hit_z[jhit];
      dist= sqrt( (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + 150*150*(zi-zj)*(zi-zj) );
      if(dist < mindist) mindist = dist;
    }
  }

  return mindist;
}


void protoAnalysis::SimpleDistributionsShower(TString outputname="")
{

  // --------------
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //-----------------

    
  TH1F* SCA_shower = new TH1F("SCA_shower","SCA_shower",15,-0.5,14.5);
  TH1F* BCID_shower = new TH1F("BCID_shower","BCID_shower",60,25,3025);
  TH2F* SCA_BCID_shower = new TH2F("SCA_BCID_shower","SCA_BCID_shower",15,-0.5,14.5,60,25,3025);
  TH2F* BCID_PREV_shower = new TH2F("BCID_PREV_shower","BCID_PREV_shower",1200,2.5,6002.5,200,2.5,1002.5);



  double mipcut=0.5;
  // Signal readout
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(bcid<1250) continue;
    if(bcid>2900) continue;

    //if(nhit_slab<7) continue;


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
      if(hit_energy[ihit]>mipcut && hit_isMasked[ihit]==0 && hit_isHit[ihit]==1) {
	bool_hit_slab[z]=true;
      }
    }
 

    int nslabhitted =0;

    for(int i=0; i<7; i++ ){
      if( bool_hit_slab[i]==true) nslabhitted++;
    }

    if( nslabhitted<6) continue;

    bool fill=false;
    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      int z=hit_z[ihit];
      if(z==9) z=6;
      
      if(hit_energy[ihit]>mipcut && hit_isMasked[ihit]==0  && hit_isHit[ihit]==1) {
    	double xi=hit_x[ihit];
    	double yi=hit_y[ihit];
    	double zi=hit_z[ihit];
	
    	double mindist=1000000000.;
    	double dist=-10;
    	for(int jhit=0; jhit< nhit_chan; jhit ++) {
    	  if(hit_energy[jhit]>mipcut && hit_isMasked[jhit]==0 && hit_isHit[jhit]==1 && ihit!=jhit) {
    	    double xj=hit_x[jhit];
    	    double yj=hit_y[jhit];
    	    double zj=hit_z[jhit];
    	    dist= sqrt( (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + 150*150*(zi-zj)*(zi-zj) );
	    if(dist < mindist) mindist = dist;
    	  }
    	}

	if(mindist<7. && hit_isMasked[ihit]==0 && hit_isHit[ihit]==1 ) {
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


void protoAnalysis::PedestalAnalysis_showers(TString folder="y", TString configuration="conf1", TString energy_string="3GeV", TString gridpoint ="grid20", double mipcut=0.5, int bcid_max=2900, int nslabs_selection=6, bool RMIsolatedHits=true, int analyze_chip=12, int analyze_layer=2)
{

  ReadMap("../fev10_chip_channel_x_y_mapping.txt");
  ReadCalibrated("../mip_calib/MIP_");

  // --------------
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //-----------------

  int number_hits=0;
  // -----------------------------------------------------------------------------------------------------

  // -----------------------------------------------------------------------------------------------------
  //#******************************************************************************************************
  // PEDESTAL PART:
  // th1f, i,j,k,sca
  TH2F* totalcharge_vs_pedcharge[15];
  TH2F* chargeinanotherchip_vs_pedcharge[15];
  TH2F* averagecharge_vs_pedcharge[15];
  TH2F* nhits_vs_pedcharge[15];
  TH2F* nchannelsElarger5_vs_pedcharge[15];
  TH2F* nchannelsElarger2_vs_pedcharge[15];
  TH2F* nchannelsElarger1_vs_pedcharge[15];
  TH2F* nchannelsElarger05_vs_pedcharge[15];

  TH2F* totalcharge_vs_pedcharge_prevsca[15];
  TH2F* chargeinanotherchip_vs_pedcharge_prevsca[15];
  TH2F* averagecharge_vs_pedcharge_prevsca[15];
  TH2F* nhits_vs_pedcharge_prevsca[15];
  TH2F* nchannelsElarger5_vs_pedcharge_prevsca[15];
  TH2F* nchannelsElarger2_vs_pedcharge_prevsca[15];
  TH2F* nchannelsElarger1_vs_pedcharge_prevsca[15];
  TH2F* nchannelsElarger05_vs_pedcharge_prevsca[15];

  for(int isca=0; isca<15; isca++) {
    totalcharge_vs_pedcharge[isca] = new TH2F(TString::Format("totalcharge_vs_pedcharge_sca%i",isca),TString::Format("totalcharge_vs_pedcharge_sca%i",isca),101,-0.5,100.5,80,-1,1);
    chargeinanotherchip_vs_pedcharge[isca] = new TH2F(TString::Format("chargeinanotherchip_vs_pedcharge_sca%i",isca),TString::Format("chargeinanotherchip_vs_pedcharge_sca%i",isca),101,-0.5,100.5,80,-1,1);
    
    averagecharge_vs_pedcharge[isca] = new TH2F(TString::Format("averagecharge_vs_pedcharge_sca%i",isca),TString::Format("averagecharge_vs_pedcharge_sca%i",isca),20,0,10,80,-1,1);
    nhits_vs_pedcharge[isca] = new TH2F(TString::Format("nhits_vs_pedcharge_sca%i",isca),TString::Format("nhits_vs_pedcharge_sca%i",isca),20,0.5,20.5,80,-1,1);

    nchannelsElarger5_vs_pedcharge[isca] = new TH2F(TString::Format("nchannelsElarger5_vs_pedcharge_sca%i",isca),TString::Format("nchannelsElarger5_vs_pedcharge_sca%i",isca),65,-0.5,64.5,80,-1,1);
    nchannelsElarger2_vs_pedcharge[isca] = new TH2F(TString::Format("nchannelsElarger2_vs_pedcharge_sca%i",isca),TString::Format("nchannelsElarger2_vs_pedcharge_sca%i",isca),65,-0.5,64.5,80,-1,1);
    nchannelsElarger1_vs_pedcharge[isca] = new TH2F(TString::Format("nchannelsElarger1_vs_pedcharge_sca%i",isca),TString::Format("nchannelsElarger1_vs_pedcharge_sca%i",isca),65,-0.5,64.5,80,-1,1);
    nchannelsElarger05_vs_pedcharge[isca] = new TH2F(TString::Format("nchannelsElarger05_vs_pedcharge_sca%i",isca),TString::Format("nchannelsElarger05_vs_pedcharge_sca%i",isca),65,-0.5,64.5,80,-1,1);

  }

  for(int isca=0; isca<15; isca++) {
    totalcharge_vs_pedcharge_prevsca[isca] = new TH2F(TString::Format("totalcharge_vs_pedcharge_prevsca_sca%i",isca),TString::Format("totalcharge_vs_pedcharge_prevsca_sca%i",isca),101,-0.5,100.5,80,-1,1);
    chargeinanotherchip_vs_pedcharge_prevsca[isca] = new TH2F(TString::Format("chargeinanotherchip_vs_pedcharge_prevsca_sca%i",isca),TString::Format("chargeinanotherchip_vs_pedcharge_prevsca_sca%i",isca),101,-0.5,100.5,80,-1,1);
    
    averagecharge_vs_pedcharge_prevsca[isca] = new TH2F(TString::Format("averagecharge_vs_pedcharge_prevsca_sca%i",isca),TString::Format("averagecharge_vs_pedcharge_prevsca_sca%i",isca),20,0,10,80,-1,1);
    nhits_vs_pedcharge_prevsca[isca] = new TH2F(TString::Format("nhits_vs_pedcharge_prevsca_sca%i",isca),TString::Format("nhits_vs_pedcharge_prevsca_sca%i",isca),20,0.5,20.5,80,-1,1);

    nchannelsElarger5_vs_pedcharge_prevsca[isca] = new TH2F(TString::Format("nchannelsElarger5_vs_pedcharge_prevsca_sca%i",isca),TString::Format("nchannelsElarger5_vs_pedcharge_prevsca_sca%i",isca),65,-0.5,64.5,80,-1,1);
    nchannelsElarger2_vs_pedcharge_prevsca[isca] = new TH2F(TString::Format("nchannelsElarger2_vs_pedcharge_prevsca_sca%i",isca),TString::Format("nchannelsElarger2_vs_pedcharge_prevsca_sca%i",isca),65,-0.5,64.5,80,-1,1);
    nchannelsElarger1_vs_pedcharge_prevsca[isca] = new TH2F(TString::Format("nchannelsElarger1_vs_pedcharge_prevsca_sca%i",isca),TString::Format("nchannelsElarger1_vs_pedcharge_prevsca_sca%i",isca),65,-0.5,64.5,80,-1,1);
    nchannelsElarger05_vs_pedcharge_prevsca[isca] = new TH2F(TString::Format("nchannelsElarger05_vs_pedcharge_prevsca_sca%i",isca),TString::Format("nchannelsElarger05_vs_pedcharge_prevsca_sca%i",isca),65,-0.5,64.5,80,-1,1);

    
  }
  
  
  TH1F *pull_beamspot[15][7];
  TH1F *stdev_beamspot[15][7];
  TH1F *mean_beamspot[15][7];
  TH1F *mean_bcid_beamspot[15][7];

  for(int isca=0; isca<15; isca++) {
    for(int ilayer=0; ilayer<7; ilayer++) {
      pull_beamspot[isca][ilayer] = new TH1F(TString::Format("pull_sca%i_layer%i_bs",isca,ilayer),TString::Format("pull_sca%i_layer%i_bs",isca,ilayer),21,-1.05,1.05);
      stdev_beamspot[isca][ilayer] = new TH1F(TString::Format("stdev_sca%i_layer%i_bs",isca,ilayer),TString::Format("stdev_sca%i_layer%i_bs",isca,ilayer),21,-1.05,1.05);
      mean_beamspot[isca][ilayer] = new TH1F(TString::Format("mean_sca%i_layer%i_bs",isca,ilayer),TString::Format("mean_sca%i_layer%i_bs",isca,ilayer),21,-1.05,1.05);
      mean_bcid_beamspot[isca][ilayer] = new TH1F(TString::Format("mean_bcid%i_layer%i_bs",isca,ilayer),TString::Format("mean_bcid%i_layer%i_bs",isca,ilayer),21,-1.05,1.05);
    }
  }

  TH1F *pull_outbeamspot[15][7];
  TH1F *stdev_outbeamspot[15][7];
  TH1F *mean_outbeamspot[15][7];
  TH1F *mean_bcid_outbeamspot[15][7];

  for(int isca=0; isca<15; isca++) {
    for(int ilayer=0; ilayer<7; ilayer++) {
      pull_outbeamspot[isca][ilayer] = new TH1F(TString::Format("pull_sca%i_layer%i",isca,ilayer),TString::Format("pull_sca%i_layer%i",isca,ilayer),21,-1.05,1.05);
      stdev_outbeamspot[isca][ilayer] = new TH1F(TString::Format("stdev_sca%i_layer%i",isca,ilayer),TString::Format("stdev_sca%i_layer%i",isca,ilayer),21,-1.05,1.05);
      mean_outbeamspot[isca][ilayer] = new TH1F(TString::Format("mean_sca%i_layer%i",isca,ilayer),TString::Format("mean_sca%i_layer%i",isca,ilayer),21,-1.05,1.05);
      mean_bcid_outbeamspot[isca][ilayer] = new TH1F(TString::Format("mean_bcid_sca%i_layer%i",isca,ilayer),TString::Format("mean_bcid_sca%i_layer%i",isca,ilayer),21,-1.05,1.05);
    }
  }

  
  TH1F* pedestal_closefromhit[15];
  TH1F* pedestal_farfromhit[15];
  for(int isca=0; isca<15; isca++) {
    pedestal_closefromhit[isca]   = new TH1F(TString::Format("pedestal_closefromhit_sca%i",isca),TString::Format("pedestal_closefromhit_sca%i",isca),80,-1,1);
    pedestal_farfromhit[isca]   = new TH1F(TString::Format("pedestal_farfromhit_sca%i",isca),TString::Format("pedestal_farfromhit_sca%i",isca),80,-1,1);
  }
  TH2F* mean_vs_entries = new TH2F("pedestal_mean_vs_entries","pedestal_mean_vs_entries",2001,-10.005,10.005,500,5,5005);
  TH2F* rms_vs_entries = new TH2F("pedestal_rms_vs_entries","pedestal_rms_vs_entries",2001,-10.005,10.005,500,5,5005);
  TH2F* mean_vs_rms = new TH2F("pedestal_mean_vs_rms","pedestal_mean_vs_rms",2001,-10.005,10.005,2001,-10.005,10.005);

  TH2F* pull_map_xy = new TH2F("pedestal_pull_map_xy","pedestal_pull_map_xy",32,-90,90,32,-90,90);
  TH2F* pull_map_xz = new TH2F("pedestal_pull_map_xz","pedestal_pull_map_xz",32,-90,90,14,-7.5,6.5);

  TH2F* pull_map_bad_xy = new TH2F("pedestal_pull_map_bad_xy","pedestal_pull_map_bad_xy",32,-90,90,32,-90,90);
  TH2F* pull_map_bad_xz = new TH2F("pedestal_pull_map_bad_xz","pedestal_pull_map_bad_xz",32,-90,90,14,-7.5,6.5);
  
  TH1F *wrongsca[16];
  TH1F *goodsca[16]; 
  TH1F *not_recalculated_sca[16];
  for(int ichip=0; ichip<16; ichip++) {
    wrongsca[ichip] = new TH1F(TString::Format("pedestal_wrong_sca_chip%i",ichip),TString::Format("pedestal_wrong_sca_chip%i",ichip),15,-0.5,14.5);
    goodsca[ichip] = new TH1F(TString::Format("pedestal_good_sca_chip%i",ichip),TString::Format("pedestal_good_sca_chip%i",ichip),15,-0.5,14.5);
    not_recalculated_sca[ichip] = new TH1F(TString::Format("pedestal_not_recalculated_sca_chip%i",ichip),TString::Format("pedestal_not_recalculated_sca_chip%i",ichip),15,-0.5,14.5);
  }
  
  TH1F* pedestal_histo[7][16][64][15];
  TH1F* charge_histo[7][16][64][15];
  TH1F* pedestal_histo_bcid[7][16][64][15];

  Double_t pedestal[7][16][64][15];
  Double_t pedestal_bcid[7][16][64][15];

  for(int iz=0; iz<7; iz++) {
    for(int ichip=0; ichip<16; ichip++) {
      for(int ichan=0; ichan<64; ichan++) {
  	for(int isca=0; isca<15;isca++) {
  	  pedestal_histo[iz][ichip][ichan][isca] = new TH1F(TString::Format("pedestal_z%i_chip%i_chan%i_sca%i",iz,ichip,ichan,isca), TString::Format("pedestal_z%i_chip%i_chan%i_sca%i",iz,ichip,ichan,isca),2000,-10,10);
	  charge_histo[iz][ichip][ichan][isca] = new TH1F(TString::Format("charge_z%i_chip%i_chan%i_sca%i",iz,ichip,ichan,isca), TString::Format("charge_z%i_chip%i_chan%i_sca%i",iz,ichip,ichan,isca),2000,-10,10);
  	  pedestal[iz][ichip][ichan][isca] = -9999.;

	  pedestal_histo_bcid[iz][ichip][ichan][isca] = new TH1F(TString::Format("pedestal_z%i_chip%i_chan%i_bcid%i",iz,ichip,ichan,isca), TString::Format("pedestal_z%i_chip%i_chan%i_bcid%i",iz,ichip,ichan,isca),2000,-10,10);
  	  pedestal_bcid[iz][ichip][ichan][isca] = -9999.;
  	}
      }
    }
  }

  TH2F* pedestal_vs_prevbcid[7][16][15];
    TH2F* pedestal_vs_nextbcid[7][16][15];

  for(int iz=0; iz<7; iz++) {
    for(int ichip=0; ichip<16; ichip++) {
      for(int isca=0; isca<15;isca++) {
	pedestal_vs_prevbcid[iz][ichip][isca] = new TH2F(TString::Format("pedestal_vs_prevbcid_z%i_chip%i_sca%i",iz,ichip,isca),TString::Format("pedestal_vs_prevbcid_z%i_chip%i_sca%i",iz,ichip,isca),210,-1.05,1.05,50,-10,990);

	pedestal_vs_nextbcid[iz][ichip][isca] = new TH2F(TString::Format("pedestal_vs_nextbcid_z%i_chip%i_sca%i",iz,ichip,isca),TString::Format("pedestal_vs_nextbcid_z%i_chip%i_sca%i",iz,ichip,isca),210,-1.05,1.05,50,-10,990);

      }
    }
  }

  double prev_spill=-1;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;

    if(ShowerBasicSelection(mipcut,bcid_max,nslabs_selection) == false) continue;

    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      int z=hit_z[ihit];
      if(z==9) z=6;
      if( IsPedestal(ihit)==true)  {
	pedestal_histo[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]->Fill(hit_energy[ihit]);
	int bcid_range=(bcid-1200)/100;
	if(bcid_range<15 && bcid>1200) pedestal_histo_bcid[z][hit_chip[ihit]][hit_chan[ihit]][bcid_range]->Fill(hit_energy[ihit]);
	
	if(prev_spill==spill ) {
	  pedestal_vs_prevbcid[z][hit_chip[ihit]][hit_sca[ihit]]->Fill(hit_energy[ihit],bcid-prev_bcid);
	  pedestal_vs_nextbcid[z][hit_chip[ihit]][hit_sca[ihit]]->Fill(hit_energy[ihit],next_bcid-bcid);
	}
      }
    }

    //  prev_bcid=bcid;
    prev_spill=spill;


  }

  cout<<"save file: "<<TString::Format("%s/Pedestals_%s_%s_%s_mipcut%0.1f_showers.txt",folder.Data(),configuration.Data(),gridpoint.Data(),energy_string.Data(),mipcut)<<endl;

  
  ofstream fout_ped(TString::Format("%s/Pedestals_%s_%s_%s_mipcut%0.1f_showers.txt",folder.Data(),configuration.Data(),gridpoint.Data(),energy_string.Data(),mipcut),ios::out);
  fout_ped<<"#pedestal results (fit to a gaussian) remove channels/sca with two pedestals peaks from the analysis"<<endl;
  fout_ped<<"#z chip channel ped0 eped0 widthped0 ped1 eped1 widthped1... ped14 eped14 widhtped14 (all SCA)"<<endl;

  /// Calculate pedestals
  // do pedestal (chip/channel/sca based) analysis
  for(int iz=0; iz<7; iz++) {
    for(int ichip=0; ichip<16; ichip++) {
      for(int ichn=0; ichn<64; ichn++) {
      	fout_ped << iz <<" "<< ichip <<" " <<ichn<< " ";
  	//	cout<< iz <<" "<< ichip <<" " <<ichn<< " "; 

  	for(int isca=0; isca<15; isca++) {
	  double entries_v = pedestal_histo[iz][ichip][ichn][isca]->GetEntries();
  	  if( entries_v> 100 ){ 
	    // pedestal_histo[iz][ichip][ichn][isca]->GetXaxis()->SetRangeUser(-1,1);
	    pedestal[iz][ichip][ichn][isca]=pedestal_histo[iz][ichip][ichn][isca]->GetMean();
	    double mean_v = pedestal_histo[iz][ichip][ichn][isca]->GetMean();
	    double rms_v = pedestal_histo[iz][ichip][ichn][isca]->GetRMS();
	    double error_v = rms_v/sqrt(entries_v);
	    double pull_v=mean_v/rms_v;

	    if(  (map_pointX[ichip][ichn]<-20 & map_pointX[ichip][ichn]>-50 ) && (map_pointY[ichip][ichn]<-30 && map_pointY[ichip][ichn]>-60 ) ) {
	      pull_beamspot[isca][iz]->Fill(pull_v);
	      stdev_beamspot[isca][iz]->Fill(rms_v);
	      mean_beamspot[isca][iz]->Fill(pedestal_histo[iz][ichip][ichn][isca]->GetMean());
	      mean_vs_entries->Fill(mean_v,entries_v);
	      rms_vs_entries->Fill(rms_v,entries_v);
	      mean_vs_rms->Fill(mean_v,rms_v);
	    } else {
	      pull_outbeamspot[isca][iz]->Fill(pull_v);
	      stdev_outbeamspot[isca][iz]->Fill(rms_v);
	      mean_outbeamspot[isca][iz]->Fill(pedestal_histo[iz][ichip][ichn][isca]->GetMean());
	    }
	    
	    // pedestal_histo[iz][ichip][ichn][isca]->GetXaxis()->SetRangeUser(-10,10);

	    if(fabs(mean_v)<0.05) {
	      goodsca[ichip]->Fill(isca); 
	      pull_map_xy->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn]);
	      pull_map_xz->Fill(map_pointX[ichip][ichn],iz);
	    }
	    
	    if(fabs(mean_v)>0.05) {
	      wrongsca[ichip]->Fill(isca);
	      pull_map_bad_xy->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn]);
	      pull_map_bad_xz->Fill(map_pointX[ichip][ichn],iz);
	      if(fabs(mean_v)>1.0)cout<<iz<<" "<<ichip<< " " <<ichn<<" " <<isca <<endl;
	    }

	    fout_ped<<mean_v<< " " << error_v<<" "<<rms_v<<" ";

  	  } else {
	    not_recalculated_sca[ichip]->Fill(isca); 
	    fout_ped<<-200<< " " << -200<<" "<<-200<<" ";
  	  }
	  fout_ped<<endl;

  	}//end for sca
	
	for(int ibcid=0; ibcid<15; ibcid++) {
	  double entries_v = pedestal_histo_bcid[iz][ichip][ichn][ibcid]->GetEntries();
	  
	  if( entries_v> 50 ){ 
	    // pedestal_histo[iz][ichip][ichn][ibcid]->GetXaxis()->SetRangeUser(-1,1);
	    pedestal_bcid[iz][ichip][ichn][ibcid]=pedestal_histo[iz][ichip][ichn][ibcid]->GetMean();
	    
	    if(  (map_pointX[ichip][ichn]<-20 & map_pointX[ichip][ichn]>-50 ) && (map_pointY[ichip][ichn]<-30 && map_pointY[ichip][ichn]>-60 ) ) {
	      mean_bcid_beamspot[ibcid][iz]->Fill(pedestal_histo_bcid[iz][ichip][ichn][ibcid]->GetMean());
	    } else {
	      mean_bcid_outbeamspot[ibcid][iz]->Fill(pedestal_histo[iz][ichip][ichn][ibcid]->GetMean());
	    }
	  }//end for bcid range
	}
      }// end if savepedestal
    }
  }

  
  // -----------------------------------------------------------------------------------------------------
  //#******************************************************************************************************
  
  // Signal readout
  nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    //if(ShowerBasicSelection_pedestal(mipcut,nslabs_selection,bcid_max, pedestal) == false)  continue;

    //-------------------------------------------------
    //basic selection
    if(bcid<1250)  continue;
    if(bcid>bcid_max) continue;
    
    
    bool bool_hit_slab[7];
    bool_hit_slab[0] = false;
    bool_hit_slab[1] = false;
    bool_hit_slab[2] = false;
    bool_hit_slab[3] = false;
    bool_hit_slab[4] = false;
    bool_hit_slab[5] = false;
    bool_hit_slab[6] = false;

    bool allpedestalcalculated=true;
    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      int z=hit_z[ihit];
      if(z==9) z=6;
      if((hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]])>mipcut && IsHit(ihit)==true ) {
	bool_hit_slab[z]=true;
      }
      if(pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]<-5000) allpedestalcalculated=false;
    }

    //    if(allpedestalcalculated==false) continue;
    
    int nslabhitted =0;
    
    for(int i=0; i<7; i++ ){
      if( bool_hit_slab[i]==true) nslabhitted++;
    }
    
    if( nslabhitted < nslabs_selection) continue;
  
   
    //    if(energy_sum_tmp<0.5) continue;
    //  number_hits++;

    int nchannelswithhits_sca[15];
    int nchannelswithElarger5_sca[15];
    int nchannelswithElarger2_sca[15];
    int nchannelswithElarger1_sca[15];
    int nchannelswithElarger05_sca[15];
    double totalcharge_sca[15];
    double chargeinanotherchip_sca[15];
    
    for(int isca=0; isca< 15; isca ++) {
      nchannelswithhits_sca[isca]=0;
      totalcharge_sca[isca]=0;
      chargeinanotherchip_sca[isca]=0;
      nchannelswithElarger5_sca[isca]=0;
      nchannelswithElarger2_sca[isca]=0;
      nchannelswithElarger1_sca[isca]=0;
      nchannelswithElarger05_sca[isca]=0;
    }


    int spill0=spill;
    int bcid0=bcid;
    int sca0=-1;

    double x_hit[64];
    double y_hit[64];
    int nhitsread=0;
    
    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      if(IsHit(ihit)==true && hit_chip[ihit]==analyze_chip && hit_z[ihit]==analyze_layer) {
	nchannelswithhits_sca[hit_sca[ihit]]++;
	totalcharge_sca[hit_sca[ihit]]+=hit_energy[ihit];
	if(hit_energy[ihit]>0.5) nchannelswithElarger05_sca[hit_sca[ihit]]++;
	if(hit_energy[ihit]>5) nchannelswithElarger5_sca[hit_sca[ihit]]++;
	if(hit_energy[ihit]>2) nchannelswithElarger2_sca[hit_sca[ihit]]++;
	if(hit_energy[ihit]>1) nchannelswithElarger1_sca[hit_sca[ihit]]++;
	sca0=hit_sca[ihit];
	x_hit[nhitsread]=hit_x[ihit];
	y_hit[nhitsread]=hit_y[ihit];
	nhitsread++;
      }
      if(IsHit(ihit)==true && hit_chip[ihit]!=analyze_chip && hit_z[ihit]==analyze_layer) {
	chargeinanotherchip_sca[hit_sca[ihit]]+=hit_energy[ihit];
      }
	
    }


    double x_pedestal[64];
    double y_pedestal[64];
    double energy_pedestal[64];
    int sca_pedestal[64];
    int npedestalsread=0;
    
    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      if(IsPedestal(ihit)==true && nchannelswithhits_sca[hit_sca[ihit]]>0 && hit_chip[ihit]==analyze_chip && hit_z[ihit]==analyze_layer && sca0==hit_sca[ihit]) {
	totalcharge_vs_pedcharge[hit_sca[ihit]]->Fill(totalcharge_sca[hit_sca[ihit]],hit_energy[ihit]);
	averagecharge_vs_pedcharge[hit_sca[ihit]]->Fill(totalcharge_sca[hit_sca[ihit]]/nchannelswithhits_sca[hit_sca[ihit]],hit_energy[ihit]);
	nhits_vs_pedcharge[hit_sca[ihit]]->Fill(nchannelswithhits_sca[hit_sca[ihit]],hit_energy[ihit]);

	nchannelsElarger5_vs_pedcharge[hit_sca[ihit]]->Fill(nchannelswithElarger5_sca[hit_sca[ihit]],hit_energy[ihit]);
	nchannelsElarger2_vs_pedcharge[hit_sca[ihit]]->Fill(nchannelswithElarger2_sca[hit_sca[ihit]],hit_energy[ihit]);
	nchannelsElarger1_vs_pedcharge[hit_sca[ihit]]->Fill(nchannelswithElarger1_sca[hit_sca[ihit]],hit_energy[ihit]);
	nchannelsElarger05_vs_pedcharge[hit_sca[ihit]]->Fill(nchannelswithElarger05_sca[hit_sca[ihit]],hit_energy[ihit]);
	x_pedestal[npedestalsread]=hit_x[ihit];
	y_pedestal[npedestalsread]=hit_y[ihit];
	energy_pedestal[npedestalsread]=hit_energy[ihit];
	sca_pedestal[npedestalsread]=hit_sca[ihit];
	npedestalsread++;
      }

      if(IsPedestal(ihit)==true && nchannelswithhits_sca[hit_sca[ihit]]>0 && hit_chip[ihit]!=analyze_chip && hit_z[ihit]==analyze_layer && sca0==hit_sca[ihit]) {
	chargeinanotherchip_vs_pedcharge[hit_sca[ihit]]->Fill(chargeinanotherchip_sca[hit_sca[ihit]],hit_energy[ihit]);
      }
    }

    for(int iped=0; iped< npedestalsread; iped ++) {
      double mindist=1000;
      
      for(int ihit=0; ihit< nhitsread; ihit ++) {
	double dist=sqrt(TMath::Power(x_pedestal[iped]-x_hit[ihit],2) + TMath::Power(y_pedestal[iped]-y_hit[ihit],2));
	if(dist<mindist) mindist=dist;
      }
      if(mindist <= 2*sqrt(50) ) pedestal_closefromhit[sca_pedestal[iped]]->Fill(energy_pedestal[iped]);
      if(mindist > 5*sqrt(50) ) pedestal_farfromhit[sca_pedestal[iped]]->Fill(energy_pedestal[iped]);
    }
	    

    
    
    int nchannelswithhits_prevsca=0;
    int nchannelswithElarger5_prevsca=0;
    int nchannelswithElarger2_prevsca=0;
    int nchannelswithElarger1_prevsca=0;
    int nchannelswithElarger05_prevsca=0;
    double totalcharge_prevsca=0;
    double chargeinanotherchip_prevsca=0;
    
     bool exitloop=false;

    if(sca0>0 && jentry>101 && nchannelswithhits_sca[sca0]>0) {
      for (Long64_t jentry2=jentry-1; jentry2>jentry-500;jentry2--) {
	Long64_t ientry2 = LoadTree(jentry2);
	if (ientry2 < 0) break;
	Long64_t nbytes2 = 0, nb2 = 0;
	nb2 = fChain->GetEntry(jentry2);   nbytes2 += nb2;
	if(spill<spill0) continue;
	if(exitloop==true) continue;
	if(bcid<bcid0) {
	   bool exitloop0=false;
	  for(int ihit2=0; ihit2< nhit_chan; ihit2 ++) {
	    if(IsHit(ihit2)==true && hit_chip[ihit2]==analyze_chip && hit_z[ihit2]==analyze_layer && (sca0-hit_sca[ihit2])==1) {
	      nchannelswithhits_prevsca++;
	      totalcharge_prevsca+=hit_energy[ihit2];
	      
	      if(hit_energy[ihit2]>0.5) nchannelswithElarger05_prevsca++;
	      if(hit_energy[ihit2]>5) nchannelswithElarger5_prevsca++;
	      if(hit_energy[ihit2]>2) nchannelswithElarger2_prevsca++;
	      if(hit_energy[ihit2]>1) nchannelswithElarger1_prevsca++;
	       exitloop0=true;
	    }
	    if(IsHit(ihit2)==true && hit_chip[ihit2]!=analyze_chip && hit_z[ihit2]==analyze_layer && (sca0-hit_sca[ihit2])==1) {
	      chargeinanotherchip_prevsca+=hit_energy[ihit2];
	    }
	  }
	  if(exitloop0==true) exitloop=true;
	}
      }
    }

    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(sca0>0 && jentry>101 && nchannelswithhits_prevsca>0) {

      for(int ihit=0; ihit< nhit_chan; ihit ++) {
	if(IsPedestal(ihit)==true && hit_chip[ihit]==analyze_chip && hit_z[ihit]==analyze_layer && sca0==hit_sca[ihit]) {
	  totalcharge_vs_pedcharge_prevsca[hit_sca[ihit]]->Fill(totalcharge_prevsca,hit_energy[ihit]);
	  averagecharge_vs_pedcharge_prevsca[hit_sca[ihit]]->Fill(totalcharge_prevsca/nchannelswithhits_prevsca,hit_energy[ihit]);
	  nhits_vs_pedcharge_prevsca[hit_sca[ihit]]->Fill(nchannelswithhits_prevsca,hit_energy[ihit]);
	  
	  nchannelsElarger5_vs_pedcharge_prevsca[hit_sca[ihit]]->Fill(nchannelswithElarger5_prevsca,hit_energy[ihit]);
	  nchannelsElarger2_vs_pedcharge_prevsca[hit_sca[ihit]]->Fill(nchannelswithElarger2_prevsca,hit_energy[ihit]);
	  nchannelsElarger1_vs_pedcharge_prevsca[hit_sca[ihit]]->Fill(nchannelswithElarger1_prevsca,hit_energy[ihit]);
	  nchannelsElarger05_vs_pedcharge_prevsca[hit_sca[ihit]]->Fill(nchannelswithElarger05_prevsca,hit_energy[ihit]);
	}
	if(IsPedestal(ihit)==true && hit_chip[ihit]!=analyze_chip && hit_z[ihit]==analyze_layer && sca0==hit_sca[ihit]) {
	  chargeinanotherchip_vs_pedcharge_prevsca[hit_sca[ihit]]->Fill(chargeinanotherchip_prevsca,hit_energy[ihit]);
	}
      }
    }
    
    
  }//if jentry

  
  // Signal analysis
  TFile *file = new TFile(TString::Format("%s%s_%s_%s_mipcut%0.1f_showers.root",folder.Data(),configuration.Data(),gridpoint.Data(),energy_string.Data(),mipcut) , "RECREATE");
  file->cd();

  for(int iz=0; iz<7; iz++) {
    for(int ichip=0; ichip<16; ichip++) {
      for(int isca=0; isca<15; isca++){
	pedestal_vs_prevbcid[iz][ichip][isca]->Write();
	pedestal_vs_nextbcid[iz][ichip][isca]->Write();
      }
    }
  }

  for(int isca=0; isca<15; isca++) {
    for(int iz=0; iz<7; iz++) {
      pull_beamspot[isca][iz]->SetName(TString::Format("bs_pedestal_pull_distribution_isca%i_ilayer%i",isca,iz));
      pull_beamspot[isca][iz]->Write();
      stdev_beamspot[isca][iz]->SetName(TString::Format("bs_pedestal_stdev_distribution_isca%i_ilayer%i",isca,iz));
      stdev_beamspot[isca][iz]->Write();
      mean_beamspot[isca][iz]->SetName(TString::Format("bs_pedestal_mean_distribution_isca%i_ilayer%i",isca,iz));
      mean_beamspot[isca][iz]->Write();
      mean_bcid_beamspot[isca][iz]->SetName(TString::Format("bs_pedestal_mean_distribution_ibcid%i_ilayer%i",isca,iz));
      mean_bcid_beamspot[isca][iz]->Write();
      
      pull_outbeamspot[isca][iz]->SetName(TString::Format("pedestal_pull_distribution_isca%i_ilayer%i",isca,iz));
      pull_outbeamspot[isca][iz]->Write();
      stdev_outbeamspot[isca][iz]->SetName(TString::Format("pedestal_stdev_distribution_isca%i_ilayer%i",isca,iz));
      stdev_outbeamspot[isca][iz]->Write();
      mean_outbeamspot[isca][iz]->SetName(TString::Format("pedestal_mean_distribution_isca%i_ilayer%i",isca,iz));
      mean_outbeamspot[isca][iz]->Write();
      mean_bcid_outbeamspot[isca][iz]->SetName(TString::Format("pedestal_mean_distribution_ibcid%i_ilayer%i",isca,iz));
      mean_bcid_outbeamspot[isca][iz]->Write();
    }
  }

  for(int isca=0; isca<15; isca++) {
    pedestal_closefromhit[isca]->Write();
    pedestal_farfromhit[isca]->Write();
  }
  mean_vs_entries->Write();
  rms_vs_entries->Write();
  mean_vs_rms->Write();


  pull_map_xy->Write();
  pull_map_xz->Write();

  pull_map_bad_xy->Write();
  pull_map_bad_xz->Write();
  
  for(int ichip=0; ichip<16; ichip++) {
    not_recalculated_sca[ichip]->Write();
    goodsca[ichip]->Write();
    wrongsca[ichip]->Write();
  }

  
  for(int isca=0; isca<15; isca++) {

    //standard maps
    totalcharge_vs_pedcharge[isca]->Write();
    chargeinanotherchip_vs_pedcharge[isca]->Write();
    
    averagecharge_vs_pedcharge[isca]->Write();
    nhits_vs_pedcharge[isca]->Write();
    
    nchannelsElarger5_vs_pedcharge[isca]->Write();
    nchannelsElarger2_vs_pedcharge[isca]->Write();
    nchannelsElarger1_vs_pedcharge[isca]->Write();
    nchannelsElarger05_vs_pedcharge[isca]->Write();
    /* // normalized plots
    TH2F *totalcharge_vs_pedcharge_norm_temp = (TH2F*)totalcharge_vs_pedcharge[isca]->Clone("totalcharge_vs_pedcharge_norm_temp");
    TH2F *averagecharge_vs_pedcharge_norm_temp = (TH2F*)averagecharge_vs_pedcharge[isca]->Clone("averagecharge_vs_pedcharge_norm_temp");
    TH2F *nhits_vs_pedcharge_norm_temp = (TH2F*)nhits_vs_pedcharge[isca]->Clone("nhits_vs_pedcharge_norm_temp");

    TH2F *nchannelsElarger5_vs_pedcharge_norm_temp = (TH2F*)nchannelsElarger5_vs_pedcharge[isca]->Clone("nchannelsElarger5_vs_pedcharge_norm_temp");
    TH2F *nchannelsElarger2_vs_pedcharge_norm_temp = (TH2F*)nchannelsElarger2_vs_pedcharge[isca]->Clone("nchannelsElarger2_vs_pedcharge_norm_temp");
    TH2F *nchannelsElarger1_vs_pedcharge_norm_temp = (TH2F*)nchannelsElarger1_vs_pedcharge[isca]->Clone("nchannelsElarger1_vs_pedcharge_norm_temp");
    TH2F *nchannelsElarger05_vs_pedcharge_norm_temp = (TH2F*)nchannelsElarger05_vs_pedcharge[isca]->Clone("nchannelsElarger05_vs_pedcharge_norm_temp");

    for(int ix=0; ix<totalcharge_vs_pedcharge_norm_temp->GetNbinsX(); ix++) {
      double entries_x=0;
      for(int iy=0; iy<totalcharge_vs_pedcharge_norm_temp->GetNbinsY(); iy++) 
	entries_x+=totalcharge_vs_pedcharge_norm_temp->GetBinContent(ix,iy);

      for(int iy=0; iy<totalcharge_vs_pedcharge_norm_temp->GetNbinsY(); iy++) 
	if(entries_x>0) totalcharge_vs_pedcharge_norm_temp->SetBinContent(ix,iy,totalcharge_vs_pedcharge_norm_temp->GetBinContent(ix,iy)/entries_x);
    }

    totalcharge_vs_pedcharge_norm_temp->SetName(TString::Format("totalcharge_vs_pedcharge_sca%i_norm",isca));
    totalcharge_vs_pedcharge_norm_temp->Write();*/

    //---------------------------------------------------------------------------------
    // prevsca plots

    totalcharge_vs_pedcharge_prevsca[isca]->Write();
    chargeinanotherchip_vs_pedcharge_prevsca[isca]->Write();
    
    averagecharge_vs_pedcharge_prevsca[isca]->Write();
    nhits_vs_pedcharge_prevsca[isca]->Write();
    
    nchannelsElarger5_vs_pedcharge_prevsca[isca]->Write();
    nchannelsElarger2_vs_pedcharge_prevsca[isca]->Write();
    nchannelsElarger1_vs_pedcharge_prevsca[isca]->Write();
    nchannelsElarger05_vs_pedcharge_prevsca[isca]->Write();

  }

  
  TDirectory *cdhisto = file->mkdir("histograms");
  cdhisto->cd();
 
  for(int iz=0; iz<7; iz++) {
    for(int ichip=0; ichip<16; ichip++) {
      for(int ichn=0; ichn<64; ichn++) {
	for(int isca=0; isca<15; isca++){
	  if( pedestal_histo[iz][ichip][ichn][isca]->GetEntries()>50) pedestal_histo[iz][ichip][ichn][isca]->Write();
	  charge_histo[iz][ichip][ichn][isca]->Write();
	}
      }
    }
  }

  file->Close();
  
  
}


void protoAnalysis::ShowerDistributions_pedestal(TString folder="y", TString configuration="conf1", TString energy_string="3GeV", TString gridpoint ="grid20", double mipcut=0.5, int bcid_max=2900, int nslabs_selection=6, bool RMIsolatedHits=true)
{

  ReadMap("../fev10_chip_channel_x_y_mapping.txt");
  ReadCalibrated("../mip_calib/MIP_");

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
  TH1F* energy_center_05 = new TH1F("energy_center_zm05","energy_center_zm05",500,1,1001);
  TH1F* energy_center_1 = new TH1F("energy_center_zm1","energy_center_zm1",500,1,1001);
  TH1F* energy_center_2 = new TH1F("energy_center_zm2","energy_center_zm2",500,1,1001);
    
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
  TH1F* nhitstotal_center_05 = new TH1F("nhitstotal_center_zm05","nhitstotal_center_zm05",500,1,1001);
  TH1F* nhitstotal_center_1 = new TH1F("nhitstotal_center_zm1","nhitstotal_center_zm1",500,1,1001);
  TH1F* nhitstotal_center_2 = new TH1F("nhitstotal_center_zm2","nhitstotal_center_zm2",500,1,1001);
      
  TH2F* nhits_vs_energy= new TH2F("nhits_vs_energy","nhits_vs_energy",51,-2.5,252.5,101,-5,1005);
  TH2F* min_dist_vs_energy= new TH2F("min_dist_vs_energy","min_dist_vs_energy",250,0,250,100,-5,995);

  int number_hits=0;
  // -----------------------------------------------------------------------------------------------------

  // -----------------------------------------------------------------------------------------------------
  //#******************************************************************************************************
  // PEDESTAL PART:
  // th1f, i,j,k,sca
  TH2F* totalcharge_vs_pedcharge[15];
  TH2F* averagecharge_vs_pedcharge[15];
  TH2F* nhits_vs_pedcharge[15];
  TH2F* nchannelsElarger5_vs_pedcharge[15];
  TH2F* nchannelsElarger10_vs_pedcharge[15];
  TH2F* nchannelsElarger20_vs_pedcharge[15];
  TH2F* nchannelsEsmaller05_vs_pedcharge[15];

  TH2F* totalcharge_vs_pedcharge_prevsca[14];
  TH2F* averagecharge_vs_pedcharge_prevsca[14];
  TH2F* nhits_vs_pedcharge_prevsca[14];
  TH2F* nchannelsElarger5_vs_pedcharge_prevsca[15];
  TH2F* nchannelsElarger10_vs_pedcharge_prevsca[15];
  TH2F* nchannelsElarger20_vs_pedcharge_prevsca[15];
  TH2F* nchannelsEsmaller05_vs_pedcharge_prevsca[15];

  for(int isca=0; isca<15; isca++) {
    totalcharge_vs_pedcharge[isca] = new TH2F(TString::Format("totalcharge_vs_pedcharge_sca%i",isca),TString::Format("totalcharge_vs_pedcharge_sca%i",isca),100,5,1005,1500,-4.995,10.005);
    averagecharge_vs_pedcharge[isca] = new TH2F(TString::Format("averagecharge_vs_pedcharge_sca%i",isca),TString::Format("averagecharge_vs_pedcharge_sca%i",isca),300,0,30,1500,-4.995,10.005);
    nhits_vs_pedcharge[isca] = new TH2F(TString::Format("nhits_vs_pedcharge_sca%i",isca),TString::Format("nhits_vs_pedcharge_sca%i",isca),200,0.5,200.5,1500,-4.995,10.005);

    nchannelsElarger5_vs_pedcharge[isca] = new TH2F(TString::Format("nchannelsElarger5_vs_pedcharge_sca%i",isca),TString::Format("nchannelsElarger5_vs_pedcharge_sca%i",isca),200,0.5,200.5,1500,-4.995,10.005);
    nchannelsElarger10_vs_pedcharge[isca] = new TH2F(TString::Format("nchannelsElarger10_vs_pedcharge_sca%i",isca),TString::Format("nchannelsElarger10_vs_pedcharge_sca%i",isca),200,0.5,200.5,1500,-4.995,10.005);
    nchannelsElarger20_vs_pedcharge[isca] = new TH2F(TString::Format("nchannelsElarger20_vs_pedcharge_sca%i",isca),TString::Format("nchannelsElarger20_vs_pedcharge_sca%i",isca),200,0.5,200.5,1500,-4.995,10.005);
    nchannelsEsmaller05_vs_pedcharge[isca] = new TH2F(TString::Format("nchannelsEsmaller05_vs_pedcharge_sca%i",isca),TString::Format("nchannelsEsmaller05_vs_pedcharge_sca%i",isca),200,0.5,200.5,1500,-4.995,10.005);

  }

  for(int isca=0; isca<14; isca++) {
    totalcharge_vs_pedcharge_prevsca[isca] = new TH2F(TString::Format("totalcharge_vs_pedcharge_prevsca_sca%i",isca+1),TString::Format("totalcharge_vs_pedcharge_prevsca_sca%i",isca+1),100,5,1005,1500,-4.995,10.005);
    averagecharge_vs_pedcharge_prevsca[isca] = new TH2F(TString::Format("averagecharge_vs_pedcharge_prevsca_sca%i",isca+1),TString::Format("averagecharge_vs_pedcharge_prevsca_sca%i",isca+1),150,0,30,1500,-4.995,10.005);
    nhits_vs_pedcharge_prevsca[isca] = new TH2F(TString::Format("nhits_vs_pedcharge_prevsca_sca%i",isca+1),TString::Format("nhits_vs_pedcharge_prevsca_sca%i",isca+1),200,0.5,200.5,1500,-4.995,10.005);

    nchannelsElarger5_vs_pedcharge_prevsca[isca] = new TH2F(TString::Format("nchannelsElarger5_vs_pedcharge_prevsca_sca%i",isca+1),TString::Format("nchannelsElarger5_vs_pedcharge_prevsca_sca%i",isca+1),200,0.5,200.5,1500,-4.995,10.005);
    nchannelsElarger10_vs_pedcharge_prevsca[isca] = new TH2F(TString::Format("nchannelsElarger10_vs_pedcharge_prevsca_sca%i",isca+1),TString::Format("nchannelsElarger10_vs_pedcharge_prevsca_sca%i",isca+1),200,0.5,200.5,1500,-4.995,10.005);
    nchannelsElarger20_vs_pedcharge_prevsca[isca] = new TH2F(TString::Format("nchannelsElarger20_vs_pedcharge_prevsca_sca%i",isca+1),TString::Format("nchannelsElarger20_vs_pedcharge_prevsca_sca%i",isca+1),200,0.5,200.5,1500,-4.995,10.005);
    nchannelsEsmaller05_vs_pedcharge_prevsca[isca] = new TH2F(TString::Format("nchannelsEsmaller05_vs_pedcharge_prevsca_sca%i",isca+1),TString::Format("nchannelsEsmaller05_vs_pedcharge_prevsca_sca%i",isca+1),200,0.5,200.5,1500,-4.995,10.005);

    
  }


  
  
  TH1F *pull_beamspot[15];
  TH1F *stdev_beamspot[15];
  TH1F *mean_beamspot[15];

  for(int isca=0; isca<15; isca++) {
    pull_beamspot[isca] = new TH1F(TString::Format("pull_sca%i_bs",isca),TString::Format("pull_sca%i_bs",isca),201,-10.05,10.05);
    stdev_beamspot[isca] = new TH1F(TString::Format("stdev_sca%i_bs",isca),TString::Format("stdev_sca%i_bs",isca),2001,-10.005,10.005);
    mean_beamspot[isca] = new TH1F(TString::Format("mean_sca%i_bs",isca),TString::Format("mean_sca%i_bs",isca),2001,-10.005,10.005);
  }

  TH1F *pull_outbeamspot[15];
  TH1F *stdev_outbeamspot[15];
  TH1F *mean_outbeamspot[15];

  for(int isca=0; isca<15; isca++) {
    pull_outbeamspot[isca] = new TH1F(TString::Format("pull_sca%i",isca),TString::Format("pull_sca%i",isca),201,-10.05,10.05);
    stdev_outbeamspot[isca] = new TH1F(TString::Format("stdev_sca%i",isca),TString::Format("stdev_sca%i",isca),2001,-10.005,10.005);
    mean_outbeamspot[isca] = new TH1F(TString::Format("mean_sca%i",isca),TString::Format("mean_sca%i",isca),2001,-10.005,10.005);
  }
  
  
  TH2F* mean_vs_entries = new TH2F("pedestal_mean_vs_entries","pedestal_mean_vs_entries",2001,-10.005,10.005,500,5,5005);
  TH2F* rms_vs_entries = new TH2F("pedestal_rms_vs_entries","pedestal_rms_vs_entries",2001,-10.005,10.005,500,5,5005);
  TH2F* mean_vs_rms = new TH2F("pedestal_mean_vs_rms","pedestal_mean_vs_rms",2001,-10.005,10.005,2001,-10.005,10.005);

  TH2F* pull_map_xy = new TH2F("pedestal_pull_map_xy","pedestal_pull_map_xy",32,-90,90,32,-90,90);
  TH2F* pull_map_xz = new TH2F("pedestal_pull_map_xz","pedestal_pull_map_xz",32,-90,90,14,-7.5,6.5);

  TH2F* pull_map_bad_xy = new TH2F("pedestal_pull_map_bad_xy","pedestal_pull_map_bad_xy",32,-90,90,32,-90,90);
  TH2F* pull_map_bad_xz = new TH2F("pedestal_pull_map_bad_xz","pedestal_pull_map_bad_xz",32,-90,90,14,-7.5,6.5);
  
  TH1F *wrongsca[16];
  TH1F *goodsca[16]; 
  TH1F *not_recalculated_sca[16];
  for(int ichip=0; ichip<16; ichip++) {
    wrongsca[ichip] = new TH1F(TString::Format("pedestal_wrong_sca_chip%i",ichip),TString::Format("pedestal_wrong_sca_chip%i",ichip),15,-0.5,14.5);
    goodsca[ichip] = new TH1F(TString::Format("pedestal_good_sca_chip%i",ichip),TString::Format("pedestal_good_sca_chip%i",ichip),15,-0.5,14.5);
    not_recalculated_sca[ichip] = new TH1F(TString::Format("pedestal_not_recalculated_sca_chip%i",ichip),TString::Format("pedestal_not_recalculated_sca_chip%i",ichip),15,-0.5,14.5);
  }
  
  TH1F* pedestal_histo[7][16][64][15];
  TH1F* charge_histo[7][16][64][15];

  Double_t pedestal[7][16][64][15];

  for(int iz=0; iz<7; iz++) {
    for(int ichip=0; ichip<16; ichip++) {
      for(int ichan=0; ichan<64; ichan++) {
  	for(int isca=0; isca<15;isca++) {
  	  pedestal_histo[iz][ichip][ichan][isca] = new TH1F(TString::Format("pedestal_z%i_chip%i_chan%i_sca%i",iz,ichip,ichan,isca), TString::Format("pedestal_z%i_chip%i_chan%i_sca%i",iz,ichip,ichan,isca),2000,-10,10);
	  charge_histo[iz][ichip][ichan][isca] = new TH1F(TString::Format("charge_z%i_chip%i_chan%i_sca%i",iz,ichip,ichan,isca), TString::Format("charge_z%i_chip%i_chan%i_sca%i",iz,ichip,ichan,isca),2000,-10,10);
  	  pedestal[iz][ichip][ichan][isca] = 0.;
  	}
      }
    }
  }

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;

    if(ShowerBasicSelection(mipcut,bcid_max,nslabs_selection) == false) continue;

    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      int z=hit_z[ihit];
      if(z==9) z=6;
      if( IsPedestal(ihit)==true)  {
	pedestal_histo[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]->Fill(hit_energy[ihit]);
      }
    }					       

  }

  cout<<"save file: "<<TString::Format("%s/Pedestals_%s_%s_%s_mipcut%0.1f_showers.txt",folder.Data(),configuration.Data(),gridpoint.Data(),energy_string.Data(),mipcut)<<endl;

  
  ofstream fout_ped(TString::Format("%s/Pedestals_%s_%s_%s_mipcut%0.1f_showers.txt",folder.Data(),configuration.Data(),gridpoint.Data(),energy_string.Data(),mipcut),ios::out);
  fout_ped<<"#pedestal results (fit to a gaussian) remove channels/sca with two pedestals peaks from the analysis"<<endl;
  fout_ped<<"#z chip channel ped0 eped0 widthped0 ped1 eped1 widthped1... ped14 eped14 widhtped14 (all SCA)"<<endl;

  /// Calculate pedestals
  // do pedestal (chip/channel/sca based) analysis
  for(int iz=0; iz<7; iz++) {
    for(int ichip=0; ichip<16; ichip++) {
      for(int ichn=0; ichn<64; ichn++) {
      	fout_ped << iz <<" "<< ichip <<" " <<ichn<< " ";
  	//	cout<< iz <<" "<< ichip <<" " <<ichn<< " "; 

  	for(int isca=0; isca<15; isca++) {

	  double entries_v = pedestal_histo[iz][ichip][ichn][isca]->GetEntries();

  	  if( entries_v> 50 ){ 


	    // pedestal_histo[iz][ichip][ichn][isca]->GetXaxis()->SetRangeUser(-1,1);
	    pedestal[iz][ichip][ichn][isca]=pedestal_histo[iz][ichip][ichn][isca]->GetMean();
	    double mean_v = pedestal_histo[iz][ichip][ichn][isca]->GetMean();
	    double rms_v = pedestal_histo[iz][ichip][ichn][isca]->GetRMS();
	    double error_v = rms_v/sqrt(entries_v);
	    double pull_v=mean_v/rms_v;

	    if(  (map_pointX[ichip][ichn]<-20 & map_pointX[ichip][ichn]>-50 ) && (map_pointY[ichip][ichn]<-30 && map_pointY[ichip][ichn]>-60 ) ) {
	      pull_beamspot[isca]->Fill(pull_v);
	      stdev_beamspot[isca]->Fill(rms_v);
	      mean_beamspot[isca]->Fill(pedestal_histo[iz][ichip][ichn][isca]->GetMean());
	      
	      mean_vs_entries->Fill(mean_v,entries_v);
	      rms_vs_entries->Fill(rms_v,entries_v);
	      mean_vs_rms->Fill(mean_v,rms_v);
	    } else {
	      pull_outbeamspot[isca]->Fill(pull_v);
	      stdev_outbeamspot[isca]->Fill(rms_v);
	      mean_outbeamspot[isca]->Fill(pedestal_histo[iz][ichip][ichn][isca]->GetMean());
	    }
	    
	    // pedestal_histo[iz][ichip][ichn][isca]->GetXaxis()->SetRangeUser(-10,10);

	    if(fabs(mean_v)<0.05) {
	      goodsca[ichip]->Fill(isca); 
	      pull_map_xy->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn]);
	      pull_map_xz->Fill(map_pointX[ichip][ichn],iz);
	    }
	    
	    if(fabs(mean_v)>0.05) {
	      wrongsca[ichip]->Fill(isca);
	      pull_map_bad_xy->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn]);
	      pull_map_bad_xz->Fill(map_pointX[ichip][ichn],iz);
	      if(fabs(mean_v)>1.0)cout<<iz<<" "<<ichip<< " " <<ichn<<" " <<isca <<endl;
	    }

	    fout_ped<<mean_v<< " " << error_v<<" "<<rms_v<<" ";

  	  } else {
	    not_recalculated_sca[ichip]->Fill(isca); 
	    fout_ped<<-200<< " " << -200<<" "<<-200<<" ";
  	  }	  

  	}//end for sca
  	fout_ped<<endl;
  	//cout<<endl;
      }// end if savepedestal
    }
  }

  
  

  
  // -----------------------------------------------------------------------------------------------------
  //#******************************************************************************************************
  
  // Signal readout
  nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    //if(ShowerBasicSelection_pedestal(mipcut,nslabs_selection,bcid_max, pedestal) == false)  continue;

    //-------------------------------------------------
    //basic selection
    if(bcid<1250)  continue;
    if(bcid>bcid_max) continue;
    
    
    bool bool_hit_slab[7];
    bool_hit_slab[0] = false;
    bool_hit_slab[1] = false;
    bool_hit_slab[2] = false;
    bool_hit_slab[3] = false;
    bool_hit_slab[4] = false;
    bool_hit_slab[5] = false;
    bool_hit_slab[6] = false;
    
    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      int z=hit_z[ihit];
      if(z==9) z=6;
      if((hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]])>mipcut && IsHit(ihit)==true ) {
	bool_hit_slab[z]=true;
      }
    }
    
    
    int nslabhitted =0;
    
    for(int i=0; i<7; i++ ){
      if( bool_hit_slab[i]==true) nslabhitted++;
    }
    
    if( nslabhitted < nslabs_selection) continue;
  
    //-------------------------------------------------
    

    // ---------------------------------------------
    //   the energy and other variables,
    double energy_sum_tmp=0;
    double weight_X0_sum_tmp=0;
    double energy_X0_sum_tmp=0, energy_X0_sum_tmp_dist=0;
    double xm=0, ym=0., zm=0.;

    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      int z=hit_z[ihit];
      if(z==9) z=6;
      
      energy_X0_sum_tmp_dist += w[int(hit_z[ihit])] * (hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]);
    }
    
    int n_hitstotal=0;
    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      int z=hit_z[ihit];
      if(z==9) z=6;
      
      if((hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]])>mipcut && IsHit(ihit)==true) {
	charge_histo[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]->Fill(hit_energy[ihit]);
      }

      double mindist=MinDistCalc(ihit,mipcut);
      min_dist_vs_energy->Fill(mindist,energy_X0_sum_tmp_dist);

      double mindist_th=10000.;
      if(RMIsolatedHits==true) mindist_th=7.0;

      if(mindist< mindist_th && IsHit(ihit)==true && (hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]])>mipcut) {
	n_hitstotal++;
	energy_sum_tmp += (hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]);
	energy_X0_sum_tmp += w[int(hit_z[ihit])] * (hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]);
	weight_X0_sum_tmp += w[int(hit_z[ihit])];

	xm +=  - hit_x[ihit] * (hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]) * w[int(hit_z[ihit])];
	ym +=  - hit_y[ihit] * (hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]) * w[int(hit_z[ihit])];
	zm +=  hit_z[ihit] * (hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]) * w[int(hit_z[ihit])];
      }
	
    }

    
    if(energy_sum_tmp<0.5) continue;
    number_hits++;

    nhits_vs_energy->Fill(n_hitstotal,energy_X0_sum_tmp);
    nhitstotal->Fill(n_hitstotal);


    xm = xm / energy_X0_sum_tmp;
    ym = ym / energy_X0_sum_tmp;
    zm = zm / energy_X0_sum_tmp;

    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      int z=hit_z[ihit];
      if(z==9) z=6;
      if((hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]])>mipcut && IsHit(ihit)==true) {
	double zlayer=hit_z[ihit];
	if(hit_z[ihit]>8) zlayer=6;
	energy_profile_z->Fill(zlayer,(hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]])/energy_sum_tmp);
	energy_profile_X0->Fill(hit_x0[ihit],(hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]])/energy_sum_tmp);
	
	energy_xy->Fill(-hit_x[ihit],-hit_y[ihit],(hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]])/energy_sum_tmp );
	energy_xz->Fill(-hit_x[ihit],hit_z[ihit],(hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]])/energy_sum_tmp );
	energy_yz->Fill(-hit_y[ihit],hit_z[ihit],(hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]])/energy_sum_tmp );     
	
	energy_xX0->Fill(-hit_x[ihit],hit_x0[ihit],(hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]])/energy_sum_tmp);
	energy_yX0->Fill(-hit_y[ihit],hit_x0[ihit],(hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]])/energy_sum_tmp);
      }
    }
    int nchannelswithhits=0;
    double totalcharge=0;

    int nchannelswithhits_sca[15];
    int nchannelswithElarger5_sca[15];
    int nchannelswithElarger10_sca[15];
    int nchannelswithElarger20_sca[15];
    int nchannelswithEsmaller05_sca[15];
    double totalcharge_sca[15];
    
    for(int isca=0; isca< 15; isca ++) {
      nchannelswithhits_sca[isca]=0;
      totalcharge_sca[isca]=0;
      nchannelswithElarger5_sca[isca]=0;
      nchannelswithElarger10_sca[isca]=0;
      nchannelswithElarger20_sca[isca]=0;
      nchannelswithEsmaller05_sca[isca]=0;
    }

    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      if(IsHit(ihit)==true && hit_chip[ihit]==12 && hit_z[ihit]==9) {
	nchannelswithhits++;
	totalcharge+=hit_energy[ihit];
	nchannelswithhits_sca[hit_sca[ihit]]++;
	totalcharge_sca[hit_sca[ihit]]+=hit_energy[ihit];
	if(hit_energy[ihit]<0.5) nchannelswithEsmaller05_sca[hit_sca[ihit]]++;
	if(hit_energy[ihit]>5) nchannelswithElarger5_sca[hit_sca[ihit]]++;
	if(hit_energy[ihit]>10) nchannelswithElarger10_sca[hit_sca[ihit]]++;
	if(hit_energy[ihit]>20) nchannelswithElarger20_sca[hit_sca[ihit]]++;
	
      }
    }

	
    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      if(IsPedestal(ihit)==true && nchannelswithhits_sca[hit_sca[ihit]]>0 && hit_chip[ihit]==12 && hit_z[ihit]==9) {
	totalcharge_vs_pedcharge[hit_sca[ihit]]->Fill(totalcharge_sca[hit_sca[ihit]],hit_energy[ihit]);
	averagecharge_vs_pedcharge[hit_sca[ihit]]->Fill(totalcharge_sca[hit_sca[ihit]]/nchannelswithhits_sca[hit_sca[ihit]],hit_energy[ihit]);
	nhits_vs_pedcharge[hit_sca[ihit]]->Fill(nchannelswithhits_sca[hit_sca[ihit]],hit_energy[ihit]);

	nchannelsElarger5_vs_pedcharge[hit_sca[ihit]]->Fill(nchannelswithElarger5_sca[hit_sca[ihit]],hit_energy[ihit]);
	nchannelsElarger10_vs_pedcharge[hit_sca[ihit]]->Fill(nchannelswithElarger10_sca[hit_sca[ihit]],hit_energy[ihit]);
	nchannelsElarger20_vs_pedcharge[hit_sca[ihit]]->Fill(nchannelswithElarger20_sca[hit_sca[ihit]],hit_energy[ihit]);
	nchannelsEsmaller05_vs_pedcharge[hit_sca[ihit]]->Fill(nchannelswithEsmaller05_sca[hit_sca[ihit]],hit_energy[ihit]);

      }

      for (Long64_t jentry2=jentry-20; jentry2<jentry+20;jentry2++) {
	Long64_t ientry2 = LoadTree(jentry2);
	if (ientry2 < 0) break;
	Long64_t nbytes2 = 0, nb2 = 0;
	nb2 = fChain->GetEntry(jentry2);   nbytes2 += nb2;

	//basic selection
	if(bcid<1250)  continue;
	if(bcid>bcid_max) continue;
      
	if( hit_sca[ihit]> 0) {
	  if(IsPedestal(ihit)==true && nchannelswithhits_sca[hit_sca[ihit]-1]>0 && hit_chip[ihit]==12 ) cout<<hit_z[ihit]<<endl;
	  if(IsPedestal(ihit)==true && nchannelswithhits_sca[hit_sca[ihit]-1]>0 && hit_chip[ihit]==12 && hit_z[ihit]==9) {
	    cout<<"--- "<<endl;
	    totalcharge_vs_pedcharge_prevsca[hit_sca[ihit]-1]->Fill(totalcharge_sca[hit_sca[ihit]-1],hit_energy[ihit]);
	    averagecharge_vs_pedcharge_prevsca[hit_sca[ihit]-1]->Fill(totalcharge_sca[hit_sca[ihit]-1]/nchannelswithhits_sca[hit_sca[ihit]-1],hit_energy[ihit]);
	    nhits_vs_pedcharge_prevsca[hit_sca[ihit]-1]->Fill(nchannelswithhits_sca[hit_sca[ihit]-1],hit_energy[ihit]);

	    nchannelsElarger5_vs_pedcharge_prevsca[hit_sca[ihit]-1]->Fill(nchannelswithElarger5_sca[hit_sca[ihit]-1],hit_energy[ihit]);
	    nchannelsElarger10_vs_pedcharge_prevsca[hit_sca[ihit]-1]->Fill(nchannelswithElarger10_sca[hit_sca[ihit]-1],hit_energy[ihit]);
	    nchannelsElarger20_vs_pedcharge_prevsca[hit_sca[ihit]-1]->Fill(nchannelswithElarger20_sca[hit_sca[ihit]-1],hit_energy[ihit]);
	    nchannelsEsmaller05_vs_pedcharge_prevsca[hit_sca[ihit]-1]->Fill(nchannelswithEsmaller05_sca[hit_sca[ihit]-1],hit_energy[ihit]);

	  }
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
    xbarycenter->Fit(fit_xbarycenter,"MQE");
  
    double xm_max=fit_xbarycenter->GetParameter(1)+2.*fit_xbarycenter->GetParameter(2);
    double xm_min=fit_xbarycenter->GetParameter(1)-2.*fit_xbarycenter->GetParameter(2);

    TF1 *fit_ybarycenter = new TF1("fit_ybarycenter","gaus");
    ybarycenter->Fit(fit_ybarycenter,"MQE");
  
    double ym_max=fit_ybarycenter->GetParameter(1)+2.*fit_ybarycenter->GetParameter(2);
    double ym_min=fit_ybarycenter->GetParameter(1)-2.*fit_ybarycenter->GetParameter(2);

    TF1 *fit_zbarycenter = new TF1("fit_zbarycenter","gaus");
    zbarycenter->Fit(fit_zbarycenter,"MEQ");
  
    double zm_max_05=fit_zbarycenter->GetParameter(1)+0.5*fit_zbarycenter->GetParameter(2);
    double zm_min_05=fit_zbarycenter->GetParameter(1)-0.5*fit_zbarycenter->GetParameter(2);

    double zm_max_1=fit_zbarycenter->GetParameter(1)+1.*fit_zbarycenter->GetParameter(2);
    double zm_min_1=fit_zbarycenter->GetParameter(1)-1.*fit_zbarycenter->GetParameter(2);

    double zm_max_2=fit_zbarycenter->GetParameter(1)+2.*fit_zbarycenter->GetParameter(2);
    double zm_min_2=fit_zbarycenter->GetParameter(1)-2.*fit_zbarycenter->GetParameter(2);

    // -----------------------------------------------------------------------------------------------------   
  
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
      double xm=0, ym=0., zm=0.;

      for(int ihit=0; ihit< nhit_chan; ihit ++) {
	int z=hit_z[ihit];
	if(z==9) z=6;
	if((hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]])>mipcut && IsHit(ihit)==true ) {
	  energy_sum_tmp += (hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]);
	  energy_X0_sum_tmp += w[int(hit_z[ihit])] * (hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]);
	  weight_X0_sum_tmp += w[int(hit_z[ihit])];

	  xm +=  - hit_x[ihit] * (hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]) * w[int(hit_z[ihit])];
	  ym +=  - hit_y[ihit] * (hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]) * w[int(hit_z[ihit])];
	  zm +=  hit_z[ihit] * (hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]) * w[int(hit_z[ihit])];
	}
      }


      // ---------------------------------------------
      //  recalculate the energy and other variables, with a cut in the distance (avoid isolated cells)
      energy_sum_tmp=0;
      weight_X0_sum_tmp=0;
      energy_X0_sum_tmp=0;
      xm=0, ym=0., zm=0.;

      int n_hitstotal=0;

      for(int ihit=0; ihit< nhit_chan; ihit ++) {
	int z=hit_z[ihit];
	if(z==9) z=6;
      
	if((hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]])>mipcut && IsHit(ihit)==true) {
	  charge_histo[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]->Fill(hit_energy[ihit]);
	}

	double mindist=MinDistCalc(ihit,mipcut);
    
	double mindist_th=10000.;
	if(RMIsolatedHits==true) mindist_th=7.0;

	if(mindist< mindist_th && IsHit(ihit)==true && (hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]])>mipcut ) {
	  n_hitstotal++;
	  energy_sum_tmp += (hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]);
	  energy_X0_sum_tmp += w[int(hit_z[ihit])] * (hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]);
	  weight_X0_sum_tmp += w[int(hit_z[ihit])];

	  xm +=  - hit_x[ihit] * (hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]) * w[int(hit_z[ihit])];
	  ym +=  - hit_y[ihit] * (hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]) * w[int(hit_z[ihit])];
	  zm +=  hit_z[ihit] * (hit_energy[ihit]-pedestal[z][hit_chip[ihit]][hit_chan[ihit]][hit_sca[ihit]]) * w[int(hit_z[ihit])];
	}
	
      }

      if(energy_sum_tmp<0.5) continue;

      xm = xm / energy_X0_sum_tmp;
      ym = ym / energy_X0_sum_tmp;
      zm = zm / energy_X0_sum_tmp;

      if( zm > zm_min_05 && zm < zm_max_05 ) {
	energy_center_05->Fill(energy_X0_sum_tmp);
	nhitstotal_center_05->Fill(n_hitstotal);
      }

      if( zm > zm_min_1 && zm < zm_max_1 ) {
	energy_center_1->Fill(energy_X0_sum_tmp);
	nhitstotal_center_1->Fill(n_hitstotal);
      }

      if( zm > zm_min_2 && zm < zm_max_2 ) {
	energy_center_2->Fill(energy_X0_sum_tmp);
	nhitstotal_center_2->Fill(n_hitstotal);
      }
    


    }

    // Signal analysis
    TFile *file = new TFile(TString::Format("%s%s_%s_%s_mipcut%0.1f_showers.root",folder.Data(),configuration.Data(),gridpoint.Data(),energy_string.Data(),mipcut) , "RECREATE");
    file->cd();

    for(int isca=0; isca<15; isca++) {
      pull_beamspot[isca]->SetName(TString::Format("bs_pedestal_pull_distribution_isca%i",isca));
      pull_beamspot[isca]->Write();
      stdev_beamspot[isca]->SetName(TString::Format("bs_pedestal_stdev_distribution_isca%i",isca));
      stdev_beamspot[isca]->Write();
      mean_beamspot[isca]->SetName(TString::Format("bs_pedestal_mean_distribution_isca%i",isca));
      mean_beamspot[isca]->Write();

      pull_outbeamspot[isca]->SetName(TString::Format("pedestal_pull_distribution_isca%i",isca));
      pull_outbeamspot[isca]->Write();
      stdev_outbeamspot[isca]->SetName(TString::Format("pedestal_stdev_distribution_isca%i",isca));
      stdev_outbeamspot[isca]->Write();
      mean_outbeamspot[isca]->SetName(TString::Format("pedestal_mean_distribution_isca%i",isca));
      mean_outbeamspot[isca]->Write();
    
    }
  
    mean_vs_entries->Write();
    rms_vs_entries->Write();
    mean_vs_rms->Write();


    pull_map_xy->Write();
    pull_map_xz->Write();

    pull_map_bad_xy->Write();
    pull_map_bad_xz->Write();
  
    for(int ichip=0; ichip<16; ichip++) {
      not_recalculated_sca[ichip]->Write();
      goodsca[ichip]->Write();
      wrongsca[ichip]->Write();
    }

    energy->Write();
  
    energy_center_05->Write();
    energy_center_1->Write();
    energy_center_2->Write();
  
    for(int isca=0; isca<15; isca++) {
      totalcharge_vs_pedcharge[isca]->Write();
      averagecharge_vs_pedcharge[isca]->Write();
      nhits_vs_pedcharge[isca]->Write();
    
      nchannelsElarger5_vs_pedcharge[isca]->Write();
      nchannelsElarger10_vs_pedcharge[isca]->Write();
      nchannelsElarger20_vs_pedcharge[isca]->Write();
      nchannelsEsmaller05_vs_pedcharge[isca]->Write();

      if(isca>0){
	totalcharge_vs_pedcharge_prevsca[isca-1]->Write();
	averagecharge_vs_pedcharge_prevsca[isca-1]->Write();
	nhits_vs_pedcharge_prevsca[isca-1]->Write();

	nchannelsElarger5_vs_pedcharge_prevsca[isca-1]->Write();
	nchannelsElarger10_vs_pedcharge_prevsca[isca-1]->Write();
	nchannelsElarger20_vs_pedcharge_prevsca[isca-1]->Write();
	nchannelsEsmaller05_vs_pedcharge_prevsca[isca-1]->Write();

      }
    }


  

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
    nhitstotal_center_05->Write();
    nhitstotal_center_1->Write();
    nhitstotal_center_2->Write();
  
    nhits_vs_energy->Write();
    min_dist_vs_energy->Write();

  
    TDirectory *cdhisto = file->mkdir("histograms");
    cdhisto->cd();
 
    for(int iz=0; iz<7; iz++) {
      for(int ichip=0; ichip<16; ichip++) {
	for(int ichn=0; ichn<64; ichn++) {
	  for(int isca=0; isca<15; isca++){
	    if( pedestal_histo[iz][ichip][ichn][isca]->GetEntries()>50) pedestal_histo[iz][ichip][ichn][isca]->Write();
	    charge_histo[iz][ichip][ichn][isca]->Write();
	  }
	}
      }
    }

    file->Close();
  
  
  }

}

void protoAnalysis::ShowerNoiseDistributions(TString folder="", TString configuration="conf1", double mipcut=0.5, int bcid_max=2800, int nslabs_selection=6, bool RMIsolatedHits=true)
{


  ReadCalibrated("../mip_calib/MIP_");

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

  TH1F* bcid_vs_next_1slab= new TH1F("bcid_vs_next_1slab","bcid_vs_next_1slab",300,-497.5,1002.5);
  TH2F* bcid_vs_next_map_1slab= new TH2F("bcid_vs_next_map_1slab","bcid_vs_next_map_1slab",700,2.5,3502.5,300,-497.5,1002.5);
  TH1F* bcid_vs_next_1slab_different= new TH1F("bcid_vs_next_1slab_different","bcid_vs_next_1slab_different",300,-497.5,1002.5);
  TH2F* bcid_vs_next_map_1slab_different= new TH2F("bcid_vs_next_map_1slab_different","bcid_vs_next_map_1slab_different",700,2.5,3502.5,300,-497.5,1002.5);

  TH2F* chip_correl_1slab= new TH2F("chip_correl_1slab","chip_correl_1slab",16,-0.5,15.5,16,-0.5,15.5);
  TH2F* chip_slab_1slab= new TH2F("chip_slab_1slab","chip_slab_1slab",16,-0.5,15.5,7,-0.5,6.5);

    
  TH1F* bcid_vs_next_2slab= new TH1F("bcid_vs_next_2slab","bcid_vs_next_2slab",300,-497.5,1002.5);
  TH2F* bcid_vs_next_map_2slab= new TH2F("bcid_vs_next_map_2slab","bcid_vs_next_map_2slab",700,2.5,3502.5,300,-497.5,1002.5);

  TH1F* bcid_vs_next_3slab= new TH1F("bcid_vs_next_3slab","bcid_vs_next_3slab",300,-497.5,1002.5);
  TH2F* bcid_vs_next_map_3slab= new TH2F("bcid_vs_next_map_3slab","bcid_vs_next_map_3slab",700,2.5,3502.5,300,-497.5,1002.5);

  TH1F* bcid_vs_next_4slab= new TH1F("bcid_vs_next_4slab","bcid_vs_next_4slab",300,-497.5,1002.5);
  TH2F* bcid_vs_next_map_4slab= new TH2F("bcid_vs_next_map_4slab","bcid_vs_next_map_4slab",700,2.5,3502.5,300,-497.5,1002.5);

  TH1F* bcid_vs_next_5slab= new TH1F("bcid_vs_next_5slab","bcid_vs_next_5slab",300,-497.5,1002.5);
  TH2F* bcid_vs_next_map_5slab= new TH2F("bcid_vs_next_map_5slab","bcid_vs_next_map_5slab",700,2.5,3502.5,300,-497.5,1002.5);

  TH2F* nhits_vs_energy= new TH2F("nhits_vs_energy","nhits_vs_energy",51,-2.5,252.5,101,-5,1005);
  TH2F* min_dist_vs_energy= new TH2F("min_dist_vs_energy","min_dist_vs_energy",250,0,250,100,-5,995);

  int number_hits=0;
  // -----------------------------------------------------------------------------------------------------
  Long64_t nbytes = 0, nb = 0;

 
  // Signal readout
  nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries-1;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(ShowerAntiSelection(mipcut,bcid_max,nslabs_selection)==false) continue;

    int islab_withhit=-1;
    int chip_1[1000];
    for(int ichip=0; ichip<1000; ichip++) chip_1[ichip]=-1;
    // ---------------------------------------------
    //   the energy and other variables,
    double energy_sum_tmp=0;
    double weight_X0_sum_tmp=0, energy_X0_sum_tmp_dist=0;
    double energy_X0_sum_tmp=0;
    double xm=0, ym=0., zm=0.;

    int ichip=0;
    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      if(hit_energy[ihit]>mipcut && IsHit(ihit)==true) {
	energy_X0_sum_tmp_dist += w[int(hit_z[ihit])] * (hit_energy[ihit]);
	if(nslabs_selection==1) {
	  islab_withhit=hit_z[ihit];
	  if(islab_withhit==9) islab_withhit=6;
	  ichip++;
	  chip_1[ichip]=hit_chip[ihit];
	}
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
	  
	  xm +=  - hit_x[ihit] * (hit_energy[ihit]) * w[int(hit_z[ihit])];
	  ym +=  - hit_y[ihit] * (hit_energy[ihit]) * w[int(hit_z[ihit])];
	  zm +=  hit_z[ihit] * (hit_energy[ihit]) * w[int(hit_z[ihit])];
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

    for(int ihit=0; ihit< nhit_chan; ihit ++) {

      if(hit_energy[ihit]>mipcut && IsHit(ihit)==true ) {
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

    /// correlation of good hits with noise hits
    int spill_1=spill;
    int spill_2=0;
    int bcid_1=bcid;
    cout<< jentry<< " ------ Spill= "<<spill<<" bcid = "<<bcid<< " ";
    bool event_selected[5];
    for(int i=0; i<5; i++) event_selected[i]=false;
    
    for (Long64_t jentry2=jentry+1; jentry2<nentries;jentry2++) {
      Long64_t ientry2 = LoadTree(jentry2);
      if (ientry2 < 0) break;
      Long64_t nbytes2 = 0, nb2 = 0;
      nb2 = fChain->GetEntry(jentry2);   nbytes2 += nb2;

      spill_2=spill;
      if(spill_2==spill_1) { 


	
	if(ShowerAntiSelection(mipcut,bcid_max,5)==true && event_selected[4]==false) {
	  bcid_vs_next_5slab->Fill(bcid-bcid_1);
	  bcid_vs_next_map_5slab->Fill(bcid,bcid-bcid_1);
	  event_selected[4]=true;
	}
	if(ShowerAntiSelection(mipcut,bcid_max,4)==true && event_selected[3]==false) {
	  bcid_vs_next_4slab->Fill(bcid-bcid_1);
	  bcid_vs_next_map_4slab->Fill(bcid,bcid-bcid_1);
	  event_selected[3]=true;
	}
	if(ShowerAntiSelection(mipcut,bcid_max,3)==true && event_selected[2]==false) {
	  bcid_vs_next_3slab->Fill(bcid-bcid_1);
	  bcid_vs_next_map_3slab->Fill(bcid,bcid-bcid_1);
	  event_selected[2]=true;
	}
	if(ShowerAntiSelection(mipcut,bcid_max,2)==true && event_selected[1]==false) {
	  bcid_vs_next_2slab->Fill(bcid-bcid_1);
	  bcid_vs_next_map_2slab->Fill(bcid,bcid-bcid_1);
	  event_selected[1]=true;
	}
	if(ShowerAntiSelection(mipcut,bcid_max,1)==true && event_selected[0]==false) {
	  cout<<bcid<<" ";
	  event_selected[0]=true;

	  if(nslabs_selection==1) {
	    int z=-1;
	    int chip_2=-1;
	    for(int ihit=0; ihit< nhit_chan; ihit ++) {
	      if( hit_energy[ihit]>mipcut && IsHit(ihit)==true ) {
		z=hit_z[ihit];
		if(z==9) z=6;
		chip_2=hit_chip[ihit];
	      }
	    }
	    if(z==islab_withhit) {
	      cout<<" z, islab ="<<z<<" "<<islab_withhit<<endl;
	      bcid_vs_next_1slab->Fill(bcid-bcid_1);
	      bcid_vs_next_map_1slab->Fill(bcid,bcid-bcid_1);
	      for(int ichip=0; ichip<1000; ichip++) {
		if(chip_1[ichip]>-0.5)  chip_correl_1slab->Fill(chip_1[ichip],chip_2);
		chip_slab_1slab->Fill(chip_2,z);
	      }
	    } else {
	      bcid_vs_next_1slab_different->Fill(bcid-bcid_1);
	      bcid_vs_next_map_1slab_different->Fill(bcid,bcid-bcid_1);
	    }
	  } else {
	    bcid_vs_next_1slab->Fill(bcid-bcid_1);
	    bcid_vs_next_map_1slab->Fill(bcid,bcid-bcid_1);
	  }
	    
	}
	  
      } else break;

    }
	
    cout<<endl;


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



  
  TFile *file = new TFile(TString::Format("%s%s_mipcut%0.1f_showers.root",folder.Data(),configuration.Data(),mipcut) , "RECREATE");
  file->cd();

  bcid_vs_next_map_1slab->Write();
  bcid_vs_next_1slab->Write();
  bcid_vs_next_map_1slab_different->Write();
  bcid_vs_next_1slab_different->Write();
  bcid_vs_next_map_2slab->Write();
  chip_correl_1slab->Write();
  chip_slab_1slab->Write();
  bcid_vs_next_2slab->Write();
  bcid_vs_next_map_3slab->Write();
  bcid_vs_next_3slab->Write();
  bcid_vs_next_map_4slab->Write();
  bcid_vs_next_4slab->Write();
  bcid_vs_next_map_5slab->Write();
  bcid_vs_next_5slab->Write();
  
  energy->Write();
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
	if( (fabs(hit_x[ihit]-hit_x[0]) > 0 || fabs(hit_y[ihit]-hit_y[0]) > 0) && hit_isHit[ihit]==1 ) nout++;
	if( hit_energy[ihit]> 0.5 && hit_isHit[ihit]==1 ) ngoodhit++;

      }

      if(nout > 0 && ngoodhit!=7 ) track=false;
	
      if(track==true && bcid>1250 && bcid<2900) {
	BCID_PREV_track->Fill(bcid,bcid-prev_bcid);
	BCID_track->Fill(bcid);
	for(int ihit=0; ihit< nhit_chan; ihit ++) {
	  if( hit_isHit[ihit]==1 ) {
	    SCA_track->Fill(hit_sca[ihit]);
	    SCA_BCID_track->Fill(hit_sca[ihit],bcid);
	  }
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





void protoAnalysis::ShowerDistributions(TString folder="", TString configuration="conf1", TString energy_string="3GeV", TString gridpoint ="grid20", double mipcut=0.5, int bcid_max=2900, int nslabs_selection=6, bool RMIsolatedHits=true)
{


  ReadCalibrated("../mip_calib/MIP_");

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
  TH1F* energy_center_05 = new TH1F("energy_center_zm05","energy_center_zm05",500,1,1001);
  TH1F* energy_center_1 = new TH1F("energy_center_zm1","energy_center_zm1",500,1,1001);
  TH1F* energy_center_2 = new TH1F("energy_center_zm2","energy_center_zm2",500,1,1001);
  
  TH1F* energy_profile_z = new TH1F("energy_profile_z","energy_profile_z",7,-0.5,6.5);//binnumz,binsz);
  TH1F* energy_profile_z_center_05 = new TH1F("energy_profile_z_center_05","energy_profile_z_center_05",7,-0.5,6.5);//binnumz,binsz);
  TH1F* energy_profile_z_center_1 = new TH1F("energy_profile_z_center_1","energy_profile_z_center_1",7,-0.5,6.5);//binnumz,binsz);
  TH1F* energy_profile_z_center_2 = new TH1F("energy_profile_z_center_2","energy_profile_z_center_2",7,-0.5,6.5);//binnumz,binsz);

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
  TH1F* nhitstotal_center_05 = new TH1F("nhitstotal_center_zm05","nhitstotal_center_zm05",500,1,1001);
  TH1F* nhitstotal_center_1 = new TH1F("nhitstotal_center_zm1","nhitstotal_center_zm1",500,1,1001);
  TH1F* nhitstotal_center_2 = new TH1F("nhitstotal_center_zm2","nhitstotal_center_zm2",500,1,1001);

  

  TH2F* nhits_vs_energy= new TH2F("nhits_vs_energy","nhits_vs_energy",51,-2.5,252.5,101,-5,1005);
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
	  
	  xm +=  - hit_x[ihit] * (hit_energy[ihit]) * w[int(hit_z[ihit])];
	  ym +=  - hit_y[ihit] * (hit_energy[ihit]) * w[int(hit_z[ihit])];
	  zm +=  hit_z[ihit] * (hit_energy[ihit]) * w[int(hit_z[ihit])];
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

    for(int ihit=0; ihit< nhit_chan; ihit ++) {

      if(hit_energy[ihit]>mipcut && IsHit(ihit)==true ) {
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

   
    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      if(hit_energy[ihit]>mipcut && IsHit(ihit)==true ) {
	double zlayer=hit_z[ihit];
	if(hit_z[ihit]>8) zlayer=6;
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
  TFile *file = new TFile(TString::Format("%s%s_%s_%s_mipcut%0.1f_showers.root",folder.Data(),configuration.Data(),gridpoint.Data(),energy_string.Data(),mipcut) , "RECREATE");
  file->cd();

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
  min_dist_vs_energy->Write();
  file->Close();
  
}

void protoAnalysis::SimpleMIPAnalysis(TString outputname="", int nin=5, int nout=10, int bcid_max=10000)
{


  ReadCalibrated("../mip_calib/MIP_");

  ofstream fout_mip("results_mipcalibration/eff_files/MIP"+outputname+".txt",ios::out);
  fout_mip<<"#mip results, all channels together SimpleMIPAnalysis using tracks"<<endl;
  fout_mip<<"#mpv empv widthmpv widthgauss chi2ndf nentries"<<endl;

  // --------------
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //-----------------

  //Objects to store the chip/channel results:
  // numbers 
  std::vector<std::vector<std::vector<TH1F*> > > mip_histo;
  TH1F* mip_histo_all = new TH1F("mip_histo_all","mip_histo_all",250,0.1,10.1);
  TH1F* ped_histo_all = new TH1F("ped_histo_all","ped_histo_all",500,-10,10);

  //inefficiencies
  std::vector<TH2F*> mip_ineff_chip_chn;
  std::vector<TH2F*> mip_ineff;
  std::vector<TH2F*> mip_tracks_chip_chn;
  std::vector<TH2F*> mip_tracks;
  std::vector<TH2F*> mip_ineff_full_chip_chn;
  std::vector<TH2F*> mip_ineff_full;
  
  // sca, bcid, bcid-bcid_prev plots
  TH1F* SCA_track = new TH1F("SCA_track","SCA_track",15,-0.5,14.5);
  TH1F* BCID_track = new TH1F("BCID_track","BCID_track",60,25,3025);
  TH2F* SCA_BCID_track = new TH2F("SCA_BCID_track","SCA_BCID_track",15,-0.5,14.5,60,25,3025);
  TH2F* BCID_PREV_track = new TH2F("BCID_PREV_track","BCID_PREV_track",1200,2.5,6002.5,200,2.5,1002.5);
  TH2F* BCID_NEXT_track = new TH2F("BCID_NEXT_track","BCID_NEXT_track",1200,2.5,6002.5,200,2.5,1002.5);


  Float_t binsx[33];
  binsx[0]=-86.;
  for(int ibins=1;ibins<33;ibins++) binsx[ibins]=binsx[ibins-1]+5.5;
  Int_t  binnumx = 32; 

  for(int islab=0; islab<7; islab++) {

    TString histo_title_2=TString::Format("mip_tracks_slab%i_chipchn",islab);
    TH2F *mip_tracks_chip_chn_tmp = new TH2F(histo_title_2,histo_title_2,16,-0.5,15.5,64,-0.5,63.5);
    histo_title_2=TString::Format("mip_tracks_slab%i_xy",islab);
    TH2F *mip_tracks_tmp = new TH2F(histo_title_2,histo_title_2,binnumx,binsx,binnumx,binsx);
    mip_tracks.push_back(mip_tracks_tmp);
    mip_tracks_chip_chn.push_back(mip_tracks_chip_chn_tmp);

    histo_title_2=TString::Format("mip_ineff_slab%i_chipchn",islab);
    TH2F *mip_ineff_chip_chn_tmp = new TH2F(histo_title_2,histo_title_2,16,-0.5,15.5,64,-0.5,63.5);
    histo_title_2=TString::Format("mip_ineff_slab%i_xy",islab);
    TH2F *mip_ineff_tmp = new TH2F(histo_title_2,histo_title_2,binnumx,binsx,binnumx,binsx);
    mip_ineff.push_back(mip_ineff_tmp);
    mip_ineff_chip_chn.push_back(mip_ineff_chip_chn_tmp);

    histo_title_2=TString::Format("mip_ineff_full_slab%i_chipchn",islab);
    TH2F *mip_ineff_full_chip_chn_tmp = new TH2F(histo_title_2,histo_title_2,16,-0.5,15.5,64,-0.5,63.5);
    histo_title_2=TString::Format("mip_ineff_full_slab%i_xy",islab);
    TH2F *mip_ineff_full_tmp = new TH2F(histo_title_2,histo_title_2,binnumx,binsx,binnumx,binsx);
    mip_ineff_full.push_back(mip_ineff_full_tmp);
    mip_ineff_full_chip_chn.push_back(mip_ineff_full_chip_chn_tmp);

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


  int ntracks=0;
 
  // -----------------------------------------------------------------------------------------------------   
  // Signal readout
  Long64_t nbytes = 0, nb = 0;
  //  nentries=500000;
  for (Long64_t jentry=100; jentry<nentries-100;jentry++){//nentries-100;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    int fivepercent=nentries/20;
    if( (jentry % fivepercent) == 0 ) cout<< "Progress: "<<100.*jentry/nentries<<" %" <<endl;

    if( bcid<1290 || bcid>bcid_max) continue;


    bool track=true; //track candidate
    
    //get the hit with largest energy as the reference for track identification
    int imaxenergy=0;
    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      int z=hit_z[ihit];
      if(hit_z[ihit]==9) z=6;
      if( ( hit_isMasked[ihit] == 0 && calibrated[z][hit_chip[ihit]][hit_chan[ihit]]==0 ) && hit_isHit[ihit]==1 && hit_energy[ihit]>0.5&& hit_energy[ihit]>hit_energy[imaxenergy]) {
	imaxenergy=ihit;
      }
    }
    // count hits in the perpendicular track
    int nintrack=0;
    int nouttrack=0;
    
    for(int ihit=0; ihit< nhit_chan; ihit ++) {
      int z=hit_z[ihit];
      if(hit_z[ihit]==9) z=6;

      if( fabs(hit_x[ihit]-hit_x[imaxenergy])< 8.6 && fabs(hit_y[ihit]-hit_y[imaxenergy])< 8.6  && ( hit_isMasked[ihit] == 0 && calibrated[z][hit_chip[ihit]][hit_chan[ihit]]==0 ) && hit_isHit[ihit]==1 && hit_energy[ihit]>0.15) nintrack++;
    
      if( (fabs(hit_x[ihit]-hit_x[imaxenergy]) > 8.6 || fabs(hit_y[ihit]-hit_y[imaxenergy]) > 8.6 )  && ( hit_isMasked[ihit] == 0 && calibrated[z][hit_chip[ihit]][hit_chan[ihit]]==0 ) && hit_isHit[ihit]==1 && hit_energy[ihit]>0.15) nouttrack++;
    }
    
    if(nintrack < nin ) track=false;
    if(nouttrack > nout ) track=false;
    
    if(track==true) {

      BCID_PREV_track->Fill(bcid,bcid-prev_bcid);
      BCID_NEXT_track->Fill(bcid,bcid-next_bcid);
      BCID_track->Fill(bcid);
      
      ntracks++;

      bool workinglayer[7];
      bool infolayer[7];
      bool layerwithhitintrack[7];
      
      
      for(int iz=0; iz<7;iz++) {
	workinglayer[iz]=false;
	infolayer[iz]=false;
	layerwithhitintrack[iz]=false;
	if(  calibrated[iz][hit_chip[imaxenergy]][hit_chan[imaxenergy]]==0  )  workinglayer[iz]=true;
      }


      for(int ihit=0; ihit< nhit_chan; ihit ++) {
	int z=hit_z[ihit];
	if(hit_z[ihit]==9) z=6;

	if( (fabs(hit_x[ihit]-hit_x[imaxenergy])< 8.6 && fabs(hit_y[ihit]-hit_y[imaxenergy])< 8.6 )  ) {

	  infolayer[z]=true;

	  if( ( hit_isMasked[ihit] == 0 && calibrated[z][hit_chip[ihit]][hit_chan[ihit]]==0 ) && hit_isHit[ihit]==1 && hit_energy[ihit]>0.15) {
	    mip_histo_all->Fill(hit_energy[ihit]);
	    SCA_track->Fill(hit_sca[ihit]);
	    SCA_BCID_track->Fill(hit_sca[ihit],bcid);
	    mip_tracks.at(z)->Fill(hit_x[ihit],hit_y[ihit]);
	    mip_tracks_chip_chn.at(z)->Fill(hit_chip[ihit],hit_chan[ihit]);
	    layerwithhitintrack[z]=true;
	  }

	}
      }
      
      for(int ihit=0; ihit< nhit_chan; ihit ++) {
	int z=hit_z[ihit];
	if(hit_z[ihit]==9) z=6;
	
	if( (fabs(hit_x[ihit]-hit_x[imaxenergy])< 8.6 && fabs(hit_y[ihit]-hit_y[imaxenergy])< 8.6 )  ) {
	  
	  if( ( hit_isMasked[ihit] == 0 && calibrated[z][hit_chip[ihit]][hit_chan[ihit]]==0 ) && hit_isHit[ihit]==0 && layerwithhitintrack[z]==false) {
	    int z=hit_z[ihit];
	    if(hit_z[ihit]==9) z=6;
	    
	    if( hit_energy[ihit]>0.15) {
	      mip_ineff.at(z)->Fill(hit_x[ihit],hit_y[ihit]);
	      mip_ineff_chip_chn.at(z)->Fill(hit_chip[ihit],hit_chan[ihit]);
	    }

	    ped_histo_all->Fill(hit_energy[ihit]);    
	  }
	  
	}// x=x, y=y
	
      }//for nhits

      for(int z=0; z<7;z++) {
	if(workinglayer[z]==true && infolayer[z]==false && layerwithhitintrack[z]==false) {
	  double ix=hit_x[imaxenergy];
	  double iy=hit_y[imaxenergy];
	  double iz=hit_z[imaxenergy];
	  int spill_ref=spill;
	  int spill_new=-1;
	  int last_sca=-1;
	  for (Long64_t jentry2=jentry-200; jentry2<jentry+200;jentry2++) {
	    Long64_t ientry2 = LoadTree(jentry2);
	    if (ientry2 < 0) break;
	    Long64_t nbytes2 = 0, nb2 = 0;
	    nb2 = fChain->GetEntry(jentry2);   nbytes2 += nb2;
	    
	    spill_new=spill;
	    if(spill_new<spill_ref) continue;
	    if(spill_new>spill_ref) break;

	    for(int ihit=0; ihit< nhit_chan; ihit ++) {
	      if( spill_new==spill_ref && fabs( ix - hit_x[ihit]) < 8.6 && fabs(iy - hit_y[ihit])< 8.6 && iz==hit_z[ihit]) last_sca=hit_sca[ihit];
	    }
	  }

	  if(last_sca<13) {
	    mip_ineff.at(z)->Fill(hit_x[imaxenergy],hit_y[imaxenergy]);
	    mip_ineff_chip_chn.at(z)->Fill(hit_chip[imaxenergy],hit_chan[imaxenergy]);
	  } else {
	    mip_ineff_full.at(z)->Fill(hit_x[imaxenergy],hit_y[imaxenergy]);
	    mip_ineff_full_chip_chn.at(z)->Fill(hit_chip[imaxenergy],hit_chan[imaxenergy]);
	  }
	  
	  
	}
      }

    }//if track
    
  }

  // -----------------------------------------------------------------------------------------------------   
  // Signal analysis
  TFile *signalfile_summary = new TFile("results_mipcalibration/eff_files/Signal_summary"+outputname+".root" , "RECREATE");
  signalfile_summary->cd();

  ped_histo_all->GetXaxis()->SetTitle("Energy [MIP]");
  ped_histo_all->GetYaxis()->SetTitle("# entries");
  ped_histo_all->SetTitle("Single cell pedestal value");
  ped_histo_all->SetName("pedestal_distribution");
  ped_histo_all->Write();
  
  mip_histo_all->GetXaxis()->SetTitle("Energy [MIP]");
  mip_histo_all->GetYaxis()->SetTitle("# entries");
  mip_histo_all->SetTitle("Single cell hit energy in 3GeV e^{+} beam");
  mip_histo_all->SetName("energy_distribution");
  mip_histo_all->Write();

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

  BCID_NEXT_track->GetXaxis()->SetTitle("bcid");
  BCID_NEXT_track->GetYaxis()->SetTitle("bcid-prev_bcid");
  BCID_NEXT_track->SetTitle("bcid vs bcid-next_bcid distribution for track-like events");
  BCID_NEXT_track->SetName("BCID_NEXT_track");
  BCID_NEXT_track->Write();
  
  for(int i=0; i<7; i++) {
    mip_tracks.at(i)->GetXaxis()->SetTitle("chip");
    mip_tracks.at(i)->GetYaxis()->SetTitle("channel");
    mip_tracks.at(i)->SetTitle(TString::Format("Total hits per channel in tracks, layer %i ",i));
    mip_tracks.at(i)->SetName(TString::Format("mip_tracks_layer%i",i));
    mip_tracks.at(i)->Write();
  }
  for(int i=0; i<7; i++) {
    mip_tracks_chip_chn.at(i)->GetXaxis()->SetTitle("chip");
    mip_tracks_chip_chn.at(i)->GetYaxis()->SetTitle("channel");
    mip_tracks_chip_chn.at(i)->SetTitle(TString::Format("Total hits per channel in tracks, layer %i ",i));
    mip_tracks_chip_chn.at(i)->SetName(TString::Format("mip_tracks_layer%i_chip_chn",i));
    mip_tracks_chip_chn.at(i)->Write();
  }
  
  for(int i=0; i<7; i++) {
    for(int ix=0; ix< mip_ineff.at(i)->GetNbinsX(); ix++) {
      for(int iy=0; iy< mip_ineff.at(i)->GetNbinsY(); iy++) {
	double norm=mip_tracks.at(i)->GetBinContent(ix+1,iy+1)/0.01;
	if(norm>0) mip_ineff.at(i)->SetBinContent(ix+1,iy+1,mip_ineff.at(i)->GetBinContent(ix+1,iy+1)/norm);
	else mip_ineff.at(i)->SetBinContent(ix+1,iy+1,0);
      }
    }
    mip_ineff.at(i)->GetXaxis()->SetTitle("x");
    mip_ineff.at(i)->GetYaxis()->SetTitle("y");
    mip_ineff.at(i)->SetTitle(TString::Format("MIP inefficiency layer %i ",i)+"[%]");
    mip_ineff.at(i)->SetName(TString::Format("mip_ineff_layer%i",i));
    mip_ineff.at(i)->Write();
  }
  for(int i=0; i<7; i++) {
    for(int ix=0; ix< mip_ineff_chip_chn.at(i)->GetNbinsX(); ix++) {
      for(int iy=0; iy< mip_ineff_chip_chn.at(i)->GetNbinsY(); iy++) {
	double norm=mip_tracks_chip_chn.at(i)->GetBinContent(ix+1,iy+1)/0.01;
	if(norm>0) mip_ineff_chip_chn.at(i)->SetBinContent(ix+1,iy+1,mip_ineff_chip_chn.at(i)->GetBinContent(ix+1,iy+1)/norm);
	else mip_ineff_chip_chn.at(i)->SetBinContent(ix+1,iy+1,0);
      }
    }
    mip_ineff_chip_chn.at(i)->GetXaxis()->SetTitle("chip");
    mip_ineff_chip_chn.at(i)->GetYaxis()->SetTitle("channel");
    mip_ineff_chip_chn.at(i)->SetTitle(TString::Format("MIP inefficiency layer %i ",i)+"[%]");
    mip_ineff_chip_chn.at(i)->SetName(TString::Format("mip_ineff_layer%i_chip_chn",i));
    mip_ineff_chip_chn.at(i)->Write();
  }

  for(int i=0; i<7; i++) {
    for(int ix=0; ix< mip_ineff_full.at(i)->GetNbinsX(); ix++) {
      for(int iy=0; iy< mip_ineff_full.at(i)->GetNbinsY(); iy++) {
	double norm=mip_tracks.at(i)->GetBinContent(ix+1,iy+1)/0.01;
	if(norm>0) mip_ineff_full.at(i)->SetBinContent(ix+1,iy+1,mip_ineff_full.at(i)->GetBinContent(ix+1,iy+1)/norm);
	else mip_ineff_full.at(i)->SetBinContent(ix+1,iy+1,0);
      }
    }
    mip_ineff_full.at(i)->GetXaxis()->SetTitle("x");
    mip_ineff_full.at(i)->GetYaxis()->SetTitle("y");
    mip_ineff_full.at(i)->SetTitle(TString::Format("MIP inefficiency for blindness layer %i ",i)+"[%]");
    mip_ineff_full.at(i)->SetName(TString::Format("mip_ineff_full_layer%i",i));
    mip_ineff_full.at(i)->Write();
  }
  for(int i=0; i<7; i++) {
    for(int ix=0; ix< mip_ineff_full_chip_chn.at(i)->GetNbinsX(); ix++) {
      for(int iy=0; iy< mip_ineff_full_chip_chn.at(i)->GetNbinsY(); iy++) {
	double norm=mip_tracks_chip_chn.at(i)->GetBinContent(ix+1,iy+1)/0.01;
	if(norm>0) mip_ineff_full_chip_chn.at(i)->SetBinContent(ix+1,iy+1,mip_ineff_full_chip_chn.at(i)->GetBinContent(ix+1,iy+1)/norm);
	else mip_ineff_full_chip_chn.at(i)->SetBinContent(ix+1,iy+1,0);
      }
    }
    mip_ineff_full_chip_chn.at(i)->GetXaxis()->SetTitle("chip");
    mip_ineff_full_chip_chn.at(i)->GetYaxis()->SetTitle("channel");
    mip_ineff_full_chip_chn.at(i)->SetTitle(TString::Format("MIP inefficiency for blindness layer %i ",i)+"[%]");
    mip_ineff_full_chip_chn.at(i)->SetName(TString::Format("mip_ineff_full_layer%i_chip_chn",i));
    mip_ineff_full_chip_chn.at(i)->Write();
  }

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

