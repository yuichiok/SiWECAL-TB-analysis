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



void protoAnalysis::SimpleMIPAnalysis(TString outputname="", TString map_filename="../fev10_chip_channel_x_y_mapping.txt")
{

  int maxnhit=5; // plane event threshold

  //Read the channel/chip -- x/y mapping
  ReadMap(map_filename);
  
  //Read Masked channels list
  ReadMasked();

  //Read the list of pedestals (this information contains, implicitily, the masked channels information )

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
	
	  if(track==true && bcid>1260 && bcid<2850) {
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

  /*

  // -----------------------------------------------------------------------------------
  // iteration 2
  // ***********************************************************************************

  //-----------------------------------
  // Peak 1
  //------------------------------------
  TH1F *hist_peak1_iter2 = (TH1F*)mip_histo_all->Clone();
  //hist_peak1_iter2->Add(fpeak2_iter1,-1);
  //hist_peak1_iter2->Add(fpeak3_iter1,-1);

  pllo[0]=fpeak1_iter1->GetParameter(0)-10*fpeak1_iter1->GetParameter(0); 
  pllo[1]=fpeak1_iter1->GetParameter(1)-10*fpeak1_iter1->GetParameter(1); 
  pllo[2]=fpeak1_iter1->GetParameter(2)-10*fpeak1_iter1->GetParameter(2); 
  pllo[3]=fpeak1_iter1->GetParameter(3)-10*fpeak1_iter1->GetParameter(3);

  plhi[0]=fpeak1_iter1->GetParameter(0)+10*fpeak1_iter1->GetParameter(0); 
  plhi[1]=fpeak1_iter1->GetParameter(1)+10*fpeak1_iter1->GetParameter(1); 
  plhi[2]=fpeak1_iter1->GetParameter(2)+10*fpeak1_iter1->GetParameter(2); 
  plhi[3]=fpeak1_iter1->GetParameter(3)+10*fpeak1_iter1->GetParameter(3);


  fr[0]=fpeak1_iter1->GetParameter(1)-1.5*fpeak1_iter1->GetParameter(0);
  fr[1]=fpeak1_iter1->GetParameter(1)+3*fpeak1_iter1->GetParameter(0);
  sv[0]=fpeak1_iter1->GetParameter(0);
  sv[1]=fpeak1_iter1->GetParameter(1);
  sv[2]=fpeak1_iter1->GetParameter(2);
  sv[3]=fpeak1_iter1->GetParameter(3);
  
  TF1 *fit_2nd_iter_temp=langaufit(hist_peak1_iter2,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  
  hist_peak1_iter2->GetXaxis()->SetTitle("Energy [MIP]");
  hist_peak1_iter2->GetYaxis()->SetTitle("# entries");
  hist_peak1_iter2->SetTitle("First MIP");
  hist_peak1_iter2->Write();

  mpv=fit_2nd_iter_temp->GetParameter(1);
  empv=fit_2nd_iter_temp->GetParError(1);
  wmpv=fit_2nd_iter_temp->GetParameter(0);
  wg=fit_2nd_iter_temp->GetParameter(3);
  chi2ndf=0;
  if(ndf>0) chi2ndf=chisqr/ndf;
  mipentries=mip_histo_all->GetEntries();

  TF1 *fpeak1_iter2 = new TF1("fpeak1_iter2",langaufun,0,5,4);
  fpeak1_iter2->SetParameter(0,fit_2nd_iter_temp->GetParameter(0));
  fpeak1_iter2->SetParameter(1,fit_2nd_iter_temp->GetParameter(1));
  fpeak1_iter2->SetParameter(2,fit_2nd_iter_temp->GetParameter(2));
  fpeak1_iter2->SetParameter(3,fit_2nd_iter_temp->GetParameter(3));
 
  fout_mip<<" "<<mpv<<" "<<empv<<" "<<wmpv<<" "<<wg<<" "<<chi2ndf<<" "<<mipentries<<endl;
  cout<<" "<<mpv<<" "<<empv<<" "<<wmpv<<" "<<wg<<" "<<chi2ndf<<" "<<mipentries<<endl;

  
  // ------------------------------------------------------
  // Peak 2
  // -----------------------------------------------------
  TH1F *hist_peak2_iter2 = (TH1F*)mip_histo_all->Clone();
  hist_peak2_iter2->Add(fpeak1_iter2,-1);

  TSpectrum *s1 = new TSpectrum(1);
  int npeaks_iter2_1 = s1->Search(hist_peak2_iter2,0.5,"",0.8); 
  Double_t *mean_peak_iter2_1=s1->GetPositionX();
  fr[0]=mean_peak_iter2_1[0]-0.2*mip_histo_all->GetRMS();//-0.65*mip_histo_all->GetRMS();
  fr[1]=mean_peak_iter2_1[0]+0.3*mip_histo_all->GetRMS();//mip_histo_all->GetMean();
  sv[0]=wmpv;
  sv[1]=2*mpv;
  sv[2]=hist_peak2_iter2->Integral("width");
  sv[3]=hist_peak2_iter2->GetRMS()/10.;

  pllo[1]=2*mpv - 15*wmpv;
  plhi[1]=2*mpv + 15*wmpv;

  TF1 *fit_2nd_iter_temp2=langaufit(hist_peak2_iter2,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  mpv=fit_2nd_iter_temp2->GetParameter(1);
  empv=fit_2nd_iter_temp2->GetParError(1);
  wmpv=fit_2nd_iter_temp2->GetParameter(0);
  wg=fit_2nd_iter_temp2->GetParameter(3);
  chi2ndf=0;
  if(ndf>0) chi2ndf=chisqr/ndf;

  fout_mip<<" "<<mpv<<" "<<empv<<" "<<wmpv<<" "<<wg<<" "<<chi2ndf<<" "<<mipentries<<endl;
  cout<<npeaks_iter2_1<<" "<<mean_peak_iter2_1[0]<<endl;
  cout<<" "<<mpv<<" "<<empv<<" "<<wmpv<<" "<<wg<<" "<<chi2ndf<<" "<<mipentries<<endl;

  hist_peak2_iter2->GetXaxis()->SetTitle("Energy [MIP]");
  hist_peak2_iter2->GetYaxis()->SetTitle("# entries");
  hist_peak2_iter2->SetTitle("Second MIP (first subtracted)");
  hist_peak2_iter2->Write();

  TF1 *fpeak2_iter2 = new TF1("fpeak2_iter2",langaufun,0,5,4);
  fpeak2_iter2->SetParameter(0,fit_2nd_iter_temp2->GetParameter(0));
  fpeak2_iter2->SetParameter(1,fit_2nd_iter_temp2->GetParameter(1));
  fpeak2_iter2->SetParameter(2,fit_2nd_iter_temp2->GetParameter(2));
  fpeak2_iter2->SetParameter(3,fit_2nd_iter_temp2->GetParameter(3));


  // ------------------------------------------------------
  // Peak 2
  // -----------------------------------------------------
  TH1F *hist_peak3_iter2 = (TH1F*)mip_histo_all->Clone();
  hist_peak3_iter2->Add(fpeak1_iter2,-1);
  hist_peak3_iter2->Add(fpeak2_iter2,-1);

  TSpectrum *s2 = new TSpectrum(1);
  int npeaks_iter2_2 = s2->Search(hist_peak3_iter2,0.5,"",0.8); 
  Double_t *mean_peak_iter2_2=s2->GetPositionX();
  fr[0]=mean_peak_iter2_2[0]-0.4*mip_histo_all->GetRMS();//-0.65*mip_histo_all->GetRMS();
  fr[1]=mean_peak_iter2_2[0]+0.4*mip_histo_all->GetRMS();//mip_histo_all->GetMean();
  sv[0]=wmpv;
  sv[1]=3/2*mpv;
  sv[2]=hist_peak3_iter2->Integral("width");
  sv[3]=hist_peak3_iter2->GetRMS()/10.;

  pllo[1]=3/2*mpv - 5*wmpv;
  plhi[1]=3/2*mpv + 5*wmpv;

  TF1 *fit_2nd_iter_temp3=langaufit(hist_peak3_iter2,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  mpv=fit_2nd_iter_temp3->GetParameter(1);
  empv=fit_2nd_iter_temp3->GetParError(1);
  wmpv=fit_2nd_iter_temp3->GetParameter(0);
  wg=fit_2nd_iter_temp3->GetParameter(3);
  chi2ndf=0;
  if(ndf>0) chi2ndf=chisqr/ndf;

  fout_mip<<" "<<mpv<<" "<<empv<<" "<<wmpv<<" "<<wg<<" "<<chi2ndf<<" "<<mipentries<<endl;
  cout<<npeaks_iter2_2<<" "<<mean_peak_iter2_2[0]<<endl;
  cout<<" "<<mpv<<" "<<empv<<" "<<wmpv<<" "<<wg<<" "<<chi2ndf<<" "<<mipentries<<endl;

  hist_peak3_iter2->GetXaxis()->SetTitle("Energy [MIP]");
  hist_peak3_iter2->GetYaxis()->SetTitle("# entries");
  hist_peak3_iter2->SetTitle("Third MIP (first+second subtracted)");
  hist_peak3_iter2->Write();

  TF1 *fpeak3_iter2 = new TF1("fpeak3_iter2",langaufun,0,5,4);
  fpeak3_iter2->SetParameter(0,fit_2nd_iter_temp3->GetParameter(0));
  fpeak3_iter2->SetParameter(1,fit_2nd_iter_temp3->GetParameter(1));
  fpeak3_iter2->SetParameter(2,fit_2nd_iter_temp3->GetParameter(2));
  fpeak3_iter2->SetParameter(3,fit_2nd_iter_temp3->GetParameter(3));

  // ------------------------------------------------------
  //write TF1s
  fpeak1_iter2->SetName("FirstPeak_iter2");
  fpeak1_iter2->Write();

  fpeak2_iter2->SetName("SecondPeak_iter2");
  fpeak2_iter2->Write();

  fpeak3_iter2->SetName("ThirdPeak_iter2");
  fpeak3_iter2->Write();


  TH1F* second_iteration_peak1 = new TH1F("second_iteration_peak1","second_iteration_peak1",925,0.005,1.855);
  TH1F* second_iteration_peak2 = new TH1F("second_iteration_peak2","second_iteration_peak2",500,1.855,2.855);
  TH1F* second_iteration_peak3 = new TH1F("second_iteration_peak3","second_iteration_peak3",500,2.855,3.855);
  second_iteration_peak1->Add(fpeak1_iter2);
  second_iteration_peak1->Add(fpeak2_iter2);
  second_iteration_peak1->Add(fpeak3_iter2);

  second_iteration_peak2->Add(fpeak1_iter2);
  second_iteration_peak2->Add(fpeak2_iter2);
  second_iteration_peak2->Add(fpeak3_iter2);

  second_iteration_peak3->Add(fpeak1_iter2);
  second_iteration_peak3->Add(fpeak2_iter2);
  second_iteration_peak3->Add(fpeak3_iter2);

  second_iteration_peak1->Write();
  second_iteration_peak2->Write();
  second_iteration_peak3->Write();
  */

  signalfile_summary->cd();

  signalfile_summary->Close();

}

void protoAnalysis::ReadMap(TString filename) 
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

void protoAnalysis::ReadMasked() 
{

  for(int islab=0; islab<7; islab++) {
    
    TString filename = "../masked/masked_dif_1_1_1.txt";
    if(islab==1) filename = "../masked/masked_dif_1_1_2.txt";
    if(islab==2) filename = "../masked/masked_dif_1_1_3.txt";
    if(islab==3) filename = "../masked/masked_dif_1_1_4.txt";
    if(islab==4) filename = "../masked/masked_dif_1_1_5.txt";
    if(islab==5) filename = "../masked/masked_dif_1_2_1.txt";
    if(islab==6) filename = "../masked/masked_dif_1_2_2.txt";

    std::ifstream reading_file(filename);
    if(!reading_file){
      cout<<" dameyo - damedame"<<endl;
    }
    
    for(int i=0; i<16; i++) {
      for(int j=0; j<64; j++) {
	masked[islab][i][j] = 0;
      }
    }
    
    Int_t tmp_chip = 0,tmp_channel = 0;
    Int_t tmp_masked = 0;
    TString tmpst;
    reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst ;
    reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst  ;
    
    cout<<"Read Masked: "<<filename<<endl;
    while(reading_file){
      reading_file >> tmp_chip >> tmp_channel >> tmp_masked ;
      masked[islab][tmp_chip][tmp_channel] = tmp_masked ;
    }
  }

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

