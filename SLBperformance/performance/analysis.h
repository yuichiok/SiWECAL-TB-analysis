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

void pedanalysis(TString run="Run_ILC_cosmic_test_11222019"){
  
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


  TFile *_file0 = TFile::Open(TString::Format("../results_calib/%s.root",run.Data()));

  TH2F* ped_all=new TH2F("ped_all","pedestal pos. ; Layer*20+Chip  ; SCA*100 +chn",300,-0.5,299.5,1500,-0.5,1499.5);
  TH2F* rms_all=new TH2F("rms_all","width of pedestal ; Layer*20+Chip  ; SCA*100 +chn",300,-0.5,299.5,1500,-0.5,1499.5);

  ofstream fout_ped(TString::Format("../results_calib/%s.txt",run.Data()).Data(),ios::out);
  fout_ped<<"#pedestal results (fit to a gaussian) remove channels/sca with two pedestals peaks from the analysis : PROTO15"<<endl;
  fout_ped<<"#layer chip channel ped0 eped0 widthped0 ped1 eped1 widthped1... ped14 eped14 widhtped14 (all SCA)"<<endl;

  for(int layer=0; layer<15; layer++) {

    
    // Comparing nbr entries in tag or not tag // GetWidth and Mean
    TString map="../../mapping/fev10_chip_channel_x_y_mapping.txt";
    //  if(layer==0 || layer==2)  map="../../mapping/fev11_cob_chip_channel_x_y_mapping.txt";
    ReadMap(map);

  
    TFile *_file0 = TFile::Open(TString::Format("../results_calib/%s.root",run.Data()));

 
    for(int i=0;i<nchips;i++){
    
      double nch=0;
      float avmean=0, avwidth=0;
    
      for(int j=0; j<64; j++) {
	fout_ped << layer<<" "<<i <<" " <<j<< " ";
      
	for(int isca=0; isca<15; isca++) {
	  TH1F* temp=(TH1F*)_file0->Get(TString::Format("layer_%i/ped_chip%i_chn%i_sca%i",layer,i,j,isca));
        
	  if(temp->GetEntries()>50) {
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
	      
	      float mean_start_fit=mean_peak_higher;
	      if(npeaks ==0 ) mean_start_fit=temp->GetMean();
              if(npeaks ==1 ) {
                Double_t *mean_peak=s->GetPositionX();
                mean_peak[0]=mean_peak_higher;
                
                TF1 *f0 = new TF1("f0","gaus",mean_peak[0]-temp->GetRMS(),mean_peak[0]+temp->GetRMS());
                temp->Fit("f0","RQNOC");
                
                TF1 *f1 = new TF1("f1","gaus",f0->GetParameter(1)-f0->GetParameter(2),f0->GetParameter(1)+f0->GetParameter(2));
                temp->Fit("f1","RQME");
                fout_ped<<f1->GetParameter(1) << " " << f1->GetParError(1)<<" "<<f1->GetParameter(2)<< " ";
	
		ped_all->Fill(20*layer+i,100*isca+j,f1->GetParameter(1));
		rms_all->Fill(20*layer+i,100*isca+j,f1->GetParameter(2));

		
              } else {
                fout_ped<<mean_peak_higher<< " " << 10 <<" "<<0<<" ";
              }
	    } else {
	      fout_ped<<0<< " " << 0<<" "<<0<<" ";
	    } 
	    delete s;
	  } else {
	    fout_ped<<0<< " " << 0<<" "<<0<<" ";
	  }
	
	
	}//isca
	fout_ped<<endl;
      }//j
    }
  }


  TFile *_file1 = new TFile(TString::Format("plots/PedSummary_%s.root",run.Data()),"RECREATE");
  TCanvas* canvas0= new TCanvas("PedAna_Summary","Pedestal Summary",1800,400);   
  canvas0->Divide(2,1);
  canvas0->cd(1);
  ped_all->GetZaxis()->SetRangeUser(150,350);
  ped_all->Draw("colz");
  canvas0->cd(2);
  rms_all->GetZaxis()->SetRangeUser(0,10);
  rms_all->Draw("colz");
  canvas0->Print(TString::Format("plots/PedSummart_%s.png",run.Data()));
  canvas0->Print(TString::Format("plots/PedSummart_%s.root",run.Data()));
  _file1->cd();
  ped_all->Write();
  rms_all->Write();
  canvas0->Write();
  _file1 ->Close();

  

}


void mipanalysis(TString run="Run_ILC_cosmic_test_11222019"){
  
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

  TFile *_file0 = TFile::Open(TString::Format("../results_calib/%s.root",run.Data()));

  TH2F* MIPM_all=new TH2F("MIPM_all","average of MPVs ; Layer  ; Chip",15,-0.5,14.5,16,-0.5,15.5);
  TH2F* MIPN_all=new TH2F("MIPN_all","channels fitted ; Layer  ; Chip",15,-0.5,14.5,16,-0.5,15.5);
  TH2F* MIPrms_all=new TH2F("MIPrms_all","rms  of MPVs ; Layer  ; Chip",15,-0.5,14.5,16,-0.5,15.5);

  ofstream fout_mip(TString::Format("../results_calib/%s.txt",run.Data()).Data(),ios::out);
  fout_mip<<"#mip results PROTO15-TB2021-11"<<endl;
  fout_mip<<"#layer chip channel mpv empv widthmpv chi2ndf nentries"<<endl;

  for(int layer=0; layer<15; layer++) {

    float avmpv[16]={0};
    float nmpv[16]={0};
    float rmsmpv[16]={0};

    TFile *file = new TFile(TString::Format("plots/%s_layer_%i.root",run.Data(),layer) , "RECREATE");
    file->cd();
    TDirectory *cdhisto = file->GetDirectory("histograms_mip");
    if(!cdhisto) cdhisto = file->mkdir("histograms_mip");
    cdhisto->cd();
    
    // Comparing nbr entries in tag or not tag // GetWidth and Mean
    TString map="../../mapping/fev10_chip_channel_x_y_mapping.txt";
    //  if(layer==0 || layer==2)  map="../../mapping/fev11_cob_chip_channel_x_y_mapping.txt";
    ReadMap(map);
    
    
    
    TH2F* MIPW =new TH2F(TString::Format("MIPW layer%i",layer),TString::Format("MIPW layer%i",layer),32,-90,90,32,-90,90);
    TH2F* MIPM=new TH2F(TString::Format("MIPM, layer%i",layer),TString::Format("MIPM, layer%i",layer),32,-90,90,32,-90,90);
    TH2F* MIPRMS=new TH2F(TString::Format("MIPRMS layer%i",layer),TString::Format("MIPRMS layer%i",layer),32,-90,90,32,-90,90);
    TH2F* MIPN=new TH2F(TString::Format("MIPN layer%i",layer),TString::Format("MIPN layer%i",layer),32,-90,90,32,-90,90);
    
    for(int i=0;i<nchips;i++){

      TCanvas* canvas_mip;

      canvas_mip= new TCanvas(TString::Format("MIPs_layer_%i_chip%i",layer,i),TString::Format("MIPs_layer_%i_chip%i",layer,i),1600,1600);   
      canvas_mip->Divide(8,8);
      
      for(int j=0; j<64; j++) {
	TCanvas* canvastemp= new TCanvas("temp","temp",100,100);   
	canvastemp->cd();
	//GetGoodEntries
	TH1F *temp=(TH1F*)_file0->Get(TString::Format("layer_%i/charge_layer%i_chip%i_chn%i",layer,layer,i,j));
	cout<<TString::Format("layer_%i/charge_layer%i_chip%i_chn%i",layer,layer,i,j)<<endl;
	if(temp==NULL){
	  fout_mip<<layer<<" "<<i<<" "<<j<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<"\n";
	  delete canvastemp; continue;
	}	
	MIPN->Fill(map_pointX[i][j],map_pointY[i][j],temp->GetEntries());
	
	Double_t fr[2];
	Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];   
	pllo[0]=0.5; pllo[1]=10; pllo[2]=1.0; pllo[3]=0;
	plhi[0]=100.0; plhi[1]=300.0; plhi[2]=100000000.0; plhi[3]=50.0;
	sv[0]=15.0;
	Double_t chisqr;
	Int_t    ndf;
	
	if(temp->Integral()>100){
	  fr[0]=TMath::Max(temp->GetMean()-0.75*temp->GetRMS(),10.);
	  fr[1]=temp->GetMean()+2*temp->GetRMS();
	  temp->GetXaxis()->SetRangeUser(fr[0], fr[1]);
	  temp->GetXaxis()->SetRangeUser(fr[0], fr[1]);
	  sv[0] = temp->GetRMS() * 0.25;
	  sv[1] = temp->GetMean() * 0.67;
	  sv[2] = temp->Integral("width") * 1.2;
	  sv[3] = temp->GetRMS()* 0.1;
	  TF1 *fitsnr_temp=langaufit(temp,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
	  for (int k = 0; k < 4; k++) {
	    if ((sv[k] >= pllo[k]) && (sv[k] <= plhi[k])) continue;
	    std::cout << "Langaus fit parameter " << k << " has starting value " << sv[k] 
		      <<  " outside the range [" << pllo[k] << ", " << plhi[k] << "]!" << std::endl;
	  }
	  double mpv=fitsnr_temp->GetParameter(1);
	  double empv=fitsnr_temp->GetParError(1);
	  double wmpv=fitsnr_temp->GetParameter(0);
	  double chi2ndf=0;
	  if(ndf>0) chi2ndf=chisqr/ndf;
	  if(chi2ndf>0 && mpv>10) {
	    avmpv[i]+=mpv;
	    rmsmpv[i]+=mpv*mpv;
	    nmpv[i]++;
	  }
	  fout_mip<<layer<<" "<<i<<" "<<j<<" "<<mpv<<" "<<empv<<" "<<wmpv<<" "<<chi2ndf<<" "<<temp->GetEntries()<<"\n";
	  
	  MIPM->Fill(map_pointX[i][j],map_pointY[i][j],mpv);
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
      canvas_mip->Print(TString::Format("plots/%s_layer_%i_chip%i.png",run.Data(),layer,i));
      canvas_mip->Write();
      delete canvas_mip;
    }
    for(int i=0;i<nchips;i++){
      if(nmpv[i]!=0 ) {
	avmpv[i]/=nmpv[i];
	rmsmpv[i]=sqrt( (rmsmpv[i]-nmpv[i]*avmpv[i]*avmpv[i]) / (nmpv[i]-1) );
	MIPM_all->Fill(layer,i,avmpv[i]);
	MIPrms_all->Fill(layer,i,rmsmpv[i]);
	MIPN_all->Fill(layer,i,nmpv[i]);
      }

   
    }


    TCanvas* canvas= new TCanvas(TString::Format("MIPAna_layer_%i",layer),TString::Format("MIPAna_layer_%i",layer),1600,1600);   
    canvas->Divide(2,2);
  
    canvas->cd(1);
    MIPM->GetXaxis()->SetTitle("x");
    MIPM->GetYaxis()->SetTitle("y");
    MIPM->GetZaxis()->SetRangeUser(10,150);
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
    MIPN->GetZaxis()->SetRangeUser(100,100000);
    MIPN->Draw("COLZ");

    canvas->cd(4);
    MIPRMS->GetXaxis()->SetTitle("x");
    MIPRMS->GetYaxis()->SetTitle("y");
    MIPRMS->GetZaxis()->SetRangeUser(5,200);
    MIPRMS->Draw("COLZ");

    canvas->Print(TString::Format("plots/MIPAna_%s_layer_%i.png",run.Data(),layer));

    file->cd();
    MIPM->Write();
    MIPW->Write();
    MIPN->Write();
    MIPRMS->Write();
    canvas->Write();
    canvas->Print(TString::Format("MIPAna_layer_%i.png",layer));

  }

  
    _file0 ->Close();

  TCanvas* canvassummary= new TCanvas("MIPAna","MIPAna",1600,800);   
  canvassummary->Divide(3,1);
    
  canvassummary->cd(1);
  // MIPM_all->GetXaxis()->SetTitle("x");
  // MIPM_all->GetYaxis()->SetTitle("y");
  MIPM_all->GetZaxis()->SetRangeUser(10,150);
  MIPM_all->Draw("colz");

  canvassummary->cd(2);
  // MIPrms_all->GetXaxis()->SetTitle("x");
  // MIPrms_all->GetYaxis()->SetTitle("y");
  MIPrms_all->GetZaxis()->SetRangeUser(1,50);
  MIPrms_all->Draw("COLZ");
    
  canvassummary->cd(3);
  //gPad->SetLogz();
  // MIPN_all->GetXaxis()->SetTitle("x");
  // MIPN_all->GetYaxis()->SetTitle("y");
  MIPN_all->GetZaxis()->SetRangeUser(0,64.5);
  MIPN_all->Draw("COLZ");
  canvassummary->Print("summary.root");
  canvassummary->Print("summary.png");
      
}



 void mipanalysis_summary(TString run="Run_ILC_cosmic_test_11222019"){ 
  
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


  TH2F* MIPM_all=new TH2F("MIPM_all","average of MPVs ; Layer  ; Chip",15,-0.5,14.5,16,-0.5,15.5);
  TH2F* MIPN_all=new TH2F("MIPN_all","channels fitted ; Layer  ; Chip",15,-0.5,14.5,16,-0.5,15.5);
  TH2F* MIPrms_all=new TH2F("MIPrms_all","rms  of MPVs ; Layer  ; Chip",15,-0.5,14.5,16,-0.5,15.5);

  ofstream fout_mip(TString::Format("../results_calib/%s.txt",run.Data()).Data(),ios::out);
  fout_mip<<"#mip results PROTO15-TB2021-11"<<endl;
  fout_mip<<"#layer chip channel mpv empv widthmpv chi2ndf nentries"<<endl;

  for(int layer=0; layer<15; layer++) {

    TFile *_file0 = TFile::Open(TString::Format("../results_calib/%s.root",run.Data()));
     
    float avmpv[16]={0};
    float nmpv[16]={0};
    float rmsmpv[16]={0};

    TString map="../../mapping/fev10_chip_channel_x_y_mapping.txt";
    ReadMap(map);

    cout<<"   LAYER : "<<layer<<endl;
        
    for(int i=0;i<nchips;i++){

        
      for(int j=0; j<64; j++) {
	TCanvas* canvastemp= new TCanvas("temp","temp",100,100);
	canvastemp->cd();
	//GetGoodEntries
	TH1F *temp=(TH1F*)_file0->Get(TString::Format("layer_%i/charge_layer%i_chip%i_chn%i",layer,layer,i,j));
	if(temp==NULL){
	  fout_mip<<layer<<" "<<i<<" "<<j<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<"\n";
	  delete canvastemp; continue;
	}
	//	MIPN->Fill(map_pointX[i][j],map_pointY[i][j],temp->GetEntries());
	
	Double_t fr[2];
	Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
	pllo[0]=0.5; pllo[1]=10; pllo[2]=1.0; pllo[3]=0;
	plhi[0]=100.0; plhi[1]=300.0; plhi[2]=100000000.0; plhi[3]=50.0;
	sv[0]=15.0;
	Double_t chisqr;
	Int_t    ndf;
	
	if(temp->Integral()>100){
	  fr[0]=TMath::Max(temp->GetMean()-0.75*temp->GetRMS(),10.);
	  fr[1]=temp->GetMean()+2*temp->GetRMS();
	  temp->GetXaxis()->SetRangeUser(fr[0], fr[1]);
	  temp->GetXaxis()->SetRangeUser(fr[0], fr[1]);
	  sv[0] = temp->GetRMS() * 0.25;
	  sv[1] = temp->GetMean() * 0.67;
	  sv[2] = temp->Integral("width") * 1.2;
	  sv[3] = temp->GetRMS()* 0.1;
	  TF1 *fitsnr_temp=langaufit(temp,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
	  for (int k = 0; k < 4; k++) {
	    if ((sv[k] >= pllo[k]) && (sv[k] <= plhi[k])) continue;
	    std::cout << "Langaus fit parameter " << k << " has starting value " << sv[k]
		      <<  " outside the range [" << pllo[k] << ", " << plhi[k] << "]!" << std::endl;
	  }
	  double mpv=fitsnr_temp->GetParameter(1);
	  double empv=fitsnr_temp->GetParError(1);
	  double wmpv=fitsnr_temp->GetParameter(0);
	  double chi2ndf=0;
	  if(ndf>0) chi2ndf=chisqr/ndf;
	  if(chi2ndf>0 && mpv>10) {
	    avmpv[i]+=mpv;
	    rmsmpv[i]+=mpv*mpv;
	    nmpv[i]++;
	  }
	  fout_mip<<layer<<" "<<i<<" "<<j<<" "<<mpv<<" "<<empv<<" "<<wmpv<<" "<<chi2ndf<<" "<<temp->GetEntries()<<"\n";
	  
	
	} else {
	  fout_mip<<layer<<" "<<i<<" "<<j<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<"\n";
	}
	delete canvastemp; 
      }
      //  file->cd();
      //canvas_mip->Print(TString::Format("plots/%s_layer_%i_chip%i.png",run.Data(),layer,i));
      //  canvas_mip->Write();
      //delete canvas_mip;
    }
    for(int i=0;i<nchips;i++){
      if(nmpv[i]!=0 ) {
	avmpv[i]/=nmpv[i];
	rmsmpv[i]=sqrt( (rmsmpv[i]-nmpv[i]*avmpv[i]*avmpv[i]) / (nmpv[i]-1) );
	MIPM_all->Fill(layer,i,avmpv[i]);
	MIPrms_all->Fill(layer,i,rmsmpv[i]);
	MIPN_all->Fill(layer,i,nmpv[i]);
      }

   
    }


  }

 
  TCanvas* canvassummary= new TCanvas("MIPAna","MIPAna",1600,800);
  canvassummary->Divide(3,1);
    
  canvassummary->cd(1);
  // MIPM_all->GetXaxis()->SetTitle("x");
  // MIPM_all->GetYaxis()->SetTitle("y");
  MIPM_all->GetZaxis()->SetRangeUser(10,150);
  MIPM_all->Draw("colz");

  canvassummary->cd(2);
  // MIPrms_all->GetXaxis()->SetTitle("x");
  // MIPrms_all->GetYaxis()->SetTitle("y");
  MIPrms_all->GetZaxis()->SetRangeUser(1,50);
  MIPrms_all->Draw("COLZ");
    
  canvassummary->cd(3);
  //gPad->SetLogz();
  // MIPN_all->GetXaxis()->SetTitle("x");
  // MIPN_all->GetYaxis()->SetTitle("y");
  MIPN_all->GetZaxis()->SetRangeUser(0,64.5);
  MIPN_all->Draw("COLZ");
  canvassummary->Print("summary.root");
  canvassummary->Print("summary.png");
      
}





void retriggeranalysis(TFile* file, TString run="Run_ILC_cosmic_test_11222019", int layer=2){
  
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
  TString map="../../mapping/fev10_chip_channel_x_y_mapping.txt";
  ReadMap(map);

  
  TFile *_file0 = TFile::Open(TString::Format("../results_retriggers/Retriggers_layer_%i_%s.root",layer,run.Data()));
   
  TH2F* f_ret_xy= (TH2F*)_file0->Get("first_retriggering_xy");
  TH2F* ret_xy= (TH2F*)_file0->Get("all_retriggering_xy");
  TH2F* trig_xy= (TH2F*)_file0->Get("trig_xy");

  TH2F* f_ret= (TH2F*)_file0->Get("first_retriggering");
  TH2F* ret= (TH2F*)_file0->Get("all_retriggering");
  TH2F* trig= (TH2F*)_file0->Get("trig");

  TCanvas* canvas= new TCanvas(TString::Format("Ret_layer_%i",layer),TString::Format("Ret_layer_%i",layer),1200,1600);   
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

  canvas->Print(TString::Format("plots/Ret_%s_layer_%i.png",run.Data(),layer));
  
  file->cd();
  trig_xy->Write();
  trig->Write();
  f_ret_xy->Write();
  f_ret->Write();
  ret_xy->Write();
  ret->Write();
  canvas->Write();
  
}
