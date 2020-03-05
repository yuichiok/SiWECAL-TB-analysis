//# Copyright 2020  Adrián Irles

#include "TROOT.h"
#include "TFile.h"
#include "DecodedSLBAnalysis.cc"
#include "TGraphErrors.h"
#include "TLegend.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"

void FitLinearity4Channels(TString run="calib_02272020_pa1.2fb_trig230_chn0.8.16.24_calib_small_Ascii", int nslboards=6, int channel1=0, int channel2=8, int channel3=16, int channel4=24){
  
  cout<<" Graphs file: "<<run<<endl;
 
  TGraphErrors* injection_high[15][16][64][3];
  TGraphErrors* injection_low[15][16][64][3];

  TFile *file_read = new TFile("results/graphs_"+run+".root" , "READ");
  file_read->cd();
  
  for(int i=0; i<nslboards; i++) {
    for(int j=0; j<16; j++) {
      for(int k=0; k<64; k++) {
	
	injection_high[i][j][k][0]=(TGraphErrors*)file_read->Get(TString::Format("high_gain_layer%i_chip%i_chn%i_sca0",i,j,k));
	injection_high[i][j][k][1]=(TGraphErrors*)file_read->Get(TString::Format("high_gain_layer%i_chip%i_chn%i_sca1",i,j,k));
	injection_high[i][j][k][2]=(TGraphErrors*)file_read->Get(TString::Format("high_gain_layer%i_chip%i_chn%i_sca2",i,j,k));

	injection_low[i][j][k][0]=(TGraphErrors*)file_read->Get(TString::Format("low_gain_layer%i_chip%i_chn%i_sca0",i,j,k));
	injection_low[i][j][k][1]=(TGraphErrors*)file_read->Get(TString::Format("low_gain_layer%i_chip%i_chn%i_sca1",i,j,k));
	injection_low[i][j][k][2]=(TGraphErrors*)file_read->Get(TString::Format("low_gain_layer%i_chip%i_chn%i_sca2",i,j,k));
      }
    }
  }

  file_read->Close();

  TFile *rejection_write = new TFile("results/rejectedFits_"+run+".root", "RECREATE");
  
  TFile *file_write = new TFile("results/fits_"+run+".root" , "RECREATE");
  file_write->cd();
  
  for(int i=0; i<nslboards; i++) {

    TString fileSuffix;
    fileSuffix.Form("layer%i_" + run + ".txt", i);
   
    fstream high_pedestalsFile;
    high_pedestalsFile.open("results/pedestals/pedestals_hg_" + fileSuffix,fstream::out);
    
    fstream low_pedestalsFile;
    low_pedestalsFile.open("results/pedestals/pedestals_lg_"+ fileSuffix,fstream::out);
    
    TString header = "#pedestal results (injection fits): " + run + "\n";
    header += "#chip channel ped0 eped0 accept0 ped1 eped1 accept1 ... pedn epedn acceptn (All sca)\n";
    
    high_pedestalsFile << header;
    low_pedestalsFile << header;
    
    for(int j=0; j<16; j++) {
      
      vector<TCanvas*> vcanvas;
      vector<Bool_t> rej_vcanvas_hg;
      vector<Bool_t> rej_vcanvas_lg;
      
      for(int ichan=0; ichan<4; ichan++) {
	
	int chan=0;
	if(ichan==0) chan=channel1;
	if(ichan==1) chan=channel2;
	if(ichan==2) chan=channel3;
	if(ichan==3) chan=channel4;
	
	high_pedestalsFile << j << " " << chan << " ";
	low_pedestalsFile << j << " " << chan << " ";
	for(int isca=0; isca<15; isca++) {
	  
	  if(!ichan)
	    {    
	      TCanvas* canvas = new TCanvas(TString::Format("canvas_layer%i_chip%i_sca%i",i,j,isca),TString::Format("canvas_layer%i_chip%i_sca%i",i,j,isca),1200,1000);
	      canvas->Divide(2,2);
	      vcanvas.push_back(canvas);
	      rej_vcanvas_hg.push_back(false);
	      rej_vcanvas_lg.push_back(false);
	    }

	  Bool_t acceptFit_hg = true;
	  Bool_t acceptFit_lg = true;
	  
	  TFitResultPtr fit_result_high = injection_high[i][j][chan][isca]->Fit("pol1", "QS", "", 0., 3.4);
	  TFitResultPtr fit_result_low = injection_low[i][j][chan][isca]->Fit("pol1", "QS", "", 0., 3.4);

	  TF1* fit_func_high = (TF1*)injection_high[i][j][chan][isca]->GetListOfFunctions()->FindObject("pol1");
	  TF1* fit_func_low = (TF1*)injection_low[i][j][chan][isca]->GetListOfFunctions()->FindObject("pol1");
	    
	  Double_t* fit_high = fit_func_high->GetParameters();
	  const Double_t* fit_high_err = fit_func_high->GetParErrors();
	  Double_t chi_ndf_high = fit_func_high->GetChisquare()/fit_func_high->GetNDF();

	  if(chi_ndf_high > 2.5 || fit_func_high->GetNumberFitPoints() < 5 || !fit_result_high->IsValid())
	    {
	      acceptFit_hg = false;
	      rej_vcanvas_hg[isca] = !acceptFit_hg;
	    } 
	  
	  Double_t* fit_low = fit_func_low->GetParameters();
	  const Double_t* fit_low_err = fit_func_low->GetParErrors();
	  Double_t chi_ndf_low = fit_func_low->GetChisquare()/fit_func_low->GetNDF();

	   if(chi_ndf_low > 2.5 || fit_func_low->GetNumberFitPoints() < 5 || !fit_result_low->IsValid())
	    {
	      acceptFit_lg = false;
	      rej_vcanvas_lg[isca] = !acceptFit_lg;
	    }
	  
	  high_pedestalsFile << fit_high[0] << " " << fit_high_err[0] << " " << acceptFit_hg << " ";
	  low_pedestalsFile << fit_low[0] << " " << fit_low_err[0] << " " << acceptFit_lg << " ";

	  vcanvas[isca]->cd(1+ichan);
	  TLegend *leg = new TLegend(0.1,0.5,0.25,0.85);
	  injection_high[i][j][chan][isca]->SetLineColor(1);
	  injection_high[i][j][chan][isca]->SetLineStyle(1);
	  injection_high[i][j][chan][isca]->SetMarkerColor(1);
	  injection_high[i][j][chan][isca]->SetMarkerStyle(21);
	  injection_high[i][j][chan][isca]->SetLineWidth(2);
	  injection_high[i][j][chan][isca]->GetYaxis()->SetRangeUser(100,500);

	  injection_high[i][j][chan][isca]->SetTitle(TString::Format("Layer %i, Chip %i, channel: %i, SCA %i",i,j,chan,isca));
	  injection_high[i][j][chan][isca]->DrawClone("alp");

	  TString legendTitle = "High Gain ";
	  if(!acceptFit_hg) legendTitle += "(Rejected) "; 
	  leg->AddEntry(injection_high[i][j][chan][isca], legendTitle,"lep");

	  injection_low[i][j][chan][isca]->SetLineColor(2);
	  injection_low[i][j][chan][isca]->SetLineStyle(2);
	  injection_low[i][j][chan][isca]->SetMarkerColor(2);
	  injection_low[i][j][chan][isca]->SetMarkerStyle(22);
	  injection_low[i][j][chan][isca]->SetLineWidth(2);
	  injection_low[i][j][chan][isca]->DrawClone("lp");
	  
	  legendTitle = "Low Gain ";
	  if(!acceptFit_lg) legendTitle += "(Rejected) "; 
	  leg->AddEntry(injection_low[i][j][chan][isca], legendTitle,"lep");
	  leg->DrawClone();
	}

	high_pedestalsFile << endl;
	low_pedestalsFile << endl;
	
      }

      for(int iCanvas = 0; iCanvas < vcanvas.size(); iCanvas++)	{
	file_write->cd();
	vcanvas[iCanvas]->Write();
	
	if(rej_vcanvas_lg[iCanvas] || rej_vcanvas_hg[iCanvas])
	  {
	    rejection_write->cd();
	    vcanvas[iCanvas]->Write();
	  }
	delete vcanvas[iCanvas];
      }

    }

    high_pedestalsFile.close();
    low_pedestalsFile.close();
    
  }
}

