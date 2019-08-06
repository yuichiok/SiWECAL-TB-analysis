#include <iomanip>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <bitset>
#include <math.h>
#include <vector>
#include <string>
#include <algorithm>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <sstream>
#include <string>


int time_correlation_all_files() {
  TFile *file_writing = new TFile("all_runs_coincidences_studies.root", "recreate");
  
  Int_t NFILE=4;
  TFile *file_histo[NFILE];

  Int_t NSLAB=9;
  Int_t NCHIP=4;
  Double_t N_SLAB[NCHIP][NSLAB];
  Double_t X[NFILE][NCHIP][NSLAB];
  Double_t Y[NFILE][NCHIP][NSLAB];
  Double_t e_X[NFILE][NCHIP][NSLAB];
  Double_t err_null[NFILE][NCHIP][NSLAB];
  TGraphErrors *g_bcid_sep[NFILE][NCHIP];
  TCanvas *canvas[NCHIP];
  // TCanvas *canvas2[NCHIP];
  TMultiGraph *mg[NCHIP];
  
  for(Int_t chip=0; chip<NCHIP; chip++) {
    mg[chip]= new TMultiGraph();

    canvas[chip]=new TCanvas(TString::Format("slab_5_chip_%i_time_correl_all_files", chip),TString::Format("slab_5_chip%i_time_correl_all_files", chip), 1800,1400);
    // canvas2[chip]=new TCanvas(TString::Format("slab_5_chip_%i_time_correl", chip),TString::Format("slab_5_chip%i_time_correl", chip), 1800,1400);
    // canvas2[chip]->Divide(2,2);
    
    for (Int_t file=0; file<NFILE; file++) {
      TString name = "";
      if (file==0) {
	name="32001";
      }
      if (file==1) {
	name="32002";
      }
      if (file==2) {
	name="32004";
      }
      if (file==3) {
	name="32015";
      }
      TString filename="run_";
      filename += name;
      filename +="_histograms.root";
      
      file_histo[file] = new TFile(filename, "read");
      
      
      //TCanvas* canvas2[NCHIP];
    
            
      //canvas[chip]->Divide(1,2);
      for (Int_t slab=0; slab<NSLAB; slab++) {
	err_null[file][chip][slab]=0;
	N_SLAB[chip][slab]=slab;
	Int_t chip0=chip;
	if (slab<5) {
	  chip0=15-chip;
	}
	if (slab>4) {
	  chip0=chip;
	}
	TH1F *histo_1D = (TH1F*)file_histo[file]->Get(TString::Format("number_coincidences_slab_5_chip_%i_slab_%i_chip_%i", chip, slab, chip0));
	cout<<TString::Format("number_coincidences_slab_5_chip_%i_slab_%i_chip_%i", chip, slab, chip0)<<endl;
	TSpectrum *s = new TSpectrum(1);
	histo_1D->GetXaxis()->SetRangeUser(-200,200);
	int npeaks = s->Search(histo_1D, 2, "", 0.8);
	
	Double_t *mean_X = s->GetPositionX();
	Double_t *mean_Y = s->GetPositionY();
	X[file][chip][slab]=mean_X[0];
	Y[file][chip][slab]=mean_Y[0];
	histo_1D->GetXaxis()->SetRangeUser(X[file][chip][slab]-3, X[file][chip][slab]+3);
	e_X[file][chip][slab]=histo_1D->GetRMS();
      
	// graph_bcid_separation[chip]->Fill();
	// coincidences[chip];
      }
      g_bcid_sep[file][chip]=new TGraphErrors(NSLAB,N_SLAB[chip],X[file][chip],0, e_X[file][chip]);
      


      // canvas2[chip]->cd(file+1);
      // g_bcid_sep[file][chip]->SetLineColor(0);
      // g_bcid_sep[file][chip]->SetMarkerColor(2);
      // g_bcid_sep[file][chip]->SetMarkerStyle(3);
      // g_bcid_sep[file][chip]->GetXaxis()->SetTitle("slab");
      // g_bcid_sep[file][chip]->GetYaxis()->SetTitle("bcid_sep");
      // g_bcid_sep[file][chip]->Draw("alp");

      g_bcid_sep[file][chip]->SetLineColor(2*(file+1));
      g_bcid_sep[file][chip]->SetMarkerColor(2*(file+1));
      g_bcid_sep[file][chip]->SetMarkerStyle(file+1);
      g_bcid_sep[file][chip]->SetMarkerSize(2);
      g_bcid_sep[file][chip]->Draw("alp");
      mg[chip]->Add(g_bcid_sep[file][chip]);
      mg[chip]->Draw("*");

    }
  
    
    
    canvas[chip];//->cd(2);

    // mg[chip]->GetXaxis()->SetTitle("slab");
    // mg[chip]->GetYaxis()->SetTitle("bcid_sep");
    //mg[chip]->SetTitle(TString::Format("bcid separation with slab 5 at chip %i", chip));
    mg[chip]->Draw("alp");
    
    mg[chip]->SetName(TString::Format("bcid_sep_slab_5_chip_%i", chip));
    file_writing->cd();

    mg[chip]->Write();//"all_runs_coincidences_studies.root");

    // g_bcid_sep[chip]->Write();
  }
  file_writing->Close();
  return 1;
  
}
