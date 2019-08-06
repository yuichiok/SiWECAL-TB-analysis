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


int time_correlation(TString run="32001") {
  
  

  TFile *file_histo = new TFile("run_"+run+"_histograms.root", "read");
  cout<<"run_"+run+"_histograms.root"<<endl;
  TFile *file_writing = new TFile("run_"+run+"_coincidences_studies.root", "update");

    for (int ref=5; ref<9; ref++ ) {

  Int_t NSLAB=9;
  Int_t NCHIP=4;
  Double_t N_SLAB[NCHIP][NSLAB];
  Double_t X[NCHIP][NSLAB];
  Double_t Y[NCHIP][NSLAB];
  Double_t e_X[NCHIP][NSLAB];
  Double_t err_null[NCHIP][NSLAB];
  TGraphErrors *g_coinc[NCHIP];
  TGraphErrors *g_bcid_sep[NCHIP];
  TCanvas* canvas[NCHIP];
  //TCanvas* canvas2[NCHIP];
  
  for(Int_t chip=0; chip<NCHIP; chip++) {

    canvas[chip]=new TCanvas(TString::Format("slab_%i_chip_%i_time_correl", ref,chip),TString::Format("slab_%i_chip%i_time_correl", ref,chip), 1800,1400);

    //canvas[chip]->Divide(1,2);
    for (Int_t slab=0; slab<NSLAB; slab++) {
      err_null[chip][slab]=0;
      N_SLAB[chip][slab]=slab;
      Int_t chip0=chip;
      if (slab<5) {
    	chip0=15-chip;
      }
      if (slab>4) {
    	chip0=chip;
      }
      TH1F *histo_1D = (TH1F*)file_histo->Get(TString::Format("number_coincidences_slab_%i_chip_%i_slab_%i_chip_%i", ref,chip, slab, chip0));
      TSpectrum *s = new TSpectrum(1);
      histo_1D->GetXaxis()->SetRangeUser(-200,200);
      int npeaks = s->Search(histo_1D, 2, "", 0.8);

      Double_t *mean_X = s->GetPositionX();
      Double_t *mean_Y = s->GetPositionY();
      X[chip][slab]=mean_X[0];
      Y[chip][slab]=mean_Y[0];
      histo_1D->GetXaxis()->SetRangeUser(X[chip][slab]-3, X[chip][slab]+3);
      e_X[chip][slab]=histo_1D->GetRMS();
      
      // graph_bcid_separation[chip]->Fill();
      // coincidences[chip];
    }
    
    g_coinc[chip]=new TGraphErrors(NSLAB,N_SLAB[chip],Y[chip],0, 0);
    g_bcid_sep[chip]=new TGraphErrors(NSLAB,N_SLAB[chip],X[chip],err_null[chip], e_X[chip]);

    canvas[chip];//->cd(2);
    g_bcid_sep[chip]->SetMarkerColor(1);
    g_bcid_sep[chip]->SetMarkerStyle(1);
    g_bcid_sep[chip]->GetXaxis()->SetTitle("slab");
    g_bcid_sep[chip]->GetYaxis()->SetTitle("bcid_sep");
    g_bcid_sep[chip]->SetTitle(TString::Format("bcid separation with slab %i",ref));
    g_bcid_sep[chip]->Draw("alp");
    
    // canvas[chip+4];//->cd(1);
    // g_coinc[chip]->SetMarkerColor(1);
    // g_coinc[chip]->SetMarkerStyle(1);
    // g_coinc[chip]->GetXaxis()->SetTitle("slab");
    // g_coinc[chip]->GetYaxis()->SetTitle("nb_coincidences");
    // g_coinc[chip]->Draw("alp");

    g_bcid_sep[chip]->SetName(TString::Format("slab_%i_chip_%i_time_correl", ref, chip));
    
    g_bcid_sep[chip]->Write();
  }
    }
  file_writing->Close();
  return 1;
  
}
