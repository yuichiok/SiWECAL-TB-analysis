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
#include "../../include/utils.h"

//#include "../../style/Style.C"
//#include "../../style/Labels.C"

using namespace std;

int nchips=16;
int nlayers=15;

void pedanalysis(TString run="PedestalMIP_3GeVMIPscan", TString gain="low"){
  
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

  ofstream fout_ped(TString::Format("../results_calib/%s_%sgain.txt",run.Data(),gain.Data()).Data(),ios::out);
  fout_ped<<"#pedestal results (fit to a gaussian) remove channels/sca with two pedestals peaks from the analysis : PROTO15"<<endl;
  fout_ped<<"#layer chip channel ped0 eped0 widthped0 ped1 eped1 widthped1... ped14 eped14 widhtped14 (all SCA)"<<endl;

  for(int layer=0; layer<nlayers; layer++) {

    
    // Comparing nbr entries in tag or not tag // GetWidth and Mean
    //TString map="../../mapping/fev10_chip_channel_x_y_mapping.txt";
    //  if(layer==0 || layer==2)  map="../../mapping/fev11_cob_chip_channel_x_y_mapping.txt";
    //  ReadMap(map);

 
    for(int i=0;i<nchips;i++){
    
      double nch=0;
      float avmean=0, avwidth=0;
    
      for(int j=0; j<64; j++) {
	fout_ped << layer<<" "<<i <<" " <<j<< " ";

	double ped[15]={0};
	double pederror[15]={0};
	double pedwidth[15]={0};
	double pedwidtherror[15]={0};


      
	for(int isca=0; isca<15; isca++) {
	  TH1F* temp=(TH1F*)_file0->Get(TString::Format("layer_%i/ped_%s_chip%i_chn%i_sca%i",layer,gain.Data(),i,j,isca));
        
	  if(temp->GetEntries()>500) {
	    TSpectrum *s = new TSpectrum();
	    int npeaks = s->Search(temp,2,"",0.8);
	    if(npeaks == 1) {
	      Double_t *mean_peak=s->GetPositionX();
	      Double_t *mean_high=s->GetPositionY();
	      double mean_peak_higher=0;
	      double mean_high_higher=0;
	      int npeak_max=0;
              for(int ipeak=0; ipeak<npeaks; ipeak ++) {
                if(mean_high[ipeak]>mean_high_higher && mean_peak[ipeak]>180) {
                  mean_high_higher=mean_high[ipeak];
                  mean_peak_higher=mean_peak[ipeak];
                  npeak_max=ipeak;
                }
              }
	      	   
	      float mean_start_fit=mean_peak_higher;
	      
	      TF1 *f0 = new TF1("f0","gaus",mean_peak_higher-2,mean_peak_higher+2);
	      TF1 *f1 = new TF1("f1","gaus",mean_peak_higher-20,mean_peak_higher+20);
	      temp->Fit("f0","RQE");
	      temp->Fit("f1","RQE");

	      if(f0->GetParameter(2)<f1->GetParameter(2))  {
		ped[isca]=f0->GetParameter(1);
		pedwidth[isca]=f0->GetParameter(2);
		pederror[isca]=f0->GetParError(1);
		pedwidtherror[isca]=f0->GetParError(2);
	      }  else  {
		ped[isca]=f1->GetParameter(1);
		pedwidth[isca]=f1->GetParameter(2);
		pederror[isca]=f1->GetParError(1);
		pedwidtherror[isca]=f1->GetParError(2);
	      }
	      
	      delete s;
	    }
	  } else pedwidth[isca]=-1;
	  
	}//isca

	double ped_av=weigthedAv(ped, pederror, 15);
	double pedwidth_av=weigthedAv(pedwidth, pedwidtherror, 15);
	for(int isca=0; isca<15; isca++) {
	  
	  if(pedwidth[isca]>0.5 && pedwidth[isca]<8 ) {
	    fout_ped<<ped[isca]<< " " << pederror[isca]<<" "<<pedwidth[isca]<<" "; //no peaks found
	    ped_all->Fill(20*layer+i,100*isca+j,ped[isca]);
	    rms_all->Fill(20*layer+i,100*isca+j,pedwidth[isca]);
	  } else  {
	    if(pedwidth[isca]>-0.5) {
	      fout_ped<<ped_av<< " " << -5<<" "<<pedwidth_av<<" "; //no peaks found
	      ped_all->Fill(20*layer+i,100*isca+j,ped_av);
	      rms_all->Fill(20*layer+i,100*isca+j,pedwidth_av);
	    }
	    if(pedwidth[isca]<0) {
	      fout_ped<<-5<< " " << -5<<" "<<-5<<" "; //no peaks found
	      ped_all->Fill(20*layer+i,100*isca+j,100);
	      rms_all->Fill(20*layer+i,100*isca+j,100);
	    }
	  }
	}
	fout_ped<<endl;
      }//j
    }
  }


  TFile *_file1 = new TFile(TString::Format("%s_%sgain_PedSummary.root",run.Data(),gain.Data()),"RECREATE");
  TCanvas* canvas0= new TCanvas("PedAna_Summary","Pedestal Summary",1800,400);   
  canvas0->Divide(2,1);
  canvas0->cd(1);
  ped_all->GetZaxis()->SetRangeUser(180,300);
  ped_all->Draw("colz");
  canvas0->cd(2);
  rms_all->GetZaxis()->SetRangeUser(1.0,6.0);
  rms_all->Draw("colz");
  canvas0->Print(TString::Format("%s_%sgain_PedSummary.png",run.Data(),gain.Data()));
  _file1->cd();
  ped_all->Write();
  rms_all->Write();
  canvas0->Write();
  _file1 ->Close();

  

}


