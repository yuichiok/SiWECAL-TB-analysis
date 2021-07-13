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

void PedestalFile(TString rootfile="") {
  int nSLB=15;

  ofstream fout_ped(rootfile+".txt",ios::out);
  fout_ped<<"#pedestal results (fit to a gaussian) remove channels/sca with two pedestals peaks from the analysis : PROTO15"<<endl;
  fout_ped<<"#layer chip channel ped0 eped0 widthped0 ped1 eped1 widthped1... ped14 eped14 widhtped14 (all SCA)"<<endl;

  TFile *_file0 = TFile::Open(rootfile+".root");
  for(int ilayer=0; ilayer<nSLB; ilayer++) {
    for(int ichip=0; ichip<16; ichip++) {
      for(int ichn=0;ichn<64;ichn++) {
	fout_ped << ilayer<<" "<<ichip <<" " <<ichn<< " ";

	for(int isca=0; isca<15; isca++) {
	  TH1F* htemp=(TH1F*)_file0->Get(TString::Format("layer_%i/ped_chip%i_chn%i_sca%i",ilayer,ichip,ichn,isca));
	  
	  if(htemp->GetEntries()>50) {
	    TSpectrum *s = new TSpectrum();
	    int npeaks = s->Search(htemp,2,"",0.8);
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
	      if(npeaks ==1 ) {
		Double_t *mean_peak=s->GetPositionX();
		mean_peak[0]=mean_peak_higher;
		
		TF1 *f0 = new TF1("f0","gaus",mean_peak[0]-htemp->GetRMS(),mean_peak[0]+htemp->GetRMS());
		htemp->Fit("f0","RQNOC");
		
		TF1 *f1 = new TF1("f1","gaus",f0->GetParameter(1)-f0->GetParameter(2),f0->GetParameter(1)+f0->GetParameter(2));
		htemp->Fit("f1","RQME");
		fout_ped<<f1->GetParameter(1) << " " << f1->GetParError(1)<<" "<<f1->GetParameter(2)<< " ";
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
      }//ichn
    }//ichip
  }//ilayer

}

