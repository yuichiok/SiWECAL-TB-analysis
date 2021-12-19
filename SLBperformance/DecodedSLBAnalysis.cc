//# Copyright 2020 Adrián Irles IJCLab (CNRS/IN2P3)
#define DecodedSLBAnalysis_cxx
#include "DecodedSLBAnalysis.h"
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


void DecodedSLBAnalysis::HitMonitoring(TString outputname="", int maxnhit=5, int nslabshit=8)
{

  //Read the channel/chip -- x/y mapping
  //  ReadMap(map_filename);

  int nSLB=15;// get it as argument or read it from ntuple or somethng ? 
  if (fChain == 0) return -1;

  Long64_t nentries = fChain->GetEntriesFast();
  //  nentries/=10;
  
  TH2F *hit_monitoring = new TH2F("hitmonitoring","hitmonitoring",18,0.5,18.5,300,-0.5,299.5);
  TH2F *hit_monitoring_layer = new TH2F("hitmonitoring_layer","hitmonitoring_layer",18,0.5,18.5,30,-0.5,29.5);

  TH2F *event_monitoring = new TH2F("eventmonitoring","eventmonitoring",18,0.5,18.5,300,-0.5,299.5);
  TH2F *event_monitoring_layer = new TH2F("eventmonitoring_layer","eventmonitoring_layer",18,0.5,18.5,30,-0.5,29.5);

  TH2F *hit_monitoring_norettrains = new TH2F("hitmonitoring_norettrains","hitmonitoring_norettrains",18,0.5,18.5,300,-0.5,299.5);
  TH2F *hit_monitoring_norettrains_layer = new TH2F("hitmonitoring_norettrains_layer","hitmonitoring_norettrains_layer",18,0.5,18.5,30,-0.5,29.5);

  TH2F *event_monitoring_norettrains = new TH2F("eventmonitoring_norettrains","eventmonitoring_norettrains",18,0.5,18.5,300,-0.5,299.5);
  TH2F *event_monitoring_norettrains_layer = new TH2F("eventmonitoring_norettrains_layer","eventmonitoring_norettrains_layer",18,0.5,18.5,30,-0.5,29.5);

  
  Long64_t nbytes = 0, nb = 0;

  // -----------------------------------------------------------------------------------------------------   
  // Basic Pedestal/MIP analysis
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
      
    if(n_slboards>nSLB) nSLB=n_slboards;
    
    if ( jentry > 100 && jentry % 100 ==0 ) 
      printProgress(100.*jentry/nentries) ;
    for(int ilayer=0; ilayer<n_slboards; ilayer++) {

      hit_monitoring_norettrains_layer->Fill(3,(14-ilayer));
      event_monitoring_norettrains_layer->Fill(3,(14-ilayer));
	
      hit_monitoring_layer->Fill(3,(14-ilayer));
      event_monitoring_layer->Fill(3,(14-ilayer));

      for(int ichip=0; ichip<16; ichip++) {
	hit_monitoring->Fill(3,(14-ilayer)*20. +ichip);
        event_monitoring->Fill(3,(14-ilayer)*20. +ichip);

	hit_monitoring_norettrains->Fill(3,(14-ilayer)*20. +ichip);
        event_monitoring_norettrains->Fill(3,(14-ilayer)*20. +ichip);
	
	int last_sca=0;
	bool retrigger=false;
	for(int isca=0; isca<15; isca++) {
	  if(bcid[ilayer][ichip][isca]>-1) last_sca=isca;

	  if( bcid[ilayer][ichip][isca]<50 || (bcid[ilayer][ichip][isca]>900 && bcid[ilayer][ichip][isca]<1000) || badbcid[ilayer][ichip][isca]==1 || badbcid[ilayer][ichip][isca]==2 ) continue;

	  int bcid_seen = SimpleCoincidenceTagger(ilayer,ichip,maxnhit,bcid[ilayer][ichip][isca]);
	  if(bcid_seen>nslabshit ) {
	    event_monitoring->Fill(1,(14-ilayer)*20 +ichip);
	    event_monitoring_layer->Fill(1,(14-ilayer));
	  } else {
	    event_monitoring->Fill(2,(14-ilayer)*20 +ichip);
	    event_monitoring_layer->Fill(2,(14-ilayer));
	  }

	  if(bcid_seen>nslabshit && (badbcid[ilayer][ichip][isca]!=3 || retrigger==false) ) {
            event_monitoring_norettrains->Fill(1,(14-ilayer)*20 +ichip);
            event_monitoring_norettrains_layer->Fill(1,(14-ilayer));
	  }
	  if(bcid_seen<(nslabshit+1) && (badbcid[ilayer][ichip][isca]!=3 || retrigger==false) ) {
            event_monitoring_norettrains->Fill(2,(14-ilayer)*20 +ichip);
            event_monitoring_norettrains_layer->Fill(2,(14-ilayer));
          }
 
	  // int ntaggedasbad = 0;
	  // for(int ichn=0; ichn<64; ichn++) {
	  //   if(charge_lowGain[ilayer][ichip][isca][ichn]<180 && charge_lowGain[ilayer][ichip][isca][ichn]>-1) {
	  //     ntaggedasbad++;
	  //   }
	  //   if(charge_hiGain[ilayer][ichip][isca][ichn]<200 && charge_hiGain[ilayer][ichip][isca][ichn]>-1) {
	  //     ntaggedasbad++;
	  //   }
	  // }//ichn 

	  for(int ichn=0; ichn<64; ichn++) {
	    if(gain_hit_high[ilayer][ichip][isca][ichn]==1) {
	      if(bcid_seen>nslabshit ) { 
		hit_monitoring->Fill(1,(14-ilayer)*20 +ichip);
		hit_monitoring_layer->Fill(1,(14-ilayer));
	      } else {
		hit_monitoring->Fill(2,(14-ilayer)*20 +ichip);
		hit_monitoring_layer->Fill(2,(14-ilayer));
	      }
	      
	      if(bcid_seen>nslabshit && (badbcid[ilayer][ichip][isca]!=3 || retrigger==false)) {
		hit_monitoring_norettrains->Fill(1,(14-ilayer)*20 +ichip);
                hit_monitoring_norettrains_layer->Fill(1,(14-ilayer));
              }
	      if(bcid_seen<(nslabshit+1) && (badbcid[ilayer][ichip][isca]!=3 || retrigger==false) ) {
                hit_monitoring_norettrains->Fill(2,(14-ilayer)*20 +ichip);
                hit_monitoring_norettrains_layer->Fill(2,(14-ilayer));
              }
	    }
	  }

	  if(badbcid[ilayer][ichip][isca]==3) retrigger=true;

	}//isca
	if(last_sca>0) {
	  hit_monitoring->Fill(3+last_sca,(14-ilayer)*20 +ichip);
	  hit_monitoring_layer->Fill(3+last_sca,(14-ilayer));
	  hit_monitoring_norettrains->Fill(3+last_sca,(14-ilayer)*20 +ichip);
          hit_monitoring_norettrains_layer->Fill(3+last_sca,(14-ilayer));
	}
      }//ichip 
      
    }// ilayer 
  }  // end first loop analysis to fill pedestal historgrams
  
   
  // --------------------------------------------------------------------------------
  // PEDESTAL ANALYSIS
  TFile *pedfile = new TFile("results_proto/HitMonitoring_"+outputname+".root" , "RECREATE");
  pedfile->cd();
  hit_monitoring->Write();
  hit_monitoring_layer->Write();
  event_monitoring->Write();
  event_monitoring_layer->Write();
  hit_monitoring_norettrains->Write();
  hit_monitoring_norettrains_layer->Write();
  event_monitoring_norettrains->Write();
  event_monitoring_norettrains_layer->Write();
  
  

}



int DecodedSLBAnalysis::NSlabsAnalysis(TString outputname="", int maxnhit=1, int nslabshit=8)
{

  //Read the channel/chip -- x/y mapping
  //  ReadMap(map_filename);

  int nSLB=15;// get it as argument or read it from ntuple or somethng ? 
  if (fChain == 0) return -1;

  Long64_t nentries = fChain->GetEntriesFast();
  //  nentries/=10;
    

  // histograms for all scas, chn, chip and layers
  std::vector<std::vector<std::vector<std::vector<TH1F*> > > > ped_low_sca;
  std::vector<std::vector<std::vector<std::vector<TH1F*> > > > ped_high_sca;
  std::vector<std::vector<std::vector<std::vector<TH1F*> > > > mip_low_sca;
  std::vector<std::vector<std::vector<std::vector<TH1F*> > > > mip_high_sca;
  
  for(int ilayer=0; ilayer<nSLB; ilayer++) {
    
    std::vector<std::vector<std::vector<TH1F*> > > ped_lowtemp_sca_layer;
    std::vector<std::vector<std::vector<TH1F*> > > ped_hightemp_sca_layer;
    std::vector<std::vector<std::vector<TH1F*> > > mip_lowtemp_sca_layer;
    std::vector<std::vector<std::vector<TH1F*> > > mip_hightemp_sca_layer;

    for(int ichip=0; ichip<16; ichip++) {
      std::vector<std::vector<TH1F*> >ped_lowtemp_sca;
      std::vector<std::vector<TH1F*> >ped_hightemp_sca;
      std::vector<std::vector<TH1F*> >mip_lowtemp_sca;
      std::vector<std::vector<TH1F*> >mip_hightemp_sca;
      
      for(int ichn=0; ichn<64; ichn++) {
	std::vector<TH1F*> ped_lowtemp_sca2;
	std::vector<TH1F*> ped_hightemp_sca2;
	std::vector<TH1F*> mip_lowtemp_sca2;
	std::vector<TH1F*> mip_hightemp_sca2;
	
	for(int isca=0; isca<15; isca++) {
	  TH1F *ped_low_sca2 = new TH1F(TString::Format("ped_low_layer%i_chip%i_chn%i_sca%i",14-ilayer,ichip,ichn,isca),TString::Format("ped_low_layer%i_chip%i_chn%i_sca%i",14-ilayer,ichip,ichn,isca),400,100.5,500.5);
	  TH1F *ped_high_sca2 = new TH1F(TString::Format("ped_high_layer%i_chip%i_chn%i_sca%i",14-ilayer,ichip,ichn,isca),TString::Format("ped_high_layer%i_chip%i_chn%i_sca%i",14-ilayer,ichip,ichn,isca),400,100.5,500.5);
	  TH1F *mip_low_sca2 = new TH1F(TString::Format("mip_low_layer%i_chip%i_chn%i_sca%i",14-ilayer,ichip,ichn,isca),TString::Format("mip_low_layer%i_chip%i_chn%i_sca%i",14-ilayer,ichip,ichn,isca),500,100.5,600.5);
	  TH1F *mip_high_sca2 = new TH1F(TString::Format("mip_high_layer%i_chip%i_chn%i_sca%i",14-ilayer,ichip,ichn,isca),TString::Format("mip_high_layer%i_chip%i_chn%i_sca%i",14-ilayer,ichip,ichn,isca),500,100.5,600.5);
	  ped_lowtemp_sca2.push_back(ped_low_sca2);
	  ped_hightemp_sca2.push_back(ped_high_sca2);
	  mip_lowtemp_sca2.push_back(mip_low_sca2);
	  mip_hightemp_sca2.push_back(mip_high_sca2);
	}
	ped_lowtemp_sca.push_back(ped_lowtemp_sca2);
	ped_hightemp_sca.push_back(ped_hightemp_sca2);
	mip_lowtemp_sca.push_back(mip_lowtemp_sca2);
	mip_hightemp_sca.push_back(mip_hightemp_sca2);
      }
      ped_lowtemp_sca_layer.push_back(ped_lowtemp_sca);
      ped_hightemp_sca_layer.push_back(ped_hightemp_sca);
      mip_lowtemp_sca_layer.push_back(mip_lowtemp_sca);
      mip_hightemp_sca_layer.push_back(mip_hightemp_sca);
    }
    ped_low_sca.push_back(ped_lowtemp_sca_layer);
    ped_high_sca.push_back(ped_hightemp_sca_layer);
    mip_low_sca.push_back(mip_lowtemp_sca_layer);
    mip_high_sca.push_back(mip_hightemp_sca_layer);
  }
  //  TCanvas *tempcanvas = new TCanvas("temp","temp",400,400);
  
  Long64_t nbytes = 0, nb = 0;

  // -----------------------------------------------------------------------------------------------------   
  // Basic Pedestal/MIP analysis
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
      
    if(n_slboards>nSLB) nSLB=n_slboards;

      
    if ( jentry > 10 && jentry % 10 ==0 ) 
      printProgress(100.*jentry/nentries) ;

    for(int ilayer=0; ilayer<n_slboards; ilayer++) {
      for(int ichip=0; ichip<16; ichip++) {

	for(int isca=0; isca<15; isca++) {

	  if( badbcid[ilayer][ichip][isca]!=0 )  continue;
	  if( bcid[ilayer][ichip][isca]<50 ) continue;
	  if( bcid[ilayer][ichip][isca]>900 && bcid[ilayer][ichip][isca]<1000 )  continue;
	  if( nhits[ilayer][ichip][isca]>maxnhit ) continue;
	  int bcid_seen = SimpleCoincidenceTagger(ilayer,ichip,maxnhit,bcid[ilayer][ichip][isca]);
	  if(bcid_seen<(nslabshit-1)) continue;

	  for(int ichn=0; ichn<64; ichn++) {
	    if(gain_hit_low[ilayer][ichip][isca][ichn]==0) 
	      ped_low_sca.at(ilayer).at(ichip).at(ichn).at(isca)->Fill(charge_lowGain[ilayer][ichip][isca][ichn]);
	    if(gain_hit_high[ilayer][ichip][isca][ichn]==0) 
	      ped_high_sca.at(ilayer).at(ichip).at(ichn).at(isca)->Fill(charge_hiGain[ilayer][ichip][isca][ichn]);
	    
	    if(gain_hit_low[ilayer][ichip][isca][ichn]==1) 
	      mip_low_sca.at(ilayer).at(ichip).at(ichn).at(isca)->Fill(charge_lowGain[ilayer][ichip][isca][ichn]);
	    if(gain_hit_high[ilayer][ichip][isca][ichn]==1) 
	      mip_high_sca.at(ilayer).at(ichip).at(ichn).at(isca)->Fill(charge_hiGain[ilayer][ichip][isca][ichn]);
	  }
	}//isca
	  
      }//ichip 
	
    }// ilayer 
  }  // end first loop analysis to fill pedestal historgrams
  
   
  // --------------------------------------------------------------------------------
  // PEDESTAL ANALYSIS
  TFile *pedfile = new TFile("results_calib/PedestalMIP"+outputname+".root" , "RECREATE");
  pedfile->cd();
  TDirectory *cdhisto[nSLB];
  for(int ilayer=0; ilayer<nSLB; ilayer++) {
    cdhisto[ilayer] = pedfile->mkdir(TString::Format("layer_%i",14-ilayer));
  }
    

  // do pedestal (layer/chip/channel/sca based) analysis
  for(int ilayer=0; ilayer<nSLB; ilayer++) {
    cdhisto[ilayer]->cd();
      
    for(int ichip=0; ichip<16; ichip++) {
      for(int ichn=0; ichn<64; ichn++) {
	  
	for(int isca=0; isca<15; isca++) {
	  ped_low_sca.at(ilayer).at(ichip).at(ichn).at(isca)->SetTitle(TString::Format("ped_low_chip%i_chn%i_sca%i",ichip,ichn,isca));
	  ped_low_sca.at(ilayer).at(ichip).at(ichn).at(isca)->SetName(TString::Format("ped_low_chip%i_chn%i_sca%i",ichip,ichn,isca));
	  ped_low_sca.at(ilayer).at(ichip).at(ichn).at(isca)->Write();  
	  ped_high_sca.at(ilayer).at(ichip).at(ichn).at(isca)->SetTitle(TString::Format("ped_high_chip%i_chn%i_sca%i",ichip,ichn,isca));
	  ped_high_sca.at(ilayer).at(ichip).at(ichn).at(isca)->SetName(TString::Format("ped_high_chip%i_chn%i_sca%i",ichip,ichn,isca));
	  ped_high_sca.at(ilayer).at(ichip).at(ichn).at(isca)->Write();

	  mip_low_sca.at(ilayer).at(ichip).at(ichn).at(isca)->SetTitle(TString::Format("mip_low_chip%i_chn%i_sca%i",ichip,ichn,isca));
	  mip_low_sca.at(ilayer).at(ichip).at(ichn).at(isca)->SetName(TString::Format("mip_low_chip%i_chn%i_sca%i",ichip,ichn,isca));
	  mip_low_sca.at(ilayer).at(ichip).at(ichn).at(isca)->Write();  
	  mip_high_sca.at(ilayer).at(ichip).at(ichn).at(isca)->SetTitle(TString::Format("mip_high_chip%i_chn%i_sca%i",ichip,ichn,isca));
	  mip_high_sca.at(ilayer).at(ichip).at(ichn).at(isca)->SetName(TString::Format("mip_high_chip%i_chn%i_sca%i",ichip,ichn,isca));
	  mip_high_sca.at(ilayer).at(ichip).at(ichn).at(isca)->Write();  
	    
	}
      }
    }
  }

  ped_low_sca.clear();
  ped_high_sca.clear();
  mip_low_sca.clear();
  mip_high_sca.clear();

  pedfile->Close();
  delete pedfile;
    
  return 1;

}
