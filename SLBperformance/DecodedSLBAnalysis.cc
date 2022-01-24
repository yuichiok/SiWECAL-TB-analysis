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

	  if( bcid[ilayer][ichip][isca]<10 || (bcid[ilayer][ichip][isca]>950 && bcid[ilayer][ichip][isca]<1000) || badbcid[ilayer][ichip][isca]==1 || badbcid[ilayer][ichip][isca]==2 ) continue;

	  
	  int bcid_seen = SimpleCoincidenceTagger(ilayer,maxnhit,bcid[ilayer][ichip][isca]);
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
  TFile *pedfile = new TFile("results_proto/stats_"+outputname+".root" , "RECREATE");
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



int DecodedSLBAnalysis::NSlabsAnalysis_checks(TString outputname="")
{

  //Read the channel/chip -- x/y mapping
  //  ReadMap(map_filename);

  int nSLB=15;// get it as argument or read it from ntuple or somethng ? 
  if (fChain == 0) return -1;

  Long64_t nentries = fChain->GetEntriesFast();
  //  nentries/=10;

  
  /* TH2F *bcid_correl[16];
  for(int ichip=0; ichip<16; ichip++) {
    bcid_correl[ichip]= new TH2F(TString::Format("bcid_correl_fullASIC%i",ichip),TString::Format("bcid_correl_fullASIC%i; Layer*20+Chip  ; bcid-bcid_last_full_mem",ichip),300,-0.5,299.5,2001,-1000.5,1000.5);
    }*/
  TH1F * nslabs_hit[10];
  for(int i=0; i<10; i++) nslabs_hit[i]= new TH1F(TString::Format("nslabs_hit_%i",i),TString::Format("nslabs_hit_%i; nslabs; events",i),17,-0.5,16.5);

  /* TH2F * hits_layer_chip[20];
  for(int icoinc=0; icoinc<20; icoinc++) {
    hits_layer_chip[icoinc]= new TH2F(TString::Format("events_%i_Coincidences",icoinc+1),TString::Format("events_%i_Coincidences; layer ; chip",icoinc+1),15,-0.5,14.5,16,-0.5,15.5);
    }*/
   
  Long64_t nbytes = 0, nb = 0;

  // -----------------------------------------------------------------------------------------------------   
  // Basic Pedestal/MIP analysis
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
      
    if(n_slboards>nSLB) nSLB=n_slboards;

      
    if ( jentry > 10 && jentry % 500 ==0 ) 
      printProgress(100.*jentry/nentries) ;

    //analysis nslabs with hit
    //	std::vector<int> list_of_bcids;
    //    for(int ichip=0; ichip<16; ichip++) {
    for(int bcid_v=50; bcid_v<4096; bcid_v+=2) {
      if(bcid_v>900 && bcid_v<1000) continue;
      std::vector<int> bcid_seen = SimpleCoincidenceGlobal(1000,bcid_v);
      for(int i=0; i<bcid_seen.size(); i++) if(bcid_seen.at(i)>2 )nslabs_hit[i]->Fill(bcid_seen.at(i));
      //      if(bcid_seen.at(1)>2 )nslabs_hit_forcing_hits_first_4slabs->Fill(bcid_seen.at(1));
      // if(bcid_seen.at(2)>2 )nslabs_hit_forcing_hits_first_4slabs_loose->Fill(bcid_seen.at(2));
    }
    
      // }

    
    /*   // for(int ichip=0; ichip<16; ichip++) {
    for(int ichip=0; ichip<16; ichip++) {
      int last_filled_sca=-1;
      for(int isca=0; isca<15; isca++) {

	int ncoincidences_filled=0;
	for(int ilayer=0; ilayer<n_slboards; ilayer++) {
	  if(bcid[ilayer][ichip][isca]>-1) last_filled_sca=isca;
	  //cout<<ilayer<<" "<<ichip<<" "<<isca<<" "<<endl;
	  if( badbcid[ilayer][ichip][isca]!=0 )  continue;
	  if( bcid[ilayer][ichip][isca]<50 ) continue;
	  if( bcid[ilayer][ichip][isca]>900 && bcid[ilayer][ichip][isca]<1000 )  continue;
	  if(isca>0) if(bcid[ilayer][ichip][isca]-bcid[ilayer][ichip][isca-1]<2) continue;
	  if(isca<14) if(bcid[ilayer][ichip][isca+1]-bcid[ilayer][ichip][isca]<2) continue;
	  if(nhits[ilayer][ichip][isca]>1000 ) continue;
	  int bcid_seen = SimpleCoincidenceTaggerGlobal(ichip,1000,bcid[ilayer][ichip][isca]);
	  
	  if(bcid_seen>1 ) {
	    if(ncoincidences_filled==0 ) {
	      // nslabs_hit->Fill(bcid_seen);
	      //  cout<<jentry <<" FIRST FILLS "<<ilayer<<" "<<ichip<<" "<<isca<<" "<<bcid_seen<<" "<<bcid[ilayer][ichip][isca]<<endl;
	      hits_layer_chip[bcid_seen]->Fill(14-ilayer,ichip);
	    } else {
	      // cout<<jentry <<" EXTRA "<<ilayer<<" "<<ichip<<" "<<isca<<" "<<bcid_seen<<" "<<bcid[ilayer][ichip][isca]<<endl;
	      hits_layer_chip[bcid_seen]->Fill(14-ilayer,ichip);
	    }
	    ncoincidences_filled++;
	  }
	
	}
      }
    }
    */
    /*	if( last_filled_sca==14) {
	  for(int isca=0; isca<15; isca++) {
	    if( last_filled_sca==14) {
	      for(int jchip=0; jchip<16; jchip++) {
		if(jchip==ichip) continue;
		for(int isca=0; isca<15; isca++) {
		  bcid_correl[ichip]->Fill((14-ilayer)*20+jchip,bcid[ilayer][jchip][isca]-bcid[ilayer][ichip][14]);
		}
	      }
	    }
	  }
	}
      }//ichip
      }//ilayer*/
  }//iloop
   
  // --------------------------------------------------------------------------------
  // PEDESTAL ANALYSIS
  TFile *pedfile = new TFile("results_calib/TestCoincidences_"+outputname+".root" , "RECREATE");
  pedfile->cd();

  for(int i=0; i<4; i++)  nslabs_hit[i]->Write();
  //  nslabs_hit_forcing_hits_first_4slabs->Write();
  //nslabs_hit_forcing_hits_first_4slabs_loose->Write();
  /*for(int ichip=0; ichip<16; ichip++) {
    bcid_correl[ichip]->Write();
  }
  for(int icoinc=0; icoinc<15; icoinc++) {
    hits_layer_chip[icoinc]->Write();
    }*/
  
  pedfile->Close();
  delete pedfile;
    
  return 1;

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

	//int last_filled_sca=-1;
	//for(int isca=0; isca<15; isca++) if(bcid[ilayer][ichip][isca]>-1) last_filled_sca=isca;
	//        if(last_filled_sca==-1 || last_filled_sca==14) continue;

	for(int isca=0; isca<15; isca++) {

	  if( badbcid[ilayer][ichip][isca]!=0 )  continue;
	  if( bcid[ilayer][ichip][isca]<10 ) continue;
	  if( bcid[ilayer][ichip][isca]>950 && bcid[ilayer][ichip][isca]<1000 )  continue;
	  if( nhits[ilayer][ichip][isca]>maxnhit ) continue;
	  int bcid_seen = SimpleCoincidenceTagger(ilayer,maxnhit,bcid[ilayer][ichip][isca]);
	  if(bcid_seen<(nslabshit-1)) continue;

          /*int triggers=0;
          for(int ichn=0; ichn<64; ichn++) {
            if(gain_hit_high[ilayer][ichip][isca][ichn]==1 ||gain_hit_low[ilayer][ichip][isca][ichn]==1) {
              triggers++;
            }
	    }*/
	  // if(triggers==0 || triggers>1) continue;

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
  TFile *pedfile = new TFile("results_calib/PedestalMIP_"+outputname+".root" , "RECREATE");
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



int DecodedSLBAnalysis::NSlabsAnalysisNoise(TString outputname="", int maxnhit=1, int nslabshit=8, int gain=1, int sca_test=15)
{

  //Read the channel/chip -- x/y mapping
  //  ReadMap(map_filename);

  int nSLB=15;// get it as argument or read it from ntuple or somethng ? 
  if (fChain == 0) return -1;

  Long64_t nentries = fChain->GetEntriesFast();
  //  nentries/=10;
    


  std::vector<std::vector<TH2F *> > cov;
  std::vector<std::vector<TH2F *> > nevents;
  std::vector<std::vector<std::vector<TH1F *> > > h_ped;
  std::vector<std::vector<std::vector<double> > > ped_value;
  for(int ilayer=0; ilayer<15; ilayer++) {                                                                                                                                                                
    std::vector<TH2F* > cov1;
    std::vector<TH2F* > nevents1;
    std::vector<std::vector<TH1F* > > h_ped1;
    std::vector<std::vector<double > > ped1;
    for(int i=0; i<16; i++) {                                                                                                                                                                             
      TH2F * cov2 = new TH2F(TString::Format("cov_unnorm_layer%i_chip%i",ilayer,i),TString::Format("cov_unnorm_layer%i_chip%i",ilayer,i),64,-0.5,63.5,64,-0.5,63.5);
      TH2F * nevents2 = new TH2F(TString::Format("nevents_layer%i_chip%i",ilayer,i),TString::Format("nevents_layer%i_chip%i",ilayer,i),64,-0.5,63.5,64,-0.5,63.5);
      cov1.push_back(cov2);
      nevents1.push_back(nevents2);
      std::vector<TH1F* > h_ped2;
      std::vector<double > ped2;
      for(int j=0; j<64; j++) {
	TH1F * h_ped3 = new TH1F(TString::Format("h_ped_layer%i_chip%i_chn%i",ilayer,i,j),TString::Format("h_ped_layer%i_chip%i_chn%i",ilayer,i,j),300,150.5,450.5);
	h_ped2.push_back(h_ped3);
	ped2.push_back(0);
      }
      h_ped1.push_back(h_ped2);
      ped1.push_back(ped2);
    }
    cov.push_back(cov1);
    nevents.push_back(nevents1);
    h_ped.push_back(h_ped1);
    ped_value.push_back(ped1);
  }

  Long64_t nbytes = 0, nb = 0;

  //  if(gain==1) ReadPedestalsProto("../pedestals/Pedestal_3GeVMIPscan_highgain.txt");
  //else ReadPedestalsProto("../pedestals/Pedestal_3GeVMIPscan_lowgain.txt");

  ReadMasked("../masked/masked_channels_"+outputname+".txt");
  
  // -----------------------------------------------------------------------------------------------------
  // Pedestal Analysis
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(n_slboards>nSLB) nSLB=n_slboards;


    if ( jentry > 10 && jentry % 10 ==0 )
      printProgress(100.*jentry/nentries) ;

    //    for(int ichip=0; ichip<16; ichip++) {
    //std::vector<int> list_of_bcids;
    for(int ilayer=0; ilayer<n_slboards; ilayer++) {
      for(int ichip=0; ichip<16; ichip++) {
	//    int last_filled_sca=-1;
	//    for(int isca=0; isca<15; isca++) if(bcid[ilayer][ichip][isca]>-1) last_filled_sca=isca;
	//    if(last_filled_sca==-1 || last_filled_sca==14) continue;

        for(int isca=0; isca<15; isca++) {
          if(sca_test<15 && isca!=sca_test) continue;
          if( badbcid[ilayer][ichip][isca]!=0 )  continue;
          if( bcid[ilayer][ichip][isca]<10 ) continue;
          if( bcid[ilayer][ichip][isca]>950 && bcid[ilayer][ichip][isca]<1000 )  continue;
          if( nhits[ilayer][ichip][isca]>maxnhit ) continue;
          int bcid_seen = SimpleCoincidenceTagger(ilayer,maxnhit,bcid[ilayer][ichip][isca]);
	  //	  cout<<jentry<<" "<<ilayer<<" "<<ichip<<" "<<isca<<" "<<bcid[ilayer][ichip][isca]<<" "<<bcid_seen<<endl;
	  if(bcid_seen<(nslabshit-1)) continue;

	  /*if(bcid_seen>nslabshit ) {
	    if(list_of_bcids.size()==0) {
	      list_of_bcids.push_back(bcid[ilayer][ichip][isca]);
	      for(int j=0; j<64; j++) {
		if(gain==1) h_ped.at(ilayer).at(ichip).at(j)->Fill(charge_hiGain[ilayer][ichip][isca][j]);
		if(gain==0) h_ped.at(ilayer).at(ichip).at(j)->Fill(charge_lowGain[ilayer][ichip][isca][j]);
	      }
	    } else {
	      int n=list_of_bcids.size();
	      for(int ibcid=0; ibcid<n; ibcid++) {
		if(fabs(bcid[ilayer][ichip][isca]-list_of_bcids.at(ibcid))>3) {
		  //  cout<<" BREAK "<<endl;
		  list_of_bcids.push_back(bcid[ilayer][ichip][isca]);
		  for(int j=0; j<64; j++) {
		    if(gain==1) h_ped.at(ilayer).at(ichip).at(j)->Fill(charge_hiGain[ilayer][ichip][isca][j]);
		    if(gain==0) h_ped.at(ilayer).at(ichip).at(j)->Fill(charge_lowGain[ilayer][ichip][isca][j]);
		  }
		  break;
		}
	      }
	    }
	    }*/
	  
	  
	  //  int triggers=0;
	  //  for(int ichn=0; ichn<64; ichn++) {
	  //    if(gain_hit_high[ilayer][ichip][isca][ichn]==1 ||gain_hit_low[ilayer][ichip][isca][ichn]==1) {
	  //      triggers++;
	  //   }
	  //  }
	  //  if(triggers==0 || triggers>1) continue;
	  for(int j=0; j<64; j++) {
	    if(gain==1 && gain_hit_high[ilayer][ichip][isca][j]==0) h_ped.at(ilayer).at(ichip).at(j)->Fill(charge_hiGain[ilayer][ichip][isca][j]);
	    if(gain==0 && gain_hit_low[ilayer][ichip][isca][j]==0 ) h_ped.at(ilayer).at(ichip).at(j)->Fill(charge_lowGain[ilayer][ichip][isca][j]);
	  }

	}
      }
    }
  }

  for(int ilayer=0; ilayer<n_slboards; ilayer++) {
    for(int ichip=0; ichip<16; ichip++) {
      for(int ichn=0; ichn<64; ichn++) {
	if(h_ped.at(ilayer).at(ichip).at(ichn)->GetEntries()>10) {
	  h_ped.at(ilayer).at(ichip).at(ichn)->GetXaxis()->SetRangeUser(h_ped.at(ilayer).at(ichip).at(ichn)->GetMean()-20,h_ped.at(ilayer).at(ichip).at(ichn)->GetMean()+20);
	  ped_value.at(ilayer).at(ichip).at(ichn)=h_ped.at(ilayer).at(ichip).at(ichn)->GetMean();
	} else ped_value.at(ilayer).at(ichip).at(ichn)=0;
      }
    }
  }
  
  // -----------------------------------------------------------------------------------------------------   
  // COVARIANCE ANALYSIS
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
      
    if(n_slboards>nSLB) nSLB=n_slboards;

      
    if ( jentry > 10 && jentry % 10 ==0 ) 
      printProgress(100.*jentry/nentries) ;

    //    for(int ichip=0; ichip<16; ichip++) {
    //      std::vector<int> list_of_bcids;
    for(int ilayer=0; ilayer<n_slboards; ilayer++) {
      for(int ichip=0; ichip<16; ichip++) {
	//	int last_filled_sca=-1;
	//for(int isca=0; isca<15; isca++) if(bcid[ilayer][ichip][isca]>-1) last_filled_sca=isca;
	//if(last_filled_sca==-1 || last_filled_sca==14) continue;

	for(int isca=0; isca<15; isca++) {

	  if(sca_test<15 && isca!=sca_test) continue;
	  if( badbcid[ilayer][ichip][isca]!=0 )  continue;
	  if( bcid[ilayer][ichip][isca]<10 ) continue;
	  if( bcid[ilayer][ichip][isca]>950 && bcid[ilayer][ichip][isca]<1000 )  continue;
	  if( nhits[ilayer][ichip][isca]>maxnhit ) continue;
	  int bcid_seen = SimpleCoincidenceTagger(ilayer,maxnhit,bcid[ilayer][ichip][isca]);
	  if(bcid_seen<(nslabshit-1)) continue;

	  for(int ichn=0; ichn<64; ichn++) {                                                                                                                                                       
	    if(ped_value.at(ilayer).at(ichip).at(ichn)<100) continue;                                                                                                                              
	    if(gain_hit_high[ilayer][ichip][isca][ichn]==1 || gain_hit_low[ilayer][ichip][isca][ichn]==1) {                                                                                        
	      continue;                                                                                                                                                                            
	    }                                                                                                                                                                                      
	    double charge=charge_hiGain[ilayer][ichip][isca][ichn]-ped_value.at(ilayer).at(ichip).at(ichn);                                                                                        
	    if(gain==0) charge=charge_lowGain[ilayer][ichip][isca][ichn]-ped_value.at(ilayer).at(ichip).at(ichn);                                                                                  
	    for(int kchn=0; kchn<64; kchn++) {
	      if(gain_hit_high[ilayer][ichip][isca][kchn]==1) continue;                                                                                                                            
	      if(ped_value.at(ilayer).at(ichip).at(kchn)<100) continue;                                                                                                                            
	      double charge_k=charge_hiGain[ilayer][ichip][isca][kchn]-ped_value.at(ilayer).at(ichip).at(kchn);                                                                                    
	      if(gain==0) charge_k=charge_lowGain[ilayer][ichip][isca][kchn]-ped_value.at(ilayer).at(ichip).at(kchn);                                                                              
	      cov.at(14-ilayer).at(ichip)->Fill(ichn,kchn,charge*charge_k);                                                                                                                        
	      nevents.at(14-ilayer).at(ichip)->Fill(ichn,kchn);
	    }//kchn                                                                                                                                                                                
	  }//ichn  

	  /*  int triggers=0;
	  for(int ichn=0; ichn<64; ichn++) {
	    if(gain_hit_high[ilayer][ichip][isca][ichn]==1 ||gain_hit_low[ilayer][ichip][isca][ichn]==1) {
	      triggers++;
	    }
	  }
	  if(triggers==0 || triggers>1) continue;*/
	  /*if(bcid_seen>nslabshit ) {
	    if(list_of_bcids.size()==0) {
	      list_of_bcids.push_back(bcid[ilayer][ichip][isca]);
	      for(int ichn=0; ichn<64; ichn++) {
		if(ped_value.at(ilayer).at(ichip).at(ichn)<100) continue;
		if(gain_hit_high[ilayer][ichip][isca][ichn]==1 || gain_hit_low[ilayer][ichip][isca][ichn]==1) {
		  continue;
		}
		double charge=charge_hiGain[ilayer][ichip][isca][ichn]-ped_value.at(ilayer).at(ichip).at(ichn);
		if(gain==0) charge=charge_lowGain[ilayer][ichip][isca][ichn]-ped_value.at(ilayer).at(ichip).at(ichn);
		for(int kchn=0; kchn<64; kchn++) {
		  if(gain_hit_high[ilayer][ichip][isca][kchn]==1) continue;
		  if(ped_value.at(ilayer).at(ichip).at(kchn)<100) continue;
		  double charge_k=charge_hiGain[ilayer][ichip][isca][kchn]-ped_value.at(ilayer).at(ichip).at(kchn);
		  if(gain==0) charge_k=charge_lowGain[ilayer][ichip][isca][kchn]-ped_value.at(ilayer).at(ichip).at(kchn);
		  cov.at(14-ilayer).at(ichip)->Fill(ichn,kchn,charge*charge_k);
		  nevents.at(14-ilayer).at(ichip)->Fill(ichn,kchn);
		}//kchn
	      }//ichn
	    } else {
	      int n=list_of_bcids.size();
	      for(int ibcid=0; ibcid<n; ibcid++) {
		if(fabs(bcid[ilayer][ichip][isca]-list_of_bcids.at(ibcid))>3) {
		  //  cout<<" BREAK "<<endl;
		  list_of_bcids.push_back(bcid[ilayer][ichip][isca]);
		  for(int ichn=0; ichn<64; ichn++) {
		    if(ped_value.at(ilayer).at(ichip).at(ichn)<100) continue;
		    if(gain_hit_high[ilayer][ichip][isca][ichn]==1 || gain_hit_low[ilayer][ichip][isca][ichn]==1) {
		      continue;
		    }
		    double charge=charge_hiGain[ilayer][ichip][isca][ichn]-ped_value.at(ilayer).at(ichip).at(ichn);
		    if(gain==0) charge=charge_lowGain[ilayer][ichip][isca][ichn]-ped_value.at(ilayer).at(ichip).at(ichn);
		    for(int kchn=0; kchn<64; kchn++) {
		      if(gain_hit_high[ilayer][ichip][isca][kchn]==1) continue;
		      if(ped_value.at(ilayer).at(ichip).at(kchn)<100) continue;
		      double charge_k=charge_hiGain[ilayer][ichip][isca][kchn]-ped_value.at(ilayer).at(ichip).at(kchn);
		      if(gain==0) charge_k=charge_lowGain[ilayer][ichip][isca][kchn]-ped_value.at(ilayer).at(ichip).at(kchn);
		      cov.at(14-ilayer).at(ichip)->Fill(ichn,kchn,charge*charge_k);
		      nevents.at(14-ilayer).at(ichip)->Fill(ichn,kchn);
		    }//kchn
		  }//ichn
		  break;
		}
	      }
	    }
	    }*/
	  

	}//isca
 
      }//ichip 
	
    }// ilayer 
  }  // end first loop analysis to fill pedestal historgrams
  
   
  // --------------------------------------------------------------------------------
  // PEDESTAL ANALYSIS
  TString gain_st="highgain";
  if(gain==0) gain_st="lowgain";
  TString title=gain_st;
  if(sca_test<15)  title=TString::Format("%s_sca%i",gain_st.Data(),sca_test);
  
  TFile *pedfile = new TFile("noise_covariance_matrix/NoiseCovariance_"+outputname+"_"+title+".root" , "RECREATE");
  pedfile->cd();
  TDirectory *cdhisto[nSLB];
  for(int ilayer=0; ilayer<nSLB; ilayer++) {
    cdhisto[ilayer] = pedfile->mkdir(TString::Format("layer_%i",ilayer));
  }
 
  // do pedestal (layer/chip/channel/sca based) analysis
  for(int ilayer=0; ilayer<nSLB; ilayer++) {
    cdhisto[ilayer]->cd();

    for(int ichip=0; ichip<16; ichip++) {
      if( cov.at(ilayer).at(ichip)->GetEntries()>400) {
	cov.at(ilayer).at(ichip)->Write();
	nevents.at(ilayer).at(ichip)->Write();
      }
      for(int ichn=0; ichn<64; ichn++) h_ped.at(14-ilayer).at(ichip).at(ichn)->Write();
    }
  }
  // noise->Write();
  //hits->Write();
  //}
	  
   
  return 1;

}
