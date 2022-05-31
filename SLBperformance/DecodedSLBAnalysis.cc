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
  
  TH2F *goodhit_monitoring = new TH2F("goodhitmonitoring","BADBCID=0 Hit/cycle; chn; Layer*20+asic",64,-0.5,63.5,300,-0.5,299.5);//-1 = total, 0= good, 1=no coinc
  TH2F *badhit_monitoring = new TH2F("badhitmonitoring","BADBCID==3 Hit/cycle rate; chn; Layer*20+asic",64,-0.5,63.5,300,-0.5,299.5);//-1 = total, 0= good, 1=no coinc
  TH2F *sca_monitoring = new TH2F("sca_monitoring","Filled sca/cycle; sca; Layer*20+asic",15,-0.5,14.5,300,-0.5,299.5);
  TH2F *good_eventrate_monitoring = new TH2F("goodeventrate_monitoring","BADBCID==0 Evt/cycle; layer; asic",15,-0.5,14.5,16,-0.5,15.5);
  TH2F *bad_eventrate_monitoring = new TH2F("badeventrate_monitoring","BADBCID==3 Evt/cycle; layer; asic",15,-0.5,14.5,16,-0.5,15.5);

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

      for(int ichip=0; ichip<16; ichip++) {
      
       for(int isca=0; isca<15; isca++) {
         if( bcid[ilayer][ichip][isca]<10 ) continue;

	 if(badbcid[ilayer][ichip][isca]>-1 && badbcid[ilayer][ichip][isca]<3) good_eventrate_monitoring->Fill(ilayer,ichip);
	 if(badbcid[ilayer][ichip][isca]==3) bad_eventrate_monitoring->Fill(ilayer,ichip);

	 sca_monitoring->Fill(isca,(ilayer)*20. +ichip);

        for(int ichn=0; ichn<64; ichn++) {
	  if(hitbit_high[ilayer][ichip][isca][ichn]==1) {
	   if(badbcid[ilayer][ichip][isca]>-1 && badbcid[ilayer][ichip][isca]<3)
	     goodhit_monitoring->Fill(ichn,(ilayer)*20. +ichip);
	   if(badbcid[ilayer][ichip][isca]==3)
	     badhit_monitoring->Fill(ichn,(ilayer)*20. +ichip);
	  }
	}//ichn
	
       }//isca

       
      }//ichip 
      
    }// ilayer 
  }  // end first loop analysis to fill pedestal historgrams
  

  // --------------------------------------------------------------------------------
  // PEDESTAL ANALYSIS
  TFile *pedfile = new TFile("results_proto/stats_"+outputname+".root" , "RECREATE");
  pedfile->cd();
  goodhit_monitoring->Scale(1./nentries);
  badhit_monitoring->Scale(1./nentries);
  goodhit_monitoring->Write();
  badhit_monitoring->Write();

  sca_monitoring->Scale(1./nentries);
  sca_monitoring->Write();

  good_eventrate_monitoring->Scale(1./nentries);
  bad_eventrate_monitoring->Scale(1./nentries);
  good_eventrate_monitoring->Write();
  bad_eventrate_monitoring->Write();

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
         TH1F *ped_low_sca2 = new TH1F(TString::Format("ped_low_layer%i_chip%i_chn%i_sca%i",ilayer,ichip,ichn,isca),TString::Format("ped_low_layer%i_chip%i_chn%i_sca%i",ilayer,ichip,ichn,isca),400,100.5,500.5);
         TH1F *ped_high_sca2 = new TH1F(TString::Format("ped_high_layer%i_chip%i_chn%i_sca%i",ilayer,ichip,ichn,isca),TString::Format("ped_high_layer%i_chip%i_chn%i_sca%i",ilayer,ichip,ichn,isca),400,100.5,500.5);
         TH1F *mip_low_sca2 = new TH1F(TString::Format("mip_low_layer%i_chip%i_chn%i_sca%i",ilayer,ichip,ichn,isca),TString::Format("mip_low_layer%i_chip%i_chn%i_sca%i",ilayer,ichip,ichn,isca),900,100.5,1000.5);
         TH1F *mip_high_sca2 = new TH1F(TString::Format("mip_high_layer%i_chip%i_chn%i_sca%i",ilayer,ichip,ichn,isca),TString::Format("mip_high_layer%i_chip%i_chn%i_sca%i",ilayer,ichip,ichn,isca),900,100.5,1000.5);
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

	  if( badbcid[ilayer][ichip][isca]!=0 || bcid[ilayer][ichip][isca]<0)  continue;
	  //if( bcid[ilayer][ichip][isca]<10 ) continue;
	  //if( bcid[ilayer][ichip][isca]>950 && bcid[ilayer][ichip][isca]<1000 )  continue;
	  if( nhits[ilayer][ichip][isca]>maxnhit ) continue;

	  int bcid_seen = SimpleCoincidenceTagger(ilayer,maxnhit,bcid[ilayer][ichip][isca]);
	 	  
          if(bcid_seen<(nslabshit-1)) continue;

            /*int triggers=0;
            for(int ichn=0; ichn<64; ichn++) {
              if(hitbit_high[ilayer][ichip][isca][ichn]==1 ||hitbit_low[ilayer][ichip][isca][ichn]==1) {
                triggers++;
              }
          }*/
          // if(triggers==0 || triggers>1) continue;

	  for(int ichn=0; ichn<64; ichn++) {
	    if(hitbit_low[ilayer][ichip][isca][ichn]==0) 
	      ped_low_sca.at(ilayer).at(ichip).at(ichn).at(isca)->Fill(adc_low[ilayer][ichip][isca][ichn]);
	    if(hitbit_high[ilayer][ichip][isca][ichn]==0) 
	      ped_high_sca.at(ilayer).at(ichip).at(ichn).at(isca)->Fill(adc_high[ilayer][ichip][isca][ichn]);
	    
	    if(hitbit_low[ilayer][ichip][isca][ichn]==1) 
	      mip_low_sca.at(ilayer).at(ichip).at(ichn).at(isca)->Fill(adc_low[ilayer][ichip][isca][ichn]);
	    if(hitbit_high[ilayer][ichip][isca][ichn]==1) 
	      mip_high_sca.at(ilayer).at(ichip).at(ichn).at(isca)->Fill(adc_high[ilayer][ichip][isca][ichn]);
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
    cdhisto[ilayer] = pedfile->mkdir(TString::Format("layer_%i",ilayer));
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



int DecodedSLBAnalysis::NSlabsAnalysisNoise(TString outputname="", int gain=1, int iscamax=4, bool frominjection=true,int maxnhit=2, int nslabshit=8)
{

  //Read the channel/chip -- x/y mapping
  //  ReadMap(map_filename);

  int nSLB=15;// get it as argument or read it from ntuple or somethng ? 
  if (fChain == 0) return -1;

  Long64_t nentries = fChain->GetEntriesFast();
  //  nentries/=20;



  std::vector<std::vector<std::vector<TH2F *> > > cov;
  std::vector<std::vector<std::vector<TH2F *> > > nevents;
  std::vector<std::vector<std::vector<std::vector<TH1F *> > > > h_ped;
  std::vector<std::vector<std::vector<std::vector<double> > > > ped_value;

  for(int isca=0; isca<iscamax; isca++) {        
    std::vector<std::vector<TH2F* > > cov_;
    std::vector<std::vector<TH2F* > > nevents_;
    std::vector<std::vector<std::vector<TH1F* > > > h_ped_;
    std::vector<std::vector<std::vector<double > > > ped_;
    
    for(int ilayer=0; ilayer<15; ilayer++) {        
      std::vector<TH2F* > cov1;
      std::vector<TH2F* > nevents1;
      std::vector<std::vector<TH1F* > > h_ped1;
      std::vector<std::vector<double > > ped1;
      for(int i=0; i<16; i++) {                     
	TH2F * cov2 = new TH2F(TString::Format("cov_unnorm_layer%i_chip%i_sca%i",ilayer,i,isca),TString::Format("cov_unnorm_layer%i_chip%i_sca%i",ilayer,i,isca),64,-0.5,63.5,64,-0.5,63.5);
	TH2F * nevents2 = new TH2F(TString::Format("nevents_layer%i_chip%i_sca%i",ilayer,i,isca),TString::Format("nevents_layer%i_chip%i_sca%i",ilayer,i,isca),64,-0.5,63.5,64,-0.5,63.5);
	cov1.push_back(cov2);
	nevents1.push_back(nevents2);
	std::vector<TH1F* >  h_ped2;
	std::vector<double >  ped2;
	for(int j=0; j<64; j++) {
	  TH1F * h_ped3 = new TH1F(TString::Format("h_ped_layer%i_chip%i_chn%i_sca%i",ilayer,i,j,isca),TString::Format("h_ped_layer%i_chip%i_chn%i_sca%i",ilayer,i,j,isca),300,150.5,450.5);
	  h_ped2.push_back(h_ped3);
	  ped2.push_back(0);
	}
	h_ped1.push_back(h_ped2);
	ped1.push_back(ped2);
      }
      cov_.push_back(cov1);
      nevents_.push_back(nevents1);
      h_ped_.push_back(h_ped1);
      ped_.push_back(ped1);
    }
    cov.push_back(cov_);
    nevents.push_back(nevents_);
    h_ped.push_back(h_ped_);
    ped_value.push_back(ped_);
  }

 Long64_t nbytes = 0, nb = 0;


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
        for(int isca=0; isca<iscamax; isca++) {
	  if(bcid[ilayer][ichip][isca]<0) continue;
	  if(badbcid[ilayer][ichip][isca]!=0) continue;
	  if(nhits[ilayer][ichip][isca]==0) continue;

	  if(frominjection==false) {
	    if( nhits[ilayer][ichip][isca]>maxnhit ) continue;
	    int bcid_seen = SimpleCoincidenceTagger(ilayer,maxnhit,bcid[ilayer][ichip][isca]);
	    if(bcid_seen<(nslabshit-1)) continue;
	  }
     	  for(int j=0; j<64; j++) {
	    if(gain==1 && hitbit_high[ilayer][ichip][isca][j]==0 ) h_ped.at(isca).at(ilayer).at(ichip).at(j)->Fill(adc_high[ilayer][ichip][isca][j]);
	    if(gain==0 && hitbit_low[ilayer][ichip][isca][j]==0 ) h_ped.at(isca).at(ilayer).at(ichip).at(j)->Fill(adc_low[ilayer][ichip][isca][j]);
	  }
	  
	}
      }
    }
  }


 for(int isca=0; isca<iscamax; isca++) {
   for(int ilayer=0; ilayer<n_slboards; ilayer++) {
     for(int ichip=0; ichip<16; ichip++) {
       for(int ichn=0; ichn<64; ichn++) {
	 if(h_ped.at(isca).at(ilayer).at(ichip).at(ichn)->GetEntries()>50) {
	   h_ped.at(isca).at(ilayer).at(ichip).at(ichn)->GetXaxis()->SetRangeUser(h_ped.at(isca).at(ilayer).at(ichip).at(ichn)->GetMean()-15,h_ped.at(isca).at(ilayer).at(ichip).at(ichn)->GetMean()+15);
	   ped_value.at(isca).at(ilayer).at(ichip).at(ichn)=h_ped.at(isca).at(ilayer).at(ichip).at(ichn)->GetMean();
	 } else ped_value.at(isca).at(ilayer).at(ichip).at(ichn)=0;
	 
	 //	cout<<endl;
       }
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

     for(int isca=0; isca<iscamax; isca++) {
       //     if(sca_test<iscamax && isca!=sca_test) continue;
	  //	  if(isca>0) if(badbcid[ilayer][ichip][isca]-) continue;
       if(bcid[ilayer][ichip][isca]<0) continue;
       if(badbcid[ilayer][ichip][isca]!=0) continue;
       if(nhits[ilayer][ichip][isca]==0) continue;

       if(frominjection==false) {
	 if( nhits[ilayer][ichip][isca]>maxnhit ) continue;
	 int bcid_seen = SimpleCoincidenceTagger(ilayer,maxnhit,bcid[ilayer][ichip][isca]);
	 if(bcid_seen<(nslabshit-1)) continue;
       }

       for(int ichn=0; ichn<64; ichn++) {                   
         if(ped_value.at(isca).at(ilayer).at(ichip).at(ichn)<1 || hitbit_high[ilayer][ichip][isca][ichn]==1 ) continue;
         if(adc_high[ilayer][ichip][isca][ichn]<0) continue;

         double charge=adc_high[ilayer][ichip][isca][ichn]-ped_value.at(isca).at(ilayer).at(ichip).at(ichn);  
         if(gain==0) charge=adc_low[ilayer][ichip][isca][ichn]-ped_value.at(isca).at(ilayer).at(ichip).at(ichn);
         for(int kchn=0; kchn<64; kchn++) {
           if(ped_value.at(isca).at(ilayer).at(ichip).at(kchn)<1 || hitbit_high[ilayer][ichip][isca][kchn]==1) continue;
           if(adc_high[ilayer][ichip][isca][kchn]<0) continue;

           double charge_k=adc_high[ilayer][ichip][isca][kchn]-ped_value.at(isca).at(ilayer).at(ichip).at(kchn);   
           if(gain==0) charge_k=adc_low[ilayer][ichip][isca][kchn]-ped_value.at(isca).at(ilayer).at(ichip).at(kchn);    
           cov.at(isca).at(ilayer).at(ichip)->Fill(ichn,kchn,charge*charge_k);       
           nevents.at(isca).at(ilayer).at(ichip)->Fill(ichn,kchn);

	    }//kchn                                          
	  }//ichn  

	}//isca

      }//ichip 

    }// ilayer 
  }  // end first loop analysis to fill pedestal historgrams
  

  // --------------------------------------------------------------------------------
  // PEDESTAL ANALYSIS
  TString gain_st="highgain";
  if(gain==0) gain_st="lowgain";
  TString title=gain_st;

  for(int isca=0; isca<iscamax; isca++) {

    title=TString::Format("%s_sca%i",gain_st.Data(),isca);
    
    TFile *pedfile = new TFile("results_noise/NoiseCovariance_"+outputname+"_"+title+".root" , "RECREATE");
    pedfile->cd();
    TDirectory *cdhisto[nSLB];
    for(int ilayer=0; ilayer<nSLB; ilayer++) {
      cdhisto[ilayer] = pedfile->mkdir(TString::Format("layer_%i",ilayer));
    }
    
    // do pedestal (layer/chip/channel/sca based) analysis
    for(int ilayer=0; ilayer<nSLB; ilayer++) {
      cdhisto[ilayer]->cd();
      
      for(int ichip=0; ichip<16; ichip++) {
	if( cov.at(isca).at(ilayer).at(ichip)->GetEntries()>10) {
	  cov.at(isca).at(ilayer).at(ichip)->SetName(TString::Format("cov_unnorm_layer%i_chip%i",ilayer,ichip));
	  cov.at(isca).at(ilayer).at(ichip)->Write();
	  nevents.at(isca).at(ilayer).at(ichip)->SetName(TString::Format("nevents_layer%i_chip%i",ilayer,ichip));
	  nevents.at(isca).at(ilayer).at(ichip)->Write();
	}
	for(int ichn=0; ichn<64; ichn++) {
	  h_ped.at(isca).at(ilayer).at(ichip).at(ichn)->SetName(TString::Format("h_ped_layer%i_chip%i_chn%i",ilayer,ichip,ichn));
	  h_ped.at(isca).at(ilayer).at(ichip).at(ichn)->Write();
	}
      }
    }
  }
  // noise->Write();
  //hits->Write();
  //}


return 1;

}
