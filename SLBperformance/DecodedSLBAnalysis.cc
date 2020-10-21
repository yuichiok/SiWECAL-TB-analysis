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

std::array<int,4096> DecodedSLBAnalysis::SimpleCoincidenceTagger(int jentry, int maxnhit)
{


  // so far we only look for coincicendes for exact same values of bcids... in principle, we could lookf for bcid+-1
  std::array<int,4096> bcid_seen={0};

  for(int islboard=0; islboard<n_slboards; islboard++) {
    int bcid_seen_slb[4096]={0};

    for(int ichip=0; ichip<16; ichip++) {
      for(int isca=0; isca<15; isca++) {
	if(bcid[islboard][ichip][isca]<0) break; //not include elements with no value (-999)

	if( badbcid[islboard][ichip][isca]==0 && nhits[islboard][ichip][isca]<(maxnhit+1)) {
	  bool gooddata=true;
	  for(int i=0; i<64; i++) {
	    if(charge_hiGain[islboard][ichip][isca][i]<150) {
	      gooddata=false;
	      break;
	    }
	  }
	  if(gooddata==true)  bcid_seen_slb[bcid[islboard][ichip][isca]]++; //save the recorded bcids, we add a count for each time that the bcid has been recorded
	}
      }//end isca
    }
    
    for(int i=0; i<4096; i++) {
      if(bcid_seen_slb[i]==1 ) {//only count bcids with signal in one chip, not in several
	bcid_seen.at(i)++;
      }
    }
    //requiring ==1 we filter events in which the same bcid has triggered in several chips at the same time in the same board
    // unfortunately we also remvoe events with that eventually happen for the same bcid after overcyclying...
  }

  return bcid_seen;
  
}

void DecodedSLBAnalysis::SlowControlMonitoring(TString outputname="SlowControlMonitoring")
{

  TGraph* temp[30];
  TGraph* avdd0[30];
  TGraph* avdd_diff[30];
  TGraph* avdd_temp[30];

  // --------------
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();

  int ngraph=nentries/1000;
  if(ngraph<1) ngraph=1;
  //-----------------
  cout<<"Total number of entries: "<< nentries<<endl;
  cout<<"Fill graph every "<< ngraph <<" cycles"<<endl;

  double x[30][1000]={0};
  double temp_y[30][1000]={0};
  double avdd_y[30][1000]={0};
  double avdd_diff_y[30][1000]={0};
  double avdd_temp_y[30][1000]={0};
  int n[30]={0};
  // -----------------------------------------------------------------------------------------------------
  // Signal readout
  Long64_t nbytes = 0, nb = 0;
  int n_SLB=0;
  for (Long64_t jentry=0; jentry<nentries;jentry+=ngraph) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if ( jentry > 1000 && jentry % 1000 ==0 ) std::cout << "Progress: " << 100.*jentry/nentries <<" %"<<endl;

    if(jentry==0) n_SLB=n_slboards;

    if ( jentry % (ngraph) !=0 ) continue;
   
    for(int islboard=0; islboard<n_slboards; islboard++) {
      if(TSD[islboard]>0) {
	x[slboard_id[islboard]][n[slboard_id[islboard]]]=jentry;
	temp_y[slboard_id[islboard]][n[slboard_id[islboard]]]=TSD[islboard];
	avdd_y[slboard_id[islboard]][n[slboard_id[islboard]]]=AVDD0[islboard];
	avdd_diff_y[slboard_id[islboard]][n[slboard_id[islboard]]]=(AVDD1[islboard]-AVDD0[islboard]);
	avdd_temp_y[slboard_id[islboard]][n[slboard_id[islboard]]]=AVDD0[islboard]/TSD[islboard];

	n[slboard_id[islboard]]++;
      }
    }
    

  }

  for(int i=0; i<30; i++) {
    for(int j=0; j<n[i]; j++) {
      temp_y[i][j]= ( temp_y[i][j] -  temp_y[i][10]);
    }
  }

  // -----------------------------------------------------------------------------------------------------   
  for(int i=0; i<30; i++) {
    if(n[i]!=0) {
      temp[i] = new TGraph(n[i],x[i],temp_y[i]);
      avdd0[i] = new TGraph(n[i],x[i],avdd_y[i]);
      avdd_diff[i] = new TGraph(n[i],x[i],avdd_diff_y[i]);
      avdd_temp[i] = new TGraph(n[i],x[i],avdd_temp_y[i]);
    }
  }

  TFile *monitoringfile_summary = new TFile("results_monitoring/SlowControl_"+outputname+".root" , "RECREATE");
  monitoringfile_summary->cd();
  
  for(int i=0; i<11; i++) {
    if(n[i]!=0) {
      temp[i]->SetName(TString::Format("temp_slboard_%i",i));
      temp[i]->Write();
      avdd0[i]->SetName(TString::Format("AVDD0_slboard_%i",i));
      avdd0[i]->Write();
      avdd_diff[i]->SetName(TString::Format("AVDDdiff_slboard_%i",i));
      avdd_diff[i]->Write();
      avdd_temp[i]->SetName(TString::Format("AVDDtemp_slboard_%i",i));
      avdd_temp[i]->Write();
    }
  }
  
  TCanvas * canvas = new TCanvas("canvas","canvas",1600,1600);
  canvas->Divide(2,2);
  canvas->cd(1);
  TLegend * leg = new TLegend(0.7,0.2,0.9,0.7);
  bool first=true;
  for(int i=0; i<30; i++) {
   if(n[i]!=0 && first==false) {
      temp[i]->SetLineColor(i+1);
      temp[i]->SetLineWidth(2);
      temp[i]->SetMarkerColor(i+1);
      temp[i]->SetMarkerStyle(20+i);
      temp[i]->SetLineStyle( i % 2 + 1);
      //   temp[i]->Draw("lp");
      leg->AddEntry(temp[i],TString::Format("SL@ %i",i),"lp");
    }
   if(n[i]!=0 && first==true) {
      temp[i]->GetXaxis()->SetTitle("# Acq. Since Start");
      temp[i]->GetYaxis()->SetTitle("Temp Variation since start [Celsius] ");
      temp[i]->SetLineColor(i+1);
      temp[i]->SetLineWidth(2);
      temp[i]->SetMarkerColor(i+1);
      temp[i]->SetMarkerStyle(20+i);
      temp[i]->SetLineStyle( i % 2 + 1);
      //  temp[i]->Draw("alp");
      first=false;
      leg->AddEntry(temp[i],TString::Format("SL@ %i",i),"lp");
    }
  }
  leg->Draw();
  
  canvas->cd(2);
  first=true;
  for(int i=0; i<30; i++) {
    if(n[i]!=0 && first==false) {
      avdd0[i]->SetLineColor(i+1);
      avdd0[i]->SetLineWidth(2);
      avdd0[i]->SetMarkerColor(i+1);
      avdd0[i]->SetMarkerStyle(20+i);
      avdd0[i]->SetLineStyle( i % 2 + 1);
      avdd0[i]->Draw("lp");
    }
    if(n[i]!=0 && first==true) {
      avdd0[i]->GetXaxis()->SetTitle("# Acq. Since Start");
      avdd0[i]->GetYaxis()->SetTitle("Avdd start ACQ [mV] ");
      avdd0[i]->GetYaxis()->SetRangeUser(3000,3700);
      avdd0[i]->SetLineColor(i+1);
      avdd0[i]->SetLineWidth(2);
      avdd0[i]->SetMarkerColor(i+1);
      avdd0[i]->SetMarkerStyle(20+i);
      avdd0[i]->SetLineStyle( i % 2 + 1);
      avdd0[i]->Draw("alp");
      first=false;
    }
  }

  canvas->cd(3);
   first=true;
  for(int i=0; i<30; i++) {
    if(n[i]!=0 && first==false) {
      avdd_diff[i]->SetLineColor(i+1);
      avdd_diff[i]->SetLineWidth(2);
      avdd_diff[i]->SetMarkerColor(i+1);
      avdd_diff[i]->SetMarkerStyle(20+i);
      avdd_diff[i]->SetLineStyle( i % 2 + 1);
      avdd_diff[i]->Draw("lp");
    }
    if(n[i]!=0 && first==true) {
      avdd_diff[i]->GetXaxis()->SetTitle("# Acq. Since Start");
      avdd_diff[i]->GetYaxis()->SetTitle("Avdd Start ACQ - Stop Acq [mV] ");
      avdd_diff[i]->GetYaxis()->SetRangeUser(-20,10);
      avdd_diff[i]->SetLineColor(i+1);
      avdd_diff[i]->SetLineWidth(2);
      avdd_diff[i]->SetMarkerColor(i+1);
      avdd_diff[i]->SetMarkerStyle(20+i);
      avdd_diff[i]->SetLineStyle( i % 2 + 1);
      avdd_diff[i]->Draw("alp");
      first=false;
    }
  }
  canvas->cd(4);
   first=true;
   for(int i=0; i<30; i++) {
   if(n[i]!=0 && first==false) {
      temp[i]->SetLineColor(i+1);
      temp[i]->SetLineWidth(2);
      temp[i]->SetMarkerColor(i+1);
      temp[i]->SetMarkerStyle(20+i);
      temp[i]->SetLineStyle( i % 2 + 1);
      temp[i]->Draw("lp");
      //leg->AddEntry(temp[i],TString::Format("SL@ %i",i),"lp");
    }
   if(n[i]!=0 && first==true) {
      temp[i]->GetXaxis()->SetTitle("# Acq. Since Start");
      temp[i]->GetYaxis()->SetTitle("Temp Variation since start [Celsius] ");
      temp[i]->SetLineColor(i+1);
      temp[i]->SetLineWidth(2);
      temp[i]->SetMarkerColor(i+1);
      temp[i]->SetMarkerStyle(20+i);
      temp[i]->SetLineStyle( i % 2 + 1);
      temp[i]->Draw("alp");
      first=false;
      // leg->AddEntry(temp[i],TString::Format("SL@ %i",i),"lp");
    }
  }
  
  canvas->Write();
  
  monitoringfile_summary->Close();

}


void DecodedSLBAnalysis::QuickDisplay(TString outputname="QuickDisplay", int nlayers_minimum=8)
{

  TH3F* event_display[1000];

  for(int j=0; j<1000; j++) {
    event_display[j]=new TH3F(TString::Format("event_display_coinc_xyz_%i",j),TString::Format("event_display_coinc_xyz_%i",j),32,-90,90,32,-90,90,15,-232.5,7.5);
  }
  
  // --------------
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //-----------------
  cout<<"Total number of entries: "<< nentries<<endl;


  // -----------------------------------------------------------------------------------------------------
  // Signal readout
  Long64_t nbytes = 0, nb = 0;
  int n_SLB=0;
  int n=0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if ( jentry > 1000 && jentry % 1000 ==0 ) std::cout << "Progress: " << 100.*jentry/nentries <<" %"<<endl;

    if(jentry==0) n_SLB=n_slboards;

    if(n>999) break;
    std::array<int,4096> bcid_seen = SimpleCoincidenceTagger(jentry,5);

    int nhistos[4096]={0};
    int counter_events=0;
    for(int i=0; i<4096; i++) {
      if(bcid_seen.at(i)>(nlayers_minimum-1) ){
	counter_events++;
	nhistos[i]=counter_events;
      }
    }

    if( n+counter_events > 999) break;

    for(int islboard=0; islboard<n_slboards; islboard++) {
      for(int ichip=0; ichip<16; ichip++) {
	for(int isca=0; isca<15; isca++) {
	  if(bcid[islboard][ichip][isca]<0) continue;
	  
	  if(badbcid[islboard][ichip][isca]==0 && bcid_seen.at(bcid[islboard][ichip][isca])>(nlayers_minimum-1)) {
	    for(int ichn=0;ichn<64; ichn++) {
	      double z=-15.*double(slot[islboard]);
	      if(gain_hit_high[islboard][ichip][isca][ichn]==1)
		event_display[n+nhistos[bcid[islboard][ichip][isca]]-1]->Fill(double(map_pointX[islboard][ichip][ichn]),double(map_pointY[islboard][ichip][ichn]),z,double(charge_hiGain[islboard][ichip][isca][ichn]));
	    }
	  }
	}
      }
    }
  
    n+=counter_events;
    
  }

  // -----------------------------------------------------------------------------------------------------   
  // Signal analysis
  TFile *monitoringfile_summary = new TFile("results_monitoring/Cosmics_"+outputname+".root" , "RECREATE");
  monitoringfile_summary->cd();
  
  for(int i=0; i<1000; i++) {
    event_display[i]->GetXaxis()->SetTitle("X [mm]");
    event_display[i]->GetYaxis()->SetTitle("Y [mm]");
    event_display[i]->GetZaxis()->SetTitle("Z (0=upper module) [mm]");
    if(event_display[i]->GetEntries()>0) event_display[i]->Write();
  }
  
  
  monitoringfile_summary->Close();

}


int DecodedSLBAnalysis::NSlabsAnalysis(TString outputname="", int maxnhit=5)
{

  //Read the channel/chip -- x/y mapping
  //  ReadMap(map_filename);

  int nSLB=1;// get it as argument or read it from ntuple or somethng ? 
  if (fChain == 0) return -1;

  Long64_t nentries = fChain->GetEntriesFast();
  //  nentries/=10;
    
 
  // histograms for all scas, chn, chip and slboards
  std::vector<std::vector<std::vector<std::vector<TH1F*> > > > ped_sca;
  std::vector<std::vector<std::vector<std::vector<TH1F*> > > > ped_sca_tagged;

  for(int islboard=0; islboard<nSLB; islboard++) {
    
    std::vector<std::vector<std::vector<TH1F*> > > pedtemp_sca_slboard;
    std::vector<std::vector<std::vector<TH1F*> > > pedtemp_sca_tagged_slboard;
    
    for(int ichip=0; ichip<16; ichip++) {
      std::vector<std::vector<TH1F*> >pedtemp_sca;
      std::vector<std::vector<TH1F*> >pedtemp_sca_tagged;
      
      for(int ichn=0; ichn<64; ichn++) {
	std::vector<TH1F*> pedtemp_sca2;
	std::vector<TH1F*> pedtemp_sca2_tagged;
	
	for(int isca=0; isca<15; isca++) {
	  TH1F *ped_sca2 = new TH1F(TString::Format("ped_slboard%i_chip%i_chn%i_sca%i",islboard,ichip,ichn,isca),TString::Format("ped_slboard%i_chip%i_chn%i_sca%i",islboard,ichip,ichn,isca),1000,0.5,1000.5);
	  TH1F *ped_sca2_tagged = new TH1F(TString::Format("ped_tagged_slboard%i_chip%i_chn%i_sca%i",islboard,ichip,ichn,isca),TString::Format("ped_tagged_slboard%i_chip%i_chn%i_sca%i",islboard,ichip,ichn,isca),1000,0.5,1000.5);
	  pedtemp_sca2.push_back(ped_sca2);
	  pedtemp_sca2_tagged.push_back(ped_sca2_tagged);
	}
	pedtemp_sca.push_back(pedtemp_sca2);
	pedtemp_sca_tagged.push_back(pedtemp_sca2_tagged);
      }
      pedtemp_sca_slboard.push_back(pedtemp_sca);
      pedtemp_sca_tagged_slboard.push_back(pedtemp_sca_tagged);
    }
    ped_sca.push_back(pedtemp_sca_slboard);
    ped_sca_tagged.push_back(pedtemp_sca_tagged_slboard); 
  }
  TCanvas *tempcanvas = new TCanvas("temp","temp",400,400);
  // -----------------------------------------------------------------------------------------------------   
  // SCA analysis
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // --------------------
    //  if(jentry==0) {
    if(n_slboards>nSLB) nSLB=n_slboards;
      //  }
 

    if ( jentry > 1000 && jentry % 1000 ==0 ) std::cout << "Progress: " << 100.*jentry/nentries <<" %"<<endl;

    for(int islboard=0; islboard<n_slboards; islboard++) {
      for(int ichip=0; ichip<16; ichip++) {

	for(int isca=0; isca<15; isca++) {

	  bool gooddata=true;
	  int ntaggedasbad = 0;
	  for(int ichn=0; ichn<64; ichn++) {
	    if(gain_hit_high[islboard][ichip][isca][ichn]==1 && charge_hiGain[islboard][ichip][isca][ichn]<15 && charge_hiGain[islboard][ichip][isca][ichn]>-1 ) 
	      ntaggedasbad++;
	  }//ichn 
	  if ( ntaggedasbad > 0) gooddata=false;
		
	  for(int ichn=0; ichn<64; ichn++) {

	    //good events
	    bool selection=false;

	    if( charge_hiGain[islboard][ichip][isca][ichn]>30 && badbcid[islboard][ichip][isca]==0 && nhits[islboard][ichip][isca]<maxnhit+1 && bcid[islboard][ichip][isca]>50) selection=true;
		  
	    // if(masked[islboard][ichip][ichn]==1) selection=false;
	    if(gain_hit_high[islboard][ichip][isca][ichn]==0 && selection==true && gooddata==true)
	      ped_sca.at(islboard).at(ichip).at(ichn).at(isca)->Fill(charge_hiGain[islboard][ichip][isca][ichn]);

	    //bad events
	    //	  selection=false;
	    //	  if( ( badbcid[islboard][ichip][isca]!=0  || nhits[islboard][ichip][isca]>maxnhit || gooddata==false) ) selection=true;
	    //	  if(masked[islboard][ichip][ichn]==1) selection=false;
	    if(gain_hit_high[islboard][ichip][isca][ichn]==0 && (selection==false || gooddata==false) )
	      ped_sca_tagged.at(islboard).at(ichip).at(ichn).at(isca)->Fill(charge_hiGain[islboard][ichip][isca][ichn]);
	  }
           
	}//isca

      }//ichip 

    }// islboard 
  }  // end first loop analysis to fill pedestal historgrams

  
  TString slboard = TString::Format("_%i_slboards_",n_slboards);
  if(outputname!="") slboard=slboard+outputname;
  
  // --------------------------------------------------------------------------------
  // PEDESTAL ANALYSIS
  TFile *pedfile = new TFile("results_proto/Pedestal"+slboard+".root" , "RECREATE");
  pedfile->cd();
  TDirectory *cdhisto[nSLB];
  for(int islboard=0; islboard<nSLB; islboard++) {
    cdhisto[islboard] = pedfile->mkdir(TString::Format("slboard_%i",islboard));
  }
    
  for(int islboard=0; islboard<nSLB; islboard++) {
    std::vector<std::vector<std::vector<Double_t> > > chip_ped_mean_slb;
    std::vector<std::vector<std::vector<Double_t> > > chip_ped_error_slb;
    std::vector<std::vector<std::vector<Double_t> > > chip_ped_width_slb;

    //initialize pedestal vectors
    for(int i=0; i<16; i++) {
      std::vector<std::vector<Double_t> > chip_ped_mean;
      std::vector<std::vector<Double_t> > chip_ped_error;
      std::vector<std::vector<Double_t> > chip_ped_width;
      
      for(int j=0; j<64; j++) {
	std::vector<Double_t> chn_ped_mean;
	std::vector<Double_t> chn_ped_error;
	std::vector<Double_t> chn_ped_width;

	for(int isca=0; isca<15; isca++) {
	  chn_ped_mean.push_back(-1);
	  chn_ped_error.push_back(-1);
	  chn_ped_width.push_back(-1);
	}
	chip_ped_mean.push_back(chn_ped_mean);
	chip_ped_error.push_back(chn_ped_error);
	chip_ped_width.push_back(chn_ped_width);
      }  
      chip_ped_mean_slb.push_back(chip_ped_mean);
      chip_ped_error_slb.push_back(chip_ped_error);
      chip_ped_width_slb.push_back(chip_ped_width);
    }
    ped_mean_slboard.push_back(chip_ped_mean_slb);
    ped_error_slboard.push_back(chip_ped_error_slb);
    ped_width_slboard.push_back(chip_ped_width_slb);
  }

  // do pedestal (slboard/chip/channel/sca based) analysis
  for(int islboard=0; islboard<nSLB; islboard++) {
    cdhisto[islboard]->cd();
    
    for(int ichip=0; ichip<16; ichip++) {
      for(int ichn=0; ichn<64; ichn++) {
	
	for(int isca=0; isca<15; isca++) {
	  ped_sca.at(islboard).at(ichip).at(ichn).at(isca)->SetTitle(TString::Format("ped_chip%i_chn%i_sca%i",ichip,ichn,isca));
	  ped_sca.at(islboard).at(ichip).at(ichn).at(isca)->SetName(TString::Format("ped_chip%i_chn%i_sca%i",ichip,ichn,isca));
	  ped_sca.at(islboard).at(ichip).at(ichn).at(isca)->Write();
	  
	  ped_sca_tagged.at(islboard).at(ichip).at(ichn).at(isca)->SetTitle(TString::Format("ped_tagged_chip%i_chn%i_sca%i",ichip,ichn,isca));
	  ped_sca_tagged.at(islboard).at(ichip).at(ichn).at(isca)->SetName(TString::Format("ped_tagged_chip%i_chn%i_sca%i",ichip,ichn,isca));
	  ped_sca_tagged.at(islboard).at(ichip).at(ichn).at(isca)->Write();


	if(ped_sca.at(islboard).at(ichip).at(ichn).at(isca)->GetEntries()> 50 ){ //max_entries/2 ) {
	  TSpectrum *s = new TSpectrum();
	  int npeaks = s->Search(ped_sca.at(islboard).at(ichip).at(ichn).at(isca),2,"",0.8); 
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
	      
	      TF1 *f0 = new TF1("f0","gaus",mean_peak[0]-ped_sca.at(islboard).at(ichip).at(ichn).at(isca)->GetRMS(),mean_peak[0]+ped_sca.at(islboard).at(ichip).at(ichn).at(isca)->GetRMS());
	      ped_sca.at(islboard).at(ichip).at(ichn).at(isca)->Fit("f0","RQNOC");
	      
	      TF1 *f1 = new TF1("f1","gaus",f0->GetParameter(1)-f0->GetParameter(2),f0->GetParameter(1)+f0->GetParameter(2));
	      ped_sca.at(islboard).at(ichip).at(ichn).at(isca)->Fit("f1","RQME");
	 
	      
	      ped_mean_slboard.at(islboard).at(ichip).at(ichn).at(isca)=f1->GetParameter(1);
	      ped_error_slboard.at(islboard).at(ichip).at(ichn).at(isca)=f1->GetParError(1);
	      ped_width_slboard.at(islboard).at(ichip).at(ichn).at(isca)=f1->GetParameter(2);
	  
	      delete f0, f1;
	    }
	  }
	  delete s;
	}
	}
      }
    }
  }

  ped_sca.clear();
  ped_sca_tagged.clear();

  pedfile->Close();
  delete pedfile;
  

  // --------------------------------------------------------------------------------------------
  //****************************************************
  // MIP ANALYSIS

   //Objects to store the chip/channel results:
  // numbers 
  std::vector<std::vector<std::vector<TH1F*> > > mip_histo;
  std::vector<std::vector<std::vector<TH1F*> > > s_n_histo;

  for(int islboard=0; islboard<nSLB; islboard++) {
    std::vector<std::vector<TH1F*> > mip_histo_slboard;
    std::vector<std::vector<TH1F*> > s_n_histo_slboard;
    for(int ichip=0; ichip<16; ichip++) {
      std::vector<TH1F*>  tmpmip_histo;
      std::vector<TH1F*>  tmps_n_histo;
      
      for(int ichn=0;ichn<65;ichn++) {
	TString histo_title=TString::Format("charge_slboard%i_chip%i_chn%i",islboard,ichip,ichn);
	if(ichn==64) histo_title=TString::Format("charge_slboard%i_chip%i",islboard,ichip);
	TH1F *chn_mip = new TH1F(histo_title,histo_title,4197,-100.5,4096.5);
	tmpmip_histo.push_back(chn_mip);
	TString s_n_title=TString::Format("s_n_slboard%i_chip%i_chn%i",islboard,ichip,ichn);
	if(ichn==64) s_n_title=TString::Format("s_n_slboard%i_chip%i",islboard,ichip);
	TH1F *chn_s_n = new TH1F(s_n_title,s_n_title,4197,-100.5,4096.5);
	tmps_n_histo.push_back(chn_s_n);
      }
      mip_histo_slboard.push_back(tmpmip_histo);
      s_n_histo_slboard.push_back(tmps_n_histo);
      
    }
    mip_histo.push_back(mip_histo_slboard);
    s_n_histo.push_back(s_n_histo_slboard);
  }

  TH2F* h_trig[15];
  TH2F* h_trig_xy[15];
  
  TH2F* h_first_retriggering[15];
  TH2F* h_all_retriggering[15];
  TH2F* h_first_retriggering_xy[15];
  TH2F* h_all_retriggering_xy[15];
  
  for(int islboard=0; islboard<nSLB; islboard++) {
    h_trig[islboard] = new TH2F(TString::Format("trig_slboard_%i",islboard),TString::Format("trig_slboard_%i",islboard),16,-0.5,15.5,64,-0.5,63.5);
    h_trig_xy[islboard] = new TH2F(TString::Format("trig_xy_slboard_%i",islboard),TString::Format("trig_xy_slboard_%i",islboard),32,-90,90,32,-90,90);
    
    h_first_retriggering[islboard] = new TH2F(TString::Format("first_retriggering_slboard_%i",islboard),TString::Format("first_retriggering_slboard_%i",islboard),16,-0.5,15.5,64,-0.5,63.5);
    h_all_retriggering[islboard] = new TH2F(TString::Format("all_retriggering_slboard_%i",islboard),TString::Format("all_retriggering_slboard_%i",islboard),16,-0.5,15.5,64,-0.5,63.5);
    h_first_retriggering_xy[islboard] = new TH2F(TString::Format("first_retriggering_xy_slboard_%i",islboard),TString::Format("first_retriggering_xy_slboard_%i",islboard),32,-90,90,32,-90,90);
    h_all_retriggering_xy[islboard] = new TH2F(TString::Format("all_retriggering_xy_slboard_%i",islboard),TString::Format("all_retriggering_xy_slboard_%i",islboard),32,-90,90,32,-90,90);
  }
      
  // -----------------------------------------------------------------------------------------------------
  // Signal readout
  nbytes = 0;
  nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if ( jentry > 1000 && jentry % 1000 ==0 ) std::cout << "Progress: " << 100.*jentry/nentries <<" %"<<endl;
   
    for(int islboard=0; islboard<n_slboards; islboard++) {
      cout<<islboard<<endl;
      
      for(int ichip=0; ichip<16; ichip++) {
	  // RETRIGGER STUFF
	bool first_retrig=false;
      
      	for(int isca=0; isca<15; isca++) {
	  if(bcid[islboard][ichip][isca]<0) continue;
	  if(bcid[islboard][ichip][isca]<50 || (bcid[islboard][ichip][isca]>890 && bcid[islboard][ichip][isca]<930)) continue;

	  if(isca==0) {
	    if(badbcid[islboard][ichip][isca]>2) first_retrig=true;
	    for(int ichn=0; ichn<64; ichn++) {
	      if(gain_hit_high[islboard][ichip][isca][ichn]==1) {
		h_first_retriggering[islboard] -> Fill(ichip,ichn);
		h_all_retriggering[islboard] -> Fill(ichip,ichn);

		h_first_retriggering_xy[islboard] -> Fill(map_pointX[islboard][ichip][ichn],map_pointY[islboard][ichip][ichn]);
		h_all_retriggering_xy [islboard]-> Fill(map_pointX[islboard][ichip][ichn],map_pointY[islboard][ichip][ichn]);     
	      }
	    }
	  }
	
	  if(isca > 0 && badbcid[islboard][ichip][isca]>2 && first_retrig==false) {
	    first_retrig=true;
	    for(int ichn=0; ichn<64; ichn++) {
	      if(gain_hit_high[islboard][ichip][isca][ichn]==1) {
		h_first_retriggering[islboard] -> Fill(ichip,ichn);
		h_all_retriggering[islboard] -> Fill(ichip,ichn);

		h_first_retriggering_xy[islboard] -> Fill(map_pointX[islboard][ichip][ichn],map_pointY[islboard][ichip][ichn]);
		h_all_retriggering_xy[islboard] -> Fill(map_pointX[islboard][ichip][ichn],map_pointY[islboard][ichip][ichn]);    
	      }
	    }
	  }

	  if(badbcid[islboard][ichip][isca]>2 && first_retrig==true){
	    for(int ichn=0; ichn<64; ichn++) {
	      if(gain_hit_high[islboard][ichip][isca][ichn]==1) {
		h_all_retriggering[islboard] -> Fill(ichip,ichn);
		h_all_retriggering_xy[islboard] -> Fill(map_pointX[islboard][ichip][ichn],map_pointY[islboard][ichip][ichn]);    

	      }
	    }
	  }
	}
      
	for(int isca=0; isca<15; isca++) {
	  if(bcid[islboard][ichip][isca]<0) continue;
	  if(bcid[islboard][ichip][isca]<50 || (bcid[islboard][ichip][isca]>890 && bcid[islboard][ichip][isca]<930)) continue;
	  for(int ichn=0; ichn<64; ichn++) {

	    if( ped_mean_slboard.at(islboard).at(ichip).at(ichn).at(isca)>50 &&  ped_width_slboard.at(islboard).at(ichip).at(ichn).at(isca)>0 ) {

	      bool selection=false;
	      if(charge_hiGain[islboard][ichip][isca][ichn]>0 && gain_hit_high[islboard][ichip][isca][ichn]==1 && badbcid[islboard][ichip][isca]==0 && nhits[islboard][ichip][isca]<maxnhit+1 ) {
		if(bcid[islboard][ichip][isca]>15 ) selection=true;
	      }

	      if(selection==true) {
		h_trig[islboard] -> Fill(ichip,ichn);
		h_trig_xy[islboard] -> Fill(map_pointX[islboard][ichip][ichn],map_pointY[islboard][ichip][ichn]);   
		mip_histo.at(islboard).at(ichip).at(ichn)->Fill(charge_hiGain[islboard][ichip][isca][ichn]-ped_mean_slboard.at(islboard).at(ichip).at(ichn).at(isca));
		s_n_histo.at(islboard).at(ichip).at(ichn)->Fill( (charge_hiGain[islboard][ichip][isca][ichn]-ped_mean_slboard.at(islboard).at(ichip).at(ichn).at(isca)) / ped_width_slboard.at(islboard).at(ichip).at(ichn).at(isca));

		mip_histo.at(islboard).at(ichip).at(64)->Fill(charge_hiGain[islboard][ichip][isca][ichn]-ped_mean_slboard.at(islboard).at(ichip).at(ichn).at(isca));
		s_n_histo.at(islboard).at(ichip).at(64)->Fill( (charge_hiGain[islboard][ichip][isca][ichn]-ped_mean_slboard.at(islboard).at(ichip).at(ichn).at(isca)) / ped_width_slboard.at(islboard).at(ichip).at(ichn).at(isca));
	      }
	    }
	    
	  }//ichn
	  
	}//sca
      }//ichip
    }//islboard

  }

  ped_mean_slboard.clear();
  ped_error_slboard.clear();
  ped_width_slboard.clear();
  
  
  // -----------------------------------------------------------------------------------------------------   
  // Signal analysis
  TFile *signalfile_summary = new TFile("results_proto/MIPs"+slboard+".root" , "RECREATE");
  signalfile_summary->cd();
  TDirectory *cdhistomip[nSLB];
  for(int islboard=0; islboard<nSLB; islboard++) {
    cdhistomip[islboard] = signalfile_summary->mkdir(TString::Format("slboard_%i",islboard));
  }

  // Setting fit range and start values
  Double_t fr[2];
  Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];   
  pllo[0]=0.; pllo[1]=0.0; pllo[2]=10.0; pllo[3]=0.;
  plhi[0]=100.0; plhi[1]=300.0; plhi[2]=100000000.0; plhi[3]=20.0;
  sv[0]=15.0;
  Double_t chisqr;
  Int_t    ndf;

  //   par[0]=Width (scale) parameter of Landau density
  //   par[1]=Most Probable (MP, location) parameter of Landau density
  //   par[2]=Total area (integral -inf to inf, normalization constant)
  //   par[3]=Width (sigma) of convoluted Gaussian function
  //

  for(int islboard=0; islboard<nSLB; islboard++) {
    cout<<"  MIP ANALYSIS OF SLBOARD: "<<islboard<<endl;;
    cdhistomip[islboard]->cd();
    h_trig[islboard]->Write();
    h_trig_xy[islboard]->Write();
    
    h_first_retriggering[islboard]->Write();
    h_all_retriggering[islboard]->Write();
    h_first_retriggering_xy[islboard]->Write();
    h_all_retriggering_xy[islboard]->Write();
        
    for(int ichip=0; ichip<16; ichip++) {
      cout<<"     SK: "<<ichip << " "<<endl;;

      for(int ichn=0; ichn<65; ichn++) {
	if(mip_histo.at(islboard).at(ichip).at(ichn)->GetEntries()>10){

	  fr[0]=mip_histo.at(islboard).at(ichip).at(ichn)->GetMean()-0.8*mip_histo.at(islboard).at(ichip).at(ichn)->GetRMS();
	  if(fr[0]<30) fr[0]=30;
	  fr[1]=fr[0]+0.5*mip_histo.at(islboard).at(ichip).at(ichn)->GetRMS();
	  if(fr[0]==30) fr[1]=fr[0]+mip_histo.at(islboard).at(ichip).at(ichn)->GetRMS();
	  sv[0]=mip_histo.at(islboard).at(ichip).at(ichn)->GetRMS()*0.5;
	  sv[1]=mip_histo.at(islboard).at(ichip).at(ichn)->GetMean()*0.6;
	  sv[2]=mip_histo.at(islboard).at(ichip).at(ichn)->Integral("width");
	  sv[3]=mip_histo.at(islboard).at(ichip).at(ichn)->GetRMS()/5.;
	  
	  TF1 *fitsnr_temp=langaufit(mip_histo.at(islboard).at(ichip).at(ichn),fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
	  mip_histo.at(islboard).at(ichip).at(ichn)->Write();
	  
	  double mpv=fitsnr_temp->GetParameter(1);
	  double empv=fitsnr_temp->GetParError(1);
	  double wmpv=fitsnr_temp->GetParameter(0);
	  double chi2ndf=0;
	  if(ndf>0) chi2ndf=chisqr/ndf;
	  double mipentries=mip_histo.at(islboard).at(ichip).at(ichn)->GetEntries();
	  
	  
	  /// Signal over NOISE:
	  fr[0]=s_n_histo.at(islboard).at(ichip).at(ichn)->GetMean()-0.8*s_n_histo.at(islboard).at(ichip).at(ichn)->GetRMS();
	  fr[1]=s_n_histo.at(islboard).at(ichip).at(ichn)->GetMean();//+0.1*s_n_histo.at(islboard).at(ichip).at(ichn)->GetRMS();                                                                                                         
	  sv[0]=s_n_histo.at(islboard).at(ichip).at(ichn)->GetRMS()*0.5;
	  sv[1]=s_n_histo.at(islboard).at(ichip).at(ichn)->GetMean()*0.6;
	  sv[2]=s_n_histo.at(islboard).at(ichip).at(ichn)->Integral("width");
	  sv[3]=s_n_histo.at(islboard).at(ichip).at(ichn)->GetRMS()/5.;
	  
	  TF1 *fitsnr_temp2=langaufit(s_n_histo.at(islboard).at(ichip).at(ichn),fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
	  s_n_histo.at(islboard).at(ichip).at(ichn)->Write();


	  delete fitsnr_temp;
	  delete fitsnr_temp2;

	}//nentries>10
      }//ichn
    }//ichip
  }//islboard

  cout<<"END"<<endl;
  mip_histo.clear();
  s_n_histo.clear();

  signalfile_summary->Close();
  delete signalfile_summary;

  delete  tempcanvas;
  return 1;

}

void DecodedSLBAnalysis::SynchronizationStudies(TString outputname="testMonitoring", int freq=10, bool shifter=false)
{
  TH2F* adc_bcid[15][16];

  for(int i=0; i<15; i++) {
    for(int j=0; j<16; j++) {
      adc_bcid[i][j]=new TH2F(TString::Format("ad_bcid_slb%i_chip%i",i,j),TString::Format("ad_bcid_slb%i_chip%i",i,j),1000,-0.5,999.5,1000,-0.5,999.5);
    }
  }

  // --------------
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //-----------------
  cout<<"Start SynchronizationStudies, read only 1 event each "<< freq<<endl;
  cout<<"Total number of entries: "<< nentries<<endl;

  // -----------------------------------------------------------------------------------------------------
  // Signal readout
  Long64_t nbytes = 0, nb = 0;
  int n_SLB=0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if ( jentry > 1000 && jentry % 1000 ==0 ) std::cout << "Progress: " << 100.*jentry/nentries <<" %"<<endl;

    if(jentry==0) n_SLB=n_slboards;

    if(jentry % freq !=0 ) continue;

    for(int islboard=0; islboard<n_slboards; islboard++) {
      //channels starting retriggers
      for(int ichip=0; ichip<16; ichip++) {
	
	bool first_retrig=false;
	bool trigger=false;
	int isca_trig=-1;
	int isca_retrig=-1;
	for(int isca=0; isca<15; isca++) {
	 for(int ichn=0; ichn<64; ichn++) {
	   if(gain_hit_high[islboard][ichip][isca][ichn]==1) {
	     adc_bcid[islboard][ichip]->Fill(bcid[islboard][ichip][isca],charge_hiGain[islboard][ichip][isca][ichn]);
	   }
	 }
	}
      }
    }

  

  }
  // -----------------------------------------------------------------------------------------------------   
  // Signal analysis
  TFile *monitoringfile_summary = new TFile("results_monitoring/SynchronizationStudies_summary_"+outputname+".root" , "RECREATE");
  monitoringfile_summary->cd();

  for(int i=0; i<n_SLB; i++) {
    for(int j=0; j<16; j++) {
      adc_bcid[i][j]->GetXaxis()->SetTitle(TString::Format("adc_bcid slboard_%i_chip%i",i,i));
      adc_bcid[i][j]->GetYaxis()->SetTitle(TString::Format("adc_bcid slboard_%i_chip%i",i,j));
      adc_bcid[i][j]->Write();
    }
  } 
  monitoringfile_summary->Close();

}


void DecodedSLBAnalysis::Monitoring(TString outputname="testMonitoring", int freq=10, bool shifter=false)
{
  TH2F* bcid_diff_trig_trig[15];
  TH2F* bcid_diff_trig_retrig_start[15];
  TH2F* bcid_diff[15];

  TH2F* bcid_correl[15][15];

  TH2F* hitmap_trig_xy[15];
  TH2F* hitmap_retrig_start_xy[15];
  TH2F* hitmap_retrig_xy[15];

  TH2F* bcidmap_trig[15];
  TH2F* bcidmap_retrig_start[15];
  TH2F* bcidmap_retrig[15];

  TH2F* hitmap_trig_chipchn[15];
  TH2F* hitmap_retrig_start_chipchn[15];
  TH2F* hitmap_retrig_chipchn[15];

  TH2F* lastsca_trig[15];
  TH2F* lastsca_retrig_start[15];
  TH2F* lastsca_retrig[15];

  for(int i=0; i<15; i++) {
    bcid_diff_trig_trig[i]=new TH2F(TString::Format("bcid_diff_trig_trig_slb%i",i),TString::Format("bcid_diff_trig_trig_slb%i",i),15,-0.5,14.5,2001,-1000.5,1000.5);
    bcid_diff_trig_retrig_start[i]=new TH2F(TString::Format("bcid_diff_trig_retrig_start_slb%i",i),TString::Format("bcid_diff_trig_retrig_start_slb%i",i),15,-0.5,14.5,2001,-1000.5,1000.5);
    bcid_diff[i]=new TH2F(TString::Format("bcid_dif_slb%i",i),TString::Format("bcid_diff_slb%i",i),15,-0.5,14.5,2001,-1000.5,1000.5);

    for(int j=0; j<15; j++) {
      bcid_correl[i][j]=new TH2F(TString::Format("bcid_correl_slb%i_slb%i",i,j),TString::Format("bcid_correl_slb%i_slb%i",i,j),500,0,5000,500,0,5000);
    }

    hitmap_trig_xy[i]=new TH2F(TString::Format("hitmap_trig_xy_slb%i",i),TString::Format("hitmap_trig_xy_slb%i",i),32,-90,90,32,-90,90);
    hitmap_retrig_start_xy[i]=new TH2F(TString::Format("hitmap_retrig_start_xy_slb%i",i),TString::Format("hitmap_retrig_start_xy_slb%i",i),32,-90,90,32,-90,90);
    hitmap_retrig_xy[i]=new TH2F(TString::Format("hitmap_retrig_xy_slb%i",i),TString::Format("hitmap_retrig_xy_slb%i",i),32,-90,90,32,-90,90);

    bcidmap_trig[i]=new TH2F(TString::Format("bcidmap_trig_slb%i",i),TString::Format("bcidmap_trig_slb%i",i),820,-2.5,4097.5,820,-2.5,4097.5);
    bcidmap_retrig_start[i]=new TH2F(TString::Format("bcidmap_retrig_start_slb%i",i),TString::Format("bcidmap_retrig_start_slb%i",i),820,-2.5,4097.5,820,-2.5,4097.5);
    bcidmap_retrig[i]=new TH2F(TString::Format("bcidmap_retrig_slb%i",i),TString::Format("bcidmap_retrig_slb%i",i),820,-2.5,4097.5,820,-2.5,4097.5);

    hitmap_trig_chipchn[i]=new TH2F(TString::Format("hitmap_trig_chipchn_slb%i",i),TString::Format("hitmap_trig_chipchn_slb%i",i),16,-0.5,15.5,64,-0.5,63.5);
    hitmap_retrig_start_chipchn[i]=new TH2F(TString::Format("hitmap_retrig_start_chipchn_slb%i",i),TString::Format("hitmap_retrig_start_chipchn_slb%i",i),16,-0.5,15.5,64,-0.5,63.5);
    hitmap_retrig_chipchn[i]=new TH2F(TString::Format("hitmap_retrig_chipchn_slb%i",i),TString::Format("hitmap_retrig_chipchn_slb%i",i),16,-0.5,15.5,64,-0.5,63.5);

    lastsca_trig[i]=new TH2F(TString::Format("lastsca_trig_slb%i",i),TString::Format("lastsca_trig_slb%i",i),16,-0.5,15.5,15,-0.5,14.5);
    lastsca_retrig[i]=new TH2F(TString::Format("lastsca_retrig_slb%i",i),TString::Format("lastsca_retrig_slb%i",i),16,-0.5,15.5,15,-0.5,14.5);
    lastsca_retrig_start[i]=new TH2F(TString::Format("FIRSTsca_retrig_start_slb%i",i),TString::Format("FIRSTsca_retrig_start_slb%i",i),16,-0.5,15.5,15,-0.5,14.5);

  }
 
  TH2F* trig = new TH2F("trig","trig",14,-0.5,13.5,16,-0.5,15.5);
  TH2F* retrig_start = new TH2F("retrig_start","retrig_start",14,-0.5,13.5,16,-0.5,15.5);
  TH2F* retrig = new TH2F("retrig","retrig",14,-0.5,13.5,16,-0.5,15.5);

  /*  TH1F* trig_coinc = new TH1F("trig_coinc","trig_cinc",14,0.5,14.5);
  TH1F* retrig_start_coinc = new TH1F("retrig_start_coinc","retrig__startcinc",14,0.5,14.5);
  TH1F* retrig_coinc = new TH1F("retrig_coinc","retrig_cinc",14,0.5,14.5);*/

  // --------------
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //-----------------
  cout<<"Start MONITORING, read only 1 event each "<< freq<<endl;
  cout<<"Total number of entries: "<< nentries<<endl;

  // -----------------------------------------------------------------------------------------------------
  // Signal readout
  Long64_t nbytes = 0, nb = 0;
  int n_SLB=0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if ( jentry > 1000 && jentry % 1000 ==0 ) std::cout << "Progress: " << 100.*jentry/nentries <<" %"<<endl;

    if(jentry==0) n_SLB=n_slboards;

    if(jentry % freq !=0 ) continue;

    for(int islboard=0; islboard<n_slboards; islboard++) {
      //channels starting retriggers
      for(int ichip=0; ichip<16; ichip++) {
	
	bool first_retrig=false;
	bool trigger=false;
	int isca_trig=-1;
	int isca_retrig=-1;
	for(int isca=0; isca<15; isca++) {

	  trigger=false;
	  //first retrigger
	  if(badbcid[islboard][ichip][isca]>2 && first_retrig==false) {
	    first_retrig=true;
	    for(int ichn=0; ichn<64; ichn++) {
	     if(gain_hit_high[islboard][ichip][isca][ichn]==1) {
	       hitmap_retrig_start_xy[islboard]->Fill(map_pointX[islboard][ichip][ichn],map_pointY[islboard][ichip][ichn]);
	       hitmap_retrig_start_chipchn[islboard]->Fill(ichip,ichn);
	     }
	    }
	    retrig_start->Fill(islboard,ichip);
	    lastsca_retrig_start[islboard]->Fill(ichip,isca);
	  }

	  //all retriggers
	  if(badbcid[islboard][ichip][isca]>2 && first_retrig==true) {
	    for(int ichn=0; ichn<64; ichn++) {
	     if(gain_hit_high[islboard][ichip][isca][ichn]==1) {
	       hitmap_retrig_xy[islboard]->Fill(map_pointX[islboard][ichip][ichn],map_pointY[islboard][ichip][ichn]);
	       hitmap_retrig_chipchn[islboard]->Fill(ichip,ichn);
	     }
	    }
	    retrig->Fill(islboard,ichip);
	    isca_retrig=isca;
	  }

	  //all triggers
	  if(badbcid[islboard][ichip][isca]==0) {
	    for(int ichn=0; ichn<64; ichn++) {
	     if(gain_hit_high[islboard][ichip][isca][ichn]==1) {
	       hitmap_trig_xy[islboard]->Fill(map_pointX[islboard][ichip][ichn],map_pointY[islboard][ichip][ichn]);
	       hitmap_trig_chipchn[islboard]->Fill(ichip,ichn);
	     }
	    }
	    trig->Fill(islboard,ichip);
	    trigger=true;
	    isca_trig=isca;
	  }

	  // search coincidences
	  //-------------------------------------------------------
	  for(int islboard2=0; islboard2<n_slboards; islboard2++) {
	    //channels starting retriggers
	    for(int ichip2=0; ichip2<16; ichip2++) {
	      bool first_retrig2=false;
	      for(int isca2=0; isca2<15; isca2++) {
		if(bcid[islboard2][ichip2][isca2]<0) continue;
		bcid_correl[islboard][islboard2]->Fill(bcid[islboard][ichip][isca],bcid[islboard2][ichip2][isca2]);
		  bcid_diff[islboard]->Fill(islboard2,bcid[islboard2][ichip2][isca2]-bcid[islboard][ichip][isca]);
		  if(trigger==true) {	      
		    if(islboard2!=islboard || ichip2!=ichip || isca2>isca)  {
		      //first retrigger
		      if(badbcid[islboard2][ichip2][isca2]>2 && first_retrig2==false) {
			first_retrig2=true;
			bcid_diff_trig_retrig_start[islboard]->Fill(islboard2,bcid[islboard2][ichip2][isca2]-bcid[islboard][ichip][isca]);		      
		      }
		      if(badbcid[islboard2][ichip2][isca2]==0) {
			bcid_diff_trig_trig[islboard]->Fill(islboard2,bcid[islboard2][ichip2][isca2]-bcid[islboard][ichip][isca]);
		      }
		    }
		  }//trig==true
	      }//isca2
	    }//ichip2
	  }//islboard2
	    
      }//end isca
	lastsca_trig[islboard]->Fill(ichip,isca_trig);
	lastsca_retrig[islboard]->Fill(ichip,isca_retrig);
	
    }//end chip
  }// islboard
  

  }
  // -----------------------------------------------------------------------------------------------------   
  // Signal analysis
  TFile *monitoringfile_summary = new TFile("results_monitoring/Monitoring_summary_"+outputname+".root" , "RECREATE");
  monitoringfile_summary->cd();

 for(int i=0; i<n_SLB; i++)
   for(int j=0; j<n_SLB; j++) {
     bcid_correl[i][j]->GetXaxis()->SetTitle(TString::Format("bcid slboard_%i",i));
     bcid_correl[i][j]->GetYaxis()->SetTitle(TString::Format("bcid slboard_%i",j));
     bcid_correl[i][j]->Write();
   }
 
  for(int i=0; i<n_SLB; i++) {

    lastsca_trig[i]->GetXaxis()->SetTitle("CHIP");
    lastsca_trig[i]->GetYaxis()->SetTitle("last SCA");
    lastsca_trig[i]->Write();

    lastsca_retrig[i]->GetXaxis()->SetTitle("CHIP");
    lastsca_retrig[i]->GetYaxis()->SetTitle("last SCA");
    lastsca_retrig[i]->Write();

    lastsca_retrig_start[i]->GetXaxis()->SetTitle("CHIP");
    lastsca_retrig_start[i]->GetYaxis()->SetTitle("last SCA");
    lastsca_retrig_start[i]->Write();

    bcid_diff_trig_trig[i]->SetTitle(TString::Format("correlations with slboard_%i triggers",i));
    bcid_diff_trig_trig[i]->GetYaxis()->SetTitle(TString::Format("bcid_{trigger in slboard_%i}- bcid_{trigger in slboard_xaxis}",i));
    bcid_diff_trig_trig[i]->GetXaxis()->SetTitle("SLB");
    // bcid_diff_trig_trig[i]->GetYaxis()->SetRangeUser(-100,100);
    bcid_diff_trig_trig[i]->Write();

    bcid_diff_trig_retrig_start[i]->SetTitle(TString::Format("correlations with slboard_%i triggers",i));
    bcid_diff_trig_retrig_start[i]->GetYaxis()->SetTitle(TString::Format("bcid_{trigger in slboard_%i}- bcid_{first REtrigger in slboard_xaxis}",i));
    bcid_diff_trig_retrig_start[i]->GetXaxis()->SetTitle("SLB");
    // bcid_diff_trig_retrig_start[i]->GetYaxis()->SetRangeUser(-100,100);
    bcid_diff_trig_retrig_start[i]->Write();

    bcid_diff[i]->SetTitle(TString::Format("correlations with slboard_%i",i));
    bcid_diff[i]->GetYaxis()->SetTitle(TString::Format("bcid_{all in slboard_%i}- bcid_{all in slboard_xaxis}",i));
    bcid_diff[i]->GetXaxis()->SetTitle("SLB");
    // bcid_diff_trig_retrig_start[i]->GetYaxis()->SetRangeUser(-100,100);
    bcid_diff[i]->Write();
    
    hitmap_trig_xy[i]->GetYaxis()->SetTitle("x");
    hitmap_trig_xy[i]->GetYaxis()->SetTitle("y");
    hitmap_trig_xy[i]->Write();

    hitmap_retrig_start_xy[i]->GetYaxis()->SetTitle("x");
    hitmap_retrig_start_xy[i]->GetYaxis()->SetTitle("y");
    hitmap_retrig_start_xy[i]->Write();

    hitmap_retrig_xy[i]->GetYaxis()->SetTitle("x");
    hitmap_retrig_xy[i]->GetYaxis()->SetTitle("y");
    hitmap_retrig_xy[i]->Write();

    hitmap_trig_chipchn[i]->GetYaxis()->SetTitle("CHIP");
    hitmap_trig_chipchn[i]->GetYaxis()->SetTitle("chn");
    hitmap_trig_chipchn[i]->Write();

    hitmap_retrig_start_chipchn[i]->GetYaxis()->SetTitle("CHIP");
    hitmap_retrig_start_chipchn[i]->GetYaxis()->SetTitle("chn");
    hitmap_retrig_start_chipchn[i]->Write();

    hitmap_retrig_chipchn[i]->GetYaxis()->SetTitle("CHIP");
    hitmap_retrig_chipchn[i]->GetYaxis()->SetTitle("chn");
    hitmap_retrig_chipchn[i]->Write();
   
  }

  trig->Write();
  retrig_start->Write();
  retrig->Write();

  TCanvas *c1 = new TCanvas("hitmaps_xy","hitmaps_xy",1800,900);
  c1->Divide(4,4);
  for(int i=0; i<n_SLB; i++) {
    c1->cd(i+1);
    hitmap_trig_xy[i]->Draw("colz");
  }
  c1->Write();
  c1->Modified();
  c1->Update();
  if(shifter==true) c1->WaitPrimitive();

  TCanvas *c1_ret = new TCanvas("retrigger_start_xy","retrigger_start_xy",1800,900);
  c1_ret->Divide(4,4);
  for(int i=0; i<n_SLB; i++) {
    c1->cd(i+1);
    hitmap_retrig_start_xy[i]->Draw("colz");
  }
  c1_ret->Write();
  c1_ret->Modified();
  c1_ret->Update();
  if(shifter==true) c1_ret->WaitPrimitive();


  TCanvas *c1_ret2 = new TCanvas("retrigger_xy","retrigger_xy",1800,900);
  c1_ret2->Divide(4,4);
  for(int i=0; i<n_SLB; i++) {
    c1->cd(i+1);
    hitmap_retrig_start_xy[i]->Draw("colz");
  }
  c1_ret2->Write();
  c1_ret2->Modified();
  c1_ret2->Update();
  if(shifter==true) c1_ret2->WaitPrimitive();


  TCanvas *c2 = new TCanvas("hitmaps_chipchn","hitmaps_chipchn",1800,900);
  c2->Divide(4,4);
  for(int i=0; i<n_SLB; i++) {
    c2->cd(i+1);
    hitmap_trig_chipchn[i]->Draw("colz");
  }
  c2->Write();
  c2->Modified();
  c2->Update();
  if(shifter==true) c2->WaitPrimitive();

  TCanvas *c2_ret = new TCanvas("retrigger_start_chipchn","retrigger_start_chipchn",1800,900);
  c2_ret->Divide(4,4);
  for(int i=0; i<n_SLB; i++) {
    c2->cd(i+1);
    hitmap_retrig_start_chipchn[i]->Draw("colz");
  }
  c2_ret->Write();
  c2_ret->Modified();
  c2_ret->Update();
  if(shifter==true) c2_ret->WaitPrimitive();


  TCanvas *c2_ret2 = new TCanvas("retrigger_chipchn","retrigger_chipchn",1800,900);
  c2_ret2->Divide(4,4);
  for(int i=0; i<n_SLB; i++) {
    c2->cd(i+1);
    hitmap_retrig_start_chipchn[i]->Draw("colz");
  }
  c2_ret2->Write();
  c2_ret2->Modified();
  c2_ret2->Update();
  if(shifter==true) c2_ret2->WaitPrimitive();


  TCanvas *c1bis = new TCanvas("lastsca","lastsca",1800,900);
  c1bis->Divide(n_SLB,3);
  for(int i=0; i<n_SLB; i++) {
    c1bis->cd(i+1);
    lastsca_trig[i]->Draw("colz");
    c1bis->cd(i+1+n_SLB);
    lastsca_retrig_start[i]->Draw("colz");
    c1bis->cd(i+1+2*n_SLB);
    lastsca_retrig[i]->Draw("colz");
  }
  
  c1bis->Write();
  c1bis->Modified();
  c1bis->Update();
  if(shifter==true) c1bis->WaitPrimitive();

  TCanvas *c3 = new TCanvas("bcid_correl","bcid_correl",1800,900);
  c3->Divide(n_SLB,n_SLB);
  for(int i=0; i<n_SLB; i++) {
    for(int j=0; j<n_SLB; j++) {
      c3->cd(i+1+j*n_SLB);
      bcid_correl[i][j]->Draw("colz");
    }
  }

  c3->Write();
  c3->Modified();
  c3->Update();
  if(shifter==true) c3->WaitPrimitive();


  TCanvas *c4 = new TCanvas("bcid_correl2","bcid_correl2",1800,900);    
  c4->Divide(4,4);
  for(int i=0; i<n_SLB; i++) {
    c4->cd(i+1);
    bcid_diff_trig_trig[i]->GetYaxis()->SetRangeUser(-100,100);
    bcid_diff_trig_trig[i]->Draw("colz");
  }

  c4->Write();
  c4->Modified();
  c4->Update();
  if(shifter==true) c4->WaitPrimitive();

  TCanvas *c5 = new TCanvas("bcid_correl3","bcid_correl3",1800,900);
  c5->Divide(4,4);
  for(int i=0; i<n_SLB; i++) {
    c5->cd(i+1);
    bcid_diff_trig_retrig_start[i]->GetYaxis()->SetRangeUser(-100,100);
    bcid_diff_trig_retrig_start[i]->Draw("colz");
  }

  c5->Write();
  c5->Modified();
  c5->Update();
  if(shifter==true) c5->WaitPrimitive();
  
  monitoringfile_summary->Close();

}


void DecodedSLBAnalysis::SignalAnalysis(int i_slboard, TString outputname="", int maxnhit=1, int ncoinc=7)
{

  TString slboard = TString::Format("_slboard_%i_",slboard_array_mapping[i_slboard]);
  if(outputname!="") slboard=slboard+outputname;

  //Read the list of pedestals (this information contains, implicitily, the masked channels information )
  ReadPedestals("results_pedestal/Pedestal"+slboard+".txt");

  
  ofstream fout_mip("results_mipcalibration/MIP"+slboard+".txt",ios::out);
  fout_mip<<"#mip results "<<slboard<<endl;
  fout_mip<<"#chip channel mpv empv widthmpv chi2ndf nentries"<<endl;

  // --------------
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //-----------------

  //Objects to store the chip/channel results:
  // numbers 
  std::vector<std::vector<TH1F*> > mip_histo;
  std::vector<std::vector<TH1F*> > s_n_histo;

  for(int ichip=0; ichip<16; ichip++) {
    std::vector<TH1F*>  tmpmip_histo;
    std::vector<TH1F*>  tmps_n_histo;
	
    for(int ichn=0;ichn<65;ichn++) {
      TString histo_title=TString::Format("charge_chip%i_chn%i",ichip,ichn);
      if(ichn==64) histo_title=TString::Format("charge_chip%i",ichip);
      TH1F *chn_mip = new TH1F(histo_title,histo_title,4197,-100.5,4096.5);
      tmpmip_histo.push_back(chn_mip);
      TString s_n_title=TString::Format("s_n_chip%i_chn%i",ichip,ichn);
      if(ichn==64) s_n_title=TString::Format("s_n_chip%i",ichip);
      TH1F *chn_s_n = new TH1F(s_n_title,s_n_title,4197,-100.5,4096.5);
      tmps_n_histo.push_back(chn_s_n);
    }      
    mip_histo.push_back(tmpmip_histo);
    s_n_histo.push_back(tmps_n_histo);
  
}
  //maps 
  TH2F* MPV_2d = new TH2F("MPV_2d","MPV_2d",32,-90,90,32,-90,90);
  TH2F* eMPV_2d = new TH2F("eMPV_2d","eMPV_2d",32,-90,90,32,-90,90);
  TH2F* widthMPV_2d = new TH2F("widthMPV_2d","widthMPV_2d",32,-90,90,32,-90,90);
  TH2F* chi2NDF_2d = new TH2F("chi2NDF_2d","chi2NDF_2d",32,-90,90,32,-90,90);
  TH2F* mipEntries_2d = new TH2F("mipEntries_2d","mipEntries_2d",32,-90,90,32,-90,90);
  TH2F* S_N_2d = new TH2F("S_N_2d","S_N_2d",32,-90,90,32,-90,90);

  TH2F* MPV_chipchn = new TH2F("MPV_chipchn","MPV_chipchn",64,-0.5,63.5,16,-0.5,15.5);
  TH2F* eMPV_chipchn = new TH2F("eMPV_chipchn","eMPV_chipchn",64,-0.5,63.5,16,-0.5,15.5);
  TH2F* widthMPV_chipchn = new TH2F("widthMPV_chipchn","widthMPV_chipchn",64,-0.5,63.5,16,-0.5,15.5);
  TH2F* mipEntries_chipchn = new TH2F("mipEntries_chipchn","mipEntries_chipchn",64,-0.5,63.5,16,-0.5,15.5);

  TH2F* correl_S_N_vs_nhits = new TH2F("correl_S_N_vs_nhits","correl_S_N_vs_nhits",80,0.25,40.25,20,1225,51225);
  TH2F* correl_S_N_vs_mpv = new TH2F("correl_S_N_vs_mpv","correl_S_N_vs_mpv",80,0.25,40.25,40,40.5,80.5);
  TH2F* correl_mpv_vs_nhits = new TH2F("correl_mpv_vs_nhits","correl_mpv_vs_nhits",40,40.5,80.5,20,1225,51225);
  TH2F* correl_widthped_vs_nhits = new TH2F("correl_widthped_vs_nhits","correl_widthped_vs_nhits",50,2.05,7.05,20,1225,51225);
  TH2F* correl_S_N_vs_widthmpv = new TH2F("correl_S_N_vs_widthmpv","correl_S_N_vs_widthmpv",80,0.25,40.25,50,0.5,50.5);
  TH2F* correl_S_N_vs_widthped = new TH2F("correl_S_N_vs_widthped","correl_S_N_vs_widthped",80,0.25,40.25,50,2.05,7.05);


  //summary plots
  TH1F* MPV_chip[16];
  TH1F* S_N_chip[16];
  for(int ichip=0; ichip<16; ichip++) {
    MPV_chip[ichip]= new TH1F(TString::Format("MPV_chip%i",ichip),TString::Format("MPV_chip%i",ichip),2000,0.05,200.05);
    S_N_chip[ichip]= new TH1F(TString::Format("S_N_chip%i",ichip),TString::Format("S_N_chip%i",ichip),500,0.05,50.05);
  }

  TH1F* MPV_slab= new TH1F("MPV","MPV",2000,0.05,200.05);
  TH1F* S_N_slab= new TH1F("S_N","S_N",1000,0.05,100.05);
 
  // -----------------------------------------------------------------------------------------------------
  // Signal readout
  Long64_t nbytes = 0, nb = 0;
  //  nentries=1000;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if ( jentry > 1000 && jentry % 1000 ==0 ) std::cout << "Progress: " << 100.*jentry/nentries <<" %"<<endl;

    std::array<int,4096> bcid_seen = SimpleCoincidenceTagger(jentry,maxnhit);

    for(int islboard=0; islboard<n_slboards; islboard++) {
      if(islboard != i_slboard) continue;
      for(int ichip=0; ichip<16; ichip++) {
	for(int isca=0; isca<15; isca++) {
	  
	  if(bcid[islboard][ichip][isca]<0 ) break;
	  
	  bool gooddata=true;
	  if( bcid_seen.at(bcid[islboard][ichip][isca]) < ncoinc ) gooddata=false;
	  for(int ichn=0; ichn<64; ichn++) if(charge_hiGain[islboard][ichip][isca][ichn]<100) gooddata=false;
	  
	  for(int ichn=0; ichn<64; ichn++) {

	    if( gooddata == true && ped_mean.at(ichip).at(ichn).at(isca)>50 &&  ped_width.at(ichip).at(ichn).at(isca)>0  && badbcid[islboard][ichip][isca]==0 && nhits[islboard][ichip][isca]<maxnhit+1 && bcid[islboard][ichip][isca]>15 ) {
	      if(charge_hiGain[islboard][ichip][isca][ichn]>0 && gain_hit_high[islboard][ichip][isca][ichn]==1) {
		mip_histo.at(ichip).at(ichn)->Fill(charge_hiGain[islboard][ichip][isca][ichn]-ped_mean.at(ichip).at(ichn).at(isca));
		s_n_histo.at(ichip).at(ichn)->Fill( (charge_hiGain[islboard][ichip][isca][ichn]-ped_mean.at(ichip).at(ichn).at(isca)) / ped_width.at(ichip).at(ichn).at(isca));
		
		mip_histo.at(ichip).at(64)->Fill(charge_hiGain[islboard][ichip][isca][ichn]-ped_mean.at(ichip).at(ichn).at(isca));
		s_n_histo.at(ichip).at(64)->Fill( (charge_hiGain[islboard][ichip][isca][ichn]-ped_mean.at(ichip).at(ichn).at(isca)) / ped_width.at(ichip).at(ichn).at(isca));
	      }
	    }
    
	  }//ichn
	}//sca
      }//ichip
    }//islboard
  }

  // -----------------------------------------------------------------------------------------------------   
  // Signal analysis
  TFile *signalfile_summary = new TFile("results_mipcalibration/Signal_summary"+slboard+".root" , "RECREATE");
  signalfile_summary->cd();
  TDirectory *cdhisto = signalfile_summary->mkdir("histograms");


  // Setting fit range and start values
  Double_t fr[2];
  Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];   
  pllo[0]=0.; pllo[1]=0.0; pllo[2]=10.0; pllo[3]=0.;
  plhi[0]=100.0; plhi[1]=300.0; plhi[2]=100000000.0; plhi[3]=20.0;
  sv[0]=15.0;
  Double_t chisqr;
  Int_t    ndf;

  //   par[0]=Width (scale) parameter of Landau density
  //   par[1]=Most Probable (MP, location) parameter of Landau density
  //   par[2]=Total area (integral -inf to inf, normalization constant)
  //   par[3]=Width (sigma) of convoluted Gaussian function
  //

  cdhisto->cd(); 

  for(int ichip=0; ichip<16; ichip++) {
    for(int ichn=0; ichn<65; ichn++) {
            
      if(mip_histo.at(ichip).at(ichn)->GetEntries()>10){

	fr[0]=mip_histo.at(ichip).at(ichn)->GetMean()-0.8*mip_histo.at(ichip).at(ichn)->GetRMS();
	if(fr[0]<30) fr[0]=30;
	fr[1]=fr[0]+0.5*mip_histo.at(ichip).at(ichn)->GetRMS();
	if(fr[0]==30) fr[1]=fr[0]+mip_histo.at(ichip).at(ichn)->GetRMS();
	sv[0]=mip_histo.at(ichip).at(ichn)->GetRMS()*0.5;
	sv[1]=mip_histo.at(ichip).at(ichn)->GetMean()*0.6;
	sv[2]=mip_histo.at(ichip).at(ichn)->Integral("width");
	sv[3]=mip_histo.at(ichip).at(ichn)->GetRMS()/5.;
	
	TF1 *fitsnr_temp=langaufit(mip_histo.at(ichip).at(ichn),fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
	mip_histo.at(ichip).at(ichn)->Write();
	
	double mpv=fitsnr_temp->GetParameter(1);
	double empv=fitsnr_temp->GetParError(1);
	double wmpv=fitsnr_temp->GetParameter(0);
	double chi2ndf=0;
	if(ndf>0) chi2ndf=chisqr/ndf;
	double mipentries=mip_histo.at(ichip).at(ichn)->GetEntries();

	if(ichn<64) {
	  MPV_2d->Fill(map_pointX[i_slboard][ichip][ichn],map_pointY[i_slboard][ichip][ichn],mpv);
	  eMPV_2d->Fill(map_pointX[i_slboard][ichip][ichn],map_pointY[i_slboard][ichip][ichn],empv);
	  widthMPV_2d->Fill(map_pointX[i_slboard][ichip][ichn],map_pointY[i_slboard][ichip][ichn],wmpv);
	  chi2NDF_2d->Fill(map_pointX[i_slboard][ichip][ichn],map_pointY[i_slboard][ichip][ichn],chi2ndf);
	  mipEntries_2d->Fill(map_pointX[i_slboard][ichip][ichn],map_pointY[i_slboard][ichip][ichn],mipentries);
	  MPV_chip[ichip]->Fill(mpv);
	  MPV_slab->Fill(mpv);
	  
	  MPV_chipchn->Fill(ichn,ichip,mpv);
	  eMPV_chipchn->Fill(ichn,ichip,empv);
	  widthMPV_chipchn->Fill(ichn,ichip,wmpv);
	  mipEntries_chipchn->Fill(ichn,ichip,mipentries);
	}
	/// Signal over NOISE:
	fr[0]=s_n_histo.at(ichip).at(ichn)->GetMean()-0.8*s_n_histo.at(ichip).at(ichn)->GetRMS();
	fr[1]=s_n_histo.at(ichip).at(ichn)->GetMean();//+0.1*s_n_histo.at(ichip).at(ichn)->GetRMS();                                                                                                         
        sv[0]=s_n_histo.at(ichip).at(ichn)->GetRMS()*0.5;
        sv[1]=s_n_histo.at(ichip).at(ichn)->GetMean()*0.6;
        sv[2]=s_n_histo.at(ichip).at(ichn)->Integral("width");
        sv[3]=s_n_histo.at(ichip).at(ichn)->GetRMS()/5.;

        TF1 *fitsnr_temp2=langaufit(s_n_histo.at(ichip).at(ichn),fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
        s_n_histo.at(ichip).at(ichn)->Write();

        double s_n_temp=fitsnr_temp2->GetParameter(1);

	if(ichn<64) {
	  S_N_2d->Fill(map_pointX[i_slboard][ichip][ichn],map_pointY[i_slboard][ichip][ichn],s_n_temp);
	  S_N_chip[ichip]->Fill(s_n_temp);
	  S_N_slab->Fill(s_n_temp);
	

	  //correlations
	  correl_S_N_vs_nhits->Fill(s_n_temp,mipentries);
	  correl_S_N_vs_mpv->Fill(s_n_temp,mpv);
	  correl_S_N_vs_widthmpv->Fill(s_n_temp,wmpv);
	  correl_mpv_vs_nhits->Fill(mpv,mipentries);
	  
	  for(int isca=0; isca<15; isca++) {
	    correl_S_N_vs_widthped->Fill(s_n_temp,ped_width.at(ichip).at(ichn).at(isca));
	    correl_widthped_vs_nhits->Fill(ped_width.at(ichip).at(ichn).at(isca),mipentries);
	  }
	  
	  fout_mip<<ichip<<" "<<ichn<<" "<<mpv<<" "<<empv<<" "<<wmpv<<" "<<chi2ndf<<" "<<mipentries<<endl;
	}	
      } else {
	if(ichn<64) fout_mip<<ichip<<" "<<ichn<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<endl;

      }
    }
  }


  signalfile_summary->cd();

  TCanvas *canvas_maps = new TCanvas("Signal","Signal",1200,800);
  canvas_maps->Divide(2,2);
  TCanvas *canvas_summary_mip = new TCanvas("mip_perchip","mip_perchip",1200,800);
  canvas_summary_mip->Divide(4,4);
  TCanvas *canvas_summary_s_n = new TCanvas("s_n_perchip","s_n_perchip",1200,1200);
  canvas_summary_s_n->Divide(4,4); 
  TCanvas *canvas_summary = new TCanvas("summary","summary",1200,600);
  canvas_summary->Divide(2,1);

  TCanvas *canvas_correlations_S_N = new TCanvas("correlations_S_N","correlations_S_N",800,900);
  canvas_correlations_S_N->Divide(2,3);

  
  // SUMMARY MAPS
  eMPV_2d->SetTitle("eMIP[ADC] map");
  eMPV_2d->GetXaxis()->SetTitle("x");
  eMPV_2d->GetYaxis()->SetTitle("y");
  eMPV_2d->Write();
  widthMPV_2d->SetTitle("widthMIP[ADC] map");
  widthMPV_2d->GetXaxis()->SetTitle("x");
  widthMPV_2d->GetYaxis()->SetTitle("y");
  widthMPV_2d->Write();

  canvas_maps->cd(1);
  MPV_2d->SetStats(kFALSE);
  MPV_2d->SetTitle("MIP[ADC] map");
  MPV_2d->GetXaxis()->SetTitle("x");
  MPV_2d->GetYaxis()->SetTitle("y");
  //MPV_2d->GetZaxis()->SetRangeUser(0,100);
  MPV_2d->Draw("colz");
  MPV_2d->Write();
  canvas_maps->cd(3);
  gPad->SetLogz();
  chi2NDF_2d->SetStats(kFALSE);
  chi2NDF_2d->SetTitle("chi2NDF map");
  chi2NDF_2d->GetXaxis()->SetTitle("x");
  chi2NDF_2d->GetYaxis()->SetTitle("y");
  //chi2NDF_2d->GetZaxis()->SetRangeUser(0,20);
  chi2NDF_2d->Draw("colz");
  chi2NDF_2d->Write();
  canvas_maps->cd(4);
  gPad->SetLogz();
  mipEntries_2d->SetStats(kFALSE);
  mipEntries_2d->SetTitle("Hits map");
  mipEntries_2d->GetXaxis()->SetTitle("x");
  mipEntries_2d->GetYaxis()->SetTitle("y");
  mipEntries_2d->Draw("colz");
  mipEntries_2d->Write();
  canvas_maps->cd(2);
  S_N_2d->SetStats(kFALSE);
  S_N_2d->SetTitle("S / N map");
  S_N_2d->GetXaxis()->SetTitle("x");
  S_N_2d->GetYaxis()->SetTitle("y");
  //S_N_2d->GetZaxis()->SetRangeUser(0,50);
  S_N_2d->Draw("colz");
  S_N_2d->Write();

  MPV_chipchn->SetTitle("MIP[ADC] map");
  MPV_chipchn->GetXaxis()->SetTitle("chip");
  MPV_chipchn->GetYaxis()->SetTitle("chn");
  MPV_chipchn->Write();
  widthMPV_chipchn->SetTitle("widthMIP[ADC] map");
  widthMPV_chipchn->GetXaxis()->SetTitle("chip");
  widthMPV_chipchn->GetYaxis()->SetTitle("chn");
  widthMPV_chipchn->Write();
  eMPV_chipchn->SetTitle("eMIP[ADC] map");
  eMPV_chipchn->GetXaxis()->SetTitle("chip");
  eMPV_chipchn->GetYaxis()->SetTitle("chn");
  eMPV_chipchn->Write();
  mipEntries_chipchn->SetTitle("widthMIP[ADC] map");
  mipEntries_chipchn->GetXaxis()->SetTitle("chip");
  mipEntries_chipchn->GetYaxis()->SetTitle("chn");
  mipEntries_chipchn->Write();


  
  // SUMMARY histograms (per chip)  
  for(int ichip=0; ichip<16; ichip ++) {
    TString subtitle=TString::Format("chip%i",ichip);
    canvas_summary_mip->cd(ichip+1);
    MPV_chip[ichip]->SetStats(kFALSE);
    MPV_chip[ichip]->SetTitle("MIP[ADC]_map, "+subtitle);
    MPV_chip[ichip]->GetXaxis()->SetTitle("x");
    MPV_chip[ichip]->GetYaxis()->SetTitle("y");
    MPV_chip[ichip]->GetZaxis()->SetRangeUser(0,100);
    //noise.at(ichip)->SetLineColor(2);
    MPV_chip[ichip]->Draw("colz");
    MPV_chip[ichip]->Write();
  }

  for(int ichip=0; ichip<16; ichip ++) {
    TString subtitle=TString::Format("chip%i",ichip);
    canvas_summary_s_n->cd(ichip+1);
    S_N_chip[ichip]->SetStats(kFALSE);
    S_N_chip[ichip]->SetTitle("S_N_map, "+subtitle);
    S_N_chip[ichip]->GetXaxis()->SetTitle("x");
    S_N_chip[ichip]->GetYaxis()->SetTitle("y");
    S_N_chip[ichip]->GetZaxis()->SetRangeUser(0,100);
    //noise.at(ichip)->SetLineColor(2);
    S_N_chip[ichip]->Draw("colz");
    S_N_chip[ichip]->Write();
  }

  // SUMMARY histograms (per chip)  
  canvas_summary->cd(1);
  MPV_slab->SetTitle("MIP, "+slboard);
  MPV_slab->GetXaxis()->SetTitle("MIP[ADC] (ped. subtr)");
  MPV_slab->GetYaxis()->SetTitle("fitted channels");
  MPV_slab->Draw("h");
  MPV_slab->Write();
  canvas_summary->cd(2);
  S_N_slab->SetTitle("S / N, "+slboard);
  S_N_slab->GetXaxis()->SetTitle("no units");
  S_N_slab->GetYaxis()->SetTitle("fitted channels");
  S_N_slab->Draw("h");
  S_N_slab->Write();
  canvas_summary->Print("results_mipcalibration/Signal_summary"+slboard+".png");

  //Correlations 
  canvas_correlations_S_N->cd(1);
  correl_S_N_vs_nhits->SetStats(kFALSE);
  correl_S_N_vs_nhits->SetTitle(slboard);
  correl_S_N_vs_nhits->GetXaxis()->SetTitle("S / N");
  correl_S_N_vs_nhits->GetYaxis()->SetTitle("nhits");

  correl_S_N_vs_nhits->Draw("colz");
  correl_S_N_vs_nhits->Write();

  canvas_correlations_S_N->cd(2);
  correl_S_N_vs_mpv->SetStats(kFALSE);
  correl_S_N_vs_mpv->SetTitle(slboard);
  correl_S_N_vs_mpv->GetXaxis()->SetTitle("S / N");
  correl_S_N_vs_mpv->GetYaxis()->SetTitle("MIP[ADC] (ped. subst.)");

  correl_S_N_vs_mpv->Draw("colz");
  correl_S_N_vs_mpv->Write();


  canvas_correlations_S_N->cd(3); 
  correl_S_N_vs_widthmpv->SetStats(kFALSE);
  correl_S_N_vs_widthmpv->SetTitle(slboard);
  correl_S_N_vs_widthmpv->GetXaxis()->SetTitle("S / N");
  correl_S_N_vs_widthmpv->GetYaxis()->SetTitle("MIP width [ADC]");

  correl_S_N_vs_widthmpv->Draw("colz");
  correl_S_N_vs_widthmpv->Write();


  canvas_correlations_S_N->cd(4);
  correl_S_N_vs_widthped->SetStats(kFALSE);
  correl_S_N_vs_widthped->SetTitle(slboard);
  correl_S_N_vs_widthped->GetXaxis()->SetTitle("S / N");
  correl_S_N_vs_widthped->GetYaxis()->SetTitle("Ped. width [ADC]");

  correl_S_N_vs_widthped->Draw("colz");
  correl_S_N_vs_widthped->Write();

  canvas_correlations_S_N->cd(5);
  correl_mpv_vs_nhits->SetStats(kFALSE);
  correl_mpv_vs_nhits->SetTitle(slboard);
  correl_mpv_vs_nhits->GetXaxis()->SetTitle("MIP [ADC] (ped. subs.)");
  correl_mpv_vs_nhits->GetYaxis()->SetTitle("nhits");

  correl_mpv_vs_nhits->Draw("colz");
  correl_mpv_vs_nhits->Write();

  canvas_correlations_S_N->cd(6);
  correl_widthped_vs_nhits->SetStats(kFALSE);
  correl_widthped_vs_nhits->SetTitle(slboard);
  correl_widthped_vs_nhits->GetXaxis()->SetTitle("Ped. width [ADC]");
  correl_widthped_vs_nhits->GetYaxis()->SetTitle("nhits");

  correl_widthped_vs_nhits->Draw("colz");
  correl_widthped_vs_nhits->Write();




  canvas_maps->Write();
  canvas_summary_mip->Write();
  canvas_summary_s_n->Write(); 
  canvas_summary->Write();
  canvas_correlations_S_N->Write();

  delete canvas_maps;
  delete canvas_summary_mip;
  delete canvas_summary_s_n;
  delete canvas_summary;
  delete canvas_correlations_S_N;

  signalfile_summary->Close();

}


void DecodedSLBAnalysis::ReadPedestals(TString filename) 
{

  std::ifstream reading_file(filename);
  if(!reading_file){
    cout<<" ERROR  ----------------------------------- No pedestal file: "<<filename<<endl;
  } else {
        cout<<" Pedestal input file: "<<filename<<endl;
  }

  for(int i=0; i<16; i++) {
    std::vector<std::vector<Double_t> > chip_ped_mean;
    std::vector<std::vector<Double_t> > chip_ped_error;
    std::vector<std::vector<Double_t> > chip_ped_width;

    for(int j=0; j<64; j++) {
      std::vector<Double_t> chn_ped_mean;
      std::vector<Double_t> chn_ped_error;
      std::vector<Double_t> chn_ped_width;
      chip_ped_mean.push_back(chn_ped_mean);
      chip_ped_error.push_back(chn_ped_error);
      chip_ped_width.push_back(chn_ped_width);
    }  
    ped_mean.push_back(chip_ped_mean);
    ped_error.push_back(chip_ped_error);
    ped_width.push_back(chip_ped_width);
  }

  for(int i=0; i<16; i++) {
    for(int j=0; j<64; j++) {
      for(int isca=0; isca<15; isca++) {
	ped_mean.at(i).at(j).push_back(0.);
	ped_error.at(i).at(j).push_back(0.);
	ped_width.at(i).at(j).push_back(0.);
      }
    }
  }

  Int_t tmp_chip = 0,tmp_channel = 0;
  Double_t tmp_ped[15], tmp_error[15], tmp_width[15];
  for(int isca=0; isca<15; isca++) {
    tmp_ped[isca]=0.;
    tmp_error[isca]=0.;
    tmp_width[isca]=0.;
  }
  
  TString tmpst;
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >>  tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst  >> tmpst >> tmpst ;
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >>  tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst ;

  while(reading_file){
    reading_file >> tmp_chip >> tmp_channel >> tmp_ped[0] >> tmp_error[0] >> tmp_width[0] >> tmp_ped[1] >> tmp_error[1] >> tmp_width[1] >> tmp_ped[2] >> tmp_error[2] >> tmp_width[2] >> tmp_ped[3] >> tmp_error[3] >> tmp_width[3] >> tmp_ped[4] >> tmp_error[4] >> tmp_width[4] >> tmp_ped[5] >> tmp_error[5] >> tmp_width[5] >> tmp_ped[6] >> tmp_error[6] >> tmp_width[6] >> tmp_ped[7] >> tmp_error[7] >> tmp_width[7] >> tmp_ped[8] >> tmp_error[8] >> tmp_width[8] >> tmp_ped[9] >> tmp_error[9] >> tmp_width[9] >> tmp_ped[10] >> tmp_error[10] >> tmp_width[10] >> tmp_ped[11] >> tmp_error[11] >> tmp_width[11] >> tmp_ped[12] >> tmp_error[12] >> tmp_width[12] >> tmp_ped[13] >> tmp_error[13] >> tmp_width[13] >> tmp_ped[14] >> tmp_error[14] >> tmp_width[14];
    //    cout<<tmp_chip <<" "<< tmp_channel << " "<< tmp_ped[0]<<endl;

    for(int isca=0; isca<15; isca++) {
      if(tmp_ped[isca]>0. ){//&& (tmp_error[isca]<ped_error.at(tmp_chip).at(tmp_channel).at(isca) || ped_error.at(tmp_chip).at(tmp_channel).at(isca)==0) ){
	ped_mean.at(tmp_chip).at(tmp_channel).at(isca)=tmp_ped[isca];
	ped_error.at(tmp_chip).at(tmp_channel).at(isca)=tmp_error[isca];
	ped_width.at(tmp_chip).at(tmp_channel).at(isca)=tmp_width[isca];
	}
    }
    
  }

  /*
  TString filename2=filename+"new";
  ofstream fout_ped(filename2,ios::out);
  fout_ped<<"#pedestal results (fit to a gaussian) remove channels/sca with two pedestals peaks from the analysis "<<filename2<<endl;
  fout_ped<<"#chip channel ped0 eped0 widthped0 ped1 eped1 widthped1... ped14 eped14 widhtped14 (all SCA)"<<endl;
  for(int ichip=0; ichip<16; ichip++) {
    for(int ichn=0; ichn<64; ichn++) {
      fout_ped << ichip <<" " <<ichn<< " ";
      cout << ichip <<" " <<ichn<< " ";
      for(int isca=0; isca<15; isca++) {
    	fout_ped <<ped_mean.at(ichip).at(ichn).at(isca)<< " "<< " "<<ped_error.at(ichip).at(ichn).at(isca)<< " "<<ped_width.at(ichip).at(ichn).at(isca)<<" ";
	cout <<ped_mean.at(ichip).at(ichn).at(isca)<< " "<< " "<<ped_error.at(ichip).at(ichn).at(isca)<< " "<<ped_width.at(ichip).at(ichn).at(isca)<<" ";

      }
      fout_ped<<endl;
      cout<<endl;
    }
    }*/
  

}

void DecodedSLBAnalysis::PedestalAnalysis(int i_slboard,TString outputname="", int maxnhit=5, int ncoinc=7)
{

  //Read the channel/chip -- x/y mapping
  //  ReadMap(map_filename);

  TString slboard = TString::Format("_slboard_%i_",slboard_array_mapping[i_slboard]);
  if(outputname!="") slboard=slboard+outputname;
  
  //Read the list of masked channels
  ReadMasked("maskedchannels/masked"+slboard+".txt");
  
  ofstream fout_ped("results_pedestal/Pedestal"+slboard+".txt",ios::out);

  fout_ped<<"#pedestal results (fit to a gaussian) remove channels/sca with two pedestals peaks from the analysis : "<<slboard<<endl;
  fout_ped<<"#chip channel ped0 eped0 widthped0 ped1 eped1 widthped1... ped14 eped14 widhtped14 (all SCA)"<<endl;

  
  bool global = true;
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //  nentries=10000;
  TH2F* pedestal_map[15];
  TH2F* pedestal_width_map[15];
  TH2F* pedestal_error_map[15];
  TH2F* pedestal_chi2ndf_map[15];
  TH2F* pedestal_npeaks_map[15];
  TH2F* pedestal_entries_map[15];

  for(int isca=0; isca<15; isca++) {
    pedestal_map[isca]= new TH2F(TString::Format("pedestal_map_sca%i",isca),TString::Format("pedestal_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_width_map[isca]= new TH2F(TString::Format("pedestal_width_map_sca%i",isca),TString::Format("pedestal_width_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_error_map[isca]= new TH2F(TString::Format("pedestal_error_map_sca%i",isca),TString::Format("pedestal_error_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_chi2ndf_map[isca]= new TH2F(TString::Format("pedestal_chi2ndf_map_sca%i",isca),TString::Format("pedestal_chi2ndf_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_npeaks_map[isca]= new TH2F(TString::Format("pedestal_npeaks_map_sca%i",isca),TString::Format("pedestal_npeaks_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_entries_map[isca]= new TH2F(TString::Format("pedestal_entries_map_sca%i",isca),TString::Format("pedestal_entries_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
  }

  TH2F* pedestal_tagged_map[15];
  TH2F* pedestal_tagged_width_map[15];
  TH2F* pedestal_tagged_error_map[15];
  TH2F* pedestal_tagged_chi2ndf_map[15];
  TH2F* pedestal_tagged_npeaks_map[15];
  TH2F* pedestal_tagged_entries_map[15];

  for(int isca=0; isca<15; isca++) {
    pedestal_tagged_map[isca]= new TH2F(TString::Format("pedestal_tagged_map_sca%i",isca),TString::Format("pedestal_tagged_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_tagged_width_map[isca]= new TH2F(TString::Format("pedestal_tagged_width_map_sca%i",isca),TString::Format("pedestal_tagged_width_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_tagged_error_map[isca]= new TH2F(TString::Format("pedestal_tagged_error_map_sca%i",isca),TString::Format("pedestal_tagged_error_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_tagged_chi2ndf_map[isca]= new TH2F(TString::Format("pedestal_tagged_chi2ndf_map_sca%i",isca),TString::Format("pedestal_tagged_chi2ndf_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_tagged_npeaks_map[isca]= new TH2F(TString::Format("pedestal_tagged_npeaks_map_sca%i",isca),TString::Format("pedestal_tagged_npeaks_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_tagged_entries_map[isca]= new TH2F(TString::Format("pedestal_tagged_entries_map_sca%i",isca),TString::Format("pedestal_tagged_entries_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
  }
 
  // --------------------
  // all sca
  //  std::vector<TCanvas*> chip;
  std::vector<std::vector<std::vector<TH1F*> > > ped_sca;
  std::vector<std::vector<std::vector<TH1F*> > > ped_sca_tagged;

  std::vector<TH1F*> pedestal_chip ;
  std::vector<TH1F*> pedestal_tagged_chip ;
  std::vector<TH1F*> pedestal_slboard_chip ;
  std::vector<TH1F*> pedestal_tagged_slboard_chip ;


  for(int ichip=0; ichip<16; ichip++) {
    TH1F *ped_chip = new TH1F(TString::Format("ped_chip%i",ichip),TString::Format("ped_chip%i",ichip),1000,0.5,1000.5);
    pedestal_chip.push_back(ped_chip);

    TH1F *ped_tagged_chip = new TH1F(TString::Format("ped_tagged_chip%i",ichip),TString::Format("ped_tagged_chip%i",ichip),1000,0.5,1000.5);
    pedestal_tagged_chip.push_back(ped_tagged_chip);

    TH1F *ped_slboard_chip = new TH1F(TString::Format("ped_slboard_chip%i",ichip),TString::Format("ped_slboard_chip%i",ichip),1002,-500,500);
    pedestal_slboard_chip.push_back(ped_slboard_chip);
    
    TH1F *ped_tagged_slboard_chip = new TH1F(TString::Format("ped_tagged_slboard_chip%i",ichip),TString::Format("ped_tagged_slboard_chip%i",ichip),1002,-500,500);
    pedestal_tagged_slboard_chip.push_back(ped_tagged_slboard_chip);

    std::vector<std::vector<TH1F*> >pedtemp_sca;
    std::vector<std::vector<TH1F*> >pedtemp_sca_tagged;

    for(int ichn=0; ichn<64; ichn++) {
      std::vector<TH1F*> pedtemp_sca2;
      std::vector<TH1F*> pedtemp_sca2_tagged;

      for(int isca=0; isca<15; isca++) {
	TH1F *ped_sca2 = new TH1F(TString::Format("ped_chip%i_chn%i_sca%i",ichip,ichn,isca),TString::Format("ped_chip%i_chn%i_sca%i",ichip,ichn,isca),1000,0.5,1000.5);
	TH1F *ped_sca2_tagged = new TH1F(TString::Format("ped_tagged_chip%i_chn%i_sca%i",ichip,ichn,isca),TString::Format("ped_tagged_chip%i_chn%i_sca%i",ichip,ichn,isca),1000,0.5,1000.5);
	pedtemp_sca2.push_back(ped_sca2);
	pedtemp_sca2_tagged.push_back(ped_sca2_tagged);
      }
      pedtemp_sca.push_back(pedtemp_sca2);
      pedtemp_sca_tagged.push_back(pedtemp_sca2_tagged);
    }
    ped_sca.push_back(pedtemp_sca);
    ped_sca_tagged.push_back(pedtemp_sca_tagged);
  }


  // -----------------------------------------------------------------------------------------------------   
  // SCA analysis
  Long64_t nbytes = 0, nb = 0;
  //  nentries=1000;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if ( jentry > 1000 && jentry % 1000 ==0 ) std::cout << "Progress: " << 100.*jentry/nentries <<" %"<<endl;


    std::array<int,4096> bcid_seen = SimpleCoincidenceTagger(jentry,maxnhit);

    
    for(int islboard=0; islboard<n_slboards; islboard++) {
      //if(islboard != i_slboard) continue;
    
      for(int ichip=0; ichip<16; ichip++) {

	for(int isca=0; isca<15; isca++) {
	  if(bcid[islboard][ichip][isca]<0) break;

	  bool gooddata=true;
	  if( bcid_seen.at(bcid[islboard][ichip][isca]) < ncoinc ) gooddata=false; //only analyze events with at least ncoinc slabs with hit in such bcid
	  for(int ichn=0; ichn<64; ichn++) if(charge_hiGain[islboard][ichip][isca][ichn]<150 && charge_hiGain[islboard][ichip][isca][ichn]>-1 ) gooddata=false;
	  /*bool gooddata=true;
	  if(global == true) {
	    int ntaggedasbad = 0;
	    for(int ichn=0; ichn<64; ichn++) {
	      if(gain_hit_high[islboard][ichip][isca][ichn]==1 && charge_hiGain[islboard][ichip][isca][ichn]<15 && charge_hiGain[islboard][ichip][isca][ichn]>-1 ) 
		ntaggedasbad++;
	    }//ichn 
	    if ( ntaggedasbad > 0) gooddata=false;
	    }*/
	
	  for(int ichn=0; ichn<64; ichn++) {

	    //good events
	    bool selection=false;
	    if( charge_hiGain[islboard][ichip][isca][ichn]>30 && badbcid[islboard][ichip][isca]==0 && nhits[islboard][ichip][isca]<maxnhit+1 && bcid[islboard][ichip][isca]>50) selection=true;
		  
	    // if(masked[islboard][ichip][ichn]==1) selection=false;
	    if(gain_hit_high[islboard][ichip][isca][ichn]==0 && gooddata==true && selection==true)
	      ped_sca.at(ichip).at(ichn).at(isca)->Fill(charge_hiGain[islboard][ichip][isca][ichn]);

	    //bad events
	    //	  selection=false;
	    //	  if( ( badbcid[islboard][ichip][isca]!=0  || nhits[islboard][ichip][isca]>maxnhit || gooddata==false) ) selection=true;
	    //	  if(masked[islboard][ichip][ichn]==1) selection=false;
	    if(gain_hit_high[islboard][ichip][isca][ichn]==0 && (selection==false || gooddata==false) )
	      ped_sca_tagged.at(ichip).at(ichn).at(isca)->Fill(charge_hiGain[islboard][ichip][isca][ichn]);

	  }
           
	}//isca

      }//ichip 

    }// islboard 
  }  // end first loop analysis to fill pedestal historgrams


  TFile *pedfile = new TFile("results_pedestal/Pedestal"+slboard+".root" , "RECREATE");
  pedfile->cd();

  //initialize pedestal vectors
  for(int i=0; i<16; i++) {
    std::vector<std::vector<Double_t> > chip_ped_mean;
    std::vector<std::vector<Double_t> > chip_ped_error;
    std::vector<std::vector<Double_t> > chip_ped_width;

    for(int j=0; j<64; j++) {
      std::vector<Double_t> chn_ped_mean;
      std::vector<Double_t> chn_ped_error;
      std::vector<Double_t> chn_ped_width;
      chip_ped_mean.push_back(chn_ped_mean);
      chip_ped_error.push_back(chn_ped_error);
      chip_ped_width.push_back(chn_ped_width);
    }  
    ped_mean.push_back(chip_ped_mean);
    ped_error.push_back(chip_ped_error);
    ped_width.push_back(chip_ped_width);
  }

  // do pedestal (chip/channel/sca based) analysis
  for(int ichip=0; ichip<16; ichip++) {
    for(int ichn=0; ichn<64; ichn++) {
      
      fout_ped << ichip <<" " <<ichn<< " "; 
      for(int isca=0; isca<15; isca++) {

	ped_sca.at(ichip).at(ichn).at(isca)->SetTitle(TString::Format("ped_chip%i_chn%i_sca%i",ichip,ichn,isca));
	ped_sca.at(ichip).at(ichn).at(isca)->SetName(TString::Format("ped_chip%i_chn%i_sca%i",ichip,ichn,isca));
	ped_sca.at(ichip).at(ichn).at(isca)->Write();

	ped_sca_tagged.at(ichip).at(ichn).at(isca)->SetTitle(TString::Format("ped_tagged_chip%i_chn%i_sca%i",ichip,ichn,isca));
	ped_sca_tagged.at(ichip).at(ichn).at(isca)->SetName(TString::Format("ped_tagged_chip%i_chn%i_sca%i",ichip,ichn,isca));
	ped_sca_tagged.at(ichip).at(ichn).at(isca)->Write();

	pedestal_entries_map[isca] -> Fill( map_pointX[i_slboard][ichip][ichn] , map_pointY[i_slboard][ichip][ichn] , ped_sca.at(ichip).at(ichn).at(isca)->GetEntries());
	pedestal_tagged_entries_map[isca] -> Fill( map_pointX[i_slboard][ichip][ichn] , map_pointY[i_slboard][ichip][ichn] , ped_sca_tagged.at(ichip).at(ichn).at(isca)->GetEntries());

	
	ped_mean.at(ichip).at(ichn).push_back(0.);
	ped_error.at(ichip).at(ichn).push_back(0.);
	ped_width.at(ichip).at(ichn).push_back(0.);

	if(ped_sca.at(ichip).at(ichn).at(isca)->GetEntries()> 50 ){ //max_entries/2 ) {
	  TSpectrum *s = new TSpectrum();
	  int npeaks = s->Search(ped_sca.at(ichip).at(ichn).at(isca),2,"",0.8); 
	  pedestal_npeaks_map[isca] -> Fill( map_pointX[i_slboard][ichip][ichn] , map_pointY[i_slboard][ichip][ichn] , npeaks);
	  float mean_for_fit=0;
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
            for(int ipeak=0; ipeak<npeaks; ipeak ++) {
	      if(ipeak != npeak_max) pedestal_slboard_chip.at(ichip)->Fill( mean_peak[npeak_max] - mean_peak[ipeak]);
	    }
	    //	    if(npeaks ==1 ) {
	    //Double_t *mean_peak=s->GetPositionX();
	    mean_peak[0]=mean_peak_higher;
	    mean_for_fit=mean_peak[0];
	  } else mean_for_fit=ped_sca.at(ichip).at(ichn).at(isca)->GetMean();
	      
	  TF1 *f0 = new TF1("f0","gaus",mean_for_fit-ped_sca.at(ichip).at(ichn).at(isca)->GetRMS()/2.,mean_for_fit+ped_sca.at(ichip).at(ichn).at(isca)->GetRMS()/2.);
	  ped_sca.at(ichip).at(ichn).at(isca)->Fit("f0","RQNOC");
	  
	  TF1 *f1 = new TF1("f1","gaus",f0->GetParameter(1)-f0->GetParameter(2)/2.,f0->GetParameter(1)+f0->GetParameter(2)/2.);
	  ped_sca.at(ichip).at(ichn).at(isca)->Fit("f1","RQME");
	  fout_ped<<f1->GetParameter(1) << " " << f1->GetParError(1)<<" "<<f1->GetParameter(2)<< " ";
	  
	  ped_mean.at(ichip).at(ichn).at(isca)=f1->GetParameter(1);
	  ped_error.at(ichip).at(ichn).at(isca)=f1->GetParError(1);
	  ped_width.at(ichip).at(ichn).at(isca)=f1->GetParameter(2);
	  pedestal_chip.at(ichip)->Fill(f1->GetParameter(1));
	  pedestal_map[isca] -> Fill( map_pointX[i_slboard][ichip][ichn] , map_pointY[i_slboard][ichip][ichn] , f1->GetParameter(1));
	  pedestal_width_map[isca] -> Fill( map_pointX[i_slboard][ichip][ichn] , map_pointY[i_slboard][ichip][ichn] , f1->GetParameter(2));
	  pedestal_error_map[isca] -> Fill( map_pointX[i_slboard][ichip][ichn] , map_pointY[i_slboard][ichip][ichn] , f1->GetParError(1));
	  pedestal_chi2ndf_map[isca] -> Fill( map_pointX[i_slboard][ichip][ichn] , map_pointY[i_slboard][ichip][ichn] , f1->GetChisquare() / f1->GetNDF());
	    
	    /*} else {
	      fout_ped<<ped_sca.at(ichip).at(ichn).at(isca)->GetMean()<< " " << ped_sca.at(ichip).at(ichn).at(isca)->GetRMS()/sqrt(ped_sca.at(ichip).at(ichn).at(isca)->GetEntries())<<" "<<ped_sca.at(ichip).at(ichn).at(isca)->GetRMS()<<" ";
	      pedestal_map[isca] -> Fill( map_pointX[i_slboard][ichip][ichn] , map_pointY[i_slboard][ichip][ichn] , ped_sca.at(ichip).at(ichn).at(isca)->GetMean());
	      pedestal_width_map[isca] -> Fill( map_pointX[i_slboard][ichip][ichn] , map_pointY[i_slboard][ichip][ichn] , ped_sca.at(ichip).at(ichn).at(isca)->GetRMS() );
	      pedestal_error_map[isca] -> Fill( map_pointX[i_slboard][ichip][ichn] , map_pointY[i_slboard][ichip][ichn] , ped_sca.at(ichip).at(ichn).at(isca)->GetRMS()/sqrt(ped_sca.at(ichip).at(ichn).at(isca)->GetMean()));
	      pedestal_chi2ndf_map[isca] -> Fill( map_pointX[i_slboard][ichip][ichn] , map_pointY[i_slboard][ichip][ichn] , 0);
	      }*/
	    
	} else {
	  fout_ped<<0<< " " << 0<<" "<<0<<" "; 
	} 
     
	// analyze pedestal for tagged events
	if(ped_sca_tagged.at(ichip).at(ichn).at(isca)->GetEntries()> 250) {
	  TSpectrum *s = new TSpectrum();
	  int npeaks = s->Search(ped_sca_tagged.at(ichip).at(ichn).at(isca),2,"",0.2);
	  pedestal_tagged_npeaks_map[isca] -> Fill( map_pointX[i_slboard][ichip][ichn] , map_pointY[i_slboard][ichip][ichn] , npeaks);
	  
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
	    
	    for(int ipeak=0; ipeak<npeaks; ipeak ++) {
              if(ipeak != npeak_max) pedestal_tagged_slboard_chip.at(ichip)->Fill( mean_peak[npeak_max] - mean_peak[ipeak] );
            }
	    
	    if(mean_peak_higher>0) {
	      TF1 *f0 = new TF1("f0","gaus",mean_peak_higher-10.,mean_peak_higher+10.);
	      ped_sca_tagged.at(ichip).at(ichn).at(isca)->Fit("f0","RQNOC");
	      if(f0->GetParameter(1)>0) {
		double xmin=100;
		double xmax=500;
		if(f0->GetParameter(1)-2.*f0->GetParameter(2) > 0 ) xmin=f0->GetParameter(1)-2.*f0->GetParameter(2);
		if(f0->GetParameter(1)+2.*f0->GetParameter(2) < 500 ) xmax=f0->GetParameter(1)+2.*f0->GetParameter(2);
		
		TF1 *f1 = new TF1("f1","gaus",xmin,xmax);
		ped_sca_tagged.at(ichip).at(ichn).at(isca)->Fit("f1","RQ");
		pedestal_tagged_chip.at(ichip)->Fill(f1->GetParameter(1));
		
		pedestal_tagged_map[isca] -> Fill( map_pointX[i_slboard][ichip][ichn] , map_pointY[i_slboard][ichip][ichn] , f1->GetParameter(1));
		pedestal_tagged_width_map[isca] -> Fill( map_pointX[i_slboard][ichip][ichn] , map_pointY[i_slboard][ichip][ichn] , f1->GetParameter(2));
		pedestal_tagged_error_map[isca] -> Fill( map_pointX[i_slboard][ichip][ichn] , map_pointY[i_slboard][ichip][ichn] , f1->GetParError(1));
		pedestal_tagged_chi2ndf_map[isca] -> Fill( map_pointX[i_slboard][ichip][ichn] , map_pointY[i_slboard][ichip][ichn] , f1->GetChisquare() / f1->GetNDF());
	      }
	    }
	  }
	}


      }//isca
      fout_ped<<endl;
    }//ichn
  }//ichip

  pedfile->Close();

  TFile *pedfile_summary = new TFile("results_pedestal/Pedestal_summary"+slboard+".root" , "RECREATE");
  pedfile_summary->cd();
  
  // good pedestal events (not tagged events)
  TCanvas *canvas_pedestal_map = new TCanvas("pedestal_map", "pedestal_map",1200,1200);
  canvas_pedestal_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString subtitle=TString::Format("sca%i",isca);
    canvas_pedestal_map->cd(isca+1);
    pedestal_map[isca]->SetStats(kFALSE);
    pedestal_map[isca]->SetTitle("pedestal_map, "+subtitle);
    pedestal_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_map[isca]->GetZaxis()->SetRangeUser(200,400);

    //noise.at(ichip)->SetLineColor(2);
    pedestal_map[isca]->Draw("colz");
    pedestal_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_width_map = new TCanvas("pedestal_width_map", "pedestal_width_map",1200,1200);
  canvas_pedestal_width_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString subtitle=TString::Format("sca%i",isca);
    canvas_pedestal_width_map->cd(isca+1);
    pedestal_width_map[isca]->SetStats(kFALSE);
    pedestal_width_map[isca]->SetTitle("pedestal_width_map, "+subtitle);
    pedestal_width_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_width_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_width_map[isca]->GetZaxis()->SetRangeUser(1,20);
    //noise.at(ichip)->SetLineColor(2);
    pedestal_width_map[isca]->Draw("colz");
    pedestal_width_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_error_map = new TCanvas("pedestal_error_map", "pedestal_error_map",1200,1200);
  canvas_pedestal_error_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString subtitle=TString::Format("sca%i",isca);
    canvas_pedestal_error_map->cd(isca+1);
    pedestal_error_map[isca]->SetStats(kFALSE);
    pedestal_error_map[isca]->SetTitle("pedestal_error_map, "+subtitle);
    pedestal_error_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_error_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_error_map[isca]->Draw("colz");
    pedestal_error_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_npeaks_map = new TCanvas("pedestal_npeaks_map", "pedestal_npeaks_map",1200,1200);
  canvas_pedestal_npeaks_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString subtitle=TString::Format("sca%i",isca);
    canvas_pedestal_npeaks_map->cd(isca+1);
    pedestal_npeaks_map[isca]->SetStats(kFALSE);
    pedestal_npeaks_map[isca]->SetTitle("pedestal_npeaks_map, "+subtitle);
    pedestal_npeaks_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_npeaks_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_npeaks_map[isca]->GetZaxis()->SetRangeUser(0.5,4.5);
    //noise.at(ichip)->SetLineColor(2);
    pedestal_npeaks_map[isca]->Draw("colz");
    pedestal_npeaks_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_chi2ndf_map = new TCanvas("pedestal_chi2ndf_map", "pedestal_chi2ndf_map",1200,1200);
  canvas_pedestal_chi2ndf_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString subtitle=TString::Format("sca%i",isca);
    canvas_pedestal_chi2ndf_map->cd(isca+1);
    pedestal_chi2ndf_map[isca]->SetStats(kFALSE);
    pedestal_chi2ndf_map[isca]->SetTitle("pedestal_chi2ndf_map, "+subtitle);
    pedestal_chi2ndf_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_chi2ndf_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_chi2ndf_map[isca]->Draw("colz");
    pedestal_chi2ndf_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_entries_map = new TCanvas("pedestal_entries_map", "pedestal_entries_map",1200,1200);
  canvas_pedestal_entries_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString subtitle=TString::Format("_sca%i",isca);
    canvas_pedestal_entries_map->cd(isca+1);
    pedestal_entries_map[isca]->SetStats(kFALSE);
    pedestal_entries_map[isca]->SetTitle("pedestal_entries_map, "+subtitle);
    pedestal_entries_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_entries_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_entries_map[isca]->Draw("colz");
    pedestal_entries_map[isca]->Write();
  }
  

  canvas_pedestal_map->Write();
  canvas_pedestal_width_map->Write();
  canvas_pedestal_error_map->Write();
  canvas_pedestal_npeaks_map->Write();
  canvas_pedestal_entries_map->Write();
  canvas_pedestal_chi2ndf_map->Write();

   
  TCanvas *canvas_pedestal_pedestal = new TCanvas("pedestal_average", "pedestal_average",1200,1200);
  canvas_pedestal_pedestal->Divide(4,4);
  for(int ichip=0; ichip<16; ichip++) {
    canvas_pedestal_pedestal->cd(ichip+1);
    //gPad->SetLogy();
    pedestal_chip.at(ichip)->GetXaxis()->SetRangeUser(0,500);
    pedestal_chip.at(ichip)->SetTitle(TString::Format("Average pedestal, chip-%i",ichip));
    pedestal_chip.at(ichip)->GetXaxis()->SetTitle("ADC");
    pedestal_chip.at(ichip)->GetYaxis()->SetTitle("#");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_chip.at(ichip)->Draw("hs");
    //pedestal_chip.at(ichip)->Write();
  }

  canvas_pedestal_pedestal->Write();

  TCanvas *canvas_pedestal_slboard = new TCanvas("pedestal_slboard", "pedestal_slboard",1200,1200);
  canvas_pedestal_slboard->Divide(4,4);
  for(int ichip=0; ichip<16; ichip++) {
    canvas_pedestal_slboard->cd(ichip+1);

    pedestal_slboard_chip.at(ichip)->GetXaxis()->SetRangeUser(-100,100);
    pedestal_slboard_chip.at(ichip)->SetTitle(TString::Format("Pedestal slboard, chip-%i",ichip));
    pedestal_slboard_chip.at(ichip)->GetXaxis()->SetTitle("ADC");
    pedestal_slboard_chip.at(ichip)->GetYaxis()->SetTitle("#");
    pedestal_slboard_chip.at(ichip)->Draw("hs");
  }

  canvas_pedestal_slboard->Write();


  // Tagged events
  TCanvas *canvas_pedestal_tagged_map = new TCanvas("pedestal_tagged_map", "pedestal_tagged_map",1200,1200);
  canvas_pedestal_tagged_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString subtitle=slboard+TString::Format("_sca%i",isca);
    canvas_pedestal_tagged_map->cd(isca+1);
    pedestal_tagged_map[isca]->SetStats(kFALSE);
    pedestal_tagged_map[isca]->SetTitle("pedestal_tagged_map, "+subtitle);
    pedestal_tagged_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_tagged_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_tagged_map[isca]->GetZaxis()->SetRangeUser(200,400);

    //noise.at(ichip)->SetLineColor(2);
    pedestal_tagged_map[isca]->Draw("colz");
    pedestal_tagged_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_tagged_width_map = new TCanvas("pedestal_tagged_width_map", "pedestal_tagged_width_map",1200,1200);
  canvas_pedestal_tagged_width_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString subtitle=TString::Format("sca%i",isca);
    canvas_pedestal_tagged_width_map->cd(isca+1);
    pedestal_tagged_width_map[isca]->SetStats(kFALSE);
    pedestal_tagged_width_map[isca]->SetTitle("pedestal_tagged_width_map, "+subtitle);
    pedestal_tagged_width_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_tagged_width_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_tagged_width_map[isca]->GetZaxis()->SetRangeUser(0,7);
    //noise.at(ichip)->SetLineColor(2);
    pedestal_tagged_width_map[isca]->Draw("colz");
    pedestal_tagged_width_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_tagged_error_map = new TCanvas("pedestal_tagged_error_map", "pedestal_tagged_error_map",1200,1200);
  canvas_pedestal_tagged_error_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString subtitle=TString::Format("sca%i",isca);
    canvas_pedestal_tagged_error_map->cd(isca+1);
    pedestal_tagged_error_map[isca]->SetStats(kFALSE);
    pedestal_tagged_error_map[isca]->SetTitle("pedestal_tagged_error_map, "+subtitle);
    pedestal_tagged_error_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_tagged_error_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_tagged_error_map[isca]->Draw("colz");
    pedestal_tagged_error_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_tagged_npeaks_map = new TCanvas("pedestal_tagged_npeaks_map", "pedestal_tagged_npeaks_map",1200,1200);
  canvas_pedestal_tagged_npeaks_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString subtitle=TString::Format("sca%i",isca);
    canvas_pedestal_tagged_npeaks_map->cd(isca+1);
    pedestal_tagged_npeaks_map[isca]->SetStats(kFALSE);
    pedestal_tagged_npeaks_map[isca]->SetTitle("pedestal_tagged_npeaks_map, "+subtitle);
    pedestal_tagged_npeaks_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_tagged_npeaks_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_tagged_npeaks_map[isca]->GetZaxis()->SetRangeUser(0.5,4.5);
    //noise.at(ichip)->SetLineColor(2);
    pedestal_tagged_npeaks_map[isca]->Draw("colz");
    pedestal_tagged_npeaks_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_tagged_chi2ndf_map = new TCanvas("pedestal_tagged_chi2ndf_map", "pedestal_tagged_chi2ndf_map",1200,1200);
  canvas_pedestal_tagged_chi2ndf_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString subtitle=TString::Format("sca%i",isca);
    canvas_pedestal_tagged_chi2ndf_map->cd(isca+1);
    pedestal_tagged_chi2ndf_map[isca]->SetStats(kFALSE);
    pedestal_tagged_chi2ndf_map[isca]->SetTitle("pedestal_tagged_chi2ndf_map, "+subtitle);
    pedestal_tagged_chi2ndf_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_tagged_chi2ndf_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_tagged_chi2ndf_map[isca]->Draw("colz");
    pedestal_tagged_chi2ndf_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_tagged_entries_map = new TCanvas("pedestal_tagged_entries_map", "pedestal_tagged_entries_map",1200,1200);
  canvas_pedestal_tagged_entries_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString subtitle=TString::Format("sca%i",isca);
    canvas_pedestal_tagged_entries_map->cd(isca+1);
    pedestal_tagged_entries_map[isca]->SetStats(kFALSE);
    pedestal_tagged_entries_map[isca]->SetTitle("pedestal_tagged_entries_map, "+subtitle);
    pedestal_tagged_entries_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_tagged_entries_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_tagged_entries_map[isca]->Draw("colz");
    pedestal_tagged_entries_map[isca]->Write();
  }
  

  canvas_pedestal_tagged_map->Write();
  canvas_pedestal_tagged_width_map->Write();
  canvas_pedestal_tagged_error_map->Write();
  canvas_pedestal_tagged_npeaks_map->Write();
  canvas_pedestal_tagged_entries_map->Write();
  canvas_pedestal_tagged_chi2ndf_map->Write();

   
  TCanvas *canvas_pedestal_tagged_pedestal = new TCanvas("pedestal_tagged_average", "pedestal_tagged_average",1200,1200);
  canvas_pedestal_tagged_pedestal->Divide(4,4);
  for(int ichip=0; ichip<16; ichip++) {
    canvas_pedestal_tagged_pedestal->cd(ichip+1);
    pedestal_tagged_chip.at(ichip)->GetXaxis()->SetRangeUser(0,500);
    pedestal_tagged_chip.at(ichip)->SetTitle(TString::Format("Average pedestal_tagged, chip-%i",ichip));
    pedestal_tagged_chip.at(ichip)->GetXaxis()->SetTitle("ADC");
    pedestal_tagged_chip.at(ichip)->GetYaxis()->SetTitle("#");
    pedestal_tagged_chip.at(ichip)->Draw("hs");
  }

  canvas_pedestal_tagged_pedestal->Write();


  TCanvas *canvas_pedestal_tagged_slboard = new TCanvas("pedestal_tagged_slboard","pedestal_tagged_slboard",1200,1200);
  canvas_pedestal_tagged_slboard->Divide(4,4);
  for(int ichip=0; ichip<16; ichip++) {
    canvas_pedestal_tagged_slboard->cd(ichip+1);

    pedestal_tagged_slboard_chip.at(ichip)->GetXaxis()->SetRangeUser(-100,100);
    pedestal_tagged_slboard_chip.at(ichip)->SetTitle(TString::Format("Pedestal slboard, chip-%i",ichip));
    pedestal_tagged_slboard_chip.at(ichip)->GetXaxis()->SetTitle("ADC");
    pedestal_tagged_slboard_chip.at(ichip)->GetYaxis()->SetTitle("#");
    pedestal_tagged_slboard_chip.at(ichip)->Draw("hs");
  }

  canvas_pedestal_tagged_slboard->Write();

  pedfile_summary->Close();



}
 
void DecodedSLBAnalysis::Retriggers(int i_slboard, TString outputname="",int maxnhit=10)
{

  //function that reads a root file and check which channels have or not signal (hit bit ==1).
  //should be used for root files that contain a full position scan, in order to provide meaninful list of masked channels.
 
  TString slboard = TString::Format("_slboard_%i_",slboard_array_mapping[i_slboard]);
  if(outputname!="") slboard=slboard+outputname;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();

  //summary
  TH2F* h_total_planeE = new TH2F("total_planeE","total_planeE",16,-0.5,15.5,15,-0.5,14.5);
  TH2F* h_total_negE = new TH2F("total_negE","total_negE",16,-0.5,15.5,15,-0.5,14.5);
  TH2F* h_total_negE_nohit = new TH2F("total_negE_nohit","total_negE_nohit",16,-0.5,15.5,15,-0.5,14.5);
  TH2F* h_total_burst = new TH2F("total_burst","total_burst",16,-0.5,15.5,15,-0.5,14.5);
  TH2F* h_total_retrig = new TH2F("total_retrig","total_retrig",16,-0.5,15.5,15,-0.5,14.5);
  TH2F* h_total_retrig_trains = new TH2F("total_retrig_trains","total_retrig_trains",16,-0.5,15.5,15,-0.5,14.5);
  TH2F* h_total_trig = new TH2F("total_trig","total_trig",16,-0.5,15.5,15,-0.5,14.5);

  TH1F* h_signal = new TH1F("signal","signal",4000,50.5,4050.5);
  TH2F* h_trig = new TH2F("trig","trig",16,-0.5,15.5,64,-0.5,63.5);
  TH2F* h_trig_xy = new TH2F("trig_xy","trig_xy",32,-90,90,32,-90,90);
  
  TH1F* h_bad = new TH1F("bad","bad",4000,50.5,4050.5);

  TH2F* h_first_retriggering = new TH2F("first_retriggering","first_retriggering",16,-0.5,15.5,64,-0.5,63.5);
  TH2F* h_all_retriggering = new TH2F("all_retriggering","all_retriggering",16,-0.5,15.5,64,-0.5,63.5);
  TH2F* h_first_retriggering_xy = new TH2F("first_retriggering_xy","first_retriggering_xy",32,-90,90,32,-90,90);
  TH2F* h_all_retriggering_xy = new TH2F("all_retriggering_xy","all_retriggering_xy",32,-90,90,32,-90,90);

  TH2F* h_dist_sca15[16];
  for(int ichip=0; ichip<16; ichip++) h_dist_sca15[ichip] = new TH2F(TString::Format("h_dist_sca15_chip%i",ichip),TString::Format("h_dist_sca15_chip%i",ichip),16,-0.5,15.5,401,-200.5,200.5);

  TH2F* h_dist_trig_retrig[16];
  for(int ichip=0; ichip<16; ichip++) h_dist_trig_retrig[ichip] = new TH2F(TString::Format("h_dist_trig_retrig_chip%i",ichip),TString::Format("h_dist_trig_retrig_chip%i",ichip),16,-0.5,15.5,401,-200.5,200.5);

  TH2F* h_dist_retrig_all[16];
  for(int ichip=0; ichip<16; ichip++) h_dist_retrig_all[ichip] = new TH2F(TString::Format("h_dist_retrig_all_chip%i",ichip),TString::Format("h_dist_retrig_all_chip%i",ichip),16,-0.5,15.5,401,-200.5,200.5); 
  
  // -----------------------------------------------------------------------------------------------------   
  // SCA analysis
  Long64_t nbytes = 0, nb = 0;
  // nentries=50000;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if ( jentry > 1000 && jentry % 1000 ==0 ) std::cout << "Progress: " << 100.*jentry/nentries <<" %"<<endl;

    for(int islboard=0; islboard<n_slboards; islboard++) {
      if(islboard != i_slboard) continue;
      
      for(int ichip=0; ichip<16; ichip++) {

	bool first_retrig=false;
      
	for(int isca=0; isca<15; isca++) {
	  bool retrig=false;
	  bool burst=false;
	  if(bcid[islboard][ichip][isca]<0) continue;
	
	  if(badbcid[islboard][ichip][isca]>2 && bcid[islboard][ichip][isca]>50 ) {
	    retrig=true;
	    if(first_retrig==false) {
	      first_retrig=true;
	      h_total_retrig_trains->Fill(ichip,isca);
	    }
	    h_total_retrig->Fill(ichip,isca);

	  }
	  if(bcid[islboard][ichip][isca]<50  ) {
	    burst=true;
	    h_total_burst->Fill(ichip,isca);
	  }
	
	  int hits_plane=0;
	  int hits_negative=0;
	  int nohits_negative=0;
	  for(int ichn=0;ichn<64;ichn++) {
	    if(gain_hit_high[islboard][ichip][isca][ichn]==1) hits_plane++; 
	    if(gain_hit_high[islboard][ichip][isca][ichn]==1 && charge_hiGain[islboard][ichip][isca][ichn]<100) hits_negative++;
	    if(gain_hit_high[islboard][ichip][isca][ichn]==0 && charge_hiGain[islboard][ichip][isca][ichn]<100)  nohits_negative++;
	  }
	
	  if( burst == false && retrig == false) {
	    bool bad = false;
	  
	    if(hits_plane>(maxnhit-1)) { bad=true; h_total_planeE->Fill(ichip,isca);}
	    if(hits_negative>(maxnhit-1)) { bad=true; h_total_negE->Fill(ichip,isca);}
	    if(nohits_negative>(maxnhit-1)) { bad=true; h_total_negE_nohit->Fill(ichip,isca);}
	

	    if(bad==false) {
	      h_total_trig->Fill(ichip,isca);
	      for(int ichn=0;ichn<64;ichn++)
		if(gain_hit_high[islboard][ichip][isca][ichn]==1) {
		  h_signal->Fill(charge_hiGain[islboard][ichip][isca][ichn]);
		  h_trig -> Fill(ichip,ichn);
		  h_trig_xy -> Fill(map_pointX[islboard][ichip][ichn],map_pointY[islboard][ichip][ichn]);    

		}
	    } else {
	      for(int ichn=0;ichn<64;ichn++) if(gain_hit_high[islboard][ichip][isca][ichn]==1) h_bad->Fill(charge_hiGain[islboard][ichip][isca][ichn]);
	    }
	  
	  } else {
	    for(int ichn=0;ichn<64;ichn++) if(gain_hit_high[islboard][ichip][isca][ichn]==1) h_bad->Fill(charge_hiGain[islboard][ichip][isca][ichn]);
	  }
	}
      }

    }// for islboards

    for(int islboard=0; islboard<n_slboards; islboard++) {
      if(islboard != i_slboard) continue;
      //channels starting retriggers
      for(int ichip=0; ichip<16; ichip++) {
	
	bool first_retrig=false;
      
	for(int isca=0; isca<15; isca++) {
	  if(bcid[islboard][ichip][isca]<0) continue;
	  if(bcid[islboard][ichip][isca]<50 || (bcid[islboard][ichip][isca]>890 && bcid[islboard][ichip][isca]<930)) continue;

	  if(isca==0) {
	    if(badbcid[islboard][ichip][isca]>2) first_retrig=true;
	    for(int ichn=0; ichn<64; ichn++) {
	      if(gain_hit_high[islboard][ichip][isca][ichn]==1) {
		h_first_retriggering -> Fill(ichip,ichn);
		h_all_retriggering -> Fill(ichip,ichn);

		h_first_retriggering_xy -> Fill(map_pointX[islboard][ichip][ichn],map_pointY[islboard][ichip][ichn]);
		h_all_retriggering_xy -> Fill(map_pointX[islboard][ichip][ichn],map_pointY[islboard][ichip][ichn]);     
	      }
	    }
	  }
	
	  if(isca > 0 && badbcid[islboard][ichip][isca]>2 && first_retrig==false) {
	    first_retrig=true;
	    for(int ichn=0; ichn<64; ichn++) {
	      if(gain_hit_high[islboard][ichip][isca][ichn]==1) {
		h_first_retriggering -> Fill(ichip,ichn);
		h_all_retriggering -> Fill(ichip,ichn);

		h_first_retriggering_xy -> Fill(map_pointX[islboard][ichip][ichn],map_pointY[islboard][ichip][ichn]);
		h_all_retriggering_xy -> Fill(map_pointX[islboard][ichip][ichn],map_pointY[islboard][ichip][ichn]);    
	      }
	    }
	  }

	  if(badbcid[islboard][ichip][isca]>2 && first_retrig==true){
	    for(int ichn=0; ichn<64; ichn++) {
	      if(gain_hit_high[islboard][ichip][isca][ichn]==1) {
		h_all_retriggering -> Fill(ichip,ichn);
		h_all_retriggering_xy -> Fill(map_pointX[islboard][ichip][ichn],map_pointY[islboard][ichip][ichn]);    

	      }
	    }
	  }

	}//end isca
      }//end chip
    }// islboard

    for(int islboard=0; islboard<n_slboards; islboard++) {
      if(islboard != i_slboard) continue;
      //if sca = 15, retrigger distance?
      for(int ichip=0; ichip<16; ichip++) {
	int sca=14;
	if(bcid[islboard][ichip][sca]<0) continue;

	for(int ichip2=0; ichip2<16; ichip2++) {
	  for(int isca=0; isca<15; isca++) {
	    if(badbcid[islboard][ichip2][isca]>2 && ichip2!=ichip) h_dist_sca15[ichip]->Fill(ichip2,bcid[islboard][ichip][sca]-bcid[islboard][ichip2][isca]);
	  }
	}
      }

      //if a trigger is recorded... when appears the next retrigger?
      for(int ichip=0; ichip<16; ichip++) {
	for(int isca=0; isca<15; isca++) {
	  if(bcid[islboard][ichip][isca]<0) continue;
	  if(bcid[islboard][ichip][isca]<50)  continue;

	  if(badbcid[islboard][ichip][isca]==0) {
	  
	    for(int ichip2=0; ichip2<16; ichip2++) {
	      for(int isca2=0; isca2<15; isca2++) {
		if(bcid[islboard][ichip2][isca2]<0) continue;
		if(bcid[islboard][ichip2][isca2]<50 ) continue;
		if(badbcid[islboard][ichip2][isca]>2 && ichip2!=ichip) h_dist_trig_retrig[ichip]->Fill(ichip2,bcid[islboard][ichip][isca]-bcid[islboard][ichip2][isca]);
	      }
	    }
	  
	  }
	}
      }
    }//islboard
	  
      
  }
  //end loop
  TFile *summary = new TFile("results_retriggers/Retriggers"+slboard+".root" , "RECREATE");
  summary->cd();

  h_total_planeE->Write();
  h_total_negE->Write();
  h_total_negE_nohit->Write();
  h_total_burst->Write();
  h_total_retrig->Write();
  h_total_retrig_trains->Write();
  h_total_trig->Write();

  
  h_first_retriggering->Write();
  h_all_retriggering->Write();
  h_first_retriggering_xy->Write();
  h_all_retriggering_xy->Write();
  h_signal->Write();
  h_trig->Write();
  h_trig_xy->Write();
  h_bad->Write();

  for(int ichip=0; ichip<16; ichip++) h_dist_sca15[ichip] ->Write();
  for(int ichip=0; ichip<16; ichip++) h_dist_trig_retrig[ichip] ->Write();

}


void DecodedSLBAnalysis::ReadMasked(TString filename) 
{

  std::ifstream reading_file(filename);
  if(!reading_file){
    cout<<" dameyo - damedame"<<endl;
  }
  for(int i=0; i<16; i++) {
    for(int j=0; j<64; j++) {
      masked[i][j] = 0;
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
    //    int masked_=0;
    //if(tmp_masked==0) masked_=1;
    //if(tmp_masked==1) masked_=0;
    masked[tmp_chip][tmp_channel] = tmp_masked;
  }

  double nmasked=0.;
  for(int i=0; i<16; i++) {
    for(int j=0; j<64; j++) {
      if(masked[i][j]==1) nmasked++;
    }
  }
  nmasked=100.*nmasked/1024;
  cout<< "In file " <<filename << " we read that "<<nmasked<<"% of channels are masked"<<endl;
  
}


void DecodedSLBAnalysis::ReadMap(TString filename, int slboard) 
{

  std::ifstream reading_file(filename);
  if(!reading_file){
    cout<<" dameyo - damedame"<<endl;
  }
  for(int i=0; i<16; i++) {
    for(int j=0; j<64; j++) {
      map_pointX[slboard][i][j] = -1000.;
      map_pointY[slboard][i][j] = -1000.;
    }
  }

  Int_t tmp_chip = 0,tmp_channel = 0;
  Float_t tmp_x0 = 0 ,tmp_y0 = 0 , tmp_x = 0 , tmp_y = 0 ;
  TString tmpst;
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst ;
  while(reading_file){
    reading_file >> tmp_chip >> tmp_x0 >> tmp_y0 >> tmp_channel >> tmp_x >> tmp_y ;
    map_pointX[slboard][tmp_chip][tmp_channel] = -tmp_x ;
    map_pointY[slboard][tmp_chip][tmp_channel] = -tmp_y ;
  }

}


// LANDAU STUFF
//---------------------------------------------------------------------------------------

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



TF1* DecodedSLBAnalysis::langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
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

