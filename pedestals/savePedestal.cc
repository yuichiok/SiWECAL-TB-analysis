#define savePedestal_cxx
#include "savePedestal.h"
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


void savePedestal::FindMasked(TString dif)
{

  //function that reads a root file and check which channels have or not signal (hit bit ==1).
  //should be used for root files that contain a full position scan, in order to provide meaninful list of masked channels.

  ofstream fout_masked("masked_"+dif+".txt",ios::out);

  fout_masked<<"#list of masked channels, per dif: "<<dif<<endl;
  fout_masked<<"#chip channel mask (0=not masked, 1=masked)"<<endl;

  
  bool global = true;
  int maxnhit=5;

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  // --------------------
  // all sca
  //  std::vector<TCanvas*> chip;
  std::vector<std::vector<TH1F*> > h_charge_channel;

  for(int ichip=0; ichip<16; ichip++) {
    
    std::vector<TH1F*> chargetemp_sca;
    for(int ichn=0; ichn<64; ichn++) {
      TH1F *tmp_charge_chn = new TH1F(TString::Format("tmp_charge_chip%i_chn%i",ichip,ichn),TString::Format("tmp_charge_chip%i_chn%i",ichip,ichn),4096,0.5,4096.5);
      chargetemp_sca.push_back(tmp_charge_chn);
    }
    h_charge_channel.push_back(chargetemp_sca);
  }


  
  // -----------------------------------------------------------------------------------------------------   
  // SCA analysis
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    for(int ichip=0; ichip<16; ichip++) {
      for(int isca=0; isca<15; isca++) {
	bool gooddata=true;
	if(global == true) {
	  int ntaggedasbad = 0;
	  for(int ichn=0; ichn<64; ichn++) {
	    if(charge_hiGain[ichip][isca][ichn]<10 && charge_hiGain[ichip][isca][ichn]>-1 ) 
	      ntaggedasbad++;
	  }//ichn 
	  if ( ntaggedasbad > 0) gooddata=false;
	}
	
	for(int ichn=0; ichn<64; ichn++) {

	  bool selection=false;
	  if(charge_hiGain[ichip][isca][ichn]>10 && badbcid[ichip][isca]==0 && nhits[ichip][isca]<maxnhit+1 && corrected_bcid[ichip][isca]>1250 && corrected_bcid[ichip][isca]<2900 ) selection=true;    
	  if(gain_hit_high[ichip][isca][ichn]==1 && selection==true && gooddata==true)
	    h_charge_channel.at(ichip).at(ichn)->Fill(charge_hiGain[ichip][isca][ichn]);
	}//ichn
           
      }//isca
    }//ichip 
  }  // end first loop analysis to fill pedestal historgrams


  //------------------------------------------------------------------
  // do signal analysis (chip/channel/sca based) analysis

  //first we calculate a first estimation of total number of hits per chip and total of channels non masked per chip
  double total_entries_chip[16];
  double total_channels_chip[16];
  for(int ichip=0; ichip<16; ichip++) {
    total_entries_chip[ichip]=0.;
    total_channels_chip[ichip]=0.;
  }
  for(int ichip=0; ichip<16; ichip++) {
    for(int ichn=0; ichn<64; ichn++) {
      total_entries_chip[ichip]+=h_charge_channel.at(ichip).at(ichn)->GetEntries();
      if(h_charge_channel.at(ichip).at(ichn)->GetEntries()>200) total_channels_chip[ichip]++;

    }
  }

  // the real estimation of masked channels is done by requiring a 0.3 of the estimated Nhits per channel for the whole data taking period run (per chip)
  // this is done because it has been observed that some times the masked channels in the xml slightly disagree within a set of runs.
  // possible reasons :
  // 1- a remasking is done during the mip runs, so we will have some channels that were not masked at the begining. we do not include them in the pedestal analysis since they are masked for following runs.
  // 2- human mistake in the masking for few runs.
  // 3- error in the slow control loading ?

  for(int ichip=0; ichip<16; ichip++) {
    for(int ichn=0; ichn<64; ichn++) {
      fout_masked << ichip <<" " <<ichn<< " "; 
      //minimum of 100 entries per SCA
      if(total_channels_chip[ichip]>0) {
	if(h_charge_channel.at(ichip).at(ichn)->GetEntries()>0.3*(total_entries_chip[ichip]/total_channels_chip[ichip])) fout_masked<<"0"<<endl;
	else fout_masked<<"1"<<endl;
      }
      else fout_masked<<"1"<<endl;    
    }
  }


}

void savePedestal::ReadMap(TString filename) 
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

void savePedestal::ReadMasked(TString filename) 
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

  while(reading_file){
    reading_file >> tmp_chip >> tmp_channel >> tmp_masked ;
    masked[tmp_chip][tmp_channel] = tmp_masked ;
   }

}

void savePedestal::PedestalAnalysis(TString dif,TString grid="",TString map_filename="../fev10_chip_channel_x_y_mapping.txt")
{

  //Read the channel/chip -- x/y mapping
  ReadMap(map_filename);
  //Read the list of masked channels
  ReadMasked("masked_"+dif+".txt");

  if(grid!="") dif=dif+"_"+grid;

  ofstream fout_ped("Pedestal_"+dif+".txt",ios::out);

  fout_ped<<"#pedestal results (fit to a gaussian) remove channels/sca with two pedestals peaks from the analysis : "<<dif<<endl;
  fout_ped<<"#chip channel ped0 eped1 ped1 eped1 ... ped14 ped14 (all SCA)"<<endl;

  
  bool global = true;
  int maxnhit=5;

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  // std::vector<std::vector<std::vector<float> > > ped_mean;
  // std::vector<std::vector<std::vector<float> > > ped_rms;

  // std::vector<std::vector<std::vector<float> > > ped_tagged_mean;
  // std::vector<std::vector<std::vector<float> > > ped_tagged_rms;
  // for(int ichip=0; ichip<16; ichip++) {
  //   ped_mean.push_back(std::vector<std::vector<float> >() );
  //   ped_rms.push_back(std::vector<std::vector<float> >() );

  //   ped_tagged_mean.push_back(std::vector<std::vector<float> >() );
  //   ped_tagged_rms.push_back(std::vector<std::vector<float> >() );

  //   for(int ichn=0; ichn<64; ichn++) {
  //     ped_mean.at(ichip).push_back(std::vector<float> () );
  //     ped_rms.at(ichip).push_back(std::vector<float> () );

  //     ped_tagged_mean.at(ichip).push_back(std::vector<float> () );
  //     ped_tagged_rms.at(ichip).push_back(std::vector<float> () );
  //   }
  // }
  //


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
  std::vector<TH1F*> pedestal_diff_chip ;
  std::vector<TH1F*> pedestal_tagged_diff_chip ;


  for(int ichip=0; ichip<16; ichip++) {
    TH1F *ped_chip = new TH1F(TString::Format("ped_chip%i",ichip),TString::Format("ped_chip%i",ichip),1000,0.5,1000.5);
    pedestal_chip.push_back(ped_chip);

    TH1F *ped_tagged_chip = new TH1F(TString::Format("ped_tagged_chip%i",ichip),TString::Format("ped_tagged_chip%i",ichip),1000,0.5,1000.5);
    pedestal_tagged_chip.push_back(ped_tagged_chip);

    TH1F *ped_diff_chip = new TH1F(TString::Format("ped_diff_chip%i",ichip),TString::Format("ped_diff_chip%i",ichip),1002,-500,500);
    pedestal_diff_chip.push_back(ped_diff_chip);
    
    TH1F *ped_tagged_diff_chip = new TH1F(TString::Format("ped_tagged_diff_chip%i",ichip),TString::Format("ped_tagged_diff_chip%i",ichip),1002,-500,500);
    pedestal_tagged_diff_chip.push_back(ped_tagged_diff_chip);

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
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    for(int ichip=0; ichip<16; ichip++) {

      for(int isca=0; isca<15; isca++) {

	bool gooddata=true;
	if(global == true) {
	  int ntaggedasbad = 0;
	  for(int ichn=0; ichn<64; ichn++) {
	    if(charge_hiGain[ichip][isca][ichn]<10 && charge_hiGain[ichip][isca][ichn]>-1 ) 
	      ntaggedasbad++;
	  }//ichn 
	  if ( ntaggedasbad > 0) gooddata=false;
	}
	
	for(int ichn=0; ichn<64; ichn++) {

	  //good events
	  bool selection=false;
	  if(charge_hiGain[ichip][isca][ichn]>10 && badbcid[ichip][isca]==0 && nhits[ichip][isca]<maxnhit+1 && corrected_bcid[ichip][isca]>1247 && corrected_bcid[ichip][isca]<2900 ) selection=true;
	  if(masked[ichip][ichn]==1) selection=false;
 	  if(gain_hit_high[ichip][isca][ichn]==0 && selection==true && gooddata==true)
	    ped_sca.at(ichip).at(ichn).at(isca)->Fill(charge_hiGain[ichip][isca][ichn]);

	  //bad events
	  selection=false;
	  if( ( badbcid[ichip][isca]>0 || nhits[ichip][isca]>maxnhit || gooddata==false) && (corrected_bcid[ichip][isca]>1247 && corrected_bcid[ichip][isca]<2900)  ) selection=true;
	  if(masked[ichip][ichn]==1) selection=false;
 	  if(gain_hit_high[ichip][isca][ichn]==0 && selection==true )
	    ped_sca_tagged.at(ichip).at(ichn).at(isca)->Fill(charge_hiGain[ichip][isca][ichn]);

	}
           
      }//isca

    }//ichip 
   
  }  // end first loop analysis to fill pedestal historgrams


  TFile *pedfile = new TFile("Pedestal_"+dif+".root" , "RECREATE");
  pedfile->cd();

  
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

	pedestal_entries_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , ped_sca.at(ichip).at(ichn).at(isca)->GetEntries());
	pedestal_tagged_entries_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , ped_sca_tagged.at(ichip).at(ichn).at(isca)->GetEntries());

	
	if(ped_sca.at(ichip).at(ichn).at(isca)->GetEntries()> 300 ){ //max_entries/2 ) {
	  TSpectrum *s = new TSpectrum();
	  int npeaks = s->Search(ped_sca.at(ichip).at(ichn).at(isca),2,"",0.2); 
	  pedestal_npeaks_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , npeaks);


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
	      if(ipeak != npeak_max) pedestal_diff_chip.at(ichip)->Fill( mean_peak[npeak_max] - mean_peak[ipeak]);
	    }

	  }

	  if(npeaks == 1) {
	    Double_t *mean_peak=s->GetPositionX();
	    
	    TF1 *f0 = new TF1("f0","gaus",mean_peak[0]-2*ped_sca.at(ichip).at(ichn).at(isca)->GetRMS(),mean_peak[0]+2*ped_sca.at(ichip).at(ichn).at(isca)->GetRMS());
	    ped_sca.at(ichip).at(ichn).at(isca)->Fit("f0","RQNOC");

	    TF1 *f1 = new TF1("f1","gaus",f0->GetParameter(1)-2.*f0->GetParameter(2),f0->GetParameter(1)+2.*f0->GetParameter(2));
	    ped_sca.at(ichip).at(ichn).at(isca)->Fit("f1","RQME");
	    // ped_mean.at(ichip).at(ichn).push_back(f1->GetParameter(1));
	    //  ped_rms.at(ichip).at(ichn).push_back(f1->GetParameter(2));
	      fout_ped<<f1->GetParameter(1) << " " << f1->GetParError(1)<<" ";
	      pedestal_chip.at(ichip)->Fill(f1->GetParameter(1));

	      pedestal_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParameter(1));
	      pedestal_width_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParameter(2));
	      pedestal_error_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParError(1));
	      pedestal_chi2ndf_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetChisquare() / f1->GetNDF());
	      
	  } else {
	    // ped_mean.at(ichip).at(ichn).push_back(0);
	    //  ped_rms.at(ichip).at(ichn).push_back(0);
	    fout_ped<<0<< " " << 0<<" ";
	  }
	} else {
	  // ped_mean.at(ichip).at(ichn).push_back(0);
          //ped_rms.at(ichip).at(ichn).push_back(0);
	  fout_ped<<0<< " " << 0<<" ";
	}

	// analyze pedestal for tagged events
	if(ped_sca_tagged.at(ichip).at(ichn).at(isca)->GetEntries()> 150 ){
	  TSpectrum *s = new TSpectrum();
	  int npeaks = s->Search(ped_sca_tagged.at(ichip).at(ichn).at(isca));
	  pedestal_tagged_npeaks_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , npeaks);

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
              if(ipeak != npeak_max) pedestal_tagged_diff_chip.at(ichip)->Fill( mean_peak[npeak_max] - mean_peak[ipeak] );
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
		
		pedestal_tagged_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParameter(1));
		pedestal_tagged_width_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParameter(2));
		pedestal_tagged_error_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParError(1));
		pedestal_tagged_chi2ndf_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetChisquare() / f1->GetNDF());
	      }
	    }
	  }
	}


      }
      fout_ped<<endl;
    }
  }

  pedfile->Close();

  TFile *pedfile_summary = new TFile("Pedestal_summary_"+dif+".root" , "RECREATE");
  pedfile_summary->cd();
  
  // good pedestal events (not tagged events)
  TCanvas *canvas_pedestal_map = new TCanvas("pedestal_map_"+dif, "pedestal_map_"+dif,1200,1200);
  canvas_pedestal_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString dif0=dif+TString::Format("_sca%i",isca);
    canvas_pedestal_map->cd(isca+1);
    pedestal_map[isca]->SetStats(kFALSE);
    pedestal_map[isca]->SetTitle("pedestal_map, "+dif0);
    pedestal_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_map[isca]->GetZaxis()->SetRangeUser(200,400);

    //noise.at(ichip)->SetLineColor(2);
    pedestal_map[isca]->Draw("colz");
    pedestal_map[isca]->Write();
   }

  TCanvas *canvas_pedestal_width_map = new TCanvas("pedestal_width_map_"+dif, "pedestal_width_map_"+dif,1200,1200);
  canvas_pedestal_width_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString dif0=dif+TString::Format("_sca%i",isca);
    canvas_pedestal_width_map->cd(isca+1);
    pedestal_width_map[isca]->SetStats(kFALSE);
    pedestal_width_map[isca]->SetTitle("pedestal_width_map, "+dif0);
    pedestal_width_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_width_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_width_map[isca]->GetZaxis()->SetRangeUser(0,7);
    //noise.at(ichip)->SetLineColor(2);
    pedestal_width_map[isca]->Draw("colz");
    pedestal_width_map[isca]->Write();
   }

  TCanvas *canvas_pedestal_error_map = new TCanvas("pedestal_error_map_"+dif, "pedestal_error_map_"+dif,1200,1200);
  canvas_pedestal_error_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString dif0=dif+TString::Format("_sca%i",isca);
    canvas_pedestal_error_map->cd(isca+1);
    pedestal_error_map[isca]->SetStats(kFALSE);
    pedestal_error_map[isca]->SetTitle("pedestal_error_map, "+dif0);
    pedestal_error_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_error_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_error_map[isca]->Draw("colz");
    pedestal_error_map[isca]->Write();
   }

  TCanvas *canvas_pedestal_npeaks_map = new TCanvas("pedestal_npeaks_map_"+dif, "pedestal_npeaks_map_"+dif,1200,1200);
  canvas_pedestal_npeaks_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString dif0=dif+TString::Format("_sca%i",isca);
    canvas_pedestal_npeaks_map->cd(isca+1);
    pedestal_npeaks_map[isca]->SetStats(kFALSE);
    pedestal_npeaks_map[isca]->SetTitle("pedestal_npeaks_map, "+dif0);
    pedestal_npeaks_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_npeaks_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_npeaks_map[isca]->GetZaxis()->SetRangeUser(0.5,4.5);
    //noise.at(ichip)->SetLineColor(2);
    pedestal_npeaks_map[isca]->Draw("colz");
    pedestal_npeaks_map[isca]->Write();
   }

  TCanvas *canvas_pedestal_chi2ndf_map = new TCanvas("pedestal_chi2ndf_map_"+dif, "pedestal_chi2ndf_map_"+dif,1200,1200);
  canvas_pedestal_chi2ndf_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString dif0=dif+TString::Format("_sca%i",isca);
    canvas_pedestal_chi2ndf_map->cd(isca+1);
    pedestal_chi2ndf_map[isca]->SetStats(kFALSE);
    pedestal_chi2ndf_map[isca]->SetTitle("pedestal_chi2ndf_map, "+dif0);
    pedestal_chi2ndf_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_chi2ndf_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_chi2ndf_map[isca]->Draw("colz");
    pedestal_chi2ndf_map[isca]->Write();
   }

  TCanvas *canvas_pedestal_entries_map = new TCanvas("pedestal_entries_map_"+dif, "pedestal_entries_map_"+dif,1200,1200);
  canvas_pedestal_entries_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString dif0=dif+TString::Format("_sca%i",isca);
    canvas_pedestal_entries_map->cd(isca+1);
    pedestal_entries_map[isca]->SetStats(kFALSE);
    pedestal_entries_map[isca]->SetTitle("pedestal_entries_map, "+dif0);
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

   
  TCanvas *canvas_pedestal_pedestal = new TCanvas("pedestal_average_"+dif, "pedestal_average_"+dif,1200,1200);
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

  TCanvas *canvas_pedestal_diff = new TCanvas("pedestal_diff_"+dif, "pedestal_diff_"+dif,1200,1200);
  canvas_pedestal_diff->Divide(4,4);
  for(int ichip=0; ichip<16; ichip++) {
    canvas_pedestal_diff->cd(ichip+1);

    pedestal_diff_chip.at(ichip)->GetXaxis()->SetRangeUser(-100,100);
    pedestal_diff_chip.at(ichip)->SetTitle(TString::Format("Pedestal diff, chip-%i",ichip));
    pedestal_diff_chip.at(ichip)->GetXaxis()->SetTitle("ADC");
    pedestal_diff_chip.at(ichip)->GetYaxis()->SetTitle("#");
    pedestal_diff_chip.at(ichip)->Draw("hs");
  }

  canvas_pedestal_diff->Write();


  // Tagged events
  TCanvas *canvas_pedestal_tagged_map = new TCanvas("pedestal_tagged_map_"+dif, "pedestal_tagged_map_"+dif,1200,1200);
  canvas_pedestal_tagged_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString dif0=dif+TString::Format("_sca%i",isca);
    canvas_pedestal_tagged_map->cd(isca+1);
    pedestal_tagged_map[isca]->SetStats(kFALSE);
    pedestal_tagged_map[isca]->SetTitle("pedestal_tagged_map, "+dif0);
    pedestal_tagged_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_tagged_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_tagged_map[isca]->GetZaxis()->SetRangeUser(200,400);

    //noise.at(ichip)->SetLineColor(2);
    pedestal_tagged_map[isca]->Draw("colz");
    pedestal_tagged_map[isca]->Write();
   }

  TCanvas *canvas_pedestal_tagged_width_map = new TCanvas("pedestal_tagged_width_map_"+dif, "pedestal_tagged_width_map_"+dif,1200,1200);
  canvas_pedestal_tagged_width_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString dif0=dif+TString::Format("_sca%i",isca);
    canvas_pedestal_tagged_width_map->cd(isca+1);
    pedestal_tagged_width_map[isca]->SetStats(kFALSE);
    pedestal_tagged_width_map[isca]->SetTitle("pedestal_tagged_width_map, "+dif0);
    pedestal_tagged_width_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_tagged_width_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_tagged_width_map[isca]->GetZaxis()->SetRangeUser(0,7);
    //noise.at(ichip)->SetLineColor(2);
    pedestal_tagged_width_map[isca]->Draw("colz");
    pedestal_tagged_width_map[isca]->Write();
   }

  TCanvas *canvas_pedestal_tagged_error_map = new TCanvas("pedestal_tagged_error_map_"+dif, "pedestal_tagged_error_map_"+dif,1200,1200);
  canvas_pedestal_tagged_error_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString dif0=dif+TString::Format("_sca%i",isca);
    canvas_pedestal_tagged_error_map->cd(isca+1);
    pedestal_tagged_error_map[isca]->SetStats(kFALSE);
    pedestal_tagged_error_map[isca]->SetTitle("pedestal_tagged_error_map, "+dif0);
    pedestal_tagged_error_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_tagged_error_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_tagged_error_map[isca]->Draw("colz");
    pedestal_tagged_error_map[isca]->Write();
   }

  TCanvas *canvas_pedestal_tagged_npeaks_map = new TCanvas("pedestal_tagged_npeaks_map_"+dif, "pedestal_tagged_npeaks_map_"+dif,1200,1200);
  canvas_pedestal_tagged_npeaks_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString dif0=dif+TString::Format("_sca%i",isca);
    canvas_pedestal_tagged_npeaks_map->cd(isca+1);
    pedestal_tagged_npeaks_map[isca]->SetStats(kFALSE);
    pedestal_tagged_npeaks_map[isca]->SetTitle("pedestal_tagged_npeaks_map, "+dif0);
    pedestal_tagged_npeaks_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_tagged_npeaks_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_tagged_npeaks_map[isca]->GetZaxis()->SetRangeUser(0.5,4.5);
    //noise.at(ichip)->SetLineColor(2);
    pedestal_tagged_npeaks_map[isca]->Draw("colz");
    pedestal_tagged_npeaks_map[isca]->Write();
   }

  TCanvas *canvas_pedestal_tagged_chi2ndf_map = new TCanvas("pedestal_tagged_chi2ndf_map_"+dif, "pedestal_tagged_chi2ndf_map_"+dif,1200,1200);
  canvas_pedestal_tagged_chi2ndf_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString dif0=dif+TString::Format("_sca%i",isca);
    canvas_pedestal_tagged_chi2ndf_map->cd(isca+1);
    pedestal_tagged_chi2ndf_map[isca]->SetStats(kFALSE);
    pedestal_tagged_chi2ndf_map[isca]->SetTitle("pedestal_tagged_chi2ndf_map, "+dif0);
    pedestal_tagged_chi2ndf_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_tagged_chi2ndf_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_tagged_chi2ndf_map[isca]->Draw("colz");
    pedestal_tagged_chi2ndf_map[isca]->Write();
   }

  TCanvas *canvas_pedestal_tagged_entries_map = new TCanvas("pedestal_tagged_entries_map_"+dif, "pedestal_tagged_entries_map_"+dif,1200,1200);
  canvas_pedestal_tagged_entries_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString dif0=dif+TString::Format("_sca%i",isca);
    canvas_pedestal_tagged_entries_map->cd(isca+1);
    pedestal_tagged_entries_map[isca]->SetStats(kFALSE);
    pedestal_tagged_entries_map[isca]->SetTitle("pedestal_tagged_entries_map, "+dif0);
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

   
  TCanvas *canvas_pedestal_tagged_pedestal = new TCanvas("pedestal_tagged_average_"+dif, "pedestal_tagged_average_"+dif,1200,1200);
  canvas_pedestal_tagged_pedestal->Divide(4,4);
  for(int ichip=0; ichip<16; ichip++) {
    canvas_pedestal_tagged_pedestal->cd(ichip+1);
    //gPad->SetLogy();
    pedestal_tagged_chip.at(ichip)->GetXaxis()->SetRangeUser(0,500);
    pedestal_tagged_chip.at(ichip)->SetTitle(TString::Format("Average pedestal_tagged, chip-%i",ichip));
    pedestal_tagged_chip.at(ichip)->GetXaxis()->SetTitle("ADC");
    pedestal_tagged_chip.at(ichip)->GetYaxis()->SetTitle("#");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_tagged_chip.at(ichip)->Draw("hs");
    //pedestal_tagged_chip.at(ichip)->Write();
  }

  canvas_pedestal_tagged_pedestal->Write();


  TCanvas *canvas_pedestal_tagged_diff = new TCanvas("pedestal_tagged_diff_"+dif, "pedestal_tagged_diff_"+dif,1200,1200);
  canvas_pedestal_tagged_diff->Divide(4,4);
  for(int ichip=0; ichip<16; ichip++) {
    canvas_pedestal_tagged_diff->cd(ichip+1);

    pedestal_tagged_diff_chip.at(ichip)->GetXaxis()->SetRangeUser(-100,100);
    pedestal_tagged_diff_chip.at(ichip)->SetTitle(TString::Format("Pedestal diff, chip-%i",ichip));
    pedestal_tagged_diff_chip.at(ichip)->GetXaxis()->SetTitle("ADC");
    pedestal_tagged_diff_chip.at(ichip)->GetYaxis()->SetTitle("#");
    pedestal_tagged_diff_chip.at(ichip)->Draw("hs");
  }

  canvas_pedestal_tagged_diff->Write();

  pedfile_summary->Close();



}


void savePedestal::BcidCorrelations(TString dif)
{

  ReadMasked("masked_"+dif+".txt");

  //function that reads a root file and check which channels have or not signal (hit bit ==1).
  //should be used for root files that contain a full position scan, in order to provide meaninful list of masked channels.
 
  bool global = true;
  int maxnhit=65;

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  // --------------------
  // all sca
  //  std::vector<TCanvas*> chip;
  TH1F* h_bcid_dif_retrig = new TH1F("bcid_dif_retrig_dif_"+dif,"bcid_dif_retrig_dif_"+dif,8000,0,800000);
  TH2F* h_bcid_corr_retrig = new TH2F("bcid_correl_retrig_dif_"+dif,"bcid_correl_retrig_dif_"+dif,800,0,80000,800,0,80000);
  TH1F* h_bcid_dif_total_retrig = new TH1F("bcid_retrig_dif_"+dif,"bcid_retrig_dif_"+dif,8002,-4000,4000);
  TH1F* h_bcid_dif_plane = new TH1F("bcid_dif_plane_dif_"+dif,"bcid_dif_plane_dif_"+dif,8000,0,800000);
  TH2F* h_bcid_corr_plane = new TH2F("bcid_correl_plane_dif_"+dif,"bcid_correl_plane_dif_"+dif,800,0,80000,800,0,80000);
  TH1F* h_bcid_dif_total_plane = new TH1F("bcid_plane_dif_"+dif,"bcid_plane_dif_"+dif,8002,-4000,4000);
  TH1F* h_bcid_dif_good = new TH1F("bcid_good_dif_"+dif,"bcid_good_dif_"+dif,8002,-4000,4000);

  // -----------------------------------------------------------------------------------------------------   
  // SCA analysis
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    for(int ichip=0; ichip<16; ichip++) {
      for(int isca=0; isca<15; isca++) {
	bool gooddata=true;
	if(global == true) {
	  int ntaggedasbad = 0;
	  for(int ichn=0; ichn<64; ichn++) {
	    if(charge_hiGain[ichip][isca][ichn]<10 && charge_hiGain[ichip][isca][ichn]>-1 ) 
	      ntaggedasbad++;
	  }//ichn 
	  if ( ntaggedasbad > 0) gooddata=false;
	}

	for(int ichn=0; ichn<64; ichn++) {

	  
	  bool selection=false;
	  if(gain_hit_high[ichip][isca][ichn]==1 && charge_hiGain[ichip][isca][ichn]>10 && badbcid[ichip][isca]==0 && nhits[ichip][isca]<maxnhit+1 && corrected_bcid[ichip][isca]>1247 && corrected_bcid[ichip][isca]<2900 ) selection=true;


	  
	  // if(masked[ichip][ichn]==1) selection=false;
	  
	  double chip_bcid=ichip*5000.+bcid[ichip][isca];
	  
	  if(selection == true) {
	    for(int ichip2=0; ichip2<16; ichip2++) {
	      for(int isca2=0; isca2<15; isca2++) {
		double chip_bcid2=ichip2*5000.+bcid[ichip2][isca2];

		for(int ichn2=0; ichn2<64; ichn2++) {
		  bool selection2=false;
		  if(gain_hit_high[ichip2][isca2][ichn2]==1 && charge_hiGain[ichip2][isca2][ichn2]>10 && badbcid[ichip2][isca2]==3 && nhits[ichip2][isca2]<maxnhit+1) selection2=true;
	  
		  //  if(masked[ichip2][ichn2]==1) selection2=false;
		  if(selection2==true) {
		    h_bcid_dif_retrig->Fill(ichip2*50000.+(bcid[ichip][isca]-bcid[ichip2][isca2]));
		    h_bcid_dif_total_retrig->Fill(bcid[ichip][isca]-bcid[ichip2][isca2]);
		    h_bcid_corr_retrig->Fill(chip_bcid,chip_bcid2);
		  }

		  bool selection3=false;
		  if(gain_hit_high[ichip2][isca2][ichn2]==1 && badbcid[ichip2][isca2]>0 && nhits[ichip2][isca2]>maxnhit) selection3=true;
		  if(selection3==true) {
		    h_bcid_dif_plane->Fill(ichip2*50000.+(bcid[ichip][isca]-bcid[ichip2][isca2]));
		    h_bcid_dif_total_plane->Fill(bcid[ichip][isca]-bcid[ichip2][isca2]);
		    h_bcid_corr_plane->Fill(chip_bcid,chip_bcid2);
		  }

		    if(ichip2==ichip &&  isca2>isca && ichn==ichn2 && gain_hit_high[ichip2][isca2][ichn2]==1 && charge_hiGain[ichip2][isca2][ichn2]>10 && badbcid[ichip2][isca2]==0 && nhits[ichip2][isca2]<maxnhit+1 && corrected_bcid[ichip2][isca2]>1247 && corrected_bcid[ichip2][isca2]<2900 ) 	h_bcid_dif_good->Fill(corrected_bcid[ichip][isca]-corrected_bcid[ichip2][isca2]);
		}
		
	      }
	    }
	  }
	}
           
      }//isca

    }//ichip 
  }  // end first loop analysis to fill pedestal historgrams

  
 
  TCanvas *canvas= new TCanvas("bcid_correl_"+dif, "bcid_correl_"+dif,1600,1000);
  canvas->Divide(3,2);
  canvas->cd(1);
  h_bcid_corr_retrig->GetXaxis()->SetTitle("bcid hit (+ 5000*ichip)");
  h_bcid_corr_retrig->GetYaxis()->SetTitle("bcid retrig (+ 5000*ichip)");   
  h_bcid_corr_retrig->Draw("colz");

  canvas->cd(2);
  h_bcid_dif_retrig->GetXaxis()->SetTitle("50000*ichip + (bcid_hit-bcid_retrig)");
  h_bcid_dif_retrig->Draw("h");

  canvas->cd(3);
  h_bcid_dif_total_retrig->GetXaxis()->SetTitle("bcid_hit-bcid_retrig (all)");
  h_bcid_dif_total_retrig->GetXaxis()->SetRangeUser(-40,40);
  h_bcid_dif_total_retrig->Draw("h");

  canvas->cd(4);
  h_bcid_corr_plane->GetXaxis()->SetTitle("bcid hit (+ 5000*ichip)");
  h_bcid_corr_plane->GetYaxis()->SetTitle("bcid plane (+ 5000*ichip)");   
  h_bcid_corr_plane->Draw("colz");

  canvas->cd(5);
  h_bcid_dif_plane->GetXaxis()->SetTitle("50000*ichip + (bcid_hit-bcid_plane)");
  h_bcid_dif_plane->Draw("h");

  
  canvas->cd(6);
  h_bcid_dif_total_plane->GetXaxis()->SetTitle("bcid_hit-bcid_plane (all)");
  h_bcid_dif_total_plane->GetXaxis()->SetRangeUser(-40,40);
  h_bcid_dif_total_plane->Draw("h");

  canvas->Print("bcid_correlations_"+dif+".root");
  
  
}
