#define singleSlabAnalysis_cxx
#include "singleSlabAnalysis.h"
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


void singleSlabAnalysis::SignalAnalysis(TString dif="dif_1_1_1", TString outputname="", bool readpedestal=false, TString map_filename="/home/calice/SiWECAL-TB-analysis/fev10_chip_channel_x_y_mapping.txt")
{

  int maxnhit=2; // plane event threshold

  //Read the channel/chip -- x/y mapping
  ReadMap(map_filename);
  //Read the list of pedestals (this information contains, implicitily, the masked channels information )
  if(readpedestal==true) ReadPedestals("results_pedestal/pedestals/Pedestal_"+dif+".txt");

  ofstream fout_mip("results_mipcalibration/MIP_"+dif+outputname+".txt",ios::out);
  fout_mip<<"#mip results "<<dif<<endl;
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
	
    for(int ichn=0;ichn<64;ichn++) {
      TString histo_title=TString::Format("charge_%s_chip%i_chn%i",dif.Data(),ichip,ichn);
      TH1F *chn_mip = new TH1F(histo_title,histo_title,4197,-100.5,4096.5);
      tmpmip_histo.push_back(chn_mip);
      TString s_n_title=TString::Format("s_n_%s_chip%i_chn%i",dif.Data(),ichip,ichn);
      TH1F *chn_s_n = new TH1F(s_n_title,s_n_title,4197,-100.5,4096.5);
      tmps_n_histo.push_back(chn_s_n);
    }      
    mip_histo.push_back(tmpmip_histo);
    s_n_histo.push_back(tmps_n_histo);
  }

  //maps 
  TH2F* MPV_2d = new TH2F("MPV_2d_"+dif,"MPV_2d_"+dif,32,-90,90,32,-90,90);
  TH2F* eMPV_2d = new TH2F("eMPV_2d_"+dif,"eMPV_2d_"+dif,32,-90,90,32,-90,90);
  TH2F* widthMPV_2d = new TH2F("widthMPV_2d_"+dif,"widthMPV_2d_"+dif,32,-90,90,32,-90,90);
  TH2F* chi2NDF_2d = new TH2F("chi2NDF_2d_"+dif,"chi2NDF_2d_"+dif,32,-90,90,32,-90,90);
  TH2F* mipEntries_2d = new TH2F("mipEntries_2d_"+dif,"mipEntries_2d_"+dif,32,-90,90,32,-90,90);
  TH2F* S_N_2d = new TH2F("S_N_2d_"+dif,"S_N_2d_"+dif,32,-90,90,32,-90,90);

  TH2F* correl_S_N_vs_nhits = new TH2F("correl_S_N_vs_nhits"+dif,"correl_S_N_vs_nhits"+dif,80,0.25,40.25,20,1225,51225);
  TH2F* correl_S_N_vs_mpv = new TH2F("correl_S_N_vs_mpv"+dif,"correl_S_N_vs_mpv"+dif,80,0.25,40.25,40,40.5,80.5);
  TH2F* correl_mpv_vs_nhits = new TH2F("correl_mpv_vs_nhits"+dif,"correl_mpv_vs_nhits"+dif,40,40.5,80.5,20,1225,51225);
  TH2F* correl_widthped_vs_nhits = new TH2F("correl_widthped_vs_nhits"+dif,"correl_widthped_vs_nhits"+dif,50,2.05,7.05,20,1225,51225);
  TH2F* correl_S_N_vs_widthmpv = new TH2F("correl_S_N_vs_widthmpv"+dif,"correl_S_N_vs_widthmpv"+dif,80,0.25,40.25,50,0.5,50.5);
  TH2F* correl_S_N_vs_widthped = new TH2F("correl_S_N_vs_widthped"+dif,"correl_S_N_vs_widthped"+dif,80,0.25,40.25,50,2.05,7.05);


  //summary plots
  TH1F* MPV_chip[16];
  TH1F* S_N_chip[16];
  for(int ichip=0; ichip<16; ichip++) {
    MPV_chip[ichip]= new TH1F(TString::Format("MPV_%s_chip%i",dif.Data(),ichip),TString::Format("MPV_%s_chip%i",dif.Data(),ichip),1500,0.05,150.05);
    S_N_chip[ichip]= new TH1F(TString::Format("S_N_%s_chip%i",dif.Data(),ichip),TString::Format("S_N_%s_chip%i",dif.Data(),ichip),500,0.05,50.05);
  }

  TH1F* MPV_slab= new TH1F(TString::Format("MPV_%s",dif.Data()),TString::Format("MPV_%s",dif.Data()),1500,0.05,150.05);
  TH1F* S_N_slab= new TH1F(TString::Format("S_N_%s",dif.Data()),TString::Format("S_N_%s",dif.Data()),1000,0.05,100.05);
 
  // -----------------------------------------------------------------------------------------------------   
  // Signal readout
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
     

    for(int ichip=0; ichip<16; ichip++) {
      for(int ichn=0; ichn<64; ichn++) {
	for(int isca=0; isca<15; isca++) {

	  if( ped_mean.at(ichip).at(ichn).at(isca)>50 &&  ped_width.at(ichip).at(ichn).at(isca)>0 ) {

	    bool selection=false;
	    if(charge_hiGain[ichip][isca][ichn]>0 && gain_hit_high[ichip][isca][ichn]==1 && badbcid[ichip][isca]==0  && nhits[ichip][isca]<maxnhit+1) selection=true;
	      
	    if(selection==true) {
	      mip_histo.at(ichip).at(ichn)->Fill(charge_hiGain[ichip][isca][ichn]-ped_mean.at(ichip).at(ichn).at(isca));
	      s_n_histo.at(ichip).at(ichn)->Fill( (charge_hiGain[ichip][isca][ichn]-ped_mean.at(ichip).at(ichn).at(isca)) / ped_width.at(ichip).at(ichn).at(isca));
	    }
	  }
    
    	}//ichn
      }//sca
    }//ichip
    
  }

  // -----------------------------------------------------------------------------------------------------   
  // Signal analysis
  TFile *signalfile_summary = new TFile("results_mipcalibration/Signal_summary_"+dif+outputname+".root" , "RECREATE");
  signalfile_summary->cd();
  TDirectory *cdhisto = signalfile_summary->mkdir("histograms");


  // Setting fit range and start values
  Double_t fr[2];
  Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];   
  pllo[0]=0.; pllo[1]=0.0; pllo[2]=10.0; pllo[3]=0.;
  plhi[0]=100.0; plhi[1]=200.0; plhi[2]=100000000.0; plhi[3]=20.0;
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
    for(int ichn=0; ichn<64; ichn++) {
            
      if(mip_histo.at(ichip).at(ichn)->GetEntries()> 500 && mip_histo.at(ichip).at(ichn)->GetMean()>20){

	fr[0]=mip_histo.at(ichip).at(ichn)->GetMean()-0.8*mip_histo.at(ichip).at(ichn)->GetRMS();
	fr[1]=mip_histo.at(ichip).at(ichn)->GetMean();//+0.1*mip_histo.at(ichip).at(ichn)->GetRMS();
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

	MPV_2d->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn],mpv);
	eMPV_2d->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn],empv);
        widthMPV_2d->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn],wmpv);
	chi2NDF_2d->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn],chi2ndf);
	mipEntries_2d->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn],mipentries);
	MPV_chip[ichip]->Fill(mpv);
	MPV_slab->Fill(mpv);

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

        S_N_2d->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn],s_n_temp);
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
      } else {
	fout_mip<<ichip<<" "<<ichn<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<endl;

      }
    }
  }


  signalfile_summary->cd();

  TCanvas *canvas_maps = new TCanvas(TString::Format("Signal_%s",dif.Data()),TString::Format("Signal_%s",dif.Data()),1200,800);
  canvas_maps->Divide(2,2);
  TCanvas *canvas_summary_mip = new TCanvas(TString::Format("mip_perchip_%s",dif.Data()),TString::Format("mip_perchip_%s",dif.Data()),1200,800);
  canvas_summary_mip->Divide(4,4);
  TCanvas *canvas_summary_s_n = new TCanvas(TString::Format("s_n_perchip_%s",dif.Data()),TString::Format("s_n_perchip_%s",dif.Data()),1200,1200);
  canvas_summary_s_n->Divide(4,4); 
  TCanvas *canvas_summary = new TCanvas(TString::Format("summary_%s",dif.Data()),TString::Format("summary_%s",dif.Data()),1200,600);
  canvas_summary->Divide(2,1);

  TCanvas *canvas_correlations_S_N = new TCanvas(TString::Format("correlations_S_N_%s",dif.Data()),TString::Format("correlations_S_N_%s",dif.Data()),800,900);
  canvas_correlations_S_N->Divide(2,3);

  
  // SUMMARY MAPS
  eMPV_2d->SetTitle("eMIP[ADC] map, "+dif);
  eMPV_2d->GetXaxis()->SetTitle("x");
  eMPV_2d->GetYaxis()->SetTitle("y");
  eMPV_2d->Write();
  widthMPV_2d->SetTitle("widthMIP[ADC] map, "+dif);
  widthMPV_2d->GetXaxis()->SetTitle("x");
  widthMPV_2d->GetYaxis()->SetTitle("y");
  widthMPV_2d->Write();

  canvas_maps->cd(1);
  MPV_2d->SetStats(kFALSE);
  MPV_2d->SetTitle("MIP[ADC] map, "+dif);
  MPV_2d->GetXaxis()->SetTitle("x");
  MPV_2d->GetYaxis()->SetTitle("y");
  //MPV_2d->GetZaxis()->SetRangeUser(0,100);
  MPV_2d->Draw("colz");
  MPV_2d->Write();
  canvas_maps->cd(3);
  gPad->SetLogz();
  chi2NDF_2d->SetStats(kFALSE);
  chi2NDF_2d->SetTitle("chi2NDF map, "+dif);
  chi2NDF_2d->GetXaxis()->SetTitle("x");
  chi2NDF_2d->GetYaxis()->SetTitle("y");
  //chi2NDF_2d->GetZaxis()->SetRangeUser(0,20);
  chi2NDF_2d->Draw("colz");
  chi2NDF_2d->Write();
  canvas_maps->cd(4);
  gPad->SetLogz();
  mipEntries_2d->SetStats(kFALSE);
  mipEntries_2d->SetTitle("Hits map, "+dif);
  mipEntries_2d->GetXaxis()->SetTitle("x");
  mipEntries_2d->GetYaxis()->SetTitle("y");
  mipEntries_2d->Draw("colz");
  mipEntries_2d->Write();
  canvas_maps->cd(2);
  S_N_2d->SetStats(kFALSE);
  S_N_2d->SetTitle("S / N map, "+dif);
  S_N_2d->GetXaxis()->SetTitle("x");
  S_N_2d->GetYaxis()->SetTitle("y");
  //S_N_2d->GetZaxis()->SetRangeUser(0,50);
  S_N_2d->Draw("colz");
  S_N_2d->Write();
  
  // SUMMARY histograms (per chip)  
  for(int ichip=0; ichip<16; ichip ++) {
    TString dif0=dif+TString::Format("_chip%i",ichip);
    canvas_summary_mip->cd(ichip+1);
    MPV_chip[ichip]->SetStats(kFALSE);
    MPV_chip[ichip]->SetTitle("MIP[ADC]_map, "+dif0);
    MPV_chip[ichip]->GetXaxis()->SetTitle("x");
    MPV_chip[ichip]->GetYaxis()->SetTitle("y");
    MPV_chip[ichip]->GetZaxis()->SetRangeUser(0,100);
    //noise.at(ichip)->SetLineColor(2);
    MPV_chip[ichip]->Draw("colz");
    MPV_chip[ichip]->Write();
  }
  for(int ichip=0; ichip<16; ichip ++) {
    TString dif0=dif+TString::Format("_chip%i",ichip);
    canvas_summary_s_n->cd(ichip+1);
    S_N_chip[ichip]->SetStats(kFALSE);
    S_N_chip[ichip]->SetTitle("S_N_map, "+dif0);
    S_N_chip[ichip]->GetXaxis()->SetTitle("x");
    S_N_chip[ichip]->GetYaxis()->SetTitle("y");
    S_N_chip[ichip]->GetZaxis()->SetRangeUser(0,100);
    //noise.at(ichip)->SetLineColor(2);
    S_N_chip[ichip]->Draw("colz");
    S_N_chip[ichip]->Write();
  }

  // SUMMARY histograms (per chip)  
  canvas_summary->cd(1);
  MPV_slab->SetTitle("MIP, "+dif);
  MPV_slab->GetXaxis()->SetTitle("MIP[ADC] (ped. subtr)");
  MPV_slab->GetYaxis()->SetTitle("fitted channels");
  MPV_slab->Draw("h");
  MPV_slab->Write();
  canvas_summary->cd(2);
  S_N_slab->SetTitle("S / N, "+dif);
  S_N_slab->GetXaxis()->SetTitle("no units");
  S_N_slab->GetYaxis()->SetTitle("fitted channels");
  S_N_slab->Draw("h");
  S_N_slab->Write();
  canvas_summary->Print("results_mipcalibration/Signal_summary_"+dif+outputname+".png");


  //Correlations 
  canvas_correlations_S_N->cd(1);
  correl_S_N_vs_nhits->SetStats(kFALSE);
  correl_S_N_vs_nhits->SetTitle(dif);
  correl_S_N_vs_nhits->GetXaxis()->SetTitle("S / N");
  correl_S_N_vs_nhits->GetYaxis()->SetTitle("nhits");

  correl_S_N_vs_nhits->Draw("colz");
  correl_S_N_vs_nhits->Write();

  canvas_correlations_S_N->cd(2);
  correl_S_N_vs_mpv->SetStats(kFALSE);
  correl_S_N_vs_mpv->SetTitle(dif);
  correl_S_N_vs_mpv->GetXaxis()->SetTitle("S / N");
  correl_S_N_vs_mpv->GetYaxis()->SetTitle("MIP[ADC] (ped. subst.)");

  correl_S_N_vs_mpv->Draw("colz");
  correl_S_N_vs_mpv->Write();


  canvas_correlations_S_N->cd(3); 
  correl_S_N_vs_widthmpv->SetStats(kFALSE);
  correl_S_N_vs_widthmpv->SetTitle(dif);
  correl_S_N_vs_widthmpv->GetXaxis()->SetTitle("S / N");
  correl_S_N_vs_widthmpv->GetYaxis()->SetTitle("MIP width [ADC]");

  correl_S_N_vs_widthmpv->Draw("colz");
  correl_S_N_vs_widthmpv->Write();


  canvas_correlations_S_N->cd(4);
  correl_S_N_vs_widthped->SetStats(kFALSE);
  correl_S_N_vs_widthped->SetTitle(dif);
  correl_S_N_vs_widthped->GetXaxis()->SetTitle("S / N");
  correl_S_N_vs_widthped->GetYaxis()->SetTitle("Ped. width [ADC]");

  correl_S_N_vs_widthped->Draw("colz");
  correl_S_N_vs_widthped->Write();

  canvas_correlations_S_N->cd(5);
  correl_mpv_vs_nhits->SetStats(kFALSE);
  correl_mpv_vs_nhits->SetTitle(dif);
  correl_mpv_vs_nhits->GetXaxis()->SetTitle("MIP [ADC] (ped. subs.)");
  correl_mpv_vs_nhits->GetYaxis()->SetTitle("nhits");

  correl_mpv_vs_nhits->Draw("colz");
  correl_mpv_vs_nhits->Write();

  canvas_correlations_S_N->cd(6);
  correl_widthped_vs_nhits->SetStats(kFALSE);
  correl_widthped_vs_nhits->SetTitle(dif);
  correl_widthped_vs_nhits->GetXaxis()->SetTitle("Ped. width [ADC]");
  correl_widthped_vs_nhits->GetYaxis()->SetTitle("nhits");

  correl_widthped_vs_nhits->Draw("colz");
  correl_widthped_vs_nhits->Write();




  canvas_maps->Write();
  canvas_summary_mip->Write();
  canvas_summary_s_n->Write(); 
  canvas_summary->Write();
  canvas_correlations_S_N->Write();
  
  signalfile_summary->Close();

}



void singleSlabAnalysis::FindMasked(TString dif)
{

  //function that reads a root file and check which channels have or not signal (hit bit ==1).
  //should be used for root files that contain a full position scan, in order to provide meaninful list of masked channels.

  ofstream fout_masked("maskedchannels/masked_"+dif+".txt",ios::out);

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

void singleSlabAnalysis::ReadMap(TString filename) 
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

void singleSlabAnalysis::ReadMasked(TString filename) 
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


void singleSlabAnalysis::ReadPedestals(TString filename) 
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
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >>  tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst ;
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >>  tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst ;

  while(reading_file){
    reading_file >> tmp_chip >> tmp_channel >> tmp_ped[0] >> tmp_error[0] >> tmp_width[0] >> tmp_ped[1] >> tmp_error[1] >> tmp_width[1] >> tmp_ped[2] >> tmp_error[2] >> tmp_width[2] >> tmp_ped[3] >> tmp_error[3] >> tmp_width[3] >> tmp_ped[4] >> tmp_error[4] >> tmp_width[4] >> tmp_ped[5] >> tmp_error[5] >> tmp_width[5] >> tmp_ped[6] >> tmp_error[6] >> tmp_width[6] >> tmp_ped[7] >> tmp_error[7] >> tmp_width[7] >> tmp_ped[8] >> tmp_error[8] >> tmp_width[8] >> tmp_ped[9] >> tmp_error[9] >> tmp_width[9] >> tmp_ped[10] >> tmp_error[10] >> tmp_width[10] >> tmp_ped[11] >> tmp_error[11] >> tmp_width[11] >> tmp_ped[12] >> tmp_error[12] >> tmp_width[12] >> tmp_ped[13] >> tmp_error[13] >> tmp_width[13] >> tmp_ped[14] >> tmp_error[14] >> tmp_width[14];

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

void singleSlabAnalysis::PedestalAnalysis(TString dif,TString grid="",TString map_filename="/home/calice/SiWECAL-TB-analysis/fev10_chip_channel_x_y_mapping.txt")
{

  //Read the channel/chip -- x/y mapping
  ReadMap(map_filename);
  //Read the list of masked channels
  ReadMasked("maskedchannels/masked_"+dif+".txt");

  if(grid!="") dif=dif+grid;

  ofstream fout_ped("results_pedestal/Pedestal_"+dif+".txt",ios::out);

  fout_ped<<"#pedestal results (fit to a gaussian) remove channels/sca with two pedestals peaks from the analysis : "<<dif<<endl;
  fout_ped<<"#chip channel ped0 eped0 widthped0 ped1 eped1 widthped1... ped14 eped14 widhtped14 (all SCA)"<<endl;

  
  bool global = true;
  int maxnhit=5;

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

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


  TFile *pedfile = new TFile("results_pedestal/Pedestal_"+dif+".root" , "RECREATE");
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

	pedestal_entries_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , ped_sca.at(ichip).at(ichn).at(isca)->GetEntries());
	pedestal_tagged_entries_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , ped_sca_tagged.at(ichip).at(ichn).at(isca)->GetEntries());

	
	ped_mean.at(ichip).at(ichn).push_back(0.);
	ped_error.at(ichip).at(ichn).push_back(0.);
	ped_width.at(ichip).at(ichn).push_back(0.);

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
	    fout_ped<<f1->GetParameter(1) << " " << f1->GetParError(1)<<" "<<f1->GetParameter(2)<< " ";

	    ped_mean.at(ichip).at(ichn).at(isca)=f1->GetParameter(1);
	    ped_error.at(ichip).at(ichn).at(isca)=f1->GetParError(1);
	    ped_width.at(ichip).at(ichn).at(isca)=f1->GetParameter(2);

	    pedestal_chip.at(ichip)->Fill(f1->GetParameter(1));

	    pedestal_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParameter(1));
	    pedestal_width_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParameter(2));
	    pedestal_error_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParError(1));
	    pedestal_chi2ndf_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetChisquare() / f1->GetNDF());
	      
	  } else {
	    fout_ped<<0<< " " << 0<<" "<<0<<" ";
	  }
	} else {
	  fout_ped<<0<< " " << 0<<" "<<0<<" ";
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

  TFile *pedfile_summary = new TFile("results_pedestal/Pedestal_summary_"+dif+".root" , "RECREATE");
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
    pedestal_tagged_chip.at(ichip)->GetXaxis()->SetRangeUser(0,500);
    pedestal_tagged_chip.at(ichip)->SetTitle(TString::Format("Average pedestal_tagged, chip-%i",ichip));
    pedestal_tagged_chip.at(ichip)->GetXaxis()->SetTitle("ADC");
    pedestal_tagged_chip.at(ichip)->GetYaxis()->SetTitle("#");
    pedestal_tagged_chip.at(ichip)->Draw("hs");
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


void singleSlabAnalysis::PedestalAnalysis_gridpoints(TString dif,TString grid="",TString map_filename="../fev10_chip_channel_x_y_mapping.txt")
{

  //Read the channel/chip -- x/y mapping
  ReadMap(map_filename);
  //Read the list of masked channels
  ReadMasked("maskedchannels/masked_"+dif+".txt");

  if(grid!="") dif=dif+"_"+grid;

  ofstream fout_ped("results_pedestal/gridpoints/Pedestal_"+dif+".txt",ios::out);

  //fout_ped<<"#pedestal results (fit to a gaussian) remove channels/sca with two pedestals peaks from the analysis : "<<dif<<endl;
  //fout_ped<<"#chip channel ped0 eped0 widthped0 ped1 eped1 widthped1... ped14 eped14 widhtped14 (all SCA)"<<endl;

  
  bool global = true;
  int maxnhit=10;

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

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
	  if(charge_hiGain[ichip][isca][ichn]>10 && badbcid[ichip][isca]==0 && nhits[ichip][isca]<maxnhit+1 && corrected_bcid[ichip][isca]>1250 && corrected_bcid[ichip][isca]<2850 ) selection=true;
	  if(masked[ichip][ichn]==1) selection=false;
 	  if(gain_hit_high[ichip][isca][ichn]==0 && selection==true && gooddata==true)
	    ped_sca.at(ichip).at(ichn).at(isca)->Fill(charge_hiGain[ichip][isca][ichn]);

	  //bad events
	  selection=false;
	  if( ( badbcid[ichip][isca]>0 || nhits[ichip][isca]>maxnhit || gooddata==false) && (corrected_bcid[ichip][isca]>1250 && corrected_bcid[ichip][isca]<2850)  ) selection=true;
	  if(masked[ichip][ichn]==1) selection=false;
 	  if(gain_hit_high[ichip][isca][ichn]==0 && selection==true )
	    ped_sca_tagged.at(ichip).at(ichn).at(isca)->Fill(charge_hiGain[ichip][isca][ichn]);

	}
           
      }//isca

    }//ichip 
   
  }  // end first loop analysis to fill pedestal historgrams


  TFile *pedfile = new TFile("results_pedestal/gridpoints/Pedestal_"+dif+".root" , "RECREATE");
  pedfile->cd();


  //// PREANALYZE PEDESTALS
  //  decide which channels are being saved
  int savepedestal[16][64];
  for(int i=0; i<16; i++) {
    for(int j=0; j<64; j++) {
      savepedestal[i][j]=0;
    }
  }

 // do pedestal (chip/channel/sca based) analysis
  for(int ichip=0; ichip<16; ichip++) {
    for(int ichn=0; ichn<64; ichn++) {
      
      int ngoodscas=0;
      for(int isca=0; isca<15; isca++) {	

	if(ped_sca.at(ichip).at(ichn).at(isca)->GetEntries()> 300 ){ //max_entries/2 ) {
	  TSpectrum *s = new TSpectrum();
	  int npeaks = s->Search(ped_sca.at(ichip).at(ichn).at(isca),2,"",0.2); 

	  if(npeaks == 1) ngoodscas++;
	}

      }
      if(ngoodscas>12) savepedestal[ichip][ichn]=1;
    }
  }

  //////////////////////////////
  /// WRITE DOWN PEDESTALS if and only if pedestal is calculated for at least 14 SCA
  
  // do pedestal (chip/channel/sca based) analysis
  for(int ichip=0; ichip<16; ichip++) {
    for(int ichn=0; ichn<64; ichn++) {
      
      if(savepedestal[ichip][ichn]>0) {
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
	      fout_ped<<f1->GetParameter(1) << " " << f1->GetParError(1)<<" "<<f1->GetParameter(2)<< " ";
	      pedestal_chip.at(ichip)->Fill(f1->GetParameter(1));
	      
	      pedestal_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParameter(1));
	      pedestal_width_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParameter(2));
	      pedestal_error_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParError(1));
	      pedestal_chi2ndf_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetChisquare() / f1->GetNDF());
	      
	    } else {
	      fout_ped<<0<< " " << 0<<" "<<0<<" ";
	    }
	  } else {
	    fout_ped<<0<< " " << 0<<" "<<0<<" ";
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
	  

	}//end for sca
	fout_ped<<endl;
      }// end if savepedestal
    }
  }

  pedfile->Close();

  TFile *pedfile_summary = new TFile("results_pedestal/gridpoints/Pedestal_summary_"+dif+".root" , "RECREATE");
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
    pedestal_tagged_chip.at(ichip)->GetXaxis()->SetRangeUser(0,500);
    pedestal_tagged_chip.at(ichip)->SetTitle(TString::Format("Average pedestal_tagged, chip-%i",ichip));
    pedestal_tagged_chip.at(ichip)->GetXaxis()->SetTitle("ADC");
    pedestal_tagged_chip.at(ichip)->GetYaxis()->SetTitle("#");
    pedestal_tagged_chip.at(ichip)->Draw("hs");
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


void singleSlabAnalysis::PedestalAnalysis_bcid(TString dif, TString grid="",TString map_filename="../fev10_chip_channel_x_y_mapping.txt")
{

  //Read the channel/chip -- x/y mapping
  ReadMap(map_filename);
  //Read the list of masked channels
  ReadMasked("maskedchannels/masked_"+dif+".txt");

  if(grid!="") dif=dif+"_"+grid;

 
  bool global = true;
  int maxnhit=5;

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  // TH1F* pedestal_normal[16][16][10];


  // for(int ichip=0; ichip<16; ichip++) {
  //   for(int isca=0; isca<16; isca++) {
  //     for(int ibcid=0; ibcid<10; ibcid++) {
  // 	TString title = TString::Format("pedestal_normal_chip%i_sca%i_bcidrange%i",ichip,isca,ibcid);
  // 	if(isca==15) title = TString::Format("pedestal_normal_chip%i_sca%i_bcidrange%i",ichip,isca,ibcid);
  // 	pedestal_normal[ichip][isca][ibcid]= new TH1F(title,title,500,0.5,500.5);
  //     }
  //   }
  // }
  
  Float_t val_evt=1245;
  Float_t bcidrange[] = {0,50,100,200,400,800,1600,3200,6400};
  Float_t timerange[] = {0,50,100,200,400,800,1600,3200,6400};
  for(int ibcid=0; ibcid<9; ibcid++) timerange[ibcid]=0.4*timerange[ibcid];
  Int_t  binnum = sizeof(bcidrange)/sizeof(Float_t) - 1; // or just = 9

  Float_t bindev[102];
  for(int ibcid=0; ibcid<102; ibcid++) bindev[ibcid]=100*(-0.105+(0.2/101)*ibcid);

  TH2F* pedestal_deviation[16];

  for(int ichip=0; ichip<16; ichip++) {
    TString title  = TString::Format("pedestal_deviation_chip%i",ichip);
    pedestal_deviation[ichip]= new TH2F(title,title,binnum,timerange,101,bindev);//101,-0.0505,0.0505);
  }

  TString title  = TString::Format("pedestal_deviation_%s",dif.Data()); 
  TH2F* pedestal_deviation_slab = new TH2F(title,title,binnum,timerange,101,bindev);
  TH2F* badpedestal_chip_channel = new TH2F("badpedestal_chip_channel","badpedestal_chip_channel",64,-0.5,63.5,16,-0.5,15.5);
  TH2F* badpedestal_chip_sca = new TH2F("badpedestal_chip_sca","badpedestal_chip_sca",15,-0.5,14.5,16,-0.5,15.5);
  TH2F* badpedestal_sca_channel = new TH2F("badpedestal_sca_channel","badpedestal_sca_channel",64,-0.5,63.5,15,-0.5,14.5);


  // --------------------
  // all sca
  //  std::vector<TCanvas*> chip;
  std::vector<std::vector<std::vector<std::vector<TH1F*> > > > ped_sca_bcid;

  for(int ichip=0; ichip<16; ichip++) {
    std::vector<std::vector<std::vector<TH1F*> > >pedtemp_sca;
    for(int ichn=0; ichn<64; ichn++) {
      std::vector<std::vector<TH1F*> >pedtemp_sca2;
      for(int isca=0; isca<15; isca++) {
	std::vector<TH1F*> pedtemp_bcid;
	for(int ibcid=0; ibcid<8; ibcid++) {
	  TH1F *ped_bcid2 = new TH1F(TString::Format("ped_chip%i_chn%i_sca%i_bcidrange%i",ichip,ichn,isca,ibcid),TString::Format("ped_chip%i_chn%i_sca%i_bcidrange%i",ichip,ichn,isca,ibcid),1000,0.5,1000.5);
	  pedtemp_bcid.push_back(ped_bcid2);
	}
	pedtemp_sca2.push_back(pedtemp_bcid);
      }
      pedtemp_sca.push_back(pedtemp_sca2);
    }
    ped_sca_bcid.push_back(pedtemp_sca);
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
	  for( int ibcid=0; ibcid<8; ibcid++) {
	    bool selection=false;
	    
	    if(charge_hiGain[ichip][isca][ichn]>10 && badbcid[ichip][isca]==0 && nhits[ichip][isca]<maxnhit+1 && corrected_bcid[ichip][isca]>(val_evt+bcidrange[ibcid]) && corrected_bcid[ichip][isca]<(val_evt+bcidrange[ibcid+1]) ) selection=true;
	    if(masked[ichip][ichn]==1) selection=false;
	    if(gain_hit_high[ichip][isca][ichn]==0 && selection==true && gooddata==true) {
	      ped_sca_bcid.at(ichip).at(ichn).at(isca).at(ibcid)->Fill(charge_hiGain[ichip][isca][ichn]);
	    }
	  }
	}
           
      }//isca

    }//ichip 
   
  }  // end first loop analysis to fill pedestal historgrams



  
  // do pedestal (chip/channel/sca based) analysis
  for(int ichip=0; ichip<16; ichip++) {
    for(int isca=0; isca<15; isca++) {
      for(int ichn=0; ichn<64; ichn++) {
	  
	double mean=0;
	double error=0;

	TH1F *temp = new TH1F(TString::Format("temp_ped_chip%i_chn%i_sca%i",ichip,ichn,isca),TString::Format("temp_ped_chip%i_chn%i_sca%i",ichip,ichn,isca),1000,0.5,1000.5);

	for(int ibcid=0; ibcid<8; ibcid++) {
	
	  if(ped_sca_bcid.at(ichip).at(ichn).at(isca).at(ibcid)->GetEntries()> 1 )	    
	    temp->Add(ped_sca_bcid.at(ichip).at(ichn).at(isca).at(ibcid));

	}

	TSpectrum *s = new TSpectrum();
	int npeaks = s->Search(temp,2,"",0.05);
	
	
	if(npeaks == 1) {
	  Double_t *mean_peak=s->GetPositionX();
	  
	  TF1 *f0 = new TF1("f0","gaus",temp->GetMean()-3*temp->GetRMS(),temp->GetMean()+3*temp->GetRMS());
	  temp->Fit("f0","RQME");
	  
	  TF1 *f1 = new TF1("f1","gaus",f0->GetParameter(1)-2.*f0->GetParameter(2),f0->GetParameter(1)+2.*f0->GetParameter(2));
	  temp->Fit("f1","RQME");
	  
	  mean = f1->GetParameter(1);
	  error = f1->GetParError(1);
	}
	  
	
	if(mean<10) continue;
	
	for(int ibcid=0; ibcid<8; ibcid++) { 
	  
	  if(ped_sca_bcid.at(ichip).at(ichn).at(isca).at(ibcid)->GetEntries()> 400 ){ 
	    TSpectrum *s = new TSpectrum();
	    int npeaks = s->Search(ped_sca_bcid.at(ichip).at(ichn).at(isca).at(ibcid),2,"",0.05);
	    
	    
	    if(npeaks == 1) {
	      Double_t *mean_peak=s->GetPositionX();
	      
	      TF1 *f0 = new TF1("f0","gaus",ped_sca_bcid.at(ichip).at(ichn).at(isca).at(ibcid)->GetMean()-3*ped_sca_bcid.at(ichip).at(ichn).at(isca).at(ibcid)->GetRMS(),ped_sca_bcid.at(ichip).at(ichn).at(isca).at(ibcid)->GetMean()+3*ped_sca_bcid.at(ichip).at(ichn).at(isca).at(ibcid)->GetRMS());
	      ped_sca_bcid.at(ichip).at(ichn).at(isca).at(ibcid)->Fit("f0","RQME");
	      
	      TF1 *f1 = new TF1("f1","gaus",f0->GetParameter(1)-2.*f0->GetParameter(2),f0->GetParameter(1)+2.*f0->GetParameter(2));
	      ped_sca_bcid.at(ichip).at(ichn).at(isca).at(ibcid)->Fit("f1","RQME");

	      if(mean >0 && error >0) {
		double deviation = 100*(f1->GetParameter(1) - mean) / 70.;
		pedestal_deviation[ichip]->Fill( timerange[ibcid]+1. , deviation );
		pedestal_deviation_slab->Fill( timerange[ibcid]+1. , deviation );
		if(fabs(deviation)>2.) {
		  badpedestal_chip_channel->Fill(ichn,ichip,fabs(deviation));
		  badpedestal_chip_sca->Fill(isca,ichip,fabs(deviation));
		  badpedestal_sca_channel->Fill(ichn,isca,fabs(deviation));
		}
	      }
	    }
	  }

	}
      }
    }
  }
  
  
  TFile *pedfile_bcid = new TFile("results_pedestal/tests/Pedestal_bcid_"+dif+".root" , "RECREATE");
  pedfile_bcid->cd();

  // good pedestal events (not tagged events)
  TCanvas *canvas_pedestal_bcid = new TCanvas(TString::Format("pedestal_bcid_chip_%s",dif.Data()),TString::Format("pedestal_bcid_chip_%s",dif.Data()),1200,1200);
  canvas_pedestal_bcid->Divide(4,4);
  for(int ichip=0; ichip<16; ichip ++) {
    canvas_pedestal_bcid->cd(ichip+1);
    pedestal_deviation[ichip]->SetTitle(TString::Format("pedestal_deviation_chip%i",ichip));
    pedestal_deviation[ichip]->GetXaxis()->SetTitle("time after val_evt (us)");
    pedestal_deviation[ichip]->GetYaxis()->SetTitle(" % deviation [MIP]");
    pedestal_deviation[ichip]->Draw("colz");
  }
  canvas_pedestal_bcid->Write();

  // good pedestal events (not tagged events)                                                                                                      
  TCanvas *canvas_pedestal_bcid2 = new TCanvas(TString::Format("pedestal_bcid_%s",dif.Data()),TString::Format("pedestal_bcid_%s",dif.Data()),1200,1200);
  canvas_pedestal_bcid2->Divide(2,2);
  canvas_pedestal_bcid2->cd(1);
  pedestal_deviation_slab->SetTitle(TString::Format("pedestal_deviation_%s",dif.Data()));
  pedestal_deviation_slab->GetXaxis()->SetTitle("time after val_evt (us)");
  pedestal_deviation_slab->GetYaxis()->SetTitle(" % deviation [MIP]");
  pedestal_deviation_slab->Draw("colz");

  canvas_pedestal_bcid2->cd(2);
  badpedestal_chip_channel->SetStats(kFALSE);
  badpedestal_chip_channel->SetTitle(TString::Format("pedestal_deviation_map_%s",dif.Data()));
  badpedestal_chip_channel->GetXaxis()->SetTitle("channel");
  badpedestal_chip_channel->GetYaxis()->SetTitle("chip");
  badpedestal_chip_channel->Draw("colz");

  canvas_pedestal_bcid2->cd(3);
  badpedestal_sca_channel->SetStats(kFALSE);
  badpedestal_sca_channel->SetTitle(TString::Format("pedestal_deviation_map_%s",dif.Data()));
  badpedestal_sca_channel->GetXaxis()->SetTitle("channel");
  badpedestal_sca_channel->GetYaxis()->SetTitle("sca");
  badpedestal_sca_channel->Draw("colz");

  canvas_pedestal_bcid2->cd(4);
  badpedestal_chip_sca->SetStats(kFALSE);
  badpedestal_chip_sca->SetTitle(TString::Format("pedestal_deviation_map_%s",dif.Data()));
  badpedestal_chip_sca->GetXaxis()->SetTitle("sca");
  badpedestal_chip_sca->GetYaxis()->SetTitle("chip");
  badpedestal_chip_sca->Draw("colz");
  canvas_pedestal_bcid2->Write();

  for(int ichip=0; ichip<16; ichip ++) pedestal_deviation[ichip]->Write();
  pedestal_deviation_slab->Write();
  badpedestal_chip_channel->Write();
  badpedestal_sca_channel->Write();
  badpedestal_chip_sca->Write();
  
  pedfile_bcid->Close();

}


void singleSlabAnalysis::BcidCorrelations(TString dif)
{

  ReadMasked("maskedchannels/masked_"+dif+".txt");
  ReadMap("../fev10_chip_channel_x_y_mapping.txt");
  //function that reads a root file and check which channels have or not signal (hit bit ==1).
  //should be used for root files that contain a full position scan, in order to provide meaninful list of masked channels.
 
  bool global = true;
  int maxnhit=5;

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  // --------------------
  // all sca
  //  std::vector<TCanvas*> chip;
  TH2F* h_badevents_samechip = new TH2F("badevents_ifhit_samechip_dif_"+dif,"badevents_ifhit_samechip_dif_"+dif,32,-90,90,32,-90,90);
  TH1F* h_bcid_dif_total_samechip = new TH1F("bcid_samechip_dif_"+dif,"bcid_samechip_dif_"+dif,16002,-4000,4000);
  TH2F* h_badevents_diffchip = new TH2F("badevents_ifhit_diffchip_dif_"+dif,"badevents_ifhit_diffchip__dif_"+dif,32,-90,90,32,-90,90);
  TH1F* h_bcid_dif_total_diffchip = new TH1F("bcid_diffchip_dif_"+dif,"bcid_diffchip_dif_"+dif,8002,-4000,4000);


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
	    
	  double chip_bcid=ichip*5000.+bcid[ichip][isca];
	  
	  if(selection == true && gooddata==true && masked[ichip][ichn] == 0) {

	    for(int ichip2=0; ichip2<16; ichip2++) {
	      for(int isca2=0; isca2<15; isca2++) {
		double chip_bcid2=ichip2*5000.+bcid[ichip2][isca2];

		bool gooddata=true;
		if(global == true) {
		  int ntaggedasbad = 0;
		  for(int ichn2=0; ichn2<64; ichn2++) {
		    if(charge_hiGain[ichip2][isca2][ichn2]<10 && charge_hiGain[ichip2][isca2][ichn2]>-1 ) 
		      ntaggedasbad++;
		  }//ichn 
		  if ( ntaggedasbad > 0) gooddata=false;
		}

		for(int ichn2=0; ichn2<64; ichn2++) {

		  bool selection2=false;
		  if(gain_hit_high[ichip2][isca2][ichn2]==1 && charge_hiGain[ichip2][isca2][ichn2]>10 && corrected_bcid[ichip2][isca2]>1247 && corrected_bcid[ichip2][isca2]<2900 && masked[ichip2][ichn2] == 0)  selection2=true;

		  if( badbcid[ichip2][isca2]>1 || (badbcid[ichip2][isca2]>-1 && nhits[ichip2][isca2]>maxnhit ) ) selection2=selection2*true;
	  
		  if(selection2==true && ichip2==ichip && ichn2!=ichn && isca2!=isca) {
		    h_bcid_dif_total_samechip->Fill(corrected_bcid[ichip][isca]-corrected_bcid[ichip2][isca2]);//ichn2);
		    h_badevents_samechip->Fill(map_pointX[ichip2][ichn2],map_pointY[ichip2][ichn2]);
		  }

		  if(selection2==true && ichip2!=ichip && ichn2!=ichn && isca2!=isca ) {
		    h_bcid_dif_total_diffchip->Fill(bcid[ichip][isca]-bcid[ichip2][isca2]);
		    h_badevents_diffchip->Fill(map_pointX[ichip2][ichn2],map_pointY[ichip2][ichn2]);
		  }

		}
		
	      }
	    }
	  }
	}
           
      }//isca

    }//ichip 
  }  // end first loop analysis to fill pedestal historgrams

  
 
  TCanvas *canvas= new TCanvas("bcid_correl_"+dif, "bcid_correl_"+dif,1600,1000);
  canvas->Divide(2,2);
  canvas->cd(1);
  h_badevents_samechip->GetXaxis()->SetTitle("x");
  h_badevents_samechip->GetYaxis()->SetTitle("y");   
  h_badevents_samechip->Draw("colz");

  canvas->cd(2);
  h_bcid_dif_total_samechip->GetXaxis()->SetTitle("bcid_hits-bcid_badevents (chip_hit==chip_badevent)");
  h_bcid_dif_total_samechip->GetXaxis()->SetRangeUser(-40,40);
  h_bcid_dif_total_samechip->Draw("h");

  canvas->cd(3);
  h_badevents_diffchip->GetXaxis()->SetTitle("x");
  h_badevents_diffchip->GetYaxis()->SetTitle("y");   
  h_badevents_diffchip->Draw("colz");
  
  canvas->cd(4);
  h_bcid_dif_total_diffchip->GetXaxis()->SetTitle("bcid_hits-bcid_badevents (chip_hit!=chip_badevent)");
  h_bcid_dif_total_diffchip->GetXaxis()->SetRangeUser(-40,40);
  h_bcid_dif_total_diffchip->Draw("h");

  canvas->Print("results_pedestal/tests/bcid_correlations_"+dif+".root");


  
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



TF1* singleSlabAnalysis::langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
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

