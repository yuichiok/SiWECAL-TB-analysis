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


void singleSlabAnalysis::SignalAnalysis(TString slboard="testname", TString outputname="", TString map_filename="/home/calice/TB201906/tpecal/mapping/tb-2019/fev11_cob_chip_channel_x_y_mapping.txt", int maxnhit=5)
{

  //Read the channel/chip -- x/y mapping
  ReadMap(map_filename);
  //Read the list of pedestals (this information contains, implicitily, the masked channels information )
  if(outputname!="") slboard=slboard+outputname;
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
      TString histo_title=TString::Format("charge_%s_chip%i_chn%i",slboard.Data(),ichip,ichn);
      if(ichn==64) histo_title=TString::Format("charge_%s_chip%i",slboard.Data(),ichip);
      TH1F *chn_mip = new TH1F(histo_title,histo_title,4197,-100.5,4096.5);
      tmpmip_histo.push_back(chn_mip);
      TString s_n_title=TString::Format("s_n_%s_chip%i_chn%i",slboard.Data(),ichip,ichn);
      if(ichn==64) s_n_title=TString::Format("s_n_%s_chip%i",slboard.Data(),ichip);
      TH1F *chn_s_n = new TH1F(s_n_title,s_n_title,4197,-100.5,4096.5);
      tmps_n_histo.push_back(chn_s_n);
    }      
    mip_histo.push_back(tmpmip_histo);
    s_n_histo.push_back(tmps_n_histo);
  }

  //maps 
  TH2F* MPV_2d = new TH2F("MPV_2d"+slboard,"MPV_2d"+slboard,32,-90,90,32,-90,90);
  TH2F* eMPV_2d = new TH2F("eMPV_2d"+slboard,"eMPV_2d"+slboard,32,-90,90,32,-90,90);
  TH2F* widthMPV_2d = new TH2F("widthMPV_2d"+slboard,"widthMPV_2d"+slboard,32,-90,90,32,-90,90);
  TH2F* chi2NDF_2d = new TH2F("chi2NDF_2d"+slboard,"chi2NDF_2d"+slboard,32,-90,90,32,-90,90);
  TH2F* mipEntries_2d = new TH2F("mipEntries_2d"+slboard,"mipEntries_2d"+slboard,32,-90,90,32,-90,90);
  TH2F* S_N_2d = new TH2F("S_N_2d"+slboard,"S_N_2d"+slboard,32,-90,90,32,-90,90);

  TH2F* MPV_chipchn = new TH2F("MPV_chipchn"+slboard,"MPV_chipchn"+slboard,64,-0.5,63.5,16,-0.5,15.5);
  TH2F* eMPV_chipchn = new TH2F("eMPV_chipchn"+slboard,"eMPV_chipchn"+slboard,64,-0.5,63.5,16,-0.5,15.5);
  TH2F* widthMPV_chipchn = new TH2F("widthMPV_chipchn"+slboard,"widthMPV_chipchn"+slboard,64,-0.5,63.5,16,-0.5,15.5);
  TH2F* mipEntries_chipchn = new TH2F("mipEntries_chipchn"+slboard,"mipEntries_chipchn"+slboard,64,-0.5,63.5,16,-0.5,15.5);

  TH2F* correl_S_N_vs_nhits = new TH2F("correl_S_N_vs_nhits"+slboard,"correl_S_N_vs_nhits"+slboard,80,0.25,40.25,20,1225,51225);
  TH2F* correl_S_N_vs_mpv = new TH2F("correl_S_N_vs_mpv"+slboard,"correl_S_N_vs_mpv"+slboard,80,0.25,40.25,40,40.5,80.5);
  TH2F* correl_mpv_vs_nhits = new TH2F("correl_mpv_vs_nhits"+slboard,"correl_mpv_vs_nhits"+slboard,40,40.5,80.5,20,1225,51225);
  TH2F* correl_widthped_vs_nhits = new TH2F("correl_widthped_vs_nhits"+slboard,"correl_widthped_vs_nhits"+slboard,50,2.05,7.05,20,1225,51225);
  TH2F* correl_S_N_vs_widthmpv = new TH2F("correl_S_N_vs_widthmpv"+slboard,"correl_S_N_vs_widthmpv"+slboard,80,0.25,40.25,50,0.5,50.5);
  TH2F* correl_S_N_vs_widthped = new TH2F("correl_S_N_vs_widthped"+slboard,"correl_S_N_vs_widthped"+slboard,80,0.25,40.25,50,2.05,7.05);


  //summary plots
  TH1F* MPV_chip[16];
  TH1F* S_N_chip[16];
  for(int ichip=0; ichip<16; ichip++) {
    MPV_chip[ichip]= new TH1F(TString::Format("MPV_%s_chip%i",slboard.Data(),ichip),TString::Format("MPV_%s_chip%i",slboard.Data(),ichip),2000,0.05,200.05);
    S_N_chip[ichip]= new TH1F(TString::Format("S_N_%s_chip%i",slboard.Data(),ichip),TString::Format("S_N_%s_chip%i",slboard.Data(),ichip),500,0.05,50.05);
  }

  TH1F* MPV_slab= new TH1F(TString::Format("MPV_%s",slboard.Data()),TString::Format("MPV_%s",slboard.Data()),2000,0.05,200.05);
  TH1F* S_N_slab= new TH1F(TString::Format("S_N_%s",slboard.Data()),TString::Format("S_N_%s",slboard.Data()),1000,0.05,100.05);
 
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
	    if(charge_hiGain[ichip][isca][ichn]>0 && gain_hit_high[ichip][isca][ichn]==1 && badbcid[ichip][isca]==0 && nhits[ichip][isca]<maxnhit+1 ) {
	      if ( (slboard=="_SLB_0" || slboard=="_SLB_1" || slboard=="_SLB_2" || slboard=="_SLB_3")  &&
		   bcid[ichip][isca]>30 && (bcid[ichip][isca]<890 ||bcid[ichip][isca]>915)  ) selection=true;
	    } else {
	      if( (slboard!="_SLB_0" && slboard!="_SLB_1" && slboard!="_SLB_2" && slboard!="_SLB_3")  &&
		  bcid[ichip][isca]>800 ) selection=true;
	    }

	    if(selection==true) {
	      mip_histo.at(ichip).at(ichn)->Fill(charge_hiGain[ichip][isca][ichn]-ped_mean.at(ichip).at(ichn).at(isca));
	      s_n_histo.at(ichip).at(ichn)->Fill( (charge_hiGain[ichip][isca][ichn]-ped_mean.at(ichip).at(ichn).at(isca)) / ped_width.at(ichip).at(ichn).at(isca));

	      mip_histo.at(ichip).at(64)->Fill(charge_hiGain[ichip][isca][ichn]-ped_mean.at(ichip).at(ichn).at(isca));
	      s_n_histo.at(ichip).at(64)->Fill( (charge_hiGain[ichip][isca][ichn]-ped_mean.at(ichip).at(ichn).at(isca)) / ped_width.at(ichip).at(ichn).at(isca));
	    }
	  }
    
    	}//ichn
      }//sca
    }//ichip
    
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
            
      if(mip_histo.at(ichip).at(ichn)->GetEntries()>50){

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
	  MPV_2d->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn],mpv);
	  eMPV_2d->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn],empv);
	  widthMPV_2d->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn],wmpv);
	  chi2NDF_2d->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn],chi2ndf);
	  mipEntries_2d->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn],mipentries);
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
	}	
      } else {
	if(ichn<64) fout_mip<<ichip<<" "<<ichn<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<endl;

      }
    }
  }


  signalfile_summary->cd();

  TCanvas *canvas_maps = new TCanvas(TString::Format("Signal_%s",slboard.Data()),TString::Format("Signal_%s",slboard.Data()),1200,800);
  canvas_maps->Divide(2,2);
  TCanvas *canvas_summary_mip = new TCanvas(TString::Format("mip_perchip_%s",slboard.Data()),TString::Format("mip_perchip_%s",slboard.Data()),1200,800);
  canvas_summary_mip->Divide(4,4);
  TCanvas *canvas_summary_s_n = new TCanvas(TString::Format("s_n_perchip_%s",slboard.Data()),TString::Format("s_n_perchip_%s",slboard.Data()),1200,1200);
  canvas_summary_s_n->Divide(4,4); 
  TCanvas *canvas_summary = new TCanvas(TString::Format("summary_%s",slboard.Data()),TString::Format("summary_%s",slboard.Data()),1200,600);
  canvas_summary->Divide(2,1);

  TCanvas *canvas_correlations_S_N = new TCanvas(TString::Format("correlations_S_N_%s",slboard.Data()),TString::Format("correlations_S_N_%s",slboard.Data()),800,900);
  canvas_correlations_S_N->Divide(2,3);

  
  // SUMMARY MAPS
  eMPV_2d->SetTitle("eMIP[ADC] map, "+slboard);
  eMPV_2d->GetXaxis()->SetTitle("x");
  eMPV_2d->GetYaxis()->SetTitle("y");
  eMPV_2d->Write();
  widthMPV_2d->SetTitle("widthMIP[ADC] map, "+slboard);
  widthMPV_2d->GetXaxis()->SetTitle("x");
  widthMPV_2d->GetYaxis()->SetTitle("y");
  widthMPV_2d->Write();

  canvas_maps->cd(1);
  MPV_2d->SetStats(kFALSE);
  MPV_2d->SetTitle("MIP[ADC] map, "+slboard);
  MPV_2d->GetXaxis()->SetTitle("x");
  MPV_2d->GetYaxis()->SetTitle("y");
  //MPV_2d->GetZaxis()->SetRangeUser(0,100);
  MPV_2d->Draw("colz");
  MPV_2d->Write();
  canvas_maps->cd(3);
  gPad->SetLogz();
  chi2NDF_2d->SetStats(kFALSE);
  chi2NDF_2d->SetTitle("chi2NDF map, "+slboard);
  chi2NDF_2d->GetXaxis()->SetTitle("x");
  chi2NDF_2d->GetYaxis()->SetTitle("y");
  //chi2NDF_2d->GetZaxis()->SetRangeUser(0,20);
  chi2NDF_2d->Draw("colz");
  chi2NDF_2d->Write();
  canvas_maps->cd(4);
  gPad->SetLogz();
  mipEntries_2d->SetStats(kFALSE);
  mipEntries_2d->SetTitle("Hits map, "+slboard);
  mipEntries_2d->GetXaxis()->SetTitle("x");
  mipEntries_2d->GetYaxis()->SetTitle("y");
  mipEntries_2d->Draw("colz");
  mipEntries_2d->Write();
  canvas_maps->cd(2);
  S_N_2d->SetStats(kFALSE);
  S_N_2d->SetTitle("S / N map, "+slboard);
  S_N_2d->GetXaxis()->SetTitle("x");
  S_N_2d->GetYaxis()->SetTitle("y");
  //S_N_2d->GetZaxis()->SetRangeUser(0,50);
  S_N_2d->Draw("colz");
  S_N_2d->Write();

  MPV_chipchn->SetTitle("MIP[ADC] map, "+slboard);
  MPV_chipchn->GetXaxis()->SetTitle("chip");
  MPV_chipchn->GetYaxis()->SetTitle("chn");
  MPV_chipchn->Write();
  widthMPV_chipchn->SetTitle("widthMIP[ADC] map, "+slboard);
  widthMPV_chipchn->GetXaxis()->SetTitle("chip");
  widthMPV_chipchn->GetYaxis()->SetTitle("chn");
  widthMPV_chipchn->Write();
  eMPV_chipchn->SetTitle("eMIP[ADC] map, "+slboard);
  eMPV_chipchn->GetXaxis()->SetTitle("chip");
  eMPV_chipchn->GetYaxis()->SetTitle("chn");
  eMPV_chipchn->Write();
  mipEntries_chipchn->SetTitle("widthMIP[ADC] map, "+slboard);
  mipEntries_chipchn->GetXaxis()->SetTitle("chip");
  mipEntries_chipchn->GetYaxis()->SetTitle("chn");
  mipEntries_chipchn->Write();


  
  // SUMMARY histograms (per chip)  
  for(int ichip=0; ichip<16; ichip ++) {
    TString slboard0=slboard+TString::Format("_chip%i",ichip);
    canvas_summary_mip->cd(ichip+1);
    MPV_chip[ichip]->SetStats(kFALSE);
    MPV_chip[ichip]->SetTitle("MIP[ADC]_map, "+slboard0);
    MPV_chip[ichip]->GetXaxis()->SetTitle("x");
    MPV_chip[ichip]->GetYaxis()->SetTitle("y");
    MPV_chip[ichip]->GetZaxis()->SetRangeUser(0,100);
    //noise.at(ichip)->SetLineColor(2);
    MPV_chip[ichip]->Draw("colz");
    MPV_chip[ichip]->Write();
  }
  for(int ichip=0; ichip<16; ichip ++) {
    TString slboard0=slboard+TString::Format("_chip%i",ichip);
    canvas_summary_s_n->cd(ichip+1);
    S_N_chip[ichip]->SetStats(kFALSE);
    S_N_chip[ichip]->SetTitle("S_N_map, "+slboard0);
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
  
  signalfile_summary->Close();

}


void singleSlabAnalysis::FindMasked(TString slboard)
{

  //function that reads a root file and check which channels have or not signal (hit bit ==1).
  //should be used for root files that contain a full position scan, in order to provide meaninful list of masked channels.

  ofstream fout_masked("maskedchannels/masked"+slboard+".txt",ios::out);

  fout_masked<<"#list of masked channels, per slboard: "<<slboard<<endl;
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
	for(int ichn=0; ichn<64; ichn++) {

	  bool selection=false;
	  if( badbcid[ichip][isca]==0 && nhits[ichip][isca]<maxnhit+1 ) {
	    if ( (slboard=="_SLB_0" || slboard=="_SLB_1" || slboard=="_SLB_2" || slboard=="_SLB_3")  &&
		 bcid[ichip][isca]>30 && (bcid[ichip][isca]<890 ||bcid[ichip][isca]>915)  ) selection=true;
	  } else {
	    if( (slboard!="_SLB_0" && slboard!="_SLB_1" && slboard!="_SLB_2" && slboard!="_SLB_3")  &&
		bcid[ichip][isca]>800 ) selection=true;
	  }
	  if(charge_hiGain[ichip][isca][ichn]>30 && gain_hit_high[ichip][isca][ichn]==1 && selection==true)
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
      if(h_charge_channel.at(ichip).at(ichn)->GetEntries()>10) total_channels_chip[ichip]++;

    }
  }

  // the real estimation of masked channels is done by requiring a 0.1 of the estimated Nhits per channel for the whole data taking period run (per chip)
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
	//	if(h_charge_channel.at(ichip).at(ichn)->GetEntries()>0.1*(total_entries_chip[ichip]/total_channels_chip[ichip])) fout_masked<<"0"<<endl;
	if(h_charge_channel.at(ichip).at(ichn)->GetEntries()>0) fout_masked<<"0"<<endl;
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

void singleSlabAnalysis::PedestalAnalysis(TString slboard,TString sufix="",TString map_filename="/home/calice/TB201906/tpecal/mapping/tb-2019/fev11_cob_chip_channel_x_y_mapping.txt", int maxnhit=5)
{

  //Read the channel/chip -- x/y mapping
  ReadMap(map_filename);
  //Read the list of masked channels
  ReadMasked("maskedchannels/masked"+slboard+".txt");

  if(sufix!="") slboard=slboard+sufix;

  ofstream fout_ped("results_pedestal/Pedestal"+slboard+".txt",ios::out);

  fout_ped<<"#pedestal results (fit to a gaussian) remove channels/sca with two pedestals peaks from the analysis : "<<slboard<<endl;
  fout_ped<<"#chip channel ped0 eped0 widthped0 ped1 eped1 widthped1... ped14 eped14 widhtped14 (all SCA)"<<endl;

  
  bool global = true;
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
	    if(gain_hit_high[ichip][isca][ichn]==1 && charge_hiGain[ichip][isca][ichn]<15 && charge_hiGain[ichip][isca][ichn]>-1 ) 
	      ntaggedasbad++;
	  }//ichn 
	  if ( ntaggedasbad > 0) gooddata=false;
	}
	
	for(int ichn=0; ichn<64; ichn++) {

	  //good events
	  bool selection=false;

	  if( charge_lowGain[ichip][isca][ichn]>30 && badbcid[ichip][isca]==0 && nhits[ichip][isca]<maxnhit+1 ) {
	    if ( (slboard=="_SLB_0" || slboard=="_SLB_1" || slboard=="_SLB_2" || slboard=="_SLB_3")  &&
		 bcid[ichip][isca]>30 && (bcid[ichip][isca]<890 ||bcid[ichip][isca]>915)  ) selection=true;
	  } else {
	    if( (slboard!="_SLB_0" && slboard!="_SLB_1" && slboard!="_SLB_2" && slboard!="_SLB_3")  &&
		bcid[ichip][isca]>800 ) selection=true;
	  }
		  
	  // if(masked[ichip][ichn]==1) selection=false;
 	  if(gain_hit_low[ichip][isca][ichn]==0 && selection==true && gooddata==true)
	    ped_sca.at(ichip).at(ichn).at(isca)->Fill(charge_lowGain[ichip][isca][ichn]);

	  //bad events
	  //	  selection=false;
	  //	  if( ( badbcid[ichip][isca]!=0  || nhits[ichip][isca]>maxnhit || gooddata==false) ) selection=true;
	  //	  if(masked[ichip][ichn]==1) selection=false;
 	  if(gain_hit_low[ichip][isca][ichn]==0 && (selection==false || gooddata==false) )
	    ped_sca_tagged.at(ichip).at(ichn).at(isca)->Fill(charge_lowGain[ichip][isca][ichn]);

	}
           
      }//isca

    }//ichip 
   
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

	pedestal_entries_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , ped_sca.at(ichip).at(ichn).at(isca)->GetEntries());
	pedestal_tagged_entries_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , ped_sca_tagged.at(ichip).at(ichn).at(isca)->GetEntries());

	
	ped_mean.at(ichip).at(ichn).push_back(0.);
	ped_error.at(ichip).at(ichn).push_back(0.);
	ped_width.at(ichip).at(ichn).push_back(0.);

	if(ped_sca.at(ichip).at(ichn).at(isca)->GetEntries()> 50 ){ //max_entries/2 ) {
	  TSpectrum *s = new TSpectrum();
	  int npeaks = s->Search(ped_sca.at(ichip).at(ichn).at(isca),2,"",0.8); 
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
	      if(ipeak != npeak_max) pedestal_slboard_chip.at(ichip)->Fill( mean_peak[npeak_max] - mean_peak[ipeak]);
	    }
	    if(npeaks ==1 ) {
	      Double_t *mean_peak=s->GetPositionX();
	      mean_peak[0]=mean_peak_higher;
	      
	      TF1 *f0 = new TF1("f0","gaus",mean_peak[0]-ped_sca.at(ichip).at(ichn).at(isca)->GetRMS(),mean_peak[0]+ped_sca.at(ichip).at(ichn).at(isca)->GetRMS());
	      ped_sca.at(ichip).at(ichn).at(isca)->Fit("f0","RQNOC");
	      
	      TF1 *f1 = new TF1("f1","gaus",f0->GetParameter(1)-f0->GetParameter(2),f0->GetParameter(1)+f0->GetParameter(2));
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
	      fout_ped<<ped_sca.at(ichip).at(ichn).at(isca)->GetMean()<< " " << ped_sca.at(ichip).at(ichn).at(isca)->GetRMS()/sqrt(ped_sca.at(ichip).at(ichn).at(isca)->GetEntries())<<" "<<ped_sca.at(ichip).at(ichn).at(isca)->GetRMS()<<" ";
	      pedestal_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , ped_sca.at(ichip).at(ichn).at(isca)->GetMean());
	      pedestal_width_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , ped_sca.at(ichip).at(ichn).at(isca)->GetRMS() );
	      pedestal_error_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , ped_sca.at(ichip).at(ichn).at(isca)->GetRMS()/sqrt(ped_sca.at(ichip).at(ichn).at(isca)->GetMean()));
	      pedestal_chi2ndf_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , 0);
	    }
	    
	  } else fout_ped<<0<< " " << 0<<" "<<0<<" ";
	  
	}else {
	  fout_ped<<0<< " " << 0<<" "<<0<<" ";
	}
     

	// analyze pedestal for tagged events
	if(ped_sca_tagged.at(ichip).at(ichn).at(isca)->GetEntries()> 250) {
	  TSpectrum *s = new TSpectrum();
	  int npeaks = s->Search(ped_sca_tagged.at(ichip).at(ichn).at(isca),2,"",0.2);
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

  TFile *pedfile_summary = new TFile("results_pedestal/Pedestal_summary"+slboard+".root" , "RECREATE");
  pedfile_summary->cd();
  
  // good pedestal events (not tagged events)
  TCanvas *canvas_pedestal_map = new TCanvas("pedestal_map"+slboard, "pedestal_map"+slboard,1200,1200);
  canvas_pedestal_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString slboard0=slboard+TString::Format("_sca%i",isca);
    canvas_pedestal_map->cd(isca+1);
    pedestal_map[isca]->SetStats(kFALSE);
    pedestal_map[isca]->SetTitle("pedestal_map, "+slboard0);
    pedestal_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_map[isca]->GetZaxis()->SetRangeUser(200,400);

    //noise.at(ichip)->SetLineColor(2);
    pedestal_map[isca]->Draw("colz");
    pedestal_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_width_map = new TCanvas("pedestal_width_map"+slboard, "pedestal_width_map"+slboard,1200,1200);
  canvas_pedestal_width_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString slboard0=slboard+TString::Format("_sca%i",isca);
    canvas_pedestal_width_map->cd(isca+1);
    pedestal_width_map[isca]->SetStats(kFALSE);
    pedestal_width_map[isca]->SetTitle("pedestal_width_map, "+slboard0);
    pedestal_width_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_width_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_width_map[isca]->GetZaxis()->SetRangeUser(1,20);
    //noise.at(ichip)->SetLineColor(2);
    pedestal_width_map[isca]->Draw("colz");
    pedestal_width_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_error_map = new TCanvas("pedestal_error_map"+slboard, "pedestal_error_map"+slboard,1200,1200);
  canvas_pedestal_error_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString slboard0=slboard+TString::Format("_sca%i",isca);
    canvas_pedestal_error_map->cd(isca+1);
    pedestal_error_map[isca]->SetStats(kFALSE);
    pedestal_error_map[isca]->SetTitle("pedestal_error_map, "+slboard0);
    pedestal_error_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_error_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_error_map[isca]->Draw("colz");
    pedestal_error_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_npeaks_map = new TCanvas("pedestal_npeaks_map"+slboard, "pedestal_npeaks_map"+slboard,1200,1200);
  canvas_pedestal_npeaks_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString slboard0=slboard+TString::Format("_sca%i",isca);
    canvas_pedestal_npeaks_map->cd(isca+1);
    pedestal_npeaks_map[isca]->SetStats(kFALSE);
    pedestal_npeaks_map[isca]->SetTitle("pedestal_npeaks_map, "+slboard0);
    pedestal_npeaks_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_npeaks_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_npeaks_map[isca]->GetZaxis()->SetRangeUser(0.5,4.5);
    //noise.at(ichip)->SetLineColor(2);
    pedestal_npeaks_map[isca]->Draw("colz");
    pedestal_npeaks_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_chi2ndf_map = new TCanvas("pedestal_chi2ndf_map"+slboard, "pedestal_chi2ndf_map"+slboard,1200,1200);
  canvas_pedestal_chi2ndf_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString slboard0=slboard+TString::Format("_sca%i",isca);
    canvas_pedestal_chi2ndf_map->cd(isca+1);
    pedestal_chi2ndf_map[isca]->SetStats(kFALSE);
    pedestal_chi2ndf_map[isca]->SetTitle("pedestal_chi2ndf_map, "+slboard0);
    pedestal_chi2ndf_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_chi2ndf_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_chi2ndf_map[isca]->Draw("colz");
    pedestal_chi2ndf_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_entries_map = new TCanvas("pedestal_entries_map"+slboard, "pedestal_entries_map"+slboard,1200,1200);
  canvas_pedestal_entries_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString slboard0=slboard+TString::Format("_sca%i",isca);
    canvas_pedestal_entries_map->cd(isca+1);
    pedestal_entries_map[isca]->SetStats(kFALSE);
    pedestal_entries_map[isca]->SetTitle("pedestal_entries_map, "+slboard0);
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

   
  TCanvas *canvas_pedestal_pedestal = new TCanvas("pedestal_average"+slboard, "pedestal_average"+slboard,1200,1200);
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

  TCanvas *canvas_pedestal_slboard = new TCanvas("pedestal_slboard"+slboard, "pedestal_slboard"+slboard,1200,1200);
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
  TCanvas *canvas_pedestal_tagged_map = new TCanvas("pedestal_tagged_map"+slboard, "pedestal_tagged_map"+slboard,1200,1200);
  canvas_pedestal_tagged_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString slboard0=slboard+TString::Format("_sca%i",isca);
    canvas_pedestal_tagged_map->cd(isca+1);
    pedestal_tagged_map[isca]->SetStats(kFALSE);
    pedestal_tagged_map[isca]->SetTitle("pedestal_tagged_map, "+slboard0);
    pedestal_tagged_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_tagged_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_tagged_map[isca]->GetZaxis()->SetRangeUser(200,400);

    //noise.at(ichip)->SetLineColor(2);
    pedestal_tagged_map[isca]->Draw("colz");
    pedestal_tagged_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_tagged_width_map = new TCanvas("pedestal_tagged_width_map"+slboard, "pedestal_tagged_width_map"+slboard,1200,1200);
  canvas_pedestal_tagged_width_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString slboard0=slboard+TString::Format("_sca%i",isca);
    canvas_pedestal_tagged_width_map->cd(isca+1);
    pedestal_tagged_width_map[isca]->SetStats(kFALSE);
    pedestal_tagged_width_map[isca]->SetTitle("pedestal_tagged_width_map, "+slboard0);
    pedestal_tagged_width_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_tagged_width_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_tagged_width_map[isca]->GetZaxis()->SetRangeUser(0,7);
    //noise.at(ichip)->SetLineColor(2);
    pedestal_tagged_width_map[isca]->Draw("colz");
    pedestal_tagged_width_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_tagged_error_map = new TCanvas("pedestal_tagged_error_map"+slboard, "pedestal_tagged_error_map"+slboard,1200,1200);
  canvas_pedestal_tagged_error_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString slboard0=slboard+TString::Format("_sca%i",isca);
    canvas_pedestal_tagged_error_map->cd(isca+1);
    pedestal_tagged_error_map[isca]->SetStats(kFALSE);
    pedestal_tagged_error_map[isca]->SetTitle("pedestal_tagged_error_map, "+slboard0);
    pedestal_tagged_error_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_tagged_error_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_tagged_error_map[isca]->Draw("colz");
    pedestal_tagged_error_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_tagged_npeaks_map = new TCanvas("pedestal_tagged_npeaks_map"+slboard, "pedestal_tagged_npeaks_map"+slboard,1200,1200);
  canvas_pedestal_tagged_npeaks_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString slboard0=slboard+TString::Format("_sca%i",isca);
    canvas_pedestal_tagged_npeaks_map->cd(isca+1);
    pedestal_tagged_npeaks_map[isca]->SetStats(kFALSE);
    pedestal_tagged_npeaks_map[isca]->SetTitle("pedestal_tagged_npeaks_map, "+slboard0);
    pedestal_tagged_npeaks_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_tagged_npeaks_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_tagged_npeaks_map[isca]->GetZaxis()->SetRangeUser(0.5,4.5);
    //noise.at(ichip)->SetLineColor(2);
    pedestal_tagged_npeaks_map[isca]->Draw("colz");
    pedestal_tagged_npeaks_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_tagged_chi2ndf_map = new TCanvas("pedestal_tagged_chi2ndf_map"+slboard, "pedestal_tagged_chi2ndf_map"+slboard,1200,1200);
  canvas_pedestal_tagged_chi2ndf_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString slboard0=slboard+TString::Format("_sca%i",isca);
    canvas_pedestal_tagged_chi2ndf_map->cd(isca+1);
    pedestal_tagged_chi2ndf_map[isca]->SetStats(kFALSE);
    pedestal_tagged_chi2ndf_map[isca]->SetTitle("pedestal_tagged_chi2ndf_map, "+slboard0);
    pedestal_tagged_chi2ndf_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_tagged_chi2ndf_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_tagged_chi2ndf_map[isca]->Draw("colz");
    pedestal_tagged_chi2ndf_map[isca]->Write();
  }

  TCanvas *canvas_pedestal_tagged_entries_map = new TCanvas("pedestal_tagged_entries_map"+slboard, "pedestal_tagged_entries_map"+slboard,1200,1200);
  canvas_pedestal_tagged_entries_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString slboard0=slboard+TString::Format("_sca%i",isca);
    canvas_pedestal_tagged_entries_map->cd(isca+1);
    pedestal_tagged_entries_map[isca]->SetStats(kFALSE);
    pedestal_tagged_entries_map[isca]->SetTitle("pedestal_tagged_entries_map, "+slboard0);
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

   
  TCanvas *canvas_pedestal_tagged_pedestal = new TCanvas("pedestal_tagged_average"+slboard, "pedestal_tagged_average"+slboard,1200,1200);
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


  TCanvas *canvas_pedestal_tagged_slboard = new TCanvas("pedestal_tagged_slboard"+slboard, "pedestal_tagged_slboard"+slboard,1200,1200);
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
 
void singleSlabAnalysis::Retriggers(TString slboard, TString sufix="",TString mapfile="/home/calice/TB201906/tpecal/mapping/tb-2019/fev11_cob_chip_channel_x_y_mapping.txt", int maxnhit=10)
{

  ReadMap(mapfile);
  //function that reads a root file and check which channels have or not signal (hit bit ==1).
  //should be used for root files that contain a full position scan, in order to provide meaninful list of masked channels.
 
  bool global = true;

  if(sufix!="") slboard=slboard+sufix;

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();

  //summary
  TH1F* h_total_triggers = new TH1F("total_triggers"+slboard,"total_triggers"+slboard,16,-0.5,15.5);
  TH1F* h_total_retriggers = new TH1F("total_retriggers"+slboard,"total_retriggers"+slboard,16,-0.5,15.5);
  TH1F* h_total_retriggers_trains = new TH1F("total_retriggers_trains"+slboard,"total_retriggers_trains"+slboard,16,-0.5,15.5);
  TH2F* h_total_triggers_2d = new TH2F("total_triggers_2d"+slboard,"total_triggers_2d"+slboard,16,-0.5,15.5,3,-0.5,2.5);
  TH1F* h_signal = new TH1F("signal"+slboard,"signal"+slboard,4000,50.5,4050.5);

  TH2F* h_total_planeE = new TH2F("total_planeE"+slboard,"total_planeE"+slboard,16,-0.5,15.5,15,-0.5,14.5);
  TH2F* h_total_negE = new TH2F("total_negE"+slboard,"total_negE"+slboard,16,-0.5,15.5,15,-0.5,14.5);
  TH2F* h_total_negE_nohit = new TH2F("total_negE_nohit"+slboard,"total_negE_nohit"+slboard,16,-0.5,15.5,15,-0.5,14.5);

  
  // --------------------
  // per sca
  TH2F* h_first_retriggering = new TH2F("first_retriggering"+slboard,"first_retriggering"+slboard,16,-0.5,15.5,64,-0.5,63.5);
  TH2F* h_all_retriggering = new TH2F("all_retriggering"+slboard,"all_retriggering"+slboard,16,-0.5,15.5,64,-0.5,63.5);

  TH2F* h_retrigger_train[16];
  TH1F* h_retrigger_train_length[16];
  TH1F* h_retrigger_n_trains[16];
  TH1F* h_signal_n_events[16];

  TH2F* h_bcid_correlation[16];
  TH2F* h_bcid2_correlation[16];

  TH1F* good[16];
  TH1F* bad[16];
  TH2F* good_vs_bad[16];

  for(int ichip=0; ichip<16; ichip++) {
    h_retrigger_train[ichip] = new TH2F(TString::Format("retrigger_train_chip%i_",ichip)+slboard,TString::Format("retrigger_train_chip%i_",ichip)+slboard,15,-0.5,14.5,64,-0.5,63.5);
    h_retrigger_train_length[ichip] = new TH1F(TString::Format("retrigger_train_length_chip%i_",ichip)+slboard,TString::Format("retrigger_train_length_chip%i_",ichip)+slboard,15,0.5,15.5);
    h_retrigger_n_trains[ichip] = new TH1F(TString::Format("retrigger_n_trains_chip%i_",ichip)+slboard,TString::Format("retrigger_n_trains_chip%i_",ichip)+slboard,5,-0.5,4.5);
    h_signal_n_events[ichip] = new TH1F(TString::Format("signal_n_events_chip%i_",ichip)+slboard,TString::Format("signal_n_events_chip%i_",ichip)+slboard,5,-0.5,4.5);

    h_bcid_correlation[ichip]
      = new TH2F(TString::Format("bcid_rettriggers_correlation_chip%i_",ichip)+slboard,TString::Format("bcid_rettriggers_correlation_chip%i_",ichip)+slboard,50,0,5000,1000,0.5,1000.5);

    h_bcid2_correlation[ichip]
      = new TH2F(TString::Format("bcid_triggers_correlation_chip%i_",ichip)+slboard,TString::Format("bcid_triggers_correlation_chip%i_",ichip)+slboard,500,5,5005,1000,0.5,1000.5);
    
    good[ichip]= new TH1F(TString::Format("good_chip%i_%s",ichip,slboard.Data()),TString::Format("good_chip%i_%s",ichip,slboard.Data()),15,0.5,15.5);
    bad[ichip]= new TH1F(TString::Format("bad_chip%i_%s",ichip,slboard.Data()),TString::Format("bad_chip%i_%s",ichip,slboard.Data()),15,0.5,15.5);
    good_vs_bad[ichip]= new TH2F(TString::Format("good_vs_bad_chip%i_%s",ichip,slboard.Data()),TString::Format("good_vs_bad_chip%i_%s",ichip,slboard.Data()),15,0.5,15.5,15,0.5,15.5);
  }

  double n_total_triggers[16];
  double n_total_retriggers[16];
  double n_total_retriggers_trains[16];

  for(int ichip=0; ichip<16; ichip++) {
    
    n_total_triggers[ichip]=0;
    n_total_retriggers[ichip]=0;
    n_total_retriggers_trains[ichip]=0;

  }

  // -----------------------------------------------------------------------------------------------------   
  // SCA analysis
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    for(int ichip=0; ichip<16; ichip++) {

      bool retrig=false;
      for(int isca=0; isca<15; isca++) {
	if(badbcid[ichip][isca]!=0 && retrig==false && bcid[ichip][isca]>20 && (bcid[ichip][isca]<890 || bcid[ichip][isca]>915)) retrig=true;
	if(badbcid[ichip][isca]!=0 && retrig==true && bcid[ichip][isca]>20 && (bcid[ichip][isca]<890 || bcid[ichip][isca]>915) ) {
	  h_bcid2_correlation[ichip]->Fill(bcid[ichip][isca-1],bcid[ichip][isca]-bcid[ichip][isca-1]);
	}

	int hits_plane=0;
	int hits_negative=0;
	int nohits_negative=0;

	for(int ichn=0;ichn<64;ichn++) {
	  if(gain_hit_high[ichip][isca][ichn]==1) hits_plane++;
	  if(gain_hit_high[ichip][isca][ichn]==1 && charge_hiGain[ichip][isca][ichn]<100) hits_negative++;
	  if(gain_hit_high[ichip][isca][ichn]==0 && charge_hiGain[ichip][isca][ichn]<100) nohits_negative++;
	}

	if(hits_plane>(maxnhit-1)) h_total_planeE->Fill(ichip,isca);
	if(hits_negative>(maxnhit-1)) h_total_negE->Fill(ichip,isca);
	if(nohits_negative>(maxnhit-1)) h_total_negE_nohit->Fill(ichip,isca);

      }
    }
    
    for(int ichip=0; ichip<16; ichip++) {

      double good_triggers=0;
      double bad_triggers=0;

      bool first_retrig=false;
      int first_sca_retrig=0;
      int length=0;
      int n_trains=0;
      int n_events=0;

      float bcid_signal=0;
      float bcid_retrigger=0;
            
      for(int isca=0; isca<15; isca++) {

	if(badbcid[ichip][isca]==0) {
	  bcid_signal=bcid[ichip][isca];
	  good_triggers++;
	  n_events++;
	}

	if( badbcid[ichip][isca]>2 && first_retrig==false && bcid[ichip][isca]>20 && (bcid[ichip][isca]<890 ||bcid[ichip][isca]>915)) n_total_retriggers_trains[ichip]++;
	if( badbcid[ichip][isca]>2 && bcid[ichip][isca]>20 && (bcid[ichip][isca]<890 ||bcid[ichip][isca]>915)) n_total_retriggers[ichip]++;
	if( badbcid[ichip][isca]==0 && bcid[ichip][isca]>20 && (bcid[ichip][isca]<890 ||bcid[ichip][isca]>915)) {
	  n_total_triggers[ichip]++;
	  for(int ichn=0;ichn<64;ichn++) if(gain_hit_high[ichip][isca][ichn]==1) h_signal->Fill(charge_hiGain[ichip][isca][ichn]);
	}

	if(badbcid[ichip][isca]>2 && first_retrig==true && badbcid[ichip][isca-1]==0) {
	  first_retrig=false;
	  length=0;
	  n_trains++;
	  bcid_retrigger=bcid[ichip][isca];
	}

	//train
	if(badbcid[ichip][isca]>2 && first_retrig==true && bcid[ichip][isca]>20 && (bcid[ichip][isca]<890 ||bcid[ichip][isca]>915)) {
	  int ntriggers=0;
	  for(int ichn=0; ichn<64; ichn++) 
	    if(gain_hit_high[ichip][isca][ichn]==1) ntriggers++;
	  
	  h_retrigger_train[ichip]->Fill(isca-first_sca_retrig,ntriggers);
	  length++;
	}
	
	if(badbcid[ichip][isca]>2 && bcid[ichip][isca]>20 && (bcid[ichip][isca]<890 ||bcid[ichip][isca]>915) ) {
	  for(int ichn=0; ichn<64; ichn++) {
	    if(gain_hit_high[ichip][isca][ichn]==1) {
	      h_all_retriggering -> Fill(ichip,ichn);
	    }
	  }
	}
	
	if(badbcid[ichip][isca]>2 && first_retrig==false && bcid[ichip][isca]>20 && (bcid[ichip][isca]<890 ||bcid[ichip][isca]>915)) {
	  //bcid_retrigger=bcid[ichip][isca];
	  bad_triggers++;
	  first_retrig=true;
	  first_sca_retrig=isca;
	  //n_trains++;

	  for(int ichn=0; ichn<64; ichn++) {
	    if(gain_hit_high[ichip][isca][ichn]==1) {
	      h_first_retriggering -> Fill(ichip,ichn);
	    }
	  }
	}// first retrigg

	if(bcid_signal>0 && bcid_retrigger>0 && bcid[ichip][isca]>20 && (bcid[ichip][isca]<890 ||bcid[ichip][isca]>915)) {
	  h_bcid_correlation[ichip]->Fill(bcid_signal,bcid_retrigger-bcid_signal);
	  bcid_signal=0;
	  bcid_retrigger=0;
	}

      }

      
      h_retrigger_train_length[ichip]->Fill(length);
      h_retrigger_n_trains[ichip]->Fill(n_trains);
      h_signal_n_events[ichip]->Fill(n_events);

      good[ichip]->Fill(good_triggers);
      bad[ichip]->Fill(bad_triggers);
      good_vs_bad[ichip]->Fill(good_triggers,bad_triggers);

    }
  }


  for(int ichip=0; ichip<16; ichip++) {

    h_total_triggers->Fill(ichip, n_total_triggers[ichip]);
    if(n_total_triggers[ichip]>0) {
      h_total_retriggers->Fill(ichip, n_total_retriggers[ichip]/n_total_triggers[ichip]);
      h_total_retriggers_trains->Fill(ichip, n_total_retriggers_trains[ichip]/n_total_triggers[ichip]);
    } else {
      if(n_total_retriggers[ichip]>0) h_total_retriggers->Fill(ichip, 1);
      else h_total_retriggers->Fill(ichip, 0);
      if(n_total_retriggers_trains[ichip]>0) h_total_retriggers_trains->Fill(ichip, 1);
      else h_total_retriggers_trains->Fill(ichip, 0);
    }

    h_total_triggers_2d->Fill(ichip,0.,n_total_triggers[ichip]);
    h_total_triggers_2d->Fill(ichip,1.,n_total_retriggers_trains[ichip]);
    h_total_triggers_2d->Fill(ichip,2.,n_total_retriggers[ichip]);
  }
  
  TFile *summary = new TFile("results_retriggers/Retriggers"+slboard+".root" , "RECREATE");
  summary->cd();

  h_total_planeE->Write();
  h_total_negE->Write();
  h_total_negE_nohit->Write();
	
  h_first_retriggering->Write();
  h_all_retriggering->Write();
    
  h_total_triggers->Write();
  
  h_total_retriggers->Write();
  h_total_retriggers_trains->Write();
  h_total_triggers_2d->Write();
  h_signal->Write();
  
  TCanvas *canvas0 = new TCanvas("retriggers_map","retriggers_map",1200,600);
  canvas0->Divide(2,3);
  canvas0->cd(1);
  gPad->SetLogz();
  h_total_planeE->SetTitle("Plane events (>4hits)" +slboard);
  h_total_planeE->GetXaxis()->SetTitle("chip");
  h_total_planeE->GetXaxis()->SetTitle("sca");
  h_total_planeE->Draw("colz");

  canvas0->cd(2);
  gPad->SetLogz();
  h_total_negE->SetTitle("Negative hits (>4hits)" +slboard);
  h_total_negE->GetXaxis()->SetTitle("chip");
  h_total_negE->GetXaxis()->SetTitle("sca");
  h_total_negE->Draw("colz");
    
  canvas0->cd(3);
  gPad->SetLogz();
  h_total_negE_nohit->SetTitle("Negative no-hits (>4)" +slboard);
  h_total_negE_nohit->GetXaxis()->SetTitle("chip");
  h_total_negE_nohit->GetXaxis()->SetTitle("sca");
  h_total_negE_nohit->Draw("colz");
  
  canvas0->cd(4);
  h_first_retriggering->SetTitle("Starting retriggers " +slboard);
  // h_first_retriggering->GetZaxis()->SetRangeUser(1,15000);
  h_first_retriggering->GetXaxis()->SetTitle("chip");
  h_first_retriggering->GetXaxis()->SetTitle("channel");
  h_first_retriggering->Draw("colz");

  canvas0->cd(5);
  gPad->SetLogz();
  h_all_retriggering->SetTitle("All retrigger cells " +slboard);
  // h_all_retriggering->GetZaxis()->SetRangeUser(1,15000);
  h_all_retriggering->GetXaxis()->SetTitle("chip");
  h_all_retriggering->GetXaxis()->SetTitle("channel");
  h_all_retriggering->Draw("colz");

  canvas0->cd(6);
  h_total_triggers->SetTitle("Good+Ill Triggers " +slboard);
  //  h_total_triggers_2d->GetZaxis()->SetRangeUser(1,15000);
  h_total_triggers->GetXaxis()->SetTitle("chip");
  h_total_triggers->Draw("histo");
  
  canvas0->Write();
  
  //  h_total_trig
  TCanvas *canvas = new TCanvas("summary","summary",600,600);
  canvas->cd();
  //  h_total_triggers->Draw("histo");
  h_total_retriggers->GetYaxis()->SetRangeUser(0,1.1);
    h_total_retriggers->GetYaxis()->SetTitle("Nretrig/Ntrig");
  h_total_retriggers->SetLineColor(2);
  h_total_retriggers->SetLineWidth(2);
  h_total_retriggers->SetLineStyle(2);
  h_total_retriggers->Draw("histo");
  h_total_retriggers_trains->SetLineColor(2);
  h_total_retriggers_trains->SetLineWidth(2);
  h_total_retriggers_trains->Draw("histosame");
  //h_total_retriggers->GetYaxis()->SetRangeUser(0,1);
  //h_total_retriggers->SetLineColor(2);
  //h_total_retriggers->SetLineWidth(2);
  //h_total_retriggers->SetLineStyle(2);
  //h_total_retriggers->Draw("same");

  canvas->Write();
				
  for(int ichip=0; ichip<16; ichip++) {
    h_bcid_correlation[ichip]->GetXaxis()->SetTitle("bcid trigger");
    h_bcid_correlation[ichip]->GetYaxis()->SetTitle("bcid retrigger - bcid trigger");
    h_bcid_correlation[ichip]->Write();

    h_bcid2_correlation[ichip]->GetXaxis()->SetTitle("bcid trigger");
    h_bcid2_correlation[ichip]->GetYaxis()->SetTitle("bcid next trigger - bcid trigger");
    h_bcid2_correlation[ichip]->Write();
    
    h_retrigger_train[ichip]->Write();
    h_retrigger_train_length[ichip]->Write();
    h_retrigger_n_trains[ichip]->Write();
    h_signal_n_events[ichip]->Write();

    good[ichip]->Write();
    bad[ichip]->Write();
    good_vs_bad[ichip]->Write();
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

