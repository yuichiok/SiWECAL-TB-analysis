//# Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)

#include "conf_struct.h"

void test_read_thresholds(TString filename="Run_Settings.txt", bool debug=true) {

  read_configuration_file("05262020/Run_Settings_it14.txt",false);

  TH2F* threshold_chip_chn[15];
  TH2F* threshold_x_y[15];

  TH2F* threshold_chip_chn_2[15];
  TH2F* threshold_x_y_2[15];
  for(int islab=0; islab<2; islab++) {
    TString map_name="../mapping/fev10_chip_channel_x_y_mapping.txt";

    // the two cobs are equipped with slboards 2.08 and 2.12 (26th May 2020)
    if(detector.slab[0][islab].add==8 || detector.slab[0][islab].add==12)
      map_name="../mapping/fev11_cob_chip_channel_x_y_mapping.txt";
    
    cout<<map_name<<endl;
    ReadMap(map_name);
    
    threshold_chip_chn[islab]= new TH2F(TString::Format("threshold_chip_chn_%i",islab),TString::Format("threshold_chip_chn_%i",islab),16,-0.5,15.5,64,-0.5,63.5);
    threshold_x_y[islab]= new TH2F(TString::Format("threshold_x_y_%i",islab),TString::Format("threshold_x_y_%i",islab),32,-90,90,32,-90,90);

    threshold_chip_chn_2[islab]= new TH2F(TString::Format("threshold_chip_chn_global_%i",islab),TString::Format("threshold_chip_chn_global_%i",islab),16,-0.5,15.5,64,-0.5,63.5);
    threshold_x_y_2[islab]= new TH2F(TString::Format("threshold_x_y_global_%i",islab),TString::Format("threshold_x_y_global_%i",islab),32,-90,90,32,-90,90);
    
    for(int ichip=0; ichip<16; ichip++) {
            
      for(int ichn=0; ichn<64; ichn++) {
	float th= -detector.slab[0][islab].asu[0].skiroc[ichip].chn_threshold_adj[ichn];
	float th_glob= detector.slab[0][islab].asu[0].skiroc[ichip].threshold_dac;

	if(detector.slab[0][islab].asu[0].skiroc[ichip].mask[ichn]==0) {
	  threshold_chip_chn[islab]->Fill(ichip,ichn,th);
	  threshold_x_y[islab]->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn],th);

	  threshold_chip_chn_2[islab]->Fill(ichip,ichn,th_glob);
	  threshold_x_y_2[islab]->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn],th_glob); 
	}
      }
    }

    TCanvas *canvas = new TCanvas(TString::Format("canvas_%i",islab),TString::Format("canvas_%i",islab),1600,600);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kLightTemperature);//RainBow);
    canvas->Divide(2,1);
    canvas->cd(1);
    threshold_chip_chn_2[islab]->GetXaxis()->SetTitle("CHIP");
    threshold_chip_chn_2[islab]->GetYaxis()->SetTitle("CHANNEL");
    threshold_chip_chn_2[islab]->GetZaxis()->SetRangeUser(220,280);
    threshold_chip_chn_2[islab]->Draw("colz");
    threshold_chip_chn[islab]->Draw("text0same");
    
    canvas->cd(2);
    threshold_x_y_2[islab]->GetXaxis()->SetTitle("x");
    threshold_x_y_2[islab]->GetYaxis()->SetTitle("y");
    threshold_x_y_2[islab]->GetZaxis()->SetRangeUser(220,280);
    threshold_x_y_2[islab]->Draw("colz");
    threshold_x_y[islab]->Draw("text0same");
    

  }

  
}
