//# Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)

#include "conf_struct.h"

void test_read_masked_channels(TString filename="Run_Settings.txt", bool debug=true) {

  read_configuration_file(filename,false);

  TH2F* mask_chip_chn[15];
  TH2F* mask_x_y[15];

  float totalmasked[15];
  float totalmasked_chip[15][16];
  
  int nslabs=15;
  for(int islab=0; islab<nslabs; islab++) {
    TString map_name="../mapping/fev10_chip_channel_x_y_mapping.txt";

    // the two cobs are equipped with slboards 2.08 and 2.12 (26th May 2020)
    if(detector.slab[0][islab].add==8 || detector.slab[0][islab].add==12)
      map_name="../mapping/fev11_cob_chip_channel_x_y_mapping.txt";
    
    cout<<" -------------------------------------------------------------------------------------------- " <<endl;
    cout<<"Slab idx" << islab<< "with slabAdd: "<<detector.slab[0][islab].add<< " and mapping file: "<< map_name<<endl;
    
    ReadMap(map_name);
    
    mask_chip_chn[islab]= new TH2F(TString::Format("mask_chip_chn_%i",islab),TString::Format("mask_chip_chn_%i",islab),16,-0.5,15.5,64,-0.5,63.5);
    mask_x_y[islab]= new TH2F(TString::Format("mask_x_y_%i",islab),TString::Format("mask_x_y_%i",islab),32,-90,90,32,-90,90);

    totalmasked[islab]=0;
    for(int ichip=0; ichip<16; ichip++) {
      totalmasked_chip[islab][ichip]=0;
      for(int ichn=0; ichn<64; ichn++) {
	if(detector.slab[0][islab].asu[0].skiroc[ichip].mask[ichn]==1) {
	  totalmasked[islab]++;
	  totalmasked_chip[islab][ichip]++;
	}
	
	float msk= detector.slab[0][islab].asu[0].skiroc[ichip].mask[ichn]+0.001;
	mask_chip_chn[islab]->Fill(ichip,ichn,msk);
	mask_x_y[islab]->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn],msk);
      }
      cout<< "Skiroc idx:"<<ichip<<"  Masked cells:"<<100.*totalmasked_chip[islab][ichip]/64.<<"%"<<endl;
    }
    cout<< "TOTAL Masked cells:"<<100.*totalmasked[islab]/1024.<<"%"<<endl;
    

    TCanvas *canvas = new TCanvas(TString::Format("canvas_%i",islab),TString::Format("canvas_%i",islab),1600,800);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kInvertedDarkBodyRadiator);//Cherry);
    //TColor::InvertPalette();
    canvas->Divide(2,1);
    canvas->cd(1);
    mask_chip_chn[islab]->GetXaxis()->SetTitle("CHIP");
    mask_chip_chn[islab]->GetYaxis()->SetTitle("CHANNEL");
    mask_chip_chn[islab]->Draw("col");
    canvas->cd(2);
    mask_x_y[islab]->GetXaxis()->SetTitle("x");
    mask_x_y[islab]->GetYaxis()->SetTitle("y");
    mask_x_y[islab]->Draw("col");
    

  }

  
}
