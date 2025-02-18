//# Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)

#include "conf_struct.h"

//void test_read_masked_channels_summary(TString filename="15102021/Run_Settings_DataTaking_TB2021_15102021_WRITENbyDAQ", bool debug=true) {
void test_read_masked_channels_summary(TString filename="/mnt/win2/Run_Data/Run_ILC_03112022_cosmic_it16_Ascii/Run_Settings", bool debug=true) {
  //void test_read_masked_channels_summary(TString filename="run_150014/Run_Settings", bool debug=true) {  

  read_configuration_file(filename+".txt",true);

  int mapping_slab[15];
  mapping_slab[0]=18;
  mapping_slab[1]=23;
  mapping_slab[2]=17;
  mapping_slab[3]=22;//
  mapping_slab[4]=25;
  mapping_slab[5]=24;
  mapping_slab[6]=31;
  mapping_slab[7]=30;//
  mapping_slab[8]=21;
  mapping_slab[9]=20;//
  mapping_slab[10]=19;//
  mapping_slab[11]=15;
  mapping_slab[12]=14;//
  mapping_slab[13]=13;
  mapping_slab[14]=16;



  TH2F* mask_chip_chn = new TH2F("mask_chip_chn","mask_chip_chn",16,-0.5,15.5,64,-0.5,63.5);
  TH2F* mask_x_y = new TH2F("mask_x_y","mask_x_y",32,-90,90,32,-90,90);

  TH2F* mask_layer_chip = new TH2F("mask_chip_chn","mask_chip_chn",15,-0.5,14.5,16,-0.5,15.5);

  float totalmasked[15];
  float totalmasked_chip[15][16];

  ofstream fout(filename+"_masked.txt",ios::out);
  fout<<"#masked_chns_list layer chip chns (0=not masked, 1=masked)"<<endl;
  int nslabs=15;
  for(int islab=0; islab<nslabs; islab++) {
    TString map_name="../mapping/fev10_chip_channel_x_y_mapping.txt";

    // the two cobs are equipped with slboards 2.08 and 2.12 (26th May 2020)
    //    if(detector.slab[0][islab].add==8 || detector.slab[0][islab].add==12)
    //  map_name="../mapping/fev11_cob_chip_channel_x_y_mapping.txt";
    
    cout<<" -------------------------------------------------------------------------------------------- " <<endl;
    cout<<"Slab idx" << islab<< "with slabAdd: "<<detector.slab[0][islab].add<< " and mapping file: "<< map_name<<endl;
    
    ReadMap(map_name);
    
    totalmasked[islab]=0;
    for(int ichip=0; ichip<16; ichip++) {
      totalmasked_chip[islab][ichip]=0;
      fout<<islab<<" "<<ichip;
      for(int ichn=0; ichn<64; ichn++) {
	if(detector.slab[0][islab].asu[0].skiroc[ichip].mask[ichn]==1) {
	  totalmasked[islab]++;
	  totalmasked_chip[islab][ichip]++;
	}
	
	float msk= detector.slab[0][islab].asu[0].skiroc[ichip].mask[ichn];
	if(msk>0) {
	  mask_chip_chn->Fill(ichip,ichn);
	  mask_x_y->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn]);
	}
	fout<<" "<<detector.slab[0][islab].asu[0].skiroc[ichip].mask[ichn];
      }
      fout<<endl;
      cout<< "Skiroc idx:"<<ichip<<"  Masked cells:"<<100.*totalmasked_chip[islab][ichip]/64.<<"%"<<endl;
      mask_layer_chip->Fill(islab,ichip,totalmasked_chip[islab][ichip]+0.01);

    }
    cout<< "TOTAL Masked cells:"<<100.*totalmasked[islab]/1024.<<"%"<<endl;
    
  }
  

  TCanvas *canvas = new TCanvas("canvas_accum","canvas_accum",2700,500);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);//Cherry);
  gStyle->SetPadRightMargin(0.16);

  //TColor::InvertPalette();
  canvas->Divide(3,1);
  canvas->cd(1);
  mask_chip_chn->GetXaxis()->SetTitle("CHIP");
  mask_chip_chn->GetYaxis()->SetTitle("CHANNEL");
  mask_chip_chn->GetZaxis()->SetTitle("N-slabs");
  mask_chip_chn->Draw("colz");

  canvas->cd(2);
  mask_x_y->GetXaxis()->SetTitle("x");
  mask_x_y->GetYaxis()->SetTitle("y");
  mask_x_y->GetZaxis()->SetTitle("N-slabs");
  mask_x_y->Draw("colz");
  // canvas->Print(TString::Format("masked_channels_slboard%i.eps",;
  canvas->cd(3);
  mask_layer_chip->GetXaxis()->SetTitle("SLAB-ID");
  mask_layer_chip->GetYaxis()->SetTitle("CHIP");
  mask_layer_chip->GetZaxis()->SetTitle("masked channels");
  mask_layer_chip->GetZaxis()->SetRangeUser(0,35);
  mask_layer_chip->Draw("colz");

  
}
