//# Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)

#include "conf_struct.h"


void test_read_masked_channels_inputlist(TString filename="Run_Settings.txt", bool debug=true) {
  read_configuration_file(filename,true);

  int masking[15][16][64]={0};

  std::ifstream reading_file("inputmaskedchannels.txt");
  if(!reading_file){
    cout<<" dameyo - damedame"<<endl;
  }
  
  TString tmpst;
  reading_file >> tmpst >> tmpst >> tmpst ;
  while(reading_file){
    Int_t tmp_layer=-1, tmp_chip = 1,tmp_channel = -1;
    reading_file >> tmp_layer >> tmp_chip >> tmp_channel;
    cout<<tmp_layer<<endl;
    if(tmp_layer>-1) masking[tmp_layer][tmp_chip][tmp_channel]=1;
  }


  TH2F* mask_chip_chn[15];
  TH2F* mask_x_y[15];

  float totalmasked[15];
  float totalmasked_chip[15][16];


  ofstream fout(filename+"_masked.txt",ios::out);
  fout<<"#masked_chns_list layer chip chns (0=not masked, 1=masked)"<<endl;

  int nslabs=15;
  for(int islab=0; islab<nslabs; islab++) {
    TString map_name="../mapping/fev10_chip_channel_x_y_mapping.txt";
    // TString map_name="../mapping/fev11_cob_chip_channel_x_y_mapping.txt";
    // TString map_name="../mapping/fev11_cob_rotate_chip_channel_x_y_mapping.txt";

    // the two cobs are equipped with slboards 2.08 and 2.12 (26th May 2020)
    if(detector.slab[0][islab].add==3 || detector.slab[0][islab].add==3)
      map_name="../mapping/fev11_cob_chip_channel_x_y_mapping.txt";
    
    cout<<" -------------------------------------------------------------------------------------------- " <<endl;
    cout<<"Slab idx" << islab<< "with slabAdd: "<<detector.slab[0][islab].add<< " and mapping file: "<< map_name<<endl;
    
    ReadMap(map_name);
    
    mask_chip_chn[islab]= new TH2F(TString::Format("mask_chip_chn_SLB%i",detector.slab[0][islab].add),TString::Format("mask_chip_chn_SLB%i",detector.slab[0][islab].add),16,-0.5,15.5,64,-0.5,63.5);
    mask_x_y[islab]= new TH2F(TString::Format("mask_x_y_SLB%i",detector.slab[0][islab].add),TString::Format("mask_x_y_SLB%i",detector.slab[0][islab].add),32,-90,90,32,-90,90);

    totalmasked[islab]=0;
    for(int ichip=0; ichip<16; ichip++) {
      totalmasked_chip[islab][ichip]=0;
      fout<<islab<<" "<<ichip;
      for(int ichn=0; ichn<64; ichn++) {

	if(masking[islab][ichip][ichn]==1) {
	  detector.slab[0][islab].asu[0].skiroc[ichip].mask[ichn]=1;
	  detector.slab[0][islab].asu[0].skiroc[ichip].preamplifier_mask[ichn]=1;
	  cout<<" masking "<<islab<<" "<<ichip<<" "<<ichn<<endl;
	}
	    
	if(detector.slab[0][islab].asu[0].skiroc[ichip].mask[ichn]==1) {
	  totalmasked[islab]++;
	  totalmasked_chip[islab][ichip]++;
	  detector.slab[0][islab].asu[0].skiroc[ichip].preamplifier_mask[ichn]=1;
	  cout<<" Fix PreaAmp "<<islab<<" "<<ichip<<" "<<ichn<<endl;

	}
	
	float msk= detector.slab[0][islab].asu[0].skiroc[ichip].mask[ichn]+0.001;
	mask_chip_chn[islab]->Fill(ichip,ichn,msk);
	mask_x_y[islab]->Fill(map_pointX[ichip][ichn],map_pointY[ichip][ichn],msk);
	fout<<" "<<detector.slab[0][islab].asu[0].skiroc[ichip].mask[ichn];
      }
      fout<<endl;
      cout<< "Skiroc idx:"<<ichip<<"  Masked cells:"<<100.*totalmasked_chip[islab][ichip]/64.<<"%"<<endl;
    }
    cout<< "TOTAL Masked cells:"<<100.*totalmasked[islab]/1024.<<"%"<<endl;
    

    TCanvas *canvas = new TCanvas(TString::Format("canvas_%i",islab),TString::Format("canvas_%i",islab),800,1600);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kInvertedDarkBodyRadiator);//Cherry);
    //TColor::InvertPalette();
    canvas->Divide(1,2);
    canvas->cd(1);
    // mask_chip_chn[islab]->SetTitle(TString::Format("SLABd%i",mapping_slab[islab]));
    mask_chip_chn[islab]->SetTitle(TString::Format("SlabAdd%i",islab));
    mask_chip_chn[islab]->GetXaxis()->SetTitle("CHIP");
    mask_chip_chn[islab]->GetYaxis()->SetTitle("CHANNEL");
    mask_chip_chn[islab]->Draw("col");
    canvas->cd(2);
    // mask_x_y[islab]->SetTitle(TString::Format("SLABd%i",mapping_slab[islab]));
    mask_x_y[islab]->SetTitle(TString::Format("SlabAdd%i",islab));
    mask_x_y[islab]->GetXaxis()->SetTitle("x");
    mask_x_y[islab]->GetYaxis()->SetTitle("y");
    mask_x_y[islab]->Draw("col");
    // canvas->Print(TString::Format("masked_channels_SLAB%i.png",mapping_slab[islab]));
    canvas->Print(TString::Format("plots/masked_channels_SlabAdd%i.png",islab));

  }

  cout << "################################################################" << endl;
  cout << "################################################################" << endl;
  cout << "####           CHECK THE MASK MAP FILE !!!!!!               ####" << endl;
  cout << "################################################################" << endl;
  cout << "################################################################" << endl;
  write_configuration_file("Run_Settings_new.txt");

  
}
