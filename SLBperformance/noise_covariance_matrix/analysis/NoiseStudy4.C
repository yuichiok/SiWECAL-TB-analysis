#include "Fit.h"

void CreateHistos(TString run="3GeVMIPscan", int igain=0, int isca=0){

  TString map="../../../mapping/fev10_chip_channel_x_y_mapping.txt";
  //  if(layer==0 || layer==2)  map="../../mapping/fev11_cob_chip_channel_x_y_mapping.txt";
  ReadMap(map,0);
  
  
  //  TH2F* cov_matrix = OpenCovMatrix("3GeVMIPscan","highgain",5,0);
  //  cov_matrix->Draw("colz");

  TString gain[2]={"lowgain","highgain"};

  
  TString name=TString::Format("%s_%s_sca%i",run.Data(),gain[igain].Data(),isca-1);
  TString newfile=TString::Format("Summary_%s.root",name.Data());
  
  TFile *_file1 = new TFile(newfile,"RECREATE");
  
  for(int layer=0; layer<1; layer++) {
    
    TH2F * h2_i   = new TH2F(TString::Format("2d noise-incoherent layer%i",layer),TString::Format("incoherent layer%i ; chip; channel ; #sigma_{i}",layer),16,-0.5,15.5,64,-0.5,63.5);
    TH2F * h2_c1  = new TH2F(TString::Format("2d noise-coherent-1 layer%i",layer),TString::Format("noise-coherent-1 layer%i ; chip ; channel ; #sigma_{c1}",layer),16,-0.5,15.5,64,-0.5,63.5);
    TH2F * h2_c2 = new TH2F(TString::Format("2d noise-coherent-2 layer%i",layer),TString::Format("noise-coherent-2 layer%i ; chip; channel ; #sigma_{c2}",layer),16,-0.5,15.5,64,-0.5,63.5);
    
    TH2F * h2xy_i   = new TH2F(TString::Format("2d noise-incoherent layer%i - xy",layer),TString::Format("incoherent layer%i ; x; y ; #sigma_{i}",layer),32,-90,90,32,-90,90);
    TH2F * h2xy_c1  = new TH2F(TString::Format("2d noise-coherent-1 layer%i -xy",layer),TString::Format("noise-coherent-1 layer%i ; x ; y ; #sigma_{c1}",layer),32,-90,90,32,-90,90);
    TH2F * h2xy_c2 = new TH2F(TString::Format("2d noise-coherent-2 layer%i -xy",layer),TString::Format("noise-coherent-2 layer%i ; x; y ; #sigma_{c2}",layer),32,-90,90,32,-90,90);
    
    TH2F * h2_e_i   = new TH2F(TString::Format("2d error-noise-incoherent layer%i",layer),TString::Format("incoherent layer%i ; chip; channel ; E#sigma_{i}",layer),16,-0.5,15.5,64,-0.5,63.5);
    TH2F * h2_e_c1  = new TH2F(TString::Format("2d error-noise-coherent-1 layer%i",layer),TString::Format("noise-coherent-1 layer%i ; chip ; channel ; E#sigma_{c1}",layer),16,-0.5,15.5,64,-0.5,63.5);
    TH2F * h2_e_c2 = new TH2F(TString::Format("2d error-noise-coherent-2 layer%i",layer),TString::Format("noise-coherent-2 layer%i ; chip; channel ; E#sigma_{c2}",layer),16,-0.5,15.5,64,-0.5,63.5);
    
    TH2F * h2_exy_i   = new TH2F(TString::Format("2d error-noise-incoherent layer%i - xy",layer),TString::Format("incoherent layer%i ; x; y ; E#sigma_{i}",layer),32,-90,90,32,-90,90);
    TH2F * h2_exy_c1  = new TH2F(TString::Format("2d error-noise-coherent-1 layer%i -xy",layer),TString::Format("noise-coherent-1 layer%i ; x ; y ; E#sigma_{c1}",layer),32,-90,90,32,-90,90);
    TH2F * h2_exy_c2 = new TH2F(TString::Format("2d error-noise-coherent-2 layer%i -xy",layer),TString::Format("noise-coherent-2 layer%i ; x; y ; E#sigma_{c2}",layer),32,-90,90,32,-90,90);
    
    TH2F * h2_ped   = new TH2F(TString::Format("pedestal_layer%i",layer),TString::Format("pedestal_layer%i ; chip; chn ; ADC",layer),16,-0.5,15.5,64,-0.5,63.5);
    TH2F * h2_ped_xy   = new TH2F(TString::Format("pedestal_layer%i_xy",layer),TString::Format("pedestal_layer%i ; x; y ; ADC",layer),32,-90,90,32,-90,90);
    
    TH2F * h2_eped   = new TH2F(TString::Format("error_pedestal_layer%i",layer),TString::Format("error_pedestal_layer%i ; chip; chn ; ADC",layer),16,-0.5,15.5,64,-0.5,63.5);
    TH2F * h2_eped_xy   = new TH2F(TString::Format("error_pedestal_layer%i_xy",layer),TString::Format("error_pedestal_layer%i ; x; y ; ADC",layer),32,-90,90,32,-90,90);
    
    
    for(int chip=12; chip<13; chip++) {
      
      int fitvalue=Fit(name,gain[igain],layer,chip);
      if(fitvalue==0) continue;
      for(int j=0; j<64; j++) {
	
	h2_ped->Fill(chip,j,pedestal[j]);
	h2_ped_xy->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]),pedestal[j]);
	
	h2_eped->Fill(chip,j,epedestal[j]);
	h2_eped_xy->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]),epedestal[j]);
	
	h2_i->Fill(chip,j,sigma_i[j]);
	h2_c1->Fill(chip,j,sigma_c1[j]);
	h2_c2->Fill(chip,j,sigma_c2[j]);
	
	h2xy_i->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]),sigma_i[j]);
	h2xy_c1->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]),sigma_c1[j]);
	h2xy_c2->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]),sigma_c2[j]);
	
	h2_e_i->Fill(chip,j,esigma_i[j]);
	h2_e_c1->Fill(chip,j,esigma_c1[j]);
	h2_e_c2->Fill(chip,j,esigma_c2[j]);
	
	h2_exy_i->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]),esigma_i[j]);
	h2_exy_c1->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]),esigma_c1[j]);
	  h2_exy_c2->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]),esigma_c2[j]);
      }
      
    }//chip
    _file1->cd();
    h2_ped->Write();
    h2_ped_xy->Write();
    
    h2_eped->Write();
    h2_eped_xy->Write();
    //canvas->Write();
    h2_i->Write();
    h2_c1->Write();
    h2_c2->Write();
    h2xy_i->Write();
    h2xy_c1->Write();
    h2xy_c2->Write();
    
    //canvas2->Write();
    h2_e_i->Write();
    h2_e_c1->Write();
    h2_e_c2->Write();
    h2_exy_i->Write();
    h2_exy_c1->Write();
    h2_exy_c2->Write();
    
    //delete canvas;
    //delete canvas2;
  }//layer
  _file1->Close();
    // }//gain


}




void AverageHistos(TString name_="3GeVMIPscan",TString st_gain="highgain", int igain=0, int run0=43, int runN=123){

  TString map="../../../mapping/fev10_chip_channel_x_y_mapping.txt";
  //  if(layer==0 || layer==2)  map="../../mapping/fev11_cob_chip_channel_x_y_mapping.txt";
  ReadMap(map,0);
  
  int nchips=16;
  int nlayers=15;
  
  TString gain="";
 

  // for(int igain=10; igain<11; igain++) {

  gain=TString::Format("%s_sca%i",st_gain.Data(),igain-1);
  if(igain==0)  gain=st_gain;

  double noise_i[15][16][64]={0};
  double noise_c1[15][16][64]={0};
  double noise_c2[15][16][64]={0};
  double N_i[15][16][64]={0};
  double N_c1[15][16][64]={0};
  double N_c2[15][16][64]={0};

  double noise2_i[15][16][64]={0};
  double noise2_c1[15][16][64]={0};
  double noise2_c2[15][16][64]={0};
  double N2_i[15][16][64]={0};
  double N2_c1[15][16][64]={0};
  double N2_c2[15][16][64]={0};

  double pedestal[15][16][64]={0};
  double error_pedestal[15][16][64]={0};
  double error_pedestal_2[15][16][64]={0};
  double n_pedestal[15][16][64]={0};


  TFile *_file1 = new TFile(TString::Format("%s_mean_c3/Summary_%s_%s.root",st_gain.Data(),name_.Data(),gain.Data()),"RECREATE");
  _file1->cd();
  ReadMasked("../../../masked/masked_channels_3GeV_22degrees_run_050126.txt");
  for(int layer=0; layer<15; layer++) {
 
    TH2F * h2masked   = new TH2F(TString::Format("masked-MIPs layer%i",layer),TString::Format("masked layer%i ; x; y ; #sigma_{i}",layer),16,-0.5,15.5,64,-0.5,63.5);
    TH2F * h2maskedxy   = new TH2F(TString::Format("masked-MIPs layer%i - xy",layer),TString::Format("maskedlayer%i ; x; y ; #sigma_{i}",layer),32,-90,90,32,-90,90);
    for(int chip=0; chip<16; chip++) {
      for(int j=0; j<64; j++) {

	if(masked[14-layer][chip][j]==1) {
	  h2masked->Fill(chip,j,100);
	  h2maskedxy->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]),100);
	}
      }
    }
    h2masked->Write();
    h2maskedxy->Write();
  }
  ReadMasked("../../../masked/masked_channels_3GeV_W_22degrees_run_050186.txt");
  for(int layer=0; layer<15; layer++) {
 
    TH2F * h2masked   = new TH2F(TString::Format("masked-W layer%i",layer),TString::Format("masked layer%i ; x; y ; #sigma_{i}",layer),16,-0.5,15.5,64,-0.5,63.5);
    TH2F * h2maskedxy   = new TH2F(TString::Format("masked-W layer%i - xy",layer),TString::Format("maskedlayer%i ; x; y ; #sigma_{i}",layer),32,-90,90,32,-90,90);
    for(int chip=0; chip<16; chip++) {
      for(int j=0; j<64; j++) {

	if(masked[14-layer][chip][j]==1) {
	  h2masked->Fill(chip,j,100);
	  h2maskedxy->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]),100);
	}
      }
    }
    h2masked->Write();
    h2maskedxy->Write();
  }
	  
  for (int iname=run0; iname<runN+1; iname ++) {
    
    TString name=TString::Format("%s_run_0500%i_%s",name_.Data(),iname,gain.Data());//.root3GeV_22degrees_run_050043_highgain";
    if(iname>99) name=TString::Format("%s_run_050%i_%s",name_.Data(),iname,gain.Data());
    TFile * file = TFile::Open(TString::Format("%s_mean_c3/Summary_%s.root",st_gain.Data(),name.Data()));

    for(int layer=0; layer<nlayers; layer++) {
      
      TH2F * h2_i=(TH2F*)file->Get(TString::Format("2d noise-incoherent layer%i",layer));
      TH2F * h2_c1=(TH2F*)file->Get(TString::Format("2d noise-coherent-1 layer%i",layer));
      TH2F * h2_c2=(TH2F*)file->Get(TString::Format("2d noise-coherent-2 layer%i",layer));

      TH2F * h2e_i=(TH2F*)file->Get(TString::Format("2d error-noise-incoherent layer%i",layer));
      TH2F * h2e_c1=(TH2F*)file->Get(TString::Format("2d error-noise-coherent-1 layer%i",layer));
      TH2F * h2e_c2=(TH2F*)file->Get(TString::Format("2d error-noise-coherent-2 layer%i",layer));

      TH2F * h2ped=(TH2F*)file->Get(TString::Format("pedestal_layer%i",layer));
      TH2F * h2eped=(TH2F*)file->Get(TString::Format("error_pedestal_layer%i",layer));


      for(int chip=0; chip<nchips; chip++) {

	for(int j=0; j<64; j++) {
	  //   if(layer==0 && chip==15) cout<<h2_i->GetBinContent(chip+1,j+1)<<" "<<h2e_i->GetBinContent(chip+1,j+1)<<endl;


	  if(h2ped->GetBinContent(chip+1,j+1)>0) {
	    pedestal[layer][chip][j]+=h2ped->GetBinContent(chip+1,j+1)/h2eped->GetBinContent(chip+1,j+1);
	    error_pedestal[layer][chip][j]+=1.0/h2eped->GetBinContent(chip+1,j+1);
	    error_pedestal_2[layer][chip][j]+=h2eped->GetBinContent(chip+1,j+1);
	    n_pedestal[layer][chip][j]++;

	    // if(j==4)  cout<<" *************** "<<iname<<" "<<gain<<" "<<h2ped->GetBinContent(chip+1,j+1)<<" "<<h2eped->GetBinContent(chip+1,j+1)<<endl;
	  }
	  
	  if(h2_i->GetBinContent(chip+1,j+1)>0) {
  	    if(h2e_i->GetBinContent(chip+1,j+1)>0.01 || j==0) {
	      N_i[layer][chip][j]+=1./h2e_i->GetBinContent(chip+1,j+1);
	      noise_i[layer][chip][j]+=h2_i->GetBinContent(chip+1,j+1)/h2e_i->GetBinContent(chip+1,j+1);
	      N_i[layer][chip][j]+=1./h2e_i->GetBinContent(chip+1,j+1);

	    }
	    if(h2e_c1->GetBinContent(chip+1,j+1)>0.001) {
	      noise_c1[layer][chip][j]+=h2_c1->GetBinContent(chip+1,j+1)/h2e_c1->GetBinContent(chip+1,j+1);
	      N_c1[layer][chip][j]+=1./h2e_c1->GetBinContent(chip+1,j+1);
	    }
	    if(h2e_c2->GetBinContent(chip+1,j+1)>0.001) {
	      noise_c2[layer][chip][j]+=h2_c2->GetBinContent(chip+1,j+1)/h2e_c2->GetBinContent(chip+1,j+1);
	      N_c2[layer][chip][j]+=1./h2e_c2->GetBinContent(chip+1,j+1);
	    }
	  }
	    
	  if(h2_i->GetBinContent(chip+1,j+1)>0. && h2_i->GetBinContent(chip+1,j+1)<5.5 ) {
	    if(h2e_i->GetBinContent(chip+1,j+1)>0.01 || j==0) {
	      noise2_i[layer][chip][j]+=h2_i->GetBinContent(chip+1,j+1)/h2e_i->GetBinContent(chip+1,j+1);
	      N2_i[layer][chip][j]+=1./h2e_i->GetBinContent(chip+1,j+1);
		
	    }
	    if(h2e_c1->GetBinContent(chip+1,j+1)>0.001) {
	      noise2_c1[layer][chip][j]+=h2_c1->GetBinContent(chip+1,j+1)/h2e_c1->GetBinContent(chip+1,j+1);
	      N2_c1[layer][chip][j]+=1./h2e_c1->GetBinContent(chip+1,j+1);
	    }
	    if(h2e_c2->GetBinContent(chip+1,j+1)>0.001) {
	      noise2_c2[layer][chip][j]+=h2_c2->GetBinContent(chip+1,j+1)/h2e_c2->GetBinContent(chip+1,j+1);
	      N2_c2[layer][chip][j]+=1./h2e_c2->GetBinContent(chip+1,j+1);
	    }
	  }

	}//chan
	  
      }//chip
    }//layer
    file->Close();
  }//iname


  TH2F * hacc_i   = new TH2F("2d noise-incoherent Bad","Layers with sigma_i>4.8ADC ; chip; chn ; NLayers",16,-0.5,15.5,64,-0.5,63.5);
  TH2F * hacc_c1   = new TH2F("2d noise-coherent-1 Bad","Layers with sigma_c1>2.5ADC ; chip; chn ; NLayers",16,-0.5,15.5,64,-0.5,63.5);
  TH2F * hacc_c2   = new TH2F("2d noise-coherent-2 Bad","Layers with sigma_c2>2.ADC ; chip; chn ; NLayers",16,-0.5,15.5,64,-0.5,63.5);

  TH2F * haccxy_i   = new TH2F("2dxy noise-incoherent Bad","Layers with sigma_i>4.8ADC ; x; y ; NLayers",32,-90,90,32,-90,90);
  TH2F * haccxy_c1   = new TH2F("2dxy noise-coherent-1 Bad","Layers with sigma_c1>2.5ADC ; x; y ; NLayers",32,-90,90,32,-90,90);
  TH2F * haccxy_c2   = new TH2F("2dxy noise-coherent-2 Bad","Layers with sigma_c2>2.ADC ; x; y ; NLayers",32,-90,90,32,-90,90);

  for(int layer=0; layer<15; layer++) {
    TH2F * h2_i   = new TH2F(TString::Format("2d noise-incoherent layer%i",layer),TString::Format("incoherent layer%i ;chip ; chn  ; #sigma_{i}",layer),16,-0.5,15.5,64,-0.5,63.5);
    TH2F * h2_c1  = new TH2F(TString::Format("2d noise-coherent-1 layer%i",layer),TString::Format("noise-coherent-1 layer%i ;chip ; chn  ; #sigma_{c1}",layer),16,-0.5,15.5,64,-0.5,63.5);
    TH2F * h2_c2 = new TH2F(TString::Format("2d noise-coherent-2 layer%i",layer),TString::Format("noise-coherent-2 layer%i ;chip ; chn  ; #sigma_{c2}",layer),16,-0.5,15.5,64,-0.5,63.5);

    TH2F * h2xy_i   = new TH2F(TString::Format("2d noise-incoherent layer%i - xy",layer),TString::Format("incoherent layer%i ; x ; y ; #sigma_{i}",layer),32,-90,90,32,-90,90);
    TH2F * h2xy_c1  = new TH2F(TString::Format("2d noise-coherent-1 layer%i -xy",layer),TString::Format("noise-coherent-1 layer%i ; x ; y ; #sigma_{c1}",layer),32,-90,90,32,-90,90);
    TH2F * h2xy_c2 = new TH2F(TString::Format("2d noise-coherent-2 layer%i -xy",layer),TString::Format("noise-coherent-2 layer%i ; x ; y ; #sigma_{c2}",layer),32,-90,90,32,-90,90);

    TH2F * h2pedestal   = new TH2F(TString::Format("pedestal_layer%i",layer),TString::Format("pedestal_layer%i ;chip ; chn  ; ADC",layer),16,-0.5,15.5,64,-0.5,63.5);
    TH2F * h2pedestal_xy   = new TH2F(TString::Format("pedestal_xy_layer%i",layer),TString::Format("pedestal_layer%i ; x ; y ; ADC",layer),32,-90,90,32,-90,90);

    TH2F * h2errorpedestal   = new TH2F(TString::Format("rms_pedestal_layer%i",layer),TString::Format("rms_pedestal_layer%i ;chip ; chn  ; ADC",layer),16,-0.5,15.5,64,-0.5,63.5);
    TH2F * h2errorpedestal_xy   = new TH2F(TString::Format("rms_pedestal_xy_layer%i",layer),TString::Format("rms_pedestal_layer%i ; x ; y ; ADC",layer),32,-90,90,32,-90,90);  

    for(int chip=0; chip<16; chip++) {
      for(int j=0; j<64; j++) {
	//  if(layer==0 && chip==15) cout<<N2_i[layer][chip][j]<<" "<<N_i[layer][chip][j]<<endl;
	//-----------
	if(error_pedestal[layer][chip][j]>0) {
	  pedestal[layer][chip][j]/=error_pedestal[layer][chip][j];
	  h2pedestal->Fill(chip,j,pedestal[layer][chip][j]);
	  h2pedestal_xy->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]),pedestal[layer][chip][j]);
	  h2errorpedestal->Fill(chip,j,error_pedestal_2[layer][chip][j]/n_pedestal[layer][chip][j]);
	  h2errorpedestal_xy->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]),error_pedestal_2[layer][chip][j]/n_pedestal[layer][chip][j]);
	}
		  
	if(N2_i[layer][chip][j]>0) {
	  noise2_i[layer][chip][j]/=N2_i[layer][chip][j];
	  h2_i->Fill(chip,j,noise2_i[layer][chip][j]);
	  h2xy_i->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]),noise2_i[layer][chip][j]);
	  if(noise2_i[layer][chip][j]>4.8) {
	    hacc_i->Fill(chip,j);
	    haccxy_i->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]));
	  }	      
	} else {
	  if(N_i[layer][chip][j]>0) {
	    noise_i[layer][chip][j]/=N_i[layer][chip][j];
	    h2_i->Fill(chip,j,noise_i[layer][chip][j]);
	    h2xy_i->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]),noise_i[layer][chip][j]);
	    if(noise_i[layer][chip][j]>4.8) {
	      hacc_i->Fill(chip,j);
	      haccxy_i->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]));
	    }
	  } 
	}
	//-----------
	if(N2_c1[layer][chip][j]>0) {
	  noise2_c1[layer][chip][j]/=N2_c1[layer][chip][j];
	  h2_c1->Fill(chip,j,noise2_c1[layer][chip][j]);
	  h2xy_c1->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]),noise2_c1[layer][chip][j]);
	  if(noise2_c1[layer][chip][j]>2.5) {
	    hacc_c1->Fill(chip,j);
	    haccxy_c1->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]));
	  }	      
	} else {
	  if(N_c1[layer][chip][j]>0) {
	    noise_c1[layer][chip][j]/=N_c1[layer][chip][j];
	    h2_c1->Fill(chip,j,noise_c1[layer][chip][j]);     
	    h2xy_c1->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]),noise_c1[layer][chip][j]);
	    if(noise_c1[layer][chip][j]>2.5) {
	      hacc_c1->Fill(chip,j);
	      haccxy_c1->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]));
	    }
	  } 
	}
	//-----------
	if(N2_c2[layer][chip][j]>0) {
	  noise2_c2[layer][chip][j]/=N2_c2[layer][chip][j];
	  h2_c2->Fill(chip,j,noise2_c2[layer][chip][j]);
	  h2xy_c2->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]),noise2_c2[layer][chip][j]);
	  if(noise2_c2[layer][chip][j]>2.0) {
	    hacc_c2->Fill(chip,j);
	    haccxy_c2->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]));
	  }	      
	} else {
	  if(N_c2[layer][chip][j]>0) {
	    noise_c2[layer][chip][j]/=N_c2[layer][chip][j];
	    h2_c2->Fill(chip,j,noise_c2[layer][chip][j]);     
	    h2xy_c2->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]),noise_c2[layer][chip][j]);
	    if(noise_c2[layer][chip][j]>2.0) {
	      hacc_c2->Fill(chip,j);
	      haccxy_c2->Fill(double(map_pointX[0][chip][j]),double(map_pointY[0][chip][j]));
	    }
	  } 
	}
	 	  
      }
    }
    _file1->cd();
    h2_i->Write();
    h2_c1->Write();
    h2_c2->Write();
      
    h2xy_i->Write();
    h2xy_c1->Write();
    h2xy_c2->Write();

    h2pedestal->Write();
    h2pedestal_xy->Write();
    h2errorpedestal->Write();
    h2errorpedestal_xy->Write();

  }
  hacc_i->Write();
  hacc_c1->Write();
  hacc_c2->Write();
    
  haccxy_i->Write();
  haccxy_c1->Write();
  haccxy_c2->Write();
  _file1->Close();
    
  //}

  
}

void NoiseStudy4(){

  TString run[8]={"test_noise_15lsabs_HV100","test_noise_15lsabs_HV125","test_noise_15lsabs_HV150"};
  for(int irun=2; irun<3; irun++) {
    for(int isca=1; isca<2;isca++) CreateHistos(run[irun],1,isca);
  }

}

