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
  
  for(int layer=0; layer<15; layer++) {
    
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
    
    
    for(int chip=0; chip<16; chip++) {
      
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

void NoiseStudy(){

  TString run[1]={"03102022_pedestal_13slabs"};
  for(int irun=0; irun<1; irun++) {
    for(int isca=1; isca<16;isca++) {
      CreateHistos(run[irun],1,isca);
      CreateHistos(run[irun],0,isca);
    }
  }

}

