#include "Fit.h"

void SummaryPedestal(TString name_="3GeVMIPscan",TString st_gain="highgain"){

  TString map="../../mapping/fev10_chip_channel_x_y_mapping.txt";
  //  if(layer==0 || layer==2)  map="../../mapping/fev11_cob_chip_channel_x_y_mapping.txt";
  ReadMap(map,0);
  int nlayers=15;
  int nchips=16;
  TH1F* width_i_layer[15];
  for(int layer=0; layer<nlayers; layer++) width_i_layer[layer] =new TH1F(TString::Format("width_i_layer%i",layer),"widths of pedestal ; ADC  ",300,0,8);
  TH1F* width_layer[15];
  for(int layer=0; layer<nlayers; layer++) width_layer[layer] =new TH1F(TString::Format("width_all_layer%i",layer),"widths of pedestal ; ADC  ",300,0,8);
  
  
  TH2F* ped=new TH2F("ped","pedestal pos. ; Layer*20+Chip  ; SCA*100 +chn; ADC",300,-0.5,299.5,1500,-0.5,1499.5);
  TH2F* ped_rms=new TH2F("ped_rms","pedestal RMS. ; Layer*20+Chip  ; SCA*100 +chn; ADC",300,-0.5,299.5,1500,-0.5,1499.5);
  TH2F* ped_i=new TH2F("ped_i","pedestal incoherent noise. ; Layer*20+Chip  ; SCA*100 +chn; ADC",300,-0.5,299.5,1500,-0.5,1499.5);
  TH2F* ped_c1=new TH2F("ped_c1","pedestal coherent-1 noise ; Layer*20+Chip  ; SCA*100 +chn; ADC",300,-0.5,299.5,1500,-0.5,1499.5);
  TH2F* ped_c2=new TH2F("ped_c2","pedestal coherent-1 noise ; Layer*20+Chip  ; SCA*100 +chn; ADC",300,-0.5,299.5,1500,-0.5,1499.5);

  TH2F* ped_nofit=new TH2F("ped_nofit","pedestal pos. ; Layer*20+Chip  ; SCA*100 +chn; ADC",300,-0.5,299.5,1500,-0.5,1499.5);
  TH2F* ped_rms_nofit=new TH2F("ped_rms_nofit","pedestal RMS. ; Layer*20+Chip  ; SCA*100 +chn; ADC",300,-0.5,299.5,1500,-0.5,1499.5);
  TH2F* ped_i_nofit=new TH2F("ped_i_nofit","pedestal incoherent noise. ; Layer*20+Chip  ; SCA*100 +chn; ADC",300,-0.5,299.5,1500,-0.5,1499.5);
  TH2F* ped_c1_nofit=new TH2F("ped_c1_nofit","pedestal coherent-1 noise ; Layer*20+Chip  ; SCA*100 +chn; ADC",300,-0.5,299.5,1500,-0.5,1499.5);
  TH2F* ped_c2_nofit=new TH2F("ped_c2_nofit","pedestal coherent-1 noise ; Layer*20+Chip  ; SCA*100 +chn; ADC",300,-0.5,299.5,1500,-0.5,1499.5);
  
  ofstream fout_ped(TString::Format("../../../pedestals/Pedestal_method2_%s_%s.txt",name_.Data(),st_gain.Data()).Data(),ios::out);
  fout_ped<<"#pedestal results (method of covariance matrix) : PROTO15"<<endl;
  fout_ped<<"#layer chip channel";
  for(int isca=0; isca<15; isca++) fout_ped<<TString::Format(" ped_mean%i ped_rms%i noise_incoherent_ped%i noise_coherent1_ped%i noise_coherent2_ped%i",isca,isca,isca,isca,isca);
  fout_ped<<endl;

  for(int layer=0; layer<nlayers; layer++) {
    for(int chip=0; chip<nchips; chip++) {
      for(int j=0; j<64; j++) {
	double pedv=0.;
	double pedrmsv=0.;
	double pediv=0.;
	double pedc1v=0.;
	double pedc2v=0.;
	double pedv_[15]={0.};
	double pedrmsv_[15]={0.};
	double pediv_[15]={0.};
	double pedc1v_[15]={0.};
	double pedc2v_[15]={0.};
	
	double nped=0;
	double nnoise=0;

	int fitted[15]={0};
	TString sca="";
	fout_ped << layer<<" "<<chip <<" " <<j<< " ";

	for(int isca=0; isca<15; isca++) {
	  sca=TString::Format("%s_sca%i",st_gain.Data(),isca);
	  TFile *file = TFile::Open(TString::Format("Summary_%s_%s.root",name_.Data(),sca.Data()));
	  TH2F * h2_ped=(TH2F*)file->Get(TString::Format("pedestal_layer%i",layer));
	  TH2F * h2_eped=(TH2F*)file->Get(TString::Format("error_pedestal_layer%i",layer));

	  TH2F * h2_i=(TH2F*)file->Get(TString::Format("2d noise-incoherent layer%i",layer));
	  TH2F * h2_c1=(TH2F*)file->Get(TString::Format("2d noise-coherent-1 layer%i",layer));
	  TH2F * h2_c2=(TH2F*)file->Get(TString::Format("2d noise-coherent-2 layer%i",layer));
	  if(h2_ped->GetBinContent(chip+1,j+1)>10) {
	    pedv+=h2_ped->GetBinContent(chip+1,j+1);
	    pedrmsv+=h2_eped->GetBinContent(chip+1,j+1);
	    pediv+=h2_i->GetBinContent(chip+1,j+1);
	    pedc1v+=h2_c1->GetBinContent(chip+1,j+1);
	    pedc2v+=h2_c2->GetBinContent(chip+1,j+1);

	    pedv_[isca]=h2_ped->GetBinContent(chip+1,j+1);
	    pedrmsv_[isca]=h2_eped->GetBinContent(chip+1,j+1);
	    pediv_[isca]=h2_i->GetBinContent(chip+1,j+1);
	    pedc1v_[isca]=h2_c1->GetBinContent(chip+1,j+1);
	    pedc2v_[isca]=h2_c2->GetBinContent(chip+1,j+1);
	    nped++;
	    fitted[isca]=1;
	    if(h2_i->GetBinContent(chip+1,j+1)>0) nnoise++;
	    /* ped->Fill(20*layer+chip,100*isca+j,h2_ped->GetBinContent(chip+1,j+1));
	    ped_i->Fill(20*layer+chip,100*isca+j,h2_i->GetBinContent(chip+1,j+1));
	    ped_c1->Fill(20*layer+chip,100*isca+j,h2_c1->GetBinContent(chip+1,j+1));
	    ped_c2->Fill(20*layer+chip,100*isca+j,h2_c2->GetBinContent(chip+1,j+1));*/
	    if(h2_i->GetBinContent(chip+1,j+1)>0) width_i_layer[layer]->Fill(h2_i->GetBinContent(chip+1,j+1));
	    if(h2_i->GetBinContent(chip+1,j+1)>0 && h2_c1->GetBinContent(chip+1,j+1)>0 && h2_c2->GetBinContent(chip+1,j+1)>0)
	      width_layer[layer]->Fill(sqrt(pow(h2_i->GetBinContent(chip+1,j+1),2)+pow(h2_c1->GetBinContent(chip+1,j+1),2)+pow(h2_c2->GetBinContent(chip+1,j+1),2)));

	  }
	  delete file;//file->Close();
	}//isca
	for(int isca=0; isca<15; isca++) {
	  if(fitted[isca]==0 && nped>0 && nnoise>0) {
	    ped->Fill(20*layer+chip,100*isca+j,pedv/nped);
	    ped_rms->Fill(20*layer+chip,100*isca+j,pedrmsv/nped);
	    ped_i->Fill(20*layer+chip,100*isca+j,pediv/nnoise);
	    ped_c1->Fill(20*layer+chip,100*isca+j,pedc1v/nnoise);
	    ped_c2->Fill(20*layer+chip,100*isca+j,pedc2v/nnoise);
	    ped_nofit->Fill(20*layer+chip,100*isca+j,pedv/nped);
	    ped_rms_nofit->Fill(20*layer+chip,100*isca+j,pedrmsv/nped);
	    ped_i_nofit->Fill(20*layer+chip,100*isca+j,pediv/nnoise);
	    ped_c1_nofit->Fill(20*layer+chip,100*isca+j,pedc1v/nnoise);
	    ped_c2_nofit->Fill(20*layer+chip,100*isca+j,pedc2v/nnoise);
	    fout_ped<<pedv/nped<< " " <<-5<< " " <<pediv/nnoise<<" "<<pedc1v/nnoise<< " " <<pedc2v/nnoise<<" "; //bad fits for all but with data for average
	  } else {
	    if(fitted[isca]==1 && pediv_[isca]>0 ) {
	      fout_ped<<pedv_[isca]<< " " <<pedrmsv_[isca]<< " " <<pediv_[isca]<< " " <<pedc1v_[isca]<< " " <<pedc2v_[isca]<<" "; //good fits for pedestal and noise
	      ped->Fill(20*layer+chip,100*isca+j,pedv_[isca]);
	      ped_rms->Fill(20*layer+chip,100*isca+j,pedrmsv_[isca]);
	      ped_i->Fill(20*layer+chip,100*isca+j,pediv_[isca]);
	      ped_c1->Fill(20*layer+chip,100*isca+j,pedc1v_[isca]);
	      ped_c2->Fill(20*layer+chip,100*isca+j,pedc2v_[isca]);
	    }  else if(fitted[isca]==1 && pediv_[isca]==0 ) {
	      if(nnoise>0) {
		fout_ped<<pedv_[isca]<< " " <<pedrmsv_[isca]<< " "<<pediv/nnoise<< " " <<pedc1v/nnoise<< " " <<pedc2v/nnoise<<" ";//bad fit of noise, but with data for average
		ped_i->Fill(20*layer+chip,100*isca+j,pediv/nnoise);
		ped_rms->Fill(20*layer+chip,100*isca+j,pedrmsv_[isca]);
		ped_c1->Fill(20*layer+chip,100*isca+j,pedc1v/nnoise);
		ped_c2->Fill(20*layer+chip,100*isca+j,pedc2v/nnoise);
		ped_i_nofit->Fill(20*layer+chip,100*isca+j,pediv/nnoise);
		ped_c1_nofit->Fill(20*layer+chip,100*isca+j,pedc1v/nnoise);
		ped_c2_nofit->Fill(20*layer+chip,100*isca+j,pedc2v/nnoise);
	      }
	      else  fout_ped<<pedv_[isca]<< " " <<-5<< " " <<-5<< " " <<-5<< " " <<-5<<" "; //bad fits
	    } else fout_ped<<pedv_[isca]<< " " <<-10<<" " <<-10<< " " <<-10<< " " <<-10<<" "; //no stat for fits
	  }
	}//isca
      	fout_ped<<endl;
      }//chn
    }//chip
  }//layer

  TFile *_file1 = new TFile(TString::Format("../../../pedestals/PedSummary_method2_%s_%s.root",name_.Data(),st_gain.Data()),"RECREATE");
  _file1->cd();
  for(int layer=0; layer<nlayers; layer++) {
    width_i_layer[layer]->Write();
    width_layer[layer]->Write();
  }
					     
  TCanvas* canvas0= new TCanvas("PedAna_Summary","Pedestal Summary",1800,400);   
  canvas0->Divide(2,1);
  canvas0->cd(1);
  ped->GetZaxis()->SetRangeUser(180,300);
  ped->Draw("colz");
  canvas0->cd(2);
  //ped_i->GetZaxis()->SetRangeUser(2.0,4.5);
  ped_i->Draw("colz");
  // canvas0->Print(TString::Format("../../pedestals/%s_%sgain_PedSummary.png",run.Data(),gain.Data()));
  ped->Write();
  ped_i->Write();
  canvas0->Write();
  TCanvas* canvas2= new TCanvas("PedAna_Summary2","Pedestal Summary2",1800,400);   
  canvas2->Divide(2,1);
  canvas2->cd(1);
  // ped_c1->GetZaxis()->SetRangeUser(0,1.5);
  ped_c1->Draw("colz");
  canvas2->cd(2);
  //ped_c2->GetZaxis()->SetRangeUser(0.0,1.0);
  ped_c2->Draw("colz");
  // canvas0->Print(TString::Format("../../pedestals/%s_%sgain_PedSummary.png",run.Data(),gain.Data()));
  ped_c1->Write();
  ped_c2->Write();
  canvas2->Write();

  TCanvas* canvas02= new TCanvas("PedAna_Summary_nofit","Pedestal Summary",1800,400);   
  canvas02->Divide(2,1);
  canvas02->cd(1);
  ped_nofit->GetZaxis()->SetRangeUser(180,300);
  ped_nofit->Draw("colz");
  canvas02->cd(2);
  //ped_i_nofit->GetZaxis()->SetRangeUser(2.0,4.5);
  ped_i_nofit->Draw("colz");
  // canvas0->Print(TString::Format("../../pedestals/%s_%sgain_PedSummary.png",run.Data(),gain.Data()));
  ped_nofit->Write();
  ped_i_nofit->Write();
  canvas02->Write();
  TCanvas* canvas22= new TCanvas("PedAna_Summary2_nofit","Pedestal Summary2",1800,400);   
  canvas22->Divide(2,1);
  canvas22->cd(1);
  // ped_c1_nofit->GetZaxis()->SetRangeUser(0,1.5);
  ped_c1_nofit->Draw("colz");
  canvas22->cd(2);
  // ped_c2_nofit->GetZaxis()->SetRangeUser(0.0,1.0);
  ped_c2_nofit->Draw("colz");
  // canvas0->Print(TString::Format("../../pedestals/%s_%sgain_PedSummary.png",run.Data(),gain.Data()));
  ped_c1_nofit->Write();
  ped_c2_nofit->Write();
  canvas22->Write();
  
}

void SummaryPlots(){
  TString runname="03102022_pedestal_13slabs";
  SummaryPedestal(runname,"highgain");
  SummaryPedestal(runname,"lowgain");
}

