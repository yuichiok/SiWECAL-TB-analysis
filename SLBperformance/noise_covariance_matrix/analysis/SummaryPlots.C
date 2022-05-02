#include "Fit.h"

//int new_layer_ord[15]={5,6,0,1,4,2,3,7,8,9,10,11,12,13,14};
int new_layer_ord[15]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};//2,3,5,6,4,0,1,7,8,9,10,11,12,13,14};

void Maps(TString name_="3GeVMIPscan",TString st_gain="highgain"){

  int nlayers=15;
  int nchips=16;
  int nsca=1;
  gStyle->SetOptStat(0);

  TFile *_file1 = new TFile(TString::Format("../../../pedestals/PedMaps_method2_%s_%s.root",name_.Data(),st_gain.Data()),"RECREATE");
  for(int layer=0; layer<nlayers; layer++) {

    for(int isca=0; isca<nsca; isca++) {
      TString sca=TString::Format("%s_sca%i",st_gain.Data(),isca);
      TFile *file = TFile::Open(TString::Format("Summary_%s_%s.root",name_.Data(),sca.Data()));
      
      TH2F * h2_ped_xy=(TH2F*)file->Get(TString::Format("pedestal_layer%i_xy",new_layer_ord[layer]));
      TH2F * h2_i_xy=(TH2F*)file->Get(TString::Format("2d noise-incoherent layer%i - xy",new_layer_ord[layer]));
      TH2F * h2_c1_xy=(TH2F*)file->Get(TString::Format("2d noise-coherent-1 layer%i -xy",new_layer_ord[layer]));
      TH2F * h2_c2_xy=(TH2F*)file->Get(TString::Format("2d noise-coherent-2 layer%i -xy",new_layer_ord[layer]));

      
      TH2F * h2_ped=(TH2F*)file->Get(TString::Format("pedestal_layer%i",new_layer_ord[layer]));
      TH2F * h2_i=(TH2F*)file->Get(TString::Format("2d noise-incoherent layer%i",new_layer_ord[layer]));
      TH2F * h2_c1=(TH2F*)file->Get(TString::Format("2d noise-coherent-1 layer%i",new_layer_ord[layer]));
      TH2F * h2_c2=(TH2F*)file->Get(TString::Format("2d noise-coherent-2 layer%i",new_layer_ord[layer]));

   
      _file1->cd();

      TCanvas* canvas0= new TCanvas(TString::Format("Noise, x-y maps, layer%i",layer),TString::Format("Noise, x-y maps, layer%i",layer),1200,1200);   
      canvas0->Divide(2,2);
      canvas0->cd(1);
      h2_ped_xy->SetTitle("Pedestal, "+name_);
      h2_ped_xy->GetXaxis()->SetTitle("x");
      h2_ped_xy->GetYaxis()->SetTitle("y");
      h2_ped_xy->GetZaxis()->SetRangeUser(180,280);
      h2_ped_xy->Draw("colz");
      canvas0->cd(2);
      h2_i_xy->SetTitle("Incoherent noise, "+name_);
      h2_i_xy->GetXaxis()->SetTitle("x");
      h2_i_xy->GetYaxis()->SetTitle("y");
      h2_i_xy->GetZaxis()->SetRangeUser(0,5);
      h2_i_xy->Draw("colz");
      canvas0->cd(3);
      h2_c1_xy->SetTitle("C1 noise, "+name_);
      h2_c1_xy->GetXaxis()->SetTitle("x");
      h2_c1_xy->GetYaxis()->SetTitle("y");
      h2_c1_xy->GetZaxis()->SetRangeUser(0,5);
      h2_c1_xy->Draw("colz");
    
      canvas0->cd(4);
      h2_c2_xy->SetTitle("C2 noise, "+name_);
      h2_c2_xy->GetXaxis()->SetTitle("x");
      h2_c2_xy->GetYaxis()->SetTitle("y");
      h2_c2_xy->GetZaxis()->SetRangeUser(0,5);
      h2_c2_xy->Draw("colz");
      canvas0->Write();
      delete canvas0;

      TCanvas* canvas1= new TCanvas(TString::Format("Noise, chip-chn maps, layer%i",layer),TString::Format("Noise, chip-chn maps, layer%i",layer),1200,1200);   
      canvas1->Divide(2,2);
      canvas1->cd(1);
      h2_ped->SetTitle("Pedestal, "+name_);
      h2_ped->GetXaxis()->SetTitle("chip");
      h2_ped->GetYaxis()->SetTitle("chn");
      h2_ped->GetZaxis()->SetRangeUser(180,280);
      h2_ped->Draw("colz");
      canvas1->cd(2);
      h2_i->SetTitle("Incoherent noise, "+name_);
      h2_i->GetXaxis()->SetTitle("chip");
      h2_i->GetYaxis()->SetTitle("chn");
      h2_i->GetZaxis()->SetRangeUser(0,5);
      h2_i->Draw("colz");
      canvas1->cd(3);
      h2_c1->SetTitle("C1 noise, "+name_);
      h2_c1->GetXaxis()->SetTitle("chip");
      h2_c1->GetYaxis()->SetTitle("chn");
      h2_c1->GetZaxis()->SetRangeUser(0,5);
      h2_c1->Draw("colz");
    
      canvas1->cd(3);
      h2_c2->SetTitle("C2 noise, "+name_);
      h2_c2->GetXaxis()->SetTitle("chip");
      h2_c2->GetYaxis()->SetTitle("chn");
      h2_c2->GetZaxis()->SetRangeUser(0,5);
      h2_c2->Draw("colz");
      canvas1->Write();
      delete canvas1;

    }
  }
  
}


void SummaryPedestal(TString name_="3GeVMIPscan",TString st_gain="highgain"){

  int nlayers=15;
  int nchips=16;
  int nsca=3;
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

  //ofstream fout_ped(TString::Format("../../../pedestals/Pedestal_method2_%s_%s.txt",name_.Data(),st_gain.Data()).Data(),ios::out);
  //fout_ped<<"#pedestal results (method of covariance matrix) : PROTO15"<<endl;
  //fout_ped<<"#layer chip channel";
  //for(int isca=0; isca<15; isca++) fout_ped<<TString::Format(" ped%i ped_rms%i noise_incoherent_ped%i noise_coherent1_ped%i noise_coherent2_ped%i",isca,isca,isca,isca,isca);
  //  fout_ped<<endl;

  for(int layer=0; layer<nlayers; layer++) {

    ofstream fout_ped(TString::Format("../../../pedestals/Pedestal_method2_%s_%s_layer%i.txt",name_.Data(),st_gain.Data(),layer).Data(),ios::out);
    //fout_ped<<"#pedestal results (method of covariance matrix) : PROTO15"<<endl;
    //fout_ped<<"#layer chip channel";
    //    for(int isca=0; isca<15; isca++) fout_ped<<TString::Format(" ped%i ped_rms%i noise_incoherent_ped%i noise_coherent1_ped%i noise_coherent2_ped%i",isca,isca,isca,isca,isca);
    //fout_ped<<endl;
    //TString map="../../../mapping/fev10_chip_channel_x_y_mapping.txt";
    //if(layer==2 || layer==3)  map="../../../mapping/fev11_cob_chip_channel_x_y_mapping.txt";
    //    ReadMap(map,layer);
    
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

	for(int isca=0; isca<nsca; isca++) {
	  sca=TString::Format("%s_sca%i",st_gain.Data(),isca);
	  TFile *file = TFile::Open(TString::Format("Summary_%s_%s.root",name_.Data(),sca.Data()));
	  TH2F * h2_ped=(TH2F*)file->Get(TString::Format("pedestal_layer%i",new_layer_ord[layer]));
	  TH2F * h2_eped=(TH2F*)file->Get(TString::Format("error_pedestal_layer%i",new_layer_ord[layer]));

	  TH2F * h2_i=(TH2F*)file->Get(TString::Format("2d noise-incoherent layer%i",new_layer_ord[layer]));
	  TH2F * h2_c1=(TH2F*)file->Get(TString::Format("2d noise-coherent-1 layer%i",new_layer_ord[layer]));
	  TH2F * h2_c2=(TH2F*)file->Get(TString::Format("2d noise-coherent-2 layer%i",new_layer_ord[layer]));

	  
	  if(h2_ped->GetBinContent(chip+1,j+1)>10) {
	    pedv+=h2_ped->GetBinContent(chip+1,j+1);
	    pedrmsv+=h2_eped->GetBinContent(chip+1,j+1);
	    pedv_[isca]=h2_ped->GetBinContent(chip+1,j+1);
	    pedrmsv_[isca]=h2_eped->GetBinContent(chip+1,j+1);
	    nped++;
	    fitted[isca]=1;

	    if(h2_i->GetBinContent(chip+1,j+1)>0.1 && h2_i->GetBinContent(chip+1,j+1)<7) {
	      pediv+=h2_i->GetBinContent(chip+1,j+1);
	      pedc1v+=h2_c1->GetBinContent(chip+1,j+1);
	      pedc2v+=h2_c2->GetBinContent(chip+1,j+1);
	      
	      pediv_[isca]=h2_i->GetBinContent(chip+1,j+1);
	      pedc1v_[isca]=h2_c1->GetBinContent(chip+1,j+1);
	      pedc2v_[isca]=h2_c2->GetBinContent(chip+1,j+1);
	      nnoise++;
	    }
	  }
	  delete file;//file->Close();

	
      	}//isca
      
	for(int isca=0; isca<15; isca++) {
	  if(fitted[isca]==0 && nped>0 && nnoise>0) {//bad fits for all but with data for average
	    ped->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pedv/nped);
	    ped_rms->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pedrmsv/nped);
	    ped_i->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pediv/nnoise);
	    ped_c1->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pedc1v/nnoise);
	    ped_c2->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pedc2v/nnoise);
	    ped_nofit->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pedv/nped);
	    ped_rms_nofit->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pedrmsv/nped);
	    ped_i_nofit->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pediv/nnoise);
	    ped_c1_nofit->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pedc1v/nnoise);
	    ped_c2_nofit->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pedc2v/nnoise);
	    fout_ped<<pedv/nped<< " " <<0<< " " <<pediv/nnoise<<" "<<pedc1v/nnoise<< " " <<pedc2v/nnoise<<" "; //bad fits for all but with data for average
	  } else {
	    if(fitted[isca]==1 && pediv_[isca]>0 ) {
	      fout_ped<<pedv_[isca]<< " " <<pedrmsv_[isca]<< " " <<pediv_[isca]<< " " <<pedc1v_[isca]<< " " <<pedc2v_[isca]<<" "; //good fits for pedestal and noise
	      ped->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pedv_[isca]);
	      ped_rms->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pedrmsv_[isca]);
	      ped_i->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pediv_[isca]);
	      ped_c1->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pedc1v_[isca]);
	      ped_c2->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pedc2v_[isca]);
	    }  else if(fitted[isca]==1 && pediv_[isca]==0 ) {
	      if(nnoise>0) {
		fout_ped<<pedv_[isca]<< " " <<pedrmsv_[isca]<< " "<<pediv/nnoise<< " " <<pedc1v/nnoise<< " " <<pedc2v/nnoise<<" ";//bad fit of noise, but with data for average
		ped_i->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pediv/nnoise);
		ped_rms->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pedrmsv_[isca]);
		ped_c1->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pedc1v/nnoise);
		ped_c2->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pedc2v/nnoise);
		ped_i_nofit->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pediv/nnoise);
		ped_c1_nofit->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pedc1v/nnoise);
		ped_c2_nofit->Fill(20*new_layer_ord[layer]+chip,100*isca+j,pedc2v/nnoise);
	      } else {
		fout_ped<<pedv_[isca]<< " " <<-5<< " " <<-5<< " " <<-5<< " " <<-5<<" "; //bad fits
	      }
	    } else fout_ped<<pedv_[isca]<< " " <<-10<<" " <<-10<< " " <<-10<< " " <<-10<<" "; //no stat for fits
	  }
	}//isca
      	fout_ped<<endl;
      }//chn

         
    }//chip
  }//layer

  TFile *_file1 = new TFile(TString::Format("../../../pedestals/PedSummary_method2_%s_%s.root",name_.Data(),st_gain.Data()),"RECREATE");
  _file1->cd();
 					     
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

void SummaryPlots(TString runname="Pedestal_run_050575_injection_merged"){

  SummaryPedestal(runname,"highgain");
  //SummaryPedestal(runname,"lowgain");
  Maps(runname,"highgain");
  //Maps(runname,"lowgain");
}

