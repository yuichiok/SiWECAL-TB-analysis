#include "../include/utils.h"

void Plots(TString name_="run_050575_injection_merged_ILCEnergyElectrons",TString st_gain="highgain"){

  int nlayers=15;
  int nchips=16;
  int nsca=3;
  
  ReadPedestalsProtoCovariance("Pedestal_method2_"+name_+"_"+st_gain+".txt");
  
  TFile *_file1 = new TFile(TString::Format("plots/PedSummaryPlots_method2_%s_%s.root",name_.Data(),st_gain.Data()),"RECREATE");
  _file1->cd();
  
  TFile *file0 = TFile::Open(TString::Format("PedSummary_method2_%s_%s.root",name_.Data(),st_gain.Data()));
  TH2F * h2ped=(TH2F*)file0->Get("ped");
  TH2F * h2ped_i=(TH2F*)file0->Get("ped_i");
  TH2F * h2ped_c1=(TH2F*)file0->Get("ped_c1");
  TH2F * h2ped_c2=(TH2F*)file0->Get("ped_c2");

  _file1->cd();
  h2ped->Write();
  h2ped_i->Write();
  h2ped_c1->Write();
  h2ped_c2->Write();
  //h2ped->Save(TString::Format("pedestal_map_%s.png",st_gain.Data()));
  //  h2ped_i->Save(TString::Format("incoherent_map_%s.png",st_gain.Data()));

  delete file0;

  for(int layer=0; layer<nlayers; layer++) {
    cout<<layer<<endl;
    // TString map="../../../mapping/fev10_chip_channel_x_y_mapping.txt";
    // if(layer==2 || layer==3)  map="../../../mapping/fev11_cob_chip_channel_x_y_mapping.txt";
    // ReadMap(map,layer);
    TGraphErrors *g_rms[16];
    TGraphErrors *g_inoise[16];
    TGraphErrors *g_totalnoise[16];
    double minimum=1.01;
    if(st_gain=="lowgain") minimum=0.1;
	
    for(int chip=0; chip<nchips; chip++) {

      double x[15]={0}, ex[15]={0};
      double yrms[15]={0}, eyrms[15]={0};
      double yinoise[15]={0}, eyinoise[15]={0};
      double ytotalnoise[15]={0}, eytotalnoise[15]={0};
    
      for(int isca=0; isca<nsca; isca++) {
	x[isca]=isca;

	double yrms_=0, eyrms_=0, nrms_=0;
	double yinoise_=0, eyinoise_=0, ninoise_=0;
	double ytotalnoise_=0, eytotalnoise_=0, ntotalnoise_=0;

	//average
	for(int j=0; j<64; j++) {
	  if(ped_error_slboard.at(layer).at(chip).at(j).at(isca)>minimum) {
	    yrms_+=ped_error_slboard.at(layer).at(chip).at(j).at(isca);
	    nrms_++;
	  }
	  if(ped_w_i_slboard.at(layer).at(chip).at(j).at(isca)>minimum) {
	    yinoise_+=ped_w_i_slboard.at(layer).at(chip).at(j).at(isca);
	    ninoise_++;
	  }
	  if(ped_w_i_slboard.at(layer).at(chip).at(j).at(isca)>minimum) {
	    ytotalnoise_+=sqrt(
			       pow(ped_w_i_slboard.at(layer).at(chip).at(j).at(isca),2)+
			       pow(ped_w_c1_slboard.at(layer).at(chip).at(j).at(isca),2)+
			       pow(ped_w_c2_slboard.at(layer).at(chip).at(j).at(isca),2));
	    ntotalnoise_++;
	  }
	}//channel

	if(nrms_>0) yrms_/=nrms_;
	else yrms_=0;
	if(ninoise_>0) yinoise_/=ninoise_;
	else yinoise_=0;
	if(ntotalnoise_>0) ytotalnoise_/=ntotalnoise_;
	else ytotalnoise_=0;

	//std dev
	for(int j=0; j<64; j++) {
	  if(ped_error_slboard.at(layer).at(chip).at(j).at(isca)>minimum) {
	    eyrms_+=pow(ped_error_slboard.at(layer).at(chip).at(j).at(isca)-yrms_,2);
	  }
	  if(ped_w_i_slboard.at(layer).at(chip).at(j).at(isca)>minimum) {
	    eyinoise_+=pow(ped_w_i_slboard.at(layer).at(chip).at(j).at(isca)-yinoise_,2);
	  }
	  if(ped_w_i_slboard.at(layer).at(chip).at(j).at(isca)>minimum) {
	    eytotalnoise_+=pow(sqrt(
				    pow(ped_w_i_slboard.at(layer).at(chip).at(j).at(isca),2)+
				    pow(ped_w_c1_slboard.at(layer).at(chip).at(j).at(isca),2)+
				    pow(ped_w_c2_slboard.at(layer).at(chip).at(j).at(isca),2))
			       -
			       ytotalnoise_,2);
	  }
	}//channel

	if(nrms_>0) eyrms_=sqrt(eyrms_/nrms_);
	else eyrms_=0;
	if(ninoise_>0) eyinoise_=sqrt(eyinoise_/ninoise_);
	else eyinoise_=0;
	if(ntotalnoise_>0) eytotalnoise_=sqrt(eytotalnoise_/ntotalnoise_);
	else eytotalnoise_=0;

	yrms[isca]=yrms_; eyrms[isca]=eyrms_;
  	yinoise[isca]=yinoise_; eyinoise[isca]=eyinoise_;
	ytotalnoise[isca]=ytotalnoise_; eytotalnoise[isca]=eytotalnoise_;

      }
      g_rms[chip]=new TGraphErrors(nsca,x,yrms,ex,eyrms);
      g_inoise[chip]=new TGraphErrors(nsca,x,yinoise,ex,eyinoise);
      g_totalnoise[chip]=new TGraphErrors(nsca,x,ytotalnoise,ex,eytotalnoise);

    }//ichip
  

    _file1->cd();
 					     
    TCanvas* canvas0= new TCanvas(TString::Format("NoiseEstimations_%s_%s-layer%i",name_.Data(),st_gain.Data(),layer),TString::Format("NoiseEstimations_%s_%s-layer%i",name_.Data(),st_gain.Data(),layer),1200,1200);   
    canvas0->Divide(2,2);
    canvas0->cd(1);
    TLegend * leg = new TLegend(0.2,0.2,0.8,0.8);
    for(int chip=0; chip<16; chip++) {
      if(chip!=9) {
	g_rms[chip]->SetLineColor(chip+1);
	g_rms[chip]->SetMarkerColor(chip+1);
      } else {
	g_rms[chip]->SetLineColor(1);
	g_rms[chip]->SetMarkerColor(1);
      }
      g_rms[chip]->SetMarkerStyle(20+chip);
      g_rms[chip]->SetMarkerSize(1.);

      
      if ( chip % 2 == 0) g_rms[chip]->SetLineStyle(2);
      else g_rms[chip]->SetLineStyle(1);
      if ( chip % 2 == 0) g_rms[chip]->SetLineWidth(2);
      else g_rms[chip]->SetLineWidth(4);
      if(chip==0) {
	g_rms[chip]->GetYaxis()->SetRangeUser(minimum,5);
	g_rms[chip]->GetYaxis()->SetTitle("ADC");
	g_rms[chip]->GetXaxis()->SetTitle("SCA");
	g_rms[chip]->SetTitle("Pedestal RMS");
	g_rms[chip]->Draw("alp");
      }
      else g_rms[chip]->Draw("lp");
      // g_rms[chip]->Write();
      leg->AddEntry(g_rms[chip],TString::Format("ASIC:%i",chip),"lpe");
    }

    canvas0->cd(2);
    for(int chip=0; chip<16; chip++) {
      if(chip!=9) {
	g_inoise[chip]->SetLineColor(chip+1);
	g_inoise[chip]->SetMarkerColor(chip+1);
      } else {
	g_inoise[chip]->SetLineColor(1);
	g_inoise[chip]->SetMarkerColor(1);
      }
      g_inoise[chip]->SetMarkerStyle(20+chip);
      g_inoise[chip]->SetMarkerSize(1.);

      
      if ( chip % 2 == 0) g_inoise[chip]->SetLineStyle(2);
      else g_inoise[chip]->SetLineStyle(1);
      if ( chip % 2 == 0) g_inoise[chip]->SetLineWidth(2);
      else g_inoise[chip]->SetLineWidth(4);
      if(chip==0) {
	g_inoise[chip]->GetYaxis()->SetRangeUser(minimum,5);
	g_inoise[chip]->GetYaxis()->SetTitle("ADC");
	g_inoise[chip]->GetXaxis()->SetTitle("SCA");
	g_inoise[chip]->SetTitle("Pedestal Incoherent Noise");
	g_inoise[chip]->Draw("alp");
      }
      else g_inoise[chip]->Draw("lp");
      //  g_inoise[chip]->Write();
    }

    canvas0->cd(3);
    for(int chip=0; chip<16; chip++) {
      if(chip!=9) {
	g_totalnoise[chip]->SetLineColor(chip+1);
	g_totalnoise[chip]->SetMarkerColor(chip+1);
      } else {
	g_totalnoise[chip]->SetLineColor(1);
	g_totalnoise[chip]->SetMarkerColor(1);
      }
      g_totalnoise[chip]->SetMarkerStyle(20+chip);
      g_totalnoise[chip]->SetMarkerSize(1.);

      
      if ( chip % 2 == 0) g_totalnoise[chip]->SetLineStyle(2);
      else g_totalnoise[chip]->SetLineStyle(1);
      if ( chip % 2 == 0) g_totalnoise[chip]->SetLineWidth(2);
      else g_totalnoise[chip]->SetLineWidth(4);
      if(chip==0) {
	g_totalnoise[chip]->GetYaxis()->SetRangeUser(minimum,5);
	g_totalnoise[chip]->GetYaxis()->SetTitle("ADC");
	g_totalnoise[chip]->GetXaxis()->SetTitle("SCA");
	g_totalnoise[chip]->SetTitle("Pedestal incoherent+coherent noise");
	g_totalnoise[chip]->Draw("alp");
      }
      else g_totalnoise[chip]->Draw("lp");
      //      g_totalnoise[chip]->Write();
    }
    
    canvas0->cd(4);
    leg->Draw();
    canvas0->Write();
    canvas0->Print(TString::Format("plots/NoiseEstimations_layer%i_%s.png",layer,st_gain.Data()));
    }
}

void PlotsPedestal(){

  TString runs[2]={"run_050571_injection_merged_LowEnergyElectrons","run_050575_injection_merged_ILCEnergyElectrons"};
  TString gain[2]={"highgain","lowgain"};
  
  for(int i=0; i<2; i++) {
    for(int j=0; j<2; j++) {
      Plots(runs[i],gain[j]);      
    }
  }

}
