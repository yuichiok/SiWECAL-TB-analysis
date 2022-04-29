#include "TFile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TH1.h"
#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TString.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "TLatex.h"
#include "../../include/utils.h"

//#include "../../style/Style.C"
//#include "../../style/Labels.C"

using namespace std;

int nchips=16;
int nlayers=15;


std::vector<double> resultfit (TH1F* hmips, TString gain) {
  Double_t fr[2];
  Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
  if(gain=="high") {
    pllo[0]=0.1; pllo[1]=14; pllo[2]=1.0; pllo[3]=0.5;
    plhi[0]=20.0; plhi[1]=220.0; plhi[2]=100000000.0; plhi[3]=10.0;
  } else {
    pllo[0]=0.1; pllo[1]=2; pllo[2]=1.0; pllo[3]=0.05;
    plhi[0]=5; plhi[1]=15.0; plhi[2]=100000000.0; plhi[3]=2.0;
  }
  
  Double_t chisqr;
  Int_t    ndf;
  
  // if(gain=="high") hmips->Rebin(2);
  
  if(gain=="high")  hmips->GetXaxis()->SetRangeUser(2, 500);
  else hmips->GetXaxis()->SetRangeUser(2., 50);
  
    
    if(gain=="high") {
      fr[0]=TMath::Max(hmips->GetMean()-0.8*hmips->GetRMS(),5.);
      hmips->GetXaxis()->SetRangeUser(fr[0], 150);
      fr[1]=hmips->GetMean()*1.3;
    } else {
      fr[0]=TMath::Max(hmips->GetMean()-2.*hmips->GetRMS(),2.);
      hmips->GetXaxis()->SetRangeUser(fr[0], 25.);
      fr[1]=hmips->GetMean()*1.1 + 2.*hmips->GetRMS();
    }
    
    TF1 *fitlandau= new TF1("fitlandau","landau",fr[0],fr[1]);
    hmips->Fit("fitlandau","QRE");
    hmips->GetXaxis()->SetRangeUser(fr[0], fr[1]);
    
    sv[0] = TMath::Min(TMath::Max(fitlandau->GetParameter(2),pllo[0]),plhi[0]);
    sv[1] = TMath::Min(TMath::Max(fitlandau->GetParameter(1),pllo[1]),plhi[1]);
    sv[2] = hmips->Integral("width") * 1.2;
    sv[3] = TMath::Min(TMath::Max(hmips->GetRMS()* 0.1,pllo[3]),plhi[3]);
    
    TF1 *fitsnr_hmips=langaufit(hmips,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
    for (int k = 0; k < 4; k++) {
      if ((sv[k] >= pllo[k]) && (sv[k] <= plhi[k])) continue;
      std::cout << "Langaus fit parameter " << k << " has starting value " << sv[k]
		<<  " outside the range [" << pllo[k] << ", " << plhi[k] << "]!" << std::endl;
    }
    double mpv=fitsnr_hmips->GetParameter(1);
    double empv=fitsnr_hmips->GetParError(1);
    double wmpv=fitsnr_hmips->GetParameter(0);
    double chi2ndf=0;
    if(ndf>0) chi2ndf=chisqr/ndf;

    std::vector<double> result;
    result.push_back(mpv);
    result.push_back(empv);
    result.push_back(wmpv);
    result.push_back(chi2ndf);
    return result;

}
  
void pedanalysis(TString run="PedestalMIP_3GeVMIPscan", TString gain="low"){
  
  gROOT->Reset();
  //SetIrlesStyle();
  //  gROOT->LoadMacro("Labels.C");
  gStyle->SetOptFit(0); 
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.2);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetMarkerSize(1.2);


  TFile *_file0 = TFile::Open(TString::Format("../results_calib/PedestalMIP_%s.root",run.Data()));

  TH1F* width_layer[15];
  for(int layer=0; layer<nlayers; layer++) width_layer[layer] =new TH1F(TString::Format("width_layer%i",layer),"widths of pedestal ; ADC  ",300,0,6);

  TH2F* ped_all=new TH2F("ped_all","pedestal pos. ; Layer*20+Chip  ; SCA*100 +chn",300,-0.5,299.5,1500,-0.5,1499.5);
  TH2F* width_all=new TH2F("width_all","width of pedestal ; Layer*20+Chip  ; SCA*100 +chn",300,-0.5,299.5,1500,-0.5,1499.5);

  TH2F* ped_morethan1peak=new TH2F("ped_morethan1peak","pedestal pos. ; Layer*20+Chip  ; SCA*100 +chn",300,-0.5,299.5,1500,-0.5,1499.5);
  TH2F* width_morethan1peak=new TH2F("width_morethan1peak","width of pedestal ; Layer*20+Chip  ; SCA*100 +chn",300,-0.5,299.5,1500,-0.5,1499.5);

  TH2F* ped_nofit=new TH2F("ped_nofit","pedestal pos. ; Layer*20+Chip  ; SCA*100 +chn",300,-0.5,299.5,1500,-0.5,1499.5);
  TH2F* width_nofit=new TH2F("width_nofit","width of pedestal ; Layer*20+Chip  ; SCA*100 +chn",300,-0.5,299.5,1500,-0.5,1499.5);

  ofstream fout_ped(TString::Format("../../pedestals/Pedestal_%s_%sgain.txt",run.Data(),gain.Data()).Data(),ios::out);
  fout_ped<<"#pedestal results (fit to a gaussian) remove channels/sca with two pedestals peaks from the analysis : PROTO15"<<endl;
  fout_ped<<"#layer chip channel";
  for(int isca=0; isca<15; isca++) fout_ped<<TString::Format(" ped_mean%i ped_error%i ped_width%i",isca,isca,isca);
  fout_ped<<endl;

  for(int layer=0; layer<nlayers; layer++) {
   
    // Comparing nbr entries in tag or not tag // GetWidth and Mean
    //TString map="../../mapping/fev10_chip_channel_x_y_mapping.txt";
    //  if(layer==0 || layer==2)  map="../../mapping/fev11_cob_chip_channel_x_y_mapping.txt";
    //  ReadMap(map);

     for(int i=0;i<nchips;i++){
    
      double nch=0;
      float avmean=0, avwidth=0;
    
      for(int j=0; j<64; j++) {
	fout_ped << layer<<" "<<i <<" " <<j<< " ";

	double ped[15]={0};
	double pederror[15]={0};
	double pedwidth[15]={0};
	double pedwidtherror[15]={0};


      	for(int isca=0; isca<15; isca++) {
	  TH1F* temp=(TH1F*)_file0->Get(TString::Format("layer_%i/ped_%s_chip%i_chn%i_sca%i",layer,gain.Data(),i,j,isca));
        
	  if(temp->GetEntries()>50) {
	    TSpectrum *s = new TSpectrum();
	    int npeaks = s->Search(temp,2,"",0.8);
	    if(npeaks == 1) {
	      Double_t *mean_peak=s->GetPositionX();
	      Double_t *mean_high=s->GetPositionY();
	      double mean_peak_higher=0;
	      double mean_high_higher=0;
	      int npeak_max=0;
              for(int ipeak=0; ipeak<npeaks; ipeak ++) {
                if(mean_high[ipeak]>mean_high_higher && mean_peak[ipeak]>180) {
                  mean_high_higher=mean_high[ipeak];
                  mean_peak_higher=mean_peak[ipeak];
                  npeak_max=ipeak;
                }
              }
	      	   
	      float mean_start_fit=mean_peak_higher;
	      
	      TF1 *f0 = new TF1("f0","gaus",mean_peak_higher-8,mean_peak_higher+8);
	      TF1 *f1 = new TF1("f1","gaus",mean_peak_higher-16,mean_peak_higher+16);
	      temp->Fit("f0","RQE");
	      temp->Fit("f1","RQE");

	      if(f0->GetParameter(2)<f1->GetParameter(2))  {
		ped[isca]=f0->GetParameter(1);
		pedwidth[isca]=f0->GetParameter(2);
		pederror[isca]=f0->GetParError(1);
		pedwidtherror[isca]=f0->GetParError(2);
	      }  else  {
		ped[isca]=f1->GetParameter(1);
		pedwidth[isca]=f1->GetParameter(2);
		pederror[isca]=f1->GetParError(1);
		pedwidtherror[isca]=f1->GetParError(2);
	      }
	      
	      delete s;
	    }
	  } else pedwidth[isca]=-1;
	  delete temp;
	}//isca

	double ped_av=weigthedAv(ped, pederror, 15);
	double pedwidth_av=weigthedAv(pedwidth, pedwidtherror, 15);
	for(int isca=0; isca<15; isca++) {  
	  if(pedwidth[isca]>0.5 && pedwidth[isca]<8 ) {
	    fout_ped<<ped[isca]<< " " << pederror[isca]<<" "<<pedwidth[isca]<<" "; //good fits
	    ped_all->Fill(20*layer+i,100*isca+j,ped[isca]);
	    width_all->Fill(20*layer+i,100*isca+j,pedwidth[isca]);
	    width_layer[layer]->Fill(pedwidth[isca]);
	  } else  {
	    if(pedwidth[isca]>-0.5) {
	      fout_ped<<ped_av<< " " << -5<<" "<<pedwidth_av<<" "; //bad fits
	      ped_all->Fill(20*layer+i,100*isca+j,ped_av);
	      width_all->Fill(20*layer+i,100*isca+j,pedwidth_av);
	      ped_morethan1peak->Fill(20*layer+i,100*isca+j,ped_av);
	      width_morethan1peak->Fill(20*layer+i,100*isca+j,pedwidth_av);
	    }
	    if(pedwidth[isca]<0) {
	      fout_ped<<ped_av<< " " << -10<<" "<<pedwidth_av<<" "; //no fits
	      ped_all->Fill(20*layer+i,100*isca+j,ped_av);
	      width_all->Fill(20*layer+i,100*isca+j,pedwidth_av);
	      ped_nofit->Fill(20*layer+i,100*isca+j,ped_av);
	      width_nofit->Fill(20*layer+i,100*isca+j,pedwidth_av);
	    }
	  }
	}
	fout_ped<<endl;
      }//j
    }
  }


  TFile *_file1 = new TFile(TString::Format("../../pedestals/%s_%sgain_PedSummary.root",run.Data(),gain.Data()),"RECREATE");
  _file1->cd();

 for(int layer=0; layer<nlayers; layer++) width_layer[layer]->Write();
  
  TCanvas* canvas0= new TCanvas("PedAna_Summary","Pedestal Summary",1800,400);   
  canvas0->Divide(2,1);
  canvas0->cd(1);
  ped_all->GetZaxis()->SetRangeUser(180,300);
  ped_all->Draw("colz");
  canvas0->cd(2);
  width_all->GetZaxis()->SetRangeUser(1.0,6.0);
  if(gain=="low") width_all->GetZaxis()->SetRangeUser(0.5,3.0);
  width_all->Draw("colz");
  // canvas0->Print(TString::Format("../../pedestals/%s_%sgain_PedSummary.png",run.Data(),gain.Data()));
  ped_all->Write();
  width_all->Write();
  canvas0->Write();

  TCanvas* canvas1= new TCanvas("PedAna_Summary_nofit","Pedestal Summary - nofit",1800,400);   
  canvas1->Divide(2,1);
  canvas1->cd(1);
  ped_nofit->GetZaxis()->SetRangeUser(180,300);
  ped_nofit->Draw("colz");
  canvas1->cd(2);
  width_nofit->GetZaxis()->SetRangeUser(1.0,6.0);
  if(gain=="low") width_nofit->GetZaxis()->SetRangeUser(0.5,3.0);
  width_nofit->Draw("colz");
  // canvas0->Print(TString::Format("../../pedestals/%s_%sgain_PedSummary.png",run.Data(),gain.Data()));
  ped_nofit->Write();
  width_nofit->Write();
  canvas1->Write();

  TCanvas* canvas2= new TCanvas("PedAna_Summary_morethan1peak","Pedestal Summary - morethan1peak",1800,400);   
  canvas2->Divide(2,1);
  canvas2->cd(1);
  ped_morethan1peak->GetZaxis()->SetRangeUser(180,300);
  ped_morethan1peak->Draw("colz");
  canvas2->cd(2);
  width_morethan1peak->GetZaxis()->SetRangeUser(1.0,6.0);
  if(gain=="low") width_morethan1peak->GetZaxis()->SetRangeUser(0.5,3.0);
  width_morethan1peak->Draw("colz");
  // canvas0->Print(TString::Format("../../pedestals/%s_%sgain_PedSummary.png",run.Data(),gain.Data()));
  ped_morethan1peak->Write();
  width_morethan1peak->Write();
  canvas2->Write();

}


void mipanalysis_summary(TString run="3GeVMIPscan", TString gain="high", int pedestal_mode=0, TString run_pedestal=""){

  if(run_pedestal=="") run_pedestal=run;
  // pedestal_mode==0 --> no subtraction
  // pedestal_mode==1 --> on-the-fly subtraction
  // pedestal_mode==2 --> subtraction from covariance file

  gROOT->Reset();
  //SetIrlesStyle();
  //  gROOT->LoadMacro("Labels.C");
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.2);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetMarkerSize(1.2);


  TH2F* MIPM_all=new TH2F("MIPM_all","average of MPVs ; Layer  ; Chip; ADC",15,-0.5,14.5,16,-0.5,15.5);
  TH2F* MIPN_all=new TH2F("MIPN_all","channels fitted ; Layer  ; Chip; ADC",15,-0.5,14.5,16,-0.5,15.5);
  TH2F* MIPrms_all=new TH2F("MIPrms_all","rms  of MPVs ; Layer  ; Chip; ADC",15,-0.5,14.5,16,-0.5,15.5);

  TH2F* MIPM2_all=new TH2F("MIPM2_all","MPV fitting all channels in one fit ; Layer  ; Chip; ADC",15,-0.5,14.5,16,-0.5,15.5);
  TH2F* MIPN2_all=new TH2F("MIPN2_all","total entries per chip ; Layer  ; Chip; ADC",15,-0.5,14.5,16,-0.5,15.5);
  TH2F* MIPrms2_all=new TH2F("MIPrms2_all","rms  of MPVs ; Layer  ; Chip; ADC",15,-0.5,14.5,16,-0.5,15.5);

  TH2F* mpv_all=new TH2F("mpv_all","mpv pos. ; Layer*20+Chip  ; chn",300,-0.5,299.5,64,-0.5,63.5);
  TH2F* width_all=new TH2F("width_all","width of Landau ; Layer*20+Chip  ; chn",300,-0.5,299.5,64,-0.5,63.5);
  TH2F* nentries_all=new TH2F("nentries_all","N entries; Layer*20+Chip  ; chn",300,-0.5,299.5,64,-0.5,63.5);

  ofstream fout_mip(TString::Format("../../mip_calib/MIP_pedestalsubmode%i_%s_%sgain.txt",pedestal_mode,run.Data(),gain.Data()).Data(),ios::out);
  fout_mip<<"#mip results PROTO15-TB2022-03"<<endl;
  fout_mip<<"#layer chip channel mpv empv widthmpv chi2ndf nentries"<<endl;

  if(pedestal_mode==2) {
    TString pedfilename=TString::Format("../../pedestals/Pedestal_method2_%s_%sgain.txt",run_pedestal.Data(),gain.Data());
    ReadPedestalsProtoCovariance(pedfilename.Data());
    cout<<pedfilename<<endl;
  }

  TFile *_file1 = new TFile(TString::Format("../../mip_calib/MIPSummary_pedestalsubmode%i_%s_%sgain.root",pedestal_mode,run.Data(),gain.Data()),"RECREATE");

  TDirectory *cdhisto[15];
  for(int ilayer=0; ilayer<nlayers; ilayer++) {
    cdhisto[ilayer] = _file1->mkdir(TString::Format("layer_%i",ilayer));
  }
  
  TH2F *mpv_layer[15];
  TH2F *entries_layer[15];
  TH2F *mpv_layer_xy[15];
  TH2F *entries_layer_xy[15];

  int new_layer_ord[15]={5,6,0,1,4,2,3,7,8,9,10,11,12,13,14};
  
  for(int layer=0; layer<nlayers; layer++) {

    
    mpv_layer_xy[layer]= new TH2F(TString::Format("mpv_layer%i_xy",layer),TString::Format("mpv_layer%i_xy; x; y; Fitted MPV",layer),32,-90,90,32,-90,90);
    entries_layer_xy[layer]= new TH2F(TString::Format("entries_layer%i_xy",layer),TString::Format("entries_layer%i_xy; x; y; Fitted MPV",layer),32,-90,90,32,-90,90);

    mpv_layer[layer]= new TH2F(TString::Format("mpv_layer%i",layer),TString::Format("mpv_layer%i; chip; chn; Fitted MPV",layer),16,-0.5,15.5,64,-0.5,63.5);
    entries_layer[layer]= new TH2F(TString::Format("entries_layer%i",layer),TString::Format("entries_layer%i; chip; chn; Fitted MPV",layer),16,-0.5,15.5,64,-0.5,63.5);

    TString map="../../mapping/fev10_chip_channel_x_y_mapping.txt";
    if(layer==5 || layer==6)  map="../../mapping/fev11_cob_chip_channel_x_y_mapping.txt";
    ReadMap(map,layer);
 

    float avmpv[16]={0};
    float nmpv[16]={0};
    float rmsmpv[16]={0};
    float avmpv2[16]={0};
    float nmpv2[16]={0};

    cout<<"   LAYER : "<<layer<<endl;
    TFile *_file0 = TFile::Open(TString::Format("../results_calib/PedestalMIP_%s.root",run.Data()));

    TH1F *hmip_chip[16];

    for(int i=0;i<nchips;i++){
    cout<<"       -  chip : "<<i<<endl;

      hmip_chip[i] = new TH1F(TString::Format("mip_%s_layer%i_chip%i",gain.Data(),layer,i),TString::Format("mip_%s_layer%i_chip%i",gain.Data(),layer,i),900,-100.5,799.5);//600,-199.5,400.5);

      for(int j=0; j<64; j++) {
	TH1F *hmips = new TH1F("hmips","hmips",900,-100.5,799.5);//600,-199.5,400.5);
	hmips->Sumw2();
	for(int isca=0; isca<15; isca++) {
	  _file0->cd();
	  // TCanvas* canvashmips=new TCanvas;
	  
	  TH1F *temp=(TH1F*)_file0->Get(TString::Format("layer_%i/mip_%s_chip%i_chn%i_sca%i",layer,gain.Data(),i,j,isca));
	  TH1F *temp2=(TH1F*)_file0->Get(TString::Format("layer_%i/ped_%s_chip%i_chn%i_sca%i",layer,gain.Data(),i,j,isca));

	  if(temp==NULL || temp2==NULL) continue;

    //	  cout<<layer<<" "<<i<<" "<<j<<" "<<isca<<" "<<temp->GetEntries()<<endl;
	  double ped_mean=0;
	  if(pedestal_mode==1) {
	    temp2->GetXaxis()->SetRangeUser(temp2->GetMean()-20,temp2->GetMean()+20);
	    ped_mean=temp2->GetMean();
	    //cout<<ped_mean<<endl;
	  }
	  if(pedestal_mode==2) ped_mean=ped_mean_slboard.at(new_layer_ord[layer]).at(i).at(j).at(isca);
	  for (int k=0;k<900;k++) {
	    double y = temp->GetBinContent(k);
	    if(y>0 && pedestal_mode==0) hmips->Fill(int(temp->GetXaxis()->GetBinCenter(k)),y);
	    if(y>0 && pedestal_mode>0 && ped_mean>0) {
	      hmips->Fill(int(temp->GetXaxis()->GetBinCenter(k)-ped_mean),y);
	    }

	  }
	  delete temp;
	  delete temp2;
	  //  delete canvashmips;
	}
	
	for(int k=0; k<900; k++) {
	  hmips->SetBinError(k,sqrt(hmips->GetBinContent(k)));
	  hmip_chip[i]->SetBinContent(k,hmip_chip[i]->GetBinContent(k)+hmips->GetBinContent(k));
	}

	MIPN2_all->Fill(layer,i,hmip_chip[i]->Integral());

	if(hmips->Integral()>500) {
	  std::vector<double> result=resultfit(hmips,gain);
	  double mpv=result.at(0);
	  double empv=result.at(1);
	  double wmpv=result.at(2);
	  double chi2ndf=result.at(3);
	  
	  //	MIPN->Fill(map_pointX[i][j],map_pointY[i][j],hmips->Integral());

	  float mpvmin=20.;
	  if(gain=="low") mpvmin=0.8;
	  float mpvmax=220;
	  if(gain=="low") mpvmax=15;
	  if(chi2ndf>0){// && mpv>mpvmin && mpv<mpvmax) {
	    avmpv[i]+=mpv/empv;
	    rmsmpv[i]+=1./(empv);
	    nmpv[i]++;
	    fout_mip<<layer<<" "<<i<<" "<<j<<" "<<mpv<<" "<<empv<<" "<<wmpv<<" "<<chi2ndf<<" "<<hmips->Integral()<<"\n";
	    mpv_all->Fill(20*layer+i,j,mpv);
	    width_all->Fill(20*layer+i,j,wmpv);
	    nentries_all->Fill(20*layer+i,j,hmips->Integral());
	    mpv_layer[layer]->Fill(i,j,mpv);
	    entries_layer[layer]->Fill(i,j,hmips->Integral());
	    mpv_layer_xy[layer]->Fill(double(map_pointX[layer][i][j]),double(map_pointY[layer][i][j]),mpv);
	    entries_layer_xy[layer]->Fill(double(map_pointX[layer][i][j]),double(map_pointY[layer][i][j]),hmips->Integral());
	  } else fout_mip<<layer<<" "<<i<<" "<<j<<" "<<0<<" "<<-5<<" "<<0<<" "<<0<<" "<<0<<"\n";

	} else {
	  fout_mip<<layer<<" "<<i<<" "<<j<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<"\n";
	}
	_file1->cd();
	cdhisto[layer]->cd();
        hmips->SetName(TString::Format("mip_%s_layer%i_chip%i_chn%i",gain.Data(),layer,i,j));
        hmips->SetTitle(TString::Format("mip_%s_layer%i_chip%i_chn%i",gain.Data(),layer,i,j));
	hmips->Write();

	delete hmips;
      }
      cdhisto[layer]->cd();
      for(int k=0; k<900; k++) {
	hmip_chip[i]->SetBinError(k,sqrt(hmip_chip[i]->GetBinContent(k)));
      }
      if(hmip_chip[i]->Integral()>2000) {
	std::vector<double> resultchip=resultfit(hmip_chip[i],gain);
	avmpv2[i]=resultchip.at(0);
      } else avmpv2[i]=0;
      hmip_chip[i]->Write();

    }
    for(int i=0;i<nchips;i++){
      if(nmpv[i]>0. ) {
    	avmpv[i]/=rmsmpv[i];//nmpv[i];
    	MIPM_all->Fill(layer,i,avmpv[i]);
    	//MIPrms_all->Fill(layer,i,rmsmpv[i]);
    	MIPN_all->Fill(layer,i,nmpv[i]);
      }
      MIPM2_all->Fill(layer,i,avmpv2[i]);
        
    }
    _file1->cd();

    int chip_vs_canvas[16]={6,4,2,0,7,5,3,1,14,12,10,8,15,13,11,9};
    TCanvas* canvassummary= new TCanvas(TString::Format("MIPSummary_layer%i",layer),TString::Format("MIPSummary_layer%i",layer),800,800);
    canvassummary->Divide(4,4);
    for(int i=0; i<nchips; i++) {
      canvassummary->cd(i+1);

      hmip_chip[chip_vs_canvas[i]]->GetXaxis()->SetRangeUser(0,300);

      hmip_chip[chip_vs_canvas[i]]->Draw();
    }
    canvassummary->Write();
    delete canvassummary;

    TCanvas* canvassummary2= new TCanvas(TString::Format("MIPSummary2_layer%i",layer),TString::Format("MIPSummary2_layer%i",layer),800,800);
    canvassummary2->Divide(2,2);
    canvassummary2->cd(1);
    mpv_layer_xy[layer]->Draw("colz");
    canvassummary2->cd(2);
    entries_layer_xy[layer]->Draw("colz");
    canvassummary2->cd(3);
    mpv_layer[layer]->Draw("colz");
    canvassummary2->cd(4);
    entries_layer[layer]->Draw("colz");
      
    canvassummary2->Write();
    delete canvassummary2;

    delete _file0;
  }

  gStyle->SetPalette(kInvertedDarkBodyRadiator);


  TCanvas* canvassummary= new TCanvas("MIPAna","MIPAna",1400,600);
  canvassummary->Divide(2,2);
    
  canvassummary->cd(1);
  // MIPM_all->GetXaxis()->SetTitle("x");
  // MIPM_all->GetYaxis()->SetTitle("y");
  if(gain=="high")  MIPM_all->GetZaxis()->SetRangeUser(0,100);
  else MIPM_all->GetZaxis()->SetRangeUser(0,15);
  MIPM_all->Draw("colz");

  canvassummary->cd(2);
  // MIPM_all->GetXaxis()->SetTitle("x");
  // MIPM_all->GetYaxis()->SetTitle("y");
  if(gain=="high")  MIPM2_all->GetZaxis()->SetRangeUser(0,100);
  else MIPM2_all->GetZaxis()->SetRangeUser(0,15);
  MIPM2_all->Draw("colz");

  /*  canvassummary->cd(2);
  // MIPrms_all->GetXaxis()->SetTitle("x");
  // MIPrms_all->GetYaxis()->SetTitle("y");
  if(gain=="high") MIPrms_all->GetZaxis()->SetRangeUser(0,10);
  else MIPrms_all->GetZaxis()->SetRangeUser(0,5);
  MIPrms_all->Draw("COLZ");*/
    
  canvassummary->cd(3);
  //gPad->SetLogz();
  // MIPN_all->GetXaxis()->SetTitle("x");
  // MIPN_all->GetYaxis()->SetTitle("y");
  MIPN_all->GetZaxis()->SetRangeUser(0,64.5);
  MIPN_all->Draw("COLZ");

  canvassummary->cd(4);
  //MIPN_all->GetZaxis()->SetRangeUser(0,64.5);
  MIPN2_all->Draw("COLZ");

  _file1->cd();
  canvassummary->Write();
  //  hmips->Draw();
  mpv_all->Write();
  width_all->Write();
  nentries_all->Write();
  //_file1->Write();
  //_file1->Close();
  //canvassummary->Print(run+"_MIPAna.root");
  //  delete canvassummary;
}



