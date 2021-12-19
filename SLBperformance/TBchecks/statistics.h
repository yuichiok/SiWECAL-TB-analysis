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
#include "../../include/utils3.h"

//#include "../../style/Style.C"
//#include "../../style/Labels.C"

using namespace std;
std::vector<TH2F*> defaultvec;

std::vector<TH2F*> statistics_calc(TString run="Run_ILC_cosmic_test_11222019", TString type="Trigger", int jmax=11){
  TString name="hit";
  if(type!="Triggers") name="event";
  TH2F* hitrate_layer_beam=new TH2F(name+"rate_layer_beam",type +"(coinc) Rates per chip  (kHz) in beam spot ; File  ; Layer ",1000,-0.5,999.5,15,-0.5,14.5);
  TH2F* hitnoiserate_layer_beam=new TH2F(name+"noiserate_layer_beam",type +"(not coinc) Rates per chip (kHz) in beam spot ; File  ; Layer ",1000,-0.5,999.5,15,-0.5,14.5);
  TH2F* hitrate_layer_nobeam=new TH2F(name+"rate_layer_nobeam",type +"(coinc) Rates per chip  (kHz) out-beam spot ; File  ; Layer ",1000,-0.5,999.5,15,-0.5,14.5);
  TH2F* hitnoiserate_layer_nobeam=new TH2F(name+"noiserate_layer_nobeam",type +"(not coinc) Rates per chip (kHz) out-beam spot ; File  ; Layer ",1000,-0.5,999.5,15,-0.5,14.5);
  TH2F* hitrate_layer=new TH2F(name+"rate_layer",type +"(coinc) Rates per layer (kHz) ; File  ; Layer ",1000,-0.5,999.5,15,-0.5,14.5);
  TH2F* hitnoiserate_layer=new TH2F(name+"noiserate_layer",type +"(not coinc) Rates per layer (kHz) ; File  ; Layer ",1000,-0.5,999.5,15,-0.5,14.5);

  TH2F* hitrate_chip=new TH2F(name+"rate_chip",type +"(coinc) Rates per layer (kHz) ; File  ; Layer*20 + chip ",1000,-0.5,999.5,300,-0.5,299.5);
  TH2F* hitnoiserate_chip=new TH2F(name+"noiserate_chip",type +"(not coinc) Rates per layer (kHz) ; File  ;Layer*20 + chip ",1000,-0.5,999.5,300,-0.5,299.5);

  TH2F* cycles_layer=new TH2F("cycles per layer","File ; layer ;cycles ",1000,-0.5,999.5,15,-0.5,14.5);
  TH2F* signal_layer=new TH2F("signal per layer","File ; layer ;signal ",1000,-0.5,999.5,15,-0.5,14.5);
  TH2F* noise_layer=new TH2F("noise per layer","File ; layer ;noise ",1000,-0.5,999.5,15,-0.5,14.5);

  /*
  int arr[] = { 1, 5, 2, 1, 3, 2, 1 };
    int n = sizeof(arr) / sizeof(arr[0]);
    cout << mostFrequent(arr, n);
  */
  int maxhits[10000] ={0};
  int n=0;
  for(int j=0; j<jmax; j++) {
    
    TString filename = TString::Format("../results_proto/%s_000%i.root",run.Data(),j);
    if(j>9) filename = TString::Format("../results_proto/%s_00%i.root",run.Data(),j);
    if(j>99) filename = TString::Format("../results_proto/%s_0%i.root",run.Data(),j);
    if(j>999) filename = TString::Format("../results_proto/%s_%i.root",run.Data(),j);

    TFile *_file0  = TFile::Open(filename);
    TH1F* temp=(TH1F*)_file0->Get(name+"monitoring_norettrains_layer");
    TH1F* temp2=(TH1F*)_file0->Get(name+"monitoring_norettrains");

    for(int layer=0; layer<15; layer++) {
      double cycles=temp->GetBinContent(3,layer+1);
      double hits=temp->GetBinContent(1,layer+1);
      double noise=temp->GetBinContent(2,layer+1);
      signal_layer->SetBinContent(j+1,layer+1,double(hits));
      noise_layer->SetBinContent(j+1,layer+1,double(noise));

      if(cycles>0) {
	hitrate_layer->SetBinContent(j+1,layer+1,double(hits/cycles));
	hitnoiserate_layer->SetBinContent(j+1,layer+1,double(noise/cycles));
      }
      cycles_layer->SetBinContent(j+1,layer+1,cycles);
      int maxhitschip=0;
      int maxchip=-1;
      for(int ichip=0; ichip<16; ichip++) {
	double hits_c=temp2->GetBinContent(1,20*layer+ichip+1);
	double noise_c=temp2->GetBinContent(2,20*layer+ichip+1);

	if(hits_c>maxhitschip) {
	  maxhitschip=hits_c;
	  maxchip=ichip;
	}
	
	if(cycles>0) {
	  hitrate_chip->SetBinContent(j+1,20*layer+ichip+1,double(hits_c/cycles));
	  hitnoiserate_chip->SetBinContent(j+1,20*layer+ichip+1,double(noise_c/cycles));
	}

      }
      maxhits[n]=maxchip;
      n++;

    }
    _file0->Close();
  }

  int chip_beam = mostFrequent(maxhits, n);
  for(int j=0; j<jmax; j++) {
    
    TString filename = TString::Format("../results_proto/%s_000%i.root",run.Data(),j);
    if(j>9) filename = TString::Format("../results_proto/%s_00%i.root",run.Data(),j);
    if(j>99) filename = TString::Format("../results_proto/%s_0%i.root",run.Data(),j);
    if(j>999) filename = TString::Format("../results_proto/%s_%i.root",run.Data(),j);
    
    cout<<filename<<endl;
    TFile *_file0  = TFile::Open(filename);
    TH1F* temp=(TH1F*)_file0->Get(name+"monitoring_norettrains_layer");
    TH1F* temp2=(TH1F*)_file0->Get(name+"monitoring_norettrains");
    
    for(int layer=0; layer<15; layer++) {
      double cycles=temp->GetBinContent(3,layer+1);
      double hits=temp->GetBinContent(1,layer+1);
      double noise=temp->GetBinContent(2,layer+1);
      
      for(int ichip=0; ichip<16; ichip++) {
	double hits_c=temp2->GetBinContent(1,20*layer+ichip+1);
	double noise_c=temp2->GetBinContent(2,20*layer+ichip+1);
	if(cycles>0) {
	  
	  if(ichip==chip_beam) {
	    hitrate_layer_beam->SetBinContent(j+1,layer+1,double(hits_c/cycles));
	    hitnoiserate_layer_beam->SetBinContent(j+1,layer+1,double(noise_c/cycles));
	  } else {
	    hitrate_layer_nobeam->SetBinContent(j+1,layer+1,hitrate_layer_nobeam->GetBinContent(j+1,layer+1)+double(hits_c/cycles)/15.);
	    hitnoiserate_layer_nobeam->SetBinContent(j+1,layer+1,hitnoiserate_layer_nobeam->GetBinContent(j+1,layer+1)+double(noise_c/cycles)/15.);
	  }
	}
      }
    }
    _file0->Close();
  }
  

  std::vector<TH2F*> result;
  result.push_back(hitrate_layer_beam);
  result.push_back(hitnoiserate_layer_beam);
  result.push_back(hitrate_layer_nobeam);
  result.push_back(hitnoiserate_layer_nobeam);
  result.push_back(hitrate_layer);
  result.push_back(hitnoiserate_layer);
  result.push_back(hitrate_chip);
  result.push_back(hitnoiserate_chip);
  result.push_back(cycles_layer);
  result.push_back(signal_layer);
  result.push_back(noise_layer);

  return result;

}

void statistics_plots(TString run="Run_ILC_cosmic_test_11222019", TString type="Trigger", std::vector<TH2F*>results=defaultvec, int jmax=63){

  
  TH2F* hitrate_layer_beam = results.at(0);
  TH2F* hitnoiserate_layer_beam = results.at(1);
  TH2F* hitrate_layer_nobeam = results.at(2);
  TH2F* hitnoiserate_layer_nobeam = results.at(3);
  TH2F* hitrate_layer = results.at(4);
  TH2F* hitnoiserate_layer = results.at(5);

  TH2F* hitrate_chip = results.at(6);
  TH2F* hitnoiserate_chip = results.at(7);
  
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleSize(0.5);

  TCanvas* canvas= new TCanvas("TotalRates_"+type,"TotalRates_"+type,1600,800);   
  canvas->Divide(2,2);
  canvas->cd(1);
  hitrate_layer->GetXaxis()->SetRangeUser(0,jmax);
  hitrate_layer->GetZaxis()->SetRangeUser(0.01,hitrate_layer->GetMaximum()*2);
  hitrate_layer->Draw("colz");

  canvas->cd(2);
  hitnoiserate_layer->GetXaxis()->SetRangeUser(0,jmax);
  hitnoiserate_layer->GetZaxis()->SetRangeUser(0.01,hitrate_layer->GetMaximum());
  hitnoiserate_layer->Draw("colz");

  canvas->cd(3);
  hitrate_chip->GetXaxis()->SetRangeUser(0,jmax);
  hitrate_chip->GetZaxis()->SetRangeUser(0.001,hitrate_chip->GetMaximum()*2);
  hitrate_chip->Draw("colz");

  canvas->cd(4);
  hitnoiserate_chip->GetXaxis()->SetRangeUser(0,jmax);
  hitnoiserate_chip->GetZaxis()->SetRangeUser(0.001,hitrate_chip->GetMaximum());
  hitnoiserate_chip->Draw("colz");

  TCanvas* canvas2= new TCanvas("PerChipRates_"+type,"PerChipRates_"+type,1600,800);   
  canvas2->Divide(2,2);
  canvas2->cd(1);
  hitrate_layer_beam->GetXaxis()->SetRangeUser(0,jmax);
  hitrate_layer_beam->GetZaxis()->SetRangeUser(0.01,hitrate_layer_beam->GetMaximum());
  hitrate_layer_beam->Draw("colz");

  canvas2->cd(2);
  hitnoiserate_layer_beam->GetXaxis()->SetRangeUser(0,jmax);
  hitnoiserate_layer_beam->GetZaxis()->SetRangeUser(0.01,hitrate_layer_beam->GetMaximum());
  hitnoiserate_layer_beam->Draw("colz");

  canvas2->cd(3);
  TH1F *hnew = (TH1F*)hitnoiserate_layer_beam->Clone("hnew");
  hnew->Divide(hitrate_layer_beam);
  hnew->Scale(100);
  hnew->GetXaxis()->SetRangeUser(0,jmax);
  hnew->GetZaxis()->SetRangeUser(0.,20);
  hnew->SetTitle(type +"B/S per chip (IN the beam spot) [%]");
  hnew->Draw("colz");

   
  canvas2->cd(4);
  TH1F *hnew2 = (TH1F*)hitnoiserate_layer_nobeam->Clone("hnew2");
  hnew2->Divide(hitrate_layer_beam);
  hnew2->Scale(100);
  hnew2->GetXaxis()->SetRangeUser(0,jmax);
  hnew2->GetZaxis()->SetRangeUser(0,20);
  hnew2->SetTitle(type +"B/S per chip (OUTSIDE the beam spot) [%]");
  hnew2->Draw("colz");
  
  TString mode="RECREATE";
  if(type!="Trigger") mode="UPDATE";
  TFile *pedfile = new TFile("../results_calib/"+run+".root" , mode);
  pedfile->cd();
  canvas->Write();
  canvas2->Write();

  hitrate_layer->Write();
  hitnoiserate_layer->Write();
  hitnoiserate_chip->Write();
  hitrate_layer_beam->Write();
  hitnoiserate_layer_beam->Write();
  TString name1="hitnoiserate_beam_over_hitratelayer";
  if(type!="Trigger") name1="eventnoiserate_beam_over_eventratelayer";
  hnew->SetName(name1);
  hnew->Write();
  TString name2="hitnoiserate_nobeam_over_hitratelayer";
  if(type!="Trigger") name2="eventnoiserate_nobeam_over_eventratelayer";
  hnew2->SetName(name2);
  hnew2->Write();
  pedfile->Close();


}



void statistics_hits(TString run="Run_ILC_cosmic_test_11222019", int jmax=63){

  std::vector<TH2F*> results = statistics_calc(run,"Trigger",jmax);
  statistics_plots(run,"Trigger",results,jmax); 

}


void statistics_evts(TString run="Run_ILC_cosmic_test_11222019", int jmax=63){

  std::vector<TH2F*> results = statistics_calc(run,"Events",jmax);
  statistics_plots(run,"Events",results,jmax);

  ofstream fout(TString::Format("../results_calib/statistics_%s.txt",run.Data()).Data(),ios::out);
  TH2F* cycles = results.at(8);
  TH2F* signal = results.at(9);
  TH2F* noise = results.at(10);

  float tot_cycles[15]={0};
  float tot_signal[15]={0};
  float tot_noise[15]={0};
  float max_cycles=0;
  fout << run <<" ";
  for(int ilayer=0; ilayer<15; ilayer++) {
    for(int j=0; j<jmax; j++) {
      tot_cycles[ilayer]+=cycles->GetBinContent(j+1,ilayer+1);
      tot_signal[ilayer]+=signal->GetBinContent(j+1,ilayer+1);
      tot_noise[ilayer]+=noise->GetBinContent(j+1,ilayer+1);
    }
    if(tot_cycles[ilayer]>max_cycles) max_cycles=tot_cycles[ilayer];
  }
  fout << "cycles: "<<max_cycles << " "<<endl;

  for(int ilayer=0; ilayer<15; ilayer++) {
    fout << " | s: "<<tot_signal[ilayer]<<" n: "<<tot_noise[ilayer];      
  }
  fout<<endl;
  
}

