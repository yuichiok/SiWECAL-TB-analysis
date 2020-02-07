#include "TFile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TH1.h"
#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TString.h"
//#include "../../style/Style.C"
//#include "../../style/Labels.C"

using namespace std;
Float_t map_pointX[16][64];
Float_t map_pointY[16][64];

void ReadMap(TString filename="/home/irles/WorkAreaECAL/2019/SiWECAL-TB-analysis/mapping/fev10_chip_channel_x_y_mapping.txt")   
{

  std::ifstream reading_file(filename);
  if(!reading_file){
    cout<<" dameyo - damedame"<<endl;
  }
  for(int i=0; i<16; i++) {
    for(int j=0; j<64; j++) {
      map_pointX[i][j] = -1000.;
      map_pointY[i][j] = -1000.;
    }
  }

  Int_t tmp_chip = 0,tmp_channel = 0;
  Float_t tmp_x0 = 0 ,tmp_y0 = 0 , tmp_x = 0 , tmp_y = 0 ;
  TString tmpst;
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst ;
  while(reading_file){
    reading_file >> tmp_chip >> tmp_x0 >> tmp_y0 >> tmp_channel >> tmp_x >> tmp_y ;
    map_pointX[tmp_chip][tmp_channel] = -tmp_x ;
    map_pointY[tmp_chip][tmp_channel] = -tmp_y ;
  }

}


void pedanalysis0(TString run="Run_ILC_cosmic_test_11222019", int slboard=2){
  
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
  
  // Comparing nbr entries in tag or not tag // GetWidth and Mean
  TString map="/home/irles/WorkAreaECAL/2019/SiWECAL-TB-analysis/mapping/fev10_chip_channel_x_y_mapping.txt";
  if(slboard==0 || slboard==2)  map="/home/irles/WorkAreaECAL/2019/SiWECAL-TB-analysis/mapping/fev11_cob_chip_channel_x_y_mapping.txt";
  ReadMap(map);

  
  TFile *_file0 = TFile::Open(TString::Format("../results_pedestal/Pedestal_SLB_%i_%s.root",slboard,run.Data()));

 
  TH2F* PedW[15];
  TH2F* PedM[15];
  TH2F* PedNPeaks[15];
  TH2F* PedRMS[15];

  for(int i=0; i<15; i++) {
    PedW[i]=new TH2F(TString::Format("PedW_sca%i",i),TString::Format("PedW_sca%i",i),32,-90,90,32,-90,90);
    PedM[i]=new TH2F(TString::Format("PedM_sca%i",i),TString::Format("PedM_sca%i",i),32,-90,90,32,-90,90);

    PedRMS[i]=new TH2F(TString::Format("PedRMS_sca%i",i),TString::Format("PedRMS_sca%i",i),32,-90,90,32,-90,90);
    PedNPeaks[i]=new TH2F(TString::Format("PedNPeaks_sca%i",i),TString::Format("PedNPeaks_sca%i",i),32,-90,90,32,-90,90);
  }

  TH2F* PedW_sca;
  TH2F* PedM_sca;
  PedW_sca=new TH2F("PedW_sca","PedW_sca",16,0,15,15,0,14);
  PedM_sca=new TH2F("PedM_sca","PedM_sca",16,0,15,15,0,14);

  for(int n=0; n<15;n++) {
    for(int i=0;i<4;i++){

      double nch=0;
      float avmean=0, avwidth=0;
      
      for(int j=0; j<64; j++) {
	//GetGoodEntries
	TH1F *temp=(TH1F*)_file0->Get(TString::Format("ped_chip%i_chn%i_sca%i",i,j,n));
	
	if(temp->GetEntries()> 200 ){ //max_entries/2 ) {
          TSpectrum *s = new TSpectrum();
          int npeaks = s->Search(temp,2,"",0.8); 
          if(npeaks > 0) {
            Double_t *mean_peak=s->GetPositionX();
            Double_t *mean_high=s->GetPositionY();
            double mean_peak_higher=0;
            double mean_high_higher=0;
            int npeak_max=0;
            for(int ipeak=0; ipeak<npeaks; ipeak ++) {
              if(mean_high[ipeak]>mean_high_higher && mean_high[ipeak]>50) {
                mean_high_higher=mean_high[ipeak];
                mean_peak_higher=mean_peak[ipeak];
                npeak_max=ipeak;
              }
            }
	    PedNPeaks[n]->Fill(map_pointX[i][j],map_pointY[i][j],npeaks);
	    
            if(npeaks ==1 ) {
	      nch++;

              Double_t *mean_peak=s->GetPositionX();
              mean_peak[0]=mean_peak_higher;
              
              TF1 *f0 = new TF1("f0","gaus",mean_peak[0]-temp->GetRMS()/2,mean_peak[0]+temp->GetRMS()/2);
              temp->Fit("f0","RQNOC");
              
              TF1 *f1 = new TF1("f1","gaus",f0->GetParameter(1)-f0->GetParameter(2)/2,f0->GetParameter(1)+f0->GetParameter(2)/2);
              temp->Fit("f1","RQME");
	      avmean+=f1->GetParameter(1);
	      avwidth+=f1->GetParameter(2);
	      PedW[n]->Fill(map_pointX[i][j],map_pointY[i][j],f1->GetParameter(2));
	      PedRMS[n]->Fill(map_pointX[i][j],map_pointY[i][j],temp->GetRMS());
	      PedM[n]->Fill(map_pointX[i][j],map_pointY[i][j],f1->GetParameter(1));
            }
	  }
       }

      }
      
      if(nch>20) {
	PedW_sca->Fill(i,n,avwidth/nch);
	PedM_sca->Fill(i,n,avmean/nch);
      }
    }

    
  }

  //  _file0 ->Close();

  TCanvas* canvas= new TCanvas(TString::Format("PedAna_%s_SLB_%i",run.Data(),slboard),TString::Format("PedAna_%s_SLB_%i",run.Data(),slboard),1600,1600);   
  canvas->Divide(4,4);
  
  for(int i=0; i<4; i++) {
    canvas->cd(1+i);
    PedM[i]->GetXaxis()->SetTitle("x");
    PedM[i]->GetYaxis()->SetTitle("y");
    PedM[i]->GetZaxis()->SetRangeUser(200,400);
    PedM[i]->Draw("colz");

    canvas->cd(5+i);
    PedW[i]->GetXaxis()->SetTitle("x");
    PedW[i]->GetYaxis()->SetTitle("y");
    PedW[i]->GetZaxis()->SetRangeUser(2,12);
    PedW[i]->Draw("COLZ");
    
    canvas->cd(9+i);
    PedNPeaks[i]->GetXaxis()->SetTitle("x");
    PedNPeaks[i]->GetYaxis()->SetTitle("y");
    PedNPeaks[i]->GetZaxis()->SetRangeUser(0,5);
    PedNPeaks[i]->Draw("COLZ");

    canvas->cd(13+i);
    PedRMS[i]->GetXaxis()->SetTitle("x");
    PedRMS[i]->GetYaxis()->SetTitle("y");
    PedRMS[i]->GetZaxis()->SetRangeUser(2,12);
    PedRMS[i]->Draw("COLZ");
    
  }
  canvas->Print(TString::Format("plots/PedAna_%s_SLB_%i.eps",run.Data(),slboard));

  TCanvas* canvas2= new TCanvas(TString::Format("PedAna2_%s_SLB_%i",run.Data(),slboard),TString::Format("PedAna2_%s_SLB_%i",run.Data(),slboard),1200,600);   
  canvas2->Divide(2,1);
  
  canvas2->cd(1);
  PedM_sca->GetXaxis()->SetTitle("x");
  PedM_sca->GetYaxis()->SetTitle("y");
  PedM_sca->GetZaxis()->SetRangeUser(200,400);
  PedM_sca->Draw("colz");

  canvas2->cd(2);
  PedW_sca->GetXaxis()->SetTitle("x");
  PedW_sca->GetYaxis()->SetTitle("y");
  PedW_sca->GetZaxis()->SetRangeUser(2,12);
  PedW_sca->Draw("colz");
   
  canvas2->Print(TString::Format("plots/PedAna2_%s_SLB_%i.eps",run.Data(),slboard));
}

