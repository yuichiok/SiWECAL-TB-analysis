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

void pedanalysis0(int run=0, int ifile=0){
  
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

  //Defining variables [chip][chn][sca]
  TH1F *G[10][100][20];
  TH1F *T[10][100][20];
  

  TCanvas* canvas;
  canvas= new TCanvas(TString::Format("PedAna_run%i_SLB_%i",run,ifile),TString::Format("PedAna_run%i_SLB_%i",run,ifile),1600,750);   
  canvas->Divide(3,2);
  Double_t Mean[10][100][20];
  Double_t Width[10][100][20];
  Double_t Ge[10][100][20];
  Double_t Te[10][100][20];
  Double_t q;
  Double_t p;
  Double_t o;
  Double_t r;
  Double_t x[400];
  Double_t y[10][100][20];
  Double_t m;
  TH2F *Gentries;
  TH2F *Tentries;
  TH2F *Ratio;
  TH2F *W;
  TH2F *M;
  
  //defining the 2Dhisto
  Gentries=new TH2F("Gentries","Gentries",400,0.5,400.5,15,0.5,15.5);
  Tentries=new TH2F("Tentries","Tentries",400,0.5,400.5,15,0.5,15.5);
  Ratio=new TH2F("Ratio","Ratio",400,0.5,400.5,15,0.5,15.5);
  W=new TH2F("Pedestal Width","Pedestal Width",400,0.5,400.5,15,0.5,15.5);
  M=new TH2F("Pedestal Mean","Pedestal Mean",400,0.5,400.5,15,0.5,15.5);
    for(int i=0;i<5;i++){
      for(int j=0; j<65; j++) {
	for(int n=0; n<17;n++) {
	  G[i][j][n]=new TH1F(TString::Format("G_chip_%i_chn_%i_sca%i",i,j,n),TString::Format("G_chip_%i_chn_%i_sca%i",i,j,n),1000,0.5,1000.5);
	  T[i][j][n]=new TH1F(TString::Format("T_chip_%i_chn_%i_sca%i",i,j,n),TString::Format("T_chip_%i_chn_%i_sca%i",i,j,n),1000,0.5,1000.5);
	}
      }
    }
    
    cout << "Reading files . . .Extracting data"  << endl;
    TFile *_file0 = TFile::Open(TString::Format("../Pedestal_SLB_%i_run_%i.root",ifile,run));
    cout<<ifile<<" "<<TString::Format("../Pedestal_SLB_%i_run_%i.root",ifile,run)<<endl;
    
    
    for(int i=0;i<4;i++){
      for(int j=0; j<64; j++) {
	for(int n=0; n<15;n++) {
	  //GetGoodEntries
	  G[i][j][n]=(TH1F*)_file0->Get(TString::Format("ped_chip%i_chn%i_sca%i",i,j,n));
	  Ge[i][j][n]=G[i][j][n]->GetEntries();
	  //GetMean
	  Mean[i][j][n]=G[i][j][n]->GetMean();
	  //GetWidth
	  m=Mean[i][j][n];
	  TF1 *f20 = new TF1("f20","gaus",m-15,m+15);
	  G[i][j][n]->Fit("f20","RQNOC");
	  Width[i][j][n]=f20->GetParameter(2);
	  //GetBadEntries
	  T[i][j][n]=(TH1F*)_file0->Get(TString::Format("ped_tagged_chip%i_chn%i_sca%i",i,j,n));
	  Te[i][j][n]=T[i][j][n]->GetEntries();
	}
      }
    }
    cout << "Filling histograms"  << endl;
    
    //filling histograms
    for(int i=0;i<15;i++){
      for(int j=0;j<64;j++){
	for(int n=0;n<4;n++){
	  p=Ge[n][j][i];
	  q=Te[n][j][i];
	  r=Width[n][j][i];
	  o=Mean[n][j][i];
	  W->Fill(j+100*n,i,r);
	  M->Fill(j+100*n,i,o);
	  Gentries->Fill(j+100*n,i, p);
	  Tentries->Fill(j+100*n,i, q);
	  Ratio->Fill(j+100*n,i, (p-q)/(p+q));
 	}
      }
    }
    
    cout << "Drawing histograms" << endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////All graphs////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    canvas->cd(1);
    Gentries->GetXaxis()->SetTitle("chn+chip*100");
    Gentries->GetYaxis()->SetTitle("SCA");
    Gentries->Draw("COLZ");
    canvas->cd(2);
    Tentries->GetXaxis()->SetTitle("chn+chip*100");
    Tentries->GetYaxis()->SetTitle("SCA");
    Tentries->Draw("COLZ");
    canvas->cd(3);
    Ratio->GetXaxis()->SetTitle("chn+chip*100");
    Ratio->GetYaxis()->SetTitle("SCA");
    Ratio->Draw("COLZ");
    canvas->cd(4);
    M->GetXaxis()->SetTitle("chn+chip*100");
    M->GetYaxis()->SetTitle("SCA");
    M->GetZaxis()->SetRangeUser(200,350);
    M->Draw("COLZ");
    canvas->cd(5);
    W->GetXaxis()->SetTitle("chn+chip*100");
    W->GetYaxis()->SetTitle("SCA");
    W->GetZaxis()->SetRangeUser(2,6);
    W->Draw("COLZ");
    canvas->Print(TString::Format("plots/PedAna_run%i_SLB_%i.C",run,ifile));
    canvas->Print(TString::Format("plots/PedAna_run%i_SLB_%i.eps",run,ifile));
    
    _file0 ->Close();
  }

