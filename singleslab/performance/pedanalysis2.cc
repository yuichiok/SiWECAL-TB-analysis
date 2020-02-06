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

void pedanalysis2(){
  
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
  
  // Comparing nbr entries in tag or not tag // Get width and Mean of the pedestals

  // Defining variables [chip][chn][sca]
  TH1F *G[10][100][20];
  TH1F *T[10][100][20];
  

  TCanvas* canvas[4];

  canvas[0]= new TCanvas("run_32001","canvas1",1200,600);
  canvas[1]= new TCanvas("run_32002","canvas2",1200,600);
  canvas[2]= new TCanvas("run_32003","canvas3",1200,600);
  canvas[3]= new TCanvas("run_32004","canvas4",1200,600);
  canvas[4]= new TCanvas("run_32005","canvas5",1200,600);
  Double_t Mean[10][100][20];
  Double_t Width[10][100][20];
  Double_t Ge[10][100][20]; //GoodEntries
  Double_t Te[10][100][20]; //BadEntries
  Double_t q; //BadEntries
  Double_t p; //GoodEntries
  Double_t o; //Mean
  Double_t r; //Width
  Double_t x[400];
  Double_t y[10][100][20];
  Double_t m; //Mean

  //GoodEntries,TaggedEntries,Ratio,Width,Mean
  TH2F *Gentries=new TH2F("Gentries","Gentries",400,0.5,400.5,15,0.5,15.5);
  TH2F *Tentries=new TH2F("Tentries","Tentries",400,0.5,400.5,15,0.5,15.5);
  TH2F *Ratio=new TH2F("Ratio","Ratio",400,0.5,400.5,15,0.5,15.5);
  TH2F *W=new TH2F("PedWidth","PedWidth",400,0.5,400.5,15,0.5,15.5);
  TH2F *M=new TH2F("PedMean","PedMean",400,0.5,400.5,15,0.5,15.5);
  
  for(int i=0;i<5;i++){
    for(int j=0; j<65; j++) {
      for(int n=0; n<17;n++) {
	G[i][j][n]=new TH1F(TString::Format("G_chip_%i_chn_%i_sca%i",i,j,n),TString::Format("G_chip_%i_chn_%i_sca%i",i,j,n),1000,0.5,1000.5);
	T[i][j][n]=new TH1F(TString::Format("T_chip_%i_chn_%i_sca%i",i,j,n),TString::Format("T_chip_%i_chn_%i_sca%i",i,j,n),1000,0.5,1000.5);
      }
    }
  }
  
  //Reading the file and extract histograms and data
  cout << "Reading the files/extracting data"<< endl;
  for(int ifile=0; ifile<1; ifile++) { // This first loop on reading different run is not used actually
    TFile *_file2 = TFile::Open(TString::Format("../Pedestal_SLB_0_run_%i.root",ifile+21007));
    cout<<ifile<<" "<<TString::Format("../Pedestal_SLB_0_run_%i.root",ifile+21007)<<endl;
    TFile *_file1 = TFile::Open(TString::Format("../Pedestal_SLB_0_run_%i.root",ifile+21012));
    cout<<ifile<<" "<<TString::Format("../Pedestal_SLB_0_run_%i.root",ifile+21012)<<endl;
    
    for(int i=0;i<4;i++){
      for(int j=0; j<64; j++) {
	for(int n=0; n<15;n++) {
	  //GetGoodEntries and GetMean
	  G[i][j][n]=(TH1F*)_file2->Get(TString::Format("ped_chip%i_chn%i_sca%i",i,j,n));
	  Ge[i][j][n]=G[i][j][n]->GetEntries();
	  Mean[i][j][n]=G[i][j][n]->GetMean();
	  //GetWidth by fitting arround the pedestal
	  m=Mean[i][j][n];
	  TF1 *f20 = new TF1("f20","gaus",m-15,m+15);
	  G[i][j][n]->Fit("f20","RQNOC");
	  Width[i][j][n]=f20->GetParameter(2);
	  //GetBadEntries
	  T[i][j][n]=(TH1F*)_file2->Get(TString::Format("ped_tagged_chip%i_chn%i_sca%i",i,j,n));
	  Te[i][j][n]=T[i][j][n]->GetEntries();
	}
      }
    }
    //Filling histograms
    cout << "Filling histograms"  << endl;
    for(int i=0;i<15;i++){
      for(int j=0;j<64;j++){
	for(int n=0;n<4;n++){
	  p=Ge[n][j][i]; //GoodEntries
	  q=Te[n][j][i]; //BadEntries
	  r=Width[n][j][i]; //Width
	  o=Mean[n][j][i]; //Mean
	  W->Fill(j+100*n,i,r);
	  M->Fill(j+100*n,i,o);
	  Gentries->Fill(j+100*n,i, p);
	  Tentries->Fill(j+100*n,i, q);
	  Ratio->Fill(j+100*n,i, (p-q)/(p+q));	  
	}
      }
    }
    
    cout << "Drawing histograms"  << endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////All graphs////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    canvas[0]->cd();
    Gentries->Draw("COLZ");
    canvas[1]->cd();
    Tentries->Draw("COLZ");
    canvas[2]->cd();
    Ratio->Draw("COLZ");
    canvas[3]->cd();
    M->GetZaxis()->SetRangeUser(200,350);
    M->Draw("COLZ");
    canvas[4]->cd();
    W->GetZaxis()->SetRangeUser(2,6);
    W->Draw("COLZ");
  }
}

