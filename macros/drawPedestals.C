#include "TFile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TH1.h"
#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include "TF1.h"
#include "Styles.hh"
#include "../SLBcommissioning/conf_struct.h"

using namespace std;

int nslabs = 15;
int nchips = 16;
int nchannels = 64;

int iSCA = 0;

TH2F* h2XY(TH2F* h2, int slab)
{
  TString hname = h2->GetName();
  TH2F *map_XY = new TH2F(hname+"_XY",hname+"_XY", 32, -90, 90, 32, -90, 90);

  TString map_name = "../mapping/fev10_chip_channel_x_y_mapping.txt";
  ReadMap(map_name);



  int nbinsX = h2->GetNbinsX();
  int nbinsY = h2->GetNbinsY();

  for( int ichip=0; ichip<nchips; ichip++ ){
    // cout << "=== chip " << ichip << " ===" << endl;
    for( int ichn=0; ichn<nchannels; ichn++ ){

      int val_x = slab*20  + ichip;
      int val_y = iSCA*100 + ichn;

      int bix_x = h2->GetXaxis()->FindBin(val_x);
      int bix_y = h2->GetYaxis()->FindBin(val_y);

      int bin = h2->GetBin(bix_x, bix_y);

      float ped = h2->GetBinContent(bin);

      float x = map_pointX[ichip][ichn];
      float y = map_pointY[ichip][ichn];

      // cout << x << " " << y << " " << ped << endl;
      map_XY->Fill(x, y, ped);

    }
  }

  map_XY->SetTitle(";x [mm];y [mm]");
  return map_XY;

}

void drawPedestals()
{
  gStyle->SetOptStat(0);

  TFile *file_PedSummary = new TFile("../pedestals/raw_siwecal_90259to90285_highgain_PedSummary.root", "READ");

  TH2F *h2_ped   = (TH2F*)file_PedSummary->Get("ped_all");
  TH2F *h2_width = (TH2F*)file_PedSummary->Get("width_all");

  TH2F* pedXY   = h2XY(h2_ped, 7);
  TH2F* widthXY = h2XY(h2_width, 7);

  TCanvas *c_ped_XY = new TCanvas("c_ped_XY", "c_ped_XY", 800, 800);
  TPad *p_ped_XY    = new TPad("p_ped_XY", "p_ped_XY", 0,0,1,1);
  StylePad(p_ped_XY,0.1,0.15,0.17,0.12);
  pedXY->GetZaxis()->SetRangeUser(180, 300);
  pedXY->GetZaxis()->SetTitle("Pedestal position [ADC]");
  pedXY->GetZaxis()->SetTitleOffset(1.4);
  pedXY->Draw("colz");

  TCanvas *c_width_XY = new TCanvas("c_width_XY", "c_width_XY", 800, 800);
  TPad *p_width_XY    = new TPad("p_width_XY", "p_width_XY", 0,0,1,1);
  StylePad(p_width_XY,0.1,0.15,0.17,0.12);
  widthXY->GetZaxis()->SetRangeUser(1, 3);
  widthXY->GetZaxis()->SetTitle("Pedestal width [ADC]");
  widthXY->GetZaxis()->SetTitleOffset(1.4);
  widthXY->Draw("colz");

  // cout << h2_ped->GetNbinsX() << " " << h2_ped->GetNbinsY() << endl;

}

