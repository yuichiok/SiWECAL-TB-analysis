// # Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)

#include "conf_struct.h"
#include "../include/style.h"

void test_read_masked_channels_TB(TString filename = "/mnt/HardDrive/beamData/ascii/1.4GeV_W_run_050163/Run_Settings.txt", TString filename_out = "test")
{

  read_configuration_file(filename, true);

  TH2F *mask_chip_chn = new TH2F("mask_chip_chn", "mask_chip_chn", 16, -0.5, 15.5, 64, -0.5, 63.5);
  TH2F *mask_x_y = new TH2F("mask_x_y", "mask_x_y", 32, -90, 90, 32, -90, 90);

  TH2F *mask_layer_chip = new TH2F("mask_layer_chip", "mask_layer_chip", 15, -0.5, 14.5, 16, -0.5, 15.5);
  TH2F *mask_layer = new TH2F("mask_layer", "mask_layer", 15, -0.5, 14.5, 20, -0.5, 20.5);
  TH1F *mask_rate = new TH1F("mask_rate", "mask_rate", 15, -0.5, 14.5);

  float totalmasked[15];
  float totalmasked_chip[15][16];

  float tmpaverage = 0;

  ofstream fout(filename_out + ".txt", ios::out);
  fout << "#masked_chns_list layer chip chns (0=not masked, 1=masked)" << endl;
  int nslabs = 15;
  for (int islab = 0; islab < nslabs; islab++)
  {
    TString map_name = "../mapping/fev10_chip_channel_x_y_mapping.txt";

    ReadMap(map_name);

    totalmasked[islab] = 0;

    cout << " ---------------------------------- "<< islab <<" ----------------------------------------- " << endl;

    for (int ichip = 0; ichip < 16; ichip++)
    {
      totalmasked_chip[islab][ichip] = 0;
      fout << islab << " " << ichip;
      for (int ichn = 0; ichn < 64; ichn++)
      {
        if (detector.slab[0][islab].asu[0].skiroc[ichip].mask[ichn] == 1)
        {
          totalmasked[islab]++;
          totalmasked_chip[islab][ichip]++;
        }

        float msk = detector.slab[0][islab].asu[0].skiroc[ichip].mask[ichn];
        if (msk > 0)
        {
          mask_chip_chn->Fill(ichip, ichn);
          mask_x_y->Fill(map_pointX[ichip][ichn], map_pointY[ichip][ichn]);
        }
        fout << " " << detector.slab[0][islab].asu[0].skiroc[ichip].mask[ichn];
      }
      fout << endl;
      cout << "Skiroc idx:" << ichip << "  Masked cells:" << 100. * totalmasked_chip[islab][ichip] / 64. << "%" << endl;
      mask_layer_chip->Fill(14 - islab, ichip, totalmasked_chip[islab][ichip] + 0.01);
    }
    cout << "TOTAL Masked cells:" << 100. * totalmasked[islab] / 1024. << "%" << endl;
    mask_layer->Fill(14 - islab, 100. * totalmasked[islab] / 1024.);
    // mask_rate->Fill(islab, 100. * totalmasked[islab] / 1024.);

    Float_t mask12_15 = 0;
    for(int ichip=12; ichip<16; ichip++){
      mask12_15 += totalmasked_chip[islab][ichip];
    }
    mask12_15 = 100. * mask12_15 / 256.;
    cout << "Masked cells in chips 12-15: " << mask12_15 << "%" << endl;
    tmpaverage += mask12_15;
    mask_rate->Fill(islab, mask12_15);

  }

  TCanvas *canvas = new TCanvas("canvas_accum", "canvas_accum", 800, 800);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow); // Cherry);
  gStyle->SetPadRightMargin(0.16);

  canvas->Divide(2, 2);
  canvas->cd(1);
  mask_chip_chn->GetXaxis()->SetTitle("CHIP");
  mask_chip_chn->GetYaxis()->SetTitle("CHANNEL");
  mask_chip_chn->GetZaxis()->SetTitle("N-slabs");
  mask_chip_chn->Draw("colz");

  canvas->cd(2);
  mask_x_y->GetXaxis()->SetTitle("x");
  mask_x_y->GetYaxis()->SetTitle("y");
  mask_x_y->GetZaxis()->SetTitle("N-slabs");
  mask_x_y->Draw("colz");

  canvas->cd(3);
  mask_layer_chip->GetXaxis()->SetTitle("Layer");
  mask_layer_chip->GetYaxis()->SetTitle("CHIP");
  mask_layer_chip->GetZaxis()->SetTitle("masked channels");
  mask_layer_chip->GetZaxis()->SetRangeUser(0, 35);
  mask_layer_chip->Draw("colz");

  canvas->cd(4);
  mask_rate->SetMarkerStyle(20);
  mask_rate->SetBarWidth(0.9);
  mask_rate->SetFillColor(49);
  mask_rate->SetTitle("Masked channels rate");
  mask_rate->GetXaxis()->SetTitle("Layer");
  mask_rate->GetYaxis()->SetTitle("[%] of channels that are masked");
  mask_rate->GetYaxis()->SetRangeUser(0, 100);
  mask_rate->Draw("HIST bar2");


  TCanvas *canvas_mask_xy = new TCanvas("canvas_mask_xy", "canvas_mask_xy", 800, 800);
  TPad *p_mask_xy = new TPad("p_mask_xy", "p_mask_xy", 0, 0, 1, 1);
  StylePad(p_mask_xy,0.1,0.1,0.15,0.12);
  mask_x_y->SetTitle("Mask cells XY positions");
  mask_x_y->Draw("colz");
  p_mask_xy->Draw();

  TCanvas *canvas_rate = new TCanvas("canvas_rate", "canvas_rate", 800, 800);
  TPad *p_rate = new TPad("p_rate", "p_rate", 0, 0, 1, 1);
  StylePad(p_rate,0.1,0.1,0.1,0.17);
  p_rate->SetGrid(0, 1);
  mask_rate->Draw("HIST bar2");
  p_rate->Draw();

  cout << "mean mask rate chips 12-15: " << tmpaverage / 15. << endl;




  // mask_layer->GetXaxis()->SetTitle("Layer");
  // mask_layer->GetYaxis()->SetTitle("[%] of channels that are masked");
  // mask_layer->GetZaxis()->SetTitle("");
  // mask_layer->GetZaxis()->SetRangeUser(0,35);
  // mask_layer->Draw("col");

  // canvas->Print(filename_out + ".root");
  // canvas->Print(filename_out + ".png");
}
