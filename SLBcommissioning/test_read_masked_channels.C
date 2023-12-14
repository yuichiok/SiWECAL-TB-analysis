//# Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)

#include "conf_struct.h"

void test_read_masked_channels(TString filename = "15102021/Run_Settings_DataTaking_TB2021_15102021", bool debug = true)
{
  read_configuration_file(filename, true);

  TH2F *mask_chip_chn[15];
  TH2F *mask_x_y[15];

  float totalmasked[15];
  float totalmasked_chip[15][16];

  int mapping_slab[15];
  mapping_slab[0] = 34;
  mapping_slab[1] = 36;
  mapping_slab[2] = 38;
  mapping_slab[3] = 35; //
  mapping_slab[4] = 19;
  mapping_slab[5] = 33;
  mapping_slab[6] = 29;
  mapping_slab[7] = 30; //
  mapping_slab[8] = 31;
  mapping_slab[9] = 24;  //
  mapping_slab[10] = 20; //
  mapping_slab[11] = 37;
  mapping_slab[12] = 18; //
  mapping_slab[13] = 23;
  mapping_slab[14] = 17;

  ofstream fout(filename + "_masked.txt", ios::out);
  fout << "#masked_chns_list layer chip chns (0=not masked, 1=masked)" << endl;

  int nslabs = 15;
  for (int islab = 0; islab < nslabs; islab++)
  {
    TString map_name = "../mapping/fev10_chip_channel_x_y_mapping.txt";
    // TString map_name="../mapping/fev11_cob_rotate_chip_channel_x_y_mapping.txt";

    // the two cobs are equipped with slboards 2.08 and 2.12 (26th May 2020)
       if(detector.slab[0][islab].add==6 || detector.slab[0][islab].add==8) 
          map_name="../mapping/fev11_cob_rotate_chip_channel_x_y_mapping.txt";

    cout << " -------------------------------------------------------------------------------------------- " << endl;
    cout << "Slab idx" << islab << "with slabAdd: " << detector.slab[0][islab].add << " and mapping file: " << map_name << endl;

    ReadMap(map_name);

    mask_chip_chn[islab] = new TH2F(TString::Format("mask_chip_chn_SLB%i", detector.slab[0][islab].add), TString::Format("mask_chip_chn_SLB%i", detector.slab[0][islab].add), 16, -0.5, 15.5, 64, -0.5, 63.5);
    mask_x_y[islab] = new TH2F(TString::Format("mask_x_y_SLB%i", detector.slab[0][islab].add), TString::Format("mask_x_y_SLB%i", detector.slab[0][islab].add), 32, -90, 90, 32, -90, 90);

    totalmasked[islab] = 0;
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

        float msk = detector.slab[0][islab].asu[0].skiroc[ichip].mask[ichn] + 0.001;
        mask_chip_chn[islab]->Fill(ichip, ichn, msk);
        mask_x_y[islab]->Fill(map_pointX[ichip][ichn], map_pointY[ichip][ichn], msk);
        fout << " " << detector.slab[0][islab].asu[0].skiroc[ichip].mask[ichn];
      }
      fout << endl;
      cout << "Skiroc idx:" << ichip << "  Masked cells:" << 100. * totalmasked_chip[islab][ichip] / 64. << "%" << endl;
    }
    cout << "TOTAL Masked cells:" << 100. * totalmasked[islab] / 1024. << "%" << endl;

    TCanvas *canvas = new TCanvas(TString::Format("canvas_%i", islab), TString::Format("canvas_%i", islab), 500, 1000);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kInvertedDarkBodyRadiator); // Cherry);
    // TColor::InvertPalette();
    canvas->Divide(1, 2);
    canvas->cd(1);
    // mask_chip_chn[islab]->SetTitle(TString::Format("SLABd%i",mapping_slab[islab]));
    mask_chip_chn[islab]->SetTitle(TString::Format("Layer %i Mask cells Chip vs. Channel", islab));
    mask_chip_chn[islab]->GetXaxis()->SetTitle("CHIP");
    mask_chip_chn[islab]->GetYaxis()->SetTitle("CHANNEL");
    mask_chip_chn[islab]->Draw("col");
    canvas->cd(2);
    mask_x_y[islab]->SetTitle(TString::Format("Layer %i Mask cells XY positions",islab));
    // mask_x_y[islab]->SetTitle(TString::Format("SLABd%i",mapping_slab[islab]));
    // mask_x_y[islab]->SetTitle(TString::Format("SlabAdd%i", islab));
    mask_x_y[islab]->GetXaxis()->SetTitle("x");
    mask_x_y[islab]->GetYaxis()->SetTitle("y");
    mask_x_y[islab]->Draw("col");
    // canvas->Print(TString::Format("masked_channels_SLAB%i.png",mapping_slab[islab]));
    canvas->Print(TString::Format("plots/masked_channels_SlabAdd%i.pdf", islab));

  }

  TCanvas *canvas_chn_chip = new TCanvas("canvas_chn_chip", "canvas_chn_chip", 800, 800);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kInvertedDarkBodyRadiator); // Cherry);
  mask_chip_chn[1]->SetTitle("Layer 1 Mask cells Chip vs. Channel");
  mask_chip_chn[1]->GetXaxis()->SetTitle("CHIP");
  mask_chip_chn[1]->GetYaxis()->SetTitle("CHANNEL");
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.09);
  mask_chip_chn[1]->Draw("col");
  canvas_chn_chip->Print("plots/mask_ChipChn_slab1.pdf");

  TCanvas *canvas_x_y = new TCanvas("canvas_x_y", "canvas_x_y", 800, 800);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kInvertedDarkBodyRadiator); // Cherry);
  mask_x_y[1]->SetTitle("Layer 1 Mask cells XY positions");
  mask_x_y[1]->GetXaxis()->SetTitle("x [mm]");
  mask_x_y[1]->GetYaxis()->SetTitle("y [mm]");
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.09);
  mask_x_y[1]->Draw("col");
  canvas_x_y->Print("plots/mask_xy_slab1.pdf");


  cout << "################################################################" << endl;
  cout << "################################################################" << endl;
  cout << "####           CHECK THE MASK MAP FILE !!!!!!               ####" << endl;
  cout << "################################################################" << endl;
  cout << "################################################################" << endl;
}
