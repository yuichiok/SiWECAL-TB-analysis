// # Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)

#include "conf_struct.h"

void test_read_thresholds_summary(TString filename = "15102021/Run_Settings_DataTaking_TB2021_15102021_WRITENbyDAQ.txt", bool debug = true)
{
  // void test_read_thresholds_summary(TString filename="15102021/Run_Settings_it10.txt", bool debug=true) {

  read_configuration_file(filename, false);

  TH2F *threshold_chip_chn[15];
  TH2F *threshold_x_y[15];

  TH2F *threshold_chip_chn_2[15];
  TH2F *threshold_x_y_2[15];

  TH2F *threshold_layer_chip = new TH2F("threshold_chip_chn", "threshold_chip_chn", 23, 12.5, 35.5, 16, -0.5, 15.5);
  int mapping_slab[15];
  mapping_slab[0] = 18;
  mapping_slab[1] = 23;
  mapping_slab[2] = 17;
  mapping_slab[3] = 22; //
  mapping_slab[4] = 25;
  mapping_slab[5] = 24;
  mapping_slab[6] = 31;
  mapping_slab[7] = 30; //
  mapping_slab[8] = 21;
  mapping_slab[9] = 20;  //
  mapping_slab[10] = 19; //
  mapping_slab[11] = 15;
  mapping_slab[12] = 14; //
  mapping_slab[13] = 13;
  mapping_slab[14] = 16;

  vector<int> thresholdVec;

  for (int islab = 0; islab < 15; islab++)
  {
    TString map_name = "../mapping/fev10_chip_channel_x_y_mapping.txt";

    // the two cobs are equipped with slbAdds 2.08 and 2.12 (26th May 2020)
    if (detector.slab[0][islab].add == 6 || detector.slab[0][islab].add == 8)
      map_name = "../mapping/fev11_cob_rotate_chip_channel_x_y_mapping.txt";

    ReadMap(map_name);

    cout << " ----------------------------------- " << endl;
    cout << "SlbAdd: " << islab << endl;

    threshold_chip_chn[islab] = new TH2F(TString::Format("threshold_chip_chn_%i", islab), TString::Format("threshold_chip_chn_%i", islab), 16, -0.5, 15.5, 64, -0.5, 63.5);
    threshold_x_y[islab] = new TH2F(TString::Format("threshold_x_y_%i", islab), TString::Format("threshold_x_y_%i", islab), 32, -90, 90, 32, -90, 90);

    threshold_chip_chn_2[islab] = new TH2F(TString::Format("threshold_chip_chn_global_%i", islab), TString::Format("threshold_chip_chn_global_%i", islab), 16, -0.5, 15.5, 64, -0.5, 63.5);
    threshold_x_y_2[islab] = new TH2F(TString::Format("threshold_x_y_global_%i", islab), TString::Format("threshold_x_y_global_%i", islab), 32, -90, 90, 32, -90, 90);


    for (int ichip = 0; ichip < 16; ichip++)
    {
      cout << "   chip: " << ichip << "  global: " << detector.slab[0][islab].asu[0].skiroc[ichip].threshold_dac << " DAC" << endl;
      cout << "   individual: ";
      thresholdVec.push_back(detector.slab[0][islab].asu[0].skiroc[ichip].threshold_dac);

      for (int ichn = 0; ichn < 64; ichn++)
      {
        float th = -detector.slab[0][islab].asu[0].skiroc[ichip].chn_threshold_adj[ichn];
        float th_glob = detector.slab[0][islab].asu[0].skiroc[ichip].threshold_dac;

        cout << th << " ";
        if (detector.slab[0][islab].asu[0].skiroc[ichip].mask[ichn] == 0)
        {
          threshold_chip_chn[islab]->Fill(ichip, ichn, th);
          threshold_x_y[islab]->Fill(map_pointX[ichip][ichn], map_pointY[ichip][ichn], th);

          threshold_chip_chn_2[islab]->Fill(ichip, ichn, th_glob);
          threshold_x_y_2[islab]->Fill(map_pointX[ichip][ichn], map_pointY[ichip][ichn], th_glob);
        }
      }
      if (mapping_slab[islab] > 23)
        threshold_layer_chip->Fill(mapping_slab[islab], ichip, detector.slab[0][islab].asu[0].skiroc[ichip].threshold_dac - 15);
      else
        threshold_layer_chip->Fill(mapping_slab[islab], ichip, detector.slab[0][islab].asu[0].skiroc[ichip].threshold_dac);
      cout << endl;
    }

    TCanvas *canvas = new TCanvas(TString::Format("canvas_%i", islab), TString::Format("canvas_%i", islab), 800, 1600);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kLightTemperature); // RainBow);
    canvas->Divide(1, 2);
    canvas->cd(1);
    threshold_chip_chn_2[islab]->GetXaxis()->SetTitle("CHIP");
    threshold_chip_chn_2[islab]->GetYaxis()->SetTitle("CHANNEL");
    threshold_chip_chn_2[islab]->GetZaxis()->SetRangeUser(190, 240);
    threshold_chip_chn_2[islab]->Draw("colz");
    threshold_chip_chn[islab]->Draw("text0same");

    canvas->cd(2);
    threshold_x_y_2[islab]->GetXaxis()->SetTitle("x");
    threshold_x_y_2[islab]->GetYaxis()->SetTitle("y");
    threshold_x_y_2[islab]->GetZaxis()->SetRangeUser(190, 240);
    threshold_x_y_2[islab]->Draw("colz");
    threshold_x_y[islab]->Draw("text0same");
    canvas->Print(TString::Format("plots/thresholds_slbAdd%i.eps", islab));

  }

  const auto [min, max] = std::minmax_element(begin(thresholdVec), end(thresholdVec));
  std::cout << "min = " << *min << ", max = " << *max << '\n';


  TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 600, 600);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow); // Cherry);
  gStyle->SetPadRightMargin(0.16);

  canvas2->cd(1);
  threshold_layer_chip->GetXaxis()->SetTitle("SLAB-ID");
  threshold_layer_chip->GetYaxis()->SetTitle("CHIP");
  threshold_layer_chip->GetZaxis()->SetTitle("Global Threshold");
  threshold_layer_chip->GetZaxis()->SetRangeUser(200, 280);
  threshold_layer_chip->Draw("colz");
}
