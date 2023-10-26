// # Copyright 2020  Adrián Irles

#include "TROOT.h"
#include "TFile.h"
#include "DecodedSLBAnalysis.cc"
#include "TGraphErrors.h"
#include "TLegend.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"

void LinearityGraphs(TString filename_in = "/home/calice/TB2020/commissioning/SiWECAL-TB-analysis_TB2020/converter_SLB/convertedfiles/",
                     TString run = "calib_02272020_pa1.2fb_trig230_chn0.8.16.24_calib_small_Ascii", int nslboards = 6)
{

  cout << " Injection folder: " << filename_in << endl;

  // slabs, chips, channels, sca, size --> yhigh, eyhigh, ylow, eylow
  std::vector<std::vector<std::array<float, 4>>> y;

  float x[100];
  float ex[100];

  for (int l = 0; l < 100; l++)
  {
    x[l] = 0;
    ex[l] = 0;
  }

  for (int i = 0; i < 15; i++)
  {
    for (int j = 0; j < 16; j++)
    {
      for (int k = 0; k < 64; k++)
      {
        for (int isca = 0; isca < 3; isca++)
        {
        std:
          vector<std::array<float, 4>> tempvec;
          for (int l = 0; l < 100; l++)
          {
            std::array<float, 4> temp;
            temp[0] = 0;
            temp[1] = 0;
            temp[2] = 0;
            temp[3] = 0;
            tempvec.push_back(temp);
          }
          y.push_back(tempvec);
        }
      }
    }
  }

  int board_id[15];

  for (int i = 0; i < 15; i++)
  {
    board_id[i] = -1;
  }

  int count = 0;
  for (int i = 0; i < 15; i++)
  {
    float pt = 0.0;
    if (i == 0)
      pt = 0.0;
    if (i == 1)
      pt = 0.1;
    if (i == 2)
      pt = 0.2;
    if (i == 3)
      pt = 0.3;
    if (i == 4)
      pt = 0.4;
    if (i == 5)
      pt = 0.5;
    if (i == 6)
      pt = 0.65;
    if (i == 7)
      pt = 0.7;
    if (i == 8)
      pt = 0.8;
    if (i == 9)
      pt = 0.9;
    if (i == 10)
      pt = 1.0;
    if (i == 11)
      pt = 1.15;
    if (i == 12)
      pt = 1.3;
    if (i == 13)
      pt = 1.5;
    if (i == 14)
      pt = 1.7;

    TString filename = filename_in + TString::Format("%s/injection_small_%.1fV.root", run.Data(), pt);
    if (i == 6 || i == 11)
      filename = filename_in + TString::Format("%s/injection_small_%.2fV.root", run.Data(), pt);

    cout << " i: " << i << endl;
    cout << " Injection file: " << filename_in << endl;

    DecodedSLBAnalysis ss(filename);
    std::vector<std::array<float, 9>> injecvalues = ss.InjectionAnalysis(-1);

    x[count] = (3.3 - pt - 1.408); // 2*(1-(pt-0.65));//(3.3-pt)/2.15;

    for (int j = 0; j < injecvalues.size(); j++)
    {
      if (board_id[int(injecvalues.at(j)[0])] == -1)
        board_id[int(injecvalues.at(j)[0])] = int(injecvalues.at(j)[1]);

      y.at(j).at(i)[0] = injecvalues.at(j)[5]; // yhigh
      y.at(j).at(i)[1] = injecvalues.at(j)[6]; // eyhigh

      y.at(j).at(i)[2] = injecvalues.at(j)[7]; // ylow
      y.at(j).at(i)[3] = injecvalues.at(j)[8]; // eylow
    }
    count++;
  }

  TGraphErrors *injection_high_sca0[15][16][64];
  TGraphErrors *injection_high_sca1[15][16][64];
  TGraphErrors *injection_high_sca2[15][16][64];

  TGraphErrors *injection_low_sca0[6][16][64];
  TGraphErrors *injection_low_sca1[6][16][64];
  TGraphErrors *injection_low_sca2[6][16][64];

  int counter = 0;
  for (int i = 0; i < nslboards; i++)
  {
    for (int j = 0; j < 16; j++)
    {
      for (int k = 0; k < 64; k++)
      {
        for (int isca = 0; isca < 3; isca++)
        {

          float y0[100], ey0[100];
          float y1[100], ey1[100];

          for (int ivec = 0; ivec < count; ivec++)
          {
            y0[ivec] = y.at(counter).at(ivec)[0];
            ey0[ivec] = y.at(counter).at(ivec)[1];
            y1[ivec] = y.at(counter).at(ivec)[2];
            ey1[ivec] = y.at(counter).at(ivec)[3];
          }

          counter++;

          if (isca == 0)
          {
            injection_high_sca0[i][j][k] = new TGraphErrors(count, x, y0, ex, ey0);
            injection_low_sca0[i][j][k] = new TGraphErrors(count, x, y1, ex, ey1);
          }
          if (isca == 1)
          {
            injection_high_sca1[i][j][k] = new TGraphErrors(count, x, y0, ex, ey0);
            injection_low_sca1[i][j][k] = new TGraphErrors(count, x, y1, ex, ey1);
          }
          if (isca == 2)
          {
            injection_high_sca2[i][j][k] = new TGraphErrors(count, x, y0, ex, ey0);
            injection_low_sca2[i][j][k] = new TGraphErrors(count, x, y1, ex, ey1);
          }
        }
      }
    }
  }

  TFile *file_summary = new TFile("results/graphs_" + run + ".root", "RECREATE");
  file_summary->cd();

  TH1F *h1 = new TH1F("layer_slboard_relation", "layer_slboard_relation", nslboards, -0.5, nslboards - 0.5);
  for (int i = 0; i < nslboards; i++)
  {
    h1->Fill(i, board_id[i]);
  }
  h1->Write();

  for (int i = 0; i < nslboards; i++)
  {
    for (int j = 0; j < 16; j++)
    {
      for (int k = 0; k < 64; k++)
      {
        injection_high_sca0[i][j][k]->SetName(TString::Format("high_gain_layer%i_slboard%i_chip%i_chn%i_sca0", i, board_id[i], j, k));
        injection_high_sca1[i][j][k]->SetName(TString::Format("high_gain_layer%i_slboard%i_chip%i_chn%i_sca1", i, board_id[i], j, k));
        injection_high_sca2[i][j][k]->SetName(TString::Format("high_gain_layer%i_slboard%i_chip%i_chn%i_sca2", i, board_id[i], j, k));
        injection_low_sca0[i][j][k]->SetName(TString::Format("low_gain_layer%i_slboard%i_chip%i_chn%i_sca0", i, board_id[i], j, k));
        injection_low_sca1[i][j][k]->SetName(TString::Format("low_gain_layer%i_slboard%i_chip%i_chn%i_sca1", i, board_id[i], j, k));
        injection_low_sca2[i][j][k]->SetName(TString::Format("low_gain_layer%i_slboard%i_chip%i_chn%i_sca2", i, board_id[i], j, k));

        injection_high_sca0[i][j][k]->Write();
        injection_high_sca1[i][j][k]->Write();
        injection_high_sca2[i][j][k]->Write();
        injection_low_sca0[i][j][k]->Write();
        injection_low_sca1[i][j][k]->Write();
        injection_low_sca2[i][j][k]->Write();
      }
    }
  }
}
