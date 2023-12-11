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
// #include "../SLBperformance/TBchecks/analysis.h"
#include "../SLBcommissioning/conf_struct.h"

using namespace std;

const Int_t NSLABS = 15;
const Int_t NCHIPS = 16;
const Int_t NCHANNELS = 64;


float FirstDrop(TGraphErrors *gr)
{

  int imax = TMath::LocMax(gr->GetN(), gr->GetY());
  Double_t xmax, ymax;
  Double_t tmpy;
  gr->GetPoint(imax, xmax, ymax);
  tmpy = ymax;

  for (int i = imax; i < gr->GetN(); i++)
  {
    gr->GetPoint(i, xmax, ymax);
    if (ymax < 0.9 * tmpy)
      return xmax;
  }

  return 0;
}

float FirstZero(TGraphErrors *gr)
{

  int imax = TMath::LocMax(gr->GetN(), gr->GetY());
  Double_t xmax, ymax;
  gr->GetPoint(imax, xmax, ymax);

  // if(imax==0) return 0;

  for (int i = imax; i < gr->GetN(); i++)
  {
    gr->GetPoint(i, xmax, ymax);
    if (ymax < 0.05)
      return xmax;
  }

  return 0;
}

// TF1 *FitScurve(TGraphErrors *gr, float xmin_ = 200, float xmax_ = 350, float x0 = 200)
TF1 *FitScurve(TGraphErrors *gr, float xmin_ = 180, float xmax_ = 350, float x0 = 200)
{

  double par1 = 0, par2 = 0, par3 = 0;

  int imax = TMath::LocMax(gr->GetN(), gr->GetY());
  int n = gr->GetN();
  double *x = gr->GetX();
  double *y = gr->GetY();

  if (n == 0 || y[imax] == 0)
    return NULL;

  double xmin = FirstZero(gr) - 15; // max(x[imax],x[n/3])-5;
  double xmax = FirstZero(gr) + 15;
  double ymax = y[imax];

  TF1 *fit1 = new TF1("fit1", "[0]*TMath::Erfc((x-[1])/(sqrt(2)*[2]))", xmin, xmax);

  // fit1->SetParLimits(0, ymax * 0.1, ymax * 10);
  fit1->SetParLimits(0, 0, ymax * 1.5);
  fit1->SetParLimits(1, xmin_, xmax_);

  fit1->SetParameter(0, ymax);
  // fit1->SetParameter(1,250);
  fit1->SetParameter(2, 3);
  gr->Fit("fit1", "QR0");

  par1 = fit1->GetParameter(0);
  par2 = fit1->GetParameter(1);
  par3 = fit1->GetParameter(2);

  xmin = max(par2 - 8 * par3, x[n / 3]);
  xmax = par2 + 20 * par3;
  TF1 *fit2 = new TF1("fit2", "[0]*TMath::Erfc((x-[1])/[2]+[3])", xmin, xmax);

  fit2->SetParLimits(0, par1 * 0.3, par1 + par1 * 0.5);
  fit2->SetParLimits(1, xmin_, xmax_);
  fit2->SetParLimits(2, 1.5, 10);

  fit2->SetParameter(0, par1);
  fit2->SetParameter(1, par2);
  fit2->SetParameter(2, par3);
  gr->Fit("fit2", "QR0");

  return fit2;
}

void drawScurve()
{
  read_configuration_file("../SLBcommissioning/thresholds/Run_Settings_090185.txt", false);

  // TFile *file_scurve = new TFile("../SLBcommissioning/thresholds/histos/scurves_06132022_1MIPinjection_chn0to7.root", "READ");
  // TFile *file_scurve = new TFile("../SLBcommissioning/thresholds/histos/scurves_injection_3p0high_chn0to8.root", "READ");
  TFile *file_scurve = new TFile("../SLBcommissioning/thresholds/histos/scurves_1.2pF.root", "READ");
  // TFile *file_scurve = new TFile("../SLBcommissioning/thresholds/histos/scurves_RateVsThresholdScan_02192020_SLBoard_test2.root", "READ");

// detector.slab[0][islab].asu[0].skiroc[ichip].mask[ichn] == 0

  TCanvas * c_scurve_sample = new TCanvas("c_scurve_sample","c_scurve_sample",800,800);
  c_scurve_sample->SetLeftMargin(0.12);
  c_scurve_sample->SetRightMargin(0.07);
  c_scurve_sample->SetGrid();

  Int_t nchannels = 8;

  TMultiGraph *mg = new TMultiGraph();
  TF1 *fit_scurve_sample[nchannels];

  Color_t color[8] = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kYellow, kBlack};
  TString legend[nchannels];

  for(int ichn=0; ichn < nchannels; ichn++){
  // for(int ichn=10; ichn < 20; ichn++){
    TGraphErrors * h_scurve_sample = (TGraphErrors*) file_scurve->Get(Form("scurve_slbAdd7_chip12_chn%i", ichn));
    // TGraphErrors * h_scurve_sample = (TGraphErrors*) file_scurve->Get(Form("scurve_layer1_chip12_chn%i", ichn));
    fit_scurve_sample[ichn] = FitScurve(h_scurve_sample);
    h_scurve_sample->SetMarkerColor(color[ichn]);
    h_scurve_sample->SetLineColor(color[ichn]);
    h_scurve_sample->SetMarkerStyle(20+ichn);
    h_scurve_sample->SetMarkerSize(1.5);
    fit_scurve_sample[ichn]->SetLineColor(color[ichn]);
    fit_scurve_sample[ichn]->SetLineWidth(2);

    legend[ichn] = Form("Channel %i", ichn);

    mg->Add(h_scurve_sample,"lp");
  }
  mg->SetTitle(";DAC;N Hits");
  gStyle->SetPalette(55);
  mg->GetXaxis()->SetRangeUser(180,310);
  mg->GetYaxis()->SetRangeUser(0,100);
  mg->Draw("a");
  for(int ichn=0; ichn < nchannels; ichn++){
    fit_scurve_sample[ichn]->Draw("same");
  }

  TLegend *leg = new TLegend(0.63,0.63,0.90,0.86);
  leg->SetHeader("#bf{Channels}");
  for(int ichn=0; ichn < nchannels; ichn++){
    leg->AddEntry(fit_scurve_sample[ichn],legend[ichn],"l");
  }
  leg->Draw();


  // TGraphErrors * h_scurve_sample = (TGraphErrors*) file_scurve->Get("scurve_slbAdd0_chip12_chn0");
  // TF1 *fit_scurve_sample = FitScurve(h_scurve_sample);

  // h_scurve_sample->SetTitle(";DAC;N Hits");
  // h_scurve_sample->Draw("");

  // TF1 *scurvefit[15][16][64];

  // for(int islab = 0; islab < NSLABS; islab++){
  //   for(int ichip = 0; ichip < NCHIPS; ichip++){
  //     for(int ichn = 0; ichn < NCHANNELS; ichn++){



  //     }

  //   }

  // }


}

